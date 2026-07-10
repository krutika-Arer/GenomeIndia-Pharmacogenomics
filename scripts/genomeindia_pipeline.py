#!/usr/bin/env python3
"""
Genome India × ClinPGx Pharmacogenomics Pipeline
=================================================
Production-grade, streaming, multiprocessing, resume-safe.

Architecture:
  - Phase 1  : Extract rsIDs from ClinPGx variants.tsv (in-memory; only 7562 rows)
  - Phase 2  : 12 chromosome workers stream VEP TSVs line-by-line
               → write matched rows to GenomeIndia_ClinPGx_Matched.tsv (append, locked)
               → accumulate novel variant counts per chrom (returned via multiprocessing)
  - Phase 3  : Write ClinPGx_NotFound.tsv from unmatched rsIDs
  - Phase 4  : Write novel VCF files per chromosome (stream original VCF, filter)
  - Phase 5  : Write Summary.txt

Usage:
  python genomeindia_pipeline.py \\
    --variants  variants.tsv \\
    --vcf-dir   vcf/ \\
    --vep-dir   vep_pgx/ \\
    --output    output/ \\
    [--workers  12] \\
    [--resume]          # skip chromosomes whose novel VCF already exists
"""

import os
import sys
import csv
import gzip
import logging
import argparse
import fcntl
import time
import multiprocessing
from pathlib import Path
from typing import Dict, Set, Tuple, List, Optional

# ─────────────────────────────────────────────────────────────────────────────
# LOGGING
# ─────────────────────────────────────────────────────────────────────────────

def setup_logging(log_path: Path) -> logging.Logger:
    logger = logging.getLogger("pgx_pipeline")
    logger.setLevel(logging.DEBUG)
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")

    fh = logging.FileHandler(log_path, mode="a")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(fmt)

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(fmt)

    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger


# ─────────────────────────────────────────────────────────────────────────────
# PHASE 1 – Load ClinPGx rsIDs
# ─────────────────────────────────────────────────────────────────────────────

def load_clinpgx(variants_tsv: Path, logger: logging.Logger) -> Dict[str, dict]:
    """
    Returns dict:  rsid (lowercase) → full ClinPGx row dict
    Only rows where Variant Name starts with 'rs' are indexed.
    """
    clinpgx: Dict[str, dict] = {}
    required = {
        "Variant ID", "Variant Name", "Gene IDs", "Gene Symbols",
        "Location", "Variant Annotation count", "Clinical Annotation count",
        "Level 1/2 Clinical Annotation count", "Guideline Annotation count",
        "Label Annotation count",
    }
    skipped = 0
    with open(variants_tsv, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not required.issubset(set(reader.fieldnames or [])):
            missing = required - set(reader.fieldnames or [])
            logger.warning(f"variants.tsv missing columns: {missing}")
        for i, row in enumerate(reader):
            rsid = (row.get("Variant Name") or "").strip().lower()
            if not rsid.startswith("rs"):
                skipped += 1
                continue
            clinpgx[rsid] = row
    logger.info(f"ClinPGx loaded: {len(clinpgx)} rsIDs  ({skipped} rows skipped/non-rs)")
    return clinpgx


# ─────────────────────────────────────────────────────────────────────────────
# VEP column helpers
# ─────────────────────────────────────────────────────────────────────────────

VEP_WANTED = [
    "#Uploaded_variation", "Location", "Allele", "Gene", "Feature",
    "Feature_type", "Consequence", "cDNA_position", "CDS_position",
    "Protein_position", "Existing_variation", "IMPACT", "DISTANCE",
    "STRAND", "VARIANT_CLASS", "SYMBOL", "SYMBOL_SOURCE", "HGNC_ID",
    "BIOTYPE", "CANONICAL", "AF", "gnomADe_AF", "gnomADe_SAS_AF",
    "CLIN_SIG", "PHENO",
]

CLINPGX_WANTED = [
    "Variant ID", "Variant Name", "Gene IDs", "Gene Symbols", "Location",
    "Variant Annotation count", "Clinical Annotation count",
    "Level 1/2 Clinical Annotation count", "Guideline Annotation count",
    "Label Annotation count",
]

OUTPUT_HEADER = VEP_WANTED + [f"ClinPGx_{c}" for c in CLINPGX_WANTED]


def parse_vep_header(line: str) -> Optional[List[str]]:
    """Return column list if this is the VEP header line, else None."""
    stripped = line.rstrip("\n")
    if stripped.startswith("#Uploaded_variation"):
        return stripped.split("\t")
    return None


def extract_rsids_from_existing(existing_variation: str) -> Set[str]:
    """Split 'rs123,rs456,COSV...' and return lowercase rsIDs."""
    return {
        tok.strip().lower()
        for tok in existing_variation.split(",")
        if tok.strip().lower().startswith("rs")
    }


# ─────────────────────────────────────────────────────────────────────────────
# PHASE 2 – Chromosome worker
# ─────────────────────────────────────────────────────────────────────────────

def process_chromosome(args: Tuple) -> dict:
    """
    Worker function (runs in a subprocess).
    Returns a summary dict for this chromosome.
    """
    (
        chrom,
        vep_path,
        clinpgx,            # dict rsid → ClinPGx row
        matched_tsv,        # Path – shared output file
        lock_path,          # Path – lockfile for append coordination
        novel_vcf_out,      # Path – novel VCF output (gz)
        vcf_input,          # Path – original VCF input (gz) for novel filter
        resume,             # bool
        log_path,           # Path
    ) = args

    logger = setup_logging(log_path)
    logger.info(f"[{chrom}] Worker started")

    stats = {
        "chrom": chrom,
        "vep_lines": 0,
        "matched_rows": 0,
        "novel_variants": 0,
        "errors": 0,
        "matched_rsids": set(),   # collected locally, merged by coordinator
        "novel_uploaded": set(),  # #Uploaded_variation values not in ClinPGx
    }

    # ── Stream VEP file ──────────────────────────────────────────────────────
    col_index: Dict[str, int] = {}
    matched_buffer: List[List[str]] = []
    FLUSH_EVERY = 5000

    open_fn = gzip.open if str(vep_path).endswith(".gz") else open

    try:
        with open_fn(vep_path, "rt", encoding="utf-8", errors="replace") as fh:
            for raw_line in fh:
                line = raw_line.rstrip("\n")

                # Skip VEP meta-comment lines
                if line.startswith("##"):
                    continue

                # Parse header
                if line.startswith("#Uploaded_variation"):
                    cols = line.split("\t")
                    col_index = {c: i for i, c in enumerate(cols)}
                    continue

                if not col_index:
                    continue  # header not yet seen

                stats["vep_lines"] += 1

                try:
                    fields = line.split("\t")

                    def get(col: str) -> str:
                        i = col_index.get(col)
                        return fields[i] if i is not None and i < len(fields) else ""

                    existing = get("Existing_variation")
                    uploaded = get("#Uploaded_variation")

                    if not existing or existing == "-":
                        stats["novel_uploaded"].add(uploaded)
                        stats["novel_variants"] += 1
                        continue

                    rsids_in_row = extract_rsids_from_existing(existing)
                    matched_rsid = rsids_in_row & clinpgx.keys()

                    if not matched_rsid:
                        stats["novel_uploaded"].add(uploaded)
                        stats["novel_variants"] += 1
                        continue

                    # Matched – build output row (one per matched rsID)
                    for rsid in matched_rsid:
                        stats["matched_rsids"].add(rsid)
                        pgx_row = clinpgx[rsid]

                        out_row = [get(c) for c in VEP_WANTED]
                        out_row += [pgx_row.get(c, "") for c in CLINPGX_WANTED]
                        matched_buffer.append(out_row)
                        stats["matched_rows"] += 1

                    if len(matched_buffer) >= FLUSH_EVERY:
                        _flush_matched(matched_buffer, matched_tsv, lock_path, logger)
                        matched_buffer.clear()

                except Exception as e:
                    stats["errors"] += 1
                    logger.debug(f"[{chrom}] Malformed line skipped: {e}")

        # Flush remainder
        if matched_buffer:
            _flush_matched(matched_buffer, matched_tsv, lock_path, logger)
            matched_buffer.clear()

    except Exception as e:
        logger.error(f"[{chrom}] Fatal VEP read error: {e}")
        stats["errors"] += 1

    logger.info(
        f"[{chrom}] VEP done: {stats['vep_lines']} lines, "
        f"{stats['matched_rows']} matched rows, "
        f"{stats['novel_variants']} novel variants"
    )

    # ── Write novel VCF ──────────────────────────────────────────────────────
    if vcf_input.exists():
        _write_novel_vcf(
            chrom, vcf_input, novel_vcf_out,
            stats["novel_uploaded"], resume, logger
        )
    else:
        logger.warning(f"[{chrom}] VCF input not found: {vcf_input}")

    # Convert set → list for pickling back to coordinator
    stats["matched_rsids"] = list(stats["matched_rsids"])
    stats["novel_uploaded"] = []  # large set; already used, don't return
    return stats


def _flush_matched(
    rows: List[List[str]],
    out_path: Path,
    lock_path: Path,
    logger: logging.Logger,
) -> None:
    """Append rows to shared matched TSV with file-lock coordination."""
    lock_path.touch(exist_ok=True)
    retries = 0
    while retries < 60:
        try:
            with open(lock_path, "r+") as lf:
                fcntl.flock(lf, fcntl.LOCK_EX)
                try:
                    with open(out_path, "a", newline="", encoding="utf-8") as fh:
                        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
                        writer.writerows(rows)
                finally:
                    fcntl.flock(lf, fcntl.LOCK_UN)
            return
        except BlockingIOError:
            retries += 1
            time.sleep(0.5)
    logger.error("Could not acquire lock after 30s – flushing anyway (risk of interleave)")
    with open(out_path, "a", newline="", encoding="utf-8") as fh:
        csv.writer(fh, delimiter="\t", lineterminator="\n").writerows(rows)


def _write_novel_vcf(
    chrom: str,
    vcf_input: Path,
    novel_vcf_out: Path,
    novel_uploaded: Set[str],
    resume: bool,
    logger: logging.Logger,
) -> None:
    """
    Stream original VCF, keep only variants whose #Uploaded_variation key
    (format: CHR_POS_REF/ALT) is in novel_uploaded set.
    Output is gzipped VCF.
    """
    if resume and novel_vcf_out.exists() and novel_vcf_out.stat().st_size > 100:
        logger.info(f"[{chrom}] Novel VCF already exists, skipping (--resume)")
        return

    novel_vcf_out.parent.mkdir(parents=True, exist_ok=True)

    # Build a fast lookup: (chrom_no_chr, pos) → True
    # VEP Uploaded_variation format: "10_10651_A/C"  (chrom without 'chr')
    pos_set: Set[Tuple[str, str]] = set()
    for uv in novel_uploaded:
        parts = uv.split("_")
        if len(parts) >= 2:
            pos_set.add((parts[0], parts[1]))

    written = 0
    open_fn = gzip.open if str(vcf_input).endswith(".gz") else open

    try:
        with open_fn(vcf_input, "rt", encoding="utf-8", errors="replace") as fin, \
             gzip.open(novel_vcf_out, "wt", encoding="utf-8") as fout:

            for raw in fin:
                line = raw  # keep original line endings for VCF
                stripped = line.rstrip("\n")

                # Header lines pass through unchanged
                if stripped.startswith("#"):
                    fout.write(line)
                    continue

                parts = stripped.split("\t")
                if len(parts) < 5:
                    continue

                chrom_field = parts[0].lstrip("chr")  # normalise 'chr1' → '1'
                pos_field   = parts[1]

                if (chrom_field, pos_field) in pos_set:
                    fout.write(line)
                    written += 1

        logger.info(f"[{chrom}] Novel VCF written: {written} variants → {novel_vcf_out.name}")

    except Exception as e:
        logger.error(f"[{chrom}] Novel VCF write error: {e}")


# ─────────────────────────────────────────────────────────────────────────────
# PHASE 3 – ClinPGx_NotFound
# ─────────────────────────────────────────────────────────────────────────────

def write_not_found(
    clinpgx: Dict[str, dict],
    all_matched_rsids: Set[str],
    out_path: Path,
    logger: logging.Logger,
) -> int:
    not_found_cols = [
        "Variant ID", "Variant Name", "Gene Symbols", "Location",
        "Clinical Annotation count", "Level 1/2 Clinical Annotation count",
        "Guideline Annotation count", "Label Annotation count",
    ]
    unmatched = {rsid: row for rsid, row in clinpgx.items() if rsid not in all_matched_rsids}
    with open(out_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=not_found_cols, delimiter="\t",
            extrasaction="ignore", lineterminator="\n"
        )
        writer.writeheader()
        for row in unmatched.values():
            writer.writerow(row)
    logger.info(f"ClinPGx_NotFound: {len(unmatched)} variants → {out_path}")
    return len(unmatched)


# ─────────────────────────────────────────────────────────────────────────────
# PHASE 5 – Summary
# ─────────────────────────────────────────────────────────────────────────────

def write_summary(
    out_path: Path,
    total_clinpgx: int,
    matched_rsids: Set[str],
    total_novel: int,
    chrom_stats: List[dict],
    logger: logging.Logger,
) -> None:
    matched = len(matched_rsids)
    unmatched = total_clinpgx - matched
    pct = (matched / total_clinpgx * 100) if total_clinpgx else 0.0

    lines = [
        "=" * 60,
        "Genome India × ClinPGx Pipeline Summary",
        "=" * 60,
        f"Total ClinPGx rsIDs             : {total_clinpgx}",
        f"Matched rsIDs                   : {matched}",
        f"Unmatched ClinPGx rsIDs         : {unmatched}",
        f"Novel Genome India variants     : {total_novel}",
        f"Match percentage                : {pct:.2f}%",
        "",
        "Per-chromosome stats:",
        f"  {'Chrom':<8} {'VEP lines':>12} {'Matched rows':>14} "
        f"{'Novel variants':>16} {'Errors':>8}",
    ]
    for s in sorted(chrom_stats, key=lambda x: x["chrom"]):
        lines.append(
            f"  {s['chrom']:<8} {s['vep_lines']:>12,} {s['matched_rows']:>14,} "
            f"{s['novel_variants']:>16,} {s['errors']:>8}"
        )

    text = "\n".join(lines) + "\n"
    out_path.write_text(text, encoding="utf-8")
    logger.info(f"Summary written → {out_path}")
    print(text)


# ─────────────────────────────────────────────────────────────────────────────
# MAIN COORDINATOR
# ─────────────────────────────────────────────────────────────────────────────

CHROMOSOMES = [str(i) for i in range(1, 23)]


def build_vep_path(vep_dir: Path, chrom: str) -> Path:
    """Try chr{N}.vep.tsv and chr{N}.vep.tsv.gz."""
    for suffix in [f"chr{chrom}.vep.tsv", f"chr{chrom}.vep.tsv.gz"]:
        p = vep_dir / suffix
        if p.exists():
            return p
    return vep_dir / f"chr{chrom}.vep.tsv"  # will be missing; worker logs it


def build_vcf_path(vcf_dir: Path, chrom: str) -> Path:
    for suffix in [f"chr{chrom}.vcf.gz", f"chr{chrom}.vcf"]:
        p = vcf_dir / suffix
        if p.exists():
            return p
    return vcf_dir / f"chr{chrom}.vcf.gz"


def main():
    parser = argparse.ArgumentParser(description="Genome India × ClinPGx Pipeline")
    parser.add_argument("--variants", required=True, type=Path, help="ClinPGx variants.tsv")
    parser.add_argument("--vcf-dir",  required=True, type=Path, help="Directory with chr*.vcf.gz")
    parser.add_argument("--vep-dir",  required=True, type=Path, help="Directory with chr*.vep.tsv")
    parser.add_argument("--output",   required=True, type=Path, help="Output directory")
    parser.add_argument("--workers",  type=int, default=12)
    parser.add_argument("--resume",   action="store_true",
                        help="Skip chromosomes whose novel VCF already exists")
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    novel_vcf_dir = args.output / "NovelVCF"
    novel_vcf_dir.mkdir(exist_ok=True)

    log_path      = args.output / "pipeline.log"
    matched_tsv   = args.output / "GenomeIndia_ClinPGx_Matched.tsv"
    not_found_tsv = args.output / "ClinPGx_NotFound.tsv"
    summary_txt   = args.output / "Summary.txt"
    lock_path     = args.output / ".matched.lock"

    logger = setup_logging(log_path)
    logger.info("=" * 60)
    logger.info("Pipeline started")
    logger.info(f"  variants : {args.variants}")
    logger.info(f"  vcf-dir  : {args.vcf_dir}")
    logger.info(f"  vep-dir  : {args.vep_dir}")
    logger.info(f"  output   : {args.output}")
    logger.info(f"  workers  : {args.workers}")
    logger.info(f"  resume   : {args.resume}")

    # ── Phase 1: Load ClinPGx ─────────────────────────────────────────────
    clinpgx = load_clinpgx(args.variants, logger)
    total_clinpgx = len(clinpgx)

    # ── Initialise matched TSV (write header once) ────────────────────────
    if not matched_tsv.exists() or matched_tsv.stat().st_size == 0:
        with open(matched_tsv, "w", newline="", encoding="utf-8") as fh:
            csv.writer(fh, delimiter="\t", lineterminator="\n").writerow(OUTPUT_HEADER)

    # ── Phase 2: Dispatch chromosome workers ─────────────────────────────
    worker_args = []
    for chrom in CHROMOSOMES:
        vep_path     = build_vep_path(args.vep_dir, chrom)
        vcf_input    = build_vcf_path(args.vcf_dir, chrom)
        novel_vcf_out = novel_vcf_dir / f"chr{chrom}_novel.vcf.gz"

        if not vep_path.exists():
            logger.warning(f"VEP file not found, skipping chromosome {chrom}: {vep_path}")
            continue

        worker_args.append((
            chrom,
            vep_path,
            clinpgx,
            matched_tsv,
            lock_path,
            novel_vcf_out,
            vcf_input,
            args.resume,
            log_path,
        ))

    logger.info(f"Dispatching {len(worker_args)} chromosome workers with {args.workers} processes")

    with multiprocessing.Pool(processes=args.workers) as pool:
        results: List[dict] = pool.map(process_chromosome, worker_args)

    # ── Aggregate results ─────────────────────────────────────────────────
    all_matched_rsids: Set[str] = set()
    total_novel = 0
    for res in results:
        all_matched_rsids.update(res.get("matched_rsids", []))
        total_novel += res.get("novel_variants", 0)

    logger.info(f"All workers done. Unique matched rsIDs: {len(all_matched_rsids)}")

    # ── Phase 3: ClinPGx_NotFound ─────────────────────────────────────────
    write_not_found(clinpgx, all_matched_rsids, not_found_tsv, logger)

    # ── Phase 5: Summary ──────────────────────────────────────────────────
    write_summary(
        summary_txt,
        total_clinpgx,
        all_matched_rsids,
        total_novel,
        results,
        logger,
    )

    logger.info("Pipeline complete.")


if __name__ == "__main__":
    # Ensure child processes can import this module cleanly on all platforms
    multiprocessing.set_start_method("fork", force=True)
    main()
