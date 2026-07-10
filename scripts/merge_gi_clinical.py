#!/usr/bin/env python3
"""
merge_gi_clinical.py
====================
Step 1: Add Genome India allele frequency (GI_AF) to GenomeIndia_ClinPGx_Matched.tsv
        Join key: CHROM + POS  (parsed from the Location column)
Step 2: Merge clinicalVariants.tsv on rsID
        Join key: rsID parsed from the Existing_variation column

Output: Master_PGx_Final.tsv

Usage:
    python merge_gi_clinical.py \\
        --matched    output/GenomeIndia_ClinPGx_Matched.tsv \\
        --vcf-dir    vcf/ \\
        --clinical   clinicalVariants.tsv \\
        --out        output/Master_PGx_Final.tsv

Notes:
  - Streams Matched.tsv line-by-line (no pandas, no full-file RAM load)
  - VCF AF index is built per-chromosome on demand (one chrom at a time)
  - clinicalVariants.tsv is small; loaded fully into a dict
  - Rows where no VCF match is found get GI_AF = "."
  - Rows where no clinical match is found get clinical columns = "."
  - Malformed lines are skipped and logged
"""

import re
import csv
import sys
import gzip
import logging
import argparse
from pathlib import Path
from typing import Dict, Optional, Tuple

# ──────────────────────────────────────────────────────────────────
# Logging
# ──────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler("merge_gi_clinical.log", mode="a"),
    ],
)
log = logging.getLogger(__name__)


# ──────────────────────────────────────────────────────────────────
# Parse VEP Location → (chrom, pos)
# Examples:  "6:1388744-1388747"  →  ("6", "1388744")
#            "10:10651"           →  ("10", "10651")
# ──────────────────────────────────────────────────────────────────
_LOC_RE = re.compile(r"^(\w+):(\d+)")

def parse_location(loc: str) -> Optional[Tuple[str, str]]:
    m = _LOC_RE.match(loc.strip())
    if not m:
        return None
    chrom = m.group(1).lstrip("chr")   # normalise → no 'chr' prefix
    pos   = m.group(2)
    return chrom, pos


# ──────────────────────────────────────────────────────────────────
# Build AF index for one chromosome from VCF
# Returns dict:  "chrom\tpos" → AF string
# ──────────────────────────────────────────────────────────────────
def build_vcf_af_index(vcf_path: Path) -> Dict[str, str]:
    index: Dict[str, str] = {}
    open_fn = gzip.open if str(vcf_path).endswith(".gz") else open
    try:
        with open_fn(vcf_path, "rt", encoding="utf-8", errors="replace") as fh:
            for raw in fh:
                line = raw.rstrip("\n")
                if line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 8:
                    continue
                chrom = parts[0].lstrip("chr")
                pos   = parts[1]
                info  = parts[7]
                # Extract AF=... from INFO field
                af = "."
                for token in info.split(";"):
                    if token.startswith("AF="):
                        af = token[3:]
                        break
                key = f"{chrom}\t{pos}"
                index[key] = af
    except Exception as e:
        log.warning(f"Error reading VCF {vcf_path}: {e}")
    log.info(f"  VCF index built: {len(index):,} positions from {vcf_path.name}")
    return index


def find_vcf(vcf_dir: Path, chrom: str) -> Optional[Path]:
    for name in [f"chr{chrom}.vcf.gz", f"chr{chrom}.vcf"]:
        p = vcf_dir / name
        if p.exists():
            return p
    return None


# ──────────────────────────────────────────────────────────────────
# Load clinicalVariants.tsv
# Returns dict:  rsid (lowercase) → row dict
# Only indexes rows whose 'variant' column starts with 'rs'
# ──────────────────────────────────────────────────────────────────
CLINICAL_COLS = ["variant", "gene", "type", "level of evidence", "chemicals", "phenotypes"]

def load_clinical(clinical_tsv: Path) -> Dict[str, dict]:
    index: Dict[str, dict] = {}
    skipped = 0
    with open(clinical_tsv, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            v = (row.get("variant") or "").strip()
            if v.lower().startswith("rs"):
                index[v.lower()] = row
            else:
                skipped += 1
    log.info(f"clinicalVariants: {len(index)} rsID rows indexed, {skipped} non-rsID rows skipped")
    return index


# ──────────────────────────────────────────────────────────────────
# Extract rsIDs from Existing_variation field
# e.g. "rs199564443" or "rs123,rs456,COSV..."
# ──────────────────────────────────────────────────────────────────
def extract_rsids(existing: str):
    return [
        tok.strip().lower()
        for tok in existing.split(",")
        if tok.strip().lower().startswith("rs")
    ]


# ──────────────────────────────────────────────────────────────────
# Main merge
# ──────────────────────────────────────────────────────────────────
def merge(matched_tsv: Path, vcf_dir: Path, clinical_tsv: Path, out_path: Path):

    # -- Load clinical lookup (small file, fine to hold in RAM) ----
    clinical = load_clinical(clinical_tsv)

    # -- Prefix for clinical columns in output ---------------------
    CLIN_PREFIX = "ClinVar_"
    clinical_out_cols = [f"{CLIN_PREFIX}{c}" for c in CLINICAL_COLS]

    # -- Stream Matched.tsv ----------------------------------------
    vcf_cache: Dict[str, Dict[str, str]] = {}   # chrom → AF index
    current_chrom: Optional[str] = None
    af_index: Dict[str, str] = {}

    total = matched = gi_af_found = clinical_found = errors = 0

    with open(matched_tsv, newline="", encoding="utf-8") as fin, \
         open(out_path, "w", newline="", encoding="utf-8") as fout:

        reader = csv.DictReader(fin, delimiter="\t")
        if reader.fieldnames is None:
            log.error("Matched.tsv appears empty or has no header")
            return

        # Output header = all input cols + GI_AF + clinical cols
        out_fields = list(reader.fieldnames) + ["GI_AF"] + clinical_out_cols
        writer = csv.DictWriter(fout, fieldnames=out_fields, delimiter="\t",
                                 extrasaction="ignore", lineterminator="\n")
        writer.writeheader()

        for row in reader:
            total += 1

            # ── Step 1: GI allele frequency ──────────────────────
            gi_af = "."
            loc = (row.get("Location") or "").strip()
            parsed = parse_location(loc)

            if parsed:
                chrom, pos = parsed

                # Load VCF index if we're on a new chromosome
                if chrom != current_chrom:
                    current_chrom = chrom
                    if chrom not in vcf_cache:
                        vcf_path = find_vcf(vcf_dir, chrom)
                        if vcf_path:
                            log.info(f"Loading VCF for chr{chrom} ...")
                            vcf_cache[chrom] = build_vcf_af_index(vcf_path)
                        else:
                            log.warning(f"No VCF found for chr{chrom}")
                            vcf_cache[chrom] = {}
                    af_index = vcf_cache[chrom]

                key = f"{chrom}\t{pos}"
                if key in af_index:
                    gi_af = af_index[key]
                    gi_af_found += 1
            else:
                errors += 1
                log.debug(f"Could not parse Location: {loc!r}")

            row["GI_AF"] = gi_af

            # ── Step 2: Clinical metadata by rsID ─────────────────
            existing = (row.get("Existing_variation") or "").strip()
            rsids = extract_rsids(existing)

            clin_row: Optional[dict] = None
            for rsid in rsids:
                if rsid in clinical:
                    clin_row = clinical[rsid]
                    clinical_found += 1
                    break   # take first match

            if clin_row:
                matched += 1
                for col in CLINICAL_COLS:
                    row[f"{CLIN_PREFIX}{col}"] = clin_row.get(col, ".")
            else:
                for col in clinical_out_cols:
                    row[col] = "."

            writer.writerow(row)

            if total % 500_000 == 0:
                log.info(f"  Processed {total:,} rows …")

    log.info("=" * 55)
    log.info(f"Total rows processed    : {total:,}")
    log.info(f"GI_AF matched           : {gi_af_found:,}")
    log.info(f"Clinical rsID matched   : {clinical_found:,}")
    log.info(f"Location parse errors   : {errors:,}")
    log.info(f"Output written to       : {out_path}")
    log.info("=" * 55)


# ──────────────────────────────────────────────────────────────────
# Entry point
# ──────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="Add GI AF + clinical metadata to Matched.tsv")
    parser.add_argument("--matched",  required=True, type=Path,
                        help="output/GenomeIndia_ClinPGx_Matched.tsv")
    parser.add_argument("--vcf-dir",  required=True, type=Path,
                        help="Directory containing chr*.vcf.gz files")
    parser.add_argument("--clinical", required=True, type=Path,
                        help="clinicalVariants.tsv")
    parser.add_argument("--out",      required=True, type=Path,
                        help="Output: Master_PGx_Final.tsv")
    args = parser.parse_args()

    for p in [args.matched, args.vcf_dir, args.clinical]:
        if not p.exists():
            log.error(f"Not found: {p}")
            sys.exit(1)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    merge(args.matched, args.vcf_dir, args.clinical, args.out)


if __name__ == "__main__":
    main()
