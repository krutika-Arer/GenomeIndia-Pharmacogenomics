#!/usr/bin/env python3

import gzip
import os
from multiprocessing import Pool

ROOT = "/mnt/Linux_storage/pharmacogenomics_9768GI_SummaryStats"

NOVEL_DIR = f"{ROOT}/output/NovelVCF"
VEP_DIR   = f"{ROOT}/vep_pgx"
OUTDIR    = f"{ROOT}/output/NovelAnnotated"

os.makedirs(OUTDIR, exist_ok=True)

WANTED = [
    "Uploaded_variation",
    "Location",
    "Allele",
    "Gene",
    "Feature",
    "Feature_type",
    "Consequence",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
    "Existing_variation",
    "DISTANCE",
    "VARIANT_CLASS",
    "SYMBOL",
    "HGNC_ID",
    "BIOTYPE",
    "CANONICAL"
]

def process_chrom(chrom):

    chr_name = f"chr{chrom}"

    print(f"[START] {chr_name}")

    novel_variants = {}

    ########################################################
    # Load Novel VCF
    ########################################################

    novel_file_gz = f"{NOVEL_DIR}/{chr_name}_novel.vcf.gz"
    novel_file    = f"{NOVEL_DIR}/{chr_name}_novel.vcf"

    if os.path.exists(novel_file_gz):
        opener = gzip.open
        infile = novel_file_gz

    elif os.path.exists(novel_file):
        opener = open
        infile = novel_file

    else:
        print(f"[SKIP] {chr_name} VCF not found")
        return

    with opener(infile, "rt") as f:

        for line in f:

            if line.startswith("#"):
                continue

            cols = line.rstrip("\n").split("\t")

            if len(cols) < 8:
                continue

            pos = cols[1]

            af = ""

            for item in cols[7].split(";"):
                if item.startswith("AF="):
                    af = item[3:]
                    break

            novel_variants[pos] = (
                cols[0],   # CHROM
                cols[1],   # POS
                cols[3],   # REF
                cols[4],   # ALT
                af
            )

    print(f"[{chr_name}] loaded {len(novel_variants):,} novel variants")

    ########################################################
    # VEP
    ########################################################

    vep_file = f"{VEP_DIR}/{chr_name}.vep.tsv"

    if not os.path.exists(vep_file):
        print(f"[SKIP] {vep_file} missing")
        return

    out_file = f"{OUTDIR}/{chr_name}_novel_annotated.tsv.gz"

    with gzip.open(out_file, "wt") as out:

        with open(vep_file, "r") as f:

            idx = {}
            keep_idx = []

            for line in f:

                if line.startswith("##"):
                    continue

                ####################################################
                # Header
                ####################################################

                if line.startswith("#Uploaded_variation"):

                    header = line.rstrip("\n").split("\t")

                    idx = {
                        c:i
                        for i,c in enumerate(header)
                    }

                    keep_idx = [
                        idx[x]
                        for x in WANTED
                        if x in idx
                    ]

                    out.write(
                        "\t".join([
                            "#CHROM",
                            "POS",
                            "REF",
                            "ALT",
                            "AF"
                        ] + WANTED)
                        + "\n"
                    )

                    continue

                cols = line.rstrip("\n").split("\t")

                if not idx:
                    continue

                try:

                    loc = cols[idx["Location"]]

                    pos = loc.split(":")[1].split("-")[0]

                except Exception:
                    continue

                if pos not in novel_variants:
                    continue

                ####################################################
                # Canonical only
                ####################################################

                if "CANONICAL" in idx:

                    try:
                        if cols[idx["CANONICAL"]] != "YES":
                            continue
                    except Exception:
                        continue

                row = list(novel_variants[pos])

                row.extend(
                    cols[i]
                    if i < len(cols)
                    else ""
                    for i in keep_idx
                )

                out.write(
                    "\t".join(map(str,row))
                    + "\n"
                )

    print(f"[DONE] {chr_name}")

if __name__ == "__main__":

    WORKERS = 12

    with Pool(WORKERS) as pool:
        pool.map(process_chrom, range(1,23))

    print("\nAll chromosomes completed.")
