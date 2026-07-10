#!/usr/bin/env python3

import os
import gzip
import re
import pandas as pd
from concurrent.futures import ThreadPoolExecutor

###########################################################################
# SETTINGS
###########################################################################

THREADS = 12

BASE = "/mnt/Linux_storage/pharmacogenomics_9768GI_SummaryStats"

VCF_DIR = os.path.join(BASE, "vcf")

OUTPUT_DIR = os.path.join(BASE, "output")

FILES_TO_UPDATE = [
    "GenomeIndia_ClinPGx_Clinical_RSMatch.tsv",
    "GenomeIndia_ClinPGx_Clinical_SelectedPhenotypes.tsv"
]

###########################################################################
# Build Genome India AF lookup
###########################################################################

def process_chr(chrom):

    lookup = {}

    vcf = os.path.join(VCF_DIR, f"chr{chrom}.vcf.gz")

    print(f"Reading chr{chrom}")

    with gzip.open(vcf, "rt") as f:

        for line in f:

            if line.startswith("#"):
                continue

            fields = line.rstrip().split("\t")

            chrom = fields[0].replace("chr","")
            pos   = fields[1]
            ref   = fields[3]
            alt   = fields[4]

            uploaded = f"{chrom}_{pos}_{ref}/{alt}"

            info = fields[7]

            af = ""

            m = re.search(r"AF=([^;]+)", info)

            if m:
                af = m.group(1)

            lookup[uploaded] = af

    print(f"Finished chr{chrom}")

    return lookup

###########################################################################
# Main
###########################################################################

print("\nBuilding Genome India AF lookup...\n")

lookup = {}

with ThreadPoolExecutor(max_workers=THREADS) as exe:

    futures = [exe.submit(process_chr, c) for c in range(1,23)]

    for f in futures:

        lookup.update(f.result())

print("\nTotal variants indexed:", len(lookup))

###########################################################################
# Update files
###########################################################################

for file in FILES_TO_UPDATE:

    infile = os.path.join(OUTPUT_DIR, file)

    outfile = infile.replace(".tsv","_AF.tsv")

    print("\nProcessing", file)

    df = pd.read_csv(
        infile,
        sep="\t",
        dtype=str,
        low_memory=False
    )

    #######################################################################
    # Rename existing AF
    #######################################################################

    if "AF" in df.columns:

        df.rename(
            columns={"AF":"Ensembl_AF"},
            inplace=True
        )

    #######################################################################
    # Add GenomeIndia AF
    #######################################################################

    df["GenomeIndia_AF"] = df["#Uploaded_variation"].map(lookup)

    #######################################################################
    # Save
    #######################################################################

    df.to_csv(
        outfile,
        sep="\t",
        index=False
    )

    matched = df["GenomeIndia_AF"].notna().sum()

    print(f"Saved : {outfile}")
    print(f"Matched AF : {matched:,} / {len(df):,}")

print("\nDone.")
