import pandas as pd
import gzip
import glob
import os
from collections import defaultdict

ROOT = "/mnt/Linux_storage/pharmacogenomics_9768GI_SummaryStats"

CLINPGX = f"{ROOT}/variants.tsv"
VEP_DIR = f"{ROOT}/vep_pgx"
VCF_DIR = f"{ROOT}/vcf"

OUTDIR = f"{ROOT}/results"

os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(f"{OUTDIR}/NovelVCF", exist_ok=True)

##############################################################################
# STEP 1
# Load ClinPGx
##############################################################################

print("Loading ClinPGx...")

clin = pd.read_csv(
    CLINPGX,
    sep="\t",
    dtype=str,
    low_memory=False
)

clin["Variant Name"] = clin["Variant Name"].astype(str)

clin_rsids = set(
    clin["Variant Name"]
    .dropna()
    .unique()
)

print("ClinPGx rsIDs:", len(clin_rsids))

##############################################################################
# STEP 2
# VEP columns to keep
##############################################################################

wanted_cols = [

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
"IMPACT",
"DISTANCE",
"STRAND",
"VARIANT_CLASS",
"SYMBOL",
"SYMBOL_SOURCE",
"HGNC_ID",
"BIOTYPE",
"CANONICAL",
"AF",
"gnomADe_AF",
"gnomADe_SAS_AF",
"CLIN_SIG",
"PHENO"

]

##############################################################################
# STEP 3
# Match ClinPGx rsIDs against VEP
##############################################################################

matched_rsids = set()

matched_file = f"{OUTDIR}/GenomeIndia_ClinPGx_Matched.tsv.gz"

first_write = True

vep_files = sorted(
    glob.glob(f"{VEP_DIR}/chr*.vep.tsv")
)

for vep in vep_files:

    print("\nScanning", os.path.basename(vep))

    chunks = pd.read_csv(
        vep,
        sep="\t",
        comment="#",
        chunksize=100000,
        dtype=str,
        low_memory=False
    )

    for chunk in chunks:

        if "Existing_variation" not in chunk.columns:
            continue

        chunk["Existing_variation"] = (
            chunk["Existing_variation"]
            .fillna("")
            .astype(str)
        )

        mask = chunk["Existing_variation"].isin(clin_rsids)

        hits = chunk.loc[mask]

        if len(hits) == 0:
            continue

        matched_rsids.update(
            hits["Existing_variation"].unique()
        )

        keep = [
            c for c in wanted_cols
            if c in hits.columns
        ]

        hits = hits[keep]

        merged = hits.merge(
            clin,
            left_on="Existing_variation",
            right_on="Variant Name",
            how="left"
        )

        merged.to_csv(
            matched_file,
            sep="\t",
            index=False,
            mode="wt" if first_write else "at",
            compression="gzip",
            header=first_write
        )

        first_write = False

##############################################################################
# STEP 4
# ClinPGx not found
##############################################################################

print("\nCreating ClinPGx_NotFound")

not_found = clin[
    ~clin["Variant Name"].isin(matched_rsids)
]

not_found.to_csv(
    f"{OUTDIR}/ClinPGx_NotFound.tsv",
    sep="\t",
    index=False
)

##############################################################################
# STEP 5
# Novel Indian annotated variants
##############################################################################

novel_file = f"{OUTDIR}/IndianNovel_Annotated.tsv.gz"

first_write = True

for vep in vep_files:

    print("\nNovel scan", os.path.basename(vep))

    chunks = pd.read_csv(
        vep,
        sep="\t",
        comment="#",
        chunksize=100000,
        dtype=str,
        low_memory=False
    )

    for chunk in chunks:

        if "Existing_variation" not in chunk.columns:
            continue

        chunk["Existing_variation"] = (
            chunk["Existing_variation"]
            .fillna("")
            .astype(str)
        )

        novel = chunk[
            ~chunk["Existing_variation"].isin(clin_rsids)
        ]

        if len(novel) == 0:
            continue

        keep = [
            c for c in wanted_cols
            if c in novel.columns
        ]

        novel = novel[keep]

        novel.to_csv(
            novel_file,
            sep="\t",
            index=False,
            compression="gzip",
            mode="wt" if first_write else "at",
            header=first_write
        )

        first_write = False

##############################################################################
# STEP 6
# Create chromosome-wise Novel VCF
##############################################################################

print("\nCreating Novel VCFs")

for chrom in range(1,23):

    chrname = f"chr{chrom}"

    vcf = f"{VCF_DIR}/{chrname}.vcf.gz"

    if not os.path.exists(vcf):
        continue

    outfile = f"{OUTDIR}/NovelVCF/{chrname}.novel.vcf"

    with gzip.open(vcf, "rt") as fin, \
         open(outfile, "w") as fout:

        for line in fin:

            if line.startswith("#"):
                fout.write(line)
                continue

            cols = line.rstrip().split("\t")

            variant_id = cols[2]

            if variant_id in clin_rsids:
                continue

            fout.write(line)

print("\nDONE")
