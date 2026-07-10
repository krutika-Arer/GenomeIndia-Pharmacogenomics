import pandas as pd
import gzip
import os
import re
from glob import glob

ROOT = "/mnt/Linux_storage/pharmacogenomics_9768GI_SummaryStats"

CLINPGX = f"{ROOT}/variants.tsv"
VEP_DIR = f"{ROOT}/vep_pgx"

OUTDIR = f"{ROOT}/results"

os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(f"{OUTDIR}/NovelVEP", exist_ok=True)

##########################################################
# LOAD CLINPGX
##########################################################

print("Loading ClinPGx...")

clin = pd.read_csv(
    CLINPGX,
    sep="\t",
    dtype=str,
    low_memory=False
)

clin["Variant Name"] = clin["Variant Name"].astype(str)

clin_dict = {
    r["Variant Name"]: r.to_dict()
    for _, r in clin.iterrows()
}

clin_rsids = set(clin_dict.keys())

print("ClinPGx rsIDs:", len(clin_rsids))

##########################################################
# OUTPUT FILES
##########################################################

matched_out = gzip.open(
    f"{OUTDIR}/GenomeIndia_ClinPGx_Matched.tsv.gz",
    "wt"
)

novel_out = gzip.open(
    f"{OUTDIR}/IndianNovel_Annotated.tsv.gz",
    "wt"
)

##########################################################
# VEP COLUMNS
##########################################################

wanted = [
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

matched_header_written = False
novel_header_written = False

matched_rsids = set()

##########################################################
# PROCESS VEP FILES
##########################################################

vep_files = sorted(glob(f"{VEP_DIR}/chr*.vep.tsv"))

for vep_file in vep_files:

    chrom = os.path.basename(
        vep_file
    ).split(".")[0]

    print(f"\nProcessing {chrom}")

    novel_chr = gzip.open(
        f"{OUTDIR}/NovelVEP/{chrom}.novel.tsv.gz",
        "wt"
    )

    with open(vep_file) as f:

        header = None

        for line in f:

            if line.startswith("##"):
                continue

            if line.startswith("#Uploaded_variation"):

                header = line.rstrip().split("\t")

                idx = {
                    c:i for i,c in enumerate(header)
                }

                keep_idx = [
                    idx[x]
                    for x in wanted
                    if x in idx
                ]

                out_header = [
                    header[i]
                    for i in keep_idx
                ]

                if not matched_header_written:

                    pgx_cols = [
                        "Variant ID",
                        "Variant Name",
                        "Gene IDs",
                        "Gene Symbols",
                        "Location",
                        "Variant Annotation count",
                        "Clinical Annotation count",
                        "Level 1/2 Clinical Annotation count",
                        "Guideline Annotation count",
                        "Label Annotation count"
                    ]

                    matched_out.write(
                        "\t".join(
                            out_header + pgx_cols
                        ) + "\n"
                    )

                    matched_header_written = True

                if not novel_header_written:

                    novel_out.write(
                        "\t".join(out_header) + "\n"
                    )

                    novel_chr.write(
                        "\t".join(out_header) + "\n"
                    )

                    novel_header_written = True

                continue

            vals = line.rstrip("\n").split("\t")

            if len(vals) < len(header):
                continue

            rsid = vals[
                idx["Existing_variation"]
            ]

            selected = [
                vals[i]
                for i in keep_idx
            ]

            ################################################
            # MATCHED
            ################################################

            if rsid in clin_rsids:

                matched_rsids.add(rsid)

                pgx = clin_dict[rsid]

                pgx_values = [

    str(pgx.get("Variant ID","")),
    str(pgx.get("Variant Name","")),
    str(pgx.get("Gene IDs","")),
    str(pgx.get("Gene Symbols","")),
    str(pgx.get("Location","")),
    str(pgx.get("Variant Annotation count","")),
    str(pgx.get("Clinical Annotation count","")),
    str(pgx.get("Level 1/2 Clinical Annotation count","")),
    str(pgx.get("Guideline Annotation count","")),
    str(pgx.get("Label Annotation count",""))

]

                matched_out.write(
                    "\t".join(selected + pgx_values)
                    + "\n"
                )

            ################################################
            # NOVEL
            ################################################

            else:

                novel_out.write(
                    "\t".join(selected) + "\n"
                )

                novel_chr.write(
                    "\t".join(selected) + "\n"
                )

    novel_chr.close()

##########################################################
# CLINPGX NOT FOUND
##########################################################

print("\nWriting ClinPGx_NotFound")

not_found = clin[
    ~clin["Variant Name"].isin(
        matched_rsids
    )
]

not_found.to_csv(
    f"{OUTDIR}/ClinPGx_NotFound.tsv",
    sep="\t",
    index=False
)

matched_out.close()
novel_out.close()

print("\nDONE")
