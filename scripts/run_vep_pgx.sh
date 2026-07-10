#!/usr/bin/env bash
set -euo pipefail

THREADS=12

VCF_DIR="/mnt/Linux_storage/pharmacogenomics_9768GI_SummaryStats/vcf_nochr"

CACHE="/mnt/Linux_storage/vep_resources/cache"

FASTA="/mnt/Linux_storage/vep_resources/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

OUTDIR="/mnt/Linux_storage/pharmacogenomics_9768GI_SummaryStats/vep_pgx"

mkdir -p "${OUTDIR}"

for VCF in ${VCF_DIR}/chr*.vcf.gz
do

    BASE=$(basename "${VCF}" .vcf.gz)

    echo "Processing ${BASE}"

    vep \
        -i "${VCF}" \
        -o "${OUTDIR}/${BASE}.vep.tsv" \
        --offline \
        --cache \
        --assembly GRCh38 \
        --dir_cache "${CACHE}" \
        --fasta "${FASTA}" \
        --fork ${THREADS} \
        --tab \
        --symbol \
        --canonical \
        --mane \
        --biotype \
        --variant_class \
        --hgvs \
        --protein \
        --sift b \
        --polyphen b \
        --check_existing \
        --af \
        --af_gnomad \
        --pubmed \
        --force_overwrite

done
