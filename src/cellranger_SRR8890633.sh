#!/usr/bin/env bash

module load sratoolkit

threads=64
mem=116

# Fetch data
prefetch SRR8890633 --max-size=50GB
fasterq-dump --split-files --include-technical --threads $threads SRR8890633

# Cleanup
# rm -rf "./SRR8890633"

bgzip -@ $threads "SRR8890633_1.fastq"
bgzip -@ $threads "SRR8890633_2.fastq"
bgzip -@ $threads "SRR8890633_3.fastq"

# Adjust as needed
mkdir -p "SRR8890633_FASTQs"
mv "SRR8890633_1.fastq.gz" "./SRR8890633_FASTQs/SRR8890633_S1_L1_R1_001.fastq.gz"
mv "SRR8890633_2.fastq.gz" "./SRR8890633_FASTQs/SRR8890633_S1_L1_R2_001.fastq.gz"
mv "SRR8890633_3.fastq.gz" "./SRR8890633_FASTQs/SRR8890633_S1_L1_I1_001.fastq.gz"

module load cellranger/8.0.1

cellranger count \
    --id SRR8890633 \
    --transcriptome /lustre/home/juicer/ExtData/10x/refdata-gex-GRCh38_and_GRCm39-2024-A \
    --fastqs "./SRR8890633_FASTQs" \
    --sample SRR8890633 \
    --create-bam false \
    --localcores $threads \
    --localmem $mem
