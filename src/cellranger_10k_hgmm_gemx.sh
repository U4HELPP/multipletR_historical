#!/usr/bin/env bash

url_5p=https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-vdj/8.0.0/10k_hgmm_5p_gemx_Multiplex/10k_hgmm_5p_gemx_Multiplex_fastqs.tar
url_3p=https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/8.0.0/10k_hgmm_3p_gemx_Multiplex/10k_hgmm_3p_gemx_Multiplex_fastqs.tar

wget $url_5p
wget $url_3p

tar -xf "10k_hgmm_5p_gemx_Multiplex_fastqs.tar"
tar -xf "10k_hgmm_3p_gemx_Multiplex_fastqs.tar"

module load cellranger/8.0.1

# 5p
cellranger count \
    --id 10k_hgmm_5p_gemx_gex1 \
    --transcriptome /lustre/home/juicer/ExtData/10x/refdata-gex-GRCh38_and_GRCm39-2024-A \
    --fastqs /lustre/home/wallbp/10x_MultipletR/10k_hgmm_5p_gemx_fastqs \
    --sample 10k_hgmm_5p_gemx_gex1 \
    --create-bam false
cellranger count \
    --id 10k_hgmm_5p_gemx_gex2 \
    --transcriptome /lustre/home/juicer/ExtData/10x/refdata-gex-GRCh38_and_GRCm39-2024-A \
    --fastqs /lustre/home/wallbp/10x_MultipletR/10k_hgmm_5p_gemx_fastqs \
    --sample 10k_hgmm_5p_gemx_gex2 \
    --create-bam false

# 3p
cellranger count \
    --id 10k_hgmm_3p_gemx_gex1 \
    --transcriptome /lustre/home/juicer/ExtData/10x/refdata-gex-GRCh38_and_GRCm39-2024-A \
    --fastqs /lustre/home/wallbp/10x_MultipletR/10k_hgmm_3p_gemx_fastqs \
    --sample 10k_hgmm_3p_gemx_gex1 \
    --create-bam false
cellranger count \
    --id 10k_hgmm_3p_gemx_gex2 \
    --transcriptome /lustre/home/juicer/ExtData/10x/refdata-gex-GRCh38_and_GRCm39-2024-A \
    --fastqs /lustre/home/wallbp/10x_MultipletR/10k_hgmm_3p_gemx_fastqs \
    --sample 10k_hgmm_3p_gemx_gex2 \
    --create-bam false
