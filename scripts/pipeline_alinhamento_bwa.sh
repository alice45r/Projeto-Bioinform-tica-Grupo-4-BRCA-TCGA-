#!/bin/bash
set -e

mkdir -p ~/dataset/mapeamento_hg38
cd ~/dataset/mapeamento_hg38

REF=~/dataset/ref_hg38/hg38.fa
FASTQ_DIR=~/dataset/sra_exoma

for SRR in SRR8619186 SRR8619133 SRR8619113 SRR8619154; do
    echo "=== A mapear amostra $SRR com BWA-MEM (8 nucleos) em $(date) ==="
    R1="${FASTQ_DIR}/${SRR}_1.fastq.gz"
    R2="${FASTQ_DIR}/${SRR}_2.fastq.gz"
    RG="@RG\tID:${SRR}\tSM:${SRR}\tPL:ILLUMINA"

    bwa mem -t 8 -R "$RG" "$REF" "$R1" "$R2" | samtools view -Sb - | samtools sort -@ 4 - "${SRR}.sorted"
    samtools index "${SRR}.sorted.bam"
    samtools flagstat "${SRR}.sorted.bam" > "${SRR}_flagstat.txt"
    echo "=== Amostra $SRR concluida em $(date) ==="
done

echo "=== Mapeamento concluido com sucesso em $(date) ==="
