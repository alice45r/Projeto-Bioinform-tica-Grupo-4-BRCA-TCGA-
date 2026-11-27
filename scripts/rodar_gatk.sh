#!/bin/bash
set -e

GATK=~/gatk/gatk
REF=~/dataset/ref_hg38/hg38.fa
BAM_DIR=~/dataset/mapeamento_hg38
OUT_DIR=~/dataset/mutacoes_somaticas

mkdir -p $OUT_DIR

# 1. Criar índice .fai e dicionário .dict do hg38.fa (exigido pelo GATK)
if [ ! -f "${REF}.fai" ]; then
    echo "Criando índice faidx..."
    samtools faidx $REF
fi

if [ ! -f "${REF%.fa}.dict" ]; then
    echo "Criando dicionário GATK..."
    $GATK CreateSequenceDictionary -R $REF -O ${REF%.fa}.dict
fi

# 2. Chamada somática (Mutect2) e filtragem (FilterMutectCalls) para as 4 linhagens
for SAMPLE in SRR8619186 SRR8619133 SRR8619113 SRR8619154; do
    echo "=========================================="
    echo "Processando amostra: $SAMPLE ($(date))"
    echo "=========================================="

    # Chamada somática em modo Tumor-Only
    $GATK --java-options "-Xmx8G" Mutect2 \
        -R $REF \
        -I ${BAM_DIR}/${SAMPLE}.sorted.bam \
        -O ${OUT_DIR}/${SAMPLE}.unfiltered.vcf.gz

    # Filtragem estatística das mutações somáticas
    $GATK --java-options "-Xmx8G" FilterMutectCalls \
        -R $REF \
        -V ${OUT_DIR}/${SAMPLE}.unfiltered.vcf.gz \
        -O ${OUT_DIR}/${SAMPLE}.filtered.vcf.gz

    # Manter apenas mutações aprovadas (PASS)
    zcat ${OUT_DIR}/${SAMPLE}.filtered.vcf.gz | grep -E "^#|PASS" > ${OUT_DIR}/${SAMPLE}.somatic.pass.vcf

    echo "Amostra $SAMPLE finalizada em $(date)"
done

echo "=== GATK CONCLUIDO COM SUCESSO! ==="
