#!/bin/bash
export MAMBA_ROOT_PREFIX=~/micromamba
GATK=~/gatk/gatk
REF=~/dataset/ref_hg38/hg38.fa
BAM_DIR=~/dataset/mapeamento_hg38
OUT_DIR=~/dataset/mutacoes_somaticas
mkdir -p $OUT_DIR

CHRS="-L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22 -L chrX -L chrY"

processar_amostra() {
    SAMPLE=$1
    echo "=== Iniciando $SAMPLE em $(date) ==="
    $GATK --java-options "-Xmx6G" Mutect2 \
        -R $REF \
        -I ${BAM_DIR}/${SAMPLE}.sorted.bam \
        $CHRS \
        --native-pair-hmm-threads 4 \
        -O ${OUT_DIR}/${SAMPLE}.unfiltered.vcf.gz

    $GATK --java-options "-Xmx6G" FilterMutectCalls \
        -R $REF \
        -V ${OUT_DIR}/${SAMPLE}.unfiltered.vcf.gz \
        -O ${OUT_DIR}/${SAMPLE}.filtered.vcf.gz

    zcat ${OUT_DIR}/${SAMPLE}.filtered.vcf.gz | grep -E "^#|PASS" > ${OUT_DIR}/${SAMPLE}.somatic.pass.vcf
    echo "=== Finalizada $SAMPLE em $(date) ==="
}

for SAMPLE in SRR8619186 SRR8619133 SRR8619113 SRR8619154; do
    processar_amostra $SAMPLE &
done

wait
echo "--- TODAS AS 4 AMOSTRAS GATK CONCLUIDAS COM SUCESSO! ---"
