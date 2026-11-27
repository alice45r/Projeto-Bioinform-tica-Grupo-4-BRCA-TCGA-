#!/bin/bash
export MAMBA_ROOT_PREFIX=~/micromamba
mkdir -p ~/dataset/mafs

for SAMPLE in SRR8619186 SRR8619133 SRR8619113 SRR8619154; do
    echo "=== Convertendo $SAMPLE para MAF ==="
    ~/bin/micromamba run -n vcf2maf_env perl ~/micromamba/envs/vcf2maf_env/bin/vcf2maf.pl \
        --input-vcf ~/dataset/mutacoes_somaticas/${SAMPLE}.somatic.pass.vcf \
        --output-maf ~/dataset/mafs/${SAMPLE}.maf \
        --tumor-id ${SAMPLE} \
        --ref-fasta ~/dataset/ref_hg38/hg38.fa \
        --vep-path ~/micromamba/envs/maf_env/bin \
        --online-vep \
        --ncbi-build GRCh38
done

echo "--- CONVERSAO MAF CONCLUIDA! ---"
curl -s -d "GATK e os 4 arquivos MAF estao prontos no servidor!" ntfy.sh/grupo4_usp_bioinfo
