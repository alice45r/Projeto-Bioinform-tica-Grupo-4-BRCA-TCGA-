set -e
mkdir -p ~/dataset/ref_hg38 ~/dataset/mapeamento_hg38
cd ~/dataset/ref_hg38

if [ ! -f "hg38.fa" ]; then
    echo "=== [1/3] A descarregar genoma de referencia hg38 em $(date) ==="
    wget --no-check-certificate -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
    echo "=== A descompactar hg38.fa.gz em $(date) ==="
    gunzip -f hg38.fa.gz
fi

if [ ! -f "hg38.fa.bwt" ]; then
    echo "=== [2/3] A criar indice BWA para hg38 (demora ~45-60 min) em $(date) ==="
    bwa index hg38.fa
fi

cd ~/dataset/mapeamento_hg38
REF=~/dataset/ref_hg38/hg38.fa
FASTQ_DIR=~/dataset/sra_exoma

for SRR in SRR8619186 SRR8619133 SRR8619113 SRR8619154; do
    echo "=== [3/3] A mapear amostra $SRR com BWA-MEM (8 nucleos) em $(date) ==="
    R1="${FASTQ_DIR}/${SRR}_1.fastq.gz"
    R2="${FASTQ_DIR}/${SRR}_2.fastq.gz"
    RG="@RG\tID:${SRR}\tSM:${SRR}\tPL:ILLUMINA"

    bwa mem -t 8 -R "$RG" "$REF" "$R1" "$R2" | samtools sort -@ 4 -o "${SRR}.sorted.bam" -
    samtools index "${SRR}.sorted.bam"
    samtools flagstat "${SRR}.sorted.bam" > "${SRR}_flagstat.txt"
    echo "=== Amostra $SRR concluida em $(date) ==="
done

echo "=== Mapeamento concluido com sucesso em $(date) ==="
