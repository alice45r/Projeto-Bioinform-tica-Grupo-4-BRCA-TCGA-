#!/bin/bash
export PATH="/home/grupo4/sratk/bin:$PATH"
mkdir -p ~/dataset/sra_exoma ~/dataset/fastqc_reports
cd ~/dataset/sra_exoma || exit 1

for SRR in SRR8619186 SRR8619133 SRR8619113 SRR8619154; do
    echo "=== [1/3] Baixando SRA bruto $SRR em $(date) ==="
    prefetch "$SRR" -O .
    echo "=== [2/3] Extraindo FASTQ de $SRR (8 nucleos) em $(date) ==="
    fasterq-dump "./$SRR" --split-files -e 8 -O .
    rm -rf "./$SRR"
    echo "=== Compactando $SRR em $(date) ==="
    gzip -f "${SRR}_1.fastq" &
    gzip -f "${SRR}_2.fastq" &
    wait
    echo "=== Amostra $SRR concluida em $(date) ==="
done

echo "=== [3/3] Gerando relatorios FastQC em $(date) ==="
fastqc -o ~/dataset/fastqc_reports -t 8 *.fastq.gz
echo "=== Pipeline concluido com sucesso em $(date) ==="
