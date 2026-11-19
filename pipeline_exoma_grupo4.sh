mkdir -p ~/dataset/sra_exoma
mkdir -p ~/dataset/fastqc_reports
cd ~/dataset/sra_exoma || exit 1

for SRR in SRR8619186 SRR8619133 SRR8619113 SRR8619154; do
    echo "=== [1/2] A transferir exoma emparelhado $SRR em $(date) ==="
    fastq-dump --split-files --gzip "$SRR"
    echo "=== Transferencia de $SRR concluida em $(date) ==="
done

echo "=== [2/2] A gerar relatorios de qualidade com FastQC em $(date) ==="
fastqc -o ~/dataset/fastqc_reports -t 2 *.fastq.gz
echo "=== Pipeline concluido com sucesso em $(date) ==="
