# Projeto de Bioinformática — Grupo 4 (TCGA-BRCA)

Repositório com os scripts, controles de qualidade e resultados de alinhamento de dados genômicos e de sequenciamento completo de exoma (WES) de carcinoma mamário (TCGA-BRCA), desenvolvidos na disciplina **RIB 0109**.

## Estrutura do Repositório

* `scripts/` e `pipeline_exoma_grupo4.sh`: Scripts em Bash para download via `sra-toolkit`, controle de qualidade (`FastQC`), indexação do genoma de referência `hg38` e alinhamento com `BWA-MEM` + `SAMtools`.
* `fastqc_reports/`: Relatórios de qualidade das leituras brutas (`FastQC`).
* `mapeamento_hg38/`: Estatísticas de alinhamento (`samtools flagstat`) das 4 amostras de exoma contra o genoma humano `GRCh38/hg38`.
* `analise_tcga_grupo4.R` e arquivos `.csv` / `.png`: Análise de prevalência de mutações por subtipo molecular do TCGA-BRCA.

## Resumo do Alinhamento de Exomas (BWA-MEM / hg38)

| Amostra (SRA) | Total de Leituras (QC-passed) | Leituras Mapeadas (`hg38`) | Pares Corretos (*Properly Paired*) | *Singletons* |
| :--- | :---: | :---: | :---: | :---: |
| **`SRR8619186`** | 117.358.731 | 114.440.841 (97,51%) | 111.481.235 (94,99%) | 1.570.565 (1,34%) |
| **`SRR8619133`** | 67.965.236  | 66.658.077 (98,08%)  | 64.585.514 (95,03%)  | 976.483 (1,44%)   |
| **`SRR8619113`** | 155.800.551 | 151.273.621 (97,09%) | 147.017.567 (94,36%) | 2.253.733 (1,45%) |
| **`SRR8619154`** | 117.319.974 | 115.034.458 (98,05%) | 112.795.779 (96,14%) | 1.224.963 (1,04%) |
