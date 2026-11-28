# 🧬 Projeto de Bioinformática — Grupo 4 (TCGA-BRCA)

Repositório completo contendo os scripts, controles de qualidade, estatísticas de alinhamento, chamada de variantes somáticas e análises comparativas de Sequenciamento Completo de Exoma (WES) de carcinoma mamário humano, desenvolvidos na disciplina **RIB 0109** (FMRP-USP) sob orientação do Prof. Dr. Daniel Guimarães Tiezzi.

---

## 🎯 Objetivo do Projeto

Caracterizar o perfil mutacional de quatro linhagens celulares humanas de carcinoma mamário (`SRR8619113`, `SRR8619133`, `SRR8619154` e `SRR8619186`) a partir de dados brutos de exoma (SRA/NCBI) e comparar a carga e o espectro de mutações somáticas com os perfis moleculares clínicos da coorte **TCGA-BRCA** (*Luminal A*, *Luminal B*, *HER2-enriquecido* e *Basal-like*).

---

## 🗂️ Estrutura do Repositório

* `scripts/` e `pipeline_exoma_grupo4.sh`: Scripts em Bash para aquisição de dados (`sra-toolkit`), controle de qualidade (`FastQC`), indexação do genoma humano e alinhamento (`BWA-MEM` + `SAMtools`).
* `fastqc_reports/`: Relatórios HTML/ZIP de qualidade das leituras brutas gerados pelo `FastQC`.
* `mapeamento_hg38/`: Estatísticas completas de alinhamento (`samtools flagstat`) das 4 amostras contra o genoma de referência `GRCh38/hg38`.
* `rodar_gatk_rapido.sh` e `converter_maf.sh`: Pipeline de chamada de variantes somáticas com **GATK4 (`Mutect2` + `FilterMutectCalls`)** e conversão automatizada dos arquivos `.vcf.gz` filtrados (`PASS`) para o formato **MAF**.
* `mafs/`: Arquivos finais no formato *Mutation Annotation Format* (`SRR8619113.maf`, `SRR8619133.maf`, `SRR8619154.maf` e `SRR8619186.maf`).
* `analise_maftools.R` e `resultados_maftools/`: Script em R (`Bioconductor/maftools`) e saídas gráficas/tabulares consolidadas das 505.570 variantes processadas.
* `analise_tcga_grupo4.R` e arquivos `.csv` / `.png`: Análise de prevalência de mutações e genes recorrentes por subtipo molecular na coorte clínica **TCGA-BRCA**.

---

## ⚙️ Fluxo Metodológico Completo (Passo a Passo)

### Etapa 1 — Aquisição de Dados Brutos e Controle de Qualidade (`SRA-Toolkit` + `FastQC`)
1. Descarregamento dos arquivos `.sra` das quatro amostras de exoma humano (`SRR8619113`, `SRR8619133`, `SRR8619154`, `SRR8619186`) e conversão para arquivos `FASTQ` *paired-end* (`fasterq-dump`).
2. Avaliação da qualidade das bases (escore Phred), conteúdo GC e presença de adaptadores com `FastQC`, confirmando qualidade elevada em todas as bibliotecas sem necessidade de *trimming* agressivo.

### Etapa 2 — Alinhamento Contra o Genoma de Referência `GRCh38/hg38` (`BWA-MEM` + `SAMtools`)
1. Indexação do genoma humano de referência (`hg38.fa`) com `bwa index` e `samtools faidx`.
2. Mapeamento das leituras *paired-end* com `bwa mem` (com atribuição de *Read Groups* `@RG` compatíveis com o GATK).
3. Conversão SAM $\rightarrow$ BAM, ordenação por coordenadas genómicas (`samtools sort`), indexação (`samtools index`) e extração de métricas de mapeamento (`samtools flagstat`).

| Amostra (SRA) | Total de Leituras (*QC-passed*) | Leituras Mapeadas (`hg38`) | Pares Corretos (*Properly Paired*) | *Singletons* |
| :--- | :---: | :---: | :---: | :---: |
| **SRR8619186** | 117.358.731 | 114.448.841 (97,51%) | 111.481.235 (94,99%) | 1.578.565 (1,34%) |
| **SRR8619133** | 87.965.238 | 86.858.877 (98,80%) | 84.585.314 (95,03%) | 978.483 (1,44%) |
| **SRR8619113** | 155.800.551 | 151.273.821 (97,09%) | 147.017.567 (94,36%) | 2.253.733 (1,45%) |
| **SRR8619154** | 117.319.974 | 115.034.458 (98,05%) | 112.795.779 (96,14%) | 1.224.963 (1,84%) |

### Etapa 3 — Chamada de Variantes Somáticas (`GATK4 Mutect2` + `FilterMutectCalls`)
1. Criação do dicionário de sequência (`CreateSequenceDictionary`) para `hg38.fa`.
2. Execução paralela no servidor da USP (`DGOVM1`) do **`GATK Mutect2`** em modo *tumor-only* restrito aos cromossomos canónicos (`chr1`–`chr22`, `chrX`, `chrY`), gerando os arquivos `.unfiltered.vcf.gz`.
3. Aplicação do filtro estatístico **`GATK FilterMutectCalls`** para distinguir variantes somáticas verdadeiras (`FILTER=PASS`) de artefactos de sequenciamento e ruído de mapeamento.

### Etapa 4 — Conversão de VCF Filtrado para MAF (`converter_maf.sh`)
* Para contornar limitações de dependências Perl/VEP do utilitário clássico `vcf2maf` no ambiente do servidor, desenvolveu-se um conversor automatizado em Python/Bash (`converter_maf.sh`) que extrai todas as variantes aprovadas (`PASS`) dos arquivos `.filtered.vcf.gz`, classifica substituições de nucleótido único e inserções/deleções (`Missense_Mutation`, `Frame_Shift_Del`, `Frame_Shift_Ins`) e estrutura as colunas padrão exigidas pelo **`maftools`** (`Hugo_Symbol`, `Chromosome`, `Start_Position`, `End_Position`, `Variant_Classification`, `Variant_Type`, `Reference_Allele`, `Tumor_Seq_Allele2`, `Tumor_Sample_Barcode`).

### Etapa 5 — Análise Mutacional Integrada (`maftools`) e Comparação com o TCGA-BRCA
1. União dos quatro arquivos `.maf` (`merge_mafs`) em ambiente R (`analise_maftools.R`), totalizando **505.570 variantes** consolidadas.
2. Geração de tabelas resumo (`carcinoma_mamario_grupo4_sampleSummary.txt`) e visualizações (`01_resumo_variantes.pdf`, `02_oncoplot_top20.pdf`, `03_titv_spectrum.pdf`).
3. Contraste dos perfis mutacionais obtidos com os dados públicos da coorte **TCGA-BRCA** (`analise_tcga_grupo4.R`).

---

## 📊 Resultados Consolidados das 4 Linhagens Celulares

| Amostra (`Tumor_Sample_Barcode`) | `Missense_Mutation` | `Frame_Shift_Del` | `Frame_Shift_Ins` | Total de Mutações |
| :--- | :---: | :---: | :---: | :---: |
| **SRR8619113** | 190.834 | 5.830 | 3.863 | **200.527** |
| **SRR8619154** | 191.563 | 5.144 | 3.618 | **200.325** |
| **SRR8619186** | 186.109 | 5.442 | 3.769 | **195.320** |
| **SRR8619133** | 156.868 | 4.633 | 3.142 | **164.643** |

---

## 🔬 Discussão Biológica e Comparação com Subtipos do TCGA-BRCA

* **Alta Eficiência de Mapeamento:** Todas as amostras apresentaram taxa de alinhamento superior a **97%** contra o `hg38` e mais de **94%** de pares devidamente alinhados (*properly paired*), garantindo cobertura robusta para a chamada de variantes.
* **Elevada Carga Mutacional em Modo *Tumor-Only*:** A contagem total entre ~164 mil e ~200 mil variantes por linhagem reflete a ausência de amostra normal pareada do mesmo doador (mantendo variantes germinativas raras não filtradas) somada ao acúmulo mutacional característico de linhagens celulares imortalizadas.
* **Predomínio de Mutações *Missense* (~95%):** Em todas as amostras, as substituições não-sinónimas (`Missense_Mutation`) correspondem a ~95% do total de alterações, enquanto pequenas inserções e deleções fora da grelha de leitura (`Frame_Shift_Ins` e `Frame_Shift_Del`) representam ~5%.
* **Estratificação por Instabilidade Genómica e Relação com o TCGA-BRCA:**
  * **Grupo de Alta Carga Mutacional (`SRR8619113`, `SRR8619154` e `SRR8619186`):** Apresentam carga elevada e homogénea (~195k–200k variantes), padrão concordante com subtipos de carcinoma mamário de maior instabilidade genómica na coorte TCGA-BRCA, como **Basal-like (Triplo-Negativo)** e **HER2-enriquecido** (frequentemente associados a mutações em *TP53* e deficiência em vias de reparo de DNA).
  * **Linhagem de Menor Carga Relativa (`SRR8619133`):** Apresenta cerca de 35 mil variantes a menos que as restantes (164.643 variantes), perfil mais próximo aos tumores luminais (**Luminal A / Luminal B**), que exibem menor carga mutacional global no TCGA.

## 📚 Referências Bibliográficas e Documentação

* **Artigos de Referência (Linhagens Celulares de Carcinoma Mamário):**
  * [Nature — Genomic landscape of breast cancer cell lines](https://www.nature.com/articles/35021093)
  * [Cell Genomics — Molecular characterization of breast cancer models](https://www.cell.com/cell-genomics/fulltext/S2666-979X(21)00083-5)
* **Documentação das Ferramentas:**
  * [FastQC — Babraham Bioinformatics](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
  * [GATK Best Practices Workflows — Broad Institute](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows)
  * [Formato MAF / vcf2maf — MSKCC](https://github.com/mskcc/vcf2maf)
