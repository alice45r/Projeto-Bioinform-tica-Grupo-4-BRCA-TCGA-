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

## 🧬 Chamada de Variantes Somáticas (`GATK Mutect2`) e Análise (`maftools`)

Após o mapeamento contra o genoma de referência `hg38`, realizou-se a chamada de variantes somáticas em modo *tumor-only* utilizando o **GATK4 (`Mutect2` + `FilterMutectCalls`)**, seguida pela conversão das variantes filtradas (`PASS`) para o formato **MAF** (*Mutation Annotation Format*) e análise exploratória no R com o pacote **`maftools`** (`analise_maftools.R`).

### 📊 Resumo de Mutações por Linhagem Celular

| Amostra (`Tumor_Sample_Barcode`) | `Missense_Mutation` | `Frame_Shift_Del` | `Frame_Shift_Ins` | Total de Mutações |
| :--- | :---: | :---: | :---: | :---: |
| **SRR8619113** | 190.834 | 5.830 | 3.863 | **200.527** |
| **SRR8619154** | 191.563 | 5.144 | 3.618 | **200.325** |
| **SRR8619186** | 186.109 | 5.442 | 3.769 | **195.320** |
| **SRR8619133** | 156.868 | 4.633 | 3.142 | **164.643** |

### 🔬 Principais Achados e Inferência de Subtipos Moleculares (Comparação com TCGA-BRCA)

* **Predomínio de Variantes *Missense*:** Cerca de **95%** das variantes somáticas identificadas nas quatro linhagens correspondem a mutações não-sinônimas (`Missense_Mutation`), enquanto inserções e deleções fora da grelha de leitura (`Frame_Shift_Ins` e `Frame_Shift_Del`) representam aproximadamente **5%**.
* **Alta Carga Mutacional (*Tumor-Only*):** A ausência de tecido normal pareado para filtragem germinativa resulta numa contagem total elevada (~164k a ~200k variantes/amostra), padrão típico de linhagens celulares imortalizadas.
* **Perfil de Instabilidade Genómica (`SRR8619113`, `SRR8619154` e `SRR8619186`):** Apresentam carga mutacional elevada e homogénea (~195k–200k variantes), comportamento compatível com subtipos de maior instabilidade cromossómica no TCGA-BRCA, como **Basal-like (Triplo-Negativo)** e **HER2-enriquecido**.
* **Perfil de Menor Carga Relativa (`SRR8619133`):** Apresenta cerca de 35 mil variantes a menos que as demais (164.643 no total), perfil mais próximo aos subtipos luminais (**Luminal A / Luminal B**).

### 📁 Arquivos Adicionados (`mafs/` e `resultados_maftools/`)

* `analise_maftools.R`: Script em R responsável pela união dos arquivos `.maf` (`merge_mafs`) e geração das métricas e visualizações comparativas.
* `carcinoma_mamario_grupo4_sampleSummary.txt`: Tabela quantitativa detalhada das variantes por amostra.
* `01_resumo_variantes.pdf`, `02_oncoplot_top20.pdf` e `03_titv_spectrum.pdf`: Gráficos de distribuição de classes de variantes, frequência mutacional e espectro de transições/transversões (`Ti/Tv`).
