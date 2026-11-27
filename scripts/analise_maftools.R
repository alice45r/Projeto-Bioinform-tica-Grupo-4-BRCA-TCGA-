# Carregar biblioteca
if (!requireNamespace("maftools", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")
  BiocManager::install("maftools", update = FALSE, ask = FALSE)
}
library(maftools)

# Criar pasta para os gráficos e tabelas finais
dir.create("~/dataset/resultados_maftools", showWarnings = FALSE, recursive = TRUE)
setwd("~/dataset/resultados_maftools")

# Listar e unir os 4 arquivos .maf gerados pelo pipeline
arquivos_maf <- list.files("~/dataset/mafs", pattern = "\\.maf$", full.names = TRUE)
cat("Arquivos MAF encontrados:\n")
print(arquivos_maf)

laml <- merge_mafs(mafs = arquivos_maf, verbose = TRUE)

# Salvar tabelas resumo (por amostra e por gene)
write.mafSummary(maf = laml, basename = "carcinoma_mamario_grupo4")

# 1. Resumo geral das variantes (SNP, INDEL, classes de mutação, top genes)
pdf("01_resumo_variantes.pdf", width = 12, height = 8)
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = "median", dashboard = TRUE, titvRaw = FALSE)
dev.off()

# 2. Oncoplot com os 20 genes mais mutados
pdf("02_oncoplot_top20.pdf", width = 11, height = 8)
oncoplot(maf = laml, top = 20, showTumorSampleBarcodes = TRUE)
dev.off()

# 3. Transições e Transversões (Ti/Tv)
pdf("03_titv_spectrum.pdf", width = 10, height = 6)
laml.titv <- titv(maf = laml, plot = FALSE, useSyn = TRUE)
plotTiTv(res = laml.titv)
dev.off()

# 4. Vias oncogênicas enriquecidas
pdf("04_vias_oncogenicas.pdf", width = 10, height = 6)
OncogenicPathways(maf = laml)
dev.off()

# 5. Interações somáticas (co-ocorrência / exclusividade mútua)
pdf("05_interacoes_somaticas.pdf", width = 9, height = 9)
somaticInteractions(maf = laml, top = 15, pvalue = c(0.05, 0.1))
dev.off()
cat("\n--- ANALISE MAFTOOLS CONCLUIDA! Resultados salvos em ~/dataset/resultados_maftools ---\n")
