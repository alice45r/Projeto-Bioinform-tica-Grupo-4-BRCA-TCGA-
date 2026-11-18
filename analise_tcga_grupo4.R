# 1. Definir pasta de trabalho do repositório local
setwd("~/Documents/trabalho_bioinfo/Projeto-Bioinform-tica-Grupo-4-BRCA-TCGA-")

# 2. Inspecionar o workspace WS.exoma_Rdata que estava no servidor
load("WS.exoma_Rdata")
cat("=== Objetos dentro de WS.exoma_Rdata ===\n")
print(ls())
for (item in ls()) {
  if (item != "item") {
    cat("\n---", item, "---\n")
    str(get(item))
  }
}

# 3. Ler tabela clínica (pulando as 4 linhas iniciais de comentário '#')
clin <- read.delim("data_clinical_patient.txt", skip = 4, header = TRUE,
                   sep = "\t", stringsAsFactors = FALSE)
clin_sub <- clin[!is.na(clin$SUBTYPE) & clin$SUBTYPE != "", c("PATIENT_ID", "SUBTYPE")]

# 4. Ler tabela de mutações somáticas
mut <- read.delim("data_mutations.txt", header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, quote = "", comment.char = "")

# Extrair PATIENT_ID (12 primeiros caracteres de Tumor_Sample_Barcode)
mut$PATIENT_ID <- substr(mut$Tumor_Sample_Barcode, 1, 12)

# Manter apenas pacientes com dados de sequenciamento
clin_seq <- clin_sub[clin_sub$PATIENT_ID %in% unique(mut$PATIENT_ID), ]

# 5. Filtrar mutações funcionais (não-sinónimas)
classes_funcionais <- c(
  "Missense_Mutation", "Nonsense_Mutation",
  "Frame_Shift_Del", "Frame_Shift_Ins",
  "Splice_Site", "In_Frame_Del", "In_Frame_Ins",
  "Translation_Start_Site", "Nonstop_Mutation"
)
mut_func <- mut[mut$Variant_Classification %in% classes_funcionais, ]

# 6. Função para integrar tabelas e calcular prevalência nos 5 subtipos
analisar_subtipos <- function(tabela_mut, tabela_clin) {
  integrado <- merge(tabela_mut[, c("PATIENT_ID", "Tumor_Sample_Barcode", "Hugo_Symbol")],
                     tabela_clin, by = "PATIENT_ID")
  
  # Contar cada gene no máximo 1 vez por amostra
  amostra_gene <- unique(integrado[, c("Tumor_Sample_Barcode", "SUBTYPE", "Hugo_Symbol")])
  
  totais <- table(tabela_clin$SUBTYPE)
  ordem_subtipos <- c("BRCA_LumA", "BRCA_LumB", "BRCA_Basal", "BRCA_Her2", "BRCA_Normal")
  nomes_limpos <- c(
    BRCA_LumA = "Luminal A",
    BRCA_LumB = "Luminal B",
    BRCA_Basal = "Basal-like",
    BRCA_Her2 = "HER2",
    BRCA_Normal = "Normal-like"
  )
  
  lista_res <- list()
  for (sub in ordem_subtipos) {
    dados_sub <- amostra_gene[amostra_gene$SUBTYPE == sub, ]
    contagem <- sort(table(dados_sub$Hugo_Symbol), decreasing = TRUE)
    n_amostras <- as.numeric(totais[sub])
    
    df_sub <- data.frame(
      Codigo_Subtipo = sub,
      Subtipo_Molecular = nomes_limpos[sub],
      Gene = names(contagem),
      Amostras_Mutadas = as.numeric(contagem),
      Total_Amostras_Subtipo = n_amostras,
      Prevalencia_Pct = round((as.numeric(contagem) / n_amostras) * 100, 2),
      stringsAsFactors = FALSE
    )
    lista_res[[sub]] <- df_sub
    
    cat(sprintf("\n=== %s (N = %d amostras) ===\n", nomes_limpos[sub], n_amostras))
    print(head(df_sub[, c("Gene", "Amostras_Mutadas", "Prevalencia_Pct")], 10), row.names = FALSE)
  }
  return(do.call(rbind, lista_res))
}

# Calcular prevalência (com filtro funcional e sem filtro)
resultado_funcional <- analisar_subtipos(mut_func, clin_seq)
resultado_todas <- analisar_subtipos(mut, clin_seq)

# 7. Criar tabela comparativa lado a lado (Top 15 genes x 5 subtipos)
top_genes <- names(sort(tapply(resultado_funcional$Amostras_Mutadas,
                               resultado_funcional$Gene, sum),
                        decreasing = TRUE))[1:15]
sub_filtrados <- resultado_funcional[resultado_funcional$Gene %in% top_genes, ]
tabela_comparativa <- reshape(
  sub_filtrados[, c("Gene", "Subtipo_Molecular", "Prevalencia_Pct")],
  idvar = "Gene", timevar = "Subtipo_Molecular", direction = "wide"
)
colnames(tabela_comparativa) <- gsub("Prevalencia_Pct\\.", "", colnames(tabela_comparativa))
tabela_comparativa[is.na(tabela_comparativa)] <- 0
tabela_comparativa <- tabela_comparativa[match(top_genes, tabela_comparativa$Gene), ]

# 8. Salvar tabelas em CSV
write.csv(resultado_funcional, "TCGA_prevalencia_por_subtipo_funcional.csv", row.names = FALSE)
write.csv(resultado_todas, "TCGA_prevalencia_por_subtipo_todas.csv", row.names = FALSE)
write.csv(tabela_comparativa, "TCGA_tabela_comparativa_top15_genes.csv", row.names = FALSE)

# 9. Gerar gráfico comparativo dos Top 10 genes por subtipo
png("grafico_top10_genes_por_subtipo.png", width = 1600, height = 1000, res = 150)
par(mfrow = c(2, 3), mar = c(4, 6, 3, 2))
cores <- c("Luminal A" = "#2b5c8f", "Luminal B" = "#4ea8de",
           "Basal-like" = "#d9534f", "HER2" = "#f0ad4e", "Normal-like" = "#5cb85c")

for (sub_nome in names(cores)) {
  df_plot <- head(resultado_funcional[resultado_funcional$Subtipo_Molecular == sub_nome, ], 10)
  df_plot <- df_plot[order(df_plot$Prevalencia_Pct), ]
  bp <- barplot(df_plot$Prevalencia_Pct, names.arg = df_plot$Gene, horiz = TRUE,
                las = 1, col = cores[sub_nome],
                xlim = c(0, max(df_plot$Prevalencia_Pct) * 1.25),
                xlab = "Prevalência (%)",
                main = paste0(sub_nome, " (N = ", df_plot$Total_Amostras_Subtipo[1], ")"))
  text(x = df_plot$Prevalencia_Pct + 1.5, y = bp,
       labels = paste0(df_plot$Prevalencia_Pct, "%"), cex = 0.85)
}
dev.off()

# Abrir tabelas na tela para visualização
View(resultado_funcional)
View(tabela_comparativa)