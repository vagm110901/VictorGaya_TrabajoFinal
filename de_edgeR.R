
# Cargar librerías necesarias
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(org.Mm.eg.db)  # Base de datos de genes humanos
library(DOSE)
library(enrichplot)
library(DT)
library(biomaRt)

colores <- c(
  "FSM" = "#6BAED6",
  "ISM" = "#FC8D59",
  "NIC" = "#78C679",
  "NNC" = "#EE6A50",
  "Genic Genomic" = "#969696",
  "Antisense" = "#66C2A4",
  "Fusion" = "goldenrod1",
  "Intergenic" = "darksalmon",
  "Genic Intron" = "#41B6C4"
)

# Leer la matriz de conteos
data <- read.csv("./isoform_counts_matrix.tsv", row.names = 1, sep = '\t')

annotations <- data[, c("structural_category", "associated_gene")] # Extraer metadatos

# Extraer los conteos y preparar el diseño experimental
counts <- data[, 1:4]  # Conteos en las 4 primeras columnas
conditions <- c("IR", "IR_control", "MI", "MI_control")  # Condiciones

# Crear objeto DGEList
group <- factor(conditions)
dge <- DGEList(counts = counts, group = group)

# Normalización TMM
dge <- calcNormFactors(dge)

dge$common.dispersion <- 0.1

# Realizar el análisis exacto para cada comparación

#logFC positivo = condición patológica
#logFC negativo = condición control

# Comparación 1: IR vs IR_control
et_IR <- exactTest(dge, pair = c("IR", "IR_control"))
results_IR <- topTags(et_IR, n = Inf)$table
write.table(results_IR, file = "DEA_IR_vs_IR_control.tsv", sep = "\t", row.names = FALSE)

results_IR <- merge(results_IR, annotations, by.x = "row.names", by.y = "row.names", all.x = TRUE)
colnames(results_IR)[1] <- "transcript_id"  # Renombrar la columna de IDs de transcritos
write.table(results_IR, file = "DEA_IR_vs_IR_control_with_annotations.tsv", sep = "\t", row.names = FALSE)

results_IR_s <- results_IR[which(results_IR$PValue < 0.01),]
genes_IR_s <- results_IR_s[,c("associated_gene","structural_category")]
table(genes_IR_s$structural_category)

# novel Genes
genes_IR_s_nG <- results_IR_s[grep("novelGene*", results_IR_s$associated_gene),]
NnovelGenes_IR <- length(grep("novelGene", results_IR_s$associated_gene))/length(results_IR_s$associated_gene)*100
print(paste0("Hay ", round(NnovelGenes_IR, 2), "% novel genes en los genes significativos de IR."))

# plot categorias estructurales
genes_IR_filtrados <- genes_IR_s %>% 
  filter(!structural_category %in% c("full-splice_match;full-splice_match", "fusion"))

ggplot(genes_IR_filtrados, aes(x = structural_category)) +
  geom_bar(fill = c("#66C2A4","#6BAED6","#969696","#41B6C4","#FC8D59","darksalmon","#78C679","#EE6A50")) +
  theme_minimal() +
  labs(title = "Distribución de Categorías Estructurales \n(reperfusión de isquemia posterior)", 
       x = "Categorías", 
       y = "Frecuencia") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Comparación 2: MI vs MI_control
et_MI <- exactTest(dge, pair = c("MI", "MI_control"))
results_MI <- topTags(et_MI, n = Inf)$table
write.table(results_MI, file = "DEA_MI_vs_MI_control.tsv", sep = "\t", row.names = FALSE)

results_MI <- merge(results_MI, annotations, by.x = "row.names", by.y = "row.names", all.x = TRUE)
colnames(results_MI)[1] <- "transcript_id"  # Renombrar la columna de IDs de transcritos
write.table(results_MI, file = "DEA_MI_vs_MI_control_with_annotations.tsv", sep = "\t", row.names = FALSE)

results_MI_s <- results_MI[which(results_MI$PValue < 0.01),]
genes_MI_s <- results_MI_s[,c("associated_gene","structural_category")]
table(genes_MI_s$structural_category)

# novel Genes
genes_MI_s_nG <- results_MI_s[grep("novelGene*", results_MI_s$associated_gene),]
NnovelGenes_MI <- length(grep("novelGene", results_MI_s$associated_gene))/length(results_MI_s$associated_gene)*100
print(paste0("Hay ", round(NnovelGenes_MI, 2), "% novel genes en los genes significativos de MI."))

# plot categorias estructurales
genes_MI_filtrados <- genes_MI_s %>% 
  filter(!structural_category %in% c("full-splice_match;full-splice_match", "fusion", "incomplete-splice_match;incomplete-splice_match"))  

ggplot(genes_MI_filtrados, aes(x = structural_category)) +
  geom_bar(fill = c("#66C2A4","#6BAED6","#969696","#41B6C4","#FC8D59","darksalmon","#78C679")) +
  theme_minimal() +
  labs(title = "Distribución de Categorías Estructurales \n(infarto de miocardio)", 
       x = "Categorías", 
       y = "Frecuencia") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

genes_IRup <- results_IR_s[which(results_IR_s$logFC > 0),"associated_gene"]
write.table(genes_IRup, file = "genes_IRup.txt", sep = " ", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
NnovelGenes_IRup <- length(
  grep("novelGene", results_IR_s[which(results_IR_s$logFC > 0),"associated_gene"])
  )/length(results_IR_s[which(results_IR_s$logFC > 0),"associated_gene"])*100
print(paste0("Hay ", round(NnovelGenes_IRup, 2), "% novel genes en los genes up-regulated de IR."))


genes_IRdown <- results_IR_s[which(results_IR_s$logFC < 0),"associated_gene"]
write.table(genes_IRdown, file = "genes_IRdown.txt", sep = " ", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
NnovelGenes_IRdown <- length(
  grep("novelGene", results_IR_s[which(results_IR_s$logFC < 0),"associated_gene"])
)/length(results_IR_s[which(results_IR_s$logFC < 0),"associated_gene"])*100
print(paste0("Hay ", round(NnovelGenes_IRdown, 2), "% novel genes en los genes down-regulated de IR."))


genes_MIup <- results_MI_s[which(results_MI_s$logFC > 0),"associated_gene"]
write.table(genes_MIup, file = "genes_MIup.txt", sep = " ", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
NnovelGenes_MIup <- length(
  grep("novelGene", results_MI_s[which(results_MI_s$logFC > 0),"associated_gene"])
)/length(results_MI_s[which(results_MI_s$logFC > 0),"associated_gene"])*100
print(paste0("Hay ", round(NnovelGenes_MIup, 2), "% novel genes en los genes up-regulated de MI."))


genes_MIdown <- results_MI_s[which(results_MI_s$logFC < 0),"associated_gene"]
write.table(genes_MIdown, file = "genes_MIdown.txt", sep = " ", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
NnovelGenes_MIdown <- length(
  grep("novelGene", results_MI_s[which(results_MI_s$logFC < 0),"associated_gene"])
)/length(results_MI_s[which(results_MI_s$logFC < 0),"associated_gene"])*100
print(paste0("Hay ", round(NnovelGenes_MIdown, 2), "% novel genes en los genes down-regulated de MI."))

