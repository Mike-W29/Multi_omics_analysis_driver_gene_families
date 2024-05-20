# Load necessary libraries
library(data.table)
library(ggvenn)
library(dplyr)

# Check and load CNV data
cnv_file <- '~/all_thresholded.by_genes.txt'  # all_thresholded.by_genes.txt from GISTIC2.0
if (!file.exists(cnv_file)) stop("CNV data file not found.")
CNV_gene <- fread(cnv_file, header = TRUE)
CNV_gene <- as.data.frame(CNV_gene)
ID_switch <- CNV_gene[, 1:3]
rownames(CNV_gene) <- CNV_gene$`Gene Symbol`
CNV_gene <- CNV_gene[, -1:-3]

# Check column frequency
cat("Column Frequency:\n")
print(table(substr(colnames(CNV_gene), 14, 16)))

# Calculate gene frequency
CNV_gene_freq <- data.frame(
  Gain = rowSums(CNV_gene >= 1),
  Loss = rowSums(CNV_gene <= -1),
  None = rowSums(CNV_gene == 0)
)
CNV_gene_freq <- CNV_gene_freq / ncol(CNV_gene)
rownames(CNV_gene_freq) <- rownames(CNV_gene)

# Filter Gain and Loss genes based on 20% threshold
CNV_gene_Gain <- rownames(CNV_gene_freq)[CNV_gene_freq$Gain > 0.2]
CNV_gene_Loss <- rownames(CNV_gene_freq)[CNV_gene_freq$Loss > 0.2]

# Load RNA data
rna_file <- '~/RNA_DEG.RData' # RNA differential analysis and expression matrix
if (!file.exists(rna_file)) stop("RNA data file not found.")
load(rna_file)

# Filter RNA data
merge_RNA_DEG <- merge_RNA_DEG[merge_RNA_DEG$Status != 'None', ]

# Align samples
co_sample <- intersect(colnames(merge_RNA_ComBat), substr(colnames(CNV_gene), 1, 15))

RNA_T <- merge_RNA_ComBat[, co_sample, drop = FALSE]
CNV_T <- CNV_gene[, substr(colnames(CNV_gene), 1, 15) %in% co_sample, drop = FALSE]

# Ensure sample order is the same
RNA_T <- RNA_T[, order(colnames(RNA_T))]
CNV_T <- CNV_T[, order(substr(colnames(CNV_T), 1, 15))]

# Initialize result data frame
CNV_VS_DMG <- data.frame(
  gene = rownames(RNA_T),
  logFC = numeric(nrow(RNA_T)),
  Spearman = numeric(nrow(RNA_T)),
  pvalue = numeric(nrow(RNA_T)),
  statusRNA = character(nrow(RNA_T)),
  stringsAsFactors = FALSE
)

# Calculate Spearman correlation
for (i in seq_len(nrow(RNA_T))) {
  cat('Processing row', i, 'of', nrow(RNA_T), '\n')
  a <- match(CNV_VS_DMG$gene[i], rownames(CNV_T))
  if (!is.na(a)) {
    p <- cor.test(as.numeric(CNV_T[a, ]), as.numeric(RNA_T[i, ]), method = 'spearman', exact = FALSE)
    CNV_VS_DMG$Spearman[i] <- p$estimate
    CNV_VS_DMG$pvalue[i] <- p$p.value
  }
}

CNV_VS_DMG <- CNV_VS_DMG[CNV_VS_DMG$pvalue != 0, ]
CNV_VS_DMG$adjustp <- p.adjust(CNV_VS_DMG$pvalue, method = 'BH')

# Fill in logFC and status information
for (i in seq_len(nrow(CNV_VS_DMG))) {
  a <- match(CNV_VS_DMG$gene[i], merge_RNA_DEG$Gene)
  if (!is.na(a)) {
    CNV_VS_DMG$logFC[i] <- merge_RNA_DEG$LogFC[a]
    CNV_VS_DMG$statusRNA[i] <- merge_RNA_DEG$Status[a]
  }
}

CNV_VS_DMG <- CNV_VS_DMG[CNV_VS_DMG$statusRNA != "", ]
CNV_VS_DMG <- CNV_VS_DMG[CNV_VS_DMG$adjustp <= 0.05, ]

# Filter Driver genes
CNV_Drivergene_Gain <- CNV_gene_Gain[CNV_gene_Gain %in% CNV_VS_DMG$gene]
CNV_Drivergene_Loss <- CNV_gene_Loss[CNV_gene_Loss %in% CNV_VS_DMG$gene]

# Save results
write.csv(CNV_VS_DMG, '~/xxxxx.csv', row.names = FALSE)
save(CNV_VS_DMG, CNV_Drivergene_Gain, CNV_Drivergene_Loss, file = '~/XXXXX.RData')

# Plotting function
plot_venn <- function(set1, set2, set3, title) {
  ggvenn(data = list(
    "CNV" = set1,
    "RNA regulation" = set2,
    "Correlation" = set3
  ), columns = c("CNV", "RNA regulation", "Correlation"),
  fill_color = c("#7F7FFF", '#FFFF7F', '#7FFF7F'),
  set_name_color = c("#1818FF", '#CACA9A', '#0C860C'),
  set_name_size = 9, text_size = 7) + ggtitle(title)
}

# CNV Loss
set1 <- CNV_Drivergene_Loss
set2 <- CNV_VS_DMG[CNV_VS_DMG$statusRNA == 'Down', 'gene']
set3 <- CNV_VS_DMG[CNV_VS_DMG$Spearman >= 0.2, 'gene']
plot_venn(set1, set2, set3, "CNV Loss")

# CNV Gain
set1 <- CNV_Drivergene_Gain
set2 <- CNV_VS_DMG[CNV_VS_DMG$statusRNA == 'Up', 'gene']
set3 <- CNV_VS_DMG[CNV_VS_DMG$Spearman >= 0.2, 'gene']
plot_venn(set1, set2, set3, "CNV Gain")
