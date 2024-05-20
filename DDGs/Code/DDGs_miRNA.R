
# Load data
# (MiRNA expression matrix, mRNA expression matrix, miRNA differential gene list, mRNA differential gene list)
load('~/RNA_DEG.RData')
load('~/miRNA_DEG.RData')

# Filter merge_RNA_DEG data
merge_RNA_DEG <- merge_RNA_DEG[merge_RNA_DEG$Status != 'None', ]
gene <- merge_RNA_DEG$Gene

# Load multiMiR package and get miRNA and mRNA interaction data
library(multiMiR)
mirna_DEG <- miRNA_DEG$Gene[miRNA_DEG$Status != 'None']

miRNA_gene_val <- get_multimir(org = 'Homo sapiens', target = gene, table = 'validated')
miRNA_gene_val <- miRNA_gene_val@data

miRNA_gene_pre <- get_multimir(org = 'Homo sapiens', target = gene, table = 'predicted', predicted.site = 'all')
miRNA_gene_pre <- miRNA_gene_pre@data

# Save intermediate results
save.image('~/miRNA_Driver.RData')
load('~/miRNA_Driver.RData')

# Select required columns
miRNA_gene_val <- dplyr::select(miRNA_gene_val, c('database', 'mature_mirna_id', 'target_symbol', 'type'))
miRNA_gene_pre <- dplyr::select(miRNA_gene_pre, c('database', 'mature_mirna_id', 'target_symbol', 'type'))

# Merge data
miRNA_VS_mRNA <- rbind(miRNA_gene_val, miRNA_gene_pre)

# Filter for differentially expressed miRNA and mRNA
miRNA_VS_mRNA <- miRNA_VS_mRNA[miRNA_VS_mRNA$mature_mirna_id %in% mirna_DEG & miRNA_VS_mRNA$target_symbol %in% gene, ]

# Log transform data
miRNA <- log2(miRNA + 1)
merge_RNA_ComBat <- log2(merge_RNA_ComBat + 1)

# Align samples
co_sample <- intersect(colnames(miRNA), colnames(merge_RNA_ComBat))

RNA_T <- merge_RNA_ComBat[, co_sample]
miRNA_T <- miRNA[, co_sample]

# Initialize columns in miRNA_VS_mRNA
miRNA_VS_mRNA <- transform(miRNA_VS_mRNA,
                           Spearman = 0, p.value = 0, p.adjust = 0,
                           miRNAlogFC = 0, mRNAlogFC = 0,
                           miRNAstatus = 0, RNAstatus = 0,
                           pair = paste(mature_mirna_id, target_symbol, sep = '-'))

# Remove duplicates
miRNA_VS_mRNA <- miRNA_VS_mRNA[!duplicated(miRNA_VS_mRNA$pair), ]

# Calculate Spearman correlation
for (i in 1:nrow(miRNA_VS_mRNA)) {
  cat('Processing row', i, 'of', nrow(miRNA_VS_mRNA), '\n')
  a <- match(miRNA_VS_mRNA$target_symbol[i], rownames(RNA_T))
  b <- match(miRNA_VS_mRNA$mature_mirna_id[i], rownames(miRNA_T))
  p <- cor.test(as.numeric(RNA_T[a, ]), as.numeric(miRNA_T[b, ]), method = "spearman", exact = FALSE)
  miRNA_VS_mRNA$Spearman[i] <- p$estimate
  miRNA_VS_mRNA$p.value[i] <- p$p.value
}

miRNA_VS_mRNA$p.adjust <- p.adjust(miRNA_VS_mRNA$p.value, method = 'BH')

# Fill in logFC and status information
for (i in 1:nrow(miRNA_VS_mRNA)) {
  cat('Processing row', i, 'of', nrow(miRNA_VS_mRNA), '\n')
  a <- match(miRNA_VS_mRNA$target_symbol[i], merge_RNA_DEG$Gene)
  b <- match(miRNA_VS_mRNA$mature_mirna_id[i], miRNA_DEG$Gene)
  if (!is.na(a) & !is.na(b)) {
    miRNA_VS_mRNA$miRNAlogFC[i] <- miRNA_DEG$LogFC[b]
    miRNA_VS_mRNA$miRNAstatus[i] <- miRNA_DEG$Status[b]
    miRNA_VS_mRNA$mRNAlogFC[i] <- merge_RNA_DEG$LogFC[a]
    miRNA_VS_mRNA$RNAstatus[i] <- merge_RNA_DEG$Status[a]
  }
}

# Filter for significant miRNA-mRNA pairs
miRNA_VS_mRNA <- miRNA_VS_mRNA[miRNA_VS_mRNA$p.adjust < 0.05, ]

# Save results
save(miRNA_VS_mRNA, merge_RNA_ComBat, file = '~/miRNA_VS_mRNA_result.RData')

# Plot Venn diagram
library(ggvenn)

# miRNA upregulation, mRNA downregulation, negative correlation
set1 <- miRNA_VS_mRNA[miRNA_VS_mRNA$miRNAstatus == 'Up', 'pair']
set2 <- miRNA_VS_mRNA[miRNA_VS_mRNA$RNAstatus == 'Down', 'pair']
set3 <- miRNA_VS_mRNA[miRNA_VS_mRNA$Spearman < -0.2, 'pair']

ggvenn(data = list("miRNA UP regulation" = set1,
                   "mRNA Down regulation" = set2, 
                   "Negative correlation" = set3),
       columns = c("miRNA UP regulation", "mRNA Down regulation", "Negative correlation"),
       fill_color = c("#7F7FFF", '#FFFF7F', '#7FFF7F'),
       set_name_color = c("#1818FF", '#CACA9A', '#0C860C'),
       set_name_size = 9,
       text_size = 7
)

# miRNA downregulation, mRNA upregulation, negative correlation
set1 <- miRNA_VS_mRNA[miRNA_VS_mRNA$miRNAstatus == 'Down', 'pair']
set2 <- miRNA_VS_mRNA[miRNA_VS_mRNA$RNAstatus == 'Up', 'pair']
set3 <- miRNA_VS_mRNA[miRNA_VS_mRNA$Spearman < -0.2, 'pair']

ggvenn(data = list("miRNA Down regulation" = set1,
                   "mRNA Up regulation" = set2, 
                   "Negative correlation" = set3),
       columns = c("miRNA Down regulation", "mRNA Up regulation", "Negative correlation"),
       fill_color = c("#7F7FFF", '#FFFF7F', '#7FFF7F'),
       set_name_color = c("#1818FF", '#CACA9A', '#0C860C'),
       set_name_size = 9,
       text_size = 7
)
