# Load required libraries
library(data.table)
library(ELMER)
library(ChAMP)
library(ggvenn)

# Load data for distal enhancers, probes, and DMP
# 'distal_DMP.RData' contains distal enhancer probes and their target information as defined by the ELMER package
load('~/distal_DMP.RData')

# Step 1: Processing distal probe data
# Retain differential probe-differential gene pairs
# 'RNA_DEG.RData' contains RNA expression matrix and differential expression results
load('~/RNA TMP/RNA_DEG.RData')
merge_RNA_DEG <- merge_RNA_DEG[merge_RNA_DEG$Status != 'None', ]
gene <- merge_RNA_DEG$Gene
a <- which(near_gene$Symbol %in% gene)
near_gene_driverDEG <- near_gene[a, ]
near_gene_driverDEG$MEstatue <- 0
for (i in 1:length(near_gene_driverDEG$ID)) {
  cat('Processing row', i, 'of', length(near_gene_driverDEG$ID), '\n')
  a <- which(rownames(DMP) == near_gene_driverDEG$ID[i])
  near_gene_driverDEG$MEstatue[i] <- DMP$staus[a]
}
a <- which(near_gene_driverDEG$MEstatue == 'None')
near_gene_driverDEG <- near_gene_driverDEG[-a, ]

# Step 2: Processing promoter probe data
promoter_probe_diff <- DMP[DMP$staus != 'None', ]
a <- which(promoter_probe_diff$gene %in% merge_RNA_DEG$Gene)
promoter_probe_diff <- promoter_probe_diff[a, ]
a <- which(promoter_probe_diff$feature == "3'UTR" | promoter_probe_diff$feature == 'Body')
promoter_probe_diff <- promoter_probe_diff[-a, ]

# Step 3: Processing body probe data
Body_probe_diff <- DMP[DMP$staus != 'None', ]
a <- which(Body_probe_diff$gene %in% merge_RNA_DEG$Gene)
Body_probe_diff <- Body_probe_diff[a, ]
a <- which(Body_probe_diff$feature == 'Body')
Body_probe_diff <- Body_probe_diff[a, ]

# Step 4: Integrate data
DMP_VS_DEG <- data.frame(
  Probe = near_gene_driverDEG$ID,
  Gene = near_gene_driverDEG$Symbol,
  Tpye = 'Distal Probe',
  Spearman = numeric(length(near_gene_driverDEG$ID)),
  P.value = numeric(length(near_gene_driverDEG$ID)),
  P.adjust = numeric(length(near_gene_driverDEG$ID)),
  DeltaB = numeric(length(near_gene_driverDEG$ID)),
  LogFC = numeric(length(near_gene_driverDEG$ID)),
  StatueME = near_gene_driverDEG$MEstatue,
  StatueRNA = character(length(near_gene_driverDEG$ID)),
  stringsAsFactors = FALSE
)

# Append promoter and body probe data
add_probe_data <- function(data, probe_data, type) {
  temp <- data.frame(
    Probe = rownames(probe_data),
    Gene = probe_data$gene,
    Tpye = type,
    Spearman = numeric(nrow(probe_data)),
    P.value = numeric(nrow(probe_data)),
    P.adjust = numeric(nrow(probe_data)),
    DeltaB = numeric(nrow(probe_data)),
    LogFC = numeric(nrow(probe_data)),
    StatueME = probe_data$staus,
    StatueRNA = character(nrow(probe_data)),
    stringsAsFactors = FALSE
  )
  rbind(data, temp)
}
DMP_VS_DEG <- add_probe_data(DMP_VS_DEG, promoter_probe_diff, 'Promoter Regions')
DMP_VS_DEG <- add_probe_data(DMP_VS_DEG, Body_probe_diff, 'Gene Body')

# Load processed data
# 'ALL_ME_DMP.RData' contains HM450k methylation beta matrix and differentially methylated probes defined by the ChAMP package
load('~/ALL_ME_DMP.RData')

# Align samples
co_sample <- intersect(colnames(merge_RNA_ComBat), substr(colnames(ME_Combat), 1, 15))

# RNA data alignment
RNA_T <- data.frame(matrix(0, nrow = nrow(merge_RNA_ComBat), ncol = length(co_sample)))
colnames(RNA_T) <- co_sample
rownames(RNA_T) <- rownames(merge_RNA_ComBat)
for (j in seq_along(co_sample)) {
  a <- which(colnames(merge_RNA_ComBat) == co_sample[j])
  RNA_T[, j] <- merge_RNA_ComBat[, a]
  cat('Aligning RNA data...', 'Processing column', j, '\n')
}

# Methylation data alignment
ME_T <- data.frame(matrix(0, nrow = nrow(ME_Combat), ncol = length(co_sample)))
colnames(ME_T) <- co_sample
rownames(ME_T) <- rownames(ME_Combat)
for (j in seq_along(co_sample)) {
  a <- which(substr(colnames(ME_Combat), 1, 15) == co_sample[j])
  ME_T[, j] <- ME_Combat[, a]
  cat('Aligning methylation data...', 'Processing column', j, '\n')
}

# Step 5: Correlation analysis
for (i in seq_along(DMP_VS_DEG$Probe)) {
  cat('Processing row', i, 'of', length(DMP_VS_DEG$Probe), '\n')
  a <- which(rownames(ME_T) == DMP_VS_DEG$Probe[i])
  b <- which(rownames(RNA_T) == DMP_VS_DEG$Gene[i])
  p <- cor.test(as.numeric(ME_T[a, ]), as.numeric(RNA_T[b, ]), method = "spearman", exact = FALSE)
  DMP_VS_DEG$Spearman[i] <- p$estimate
  DMP_VS_DEG$P.value[i] <- p$p.value
}
DMP_VS_DEG$P.adjust <- p.adjust(DMP_VS_DEG$P.value, method = 'BH')

# Add additional information
for (i in seq_along(DMP_VS_DEG$Probe)) {
  cat('Adding information for row', i, 'of', length(DMP_VS_DEG$Probe), '\n')
  a <- which(rownames(DMP) == DMP_VS_DEG$Probe[i])
  b <- which(merge_RNA_DEG$Gene == DMP_VS_DEG$Gene[i])
  DMP_VS_DEG$DeltaB[i] <- DMP$oppositedeltaBeta[a]
  DMP_VS_DEG$LogFC[i] <- merge_RNA_DEG$LogFC[b]
  DMP_VS_DEG$StatueME[i] <- DMP$staus[a]
  DMP_VS_DEG$StatueRNA[i] <- merge_RNA_DEG$Status[b]
}

# Filter significant results
significant_results <- which(DMP_VS_DEG$P.adjust <= 0.05)
DMP_VS_DEG <- DMP_VS_DEG[significant_results, ]
write.csv(DMP_VS_DEG, file = '~/DMP_VS_DEG.csv', row.names = FALSE)

# Save final results
save(DMP_VS_DEG, near_gene_driverDEG, file = '~/DMP_VS_DEG_result.RData')

# Visualization
DMP_VS_DEG$Paris <- paste(DMP_VS_DEG$Probe, DMP_VS_DEG$Gene, sep = '-')

# Define Venn diagram sets
venn_sets <- function(tpye, me_status, rna_status, correlation) {
  set1 <- DMP_VS_DEG[DMP_VS_DEG$Tpye == tpye & DMP_VS_DEG$StatueME == me_status, 'Paris']
  set2 <- DMP_VS_DEG[DMP_VS_DEG$Tpye == tpye & DMP_VS_DEG$StatueRNA == rna_status, 'Paris']
  set3 <- DMP_VS_DEG[DMP_VS_DEG$Tpye == tpye & DMP_VS_DEG$Spearman * correlation > 0.2, 'Paris']
  list(set1, set2, set3)
}

# Draw Venn diagrams
draw_venn <- function(sets, labels, colors, title) {
  ggvenn(
    data = list(labels[1] = sets[[1]], labels[2] = sets[[2]], labels[3] = sets[[3]]),
    columns = labels,
    fill_color = colors,
    set_name_color = c("#1818FF", '#CACA9A', '#0C860C'),
    set_name_size = 9,
    text_size = 7,
    main = title
  )
}

# Distal Probe - Hyper + Down + Negative correlation
sets <- venn_sets('Distal Probe', 'Hyper', 'Down', -1)
draw_venn(sets, c('Hyper', 'Down', 'Negative Correlation'), c("#1818FF", '#CACA9A', '#0C860C'), 'Distal Probe - Hyper + Down + Negative correlation')

# Distal Probe - Hypo + Up + Positive correlation
sets <- venn_sets('Distal Probe', 'Hypo', 'Up', 1)
draw_venn(sets, c('Hypo', 'Up', 'Positive Correlation'), c("#1818FF", '#CACA9A', '#0C860C'), 'Distal Probe - Hypo + Up + Positive correlation')

# Promoter Regions - Hyper + Down + Negative correlation
sets <- venn_sets('Promoter Regions', 'Hyper', 'Down', -1)
draw_venn(sets, c('Hyper', 'Down', 'Negative Correlation'), c("#1818FF", '#CACA9A', '#0C860C'), 'Promoter Regions - Hyper + Down + Negative correlation')

# Promoter Regions - Hypo + Up + Positive correlation
sets <- venn_sets('Promoter Regions', 'Hypo', 'Up', 1)
draw_venn(sets, c('Hypo', 'Up', 'Positive Correlation'), c("#1818FF", '#CACA9A', '#0C860C'), 'Promoter Regions - Hypo + Up + Positive correlation')

# Gene Body - Hyper + Down + Negative correlation
sets <- venn_sets('Gene Body', 'Hyper', 'Down', -1)
draw_venn(sets, c('Hyper', 'Down', 'Negative Correlation'), c("#1818FF", '#CACA9A', '#0C860C'), 'Gene Body - Hyper + Down + Negative correlation')

# Gene Body - Hypo + Up + Positive correlation
sets <- venn_sets('Gene Body', 'Hypo', 'Up', 1)
draw_venn(sets, c('Hypo', 'Up', 'Positive Correlation'), c("#1818FF", '#CACA9A', '#0C860C'), 'Gene Body - Hypo + Up + Positive correlation')
