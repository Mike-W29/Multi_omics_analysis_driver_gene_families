# Load necessary packages
library(dplyr)

# Define paths
path <- '~/Alignment_tsv/'

# Define paths
Domain <- read.csv('~/Hotsopt_domain.csv')
Hotsopt_MAF <- read.csv('~/Hotsopt_MAF.csv', header = TRUE)
pos_del_S_E <- read.csv('~/pos_del_S_E.csv', header = TRUE)

# Initialize Mutation_pos and Domain_length in Hotspot_MAF
Hotsopt_MAF$Mutation_pos <- 0
Hotsopt_MAF$Domain_length <- 0

# Read alignment files
alignment_files <- list.files(path = path)

# Check alignment files
check_alignment_files <- function(alignment_files, path) {
  for (file in alignment_files) {
    alignment_result <- read.table(paste(path, file, sep = '/'), header = TRUE)
    alignment_result$seq <- gsub('=', '', alignment_result$seq)  # 删除等号
    alignment_result$seq <- gsub(' ', '', alignment_result$seq)
    a <- nchar(alignment_result$seq)
    a <- a[!duplicated(a)]
    if (length(a) > 1) {
      message("Problem detected in ", file, ", needs checking", '\n')
    }
    b <- which(alignment_result$seq == '#NAME?')
    if (length(b) > 1) {
      message("Problem detected in ", file, ", needs checking", '\n')
    }
    if (length(colnames(alignment_result)) != 5) {
      message("Problem detected in ", file, ", needs checking", '\n')
    }
  }
}

check_alignment_files(alignment_files, path)

# Automatically record Domain mutation positions
record_mutation_positions <- function(alignment_files, path, Hotsopt_MAF, pos_del_S_E) {
  for (file in alignment_files) {
    cat(file, '\n')
    alignment_result <- read.table(paste(path, file, sep = '/'), header = TRUE)
    alignment_result$seq <- gsub('=', '', alignment_result$seq)
    alignment_result$seq <- gsub(' ', '', alignment_result$seq)
    alignment_result$length <- nchar(alignment_result$seq)
    alignment_result <- alignment_result[order(alignment_result$length, decreasing = TRUE),]
    
    pos <- as.data.frame(matrix(0, nrow = length(alignment_result$FASTA_NAME), ncol = alignment_result$length[1]))
    colnames(pos) <- 1:alignment_result$length[1]
    for (rows in 1:length(rownames(pos))) {
      for (cols in 1:alignment_result$length[1]) {
        pos[rows, cols] <- substr(alignment_result$seq[rows], cols, cols)
      }
    }
    
    pos_del <- as.data.frame(matrix(2, nrow = alignment_result$length[1], ncol = 3))
    colnames(pos_del) <- c('pos', 'numb_gap', 'gap_%')
    pos_del$pos <- 1:alignment_result$length[1]
    domain_index <- which(pos_del_S_E$Domain == gsub('.txt', '', file))
    
    if (pos_del_S_E$Start[domain_index] != 0) {
      b <- pos_del_S_E$Start[domain_index]
      colnames(pos)[1:b] <- 0
    }
    
    if (pos_del_S_E$End[domain_index] != alignment_result$length[1]) {
      b <- pos_del_S_E$End[domain_index]
      i <- alignment_result$length[1]
      colnames(pos)[b:i] <- 0
    }
    
    b <- pos_del_S_E$End[domain_index] - pos_del_S_E$Start[domain_index]
    i <- pos_del_S_E$Start[domain_index] + 1
    j <- pos_del_S_E$End[domain_index]
    colnames(pos)[i:j] <- 1:b
    
    pos_del$pos <- colnames(pos)
    for (rows in 1:length(rownames(pos_del))) {
      a <- which(pos[, rows] == '-')
      pos_del$numb_gap[rows] <- length(a)
      pos_del$`gap_%`[rows] <- pos_del$numb_gap[rows] / length(alignment_result$FASTA_NAME)
    }
    pos_del <- pos_del[pos_del$`gap_%` >= 0.75,]
    pos_del <- pos_del[pos_del$pos != 0,]
    alignment_result <- cbind(alignment_result, pos)
    
    for (i in 1:length(alignment_result$Uniport)) {
      sub_Domain_Pos <- alignment_result[i,]
      sub_Domain_Pos$length <- b
      a <- which(sub_Domain_Pos == '-')
      if (length(a) != 0) {
        sub_Domain_Pos <- sub_Domain_Pos[,-a]
      }
      a <- which(Hotsopt_MAF$Symbol %in% sub_Domain_Pos$Uniport)
      for (j in a) {
        if (Hotsopt_MAF$Position_aa[j] >= sub_Domain_Pos$envelope.start & Hotsopt_MAF$Position_aa[j] <= sub_Domain_Pos$envelope.end) {
          post <- Hotsopt_MAF$Position_aa[j] - sub_Domain_Pos$envelope.start
          Hotsopt_MAF$Mutation_pos[j] <- colnames(sub_Domain_Pos)[7 + post]
          Hotsopt_MAF$Domain_length[j] <- sub_Domain_Pos$length
          if (Hotsopt_MAF$Mutation_pos[j] %in% pos_del$pos) {
            Hotsopt_MAF$Mutation_pos[j] <- 'None'   # Genes with positions exceeding 75% gap rate are not recorded
          }
        }
      }
    }
  }
  
  Hotsopt_MAF$Mutation_pos[Hotsopt_MAF$Mutation_pos < 1] <- 0
  return(Hotsopt_MAF)
}

Hotsopt_MAF <- record_mutation_positions(alignment_files, path, Hotsopt_MAF, pos_del_S_E)

save(Hotsopt_MAF, file = '~/Hotsopt_MAF.RData')

# Calculate p-values
load('~/Hotsopt_MAF.RData')
Result_premutation_test <- read.csv('~/Result_premutation_test.csv', header = TRUE)
Hotsopt_MAF <- Hotsopt_MAF[Hotsopt_MAF$Mutation_pos != 0,]

# Correct the Domain_length for SHPRH gene
a <- which(Hotsopt_MAF$Symbol == 'SHPRH' & Hotsopt_MAF$pfam == 'PHD')
Hotsopt_MAF$Domain_length[a] <- 70

calculate_p_values <- function(Hotsopt_MAF, Result_premutation_test) {
  Hotspot <- data.frame()
  for (dom in 1:nrow(Result_premutation_test)) {
    Sub_Domain <- Result_premutation_test[dom,]
    Sub_Hotsopt_MAF <- Hotsopt_MAF[Hotsopt_MAF$pfam %in% Sub_Domain$Domain,]
    pos <- unique(as.numeric(Sub_Hotsopt_MAF$Mutation_pos[Sub_Hotsopt_MAF$Mutation_pos != 'None']))
    
    if (length(pos) == 0) next
    
    Sub_Hotspot <- data.frame(Domain = rep(Sub_Domain$Domain, length(pos)), Position = pos, Mut = 0, P.value = 0, P.adjust = 0)
    for (i in 1:length(Sub_Hotspot$Position)) {
      Sub_Hotspot$Mut[i] <- sum(Sub_Hotsopt_MAF$Mutation_pos == Sub_Hotspot$Position[i])
    }
    
    L <- unique(Sub_Hotsopt_MAF$Domain_length)
    px <- 1 / L
    n <- Sub_Domain$Mutation
    
    for (i in 1:length(Sub_Hotspot$Position)) {
      p.value <- 0 
      k <- Sub_Hotspot$Mut[i]
      for (m in k:n) {
        p.value <- p.value + choose(n, m) * px^m * (1 - px)^(n - m)
      }
      if (length(p.value) > 1 || p.value > 1) {
        stop('p.value error')
      }
      Sub_Hotspot$P.value[i] <- p.value
    }
    
    Sub_Hotspot$P.adjust <- p.adjust(Sub_Hotspot$P.value, method = 'BH')
    Hotspot <- rbind(Hotspot, Sub_Hotspot)
  }
  
  return(Hotspot)
}

Hotspot <- calculate_p_values(Hotsopt_MAF, Result_premutation_test)

# After calculating Hotspot data frame, add the following steps:

# Delete positions with P.adjust less than 0.05
Hotspot <- Hotspot[Hotspot$P.adjust < 0.05,]

# Initialize additional columns
Hotspot$S <- 0
Hotspot$Top_Gene1 <- ''
Hotspot$Top_Gene2 <- ''
Hotspot$Top_Gene3 <- ''

# Process Top Genes and entropy value calculation for each Domain
for (i in 1:nrow(Hotspot)) {
  Sub_MAF <- Hotsopt_MAF[Hotsopt_MAF$pfam == Hotspot$Domain[i],]
  Sub_MAF <- Sub_MAF[Sub_MAF$Mutation_pos == Hotspot$Position[i],]
  genes <- unique(Sub_MAF$Symbol)
  
  # Create Mu_Gene data frame
  Mu_Gene <- data.frame(Gene = genes, numb = 0)
  for (j in 1:nrow(Mu_Gene)) {
    a <- which(Sub_MAF$Symbol == Mu_Gene$Gene[j])
    Mu_Gene$numb[j] <- length(a)
  }
  
  # Sort by gene occurrence
  Mu_Gene <- Mu_Gene[order(Mu_Gene$numb, decreasing = TRUE),]
  
  # Fill Top_Gene1, Top_Gene2, Top_Gene3
  if (nrow(Mu_Gene) >= 3) {
    for (k in 1:3) {
      Hotspot[[paste0("Top_Gene", k)]][i] <- paste(Mu_Gene$numb[k], Mu_Gene$Gene[k], '(', sep = '')
      a <- which(Sub_MAF$Symbol == Mu_Gene$Gene[k])
      a <- a[1]
      a <- paste(Sub_MAF$Ref_aa[a], Sub_MAF$Position_aa[a], sep = '')
      a <- gsub(' ', '', a)
      Hotspot[[paste0("Top_Gene", k)]][i] <- paste(Hotspot[[paste0("Top_Gene", k)]][i], a, ')', sep = '')
    }
  } else if (nrow(Mu_Gene) == 2) {
    for (k in 1:2) {
      Hotspot[[paste0("Top_Gene", k)]][i] <- paste(Mu_Gene$numb[k], Mu_Gene$Gene[k], '(', sep = '')
      a <- which(Sub_MAF$Symbol == Mu_Gene$Gene[k])
      a <- a[1]
      a <- paste(Sub_MAF$Ref_aa[a], Sub_MAF$Position_aa[a], sep = '')
      a <- gsub(' ', '', a)
      Hotspot[[paste0("Top_Gene", k)]][i] <- paste(Hotspot[[paste0("Top_Gene", k)]][i], a, ')', sep = '')
    }
  } else if (nrow(Mu_Gene) == 1) {
    Hotspot$Top_Gene1[i] <- paste(Mu_Gene$numb[1], Mu_Gene$Gene[1], '(', sep = '')
    a <- which(Sub_MAF$Symbol == Mu_Gene$Gene[1])
    a <- a[1]
    a <- paste(Sub_MAF$Ref_aa[a], Sub_MAF$Position_aa[a], sep = '')
    a <- gsub(' ', '', a)
    Hotspot$Top_Gene1[i] <- paste(Hotspot$Top_Gene1[i], a, ')', sep = '')
  }
  
  # Calculate entropy value
  s <- 0
  L <- Hotspot$Mut[i]
  for (j in 1:nrow(Mu_Gene)) {
    s <- s - Mu_Gene$numb[j] / L * log(Mu_Gene$numb[j] / L)
  }
  if (nrow(Mu_Gene) > 1) {
    s <- s / log(nrow(Mu_Gene))
  } else {
    s <- 0
  }
  Hotspot$S[i] <- s
  
  cat('Progress: ', i / nrow(Hotspot) * 100, '%', '\n')
}

# Delete rows where Mu is greater than 1
Hotspot <- Hotspot[Hotspot$Mu > 1,]

# Output results
write.csv(Hotspot, file = '~/Domain_Hotspot.csv')




