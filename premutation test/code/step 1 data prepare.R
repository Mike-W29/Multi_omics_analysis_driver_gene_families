
#####

merge_data <- function(path = NULL, method = c('row', 'col')) {
  message('Make sure row names or column names are consistent.')
  
  if (is.null(path)) {
    stop('Please provide a file path!')
  }
  
  if (is.null(method) || !method %in% c('row', 'col')) {
    stop("Method should be 'row' or 'col'!")
  }
  
  library(data.table)
  
  filenames <- list.files(path = path, full.names = TRUE)
  mydata <- fread(filenames[1], header = TRUE)
  
  for (i in 2:length(filenames)) {
    temp_data <- fread(filenames[i], header = TRUE)
    
    if (method == 'row') {
      mydata <- rbind(mydata, temp_data)
    } else if (method == 'col') {
      mydata <- cbind(mydata, temp_data)
    }
  }
  
  print('Data merging completed.')
  return(mydata)
}

where_mutation_is <- function(MAF = NULL, pfamDomain = NULL, ID = NULL) {
  message('The first two columns of ID should be Uniport ID and Symbol.')
  
  if (is.null(MAF) || is.null(pfamDomain) || is.null(ID)) {
    stop('Please provide MAF, pfamDomain, and ID files!')
  }
  
  pfamDomain$Symbol <- ID[match(pfamDomain$UniportID, ID[,1]), 2]
  Domains <- unique(pfamDomain$`hmm name`)
  
  MAF <- MAF[MAF$Symbol %in% pfamDomain$Symbol, ]
  gene <- unique(MAF$Symbol)
  
  Domain_mutation <- data.frame(Symbol = gene, 
                                matrix(0, nrow = length(gene), ncol = length(Domains) + 1))
  colnames(Domain_mutation)[-1] <- c(Domains, "outDomain")
  
  cat('Starting mutation site analysis.\n')
  
  for (i in seq_along(gene)) {
    cat('Analyzing mutations... Progress:', i/length(gene) * 100, '%\n')
    
    Sub_MAF <- MAF[MAF$Symbol == gene[i], ]
    Sub_Domain <- pfamDomain[pfamDomain$Symbol == gene[i], ]
    
    for (mu in seq_len(nrow(Sub_MAF))) {
      domain_matched <- FALSE
      
      for (dom in seq_len(nrow(Sub_Domain))) {
        if (Sub_MAF$Position_aa[mu] >= Sub_Domain$`envelope start`[dom] &&
            Sub_MAF$Position_aa[mu] <= Sub_Domain$`envelope end`[dom]) {
          Domain_mutation[i, Sub_Domain$`hmm name`[dom] ] <- 
            Domain_mutation[i, Sub_Domain$`hmm name`[dom]] + 1
          domain_matched <- TRUE
          break
        }
      }
      
      if (!domain_matched) {
        Domain_mutation[i, "outDomain"] <- Domain_mutation[i, "outDomain"] + 1
      }
    }
  }
  
  cat('Mutation site analysis completed.\n')
  return(Domain_mutation)
}


premutation_test <- function(observe_domain = NULL, pfamDomain = NULL, n = 10000, ID = NULL, seed = 100) {
  message('The first two columns of ID should be Uniport ID and Symbol.')
  
  if (is.null(pfamDomain) || is.null(ID)) {
    stop('Please provide pfamDomain and ID files!')
  }
  
  pfamDomain$Symbol <- ID[match(pfamDomain$UniportID, ID[,1]), 2]
  Domain <- colnames(observe_domain)[!colnames(observe_domain) %in% c('Symbol', 'outDomain')]
  
  cat('Detected', length(Domain), 'domains.\n')
  
  premutation_Domain <- data.frame(matrix(0, nrow = length(Domain), ncol = n))
  rownames(premutation_Domain) <- Domain
  
  cat('Starting mutation simulation...\n')
  time_start <- Sys.time()
  
  Genes <- unique(pfamDomain$Symbol)
  
  for (dom in seq_along(Domain)) {
    cat('Analyzing records for', Domain[dom], '... Progress:', dom/length(Domain) * 100, '%\n')
    
    Sub_Domain <- pfamDomain[pfamDomain$`hmm name` == Domain[dom], ]
    subgenes <- unique(Sub_Domain$Symbol)
    
    for (gene in subgenes) {
      set.seed(seed)
      
      for (i in 1:n) {
        a <- which(observe_domain$Symbol == gene)
        
        if (length(a) == 0) {
          break
        }
        
        inDomain <- observe_domain[a, Domain[dom]]
        outDomain <- observe_domain$outDomain[a]
        
        if (inDomain + outDomain == 0) {
          break
        }
        
        a <- which(Sub_Domain$Symbol == gene)
        a <- a[1]
        
        premutation <- sample(1:Sub_Domain$ProteinLength[a], inDomain + outDomain, replace = TRUE)
        
        for (mu in seq_along(premutation)) {
          for (subdom in seq_len(nrow(Sub_Domain))) {
            if (premutation[mu] >= Sub_Domain$`envelope start`[subdom] &&
                premutation[mu] <= Sub_Domain$`envelope end`[subdom]) {
              premutation_Domain[Domain[dom], i] <- premutation_Domain[Domain[dom], i] + 1
              break
            }
          }
        }
      }
    }
    
    cat('Completed records for', Domain[dom], '.\n')
    cat(Domain[dom], 'mutation simulation completed! Time taken:', Sys.time() - time_start, '\n')
  }
  
  cat('Mutation simulation completed! Total time taken:', Sys.time() - time_start, '\n')
  return(premutation_Domain)
}

# Merge domain and length files
family_domain <- merge_data(path = '~/test data/pfam_annotation/', method = 'row')
family_length <- merge_data(path = '~/length_files/', method = 'row')

# Load MAF file
Merge_Mutation_Driverfamily <- read.csv('~/Mutation_file.csv', header = TRUE)

# Load family file
gene_family_Driver <- read.csv('~/test data/Driver_familyID.csv', header = TRUE)

# Check the number of genes and remove those without domain information
table(duplicated(family_domain$UniportID))

# Filter out non-point mutations
# a <- which(nchar(Merge_Mutation_Driverfamily$Reference_Allele) > 1)
# Merge_Mutation_Driverfamily <- Merge_Mutation_Driverfamily[-a,]

# Annotate Length file
gene_family_Driver <- gene_family_Driver[!duplicated(gene_family_Driver$uniport.ID), ]
family_length$Symbol <- ''
for (i in 1:length(family_length$From)) {
  a <- which(gene_family_Driver$uniport.ID == family_length$From[i])
  family_length$Symbol[i] <- gene_family_Driver$smybol[a]
}
family_length <- dplyr::select(family_length, c('From', 'Symbol', 'Length'))

# Filter mutations longer than total amino acid length
Merge_Mutation_Driverfamily <- Merge_Mutation_Driverfamily[Merge_Mutation_Driverfamily$Symbol %in% family_length$Symbol, ]
Merge_Mutation_Driverfamily$Length <- 0
family_length <- family_length[!duplicated(family_length$Symbol), ]
for (i in 1:length(Merge_Mutation_Driverfamily$Symbol)) {
  a <- which(family_length$Symbol == Merge_Mutation_Driverfamily$Symbol[i])
  Merge_Mutation_Driverfamily$Length[i] <- family_length$Length[a]
}
# a <- which(Merge_Mutation_Driverfamily$Position_aa > Merge_Mutation_Driverfamily$Length)
# Merge_Mutation_Driverfamily <- Merge_Mutation_Driverfamily[-a,]

# Remove extra UniportDomain entries
family_domain <- family_domain[family_domain$UniportID %in% family_length$From, ]

# Calculate Domain / Protein Length
family_domain$ProteinLength <- 0
for (i in 1:length(family_domain$UniportID)) {
  a <- which(family_length$From == family_domain$UniportID[i])
  family_domain$ProteinLength[i] <- family_length$Length[a]
}

# Calculate percentage of protein length covered by Domain
family_domain$percentage <- (family_domain$`envelope end` - family_domain$`envelope start` + 1) / family_domain$ProteinLength
family_domain <- family_domain[family_domain$percentage < 0.75, ]

# Check domain information
table(duplicated(family_domain$`hmm name`))

# Summarize domain information
Sum_domain <- data.frame(Domain = unique(family_domain$`hmm name`), Sum = 0)
for (i in 1:nrow(Sum_domain)) {
  a <- which(family_domain$`hmm name` == Sum_domain$Domain[i])
  a <- family_domain$UniportID[a]
  a <- a[!duplicated(a)]
  Sum_domain$Sum[i] <- length(a)
}

# Remove duplicates
table(duplicated(family_domain))
family_domain <- family_domain[!duplicated(family_domain), ]
table(duplicated(family_domain))

# Analyze family mutation positions
observe_domain <- where_mutation_is(MAF = Merge_Mutation_Driverfamily,
                                    pfamDomain = family_domain,
                                    ID = gene_family_Driver[, 2:3])

# Summarize mutations
observe_sum <- as.data.frame(apply(observe_domain[, 2:ncol(observe_domain)], 2, sum))
colnames(observe_sum) <- 'SumMutation'
observe_sum$Domain <- rownames(observe_sum)

# Filter out unwanted domains (less than 5 mutations or appearing only once)
observe_sum <- observe_sum[observe_sum$SumMutation >= 5, ]
b <- which(colnames(observe_domain) %in% observe_sum$Domain)
a <- observe_domain
observe_domain <- observe_domain[, b]

# Keep only domains with at least two occurrences
Sum_domain <- Sum_domain[Sum_domain$Sum > 1, ]
a <- dplyr::select(a, c('Symbol', 'outDomain'))
b <- which(colnames(observe_domain) %in% Sum_domain$Domain)
observe_domain <- observe_domain[, b]
observe_domain <- cbind(a, observe_domain)

# Recalculate sums
observe_sum <- as.data.frame(apply(observe_domain[, 2:ncol(observe_domain)], 2, sum))
colnames(observe_sum) <- 'SumMutation'

# Perform permutation test
premutation_Domain <- premutation_test(observe_domain = observe_domain,
                                       pfamDomain = family_domain,
                                       ID = gene_family_Driver[, 2:3],
                                       n = 10000)

# Save results



save.image(file = '~/result_premutation.RData')








