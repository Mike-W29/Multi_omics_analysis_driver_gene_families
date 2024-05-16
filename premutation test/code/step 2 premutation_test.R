

premutation_pvalue <- function(premutation_Domain, observe_domain, observe_sum, Sum_domain) {
  Domain <- rownames(premutation_Domain)
  Result_premutation_test <- data.frame(
    Domain = Domain,
    Gene = numeric(length(Domain)),
    p.value = numeric(length(Domain)),
    Mutation = numeric(length(Domain)),
    Top_5_Gene = character(length(Domain)),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(Domain)) {
    domain_name <- Domain[i]
    
    # Extract mutation count from observe_sum
    mutation_index <- which(rownames(observe_sum) == domain_name)
    Result_premutation_test$Mutation[i] <- observe_sum$SumMutation[mutation_index]
    
    # Extract gene count and top 5 genes from observe_domain
    domain_col <- which(colnames(observe_domain) == domain_name)
    sub_observe_domain <- observe_domain[order(observe_domain[, domain_col], decreasing = TRUE), ]
    sub_observe_domain <- sub_observe_domain[sub_observe_domain[, domain_col] != 0, ]
    
    Result_premutation_test$Gene[i] <- nrow(sub_observe_domain)
    
    if (nrow(sub_observe_domain) >= 5) {
      top_5_genes <- paste(sub_observe_domain$Symbol[1:5], sub_observe_domain[1:5, domain_col], sep = '(')
    } else {
      top_5_genes <- paste(sub_observe_domain$Symbol, sub_observe_domain[, domain_col], sep = '(')
    }
    
    Result_premutation_test$Top_5_Gene[i] <- paste(top_5_genes, collapse = ');')
    
    # Calculate p-value
    premutation_col <- which(rownames(premutation_Domain) == domain_name)
    mub <- sum(as.numeric(premutation_Domain[premutation_col, ]) >= observe_sum$SumMutation[mutation_index])
    Result_premutation_test$p.value[i] <- (mub + 1) / (10000 + 1)
  }
  
  return(Result_premutation_test)
}


Entropy <- function(Result_premutation_test, observe_domain) {
  Result_premutation_test$S <- 0
  
  for (i in seq_len(nrow(Result_premutation_test))) {
    domain_name <- Result_premutation_test$Domain[i]
    domain_col <- which(colnames(observe_domain) == domain_name)
    sub_observe_domain <- observe_domain[order(observe_domain[, domain_col], decreasing = TRUE), ]
    sub_observe_domain$sum <- rowSums(sub_observe_domain[, 2:length(sub_observe_domain)])
    sub_observe_domain <- sub_observe_domain[sub_observe_domain$sum != 0 & sub_observe_domain[, domain_col] != 0, ]
    
    n <- nrow(sub_observe_domain)
    ob_mu <- sum(sub_observe_domain[, domain_col])
    
    h <- sum(-sub_observe_domain[, domain_col] / ob_mu * log(sub_observe_domain[, domain_col] / ob_mu))
    s <- h / log(n)
    
    Result_premutation_test$S[i] <- s
  }
  
  return(Result_premutation_test)
}


# Load necessary data if not already loaded

load('~/result_premutation.RData')

# Perform mutation analysis and statistical calculations
Result_premutation_test <- premutation_pvalue(premutation_Domain, 
                                              observe_domain, observe_sum, Sum_domain)

# Filter domains with adjusted p-value < 0.05
Result_premutation_test$p.adjust <- p.adjust(Result_premutation_test$p.value, method = 'BH')
Result_premutation_test <- Result_premutation_test[Result_premutation_test$p.adjust < 0.05, ]

# Calculate entropy-based statistics
Result_premutation_test <- Entropy(Result_premutation_test, observe_domain)

# Select and reorder columns
Result_premutation_test <- Result_premutation_test[, c('Domain', 'Gene', 'p.value', 'p.adjust', 'S', 'Mutation', 'Top_5_Gene')]
Result_premutation_test$S[is.na(Result_premutation_test$S)] <- 0

# Save results to CSV file
write.csv(Result_premutation_test, file = '~/Result_premutation_test_NN.csv')

# Save objects to RData file
save(Merge_Mutation_Driverfamily, Result_premutation_test, family_domain, family_length,
     file = '~/Result_premutation_testN.RData')

