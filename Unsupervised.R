## BIC formula:
## bic <- n_parameters * log(nrow(counts)) - 2 * totallogliks

### n_parameters: the dimension of the profiles matrix, i.e. # clusters * # genes
### profiles: [960*13]

###  nrow(counts): the number of cells, i.e. “n”

### logliks: a vector giving each cell’s log likelihood under the clustering solution. 
### For each cell, this is the sum of the loglikelihoods of its individual count values.
### totallogliks: max of each row from logliks, sum across all rows

calculate_bic <- function(insitutype_output) {
        
        # Get total log likelihood----
        # Extract log-likelihood matrix
        log_liks <- insitutype_output$logliks
        # Initialize 
        log_liks_total <- 0
        
        # Loop over log_liks
        for(ii in 1:nrow(log_liks)) {
                log_liks_total <- log_liks_total + max(log_liks[ii, ])
        }
        
        # Get k = number of parameters in the model----
        # Extract profile matrix
        profiles <- insitutype_output$profiles
        # Number of genes
        num_genes <- dim(profiles)[1]
        # Number of clusters
        num_clusters <- dim(profiles)[2]
        # Number of parameters
        num_param <- num_genes * num_clusters
        
        # Get n = number of data points----
        num_data <- nrow(insitutype_output$logliks)
        
        # Calculate BIC
        bic_final <- num_param * log(num_data) - 2 * log_liks_total
        
        return(bic_final)
}

## ARI
## Sample code
#### unsupervised clustering -------------------------
# without cohorts 
unsup.nocohorts.13.10start10k = insitutype(x = counts,
                                           neg = annot$negmean,
                                           n_clusts = 13,
                                           n_starts = 10,
                                           n_phase1 = 10000)

# if cohorts
ifdata = as.matrix(annot[, paste0(c("Blue", "Red", "Yellow", "Green"), "Int")])
ifcohort = fastCohorting(mat = ifdata)

unsup.ifcohorts.13.10start10k = insitutype(x = counts,
                                           neg = annot$negmean,
                                           n_clusts = 13,
                                           n_starts = 10,
                                           n_phase1 = 10000,
                                           cohort = ifcohort)
## ARI from no cohort
mclust::adjustedRandIndex(unsup.nocohorts.13.10start10k$clust, annot$cell_line)
## ARI from IF cohort
mclust::adjustedRandIndex(unsup.ifcohorts.13.10start10k$clust, annot$cell_line)


### change different n_starts value
nstart = c(8,9,10,11,12)
ifco_output = vector("list", 5)
noco_output = vector("list", 5)

for (i in 1:length(nstart)){
        
        # without cohorts 
        noco_output[[i]] = insitutype(x = counts, 
                                      neg = annot$negmean,
                                      n_clusts = 13,
                                      n_starts = nstart[i],
                                      n_phase1 = 10000)
        
        # if cohorts
        ifdata = as.matrix(annot[, paste0(c("Blue", "Red", "Yellow", "Green"), "Int")])
        ifcohort = fastCohorting(mat = ifdata)
        
        ifco_output[[i]] = insitutype(x = counts,
                                      neg = annot$negmean, 
                                      n_clusts = 13,
                                      n_starts = nstart[i],
                                      n_phase1 = 10000,
                                      cohort = ifcohort)
       }

