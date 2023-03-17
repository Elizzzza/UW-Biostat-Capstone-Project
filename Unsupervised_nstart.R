# install.packages(devtools)
# devtools::install_github('CorradoLanera/CrossClustering', ref = 'develop')
library(tibble)
library(readr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(InSituType)
library(clue)
library(pheatmap)
library(SingleCellExperiment)
library(scater)
library(mclust)
library(scran)
library(igraph)
library(SC3)
library(Seurat)
library(CHETAH)
library(CrossClustering)
library(tidyverse)
rm(list = ls())
pal = tableau_color_pal()(10) #color palette

#### load CPA data: ---------------
load("Cell Pellet Array annotation and raw counts.RData")
load("CPA16_RNAseq.RData")
counts = t(as.matrix(raw))
rm(raw)
badprobes = read.csv("genes with efficiency 8-fold below average in old CPA panel.csv")[,2]
counts = counts[, !is.element(colnames(counts), badprobes)]
boxplot(log(annot$raw_totalCount)~annot$cell_line,las=2)
abline(h = 150, col = 2)

# remove the cell lines with failed FOVs:
failed.cell.lines = names(which(by(annot$raw_totalCount, annot$cell_line, median) < 150))
remove = is.element(annot$cell_line, failed.cell.lines)
annot = annot[!remove, ]
counts = counts[!remove, ]

#### unsupervised clustering -------------------------

##### Full Sample #####
# if cohorts
ifdata = as.matrix(annot[, paste0(c("Blue", "Red", "Yellow", "Green"), "Int")])
ifcohort = fastCohorting(mat = ifdata)

# change different n_starts value
# change nstart for other tuning parameter
nstart = c(1:20)
ifco_output = vector("list", length(nstart))
noco_output = vector("list", length(nstart))

ifco_timing = data.frame(n_starts = 1:length(nstart), 
                         ifco_time = length(nstart))
noco_timing = data.frame(n_starts = 1:length(nstart), 
                    noco_time = length(nstart))

for (i in 1:length(nstart)){
        ## No cohort
        # Record start time
        noco_start_time <- Sys.time()
        
        # without cohorts 
        noco_output[[i]] = insitutype(x = counts, 
                                      neg = annot$negmean,
                                      n_clusts = 13,
                                      n_starts = nstart[i],
                                      n_phase1 = 10000)
        
        # Record end time and elapsed time
        noco_end_time <- Sys.time()
        noco_elapsed_time <- as.numeric(noco_end_time - noco_start_time)
        
        # Store the elapsed time in the timing data frame
        noco_timing[i, "noco_time"] <- noco_elapsed_time
        
        ## IF cohort
        # Record start time
        ifco_start_time <- Sys.time()
        
        ifco_output[[i]] = insitutype(x = counts,
                                      neg = annot$negmean, 
                                      n_clusts = 13,
                                      n_starts = nstart[i],
                                      n_phase1 = 10000,
                                      cohort = ifcohort)
        # Record end time and elapsed time
        ifco_end_time <- Sys.time()
        ifco_elapsed_time <- as.numeric(ifco_end_time - ifco_start_time)
        
        # Store the elapsed time in the timing data frame
        ifco_timing[i, "ifco_time"] <- ifco_elapsed_time
        
}

# function to calculate ARI CI
ari_ci <- function(output, truth){
        temp_ari = ari(table(output$clust, truth$cell_line))
        ariCI = list(lower = attr(temp_ari$ari, "ci")["lower"],
                     upper = attr(temp_ari$ari, "ci")["upper"])
        return(ariCI)
}

# function to calculate BIC
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

# get ARI and CI
# change n_starts column in the ari_df for different parameter
get_ari <- function(noco_output, ifco_output, truth) {
        # Generate ARI for no cohort
        noco_ari <- sapply(seq_along(noco_output), function(i) {
                adjustedRandIndex(noco_output[[i]]$clust, truth$cell_line)
        })
        
        # Generate ARI for IF cohort
        ifco_ari <- sapply(seq_along(ifco_output), function(i) {
                adjustedRandIndex(ifco_output[[i]]$clust, truth$cell_line)
        })
        
        # Create data frame with ARI
        ari_df <- data.frame(
                n_starts = 1:length(noco_output),
                no_cohort_ari = noco_ari,
                if_cohort_ari = ifco_ari
        )
        
        ## Get CI
        # Initialize an empty data frame to store the results
        noco_ari_ci_df <- data.frame(lower = numeric(length(noco_output)),
                                     upper = numeric(length(noco_output)))
        ifco_ari_ci_df <- data.frame(lower = numeric(length(ifco_output)),
                                     upper = numeric(length(ifco_output)))
        
        for (i in seq_along(noco_output)) {
                noco.ariCI <- ari_ci(noco_output[[i]], truth)
                noco_ari_ci_df[i, ] <- unlist(noco.ariCI)
        }
        
        for (i in seq_along(ifco_output)) {
                ifco.ariCI <- ari_ci(ifco_output[[i]], truth)
                ifco_ari_ci_df[i, ] <- unlist(ifco.ariCI)
        }
        colnames(noco_ari_ci_df) <- c("noco_ari_ci_lower", "noco_ari_ci_upper")
        colnames(ifco_ari_ci_df) <- c("ifco_ari_ci_lower", "ifco_ari_ci_upper")
        ari_ci_df = cbind(noco_ari_ci_df, ifco_ari_ci_df)
        
        ari_df <- cbind(ari_df, ari_ci_df)
        
        return(ari_df)
}

# get BIC
# change n_starts column in the bic_df for different parameter
get_bic <- function(noco_output, ifco_output) {
        # Get BIC for no cohort
        noco_bic <- sapply(seq_along(noco_output), function(k) {
                calculate_bic(noco_output[[k]])
        })
        
        # Generate BIC for IF cohort
        ifco_bic <- sapply(seq_along(ifco_output), function(k) {
                calculate_bic(ifco_output[[k]])
        })
        
        # Create data frame with index, results from list 1, and results from list 2
        bic_df <- data.frame(
                n_starts = 1:length(noco_output),
                no_cohort_bic = noco_bic,
                if_cohort_bic = ifco_bic
        )
        
        return(bic_df)
}


# Print the timing data frame
cat("Timing data frame:\n")
print(noco_timing)
print(ifco_timing)

# Performance data frame
ari.df <- get_ari(noco_output, ifco_output, annot)
bic.df <- get_bic(noco_output, ifco_output)

# put all data frames into list
df_list <- list(ari.df, bic.df, noco_timing, ifco_timing)
# merge all data frames in list
df <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
df <- df %>% mutate_if(is.numeric, round, digits = 5)
write.csv(df, "unsup_full_df.csv", row.names = FALSE)

##### Subsample genes #####
# Half genes
subsample_half = sample.int(ncol(counts), ncol(counts)/2)
# if cohorts
ifdata = as.matrix(annot[, paste0(c("Blue", "Red", "Yellow", "Green"), "Int")])
ifcohort = fastCohorting(mat = ifdata)

# change different n_starts value
nstart = c(1:20)
half_ifco_output = vector("list", length(nstart))
half_noco_output = vector("list", length(nstart))

half_ifco_timing = data.frame(n_starts = 1:length(nstart), 
                         half_ifco_time = length(nstart))
half_noco_timing = data.frame(n_starts = 1:length(nstart), 
                         half_noco_time = length(nstart))

for (i in 1:length(nstart)){
        ## No cohort
        # Record start time
        noco_start_time <- Sys.time()
        
        # without cohorts 
        half_noco_output[[i]] = insitutype(x = counts[,subsample_half], 
                                           neg = annot$negmean,
                                           n_clusts = 13,
                                           n_starts = nstart[i],
                                           n_phase1 = 10000)
        
        # Record end time and elapsed time
        noco_end_time <- Sys.time()
        noco_elapsed_time <- as.numeric(noco_end_time - noco_start_time)
        
        # Store the elapsed time in the timing data frame
        half_noco_timing[i, "half_noco_time"] <- noco_elapsed_time
        
        ## IF cohort
        # Record start time
        ifco_start_time <- Sys.time()
        
        
        half_ifco_output[[i]] = insitutype(x = counts[,subsample_half],
                                           neg = annot$negmean, 
                                           n_clusts = 13,
                                           n_starts = nstart[i],
                                           n_phase1 = 10000,
                                           cohort = ifcohort)
        # Record end time and elapsed time
        ifco_end_time <- Sys.time()
        ifco_elapsed_time <- as.numeric(ifco_end_time - ifco_start_time)
        
        # Store the elapsed time in the timing data frame
        half_ifco_timing[i, "half_ifco_time"] <- ifco_elapsed_time
        
}

# Print the timing data frame
cat("Timing data frame:\n")
print(half_noco_timing)
print(half_ifco_timing)

# Performance dataframe
half.ari.df <- get_ari(half_noco_output, half_ifco_output, annot)
half.bic.df <- get_bic(half_noco_output, half_ifco_output)

# put all data frames into list
half_df_list <- list(half.ari.df, half.bic.df, half_noco_timing, half_ifco_timing)
# merge all data frames in list
half.df <- Reduce(function(x, y) merge(x, y, all=TRUE), half_df_list)


##### 96 HVGs #####
## Preprocess
sce_hvg = SingleCellExperiment(assays = list(counts = t(counts)), #make SCE
                               colData = annot)
sce_hvg = logNormCounts(sce_hvg) #normalize
top = getTopHVGs(sce_hvg, n = 96)
sce_hvg = sce_hvg[top, ]
sce_hvg = scater::runPCA(sce_hvg, subset_row = top)

counts_hvg = counts[,top]
hvg_filter = rowSums(counts_hvg) >= 5
annot_hvg = annot[hvg_filter,]
counts_hvg = counts_hvg[hvg_filter,]
sce_hvg = sce_hvg[,hvg_filter]
pcs_hvg = reducedDim(sce_hvg)

# if cohorts
ifdata_hvg = as.matrix(annot_hvg[, paste0(c("Blue", "Red", "Yellow", "Green"), "Int")])
ifcohort_hvg = fastCohorting(mat = ifdata_hvg)


# change different n_starts value
nstart = c(1:20)
hvg_ifco_output = vector("list", length(nstart))
hvg_noco_output = vector("list", length(nstart))

hvg_ifco_timing = data.frame(n_starts = 1:length(nstart), 
                              hvg_ifco_time = length(nstart))
hvg_noco_timing = data.frame(n_starts = 1:length(nstart), 
                              hvg_noco_time = length(nstart))

for (i in 1:length(nstart)){
        ## No cohort
        # Record start time
        noco_start_time <- Sys.time()
        
        # without cohorts 
        hvg_noco_output[[i]] = insitutype(x = counts_hvg,
                                          neg = annot_hvg$negmean,
                                          n_clusts = 13,
                                          n_starts = nstart[i],
                                          n_phase1 = 10000)
        
        # Record end time and elapsed time
        noco_end_time <- Sys.time()
        noco_elapsed_time <- as.numeric(noco_end_time - noco_start_time)
        
        # Store the elapsed time in the timing data frame
        hvg_noco_timing[i, "hvg_noco_time"] <- noco_elapsed_time
        
        ## IF cohorr
        # Record start time
        ifco_start_time <- Sys.time()
        
        
        hvg_ifco_output[[i]] = insitutype(x = counts_hvg,
                                          neg = annot_hvg$negmean,
                                          n_clusts = 13,
                                          n_starts = nstart[i],
                                          n_phase1 = 10000,
                                          cohort = ifcohort_hvg)
        # Record end time and elapsed time
        ifco_end_time <- Sys.time()
        ifco_elapsed_time <- as.numeric(ifco_end_time - ifco_start_time)
        
        # Store the elapsed time in the timing data frame
        hvg_ifco_timing[i, "hvg_ifco_time"] <- ifco_elapsed_time
        
}

# Print the timing data frame
cat("Timing data frame:\n")
print(hvg_noco_timing)
print(hvg_ifco_timing)

# Performance dataframe
hvg.ari.df <- get_ari(hvg_noco_output, hvg_ifco_output, annot_hvg)
hvg.bic.df <- get_bic(hvg_noco_output, hvg_ifco_output)

# put all data frames into list
hvg_df_list <- list(hvg.ari.df, hvg.bic.df, hvg_noco_timing, hvg_ifco_timing)
# merge all data frames in list
hvg.df <- Reduce(function(x, y) merge(x, y, all=TRUE), hvg_df_list)

write.csv(full.df, file = "upsup_full_df.csv", row.names = FALSE)
write.csv(half.df, file = "unsup_half_df.csv", row.names = FALSE)
write.csv(hvg.df, file = "unsup_hvg_df.csv", row.names = FALSE)
