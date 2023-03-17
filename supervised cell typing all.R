# Description: 
# This script contains the code for running singleR, InSituType, Seurat, CHETAH, SVM, SingleCellNet and scPred on the data 
# and outputting the results to a csv document. 
# The visualization of the results is performed in a separate script.

library(tibble)
library(readr)
library(dplyr)
library(ggplot2)
library(data.table)
library(InSituType)
library(SingleCellExperiment)
library(scater)
library(scran)
library(Seurat)

# load data
setwd("~/Desktop/Nanostring")
load("lupus nephritis SP11_1139 data.RData") # test data
HCA_reference <- readRDS("local.rds") # train data
in.glom = readRDS("in.glom.RDS") # dummy variable for glomerulus cells

# modify annot_ref and counts_ref
gene_key <- data.frame(ENSG = rownames(HCA_reference@assays[["RNA"]]@meta.features),
                       gene = as.character(HCA_reference@assays[["RNA"]]@meta.features$feature_name))
## replace gene name for CCL3L1 and RGS5
gene_key$gene[gene_key$gene == "CCL3L1"] <- "CCL3L3"
gene_key$gene[substr(gene_key$gene, 1, 4) == "RGS5"][1] <- "RGS5"
rownames(HCA_reference@assays[["RNA"]]@counts) <- gene_key$gene
counts_ref <- t(as.matrix(HCA_reference@assays[["RNA"]]@counts[colnames(counts), ]))
annot_ref <- as.character(HCA_reference@meta.data[["cell_type"]])

# pre-process test data -- all genes
sce = SingleCellExperiment(assays = list(counts = t(counts)), #make SCE
                           colData = annot)
sce = logNormCounts(sce) #normalize
sce = scater::runPCA(sce) # 61073 cells

# pre-process test data -- half genes
## subsample half of the genes
set.seed(101)
subsample_half <- sample.int(ncol(counts), ncol(counts) / 2)  #480 genes
counts_ref_half <- t(as.matrix(HCA_reference@assays[["RNA"]]@counts[colnames(counts[, subsample_half]), ]))
## pre-process test data
sce_half = SingleCellExperiment(assays = list(counts = t(counts[,subsample_half])), #make SCE
                                colData = annot)
sce_half = logNormCounts(sce_half) #normalize
sce_half = scater::runPCA(sce_half) # 61073 cells

# pre-process test data -- 10% highly variable genes
## subsample for 10% hvg
sce_hvg = SingleCellExperiment(assays = list(counts = t(counts)), #make SCE
                               colData = annot)
sce_hvg = logNormCounts(sce_hvg) #normalize
top = getTopHVGs(sce, n = 96)
sce_hvg = sce_hvg[top, ]
## pre-process
sce_hvg = scater::runPCA(sce_hvg, subset_row = top)
### Warning in (function (A, nv = 5, nu = nv, maxit = 1000, work = nv + 7, reorth = TRUE,  :
### You're computing too large a percentage of total singular values, use a standard svd instead.
counts_hvg = counts[,top]
hvg_filter = rowSums(counts_hvg) >= 5
annot_hvg = annot[hvg_filter,]
counts_hvg = counts_hvg[hvg_filter,]
sce_hvg = sce_hvg[,hvg_filter] # 61002 cells

# pre-process reference data -- all genes
sce_ref = SingleCellExperiment(assays = list(counts = t(counts_ref)),
                               colData = annot_ref)
sce_ref = sce_ref[,colSums(counts(sce_ref))>0]
sce_ref = logNormCounts(sce_ref)

# pre-process reference data -- hvg
sce_ref_hvg = sce_ref[intersect(rownames(sce_ref), top),]
sce_ref_hvg = sce_ref[,colSums(counts(sce_ref_hvg))>0]
sce_ref_hvg = logNormCounts(sce_ref_hvg)

######## Single R
library(SingleR)
# reformat reference data
meanprofiles <- InSituType:::Estep(counts = counts_ref,
                                   clust = annot_ref,
                                   neg = rep(0, length(annot_ref)))
# run SingleR
## all genes
start_time <- Sys.time()
singleR = SingleR(test = sce,
                  ref = meanprofiles,
                  labels = colnames(meanprofiles))
end_time <- Sys.time()
time_taken_all <- end_time - start_time # 1.184003 mins
## half genes
start_time <- Sys.time()
singleR_half = SingleR(test = sce_half,
                       ref = meanprofiles,
                       labels = colnames(meanprofiles))
end_time <- Sys.time()
time_taken_half <- end_time - start_time # 44.06706 secs
## hvg
start_time <- Sys.time()
singleR_hvg = SingleR(test = sce_hvg,
                      ref = meanprofiles,
                      labels = colnames(meanprofiles))
end_time <- Sys.time()
time_taken_hvg <- end_time - start_time # 3.931766 secs

# merge in.glom and singleR predictions with annot
## convert rownames to column for in.glom
in.glom <- tibble::rownames_to_column(as.data.frame(in.glom), "cell_ID")
## merge all annotations with singleR predictions
annot.pred <- merge(annot, in.glom, by = "cell_ID")
annot.pred <- merge(annot.pred, data.frame(cell_ID = rownames(singleR), singleR.labels = singleR$pruned.labels), by = "cell_ID") 
annot.pred <- merge(annot.pred, data.frame(cell_ID = rownames(singleR_half), SingleR_half.labels = singleR_half$pruned.labels), by = "cell_ID") 
annot.pred <- merge(annot.pred, data.frame(cell_ID = rownames(singleR_hvg), SingleR_hvg.labels = singleR_hvg$pruned.labels), 
                    by = "cell_ID", all.x = TRUE)

######## InSituType
library(tictoc)
tic()
# run InSituType (no cohorts)
## all genes
sup.nocohorts = insitutypeML(x = counts,
                             neg = annot$negmean,
                             reference_profiles = meanprofiles,
                             cohort = NULL)$clust
elapsed_time = toc(quiet = TRUE) # 31.249 sec
## half genes
tic()
sup.nocohorts_half = insitutypeML(x = counts[,subsample_half],
                                  neg = annot$negmean,
                                  reference_profiles = meanprofiles,
                                  cohort = NULL)$clust
elapsed_time = toc(quiet = TRUE) # 15.651 sec
## hvg
tic()
sup.nocohorts_hvg = insitutypeML(x = counts_hvg,
                                 neg = annot_hvg$negmean,
                                 reference_profiles = meanprofiles,
                                 cohort = NULL)$clust
elapsed_time = toc(quiet = TRUE) # 5.24 sec

# merge InSituType (no cohorts) predictions with annot
annot.pred <- merge(annot.pred, rownames_to_column(data.frame(sup.nocohorts), "cell_ID"), by = "cell_ID") 
annot.pred <- merge(annot.pred, rownames_to_column(data.frame(sup.nocohorts_half), "cell_ID"), by = "cell_ID") 
annot.pred <- merge(annot.pred, rownames_to_column(data.frame(sup.nocohorts_hvg), "cell_ID"), by = "cell_ID", all.x = TRUE)


# run InSituType (IF cohorts)
# create IF cohorts
## all genes and half genes
ifdata = as.matrix(annot[, c("Mean.CD298", "Mean.PanCK", "Mean.CD45", "Mean.CD20", "Mean.DAPI")])
ifcohort = fastCohorting(mat = ifdata)
## hvg
ifdata_hvg = as.matrix(annot_hvg[, c("Mean.CD298", "Mean.PanCK", "Mean.CD45", "Mean.CD20", "Mean.DAPI")])
ifcohort_hvg = fastCohorting(mat = ifdata_hvg)
## all genes
tic()
sup.ifcohorts = insitutypeML(x = counts,
                             neg = annot$negmean,
                             reference_profiles = meanprofiles, 
                             cohort = ifcohort)$clust
elapsed_time = toc(quiet = TRUE) # 32.994 sec
## half genes
tic()
sup.ifcohorts_half = insitutypeML(x = counts[,subsample_half],
                                  neg = annot$negmean,
                                  reference_profiles = meanprofiles,
                                  cohort = ifcohort)$clust
elapsed_time = toc(quiet = TRUE) # 15.729 sec
## hvg
tic()
sup.ifcohorts_hvg = insitutypeML(x = counts_hvg,
                                 neg = annot_hvg$negmean,
                                 reference_profiles = meanprofiles,
                                 cohort = ifcohort_hvg)$clust
elapsed_time = toc(quiet = TRUE) # 5.06 sec

# merge InSituType (IF cohorts) predictions with annot
annot.pred <- merge(annot.pred, rownames_to_column(data.frame(sup.ifcohorts), "cell_ID"), by = "cell_ID") 
annot.pred <- merge(annot.pred, rownames_to_column(data.frame(sup.ifcohorts_half), "cell_ID"), by = "cell_ID")
annot.pred <- merge(annot.pred, rownames_to_column(data.frame(sup.ifcohorts_hvg), "cell_ID"), by = "cell_ID", all.x = TRUE) 

######## Seurat
# create Seurat objects -- all genes
rownames(annot) = rownames(counts)
seurat.test = CreateSeuratObject(t(counts), meta.data = annot)
annot_ref <- as.data.frame(annot_ref)
rownames(annot_ref) = rownames(counts_ref)
seurat.ref = CreateSeuratObject(t(counts_ref), meta.data = annot_ref)
## normalizes the data
seurat.test = NormalizeData(seurat.test, verbose = FALSE)
seurat.ref = NormalizeData(seurat.ref, verbose = FALSE)

# create Seurat objects -- half genes (only for test data)
seurat.test_half = CreateSeuratObject(t(counts[,subsample_half]), meta.data = annot)
seurat.test_half = NormalizeData(seurat.test_half, verbose = FALSE)

# create Seurat objects -- hvg (only for test data)
seurat.test_hvg = CreateSeuratObject(t(counts_hvg), meta.data = annot_hvg)
seurat.test_hvg = NormalizeData(seurat.test_hvg, verbose = FALSE)

# find gene intersection
## all genes
gene_intersect <- intersect(colnames(counts_ref), colnames(counts))
## half genes
gene_intersect_half <- intersect(colnames(counts_ref), colnames(counts[,subsample_half]))
## hvg
gene_intersect_hvg <- intersect(colnames(counts_ref), colnames(counts_hvg))

# run Seurat
## all genes
tic()
anchors <- FindTransferAnchors(reference = seurat.ref, query = seurat.test,
                               features = gene_intersect)
seurat.predictions <- TransferData(anchorset = anchors, refdata = seurat.ref$annot_ref)
elapsed_time = toc(quiet = TRUE) # 2.987183 mins
## half genes
tic()
anchors_half <- FindTransferAnchors(reference = seurat.ref, query = seurat.test_half,
                                    features = gene_intersect_half)
seurat.predictions_half <- TransferData(anchorset = anchors_half, refdata = seurat.ref$annot_ref)
elapsed_time = toc(quiet = TRUE) # 3.98766 mins
## hvg
tic()
anchors_hvg <- FindTransferAnchors(reference = seurat.ref, query = seurat.test_hvg,
                                   features = gene_intersect_hvg)
seurat.predictions_hvg <- TransferData(anchorset = anchors_hvg, refdata = seurat.ref$annot_ref)
elapsed_time = toc(quiet = TRUE) # 5.925566 mins

# merge Seurat predictions with annot
seurat <- data.frame(cell_ID = rownames(seurat.predictions),
                     Seurat.labels = seurat.predictions$predicted.id)
seurat_half <- data.frame(cell_ID = rownames(seurat.predictions_half),
                          Seurat_half.labels = seurat.predictions_half$predicted.id) 
seurat_hvg <- data.frame(cell_ID = rownames(seurat.predictions_hvg),
                         Seurat_hvg.labels = seurat.predictions_hvg$predicted.id)
annot.pred <- Reduce(function(x, y) merge(x, y, by = "cell_ID", all.x = TRUE), 
                     list(annot.pred, seurat, seurat_half, seurat_hvg))


######## CHETAH
library(CHETAH)
# for the input we define a "counts" assay and "TSNE" reduced dimensions
## all genes
sce = runTSNE(sce)
## half genes
sce_half = runTSNE(sce_half)
## hvg
sce_hvg = runTSNE(sce_hvg)

# run CHETAH
## all genes
tic()
sce <- CHETAHclassifier(input = sce,
                        ref_cells = sce_ref,
                        ref_ct = "annot_ref")
elapsed_time = toc(quiet = TRUE) # 26.110483 mins
## half genes
tic()
sce_half <- CHETAHclassifier(input = sce_half,
                             ref_cells = sce_ref,
                             ref_ct = "annot_ref")
elapsed_time = toc(quiet = TRUE) # 34.2683 mins
## hvg
tic()
sce_hvg <- CHETAHclassifier(input = sce_hvg,
                            ref_cells = sce_ref_hvg,
                            ref_ct = "annot_ref")
elapsed_time = toc(quiet = TRUE)
# Error in ref_profiles[genes, type, drop = FALSE] : subscript out of bounds

# merge CHETAH predictions with annot
CHETAH <- data.frame(cell_ID = names(sce$celltype_CHETAH),
                     CHETAH.labels = sce$celltype_CHETAH) # 49.7% missing cell types
# table(grepl("^Node", sce$celltype_CHETAH))
# FALSE  TRUE 
# 30701 30372
CHETAH_half <- data.frame(cell_ID = names(sce_half$celltype_CHETAH),
                          CHETAH_half.labels = sce_half$celltype_CHETAH) # 53.9% missing cell types
# table(grepl("^Node", sce_half$celltype_CHETAH))
# FALSE  TRUE 
# 28157 32916

# merge CHETAH predictions with annot
annot.pred <- Reduce(function(x, y) merge(x, y, by = "cell_ID"), 
                     list(annot.pred, CHETAH, CHETAH_half))


######## SVM
library(e1071)
library(caret)

# transpose the counts
## all genes
counts_svm <- t(as.matrix(assay(sce, "logcounts")))
counts_ref_svm <- t(as.matrix(assay(sce_ref, "logcounts")))
## half gene
counts_svm_half <- t(as.matrix(assay(sce_half, "logcounts")))
counts_ref_svm_half <- counts_ref_svm[, colnames(counts_svm_half)]
## hvg
counts_svm_hvg <- t(as.matrix(assay(sce_hvg, "logcounts")))
counts_ref_svm_hvg <- counts_ref_svm[, colnames(counts_svm_hvg)]

# run svm
## all genes
tic()
svm_model <- svm(counts_ref_svm, as.factor(annot_ref$annot_ref), 
                 kernel = "linear", scale = FALSE, cost = 1)
svm.predicted <- predict(svm_model, counts_svm)
elapsed_time = toc(quiet = TRUE) # 6.478 mins
## half genes
tic()
svm_model_half <- svm(counts_ref_svm_half, as.factor(annot_ref$annot_ref), 
                      kernel = "linear", scale = FALSE, cost = 1)
svm.predicted_half <- predict(svm_model_half, counts_svm_half)
elapsed_time = toc(quiet = TRUE) # 3.5647 mins
## hvg
tic()
svm_model_hvg<- svm(counts_ref_svm_hvg, as.factor(annot_ref$annot_ref), 
                    kernel = "linear", scale = FALSE, cost = 1)
svm.predicted_hvg <- predict(svm_model_hvg, counts_svm_hvg)
elapsed_time = toc(quiet = TRUE)# 1.79163 mins

# merge SVM predictions with annot
annot.pred <- merge(annot.pred, rownames_to_column(data.frame(svm.predicted), "cell_ID"), by = "cell_ID") 
annot.pred <- merge(annot.pred, rownames_to_column(data.frame(svm.predicted_half), "cell_ID"), by = "cell_ID") 
annot.pred <- merge(annot.pred, rownames_to_column(data.frame(svm.predicted_hvg), "cell_ID"), by = "cell_ID", all.x = TRUE) 

######## SingleCellNet
library(singleCellNet)
# run SingleCellNet
## all genes
tic()
scn_model = scn_train(stTrain = annot_ref,
                      expTrain = t(counts_ref_svm),
                      nTopGenes = 10,
                      nRand = 70,
                      nTrees = 1000,
                      nTopGenePairs = 25,
                      dLevel = "annot_ref",
                      colName_samp = "row.names")
scn_predict = scn_predict(cnProc=scn_model[['cnProc']],
                          expDat=t(counts_svm),
                          nrand = 50)
elapsed_time = toc(quiet = TRUE) # ~1.5 hours
## half genes
tic()
scn_model_half = scn_train(stTrain = annot_ref,
                           expTrain = t(counts_ref_svm_half),
                           nTopGenes = 10,
                           nRand = 70,
                           nTrees = 1000,
                           nTopGenePairs = 25,
                           dLevel = "annot_ref", 
                           colName_samp = "row.names")
scn_predict_half = scn_predict(cnProc=scn_model_half[['cnProc']],
                               expDat=t(counts_svm_half),
                               nrand = 50)
elapsed_time = toc(quiet = TRUE) # 64.13271 mins
## hvg
tic()
scn_model_hvg = scn_train(stTrain = annot_ref,
                          expTrain = t(counts_ref_svm_hvg),
                          nTopGenes = 10,
                          nRand = 70,
                          nTrees = 1000,
                          nTopGenePairs = 25,
                          dLevel = "annot_ref", 
                          colName_samp = "row.names") 
# There were 29 warnings (use warnings() to see them)
scn_predict_hvg = scn_predict(cnProc=scn_model_hvg[['cnProc']],
                              expDat=t(counts_svm_hvg),
                              nrand = 50)
elapsed_time = toc(quiet = TRUE) # 39.5809 mins

# extract the predicted labels
## all genes
scn.labels = c()
scn.labels = unlist(assign_cate(scn_predict, scn.labels))
### remove the last 50 predicted labels (corresponding to the randomly sampled cells)
scn.labels = scn.labels[-c((length(scn.labels)-49):length(scn.labels))]
## half genes
scn_half.labels = c()
scn_half.labels = unlist(assign_cate(scn_predict_half, scn_half.labels))
scn_half.labels = scn_half.labels[-c((length(scn_half.labels)-49):length(scn_half.labels))]
## hvg
scn_hvg.labels = c()
scn_hvg.labels = unlist(assign_cate(scn_predict_hvg, scn_hvg.labels))
scn_hvg.labels = scn_hvg.labels[-c((length(scn_hvg.labels)-49):length(scn_hvg.labels))]

# merge singleCellNet predictions with annot
annot.pred <- cbind(annot.pred, scn.labels, scn_half.labels)
names(annot.pred)[(ncol(annot.pred)-1):ncol(annot.pred)] <- c("scn.labels", "scn_half.labels")
scn_hvg.labels <- data.frame(cell_ID = rownames(counts_svm_hvg), scn_hvg.labels)
annot.pred <- merge(annot.pred, scn_hvg.labels, by = "cell_ID", all.x = TRUE)

######## scPred (only all genes; takes too long to train the model)
library(scPred)
# process train data
scPred_ref <- seurat.ref %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:30)

tic()
# training classifiers with `scPred`
scPred_model <- getFeatureSpace(scPred_ref, "annot_ref") # 1m 54s
scPred_model <- trainModel(scPred_model, model = "svmRadial") # 04h 26m 10s
# predict cell types
scPred <- scPredict(seurat.test, scPred_model) # 3m
elapsed_time = toc(quiet = TRUE)

# merge scPred predictions with annot
scpred_prediction <- scPred$scpred_prediction
annot.pred <- cbind(annot.pred, scpred_prediction)
names(annot.pred)[ncol(annot.pred)] <- "scpred.label" # 98.52% unassigned

# output annot.pred
write.csv(annot.pred, "annotation_prediction.csv")
