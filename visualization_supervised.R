# Description: 
# This script includes the visualization of the cell typing methods.

library(binom)
library(ggplot2)
library(dplyr)
library(ggthemes)

# load results
setwd("~/Desktop/Nanostring")
annot.pred <- read.csv("annotation_prediction.csv")

# load counts
load("lupus nephritis SP11_1139 data.RData")

# count the total number of transcripts per cell
annot.pred$transcripts <- rowSums(counts)

# calculate proportion of a particular cell type in glomeruli
calculate_prop <- function(data, label_col, label_val) {
  subset_data <- data[data[[label_col]] == label_val, ]
  mean(subset_data$in.glom == TRUE)
}

# full transcirpts
## all genes
sup_all <- function(label_val) {
  data.frame(Methods = c("Insitutype (no cohorts)", "Insitutype (if cohorts)", "SingleR", "Seurat", "CHETAH", "SVM", "SingleCellNet"),
             Percent = c(calculate_prop(annot.pred, "sup.nocohorts", label_val),
                         calculate_prop(annot.pred, "sup.ifcohorts", label_val),
                         calculate_prop(annot.pred, "singleR.labels", label_val),
                         calculate_prop(annot.pred, "Seurat.labels", label_val),
                         calculate_prop(annot.pred, "CHETAH.labels", label_val),
                         calculate_prop(annot.pred, "svm.predicted", label_val),
                         calculate_prop(annot.pred, "scn.labels", label_val))) %>%
    mutate(Percent = Percent*100,
           n = nrow(annot.pred[annot.pred$in.glom == "TRUE",]),
           correct = Percent * n / 100,
           ci_low = binom.confint(correct, n, methods = "wilson")$lower*100,
           ci_high = binom.confint(correct, n, methods = "wilson")$upper*100,
           cells = "All cells",
           genes = "All genes",
           cell_type = label_val)
}

results_all <- rbind(sup_all("podocyte"), 
                     sup_all("glomerular endothelial cell"), 
                     sup_all("epithelial cell of proximal tubule"))

## half genes
sup_half <- function(label_val) {
  data.frame(Methods = c("Insitutype (no cohorts)", "Insitutype (if cohorts)", "SingleR", "Seurat", "CHETAH", "SVM", "SingleCellNet"),
             Percent = c(calculate_prop(annot.pred, "sup.nocohorts_half", label_val),
                         calculate_prop(annot.pred, "sup.ifcohorts_half", label_val),
                         calculate_prop(annot.pred, "SingleR_half.labels", label_val),
                         calculate_prop(annot.pred, "Seurat_half.labels", label_val),
                         calculate_prop(annot.pred, "CHETAH_half.labels", label_val),
                         calculate_prop(annot.pred, "svm.predicted_half", label_val),
                         calculate_prop(annot.pred, "scn_half.labels", label_val))) %>%
    mutate(Percent = Percent*100,
           n = nrow(annot.pred[annot.pred$in.glom == "TRUE",]),
           correct = Percent * n / 100,
           ci_low = binom.confint(correct, n, methods = "wilson")$lower*100,
           ci_high = binom.confint(correct, n, methods = "wilson")$upper*100,
           cells = "All cells",
           genes = "Half genes",
           cell_type = label_val)
}

results_all <- rbind(results_all,
                     sup_half("podocyte"), 
                     sup_half("glomerular endothelial cell"), 
                     sup_half("epithelial cell of proximal tubule"))

## hvg
sup_hvg <- function(label_val) {
  data.frame(Methods = c("Insitutype (no cohorts)", "Insitutype (if cohorts)", "SingleR", "Seurat", "SVM", "SingleCellNet"),
             Percent = c(calculate_prop(na.omit(annot.pred), "sup.nocohorts_hvg", label_val),
                         calculate_prop(na.omit(annot.pred), "sup.ifcohorts_hvg", label_val),
                         calculate_prop(na.omit(annot.pred), "SingleR_hvg.labels", label_val),
                         calculate_prop(na.omit(annot.pred), "Seurat_hvg.labels", label_val),
                         calculate_prop(na.omit(annot.pred), "svm.predicted_hvg", label_val),
                         calculate_prop(na.omit(annot.pred), "scn_hvg.labels", label_val))) %>%
    mutate(Percent = Percent*100,
           n = nrow(na.omit(annot.pred)[na.omit(annot.pred)$in.glom == "TRUE",]),
           correct = Percent * n / 100,
           ci_low = binom.confint(correct, n, methods = "wilson")$lower*100,
           ci_high = binom.confint(correct, n, methods = "wilson")$upper*100,
           cells = "All cells",
           genes = "hvg",
           cell_type = label_val)
}

results_all <- rbind(results_all,
                     sup_hvg("podocyte"), 
                     sup_hvg("glomerular endothelial cell"), 
                     sup_hvg("epithelial cell of proximal tubule"))


# calculate results for a given cell type and transcript range
## all genes
sup_all_subset <- function(label_val, transcripts_min, transcripts_max) {
  subset <- annot.pred %>% filter(transcripts >= transcripts_min & transcripts <= transcripts_max)
  
  results <- data.frame(Methods = c("Insitutype (no cohorts)", "Insitutype (if cohorts)", "SingleR", "Seurat", "CHETAH", "SVM", "SingleCellNet"),
                        Percent = c(calculate_prop(subset, "sup.nocohorts", label_val),
                                    calculate_prop(subset, "sup.ifcohorts", label_val),
                                    calculate_prop(subset, "singleR.labels", label_val),
                                    calculate_prop(subset, "Seurat.labels", label_val),
                                    calculate_prop(subset, "CHETAH.labels", label_val),
                                    calculate_prop(subset, "svm.predicted", label_val),
                                    calculate_prop(subset, "scn.labels", label_val))) %>%
    mutate(Percent = Percent*100,
           n = nrow(subset[subset$in.glom == "TRUE",]),
           correct = Percent * n / 100,
           ci_low = binom.confint(correct, n, methods = "wilson")$lower*100,
           ci_high = binom.confint(correct, n, methods = "wilson")$upper*100,
           cells = paste(transcripts_min, "-", transcripts_max),
           genes = "All genes", 
           cell_type = label_val)
  
  return(results)
}

results_all <- rbind(results_all,
                     sup_all_subset("podocyte", 0, 50), 
                     sup_all_subset("glomerular endothelial cell", 0, 50),
                     sup_all_subset("epithelial cell of proximal tubule", 0, 50),
                     sup_all_subset("podocyte", 51, 100), 
                     sup_all_subset("glomerular endothelial cell", 51, 100), 
                     sup_all_subset("epithelial cell of proximal tubule", 51, 100),
                     sup_all_subset("podocyte", 101, 200), 
                     sup_all_subset("glomerular endothelial cell", 101, 200), 
                     sup_all_subset("epithelial cell of proximal tubule", 101, 200),
                     sup_all_subset("podocyte", 201, Inf), 
                     sup_all_subset("glomerular endothelial cell", 201, Inf), 
                     sup_all_subset("epithelial cell of proximal tubule", 201, Inf)
                     )

## half genes
sup_half_subset <- function(label_val, transcripts_min, transcripts_max) {
  subset <- annot.pred %>% filter(transcripts >= transcripts_min & transcripts <= transcripts_max)
  
  results <- data.frame(Methods = c("Insitutype (no cohorts)", "Insitutype (if cohorts)", "SingleR", "Seurat", "CHETAH", "SVM", "SingleCellNet"),
                        Percent = c(calculate_prop(subset, "sup.nocohorts_half", label_val),
                                    calculate_prop(subset, "sup.ifcohorts_half", label_val),
                                    calculate_prop(subset, "SingleR_half.labels", label_val),
                                    calculate_prop(subset, "Seurat_half.labels", label_val),
                                    calculate_prop(subset, "CHETAH_half.labels", label_val),
                                    calculate_prop(subset, "svm.predicted_half", label_val),
                                    calculate_prop(subset, "scn_half.labels", label_val))) %>%
    mutate(Percent = Percent*100,
           n = nrow(subset[subset$in.glom == "TRUE",]),
           correct = Percent * n / 100,
           ci_low = ifelse(!is.na(correct), binom.confint(na.omit(correct), n, methods = "wilson")$lower*100, NaN),
           ci_high = ifelse(!is.na(correct), binom.confint(na.omit(correct), n, methods = "wilson")$upper*100, NaN),
           cells = paste(transcripts_min, "-", transcripts_max),
           genes = "Half genes", 
           cell_type = label_val)
  
  return(results)
}

results_all <- rbind(results_all,
                     sup_half_subset("podocyte", 0, 50), # no podocyte in scn results
                     sup_half_subset("glomerular endothelial cell", 0, 50),
                     sup_half_subset("epithelial cell of proximal tubule", 0, 50),
                     sup_half_subset("podocyte", 51, 100), 
                     sup_half_subset("glomerular endothelial cell", 51, 100), 
                     sup_half_subset("epithelial cell of proximal tubule", 51, 100),
                     sup_half_subset("podocyte", 101, 200), 
                     sup_half_subset("glomerular endothelial cell", 101, 200), 
                     sup_half_subset("epithelial cell of proximal tubule", 101, 200),
                     sup_half_subset("podocyte", 201, Inf), 
                     sup_half_subset("glomerular endothelial cell", 201, Inf), 
                     sup_half_subset("epithelial cell of proximal tubule", 201, Inf)
)

## hvg
sup_hvg_subset <- function(label_val, transcripts_min, transcripts_max) {
  subset <- annot.pred %>% filter(transcripts >= transcripts_min & transcripts <= transcripts_max)
  
  results <- data.frame(Methods = c("Insitutype (no cohorts)", "Insitutype (if cohorts)", "SingleR", "Seurat", "SVM", "SingleCellNet"),
                        Percent = c(calculate_prop(na.omit(subset), "sup.nocohorts_hvg", label_val),
                                    calculate_prop(na.omit(subset), "sup.ifcohorts_hvg", label_val),
                                    calculate_prop(na.omit(subset), "SingleR_hvg.labels", label_val),
                                    calculate_prop(na.omit(subset), "Seurat_hvg.labels", label_val),
                                    calculate_prop(na.omit(subset), "svm.predicted_hvg", label_val),
                                    calculate_prop(na.omit(subset), "scn_hvg.labels", label_val))) %>%
    mutate(Percent = Percent*100,
           n = nrow(subset[subset$in.glom == "TRUE",]),
           correct = Percent * n / 100,
           ci_low = ifelse(!is.na(correct), binom.confint(na.omit(correct), n, methods = "wilson")$lower*100, NaN),
           ci_high = ifelse(!is.na(correct), binom.confint(na.omit(correct), n, methods = "wilson")$upper*100, NaN),
           cells = paste(transcripts_min, "-", transcripts_max),
           genes = "hvg", 
           cell_type = label_val)
  
  return(results)
}

results_all <- rbind(results_all,
                     sup_hvg_subset("podocyte", 0, 50),
                     sup_hvg_subset("glomerular endothelial cell", 0, 50),
                     sup_hvg_subset("epithelial cell of proximal tubule", 0, 50),
                     sup_hvg_subset("podocyte", 51, 100), 
                     sup_hvg_subset("glomerular endothelial cell", 51, 100), 
                     sup_hvg_subset("epithelial cell of proximal tubule", 51, 100),
                     sup_hvg_subset("podocyte", 101, 200), 
                     sup_hvg_subset("glomerular endothelial cell", 101, 200), 
                     sup_hvg_subset("epithelial cell of proximal tubule", 101, 200),
                     sup_hvg_subset("podocyte", 201, Inf), 
                     sup_hvg_subset("glomerular endothelial cell", 201, Inf), 
                     sup_hvg_subset("epithelial cell of proximal tubule", 201, Inf)
)

# transform cells and methods into factor
results_all$Methods <- factor(results_all$Methods, 
                              levels = c("Insitutype (no cohorts)", "Insitutype (if cohorts)", "SingleR", "Seurat", "CHETAH",
                                         "SVM", "SingleCellNet"))
results_all$cells <- factor(results_all$cells, 
                            levels = c("All cells", "0 - 50", "51 - 100", "101 - 200", "201 - Inf"),
                            labels = c("All cells", "0-50", "51-100", "101-200", ">200"))

# create a grouped bar chart for each cell type
plot_cell_type <- function(celltype) {
  results_all %>% 
    filter(cell_type == celltype) %>%
    ggplot(aes(x = cells, y = Percent, fill = Methods)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0, position = position_dodge(0.9))+
    theme_bw() +
    ylim(if(celltype == "epithelial cell of proximal tubule") {c(0, 25)} else{c(0, 100)}) + 
    labs(x = "Transcripts", y = paste0("% of ", celltype, " in glomeruli"), fill = "Methods") +
    scale_fill_manual(values = pal[c(1, 7:9, 2, 6, 4)]) + 
    facet_grid(~genes) + 
    if (celltype == "epithelial cell of proximal tubule") {
      geom_hline(yintercept = 0, color = "red")
    } else {
      geom_hline(yintercept = 100, color = "red")
    }
  }

p_podocyte <- plot_cell_type("podocyte")
p_epithelial <- plot_cell_type("glomerular endothelial cell")
p_tubule <- plot_cell_type("epithelial cell of proximal tubule")

ggsave(filename = "p_podocyte.png", p_podocyte & theme(text = element_text(size = 20)), 
       width = 24, height = 8, device = png) 
ggsave(filename = "p_epithelial.png", p_epithelial & theme(text = element_text(size = 20)), 
       width = 24, height = 8, device = png) 
ggsave(filename = "p_tubule.png", p_tubule & theme(text = element_text(size = 20)), 
       width = 24, height = 8, device = png) 

# specify immune cells and cytokeratin cells
immune <- c("B cell", "CD4-positive, alpha-beta T cell", "CD8-positive, alpha-beta T cell", 
            "classical monocyte", "dendritic cell", "kidney resident macrophage","mature NK T cell",
            "natural killer cell", "neutrophil", "non-classical monocyte", "plasmacytoid dendritic cell")

cytokeratin <- c("epithelial cell of proximal tubule", "kidney connecting tubule epithelial cell",
                 "kidney loop of Henle thick ascending limb epithelial cell", "renal alpha-intercalated cell", 
                 "renal beta-intercalated cell", "renal intercalated cell", "renal principal cell")


# Calculates the log2-fold difference between cytokeratin+ and cytokeratin- groups
log2_diff_panCK <- function(label_col, transcripts_min, transcripts_max) {
  subset <- annot.pred %>% filter(transcripts >= transcripts_min & transcripts <= transcripts_max)
  pos_mean <- mean(subset[subset[[label_col]] %in% cytokeratin, "Mean.PanCK"])
  neg_mean <- mean(subset[!subset[[label_col]] %in% cytokeratin, "Mean.PanCK"])
  log2(pos_mean/neg_mean)
}

## all genes
panCK_all_subset <- function(transcripts_min, transcripts_max) {
  Methods <- c("Insitutype (no cohorts)", "Insitutype (if cohorts)", "SingleR", "Seurat", "CHETAH", "SVM", "SingleCellNet")
  panCK_fold <- sapply(c("sup.nocohorts", "sup.ifcohorts", "singleR.labels", "Seurat.labels", "CHETAH.labels", "svm.predicted", "scn.labels"), 
                       log2_diff_panCK, transcripts_min, transcripts_max)
  cells <- paste(transcripts_min, "-", transcripts_max)
  genes <- "All genes"
  data.frame(Methods, panCK_fold, cells, genes)
}

panCK_results_all <- rbind(panCK_all_subset(0, Inf), panCK_all_subset(0, 50), panCK_all_subset(51, 100), 
                           panCK_all_subset(101, 200), panCK_all_subset(201, Inf))

## half genes
panCK_half_subset <- function(transcripts_min, transcripts_max) {
  Methods <- c("Insitutype (no cohorts)", "Insitutype (if cohorts)", "SingleR", "Seurat", "CHETAH", "SVM", "SingleCellNet")
  panCK_fold <- sapply(c("sup.nocohorts_half", "sup.ifcohorts_half", "SingleR_half.labels", "Seurat_half.labels", 
                         "CHETAH_half.labels", "svm.predicted_half", "scn_half.labels"), 
                       log2_diff_panCK, transcripts_min, transcripts_max)
  cells <- paste(transcripts_min, "-", transcripts_max)
  genes <- "Half genes"
  data.frame(Methods, panCK_fold, cells, genes)
}

panCK_results_all <- rbind(panCK_results_all, panCK_half_subset(0, Inf), panCK_half_subset(0, 50), 
                           panCK_half_subset(51, 100), panCK_half_subset(101, 200), panCK_half_subset(201, Inf))

## hvg
panCK_hvg_subset <- function(transcripts_min, transcripts_max) {
  Methods <- c("Insitutype (no cohorts)", "Insitutype (if cohorts)", "SingleR", "Seurat", "SVM", "SingleCellNet")
  panCK_fold <- sapply(c("sup.nocohorts_hvg", "sup.ifcohorts_hvg", "SingleR_hvg.labels", "Seurat_hvg.labels", 
                         "svm.predicted_hvg", "scn_hvg.labels"), 
                       log2_diff_panCK, transcripts_min, transcripts_max)
  cells <- paste(transcripts_min, "-", transcripts_max)
  genes <- "hvg"
  data.frame(Methods, panCK_fold, cells, genes)
}

panCK_results_all <- rbind(panCK_results_all, panCK_hvg_subset(0, Inf), panCK_hvg_subset(0, 50), 
                           panCK_hvg_subset(51, 100), panCK_hvg_subset(101, 200), panCK_hvg_subset(201, Inf))

# transform cells and methods into factor
panCK_results_all$Methods <- factor(panCK_results_all$Methods, 
                                    levels = c("Insitutype (no cohorts)", "Insitutype (if cohorts)", "SingleR", "Seurat", "CHETAH", 
                                               "SVM", "SingleCellNet"))
panCK_results_all$cells <- factor(panCK_results_all$cells, 
                                  levels = c("0 - Inf", "0 - 50", "51 - 100", "101 - 200", "201 - Inf"),
                                  labels = c("All cells", "0-50", "51-100", "101-200", ">200"))

# create a grouped bar chart
panCK_results_all %>% 
  ggplot(aes(x = cells, y = panCK_fold, fill = Methods)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  labs(x = "Transcripts",  y = "Log2-fold difference in PanCK", fill = "Methods") +
  scale_fill_manual(values = pal[c(1, 7:9, 2, 6, 4)]) + 
  facet_grid(~genes) -> p_diff_stain_panCK

ggsave(filename = "p_diff_stain_panCK.png", p_diff_stain_panCK & theme(text = element_text(size = 20)), 
       width = 24, height = 8, device = png)

# Calculates the log2-fold difference between immune and non-immune groups
log2_diff_CD45 <- function(label_col, transcripts_min, transcripts_max) {
  subset <- annot.pred %>% filter(transcripts >= transcripts_min & transcripts <= transcripts_max)
  pos_mean <- mean(subset[subset[[label_col]] %in% immune, "Mean.CD45"])
  neg_mean <- mean(subset[!subset[[label_col]] %in% immune, "Mean.CD45"])
  log2(pos_mean/neg_mean)
}

## all genes
CD45_all_subset <- function(transcripts_min, transcripts_max) {
  Methods <- c("Insitutype (no cohorts)", "Insitutype (if cohorts)", "SingleR", "Seurat", "CHETAH", "SVM", "SingleCellNet")
  CD45_fold <- sapply(c("sup.nocohorts", "sup.ifcohorts", "singleR.labels", "Seurat.labels", "CHETAH.labels", "svm.predicted", "scn.labels"), 
                       log2_diff_CD45, transcripts_min, transcripts_max)
  cells <- paste(transcripts_min, "-", transcripts_max)
  genes <- "All genes"
  data.frame(Methods, CD45_fold, cells, genes)
}

CD45_results_all <- rbind(CD45_all_subset(0, Inf), CD45_all_subset(0, 50), CD45_all_subset(51, 100), 
                           CD45_all_subset(101, 200), CD45_all_subset(201, Inf))

## half genes
CD45_half_subset <- function(transcripts_min, transcripts_max) {
  Methods <- c("Insitutype (no cohorts)", "Insitutype (if cohorts)", "SingleR", "Seurat", "CHETAH", "SVM", "SingleCellNet")
  CD45_fold <- sapply(c("sup.nocohorts_half", "sup.ifcohorts_half", "SingleR_half.labels", "Seurat_half.labels", 
                         "CHETAH_half.labels", "svm.predicted_half", "scn_half.labels"), 
                       log2_diff_CD45, transcripts_min, transcripts_max)
  cells <- paste(transcripts_min, "-", transcripts_max)
  genes <- "Half genes"
  data.frame(Methods, CD45_fold, cells, genes)
}

CD45_results_all <- rbind(CD45_results_all, CD45_half_subset(0, Inf), CD45_half_subset(0, 50), 
                          CD45_half_subset(51, 100), CD45_half_subset(101, 200), CD45_half_subset(201, Inf))

## hvg
CD45_hvg_subset <- function(transcripts_min, transcripts_max) {
  Methods <- c("Insitutype (no cohorts)", "Insitutype (if cohorts)", "SingleR", "Seurat", "SVM", "SingleCellNet")
  CD45_fold <- sapply(c("sup.nocohorts_hvg", "sup.ifcohorts_hvg", "SingleR_hvg.labels", "Seurat_hvg.labels", 
                         "svm.predicted_hvg", "scn_hvg.labels"), 
                       log2_diff_CD45, transcripts_min, transcripts_max)
  cells <- paste(transcripts_min, "-", transcripts_max)
  genes <- "hvg"
  data.frame(Methods, CD45_fold, cells, genes)
}

CD45_results_all <- rbind(CD45_results_all, CD45_hvg_subset(0, Inf), CD45_hvg_subset(0, 50), 
                          CD45_hvg_subset(51, 100), CD45_hvg_subset(101, 200), CD45_hvg_subset(201, Inf))

# transform cells and methods into factor
CD45_results_all$Methods <- factor(CD45_results_all$Methods, 
                                    levels = c("Insitutype (no cohorts)", "Insitutype (if cohorts)", "SingleR", "Seurat", "CHETAH", 
                                               "SVM", "SingleCellNet"))
CD45_results_all$cells <- factor(CD45_results_all$cells, 
                                  levels = c("0 - Inf", "0 - 50", "51 - 100", "101 - 200", "201 - Inf"),
                                  labels = c("All cells", "0-50", "51-100", "101-200", ">200"))

# create a grouped bar chart
CD45_results_all %>% 
  ggplot(aes(x = cells, y = CD45_fold, fill = Methods)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  labs(x = "Transcripts",  y = "Log2-fold difference in CD45", fill = "Methods") +
  scale_fill_manual(values = pal[c(1, 7:9, 2, 6, 4)]) + 
  facet_grid(~genes) -> p_diff_stain_CD45

ggsave(filename = "p_diff_stain_CD45.png", p_diff_stain_CD45 & theme(text = element_text(size = 20)), 
       width = 24, height = 8, device = png)




