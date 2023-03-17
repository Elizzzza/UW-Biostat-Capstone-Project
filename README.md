# UW-Biostat-Capstone-Project

**Benchmark Analysis of InSituType Cell Typing Model for Spatial Transcriptomics** 

Sponsor: NanoString Technologies, Inc 

Project Members: Howard Baek, Eliza Chai, Alexis Harris, Ingrid Luo, Makayla Tang 

## Background

### Cell Typing and InSituType 

Cell typing is a fundamental aspect of single cell transcriptomics analysis, which involves categorizing cells based on their expression profiles. While unsupervised clustering is commonly used, supervised classification can also be employed in tissues with high-quality reference datasets to avoid clustering's instability and interpretation difficulties. Spatial transcriptomics platforms detect individual RNA molecules in tissue samples and maintain cells' location data, producing datasets that can be used for cell typing. In this context, the InSituType algorithm is introduced for cell typing in spatial transcriptomics data.  

NanoString’s cell typing model InSituType utilizes a likelihood model to enable supervised cell typing from reference datasets via a Bayes classifier and unsupervised or semi-supervised cell typing via an Expectation Maximization algorithm. InSituType is designed to address the challenges of sparse and multi-modal spatial transcriptomics data, and it employs an escalating subsampling scheme to handle large datasets efficiently. Compared to existing cell typing methods, InSituType offers several advantages, including its ability to handle large datasets, its use of a likelihood model to weigh the evidence from every transcript in a cell, its incorporation of alternative data types such as images and spatial context, and its ability to identify new clusters alongside reference cell types. 

## Objectives 

We have two objectives for our projects. The first objective is to understand how the performance of the unsupervised InSituType method varies due to tuning parameters of interest. We investigated the InSituType parameters relating to the number of clusters, number of iterations, subsample size for random starts, and number of cohorts. To measure the performance of the model, we used the Adjusted Rand Index, Bayesian information criterion, and recorded computation time to run the model on a consistent computing environment. The second objective is to assess the performance of the supervised InSituType method and other open-source supervised algorithms using a kidney biopsy dataset.  

## Data sources 

To benchmark InSituType, we utilized two samples of tissue datasets. To analyze the performance of the unsupervised method, we use a dataset obtained from a "cell pellet array" (CPA) where various cell lines are pelleted onto a slide. This CPA was profiled using the CosMx™ Spatial Molecular Imager. Each field of view (FOV) comprises cells solely from a distinct cell line, providing us with accurate cell identity information. This dataset contains information for 52,518 cells and has thirteen known clusters, or cell lines.  

To benchmark the performance of the InSituType supervised cell typing method, we used a dataset from a kidney biopsy taken from a lupus nephritis patient. This dataset was also profiled using the CosMx™ Spatial Molecular Imager and is a single cell RNA dataset with 61,073 mature kidney cells. We perform cell typing on this dataset using reference profiles from the Kidney Single Cell Atlas. This training dataset is from a single cell RNA study that performed an analysis of all the cell types in the human kidney and is comprised of 40,268 mature kidney cells and 33 cell types. 

## Overview of the Method 

### Benchmarking unsupervised clustering 

We examine the performance of unsupervised clustering method by tuning parameters of interest. In this analysis, we tuned the parameters for the number of clusters, the number of iterations, the subsample size for phase 1 random starts, and the number of cohorts argument in the “fastCohorting ” function (used to quickly split cells into cohorts). Benchmarking was performed on the full CPA dataset and two down sampled CPA datasets containing a random sample of 50% of the genes and the 10% most highly variable genes (HVG)​. 

To access performance, we used three evaluation metrics: Bayesian Information Criterion (BIC), Adjusted Rand Index (ARI), and computation time.  With lower BIC values, we obtain models that better fit the data. A high ARI indicates high similarity between the true cell type labels and the clustering partition. ARI equals 1 when there is perfect correspondence and takes a value near 0 for a random partition. The recorded system run time allows us to access the computation resources required when running on a consistent computing environment. 

- `Number of clusters`: The number of clusters is the number of Field of Views (FOV) in this dataset. The field of view for a microscope is the extent of the observable area in distance unit. However, this may not always give optimal fit and accuracy. In this CPA dataset used for benchmarking, we know there are thirteen cell types, however this information will not always be known when performing unsupervised cell typing on spatial transcriptomic data. To tune the number of clusters parameter, we tested a range of 10 to 25 clusters to observe the change in performance due to the initial estimate of the number of clusters.  

- `Number of iterations`: InSituType iteratively runs an Expectation Maximization (EM) algorithm on subsamples of cells in the input dataset. With more iterations, we are more likely to get a better fit at the cost of longer runtime. The default value is 10 iterations, but we run the function with values ranging from 1 to 20 iterations and evaluate the performance with the three proposed metrices. 

- `Number of cells in Phase 1`: The number of cells in phase 1 represents the subsample cell size for phase 1 random starts. The default cell sample size for phase 1 is 10,000 cells. InSituType randomly selects 10,000 cells from the input dataset and performs early iterations of EM. The early iterations analyzed in phase 1 are used for an escalating subsampling scheme to quickly approach an approximate solution in large datasets. After phase 1 is complete, phase 2 iterates the EM algorithm on a subsample of 20,000 cells to obtain a more precise solution, and then phase 3 iterates over 100,000 cells. The final phase analyzes the whole dataset once where cells are classified using profiles obtained from previous phases.  

  To tune the number of cells in phase 1, we investigated a subsample size range from 1,000 to 20,000 cells, indexing by 1,000 cells and evaluated the performance using the three proposed metrics. 

- `Number of cohorts`: 
  “Cohorting” is a pre-clustering step performed before InSituType is run. The use of the cohort parameter in InSituType incorporates alternative data types. Examples of these alternative data types include information on the spatial context of the cell types in tissue, images, and immunofluorescence stains with diverse distributions. The cohort membership defines the cell’s prior probability of belonging to each cell type, in other words informing the likelihood calculations in the EM algorithm. We have no idea how fine-grained these pre-clusters should be. As a result, we explore values ranging from 0 to 100 to see how the default value of 25 performs. Similar to the other parameters, we evaluated the performance using ARI, BIC, and system time and tested only on the full CPA dataset. 

### Benchmarking supervised cell typing

To perform supervised cell typing, it is necessary to have both a training set and a test set. Our test set is a kidney biopsy from a lupus nephritis patient, and the Kidney Single Cell Atlas serves as the training set which provides average expression profiles for each kidney cell type based on their single-cell RNA-seq studies. Since there is no definitive gold standard cell type data for this sample, we cannot determine with certainty whether a cell type call is accurate. To evaluate the quality of our cell typing results, we rely on known kidney biology and use two metrics. 

The first metric evaluates the accuracy of glomerular cell assignments to glomeruli, as kidneys contain substructures called glomeruli that comprise distinct cell types. Our expectation is that 100% of podocytes and glomerular endothelial cells will be assigned to glomeruli, while 0% of proximal tubule cells will be assigned to glomeruli. 

We use bar charts ([Figure1](https://github.com/Elizzzza/UW-Biostat-Capstone-Project/blob/main/supervised_figures/p_podocyte.png), [Figure2](https://github.com/Elizzzza/UW-Biostat-Capstone-Project/blob/main/supervised_figures/p_epithelial.png), [Figure3](https://github.com/Elizzzza/UW-Biostat-Capstone-Project/blob/main/supervised_figures/p_tubule.png)) to present the first metrics, showing the percentage of each cell type in the glomeruli, along with confidence intervals. We use a red line to indicate the target percentage for each cell type (i.e., 100% for podocytes and glomerular endothelial cells, and 0% for proximal tubule cells).  The cells are grouped based on their transcript expression level, including all cells and cells with 0-50, 51-100, 101-200, and >200 transcripts. 

We perform benchmarking on three datasets: the full dataset, and two down sampled datasets which are a random sample of 50% of the genes and the top 10% most highly variable genes. We do this in response to certain company competitors who offer products focused on analyzing the top 10% most highly variable genes. 

Second, we record the intensity of a few "marker proteins" that are highly informative of cell type for each cell. We compare the marker protein mean expression intensity between cell types that are either positive or negative for the markers, such as CD45 for immune cells versus other cell types, or PanCK for cytokeratin+ versus cytokeratin- cell types.  

Similar to the first metrics, the second metrics use bar charts ([Figure4](https://github.com/Elizzzza/UW-Biostat-Capstone-Project/blob/main/supervised_figures/p_diff_stain_panCK.png), [Figure5](https://github.com/Elizzzza/UW-Biostat-Capstone-Project/blob/main/supervised_figures/p_diff_stain_CD45.png)) to display the results. We present the logarithmic difference in mean PanCK stain intensity between cytokeratin+ and cytokeratin- cells ([Figure4](https://github.com/Elizzzza/UW-Biostat-Capstone-Project/blob/main/supervised_figures/p_diff_stain_panCK.png)), and the logarithmic difference in CD45 stain intensity between immune and non-immune cells ([Figure5](https://github.com/Elizzzza/UW-Biostat-Capstone-Project/blob/main/supervised_figures/p_diff_stain_CD45.png)). The classification of cell types as either cytokeratin-positive or cytokeratin-negative, as well as immune or non-immune, is determined through biological knowledge. We use logarithm base 2 ratio instead of a simple difference between the mean stain intensity for two cell groups because a logarithmic scale can better represent the fold difference between groups. 

To evaluate the performance of our supervised cell typing approach, we compare InSituType to SingleR, CHETAH, SeuratV3, SingleCellNet, scPred, and SVM. Additionally, we conduct all analyses in a consistent computing environment and record computation times. 

- **SingleR**: An unbiased cell typing method for scRNA-seq by leveraging reference transcriptomic datasets of pure cell types to infer the cell of origin of each single cell independently. 

- **CHETAH** (Characterization of Cell Types Aided by Hierarchical classification): A scRNA-seq classifier by hierarchical clustering of the reference data. The classification tree enables a step-wise, top-to-bottom classification.  

- **SeuratV3**: A single-cell transcriptomics classifier can anchor diverse datasets together, enabling us to integrate single-cell measurements not only across scRNA-seq, but also across 	different modalities (e.g. scATAC-seq).  

- **SingleCellNet**: A random forest classifier to learn cell type-specific gene pairs from cross-platform and cross-species datasets and thus quantitatively assesses cell identity at a single-cell resolution. 

- **ScPred**: A scRNA-seq classifier by using a combination of unbiased feature selection from a 	reduced-dimensional space (e.g. PCA), and machine-learning probability-based prediction methods.  

- **SVM** (support vector machines): Supervised learning models with associated learning algorithms that analyze data used for classification and regression analysis.  

