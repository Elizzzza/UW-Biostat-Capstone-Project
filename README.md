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

