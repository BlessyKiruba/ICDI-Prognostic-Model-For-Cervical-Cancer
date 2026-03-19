## Synergistic Profiling of Programmed Cell Death and Immune Responses Identifies a Novel Prognostic Index for Cervical Cancer


This repository contains the R scripts and workflow for the development of the **Immune–Cell Death Index (ICDI)**. By integrating Programmed Cell Death (PCD) and Immune Response (IR) genes, this study utilizes machine learning to improve risk stratification and treatment selection for Cervical Cancer (CC) patients.

## Repository Structure
* `Deseq.R`: Differential expression analysis to identify significant genes.
* `Unicox.R` & `Multicox.R`: Univariate and Multivariate Cox regression for survival analysis.
* `ML.R`: Machine learning models  
* `ConsensusClusterPlus.R`: Unsupervised clustering to identify patient subgroups.
* `ClusterProfiler.R`: Functional enrichment analysis (GO/KEGG).
* `Immune Deconvolution.R`: Analysis of immune cell infiltration (CIBERSORT, etc.).
* `Nomogram.R`: Generation of clinical predictive nomograms.

## Getting Started
### Prerequisites
You will need R (>= 4.0.0) and the following Bioconductor/CRAN packages:
```R
install.packages(c("survival", "survminer", "glmnet", "RMS"))
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "clusterProfiler", "ConsensusClusterPlus"))
