if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ConsensusClusterPlus")
library(ConsensusClusterPlus)

df <- as.list(c("FADD", "MUC4", "CLNK", "CD8B"))
df_TCGA <- as.data.frame(df_TCGA)
rownames(df_TCGA) <- df_TCGA$rownames
df1 <- df_TCGA[, colnames(df_TCGA) %in% df] #to change the count matrix

dim(df1)
df1<- as.matrix(df1)
title <- "C:/Users/Documents/ConsensusClusterPlus_1"
df1_t <- t(df1)

results = ConsensusClusterPlus(as.matrix(df1_t),
                               maxK=6,
                               reps=100,
                               pItem=0.8,
                               pFeature=1,
                               tmyPal = c('white','#008080'),
                               title='ConsensusCluster/',
                               clusterAlg="km",
                               distance="euclidean",
                               seed=123456,
                               plot = "pdf")

# Choose the optimal number of clusters (e.g., K = 3)
K <- 2

# Extract sample-to-cluster assignments
cluster_assignments <- results[[K]]$consensusClass
# Combine cluster assignments with your expression data
# Transpose so samples are rows
df_clustered <- cbind(Cluster = cluster_assignments, t(df1_t))

# Check result
View(df_clustered)
df_clustered <- as.data.frame(df_clustered)
df_clustered$rownames <- rownames(df_clustered)
library(writexl)
write_xlsx(df_clustered, "Consensusclusplus_k2.xlsx")

