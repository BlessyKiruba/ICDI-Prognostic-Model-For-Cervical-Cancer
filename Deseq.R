## Batch Correction

library(sva)
Type <- factor(meta_data$Type)
Batch <- factor(meta_data$Batch)
count_matrix <- as.matrix(count_data)
corrected_counts <- ComBat_seq(count_matrix, batch = metadata$Batch)

#############

pca_res <- prcomp(t(count_data))

pca_res <- prcomp(t(corrected_counts))

pca_df <- data.frame(
  PC1 = pca_res$x[,1],
  PC2 = pca_res$x[,2],
  Batch = metadata$batch,
  Type = metadata$type
  
)
metadata$type <- dplyr::recode(metadata$Batch, "Gtex" = "TCGA")

pca_df$Type <- as.character(pca_df$Type)

library(ggplot2)

ggplot(pca_df, aes(x = PC1, y = PC2, color = "Batch")) +
  geom_point(size = 3, alpha = 0.8) +
  theme_classic() +
  ggtitle("PCA Plot Colored by Batch") +
  theme(plot.title = element_text(hjust = 0.5))

#DESEQ2
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = corrected_counts,
                              colData = metadata,
                              design = ~definition)

## Low count genes
smallestGroupSize <- 5
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

## Adding Factor levels
colData(dds)$definition <- factor(metasataa$definition, levels = c("Normal", "Cancer"))

# Check if the levels are correctly assigned
table(colData(dds)$Type)
View(counts(dds))
# DEGs
dds <- DESeq(dds)
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing
res <- results(dds)
summary(res)
vsd <- vst(dds, blind=FALSE)