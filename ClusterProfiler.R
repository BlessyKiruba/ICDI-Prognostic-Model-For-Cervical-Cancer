library(clusterProfiler)
library(org.Hs.eg.db)
#genes$PCD <- gene list with the PCD genes
# Original enrichGO result
GO_results_BP <- enrichGO(
  gene = genes$PCD,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

GO_df <- as.data.frame(GO_results_BP)

GO_top10 <- GO_df %>%
  dplyr::arrange(desc(Count)) %>%
  head(8)

#View(GO_results_BP)
fit <- barplot(GO_results_BP,
               showCategory = GO_top10$Description,
               title = "GO BP",
               font.size = 12)
fit

#View(as.data.frame(GO_results_BP))
library(ggplot2)
ggsave("GOBP ALL.tiff", plot = fit, bg = "white", width = 10, height = 8, dpi = 600)

##################
#####KEGG#########
#################
entrez_id_ <- bitr(genes$PCD,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)


kegg_enrich <- enrichKEGG(gene         = entrez_id_$ENTREZID,
                          organism     = 'hsa',       # human
                          pvalueCutoff = 0.5)

head(kegg_enrich)


#write.csv(as.data.frame(kegg_enrich), "KEGG_Results.csv")

fit_kegg <- dotplot(kegg_enrich, showCategory = 8)
fit_kegg
ggsave("KEGG.tiff", plot = fit_kegg, bg = "white", width = 9, height = 7, dpi = 600)
