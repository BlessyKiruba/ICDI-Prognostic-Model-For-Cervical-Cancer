library(dplyr)
library(survival)

cox_res <- list()

gene_list <- c("PCDI_4", "IRGI_4", "ICDS_4", "Age","Stage",  "T", "N","M")

for (gene in gene_list) { #initiate the loop for each gene in the list
  if (gene %in% rownames(multi_tcga)) { # check if they are present in the exprssion datafram
    
    gene_expr <- expr_df[gene, , drop = FALSE] #one gene with its exp data
    
    gene_expr_long <- data.frame(
      sample = colnames(gene_expr), #extract the barcode and expression details of the one gene
      expression = as.numeric(gene_expr), #make that coloumn as numeric
      stringsAsFactors = FALSE
    )
    
    merged <- merge(gene_expr_long, df10, by = "sample")
    
    merged$ <- as.numeric(as.character(merged$days_to_last_known_alive))
    
    result <- tryCatch({
      cox_model <- coxph(Surv(days_to_last_known_alive, deceased) ~ expression, data = merged)
      cox_summary <- summary(cox_model)
      
      # Check proportional hazards assumption
      zph_result <- cox.zph(cox_model)
      
      list(
        model = cox_model,
        HR = cox_summary$coefficients[, "exp(coef)"],
        pvalue = cox_summary$coefficients[, "Pr(>|z|)"], # Extract p-value for expression
        CI = cox_summary$conf.int[, c("lower .95", "upper .95")]
      )
    }, error = function(e) {
      message(sprintf("Skipping gene %s due to error: %s", gene, e$message))
      NULL
    })
    
    cox_res[[gene]] <- result
  }
}

# in table cox result
# Initialize empty data frame
cox_df <- data.frame(
  gene = character(),
  HR = numeric(),
  pvalue = numeric(),
  CI_lower = numeric(),
  CI_upper = numeric(),
  stringsAsFactors = FALSE
)

for (gene in names(cox_res)) {
  res <- cox_res[[gene]]
  if (!is.null(res)) {
    cox_df <- rbind(
      cox_df,
      data.frame(
        gene = gene,
        HR = res$HR,
        pvalue = res$pvalue,
        CI_lower = res$CI[1],
        CI_upper = res$CI[2],
        stringsAsFactors = FALSE
      )
    )
  }
}

library(writexl)
write_xlsx(cox_df, "cox_irg.xlsx")