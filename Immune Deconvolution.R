library(readxl)
library(e1071)
library(preprocessCore)
library(future)
library(furrr)
library(parallel)
library(purrr)
library(CIBERSORT)
library(immunedeconv)
deconvolution_methods
#MCPcounter             EPIC        quanTIseq            xCell 
#"mcp_counter"           "epic"      "quantiseq"          "xcell" 
#CIBERSORT CIBERSORT (abs.)            TIMER     ConsensusTME 
#"cibersort"  "cibersort_abs"          "timer"  "consensus_tme" 
#ABIS         ESTIMATE 
#"abis"       "estimate" 


# Create the indication vector
indications <- rep("CESC", ncol(gene_exp))
# Run TIMER
timer <- immunedeconv::deconvolute(gene_exp, method = "timer", indications = indications, tumor = TRUE)

#Run Estimate
estimate <- immunedeconv::deconvolute(gene_exp, "estimate")

#EPIC
epic <- immunedeconv::deconvolute(gene_exp, "epic")

#MCP_COUNTER
mcp_counter <- immunedeconv::deconvolute(gene_exp, "mcp_counter")

#CIBERSORT
set_cibersort_binary("CIBERSORT.R")
set_cibersort_mat("LM22.txt")
cibersort <- cibersort("C:/Users/USER/Downloads/LM22.txt", gene_exp, QN = FALSE)

#cibersort_abs <- immunedeconv::deconvolute(gene_exp, "cibersort_abs")
#quantiseq <- immunedeconv::deconvolute(gene_exp, "quantiseq")
#xcell <- immunedeconv::deconvolute(gene_exp, "xcell")
#abis <- immunedeconv::deconvolute(gene_exp, "abis")
#consensus_tme <- immunedeconv::deconvolute(gene_exp, "consensus_tme", indications = indications)
