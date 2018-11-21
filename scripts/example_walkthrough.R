library(dplyr)
library(ggplot2)

# functions

source("../R/DRIMSeq_functions.R")

data <- as.data.frame(read.table("../example/aga9.expression_sites.tab.gz", sep = "\t", stringsAsFactors=FALSE, header = TRUE) )

# create support file 
WT <- select(data, contains("FUS.WT"))
D14 <- select(data, contains("FUS.D14"))

support <- data.frame( 
  sample_id = c( names(WT), names(D14)),
  group = factor(c( rep("control", 5), rep("FUS_d14", 5) ), levels = c("control", "FUS_d14"))
  )


# DRIMSeq only wants a file with gene_id, feature_id and then per sample counts

data$feature_id <- paste(data$gene_name, data$site_pos, sep = "_")

sample_counts <- select( data, gene_id, feature_id, contains("FUS.WT"), contains("FUS.D14"))

# run drimseq
drimseq_res <- DRIMrun(gr = support, df = sample_counts )


res <- DRIMSeq::results(drimseq_res, level = "feature") %>% 
  arrange(pvalue)

table(res$adj_pvalue < 0.05)

# plot top gene
top_gene <- head(res,1)$gene_id

# plot it
pdf("../example/example_plot.pdf")
plotPolyA(geneid = top_gene, data = data, DRIMSeq_res = drimseq_res)
dev.off()

# run StageR
staged_res <- StageRun(drimseq_res,case_condition = "FUS_d14")

# save output
readr::write_tsv(res,path = "../example/example_results.txt")
save( drimseq_res, staged_res, file = "../example/example_results.Rdata")


