library(DESeq2)
library(ggplot2)
library(gridExtra)


make_signedZscore <- function( log2FC, pvalue){
  signedZ <- ifelse(test = log2FC > 0,
                    yes = qnorm(1 - (pvalue / 2)),
                    no = qnorm(pvalue / 2) )
  # high pass - otherwise very low P is set to Inf
  
  signedZ[ signedZ > 6 ] <- 6
  signedZ[ signedZ < -6 ] <- -6
  return(signedZ)
}

produceZscoreTable <- function(QuantSEQ, TotalSEQ){
  QuantZ <- data.frame(row.names(QuantSEQ), make_signedZscore(QuantSEQ$log2FoldChange, QuantSEQ$pvalue))
  TotalZ <- data.frame(TotalSEQ$EnsemblID, make_signedZscore(TotalSEQ$log2FoldChange, TotalSEQ$pvalue))
  colnames(TotalZ) <- c("GeneID", "TotalZScore")
  colnames(QuantZ) <- c("GeneID", "FastQZScore")
  
  #merge and reformat both sets of z scores into a single dataframe
  Zscores <- merge(TotalZ, QuantZ, by="GeneID")
  row.names(Zscores) <- Zscores$GeneID
  Zscores <- Zscores[,-1]
  return(Zscores)
}

RunDEseq <- function (QuantSeq){
  coldata<-data.frame(sample = colnames(QuantSeq))
  coldata$Condition[1:4] <- "WT"
  coldata$Condition[5:8] <- "HOM"
  row.names(coldata) <- coldata$sample
  coldata$sample <- NULL
  coldata$Condition <- factor(coldata$Condition,levels=c("WT","HOM"))
  
  #run DeSeq on QuantSeq
  formula0 = as.formula("~ 1") 
  CDS <- DESeqDataSetFromMatrix(countData = QuantSeq, colData = coldata, design = as.formula("~ Condition"))
  CDS <- DESeq(CDS, test = "LRT", reduced = formula0, minReplicatesForReplace = 5 ) 
  deseq.res <- results(CDS)
  return(deseq.res)
}

d14TotalSEQ <- as.data.frame(read.table("deseq_Nicol_d14_differential_expression.tab", sep="\t",header=TRUE))
d14QuantRaw <- as.data.frame(read.table("jack4.expression_genes.tab", sep="\t",header=TRUE, stringsAsFactors = FALSE ))
d14QuantSeq <- data.frame(d14QuantRaw[,4], d14QuantRaw[,10:17])
rownames(d14QuantSeq) <- d14QuantSeq[,1]
d14QuantSeq <- d14QuantSeq[,-1]

d14QuantSEQ <- RunDEseq(d14QuantSeq)
d14Zscore <- produceZscoreTable(d14QuantSEQ, d14TotalSEQ)


d14PipeQuant <- as.data.frame(read.table("deseq_nicol_fus_d14_differential_expression.tab", sep="\t",header=TRUE, stringsAsFactors = FALSE ))
rownames(d14PipeQuant) <- d14PipeQuant[,1]
d14PipeQuant <- d14PipeQuant[,-1]
d14PipeZscore <- produceZscoreTable(d14PipeQuant, d14TotalSEQ)

present <- row.names(d14PipeZscore)  %in% row.names(d14Zscore)

d14 <- ggplot(d14Zscore, aes(x=TotalZScore, y=FastQZScore)) + theme_bw() + labs(title="Gregor FUS ctrl vs HOM d14 \n         ZScore (DEseq)", x="TotalRNASeq", y="QuantSeq", legend=FALSE) + geom_vline(xintercept = c(-2,2)) + geom_hline(yintercept = c(-2,2)) + coord_cartesian(xlim = c(-6, 6), ylim = c(-6,6)) + geom_point(aes(colour = TRUE), alpha=.5) + scale_colour_manual(values ="red") + guides(colour=FALSE)
d14Pipe <- ggplot(d14PipeZscore, aes(x=TotalZScore, y=(-1*FastQZScore))) + theme_bw() + labs(title="Pipeline FUS ctrl vs HOM d14 \n         ZScore (DEseq)", x="TotalRNASeq", y="QuantSeq", legend=FALSE) + geom_vline(xintercept = c(-2,2)) + geom_hline(yintercept = c(-2,2)) + coord_cartesian(xlim = c(-6, 6), ylim = c(-6,6)) + geom_point(aes(colour = present), alpha=.5) + scale_colour_manual(labels = c("Found by Gregor", "Excluded by Gregor"), values = c("red", "blue"))

grid.arrange(d14, d14Pipe, ncol=2, widths=c(1.3,2))

