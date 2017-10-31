library(limma)
library(DESeq2)
library(ggplot2)
library(gridExtra)

#function to make z score
make_signedZscore <- function( log2FC, pvalue){
  signedZ <- ifelse(test = log2FC > 0,
                    yes = qnorm(1 - (pvalue / 2)),
                    no = qnorm(pvalue / 2) )
  # high pass - otherwise very low P is set to Inf
  
  signedZ[ signedZ > 6 ] <- 6
  signedZ[ signedZ < -6 ] <- -6
  return(signedZ)
}

#function to produce table of zscores of both Quant Seq and Total RNA Seq
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

#function to run deseq
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

#importand format d14 data
d14TotalSEQ <- as.data.frame(read.table("Input/deseq_Nicol_d14_differential_expression.tab", sep="\t",header=TRUE))
d14QuantRaw <- as.data.frame(read.table("Input/jack4.expression_genes.tab", sep="\t",header=TRUE, stringsAsFactors = FALSE ))
d14QuantSeq <- data.frame(d14QuantRaw[,4], d14QuantRaw[,10:17])
rownames(d14QuantSeq) <- d14QuantSeq[,1]
d14QuantSeq <- d14QuantSeq[,-1]

#Run Deseq and create zscore table of d14
d14QuantSEQ <- RunDEseq(d14QuantSeq)
d14Zscore <- produceZscoreTable(d14QuantSEQ, d14TotalSEQ)

#importand format d14 data
WTQuantRaw <- as.data.frame(read.table("Input/jack2.expression_genes.tab", sep="\t",header=TRUE, stringsAsFactors = FALSE ))
WTTotalSEQ <- as.data.frame(read.table("Input/deseq_Nicol_KO_differential_expression.tab", sep="\t",header=TRUE))
WTQuantSeq <- data.frame(WTQuantRaw[,4], WTQuantRaw[,10:17])
rownames(WTQuantSeq) <- WTQuantSeq[,1]
WTQuantSeq <- WTQuantSeq[,-1]

#Run Deseq and create zscore table of WT
WTQuantSEQ <- RunDEseq(WTQuantSeq)
WTZscore <- produceZscoreTable(WTQuantSEQ, WTTotalSEQ)

#Create plots for d14 and KO and output using grid arrange
d14 <- ggplot(d14Zscore, aes(x=TotalZScore, y=FastQZScore)) + geom_point(alpha=.2) + theme_bw() + labs(title="FUS ctrl vs HOM d14 ZScore (DEseq)", x="TotalRNASeq", y="QuantSeq") + geom_vline(xintercept = c(-2,2)) + geom_hline(yintercept = c(-2,2)) + coord_cartesian(xlim = c(-6, 6), ylim = c(-6,6)) 
KO <- ggplot(WTZscore, aes(x=TotalZScore, y=FastQZScore)) + geom_point(alpha=.2) + theme_bw() + labs(title="FUS ctrl vs HOM KO ZScore (DEseq)", x="TotalRNASeq", y="QuantSeq") + geom_vline(xintercept = c(-2,2)) + geom_hline(yintercept = c(-2,2)) + coord_cartesian(xlim = c(-6, 6), ylim = c(-6,6)) 

grid.arrange(d14, KO, ncol=2)
