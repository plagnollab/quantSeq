library(limma)
library(DESeq2)

make_signedZscore <- function( log2FC, pvalue){
  signedZ <- ifelse(test = log2FC > 0,
                    yes = qnorm(1 - (pvalue / 2)),
                    no = qnorm(pvalue / 2) )
  # high pass - otherwise very low P is set to Inf
  
  signedZ[ signedZ > 10 ] <- 10
  signedZ[ signedZ < -10 ] <- -10
  return(signedZ)
}

produceZscoreTable <- function(QuantSEQ, TotalSEQ){
  QuantZ <- data.frame(QuantSEQ$gene_id, make_signedZscore(QuantSEQ$log2FoldChange, QuantSEQ$pvalue))
  TotalZ <- data.frame(TotalSEQ$EnsemblID, make_signedZscore(TotalSEQ$log2FoldChange, TotalSEQ$pvalue))
  colnames(TotalZ) <- c("GeneID", "TotalZScore")
  colnames(QuantZ) <- c("GeneID", "FastQZScore")
  
  #merge and reformat both sets of z scores into a single dataframe
  Zscores <- merge(TotalZ, QuantZ, by="GeneID")
  row.names(Zscores) <- Zscores$GeneID
  Zscores <- Zscores[,-1]
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
  return 
}

#######################################################################################################################
#######################################################################################################################
########################Plot two samples against each other with pre-done analysis#####################################
#######################################################################################################################
#######################################################################################################################

#read in table 
QuantSEQ <- as.data.frame(read.table("jack4.genes_de.tab", sep="\t",header=TRUE))
TotalSEQ <- as.data.frame(read.table("deseq_Nicol_d14_differential_expression.tab", sep="\t",header=TRUE))

Zscores <- produceZscoreTable(QuantSEQ, TotalSEQ)

#plot zscores against each other and add straight lines at 2 and -1 on the x and y 
plot (Zscores[complete.cases(Zscores),], xlab = "TotalRNASeq", ylab="QuantSeq", ylim = c(-6, 6), xlim=c(-6,6))
abline(h = c(2,-2), lty=3)
abline(v = c(2,-2), lty=3)

cor.test(Zscores[,1], Zscores[,2])

#######################################################################################################################
#######################################################################################################################
###################################################Run DeSeq###########################################################
#######################################################################################################################
#######################################################################################################################

#import raw reads from QuantSeq and format resulding dataframe
QuantRaw <- as.data.frame(read.table("jack4.expression_genes.tab", sep="\t",header=TRUE, stringsAsFactors = FALSE ))
QuantSeq <- data.frame(QuantRaw[,4], QuantRaw[,10:17])
rownames(QuantSeq) <- QuantSeq$QuantRaw...4.
QuantSeq <- QuantSeq[,-1]

#create reference table
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
deseq.res.df <- data.frame(deseq.res) 

#Create Z scores and merge both tables
Zscore <- produceZscoreTable()
QuantZ <- data.frame(rownames(deseq.res.df), make_signedZscore(deseq.res.df$log2FoldChange, deseq.res.df$pvalue))
TotalZ <- data.frame(TotalSEQ$EnsemblID, make_signedZscore(TotalSEQ$log2FoldChange, TotalSEQ$pvalue))
colnames(TotalZ) <- c("GeneID", "TotalRNAZScore")
colnames(QuantZ) <- c("GeneID", "QuantZScore")
Zscores <- merge(TotalZ, QuantZ, by="GeneID")
row.names(Zscores) <- Zscores$GeneID
Zscores <- Zscores[,-1]

#plot z scoresfi
grid.arrange(ncol=2)
plot (Zscores[complete.cases(Zscores),], xlab = "TotalRNASeq", ylab="QuantSeq", ylim = c(-6, 6), xlim=c(-6,6), main="FUS ctrl vs HOM d14 (DEseq)")
abline(h = c(2,-2), lty=3)
abline(v = c(2,-2), lty=3)

