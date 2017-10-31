library(limma)
library(DESeq2)

#import files and set them up correctly
#key <- as.matrix(read.csv("KeyTable.csv", stringsAsFactors=FALSE, header = FALSE))
#FastQAndTotalSEQ <- as.data.frame(read.csv("Combined.csv", stringsAsFactors=FALSE, header = TRUE))
#rownames(FastQAndTotalSEQ)<-FastQAndTotalSEQ[,1]
#FastQAndTotalSEQ <- FastQAndTotalSEQ[,-1]

#take the samples from the table combined which were done with FastQ and TotalRNASeq
#TotalSEQ <- FastQAndTotalSEQ[,(names(FastQAndTotalSEQ) %in% key[,1])]
#FastQ <- FastQAndTotalSEQ[,(names(FastQAndTotalSEQ) %in% key[,2])]

#######################################################################################################################
#######################################################################################################################
########################Plot two samples against each other with pre-done analysis#####################################
#######################################################################################################################
#######################################################################################################################

#read in table 
QuantSEQ <- as.data.frame(read.table("Input/jack2.genes_de.tab", sep="\t",header=TRUE))
TotalSEQ <- as.data.frame(read.table("Input/deseq_Nicol_KO_differential_expression.tab", sep="\t",header=TRUE))

#Function to make zscore from Jack
make_signedZscore <- function( log2FC, pvalue){
  signedZ <- ifelse(test = log2FC > 0,
                    yes = qnorm(1 - (pvalue / 2)),
                    no = qnorm(pvalue / 2) )
  # high pass - otherwise very low P is set to Inf
  
  signedZ[ signedZ > 10 ] <- 10
  signedZ[ signedZ < -10 ] <- -10
  return(signedZ)
}

#make Zscores and rename matrixes
QuantZ <- data.frame(QuantSEQ$gene_id, make_signedZscore(QuantSEQ$logFC, QuantSEQ$PValue))
TotalZ <- data.frame(TotalSEQ$EnsemblID, make_signedZscore(TotalSEQ$log2FoldChange, TotalSEQ$pvalue))
colnames(TotalZ) <- c("GeneID", "TotalZScore")
colnames(QuantZ) <- c("GeneID", "FastQZScore")

#merge and reformat both sets of z scores into a single dataframe
Zscores <- merge(TotalZ, QuantZ, by="GeneID")
row.names(Zscores) <- Zscores$GeneID
Zscores <- Zscores[,-1]

#plot zscores against each other and add straight lines at 2 and -1 on the x and y 
plot (Zscores[complete.cases(Zscores),])
abline(h = c(2,-2))
abline(v = c(2,-2))

cor.test(Zscores[,1], Zscores[,2])

#######################################################################################################################
#######################################################################################################################
###################################################Run DeSeq###########################################################
#######################################################################################################################
#######################################################################################################################

#import raw reads from QuantSeq and format resulding dataframe
QuantRaw <- as.data.frame(read.table("Input/jack2.expression_genes.tab", sep="\t",header=TRUE, stringsAsFactors = FALSE ))
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
print(head(deseq.res.df)) 

#Create Z scores and merge both tables
QuantZ <- data.frame(rownames(deseq.res.df), make_signedZscore(deseq.res.df$log2FoldChange, deseq.res.df$padj))
TotalZ <- data.frame(TotalSEQ$EnsemblID, make_signedZscore(TotalSEQ$log2FoldChange, TotalSEQ$pvalue))
colnames(TotalZ) <- c("GeneID", "TotalZScore")
colnames(QuantZ) <- c("GeneID", "QuantZScore")
Zscores <- merge(TotalZ, QuantZ, by="GeneID")
row.names(Zscores) <- Zscores$GeneID
Zscores <- Zscores[,-1]

#plot z scores using only complete rows
plot (Zscores[complete.cases(Zscores),], xlab = "TotalRNASeq", ylab="QuantSeq")
abline(h = c(2,-2))
abline(v = c(2,-2))
