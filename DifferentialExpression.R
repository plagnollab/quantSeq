library(limma)
library(DESeq2)

key <- as.matrix(read.csv("KeyTable.csv", stringsAsFactors=FALSE, header = FALSE))
FastQAndTotalSEQ <- as.data.frame(read.csv("Combined.csv", stringsAsFactors=FALSE, header = TRUE))
rownames(FastQAndTotalSEQ)<-FastQAndTotalSEQ[,1]
FastQAndTotalSEQ <- FastQAndTotalSEQ[,-1]

TotalSEQ <- FastQAndTotalSEQ[,(names(FastQAndTotalSEQ) %in% key[,1])]
FastQ <- FastQAndTotalSEQ[,(names(FastQAndTotalSEQ) %in% key[,2])]

Jack2 <- as.data.frame(read.table("jack2.genes_de.tab", sep="\t",header=TRUE))
Nicol <- as.data.frame(read.table("deseq_Nicol_KO_differential_expression.tab", sep="\t",header=TRUE))

make_signedZscore <- function( log2FC, pvalue){
  signedZ <- ifelse(test = log2FC > 0,
                    yes = qnorm(1 - (pvalue / 2)),
                    no = qnorm(pvalue / 2) )
  # high pass - otherwise very low P is set to Inf
  
  signedZ[ signedZ > 10 ] <- 10
  signedZ[ signedZ < -10 ] <- -10
  return(signedZ)
}

JackZ <- data.frame(Jack2$gene_id, make_signedZscore(Jack2$logFC, Jack2$PValue))
NicolZ <- data.frame(Nicol$EnsemblID, make_signedZscore(Nicol$log2FoldChange, Nicol$pvalue))
colnames(NicolZ) <- c("GeneID", "NicolZScore")
colnames(JackZ) <- c("GeneID", "JackZScore")

Zscores <- merge(NicolZ, JackZ, by="GeneID")
row.names(Zscores) <- Zscores$GeneID
Zscores <- Zscores[,-1]


plot (Zscores[complete.cases(Zscores),])
abline(h = c(2,-2))
abline(v = c(2,-2))
Gr2 <- subset(Zscores, JackZScore >2 & NicolZScore >2)
JackGr2 <- subset(Zscores, JackZScore >2 & NicolZScore < 2)
NicolGr2 <- subset(Zscores, NicolZScore >2 & JackZScore < 2)

le2 <- Gr2 <- subset(Zscores, JackZScore < -2 & NicolZScore < -2)
JackLeNeg2 <- subset(Zscores, JackZScore < -2 & NicolZScore > -2 )
NicolLeNeg2 <- subset(Zscores, NicolZScore < -2 & JackZScore > -2 )

cor.test(Zscores[,1], Zscores[,2])



QuantSeq <- data.frame(Jack2$gene_id, Jack2[,10:17])
rownames(QuantSeq) <- QuantSeq$Jack2.gene_id
QuantSeq <- QuantSeq[,-1]
