library(pheatmap)
library(ggplot2)

#setwd("/Users/Jack/google_drive/FUS NICOL/Seth/")
#Import Total RNA Seq data
FUS14 <- as.matrix(read.csv("Input/FUS_d14_rpkms.csv", stringsAsFactors=FALSE, header=TRUE))
FUSKO <- as.matrix(read.csv("Input/FUS_KO_rpkms.csv", stringsAsFactors=FALSE, header = TRUE))
colnames(FUS14)[2] <- "gene_id"
colnames(FUSKO)[2] <- "gene_id"

#import key table
key <- as.matrix(read.csv("Input/KeyTable.csv", stringsAsFactors=FALSE, header = FALSE))

#import quantseq data
Jack1 <- as.matrix(read.table("Input/jack1.expression_genes.tab", sep="\t",header=TRUE))
Jack2 <- as.matrix(read.table("Input/jack2.expression_genes.tab", sep="\t",header=TRUE))
Jack3 <- as.matrix(read.table("Input/jack3.expression_genes.tab", sep="\t",header=TRUE))
Jack4 <- as.matrix(read.table("Input/jack4.expression_genes.tab", sep="\t",header=TRUE))

#test <- Jack1[1,3]
#blep <- sapply(strsplit(test, ":"), "[", 2)

#this was to test the making of the RPKM calculator
#J1 <- Jack1[,c(2:4,9:17)]
#J1[,2] <- sapply(strsplit(J1[,2], ":"), "[", 2)
#J1ST <-as.numeric(sapply(strsplit(blep, "-"), "[", 1))
#J1END <-sapply(strsplit(blep, "-"), "[", 2)
#J1GeneStart <-as.numeric(sapply(strsplit(J1[,2], "-"), "[", 1))
#J1GeneEND <-as.numeric(sapply(strsplit(J1[,2], "-"), "[", 2))
#J1 <- cbind(J1, GeneStart = J1GeneStart)
#J1 <- cbind(J1, GeneEND = J1GeneEND)
#J1len <- as.numeric(J1[,14])-as.numeric(J1[,13])
#J1 <- cbind(J1, GeneLength = J1len)
#J1s <- as.numeric(J1[,5]) /((as.numeric(J1[,15])/1000)*(sum(as.numeric(J1[,5]))/1000000))
#J1 <- cbind(J1, FPKM = J1s)
#J1s <- as.numeric(J1[,6]) /((as.numeric(J1[,15])/1000)*(sum(as.numeric(J1[,5]))/1000000))
#J1 <- cbind(J1, FPKM2 = J1s)

#create a function which works out RPKM score of each sample, CPKM has full comments
RPKM <- function(Matr){
  J1 <- Jack1[,c(2:4,9:17)]
  J1[,2] <- sapply(strsplit(J1[,2], ":"), "[", 2)
  J1GeneStart <-as.numeric(sapply(strsplit(J1[,2], "-"), "[", 1))
  J1GeneEND <-as.numeric(sapply(strsplit(J1[,2], "-"), "[", 2))
  J1 <- cbind(J1, GeneStart = J1GeneStart)
  J1 <- cbind(J1, GeneEND = J1GeneEND)
  J1len <- as.numeric(J1[,14])-as.numeric(J1[,13])
  J1 <- cbind(J1, GeneLength = J1len)
  J1s <- as.numeric(J1[,5]) /((as.numeric(J1[,15])/1000)*(sum(as.numeric(J1[,5]))/1000000))
  J1 <- cbind(J1, C1 = J1s)
  J1s <- as.numeric(J1[,6]) /((as.numeric(J1[,15])/1000)*(sum(as.numeric(J1[,5]))/1000000))
  J1 <- cbind(J1, C2 = J1s)
  J1s <- as.numeric(J1[,7]) /((as.numeric(J1[,15])/1000)*(sum(as.numeric(J1[,5]))/1000000))
  J1 <- cbind(J1, C3 = J1s)
  J1s <- as.numeric(J1[,8]) /((as.numeric(J1[,15])/1000)*(sum(as.numeric(J1[,5]))/1000000))
  J1 <- cbind(J1, C4 = J1s)
  J1s <- as.numeric(J1[,9]) /((as.numeric(J1[,15])/1000)*(sum(as.numeric(J1[,5]))/1000000))
  J1 <- cbind(J1, T1 = J1s)
  J1s <- as.numeric(J1[,10]) /((as.numeric(J1[,15])/1000)*(sum(as.numeric(J1[,5]))/1000000))
  J1 <- cbind(J1, T2 = J1s)
  J1s <- as.numeric(J1[,11]) /((as.numeric(J1[,15])/1000)*(sum(as.numeric(J1[,5]))/1000000))
  J1 <- cbind(J1, T3 = J1s)
  J1s <- as.numeric(J1[,12]) /((as.numeric(J1[,15])/1000)*(sum(as.numeric(J1[,5]))/1000000))
  J1 <- cbind(J1, T4 = J1s)
  J1 <- J1[,-2:4]
  #J1 <- J1[,c(3,16:23)]
  #colnames(J1)<- colnames(Matr)[c(4, 10:17)]
  return(J1)
}

#run RPKM on each Quant Seq sample
JackFin1 <- RPKM(Jack1)
JackFin2 <- RPKM(Jack2)
JackFin3 <- RPKM(Jack3)
JackFin4 <- RPKM(Jack4)

#merge all the quant seq samples into a single matix and rename columnsand rows
comb <- merge(JackFin1, JackFin2, by="gene_id", all=TRUE)
comb <- merge(comb, JackFin3, by="gene_id", all=TRUE)
comb <- merge(comb, JackFin4, by="gene_id", all=TRUE)
comb1 <- merge(comb, FUS14[,2:14], by="gene_id")
comb1 <- merge(comb1, FUSKO[,2:14], by="gene_id")
genenames <- comb1[,1]
row.names(comb1) <- genenames
comb1 <- comb1[,-1]

#convert into a dataframe
comb2 <- as.data.frame(comb1)

#take only useful columns
comp <- comb2[,c(1,9)]
comp[,1] <- as.numeric(comp[,1])
comp[,2] <- as.numeric(comp[,2])

#ggplot(comp)
#remove repeat columns
colnames(comb2)
drop <- c("c1.d14_WT1.y", "c2.d14_WT2.y", "c3.d14_WT3.y", "c4.d14_WT4.y", "c1.KO_WT1.y",  "c2.KO_WT2.y",  "c3.KO_WT3.y",  "c4.KO_WT4.y")
comb3 <- comb2[,!(names(comb2) %in% drop)]

#write combined data frame to computer to reduce memory loss
write.csv(comb3, file="Combined.csv")

#import key and reimport combined dataframe
key <- as.matrix(read.csv("KeyTable.csv", stringsAsFactors=FALSE, header = FALSE))
comb3 <- as.data.frame(read.csv("Combined.csv", stringsAsFactors=FALSE, header = TRUE))

#plot log of RPKM scores
pdf("Log Graphs.pdf")
plot(log10(comb3[,(names(comb3) %in% key[1,])]), ylab="QuantSeq", xlab="TotalRNASeq", main="Log KO WT 1")
plot(log10(comb3[,(names(comb3) %in% key[2,])]), ylab="QuantSeq", xlab="TotalRNASeq", main="Log KO WT 2")
plot(log10(comb3[,(names(comb3) %in% key[3,])]), ylab="QuantSeq", xlab="TotalRNASeq", main="Log KO WT 3")
plot(log10(comb3[,(names(comb3) %in% key[4,])]), ylab="QuantSeq", xlab="TotalRNASeq", main="Log KO WT 4")
plot(log10(comb3[,(names(comb3) %in% key[5,])]), ylab="QuantSeq", xlab="TotalRNASeq", main="Log KO HET 1")
plot(log10(comb3[,(names(comb3) %in% key[6,])]), ylab="QuantSeq", xlab="TotalRNASeq", main="Log KO HET 2")
plot(log10(comb3[,(names(comb3) %in% key[7,])]), ylab="QuantSeq", xlab="TotalRNASeq", main="Log KO HET 3")
plot(log10(comb3[,(names(comb3) %in% key[8,])]), ylab="QuantSeq", xlab="TotalRNASeq", main="Log KO HET 4")
plot(log10(comb3[,(names(comb3) %in% key[9,])]), ylab="QuantSeq", xlab="TotalRNASeq", main="Log KO HOM 1")
plot(log10(comb3[,(names(comb3) %in% key[10,])]), ylab="QuantSeq", xlab="TotalRNASeq", main="Log KO HOM 2")
plot(log10(comb3[,(names(comb3) %in% key[11,])]), ylab="QuantSeq", xlab="TotalRNASeq", main="Log KO HOM 3")
plot(log10(comb3[,(names(comb3) %in% key[12,])]), ylab="QuantSeq", xlab="TotalRNASeq", main="Log KO HOM 4")
plot(log10(comb3[,(names(comb3) %in% key[13,])]), ylab="QuantSeq", xlab="TotalRNASeq", main="Log d14 WT 1")
plot(log10(comb3[,(names(comb3) %in% key[14,])]), ylab="QuantSeq", xlab="TotalRNASeq", main="Log d14 WT 2")
plot(log10(comb3[,(names(comb3) %in% key[15,])]), ylab="QuantSeq", xlab="TotalRNASeq", main="Log d14 WT 3")
plot(log10(comb3[,(names(comb3) %in% key[16,])]), ylab="QuantSeq", xlab="TotalRNASeq", main="Log d14 WT 4")
plot(log10(comb3[,(names(comb3) %in% key[17,])]), ylab="QuantSeq", xlab="TotalRNASeq", main="Log d14 HET 1")
plot(log10(comb3[,(names(comb3) %in% key[18,])]), ylab="QuantSeq", xlab="TotalRNASeq", main="Log d14 HET 2")
plot(log10(comb3[,(names(comb3) %in% key[19,])]), ylab="QuantSeq", xlab="TotalRNASeq", main="Log d14 HET 3")
plot(log10(comb3[,(names(comb3) %in% key[20,])]), ylab="QuantSeq", xlab="TotalRNASeq", main="Log d14 HET 4")
plot(log10(comb3[,(names(comb3) %in% key[21,])]), ylab="QuantSeq", xlab="TotalRNASeq", main="Log d14 HOM 1")
plot(log10(comb3[,(names(comb3) %in% key[22,])]), ylab="QuantSeq", xlab="TotalRNASeq", main="Log d14 HOM 2")
plot(log10(comb3[,(names(comb3) %in% key[23,])]), ylab="QuantSeq", xlab="TotalRNASeq", main="Log d14 HOM 3")
plot(log10(comb3[,(names(comb3) %in% key[24,])]), ylab="QuantSeq", xlab="TotalRNASeq", main="Log d14 HOM 4")
dev.off()

#plot RPKM scores
pdf("Graphs.pdf")
plot(comb3[,(names(comb3) %in% key[1,])], ylab="QuantSeq", xlab="TotalRNASeq", main="KO WT 1")
plot(comb3[,(names(comb3) %in% key[2,])], ylab="QuantSeq", xlab="TotalRNASeq", main="KO WT 2")
plot(comb3[,(names(comb3) %in% key[3,])], ylab="QuantSeq", xlab="TotalRNASeq", main="KO WT 3")
plot(comb3[,(names(comb3) %in% key[4,])], ylab="QuantSeq", xlab="TotalRNASeq", main="KO WT 4")
plot(comb3[,(names(comb3) %in% key[5,])], ylab="QuantSeq", xlab="TotalRNASeq", main="KO HET 1")
plot(comb3[,(names(comb3) %in% key[6,])], ylab="QuantSeq", xlab="TotalRNASeq", main="KO HET 2")
plot(comb3[,(names(comb3) %in% key[7,])], ylab="QuantSeq", xlab="TotalRNASeq", main="KO HET 3")
plot(comb3[,(names(comb3) %in% key[8,])], ylab="QuantSeq", xlab="TotalRNASeq", main="KO HET 4")
plot(comb3[,(names(comb3) %in% key[9,])], ylab="QuantSeq", xlab="TotalRNASeq", main="KO HOM 1")
plot(comb3[,(names(comb3) %in% key[10,])], ylab="QuantSeq", xlab="TotalRNASeq", main="KO HOM 2")
plot(comb3[,(names(comb3) %in% key[11,])], ylab="QuantSeq", xlab="TotalRNASeq", main="KO HOM 3")
plot(comb3[,(names(comb3) %in% key[12,])], ylab="QuantSeq", xlab="TotalRNASeq", main="KO HOM 4")
plot(comb3[,(names(comb3) %in% key[13,])], ylab="QuantSeq", xlab="TotalRNASeq", main="d14 WT 1")
plot(comb3[,(names(comb3) %in% key[14,])], ylab="QuantSeq", xlab="TotalRNASeq", main="d14 WT 2")
plot(comb3[,(names(comb3) %in% key[15,])], ylab="QuantSeq", xlab="TotalRNASeq", main="d14 WT 3")
plot(comb3[,(names(comb3) %in% key[16,])], ylab="QuantSeq", xlab="TotalRNASeq", main="d14 WT 4")
plot(comb3[,(names(comb3) %in% key[17,])], ylab="QuantSeq", xlab="TotalRNASeq", main="d14 HET 1")
plot(comb3[,(names(comb3) %in% key[18,])], ylab="QuantSeq", xlab="TotalRNASeq", main="d14 HET 2")
plot(comb3[,(names(comb3) %in% key[19,])], ylab="QuantSeq", xlab="TotalRNASeq", main="d14 HET 3")
plot(comb3[,(names(comb3) %in% key[20,])], ylab="QuantSeq", xlab="TotalRNASeq", main="d14 HET 4")
plot(comb3[,(names(comb3) %in% key[21,])], ylab="QuantSeq", xlab="TotalRNASeq", main="d14 HOM 1")
plot(comb3[,(names(comb3) %in% key[22,])], ylab="QuantSeq", xlab="TotalRNASeq", main="d14 HOM 2")
plot(comb3[,(names(comb3) %in% key[23,])], ylab="QuantSeq", xlab="TotalRNASeq", main="d14 HOM 3")
plot(comb3[,(names(comb3) %in% key[24,])], ylab="QuantSeq", xlab="TotalRNASeq", main="d14 HOM 4")
dev.off()

#Remove ensembl ID column and turn it into row names
listy <- comb3$X
rownames(comb3) <- listy
comb3 <- comb3[,-1]

#create heatmaps
pheatmap(comb3, scale = "row", show_rownames = F)
pheatmap(comb3, scale = "column")
pheatmap(comb3[,(names(comb3) %in% key[1:12,])], scale="row")
pheatmap(comb3[,(names(comb3) %in% key[1:12,])], scale="column")
pheatmap(comb3[,(names(comb3) %in% key[13:24,])], scale="row")
pheatmap(comb3[,(names(comb3) %in% key[13:24,])], scale="column")
log10(comb3)
