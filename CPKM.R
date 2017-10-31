library(pheatmap)
library(ggplot2)

#setwd("/Users/Jack/google_drive/FUS NICOL/Seth/")

#Import Total RNA SEQ data and rename columns
FUS14 <- as.matrix(read.csv("Input/FUS_d14_rpkms.csv", stringsAsFactors=FALSE, header=TRUE))
FUSKO <- as.matrix(read.csv("Input/FUS_KO_rpkms.csv", stringsAsFactors=FALSE, header = TRUE))
colnames(FUS14)[2] <- "gene_id"
colnames(FUSKO)[2] <- "gene_id"

#read in key table which shows whic sample matches which
key <- as.matrix(read.csv("Input/KeyTable.csv", stringsAsFactors=FALSE, header = FALSE))

#Import Quant Seq data
Jack1 <- as.matrix(read.table("Input/jack1.expression_genes.tab", sep="\t",header=TRUE))
Jack2 <- as.matrix(read.table("Input/jack2.expression_genes.tab", sep="\t",header=TRUE))
Jack3 <- as.matrix(read.table("Input/jack3.expression_genes.tab", sep="\t",header=TRUE))
Jack4 <- as.matrix(read.table("Input/jack4.expression_genes.tab", sep="\t",header=TRUE))

#create function which gets RPKM values
RPKM <- function(Matr){
  #take only useful columns from matrix
  J1 <- Matr[,c(2:4,9:17)]
  #Take just the actual useful information
  J1[,2] <- sapply(strsplit(J1[,2], ":"), "[", 2)
  #Split the string again and take the gene start
  J1GeneStart <-as.numeric(sapply(strsplit(J1[,2], "-"), "[", 1))
  #take the gene end point
  J1GeneEND <-as.numeric(sapply(strsplit(J1[,2], "-"), "[", 2))
  #Add gene start and end to matrix
  J1 <- cbind(J1, GeneStart = J1GeneStart)
  J1 <- cbind(J1, GeneEND = J1GeneEND)
  #calculate gene length
  J1len <- as.numeric(J1[,14])-as.numeric(J1[,13])
  J1 <- cbind(J1, GeneLength = J1len)
  #Work out FPKM score for each sample in each gene and bind it
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
  J1 <- J1[,-c(1:2, 4:14)]
  #Rename Columns
  namelist <- list()
  namelist[1] <- colnames(Matr)[4]
  namelist[2] <- "GeneLength"
  namelist[3:10] <- colnames(Matr)[10:17]
  colnames(J1)<- namelist
  return(J1)
}

#get RPKM scores for each Quant Seq Sample
JackFin1 <- RPKM(Jack1)
JackFin2 <- RPKM(Jack2)[,-2]
JackFin3 <- RPKM(Jack3)[,-2]
JackFin4 <- RPKM(Jack4)[,-2]

#merge the quant seq samples into a single dataframe and remove duplicate columns
comb <- merge(JackFin1, JackFin2, by="gene_id", all=TRUE)
comb <- merge(comb, JackFin3, by="gene_id", all=TRUE)
comb <- merge(comb, JackFin4, by="gene_id", all=TRUE)
comb <- merge(comb, FUS14[,2:14], by="gene_id")
comb <- merge(comb, FUSKO[,2:14], by="gene_id")
genenames <- comb[,1]
row.names(comb) <- genenames
comb <- comb[,-1]
comb <- as.data.frame(comb)
drop <- c("c1.d14_WT1.y", "c2.d14_WT2.y", "c3.d14_WT3.y", "c4.d14_WT4.y", "c1.KO_WT1.y",  "c2.KO_WT2.y",  "c3.KO_WT3.y",  "c4.KO_WT4.y")
comb1 <- comb[,!(names(comb) %in% drop)]

#Create empty matrix then populate it with numeric versions of the CPKM
comb3 <- matrix(, nrow = length(as.numeric(comb1[,2])), ncol = 48)
for(y in 1:49){
  comb3[,(y-1)] <- as.numeric(comb1[,y])*as.numeric(comb1[,1])
}
rownames(comb3) <- rownames(comb1)
colnames(comb3) <- colnames(comb1)[-1]
comb3 <- as.data.frame(comb3)

#plot the log of the CPKM scores in the Total RNA Seq against Total RNA Seq
pdf("Log Graphs CPKM.pdf")
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

#plot the non-log versions
pdf("Graphs CPKM.pdf")
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

#crete function to calculate correlation between the two types of analysis
correlation <- function(loc, mat, key){
  correlation<-cor.test(mat[,(names(mat) %in% key[loc,1])], mat[,(names(mat) %in% key[loc,2])])
  return(correlation)
}

core1<-correlation(1,comb3,key)
core2<-correlation(2,comb3,key)
core3<-correlation(3,comb3,key)
core4<-correlation(4,comb3,key)
core5<-correlation(5,comb3,key)
core6<-correlation(6,comb3,key)
core7<-correlation(7,comb3,key)
core8<-correlation(8,comb3,key)
core9<-correlation(9,comb3,key)
core10<-correlation(10,comb3,key)
core11<-correlation(11,comb3,key)
core12<-correlation(12,comb3,key)
core13<-correlation(13,comb3,key)
core14<-correlation(14,comb3,key)
core15<-correlation(15,comb3,key)
core16<-correlation(16,comb3,key)
core17<-correlation(17,comb3,key)
core18<-correlation(18,comb3,key)
core19<-correlation(19,comb3,key)
core20<-correlation(20,comb3,key)
core21<-correlation(21,comb3,key)
core22<-correlation(22,comb3,key)
core23<-correlation(23,comb3,key)
core24<-correlation(24,comb3,key)

#attempt to create ggplot versions of plot
ggplot(comb3, aes_string(x=key[1,1], y=key[1,2])) + geom_point(alpha=.2) + geom_smooth(method=lm) + geom_abline(slope=1, intercept = 0) +theme_bw() + labs(title="KO WT 1", x="TotalRNASeq", y="QuantSeq")
ggplot(comb3, aes_string(x=key[2,1], y=key[2,2])) + geom_point() + geom_smooth(method=lm)
ggplot(comb3, aes_string(x=key[3,1], y=key[3,2])) + geom_point()
ggplot(comb3, aes_string(x=key[4,1], y=key[4,2])) + geom_point()
ggplot(comb3, aes_string(x=key[5,1], y=key[5,2])) + geom_point()
ggplot(comb3, aes_string(x=key[6,1], y=key[6,2])) + geom_point()
ggplot(comb3, aes_string(x=key[7,1], y=key[7,2])) + geom_point()
ggplot(comb3, aes_string(x=key[8,1], y=key[8,2])) + geom_point()
ggplot(comb3, aes_string(x=key[9,1], y=key[9,2])) + geom_point()
ggplot(comb3, aes_string(x=key[10,1], y=key[10,2])) + geom_point()
ggplot(comb3, aes_string(x=key[11,1], y=key[11,2])) + geom_point()
ggplot(comb3, aes_string(x=key[12,1], y=key[12,2])) + geom_point()
ggplot(comb3, aes_string(x=key[13,1], y=key[13,2])) + geom_point()
ggplot(comb3, aes_string(x=key[14,1], y=key[14,2])) + geom_point()
ggplot(comb3, aes_string(x=key[15,1], y=key[15,2])) + geom_point()
ggplot(comb3, aes_string(x=key[16,1], y=key[16,2])) + geom_point()
ggplot(comb3, aes_string(x=key[17,1], y=key[17,2])) + geom_point()
ggplot(comb3, aes_string(x=key[18,1], y=key[18,2])) + geom_point()
ggplot(comb3, aes_string(x=key[19,1], y=key[19,2])) + geom_point()
ggplot(comb3, aes_string(x=key[20,1], y=key[20,2])) + geom_point()
ggplot(comb3, aes_string(x=key[21,1], y=key[21,2])) + geom_point()
ggplot(comb3, aes_string(x=key[22,1], y=key[22,2])) + geom_point()
ggplot(comb3, aes_string(x=key[23,1], y=key[23,2])) + geom_point()
ggplot(comb3, aes_string(x=key[24,1], y=key[24,2])) + geom_point()

#create heatmaps of the combination as as well as just in d14 and KO
pheatmap(comb3[complete.cases(comb3),], scale = "row", show_rownames = F)
pheatmap(comb3[complete.cases(comb3),(names(comb3) %in% key[1:12,])], scale="row", show_rownames = F)
pheatmap(comb3[complete.cases(comb3),(names(comb3) %in% key[13:24,])], scale="row", show_rownames = F)

#Create sample which only had complete rows
comb4 <- comb3[complete.cases(comb3),]

#Created heatmap just looking at total RNA Seq
totalRNA <- comb4[ , 25:48]
totalRNA_sd <- apply( totalRNA, MAR =1, FUN = sd)
totalRNA <- totalRNA[ totalRNA_sd > 0, ]
totalRNA <- totalRNA[ complete.cases( totalRNA),]
pheatmap( totalRNA, scale = "row", show_rownames = FALSE, main="Total RNA Seq")

#create heatmap just looking at Quant Seq
QuantSeq <- comb4[ , 1:12]
QuantSeq_sd <- apply( QuantSeq, MAR =1, FUN = sd)
QuantSeq <- QuantSeq[ QuantSeq_sd > 0, ]
QuantSeq <- QuantSeq[ complete.cases( QuantSeq),]
pheatmap( QuantSeq, scale = "row", show_rownames = FALSE, main="Quant Seq")

#created heatmap with everything
pheatmap(comb4[,], scale = "row", show_rownames = F, main = "All Samples")
# quantseq heatmap
pheatmap(comb4[,(names(comb4) %in% key[,1])], scale="row", show_rownames = F)
pheatmap(comb4[,(names(comb4) %in% key[,2])], scale="row", show_rownames = F)