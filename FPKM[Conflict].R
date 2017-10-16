library(pheatmap)
library(ggplot2)

FUS14 <- as.matrix(read.csv("FUS_d14_rpkms.csv", stringsAsFactors=FALSE, header=TRUE))
FUSKO <- as.matrix(read.csv("FUS_KO_rpkms.csv", stringsAsFactors=FALSE, header = TRUE))
colnames(FUS14)[2] <- "gene_id"
colnames(FUSKO)[2] <- "gene_id"

key <- as.matrix(read.csv("KeyTable.csv", stringsAsFactors=FALSE, header = FALSE))

Jack1 <- as.matrix(read.table("jack1.expression_genes.tab", sep="\t",header=TRUE))
Jack2 <- as.matrix(read.table("jack2.expression_genes.tab", sep="\t",header=TRUE))
Jack3 <- as.matrix(read.table("jack3.expression_genes.tab", sep="\t",header=TRUE))
Jack4 <- as.matrix(read.table("jack4.expression_genes.tab", sep="\t",header=TRUE))

#test <- Jack1[1,3]
#blep <- sapply(strsplit(test, ":"), "[", 2)

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


RPKM <- function(Matr){
  J1 <- Matr[,c(2:4,9:17)]
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
  J1 <- J1[,c(3,16:23)]
  colnames(J1)<- colnames(Matr)[c(4, 10:17)]
  return(J1)
}

JackFin1 <- RPKM(Jack1)
JackFin2 <- RPKM(Jack2)
JackFin3 <- RPKM(Jack3)
JackFin4 <- RPKM(Jack4)

comb <- merge(JackFin1, JackFin2, by="gene_id", all=TRUE)
comb <- merge(comb, JackFin3, by="gene_id", all=TRUE)
comb <- merge(comb, JackFin4, by="gene_id", all=TRUE)

comb1 <- merge(comb, FUS14[,2:14], by="gene_id")
comb1 <- merge(comb1, FUSKO[,2:14], by="gene_id")

genenames <- comb1[,1]
row.names(comb1) <- genenames
comb1 <- comb1[,-1]

comb2 <- as.data.frame(comb1)

comp <- comb2[,c(1,9)]

comp[,1] <- as.numeric(comp[,1])
comp[,2] <- as.numeric(comp[,2])

cor.test(comp$c1.KO_WT1.x, comp$c1.KO_WT1.y)

ggplot(comp)

colnames(comb2)
drop <- c("c1.d14_WT1.y", "c2.d14_WT2.y", "c3.d14_WT3.y", "c4.d14_WT4.y", "c1.KO_WT1.y",  "c2.KO_WT2.y",  "c3.KO_WT3.y",  "c4.KO_WT4.y")
comb3 <- comb2[,!(names(comb2) %in% drop)]

write.csv(comb3, file="Combined.csv")

key <- as.matrix(read.csv("KeyTable.csv", stringsAsFactors=FALSE, header = FALSE))
comb3 <- as.data.frame(read.csv("Combined.csv", stringsAsFactors=FALSE, header = TRUE))

plot(comp)

colnames(comb3)

plot(log10(comb3[,(names(comb3) %in% key[1,])]),  xlab="C1.KO_WT.x", ylab="QuantSeq")
title(main = "Log KO WT 1")
plot(log10(comb3[,(names(comb3) %in% key[2,])]))
title(Main = "Log KO WT 2")
plot(log10(comb3[,(names(comb3) %in% key[3,])]))
title(Main = "Log KO WT 3")
plot(log10(comb3[,(names(comb3) %in% key[4,])]))
title(Main = "Log KO WT 4")
plot(log10(comb3[,(names(comb3) %in% key[5,])]))
title(Main = "Log KO HET 1")
plot(log10(comb3[,(names(comb3) %in% key[6,])]))
title(Main = "Log KO HET 2")
plot(log10(comb3[,(names(comb3) %in% key[7,])]))
title(Main = "Log KO HET 3")
plot(log10(comb3[,(names(comb3) %in% key[8,])]))
title(Main = "Log KO HET 4")
plot(log10(comb3[,(names(comb3) %in% key[9,])]))
title(Main = "Log KO HOM 1")
plot(log10(comb3[,(names(comb3) %in% key[10,])]))
title(Main = "Log KO HOM 2")
plot(log10(comb3[,(names(comb3) %in% key[11,])]))
title(Main = "Log KO HOM 3")
plot(log10(comb3[,(names(comb3) %in% key[12,])]))
title(Main = "Log KO HOM 4")
plot(log10(comb3[,(names(comb3) %in% key[13,])]))
title(Main = "Log d14 WT 1")
plot(log10(comb3[,(names(comb3) %in% key[14,])]))
title(Main = "Log d14 WT 2")
plot(log10(comb3[,(names(comb3) %in% key[15,])]))
title(Main = "Log d14 WT 3")
plot(log10(comb3[,(names(comb3) %in% key[16,])]))
title(Main = "Log d14 WT 4")
plot(log10(comb3[,(names(comb3) %in% key[17,])]))
title(Main = "Log d14 HET 1")
plot(log10(comb3[,(names(comb3) %in% key[18,])]))
title(Main = "Log d14 HET 2")
plot(log10(comb3[,(names(comb3) %in% key[19,])]))
title(Main = "Log d14 HET 3")
plot(log10(comb3[,(names(comb3) %in% key[20,])]))
title(Main = "Log d14 HET 4")
plot(log10(comb3[,(names(comb3) %in% key[21,])]))
title(Main = "Log d14 HOM 1")
plot(log10(comb3[,(names(comb3) %in% key[22,])]))
title(Main = "Log d14 HOM 2")
plot(log10(comb3[,(names(comb3) %in% key[23,])]))
title(Main = "Log d14 HOM 3")
plot(log10(comb3[,(names(comb3) %in% key[24,])]))
title(Main = "Log d14 HOM 4")

plot(comb3[,(names(comb3) %in% key[1,])],  xlab="C1.KO_WT.x", ylab="QuantSeq")
title(main = " KO WT 1")
plot(comb3[,(names(comb3) %in% key[2,])])
title(Main = " KO WT 2")
plot(comb3[,(names(comb3) %in% key[3,])])
title(Main = " KO WT 3")
plot(comb3[,(names(comb3) %in% key[4,])])
title(Main = " KO WT 4")
plot(comb3[,(names(comb3) %in% key[5,])])
title(Main = " KO HET 1")
plot(comb3[,(names(comb3) %in% key[6,])])
title(Main = " KO HET 2")
plot(comb3[,(names(comb3) %in% key[7,])])
title(Main = " KO HET 3")
plot(comb3[,(names(comb3) %in% key[8,])])
title(Main = " KO HET 4")
plot(comb3[,(names(comb3) %in% key[9,])])
title(Main = " KO HOM 1")
plot(comb3[,(names(comb3) %in% key[10,])])
title(Main = " KO HOM 2")
plot(comb3[,(names(comb3) %in% key[11,])])
title(Main = " KO HOM 3")
plot(comb3[,(names(comb3) %in% key[12,])])
title(Main = " KO HOM 4")
plot(comb3[,(names(comb3) %in% key[13,])])
title(Main = " d14 WT 1")
plot(comb3[,(names(comb3) %in% key[14,])])
title(Main = " d14 WT 2")
plot(comb3[,(names(comb3) %in% key[15,])])
title(Main = " d14 WT 3")
plot(comb3[,(names(comb3) %in% key[16,])])
title(Main = " d14 WT 4")
plot(comb3[,(names(comb3) %in% key[17,])])
title(Main = " d14 HET 1")
plot(comb3[,(names(comb3) %in% key[18,])])
title(Main = " d14 HET 2")
plot(comb3[,(names(comb3) %in% key[19,])])
title(Main = " d14 HET 3")
plot(comb3[,(names(comb3) %in% key[20,])])
title(Main = " d14 HET 4")
plot(comb3[,(names(comb3) %in% key[21,])])
title(Main = " d14 HOM 1")
plot(comb3[,(names(comb3) %in% key[22,])])
title(Main = " d14 HOM 2")
plot(comb3[,(names(comb3) %in% key[23,])])
title(Main = " d14 HOM 3")
plot(comb3[,(names(comb3) %in% key[24,])])
title(Main = " d14 HOM 4")


key[1,]
names(comb3)

cor.test(comb3$c1.KO_WT1.x, comb3$k1WT)
