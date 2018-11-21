# library(DRIMSeq) 
# library(stageR)
# library(dplyr)
# library(tidyr)
# library(plyr)
# 



#make feature ID unique and reformat the data
featid <- function(df){
  df$feature_id <- paste(df$gene_name, df$site_pos, sep = "_")
  dat <- data.frame(gene_id = df$gene_id, feature_id = df$feature_id)
  # create count matrix with meta data columns
  dat[,3:22] <- df[,c(10:19, 21:30)]
  return(dat)
}


# run DRIMSeq
DRIMrun <- function(gr, df){
  d <- DRIMSeq::dmDSdata(counts = df, samples = gr) 
  d <- DRIMSeq::dmFilter(d, min_samps_gene_expr = 5, min_samps_feature_expr = 5,
                min_gene_expr = 10, min_feature_expr = 4)
  set.seed(1337)
  design_full <- model.matrix(~ group, data = d@samples)
  d <- DRIMSeq::dmPrecision(d, design = design_full)
  d <- DRIMSeq::dmFit(d, design = design_full, verbose = 1) 
  d <- DRIMSeq::dmTest(d, coef =  colnames(design_full)[2]) 
  # return d?
  return(d)
}


# for a single gene, plot the proportions of polyA sites
plotPolyA <- function(gene, data, DRIMSeq_res){
  # allow for Ensembl or regular gene names
  if( grepl("^ENS", gene)){
    ensemblID <- gene
    gene_name <- unique(filter(data, gene_id == ensemblID)$gene_name)
  }else{
    gene_name <- gene
    ensemblID <- unique(filter(data, gene_name == gene)$gene_id)
  }
  # get feature IDs for x axis
  gene_data <- dplyr::filter(data, data$gene_id == ensemblID)
  # extract per gene p value to add to plot title
  gene_res <- DRIMSeq::results(DRIMSeq_res, level = "gene") 
  gene_pvalue <- dplyr::filter(gene_res, gene_id == ensemblID)$pvalue

  if (gene_data$strand[1] == '+'){
    featid <- sort(gene_data$feature_id, decreasing = FALSE)
    sitepos <- sort(gene_data$site_pos, decreasing = FALSE)
  } else {
    featid <- sort(gene_data$feature_id, decreasing = TRUE)
    sitepos <- sort(gene_data$site_pos, decreasing = TRUE)
  }
  
  gene_pvalue <- signif(gene_pvalue, digits = 3)
  
  # make plot
  p <- DRIMSeq::plotProportions(DRIMSeq_res, gene_id = ensemblID, group_variable = "group")
  
  # add labels and title
  p <- p + scale_x_discrete(limits = featid) + xlab("") + labs(title = gene_name, subtitle = paste0("Gene P = ", gene_pvalue) )
  
  return(p)
}

# 
StageRun <- function(d, case_condition, alpha = 0.05){
  
  # get gene wide p values
  pScreen <- DRIMSeq::results(d, level = "gene")$pvalue
  
  names(pScreen) <- DRIMSeq::results(d)$gene_id
  
  # get polyA site specific p values
  pConfirmation <- matrix( DRIMSeq::results(d, level = "feature")$pvalue, ncol = 1)
  
  rownames(pConfirmation) <- DRIMSeq::results(d, level = "feature")$feature_id
  
  # map genes to features
  tx2gene <- DRIMSeq::results(d, level = "feature")[, c("feature_id", "gene_id")]
  
  # load into stageR object
  stageRObj <- stageR::stageRTx(pScreen = pScreen,
                                pConfirmation = pConfirmation,
                                pScreenAdjusted = FALSE, 
                                tx2gene = tx2gene)
  
  #perform stage-wise adjustment for differential transcript usage
  # alpha? probably should set
  stageRObj <- stageR::stageWiseAdjustment(object = stageRObj, method = "dtu",
                                   alpha = alpha)
  
  # adjust p value
  # again - all genes?
  padj <- stageR::getAdjustedPValues(stageRObj, order = TRUE,
                             onlySignificantGenes = FALSE)
  
  # 
  coeff <- DRIMSeq::coefficients(d, level="feature")
  
  # change is effect of the case condition 
  case_column <- paste0("group", case_condition)
  coe <- data.frame(txID = coeff$feature_id, change = coeff[[case_column]])
  
  results <- inner_join( padj, coe, by = "txID")
  
  # full is created from inner_join on multiple datasets - probably just have a general list of genes instead
  #full2 <- full 
  #full2 <- rename(full2, "txID" = "feature_id")
  
  # add in padj and r
  #TranscriptWPval <- merge(full2, padj, by="txID")
  
  #TranscriptWAll <- merge(TranscriptWPval, coe, by = "txID")
  
  #Temp <- dplyr::select(TranscriptWAll, gene_id, gene_name, txID, chr, strand, site_pos:change) 
  
  return( results )
  
  #return(Temp)
}
# 
# GetSites <- function(agR){
#   
#   AgaResults <- dplyr::select(agR, gene_id, gene_name, txID, gene, transcript, change)
#   AgaResults$sign <- sign(AgaResults$change)
#   
#   # unnecessary?
#   AgaResults$change <- AgaResults$change
#   
#   genes <- unique(agR$gene_name, incomparables = FALSE)
#   
#   # initialise empty data frames
#   Temp2 <- data.frame(txID=NA, PtD=NA)
#   MostDistal <- data.frame(txID=NA, PtD=NA)
#   
#   # for each gene
#   for(gene in genes){
#     # select all sites for this gene
#     tmp <- agR[agR$gene_name==gene,]
#     
#     if(nrow(tmp)==1){
#       next
#     }
#     # arrange sites by site position
#     if(tmp$strand == "-"){
#       tmp <- dplyr::arrange(tmp, desc(site_pos))
#     }
#     else{
#       tmp <- dplyr::arrange(tmp, site_pos)
#     }
#     # assign number based on strand-specific site position
#     tmp$PtD <- 1:nrow(tmp)
#     
#     tmp <- dplyr::select(tmp, txID, PtD)
#     # add site positions to big list
#     Temp2 <- rbind(Temp2, tmp)
#     # grab the most distal site and add to list of most distal sites per gene
#     MostDistal <- rbind(MostDistal, tmp[nrow(tmp),])
#   }
#   # remove first row?
#   Temp2 <- Temp2[-1,]
#   MostDistal <- MostDistal[-1,]
#   
#   ResultsFull <-  dplyr::full_join(AgaResults, Temp2)
#   MostDistalFull <-  dplyr::right_join(AgaResults, MostDistal)
#   
#   # get all the most proximal sites per gene
#   ProximalFull <-  dplyr::filter(ResultsFull, PtD == 1)
#   # get all the following sites
#   DistalFull <-  dplyr::filter(ResultsFull, PtD != 1)
#   
#   return(list(ResultsFull, MostDistalFull, ProximalFull, DistalFull))
# }
