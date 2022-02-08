source("CodeOrg/Utilities.R")
outdir <- "~/Desktop/tmp"
Dat <- Sys.Date()
dirC(outdir)
#===== From Count Processed to LFC =====
pooled_counts <- readRDS("Processed_Data/Minilib_JointCounts.rds")
mapped_info <- readRDS("Raw_Data/Barcode2Sample_Map.rds")
combs <- gsub("\\_Rep[ABC]","",mapped_info$Collapsed_identifier)%>%unique()
mapping_D_cas9 <- readRDS("Raw_Data/Minilib_Annotation.rds")
LFCs_sample <- lapply(combs,function(c){
  counts_pgDNA <- 
    pooled_counts%>%
    filter(grepl(c,Collapsed_identifier,ignore.case = T))%>%
    mutate(cell.line = case_when(cell.line=="N/A"~"pDNA",TRUE~cell.line))%>%
    mutate(Sample.Name = paste(cell.line,Collapsed_identifier,sep = "_"))%>%
    select(-Collapsed_identifier,-cell.line)%>%
    distinct()%>%
    spread(Sample.Name,pooled_counts)%>%
    replace(is.na(.), 0)%>%
    right_join(mapping_D_cas9%>%select(Label,rowname=Cas9_oligo),by=c("pairs"="Label"))%>%
    select(-pairs)%>%
    column_to_rownames("rowname")%>%
    as.matrix()%>%
    replace(is.na(.), 0)
  pDNA_name <- grep("Rep",colnames(counts_pgDNA),value=T,invert = T)
  reps_name <- grep("Rep",colnames(counts_pgDNA),value=T)%>%sort()
  counts_pgDNA_sort <- counts_pgDNA[,c(pDNA_name,reps_name)]
  reps_annot <- data.frame(
    repName = colnames(counts_pgDNA_sort),
    sampleName = gsub("\\_Rep.*","",colnames(counts_pgDNA_sort)),
    TP = as.factor(c("ETP",rep("LTP",ncol(counts_pgDNA_sort)-1))))
  lfcs <- lfc_calc(counts=counts_pgDNA_sort,replicate.map=reps_annot)%>%
    replace(is.na(.), 0)%>%
    as.data.frame()%>%
    mutate(guidepair = rownames(counts_pgDNA_sort))
  return(lfcs)
})
LFCs_sample%<>%purrr::reduce(inner_join, by = "guidepair") 
LFCs_sample%<>%column_to_rownames("guidepair")
saveRDS(LFCs_sample,"Processed_Data/Minilib_JointLFC.rds")
#===== Replicate Correlation (LFC level) =====
LFCs_repsLevel <- lapply(combs,function(c){
  counts_pgDNA <- 
    pooled_counts%>%
    filter(grepl(c,Collapsed_identifier,ignore.case = T))%>%
    mutate(cell.line = case_when(cell.line=="N/A"~"pDNA",TRUE~cell.line))%>%
    mutate(Sample.Name = paste(cell.line,Collapsed_identifier,sep = "_"))%>%
    select(-Collapsed_identifier,-cell.line)%>%
    distinct()%>%
    spread(Sample.Name,pooled_counts)%>%
    replace(is.na(.), 0)%>%
    right_join(mapping_D_cas9%>%select(Label,rowname=Cas9_oligo),by=c("pairs"="Label"))%>%
    select(-pairs)%>%
    column_to_rownames("rowname")%>%
    as.matrix()%>%
    replace(is.na(.), 0)
  pDNA_name <- grep("Rep",colnames(counts_pgDNA),value=T,invert = T)
  reps_name <- grep("Rep",colnames(counts_pgDNA),value=T)%>%sort()
  counts_pgDNA_sort <- counts_pgDNA[,c(pDNA_name,reps_name)]
  reps_annot <- data.frame(
    repName = colnames(counts_pgDNA_sort),
    sampleName = colnames(counts_pgDNA_sort),
    TP = as.factor(c("ETP",rep("LTP",ncol(counts_pgDNA_sort)-1))))
  lfcs <- lfc_calc(counts=counts_pgDNA_sort,replicate.map=reps_annot)%>%
    replace(is.na(.), 0)%>%
    as.data.frame()%>%
    mutate(guidepair = rownames(counts_pgDNA_sort))
  return(lfcs)
})
LFCs_repsLevel%<>%purrr::reduce(inner_join, by = "guidepair") 
saveRDS(LFCs_repsLevel,"Processed_Data/Minilib_JointLFC_Replicate.rds")
corLFCs <- LFCs_repsLevel%>%
  select(-guidepair)%>%
  cor(.)
heatmaply::heatmaply_cor(corLFCs,limits = c(round(min(corLFCs),1)-0.05, 1),
                         file=glue("{outdir}/MiniLibScreen_repCor.pdf"),
                         width=1400,height=1100)

