#===== Env Specification and library loading  ====
rm(list=ls())
# Set the wokring directory to a folder containing three folders
# Processed_data
# Raw_data
# CodeOrg
# If you have different names for those three folders and prefer not changing
# Please change the strings for all 3 variables specified below
setwd("~/ResubPrep/")
CodeOrg <- "CodeOrg"
Processed_Data <- "Processed_Data"
Raw_Data <- "Raw_Data"
source(glue::glue("{CodeOrg}/00_Utilities.R"))
#"The output directory you would like to put"
outdir <- "~/Desktop/SuppData"
if(!dir.exists(outdir)){dir.create(outdir)}
Dat <- Sys.Date()
sessionInfo()
#===== Supplementary Tables Output =====
Cas9_lib<- readRDS(glue("{Raw_Data}/sgRNASource_Cas9_enCas12a.rds"))$Cas9_sgRNA
enCas12a_lib  <- readRDS(glue("{Raw_Data}/sgRNASource_Cas9_enCas12a.rds"))$enCas12a_sgRNA
minilib <- readRDS(glue("{Processed_Data}/minilib_Init.rds"))
annot <- readRDS(glue("{Raw_Data}/Minilib_Annotation.rds"))%>%select(-essClass)
# no 15 copy genes+chr2+AAVS1+nonTarget 617-18=599
# 10 genes only have 3 sgRNAs, paired with AAVS1 to investigate promoter effect
# 589*6+10*3 = 3564 sgRNAs
# Suggest to remove chr2, AAVS1, nonTarget trunk since we never talked about them in the manuscript
write.xlsx(
  list(
    "Supp1(Cas9 Gene-sgRNA Map)" = Cas9_lib,
    "Supp2(Cas9 Construct Oligo)" = annot%>%select(-Cas12a_oligo),
    "Supp3(enCas12a Gene-sgRNA Map)" = enCas12a_lib,
    "Supp4(enCas12a Construct Oligo)" = annot%>%select(-Cas9_oligo),
    "Supp5(Screen Raw Count)" = readRDS(glue("{Processed_Data}/Minilib_JointCounts.rds"))%>%
      mutate(cell.line = case_when(cell.line=="N/A"~"pDNA",TRUE~cell.line))%>%
      mutate(Sample.Name = paste(cell.line,Collapsed_identifier,sep = "_"))%>%
      select(-Collapsed_identifier,-cell.line)%>%
      distinct()%>%
      spread(Sample.Name,pooled_counts)%>%
      replace(is.na(.), 0),
    "Supp6(Screen Raw LFC)" = readRDS(glue("{Processed_Data}/Minilib_JointLFC.rds"))%>%as.data.frame()%>%rownames_to_column("sgRNA_Pair"),
    "Supp7(Screen Sensitive Score)" = minilib$score%>%as.data.frame()%>%rownames_to_column("gene_Pair"),
    "Supp8(Gene-level pos-control)" = minilib$posctrl_gene,
    "Supp9(Gene-level neg-control)" = minilib$negctrl_gene,
    "Supp10(Pair-level pos-control)" = minilib$posctrl_pair,
    "Supp11(Pair-level neg-control)" = minilib$negctrl_pair),
  glue("{outdir}/Supplementary_Data_{Dat}.xlsx"),overwrite=T
)
