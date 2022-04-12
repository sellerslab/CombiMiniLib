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
outdir <- "~/Desktop/ExtendedFigure2"
if(!dir.exists(outdir)){dir.create(outdir)}
Dat <- Sys.Date()
sessionInfo()
#===== Replicate Correlation (Extended Figure 2a)  =====
LFCs_repsLevel <- readRDS(glue("Processed_Data/Minilib_JointLFC_Replicate.rds"))
corLFCs <- LFCs_repsLevel%>%
  select(-guidepair)%>%
  cor(.)
heatmaply::heatmaply_cor(corLFCs,limits = c(round(min(corLFCs),1)-0.05, 1),
                         file=glue("{outdir}/ExFigure2a_MiniLibScreen_repCor.png"),
                         width=1400,height=1100)
#===== Align with Avana v21q4  (Extended Figure 2b) =====
minilibInt <- readRDS(glue("{Processed_Data}/minilib_Init.rds"))
saspInt <- readRDS(glue("Processed_Data/SpSaLib_Init.rds"))
BalanceD_renamed_gene <- minilibInt
BalanceD_renamed_gene$avg_Single%<>%.[,grep("IPC298",colnames(BalanceD_renamed_gene$avg_Single),value=T)]
colnames(BalanceD_renamed_gene$avg_Single) <- gsub("IPC298\\_","",colnames(BalanceD_renamed_gene$avg_Single))
`[` <- function(...) base::`[`(...,drop=FALSE) 
saspInt$avg_Single%<>%.[,grep("IPC298",colnames(saspInt$avg_Single),value=T)]
colnames(saspInt$avg_Single) <- gsub("IPC298\\_","",colnames(saspInt$avg_Single))
plotOrder <- c("spCas9-saCas9", "enCas12a","VCR1-WCR3","WCR3-VCR1","WCR2-WCR3",
               "SPCR1-WCR3","SPCR1-SCR27","SPCR1-SCR43","VCR1-SCR27","VCR1-SCR43")
s2b <- align_avana(inputlist = list(BalanceD_renamed_gene,saspInt),rowno=2,cor.ordered=F,orders=plotOrder,clSingle="IPC298",
                   dataLoad=T)
ggsave(s2b,filename = glue("{outdir}/ExFigure2b_singleKO_Avana_combined_{Dat}.pdf"),width=13,height=5.8)
#===== Source data output =====
write.xlsx(list("Panel a" = corLFCs,
                "Panel b" = s2b$data),
           glue("{outdir}/ExFigure2_sourceData_{Dat}.xlsx"),
           overwrite = T)