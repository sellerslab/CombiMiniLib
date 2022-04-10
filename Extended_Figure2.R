setwd("~/ResubPrep/")
source("CodeOrg/00_Utilities.R")
#"The output directory you would like to put"
outdir <- "~/Desktop"
Dat <- Sys.Date()
reCal <- F
#===== Replicate Correlation (Extended Figure 2a)  =====
LFCs_repsLevel <- readRDS("Processed_Data/Minilib_JointLFC_Replicate.rds")
corLFCs <- LFCs_repsLevel%>%
  select(-guidepair)%>%
  cor(.)
heatmaply::heatmaply_cor(corLFCs,limits = c(round(min(corLFCs),1)-0.05, 1),
                         file=glue("{outdir}/MiniLibScreen_repCor.pdf"),
                         width=1400,height=1100)
#===== Align with Avana v21q4  (Extended Figure 2b) =====
minilibInt <- readRDS("Processed_Data/minilib_Init.rds")
saspInt <- readRDS("Processed_Data/SpSaLib_Init.rds")
BalanceD_renamed_gene <- minilibInt
BalanceD_renamed_gene$avg_Single%<>%.[,grep("IPC298",colnames(BalanceD_renamed_gene$avg_Single),value=T)]
colnames(BalanceD_renamed_gene$avg_Single) <- gsub("IPC298\\_","",colnames(BalanceD_renamed_gene$avg_Single))
`[` <- function(...) base::`[`(...,drop=FALSE) 
saspInt$avg_Single%<>%.[,grep("IPC298",colnames(saspInt$avg_Single),value=T)]
colnames(saspInt$avg_Single) <- gsub("IPC298\\_","",colnames(saspInt$avg_Single))
plotOrder <- c("spCas9-saCas9", "enCas12a","VCR1-WCR3","WCR3-VCR1","WCR2-WCR3",
               "SPCR1-WCR3","SPCR1-SCR27","SPCR1-SCR43","VCR1-SCR27","VCR1-SCR43")
align_avana(inputlist = list(BalanceD_renamed_gene,saspInt),rowno=2,cor.ordered=F,orders=plotOrder,clSingle="IPC298")
ggsave(filename = paste0(outdir,"singleKO_Avana_combined_",Sys.Date(),".pdf"),width=13,height=5.8)