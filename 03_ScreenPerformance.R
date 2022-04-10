setwd("~/ResubPrep/")
source("CodeOrg/00_Utilities.R")
#"The output directory you would like to put"
outdir <- "~/Desktop"
Dat <- Sys.Date()
reCal <- F
if(reCal){
  minilib = DiagInit(LFC=readRDS("Processed_Data/Minilib_JointLFC.rds"),
                     posctrl_gene=readRDS("Raw_Data/Ctrls/Minilib_gene_pos.rds"),
                     negctrl_gene=readRDS("Raw_Data/Ctrls/Minilib_gene_neg.rds"),
                     posctrl_pair = readRDS("Raw_Data/Ctrls/Minilib_pair_pos.rds"),
                     annot = readRDS("Raw_Data/Minilib_Annotation.rds")%>%
                       select(Cas9_oligo,symbol1,symbol2),
                     cl = c("IPC298_SKIN","MELJUSO_SKIN","PK1_PANCREAS"))
  saveRDS(minilib,"Processed_Data/minilib_Init.rds")
  SpSaLib = DiagInit(LFC=readRDS("Raw_Data/SpSa_mini3LFC.rds"),
                     posctrl_gene=readRDS("Raw_Data/Ctrls/SpSa_gene_pos.rds"),
                     negctrl_gene=readRDS("Raw_Data/Ctrls/SpSa_gene_neg.rds"),
                     posctrl_pair = readRDS("Raw_Data/Ctrls/SpSa_pair_pos.rds"),
                     annot = readRDS("Raw_Data/SpSa_Annotation.rds"),
                     cl = c("IPC298_SKIN","MELJUSO_SKIN","PK1_PANCREAS"))
  saveRDS(SpSaLib,"Processed_Data/SpSaLib_Init.rds")
}