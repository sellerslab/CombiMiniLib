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
outdir <- "~/Desktop/Figure4"
if(!dir.exists(outdir)){dir.create(outdir)}
Dat <- Sys.Date()
sessionInfo()
#===== Pair Level Separation Analysis  (Figure 4a-c) =====
minilib <- readRDS(glue("{Processed_Data}/minilib_Init.rds"))
SpSaLib <- readRDS(glue("{Processed_Data}/SpSaLib_Init.rds"))
combLib <- list("minilib"=minilib,"SpSaLib"=SpSaLib)
system3_pairSep <- EvalWrapper(combLib,gene=F)
rocList <- lapply(c("IPC298","MELJUSO","PK1"),function(l){
  if(l=="IPC298"){
    SepBar(system3_pairSep$metrics%>%filter(grepl(l,cl))%>%mutate(cl = gsub(paste0(l,"\\_"),"",cl)),colorp="Paired")%>%
      ggsave(plot = .,filename = glue("{outdir}/UnionNeg_PairSeparation_Bar_{l}_{Dat}.pdf"),width=14,height=5)
  }
  labelMap <- system3_pairSep$metrics%>%
    filter(grepl(l,cl))%>%
    filter(metric =="AUC-ROC")%>%
    arrange(-values)%>%
    mutate(fakeB =  " (",fakeE = ")")%>%
    unite("Label",c("cl","fakeB","values","fakeE"),sep="",remove=F)%>%
    mutate(Label = gsub(paste0(l,"\\_"),"",Label))
  if(l !="IPC298"){
    cols <-c("#1F78B4","#33A02C","#E31A1C","#FB9A99")
  }else{cols <- NULL}
  rocP<- SepROC(system3_pairSep$confusion%>%
                  filter(grepl(l,sample))%>%
                  inner_join(labelMap%>%select(Label,cl),by=c("sample"="cl"))%>%
                  mutate(sample = Label),
                labelLevels=labelMap$Label,ptitle=paste0("Digenic paralog - ",l),
                manucol = cols)
  return(rocP)
})
ggarrange(plotlist = rocList,nrow=1)%>%ggsave(plot=.,glue("{outdir}/Figure4ac_UnionNeg_PairSeparation_ROC_{Dat}.pdf"),width = 16.8,height=6.2)
#===== Pair-Level left-right comparison (Figure 4d) =====
minilibInt <- readRDS(glue("{Processed_Data}/minilib_Init.rds"))
saspInt <- readRDS(glue("{Processed_Data}/SpSaLib_Init.rds"))
BalanceD_renamed_pair <- minilibInt
BalanceD_renamed_pair$LFC%<>%.[,grep("IPC298",colnames(minilibInt$LFC),value=T)]
colnames(BalanceD_renamed_pair$LFC) <- gsub("IPC298\\_","",colnames(BalanceD_renamed_pair$LFC))
`[` <- function(...) base::`[`(...,drop=FALSE) 
saspInt$LFC%<>%.[,grep("IPC298",colnames(saspInt$LFC),value=T)]
colnames(saspInt$LFC) <- gsub("IPC298\\_","",colnames(saspInt$LFC))
# plotOrder <- c("VCR1-WCR3","WCR2-WCR3",
#                "enCas12a","spCas9-saCas9",
#                "WCR3-VCR1","SPCR1-WCR3",
#                "VCR1-SCR27","VCR1-SCR43",
#                "SPCR1-SCR43","SPCR1-SCR27")
plotOrder <- c("VCR1-WCR3","WCR2-WCR3",
               "enCas12a","spCas9-saCas9")
double_lf_mean <- left_right_scatter(inputlist = list(BalanceD_renamed_pair,saspInt),
                                     type="double",FUN="mean",
                                     rowno=1,cor.ordered=F,orders=plotOrder)
percD <- double_lf_mean$data%>%
  filter(side1< -1 & side2 < -1 & Essentiality=="yes")%>%
  group_by(source)%>%
  summarise(n_hit=n())%>%
  ungroup%>%
  filter(source %in% c("VCR1-WCR3","WCR2-WCR3","enCas12a","spCas9-saCas9"))%>%
  mutate(source = factor(source,levels = c("VCR1-WCR3","WCR2-WCR3","enCas12a","spCas9-saCas9")))%>%
  arrange(source)%>%
  mutate(n_pos = c(rep(length(minilibInt$posctrl_pair),3),length(saspInt$posctrl_pair)))%>%
  mutate(phits = round(n_hit*100/n_pos,1))%>%
  mutate(phits = paste0(phits,"%"))
Fig4d <- double_lf_mean+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(hjust = 0,colour = "black"))+
  xlim(c(-5.5,2.4))+
  ylim(c(-5.5,2.4))+
  geom_text(data = percD,aes(x = -4.8, y = -0.6, label=phits),
            show.legend = FALSE,size=4.5,inherit.aes = F,family="Times",color = "firebrick")
ggsave(Fig4d,filename = glue("{outdir}/Figure4d_double_left_right_comparison_mean_{Dat}.pdf"),width=12,height=3)
#===== Source data output =====
cls <- c("IPC298","MELJUSO","PK1")
Fig4a2c <- lapply(seq_along(rocList),function(s) rocList[[s]]$data%>%
                    mutate(title=glue("Digenic Paralog- {cls[s]}"))%>%
                    mutate(panel = glue("Panel {letters[s]}")))%>%
  bind_rows()
write.xlsx(list("Panel a-c" = Fig4a2c,
                "Panel d" = Fig4d$data),
           glue("{outdir}/Figure4_sourceData_{Dat}.xlsx"),
           overwrite = T)