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
source(glue("{CodeOrg}/00_Utilities.R"))
#"The output directory you would like to put"
outdir <- "~/Desktop/Figure2"
if(!dir.exists(outdir)){dir.create(outdir)}
Dat <- Sys.Date()
sessionInfo()
#===== Gene-Level ctrl separation  (Figure 2a-d) =====
minilib <- readRDS(glue("{Processed_Data}/minilib_Init.rds"))
SpSaLib <- readRDS(glue("{Processed_Data}/SpSaLib_Init.rds"))
combLib <- list("minilib"=minilib,"SpSaLib"=SpSaLib)
system3_geneSep <- EvalWrapper(combLib)
rocList <- lapply(c("IPC298","MELJUSO","PK1"),function(l){
  if(l=="IPC298"){
    SepBar(system3_geneSep$metrics%>%filter(grepl(l,cl))%>%mutate(cl = gsub(paste0(l,"\\_"),"",cl)),colorp="Paired")%>%
      ggsave(plot = .,filename = glue("{outdir}/Figure2ac_GeneSeparation_Bar_{l}_{Dat}.pdf"),width=14,height=5)
  }
  labelMap <- system3_geneSep$metrics%>%
    filter(grepl(l,cl))%>%
    filter(metric =="AUC-ROC")%>%
    arrange(-values)%>%
    mutate(fakeB =  " (",fakeE = ")")%>%
    unite("Label",c("cl","fakeB","values","fakeE"),sep="",remove=F)%>%
    mutate(Label = gsub(paste0(l,"\\_"),"",Label))
  if(l !="IPC298"){
    cols <-c("#1F78B4","#B2DF8A","#FB9A99","#E31A1C")
  }else{cols <- NULL}
  rocP<- SepROC(system3_geneSep$confusion%>%
                  filter(grepl(l,sample))%>%
                  inner_join(labelMap%>%select(Label,cl),by=c("sample"="cl"))%>%
                  mutate(sample = Label),
                labelLevels=labelMap$Label,ptitle=paste0("Single gene - ",l),
                manucol = cols)
  return(rocP)
})
ggarrange(plotlist = rocList,nrow=1)%>%ggsave(plot=.,glue("{outdir}/Figure2d_GeneSeparation_ROC_{Dat}.pdf"),width = 14,height=5.2)
#===== Gene-Level left-right comparison (Figure 2e) =====
minilibInt <- readRDS(glue("{Processed_Data}/minilib_Init.rds"))
saspInt <- readRDS(glue("{Processed_Data}/SpSaLib_Init.rds"))
BalanceD_renamed_gene <- minilibInt
BalanceD_renamed_gene$LFC%<>%.[,grep("IPC298",colnames(minilibInt$LFC),value=T)]
colnames(BalanceD_renamed_gene$LFC) <- gsub("IPC298\\_","",colnames(BalanceD_renamed_gene$LFC))
`[` <- function(...) base::`[`(...,drop=FALSE) 
saspInt$LFC%<>%.[,grep("IPC298",colnames(saspInt$LFC),value=T)]
colnames(saspInt$LFC) <- gsub("IPC298\\_","",colnames(saspInt$LFC))
plotOrder <- c("VCR1-WCR3","WCR3-VCR1",
               "WCR2-WCR3","SPCR1-WCR3",
               "enCas12a","spCas9-saCas9",
               "VCR1-SCR27","VCR1-SCR43",
               "SPCR1-SCR43","SPCR1-SCR27")
single_lf_mean <- left_right_scatter(inputlist = list(BalanceD_renamed_gene,saspInt),type="single",FUN="mean",
                                     rowno=2,cor.ordered=F,orders=plotOrder)

percD <- single_lf_mean$data%>%
  filter(side1< -1 & side2 < -1 & Essentiality=="yes")%>%
  group_by(source)%>%
  summarise(n_hit=n())%>%
  ungroup%>%
  arrange(source)%>%
  mutate(n_pos = case_when(source=="spCas9-saCas9"~length(saspInt$posctrl_gene),
                           TRUE~length(minilibInt$posctrl_gene)))

percD%<>%rbind(.,data.frame(source = setdiff(plotOrder,percD$source),
                            n_hit = 0,
                            n_pos = length(saspInt$posctrl_gene)))%>%
  mutate(phits = paste0(round(n_hit*100/n_pos,1),"%"))%>%
  mutate(source = factor(source,levels = plotOrder))

ggsave(single_lf_mean+
         theme(strip.background = element_blank(),
               strip.text.x = element_text(hjust = 0,colour = "black"))+
         geom_text(data = percD,aes(x = -5.3, y = -0.6, label=phits),
                   show.legend = FALSE,size=4.5,inherit.aes = F,family="Times",color = "firebrick"),
       filename = glue("{outdir}/Figure2e_single_left_right_comparison_mean_{Dat}.pdf"),width=13.9,height=6.5)
single_lf_mean$data%>%
  filter(side1< -1 & side2 < -1 & Essentiality=="yes")%>%
  group_by(source)%>%
  summarise(n_hit=n())%>%
  ungroup%>%
  arrange(source)%>%
  mutate(n_pos = case_when(source=="spCas9-saCas9"~length(saspInt$posctrl_gene),
                           TRUE~length(minilibInt$posctrl_gene)))%>%
  mutate(`% of hits` = round(n_hit*100/n_pos,1))%>%
  ggtexttable(., rows = NULL, theme = ttheme("lBlueWhite")) %>%
  tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 2)
ggsave(glue("{outdir}/EssHit_Table.png",width = 3,height = 3))
#===== Source data output =====
cls <- c("IPC298","MELJUSO","PK1")
Fig2d <- system3_geneSep$metrics%>%
  filter(grepl(l,cl))%>%mutate(cl = gsub(paste0(l,"\\_"),"",cl))%>%
  filter(metric=="NNMD")
Fig2a2c <- lapply(seq_along(rocList),function(s) rocList[[s]]$data%>%
                    mutate(title=glue("Single gene- {cls[s]}"))%>%
                    mutate(panel = glue("Panel {letters[s]}")))%>%
  bind_rows()
write.xlsx(list("Panel a-c" = Fig2a2c,
                "Panel d" = Fig2d,
                "Panel e" = single_lf_mean$data),
           glue("{outdir}/Figure2_sourceData_{Dat}.xlsx"),
           overwrite = T)