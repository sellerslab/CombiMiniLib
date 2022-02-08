setwd("~/ResubPrep/")
source("CodeOrg/Utilities.R")
#"The output directory you would like to put"
outdir <- "tmp/ToDiscuss_01_15_22"
Dat <- Sys.Date()
reCal <- F
#===== Supplementary Table Output =====
Cas9_lib<- readRDS("Raw_Data/sgRNASource_Cas9_enCas12a.rds")$Cas9_sgRNA
enCas12a_lib  <- readRDS("Raw_Data/sgRNASource_Cas9_enCas12a.rds")$enCas12a_sgRNA
minilib <- readRDS("Processed_Data/minilib_Init.rds")
annot <- readRDS("Raw_Data/Minilib_Annotation.rds")%>%select(-essClass)
# no 15 copy genes+chr2+AAVS1+nonTarget 617-18=599
# 10 genes only have 3 sgRNAs, paired with AAVS1 to investigate promoter effect
# 589*6+10*3 = 3564 sgRNAs
# Suggest to remove chr2, AAVS1, nonTarget trunk since we never talked about them in the manuscript
write.xlsx(
  list(
    "Supplementary table 1" = Cas9_lib,
    "Supplementary table 2" = annot%>%select(-Cas12a_oligo),
    "Supplementary table 3" = enCas12a_lib,
    "Supplementary table 4" = annot%>%select(-Cas9_oligo),
    "Supplementary table 5" = readRDS("Processed_Data/Minilib_JointCounts.rds")%>%
      mutate(cell.line = case_when(cell.line=="N/A"~"pDNA",TRUE~cell.line))%>%
      mutate(Sample.Name = paste(cell.line,Collapsed_identifier,sep = "_"))%>%
      select(-Collapsed_identifier,-cell.line)%>%
      distinct()%>%
      spread(Sample.Name,pooled_counts)%>%
      replace(is.na(.), 0),
    "Supplementary table 6" = readRDS("Processed_Data/Minilib_JointLFC.rds")%>%as.data.frame()%>%rownames_to_column("sgRNA_Pair"),
    "Supplementary table 7" = minilib$score%>%as.data.frame()%>%rownames_to_column("gene_Pair"),
    "Supplementary table 8" = minilib$posctrl_gene,
    "Supplementary table 9" = minilib$negctrl_gene,
    "Supplementary table 10" = minilib$posctrl_pair,
    "Supplementary table 11" = minilib$negctrl_pair),
  glue("{outdir}/Supplementary_Table.xlsx"),overwrite=T
)

write.xlsx(
  list(
    "Supp1(Cas9 Gene-sgRNA Map)" = Cas9_lib,
    "Supp2(Cas9 Construct Oligo)" = annot%>%select(-Cas12a_oligo),
    "Supp3(enCas12a Gene-sgRNA Map)" = enCas12a_lib,
    "Supp4(enCas12a Construct Oligo)" = annot%>%select(-Cas9_oligo),
    "Supp5(Screen Raw Count)" = readRDS("Processed_Data/Minilib_JointCounts.rds")%>%
      mutate(cell.line = case_when(cell.line=="N/A"~"pDNA",TRUE~cell.line))%>%
      mutate(Sample.Name = paste(cell.line,Collapsed_identifier,sep = "_"))%>%
      select(-Collapsed_identifier,-cell.line)%>%
      distinct()%>%
      spread(Sample.Name,pooled_counts)%>%
      replace(is.na(.), 0),
    "Supp6(Screen Raw LFC)" = readRDS("Processed_Data/Minilib_JointLFC.rds")%>%as.data.frame()%>%rownames_to_column("sgRNA_Pair"),
    "Supp7(Screen Sensitive Score)" = minilib$score%>%as.data.frame()%>%rownames_to_column("gene_Pair"),
    "Supp8(Gene-level pos-control)" = minilib$posctrl_gene,
    "Supp9(Gene-level neg-control)" = minilib$negctrl_gene,
    "Supp10(Pair-level pos-control)" = minilib$posctrl_pair,
    "Supp11(Pair-level neg-control)" = minilib$negctrl_pair),
  glue("{outdir}/Supplementary_Table_TabAnnot.xlsx"),overwrite=T
)

#===== Data Loading & Gene+Pair Level Separation Analysis  =====
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
minilib <- readRDS("Processed_Data/minilib_Init.rds")
SpSaLib <- readRDS("Processed_Data/SpSaLib_Init.rds")
combLib <- list("minilib"=minilib,"SpSaLib"=SpSaLib)
system3_geneSep <- EvalWrapper(combLib)
rocList <- lapply(c("IPC298","MELJUSO","PK1"),function(l){
  if(l=="IPC298"){
    SepBar(system3_geneSep$metrics%>%filter(grepl(l,cl))%>%mutate(cl = gsub(paste0(l,"\\_"),"",cl)),colorp="Paired")%>%
      ggsave(plot = .,filename = glue("{outdir}/GeneSeparation_Bar_{l}_{Dat}.pdf"),width=14,height=5)
  }
  labelMap <- system3_geneSep$metrics%>%
    filter(grepl(l,cl))%>%
    filter(metric =="AUC-ROC")%>%
    arrange(-values)%>%
    mutate(fakeB =  " (",fakeE = ")")%>%
    unite("Label",c("cl","fakeB","values","fakeE"),sep="",remove=F)%>%
    mutate(Label = gsub(paste0(l,"\\_"),"",Label))
  if(l !="IPC298"){
    cols <-c("#1F78B4","#33A02C","#E31A1C","#FB9A99")
  }else{cols <- NULL}
  rocP<- SepROC(system3_geneSep$confusion%>%
                  filter(grepl(l,sample))%>%
                  inner_join(labelMap%>%select(Label,cl),by=c("sample"="cl"))%>%
                  mutate(sample = Label),
                labelLevels=labelMap$Label,ptitle=paste0("Single gene - ",l),
                manucol = cols)
  return(rocP)
})
ggarrange(plotlist = rocList,nrow=1)%>%ggsave(plot=.,glue("{outdir}/GeneSeparation_ROC_{Dat}.pdf"),width = 16.8,height=6.2)
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
ggarrange(plotlist = rocList,nrow=1)%>%ggsave(plot=.,glue("{outdir}/UnionNeg_PairSeparation_ROC_{Dat}.pdf"),width = 16.8,height=6.2)
#===== Manhattan Plot (shared 296/454 pairs, 46 negative control pairs) =====
thres <- -log10(1e-3)
subset_system <- c("IPC298_VCR1-WCR3","IPC298_WCR2-WCR3","IPC298_enCas12a","IPC298_spCas9-saCas9")
rawPair_Comb <- match_matrix(list(combLib$minilib$avg_Double,combLib$SpSaLib$avg_Double),subset_system)
Synergy_Comb <- match_matrix(list(-combLib$minilib$score,-combLib$SpSaLib$score),subset_system)
ncPairs_Shared <- intersect(combLib$minilib$negctrl_pair, combLib$SpSaLib$negctrl_pair)
scores_SigTest(rawPair_Comb,Synergy_Comb,ncPairs_Shared)%>%
  mutate(sample = gsub("IPC298\\_","",sample))%>%
  mutate(sample = factor(sample,levels = gsub("IPC298\\_","",subset_system)))%>%
  ggplot(data = ., aes(x = genepair, y = GEMINI.score, color = LFC)) + 
  facet_grid(. ~ sample, scales = 'free_x', space = 'free_x', switch = 'x') + 
  scale_colour_gradientn(colours=c("red", "#CCCCCC", "#333333" )) +
  geom_point() + 
  scale_y_sqrt() +
  geom_hline(yintercept = thres, color = 'black', linetype = "dashed") +
  labs(x="",y="-log10(FDR)")+
  theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(), 
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", size = rel(0.5), fill = NA),
        strip.background = element_blank(),
        axis.title.y = element_text(size = rel(1), colour = "black",family = "Helvetica"),
        strip.text.x = element_text(size = rel(1), colour = "black", angle = 0,family = "Helvetica"))
ggsave(filename = glue("{outdir}/ManhattanShared_46Neg_IPC298_{Dat}.pdf"),width = 12,height = 6)




