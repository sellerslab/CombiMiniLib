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
outdir <- "~/Desktop/ExtendedFigure4"
if(!dir.exists(outdir)){dir.create(outdir)}
Dat <- Sys.Date()
sessionInfo()
#===== Positive Control Examples (Extended Figure 4a) =====
minilib <- readRDS(glue("{Processed_Data}/minilib_Init.rds"))
sasp <- readRDS(glue("{Processed_Data}/SpSaLib_Init.rds"))
synLeth <- fread(glue("{Raw_Data}/CRISPR-paralog.csv",data.table=F))
common.essentials <- load.from.taiga(data.name='avana-public-21q4-e1b8', data.version=4, data.file='common_essentials')
top10_synleth <- synLeth%>%
  select(Dep_Gene,Loss_Gene,genetic_ks_qvalue,genetic_hit)%>%
  filter(genetic_hit==T)%>%
  filter(!(Dep_Gene %in% common.essentials$gene| Loss_Gene %in% common.essentials$gene))%>%
  arrange(genetic_ks_qvalue)%>%
  mutate(gene1 = gsub(" \\(.*","",Dep_Gene),
         gene2 = gsub(" \\(.*","",Loss_Gene))%>%
  select(gene1,gene2,genetic_ks_qvalue)%>%
  dplyr::slice(1:10)%>%
  mutate(pair = pmap_chr(list(gene1,gene2),~paste0(sort(c(...)),collapse = ":")))%>%pull(pair)
s4a <- boxplot_stack(LFC=minilib$LFC,map=minilib$annot,
              gene_pairs=top10_synleth,
              cell_line=c("IPC298_enCas12a","IPC298_WCR2-WCR3","IPC298_VCR1-WCR3"))
ggsave(s4a,width  =10,height=7.5,filename = glue("{outdir}/ExFigure4a_AvanaRs2_comparison_boxplot_{Dat}.pdf"))
#===== Cell Line Variability (Extended Figure 4c) =====
sampleInfo <- load.from.taiga(data.name='avana-public-21q4-e1b8', 
                              data.version=4,
                              data.file="Achilles_sample_info")
sampleInfo%>%
  mutate(cl = gsub("\\_.*","",CCLE_name))%>%
  filter(cl %in% c("IPC298","MELJUSO","PK1"))%>%
  select(cl,cas9_activity)%>%
  ggtexttable(., rows = NULL,
              theme = ttheme("lBlueWhite"))
#ggsave(glue("{outdir}/DepMap_Cas9Activity.png"),width = 3,height = 3)
nnmdCL <- sampleInfo%>%
  mutate(cl = gsub("\\_.*","",CCLE_name))%>%
  select(lineage=primary_tissue,cl,nnmd = cell_line_NNMD)%>%
  na.omit()
tooSmall <- nnmdCL%>%group_by(lineage)%>%summarise(n=n())%>%filter(n<=6)
nnmdCL%<>%filter(!lineage %in% tooSmall$lineage)
orderCL <- nnmdCL%>%
  group_by(lineage)%>%
  summarise(medianNNMD=mean(nnmd))%>%
  arrange(medianNNMD)%>%
  pull(lineage)
nnmdCL%<>%mutate(lineage = factor(lineage,levels=orderCL))
nnmdCL%<>%mutate(h = case_when(cl %in% c("IPC298","MELJUSO","PK1")~"Y",TRUE~"N"))
nnmdCL%<>%mutate(ch = case_when(cl %in% c("IPC298","MELJUSO","PK1")~cl,TRUE~"others"))%>%
  mutate(ch = factor(ch,levels=c("others","PK1","IPC298","MELJUSO")))
s4c <- ggplot(nnmdCL,aes(x = lineage,y = nnmd,label=cl))+
  geom_violin(trim=T)+
  geom_jitter(size=rel(1.5),width = 0.1,aes(color=ch,shape=h))+
  scale_color_manual(values=c("grey50","blue","orange","yellow"))+
  theme(plot.title = element_text(color='black', hjust = 0.5),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        text = element_text(color='black'),
        panel.grid = element_blank(),
        axis.text = element_text(color='black'),
        strip.background = element_rect(colour="black", fill="gray81",
                                        size=1, linetype="solid"),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90,hjust=1,vjust = 1,size=rel(1.3)),
        legend.position = "none"
  )+labs(x="",y="NNMD")+
  ggrepel::geom_label_repel(data=.%>%filter(h=="Y"),
                            size = rel(5),
                            color="black",
                            family="Helvetica",
                            point.padding = 1,arrow=T)
ggsave(s4c,filename = glue("{outdir}/ExFigure4c_lineage_violin_h3_{Dat}.pdf"),width=15,height = 4.5)
#===== Source data output =====
write.xlsx(list("Panel a" = s4a$data,
                "Panel c" = s4c$data),
           glue("{outdir}/ExFigure4_sourceData_{Dat}.xlsx"),
           overwrite = T)