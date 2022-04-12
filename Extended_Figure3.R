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
outdir <- "~/Desktop/ExtendedFigure3"
if(!dir.exists(outdir)){dir.create(outdir)}
Dat <- Sys.Date()
sessionInfo()
#===== PAM investigation (Extended Figure 3a) =====
minilib <- readRDS(glue("{Processed_Data}/minilib_Init.rds"))
cas12a_pam <- readRDS(glue("{Raw_Data}/enCas12a_fullSource.rds"))%>%
  select(PAM,gRNA,gene=symbol,Tier,ontar_Rank)
cas12a_oligos <- readRDS(glue("{Raw_Data}/Minilib_Annotation.rds"))
mapping_nc_cas12a <- cas12a_oligos%>%
  select(rowname = Cas12a_oligo,Cas9_oligo,Label,symbol1,symbol2)%>%
  separate(rowname,c("gRNA_Left","gRNA_Right"),sep="_",remove=F)%>%
  filter(grepl("AAVS1",Label))%>%
  mutate(Side = case_when(symbol2=="AAVS1"~"l",TRUE~"r"),
         Gene = case_when(symbol2=="AAVS1"~symbol1,TRUE~symbol2))%>%
  select(rowname,Side,Gene,Cas9_oligo)%>%distinct()%>%filter(!Gene %in% c("AAVS1","chr2"))%>%
  mutate(gRNA = case_when(Side=="l"~gsub("\\_.*","",rowname),
                          Side=="r"~gsub(".*\\_","",rowname)))
cas12a_LFC <- minilib$LFC[,grep("enCas12a",colnames(minilib$LFC),value=T)]%>%
  as.data.frame()%>%
  rownames_to_column("Cas9_oligo")%>%
  inner_join(mapping_nc_cas12a%>%select(rowname,Cas9_oligo),
             by="Cas9_oligo")%>%
  select(-Cas9_oligo)%>%
  column_to_rownames("rowname")%>%as.matrix()
mapping_nc_cas12a%<>%cbind(cas12a_LFC[mapping_nc_cas12a$rowname,])
cas12a_single_gene <- mapping_nc_cas12a%>%
  select(Gene,gRNA,IPC298_enCas12a,MELJUSO_enCas12a,PK1_enCas12a)%>%
  mutate(gene = case_when(Gene=="VARS"~"VARS1",
                          Gene=="TARS"~"TARS1",
                          grepl("copy",Gene)~gsub("\\_copy","",Gene),
                          TRUE~Gene))%>%
  inner_join(cas12a_pam%>%distinct(),by= c("gRNA","gene"))%>%
  group_by_at(vars(-IPC298_enCas12a,-MELJUSO_enCas12a,-PK1_enCas12a))%>%
  summarize_at(vars(IPC298_enCas12a,MELJUSO_enCas12a,PK1_enCas12a),mean)%>%
  mutate(Tier = factor(case_when(grepl("TTT[ACG]",PAM)~"TTTV",TRUE~Tier),
                       levels = c("TTTV","Tier 1","Tier 2","Tier 3")))%>%
  as.data.frame()
s3a <- ggplot(cas12a_single_gene%>%
         filter(Gene %in% minilib$posctrl_gene)%>%
         select(Tier,LFC=IPC298_enCas12a),
       aes(x=Tier, y=LFC,fill=Tier)) +
  geom_boxplot()+
  scale_fill_manual(values=c("#a4bce8","#4878d0", "#ee8549","#6acc64"))+
  theme(plot.title = element_text(color='black', hjust = 0.5),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.grid = element_blank(),
        axis.text = element_text(color='black',family = "Times"),
        axis.title = element_text(color='black',family = "Times"),
        legend.key = element_rect(fill = NA),
        legend.position ="none"
  )+
  labs(x="PAM",y="LFC")+
  geom_hline(yintercept =  0,linetype = "dashed")
ggsave(s3a,width = 4,height = 4,filename = glue("{outdir}/ExFigure3a_posgene53_LFC_pam_tier1_{Dat}.pdf"))
#===== On-target Investigation (Extended Figure 3b) =====
s3b <- cas12a_single_gene%>%
  mutate(BIN_ontar = case_when(ontar_Rank<=15~"Top 15",TRUE~"Others"))%>%
  filter(Gene %in%minilib$posctrl_gene)%>%
  select(Gene,BIN_ontar,LFC=IPC298_enCas12a)%>%
  ggplot(.,aes(x=BIN_ontar, y=LFC,fill=BIN_ontar))+
  geom_boxplot()+
  scale_fill_manual(values=c("#4878d0", "#ee8549"))+
  theme(plot.title = element_text(color='black', hjust = 0.5),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.grid = element_blank(),
        axis.text = element_text(color='black',family = "Times"),
        axis.title= element_text(color='black',family = "Times"),
        legend.key = element_rect(fill = NA),
        legend.position ="none",
        axis.title.x = element_blank())
ggsave(s3b,width = 4,height = 4,
       filename = glue("{outdir}/ExFigure3b_ontarget_boxplot_{Dat}.pdf"))
#===== Avana vs RuleSet2Only (Extended Figure 3c) =====
minilib <- readRDS(glue("{Processed_Data}/minilib_Init.rds"))
sgRNA_source <- readRDS(glue("{Raw_Data}/cas9_completeAnnot.rds"))
mapping_12 <- readRDS(glue("{Processed_Data}/SideMapping_Minilib.rds"))
nctrl.sgRNA <- c("TCGATCCGCCCCGTCGTTCC","AGGGAGACATCCGTCGGAGA")
ctrl_pastten <- paste0(paste0("\\_",nctrl.sgRNA,"|",nctrl.sgRNA,"\\_"),collapse = "|")
mapping_12%<>%mutate(sgRNA = gsub(ctrl_pastten,"",rowname))
sgRNASouce <- mapping_12%>%
  mutate(symbol=gsub("\\_copy","",symbol))%>%
  inner_join(sgRNA_source,by=c("symbol"="gene","sgRNA"))%>%
  mutate(source = gsub("; ",";",source))%>%
  distinct()%>%
  mutate(RuleSource = case_when(source %in% c("Avana","Avana;pfam","Avana;GPP;pfam","Avana;GPP")~"Avana",
                                source %in% c("GPP;pfam","GPP")~"Ruleset2"))%>%
  select(RuleSource,rowname,symbol)
LFC <- minilib$LFC%>%
  as.data.frame()%>%
  rownames_to_column("rowname")%>%
  select(rowname,`IPC298_VCR1-WCR3`)%>%
  inner_join(sgRNASouce,by="rowname")
RuleComp <- lapply(list(minilib$posctrl_gene,minilib$negctrl_gene),function(s){
  p <- ggplot(LFC%>%filter(symbol %in% s),
              aes(x=RuleSource, y=`IPC298_VCR1-WCR3`,fill=RuleSource)) +
    geom_boxplot()+
    scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
    theme(plot.title = element_text(color='black', hjust = 0.5),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", size = 1, fill = NA),
          panel.grid = element_blank(),
          axis.text = element_text(color='black'),
          legend.key = element_rect(fill = NA),
          legend.position ="none")+
    geom_hline(yintercept = 0,linetype="dashed")+
    stat_compare_means(method = "wilcox.test",method.args = list(alternative="less"))+
    labs(x=NULL,y="LFC (sgRNA-AAVS1)")
  return(p)})
egg::ggarrange(plots = RuleComp,nrow = 1,top ="Essential (Left) and Non-essential (Right)")%>%
  ggsave(width  =10,height=5,filename = glue("{outdir}/ExFigure3cd_AvanaRs2_comparison_boxplot_{Dat}.pdf"))
#===== Source data output =====
write.xlsx(list("Panel a" = s3a$data,
                "Panel b" = s3b$data,
                "Panel c" = RuleComp[[1]]$data,
                "Panel d" = RuleComp[[2]]$data),
           glue("{outdir}/ExFigure3_sourceData_{Dat}.xlsx"),
           overwrite = T)