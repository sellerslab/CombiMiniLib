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
outdir <- "~/Desktop/Figure3"
if(!dir.exists(outdir)){dir.create(outdir)}
Dat <- Sys.Date()
sessionInfo()
#===== Balance Dissection (Figure 3a-b)  =====
### Fig 3a below: Dissection of VCR1 compared with WCR3 by removing influences of sgRNA choices and promoter
minilibInt <- readRDS(glue("{Processed_Data}/minilib_Init.rds"))
mapping_12 <-minilibInt$annot%>%
  filter(symbol1=="AAVS1"|symbol2=="AAVS1")%>%
  mutate(Side = case_when(symbol1=="AAVS1"~"side2",TRUE~"side1"),
         symbol = case_when(symbol1=="AAVS1"~symbol2,TRUE~symbol1))%>%
  select(rowname,Side,symbol)%>%distinct()
fig3a <- minilibInt$LFC[,c("IPC298_VCR1-WCR3","IPC298_WCR3-VCR1")]%>%
  as.data.frame()%>%
  rownames_to_column("rowname")%>%
  inner_join(mapping_12,by="rowname")%>%
  mutate(`tracrDiff_VCR1-WCR3` = case_when(Side=="side1"~`IPC298_VCR1-WCR3`-`IPC298_WCR3-VCR1`,
                                           TRUE~`IPC298_WCR3-VCR1`-`IPC298_VCR1-WCR3`))%>%
  ggplot(.,aes(x=`tracrDiff_VCR1-WCR3`,color=Side,fill=Side))+
  geom_density(alpha=0.5)+
  scale_color_manual(values=c("#009E73","#0072B2"),labels=c("U6","H1"),name="")+
  scale_fill_manual(values=c("#009E73","#0072B2"),labels=c("U6","H1"),name="")+
  theme(plot.title = element_text(color='black', hjust = 0.5),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.grid = element_blank(),
        axis.text = element_text(color='black'),
        legend.key = element_rect(fill = "white"),
        legend.position = c(0.99, 0.01),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(5, 5, 5, 5),
        legend.text=element_text(size=rel(1),family = "Helvetica",face="bold")
  )+
  geom_vline(xintercept = 0,linetype="dashed")+ 
  labs(x="Strength Difference in tracrRNA (VCR1 vs WCR3)")
ggsave(fig3a,filename = glue("{outdir}/Figure3a_tracrRNA_strengthDiff_{Dat}.pdf"),width = 6,height = 6)
### Fig 3b below: Dissection of promoter effect by removing influences of sgRNA choices and tracrRNA
nctrl.sgRNA <- c("TCGATCCGCCCCGTCGTTCC","AGGGAGACATCCGTCGGAGA")
ctrl_pastten <- paste0(paste0("\\_",nctrl.sgRNA,"|",nctrl.sgRNA,"\\_"),collapse = "|")
mapping_12%<>%mutate(sgRNA = gsub(ctrl_pastten,"",rowname))
promSub <- mapping_12%>%select(sgRNA,Side)%>%distinct()%>%
  group_by(sgRNA)%>%summarise(n=n())%>%filter(n==2)%>%pull(sgRNA)
fig3b <- minilibInt$LFC[,c("IPC298_VCR1-WCR3","IPC298_WCR3-VCR1")]%>%
  as.data.frame()%>%
  rownames_to_column("rowname")%>%
  gather(source,lfc,-rowname)%>%
  inner_join(mapping_12%>%filter(sgRNA %in% promSub),by="rowname")%>%
  mutate(source = gsub("IPC298\\_","",source))%>%
  separate(source,into=c("side1","side2"),sep="-",remove=T)%>%
  mutate(Promoter = case_when(Side=="side1"~"U6",TRUE~"H1"),
         tracrRNA = case_when(Side=="side1"~side1,TRUE~side2))%>%
  select(sgRNA,lfc,Promoter,tracrRNA)%>%
  spread(Promoter,lfc)%>%
  mutate(`Effect Size` = abs(U6-H1))%>%
  ggplot(.,aes(x=U6,y=H1,fill=tracrRNA,shape=tracrRNA,size=`Effect Size`))+
  geom_point(color="grey50")+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=c('#F16B4E','#56B4E9'))+
  theme(plot.title = element_text(color='black', hjust = 0.5),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.grid = element_blank(),
        axis.text = element_text(color='black'),
        legend.key = element_rect(fill = "white"),
        legend.position = c(0.99, 0.01),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(5, 5, 5, 5),
        legend.text=element_text(size=rel(1),family = "Helvetica",face="bold")
  )+
  labs(x="LFC from U6 promoter",y="LFC from H1 promoter")+
  geom_abline(intercept = 0,color="grey48",linetype="longdash")+ 
  xlim(c(-6.5,1))+ylim(c(-6.5,1))+
  guides(color = guide_legend(order = 1,override.aes = list(size=8,fill=NA)),
         shape = guide_legend(override.aes = list(size=8)))
ggsave(fig3b,filename = glue("{outdir}/Figure3b_PromoterEffect_scatter_{Dat}.pdf"),width = 6,height = 6)
#===== Recombination Rate (Figure 3d) =====
recombMap <- readRDS(glue("{Processed_Data}/pDNA_mappingD.rds"))
mini_mapped <- recombMap$mini_mapped
mini_mapped_tracr <- recombMap$mini_mapped_tracr
recomb_df <- mini_mapped%>%
  as.data.frame()%>%
  filter(final_note=="both")%>%
  select(Sample.Name,Barcode,sgRNA_only = perc)%>%
  inner_join(mini_mapped_tracr%>%
               as.data.frame()%>%
               filter(final_note=="both")%>%
               select(Sample.Name,Barcode,sgRNA_tracrRNA = perc),by=c("Sample.Name","Barcode"))%>%
  mutate(recomb_rate = 100-(sgRNA_tracrRNA*100/sgRNA_only))%>%
  mutate(tracrComb = gsub(".*\\_(.*\\_.*)\\_PDNA","\\1",Sample.Name))%>%
  group_by(tracrComb)%>%
  summarize(meanRecomb = mean(recomb_rate))%>%
  arrange(meanRecomb)
recomb_df%<>%
  mutate(tracrComb = gsub("_","-",tracrComb))%>%
  mutate(tracrComb = forcats::fct_reorder(tracrComb, -meanRecomb))
fig3d <- ggplot(recomb_df%>%filter(tracrComb %in% c("VCR1-WCR3","WCR2-WCR3")), 
       aes(x=tracrComb, y=100-meanRecomb, fill=tracrComb)) +
  geom_bar(stat="identity")+
  scale_fill_manual(values = c("#64b56e","#5193c4"))+
  theme(plot.title = element_text(color='black', hjust = 0.5),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 0.25, fill = NA),
        panel.grid = element_blank(),
        axis.text.x = element_text(color='black',family = "Helvetica",angle = 90),
        axis.title.y =  element_text(color='black',family = "Helvetica"),
        axis.title.x = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.position ="none")+
  labs(y="Average Mapped sgRNA/Mapped crRNA")+
  ylim(c(0,100))
ggsave(fig3d,filename = glue("{outdir}/Figure3d_pDNA_Mappability_{Dat}.pdf"),width = 5,height = 8.5)
#===== Reverse Manhattan Mappability for CR3  (Figure 3e) =====
recombLFC <- readRDS(glue("{Processed_Data}/Recomb_LFC.rds"))
mapping_12 <- readRDS(glue("{Processed_Data}/SideMapping_Minilib.rds"))
minilibInt <- readRDS(glue("{Processed_Data}/minilib_Init.rds"))
desireDesign <- mapping_12%>%group_by(symbol)%>%summarise(n=n())%>%filter(n==6)%>%pull(symbol)
mapping_12%<>%filter(symbol %in% desireDesign)
WCR3_H1 <- mapping_12%>%filter(Side=="side2")
WCR3_AlignD <- recombLFC%>%
  as.data.frame()%>%
  rownames_to_column("rowname")%>%
  inner_join(WCR3_H1,by="rowname")%>%
  select(-Side)%>%
  gather(comb,LFC,-rowname,-symbol)%>%
  mutate(WCR3_type = case_when(grepl("VCR1_WCR3",comb)~"WCR3_VCR1", TRUE~"WCR3_WCR2"),
         Recomb = case_when(grepl("Recomb",comb)~"sgRNA",TRUE~"crRNA"))%>%
  select(-comb)%>%
  spread(WCR3_type,LFC)
min_lim <- floor(apply(WCR3_AlignD%>%select(WCR3_WCR2,WCR3_VCR1),2,min,na.rm=T)%>%min())
max_lim <- ceiling(apply(WCR3_AlignD%>%select(WCR3_WCR2,WCR3_VCR1),2,max,na.rm=T)%>%max())
WCR3_AlignD_gene <- WCR3_AlignD%>%
  group_by(symbol,Recomb)%>%
  summarise_at(vars(WCR3_WCR2,WCR3_VCR1),mean)%>%
  as.data.frame()
percD <- WCR3_AlignD_gene%>%
  gather(tracr,avgLFC,-Recomb,-symbol)%>%
  filter(symbol %in% minilibInt$posctrl_gene)%>%
  group_by(tracr,Recomb)%>%
  dplyr::summarise(p_lm1 = sum(avgLFC < -1),
                   p_lm2 = sum(avgLFC < -2),
                   p_lm3 = sum(avgLFC < -3))%>%
  mutate_at(vars(p_lm1,p_lm2,p_lm3),function(s)round(s*100/length(unique(minilibInt$posctrl_gene)),1))%>%
  mutate(phits = paste0(p_lm3,"%"))%>%
  mutate(tracr = factor(tracr,levels=c("WCR3_WCR2","WCR3_VCR1")))
fig3e <- ggplot(data = WCR3_AlignD_gene%>%
         gather(tracr,avgLFC,-Recomb,-symbol)%>%
         mutate(tracr = factor(tracr,levels=c("WCR3_WCR2","WCR3_VCR1"))), 
       aes(x = symbol, y = avgLFC, color = tracr)) + 
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -3),
            fill = "#fbf7f7",linetype = "longdash",color = "tomato2")+
  facet_wrap(. ~ Recomb+tracr, scales = 'free_x', nrow=1) + 
  scale_color_manual(values = c("#64b56e","#5193c4"))+
  geom_point()+
  geom_hline(yintercept = -3, color = 'tomato4', linetype = "dashed") +
  geom_hline(yintercept = 0, color = 'black', linetype = "dashed") +
  xlab("") +
  ylab("Gene Level LFC") +
  theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(), 
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        strip.background = element_blank(),
        axis.title.y = element_text(size = 15, colour = "black",family="Times"),
        strip.text.x = element_text(size = 15, colour = "black", angle = 0,family="Times"),
        legend.position = "none")+
  geom_text(data = percD,aes(x = 300, y = -5.8, label=phits),
            show.legend = FALSE,size=4.5,inherit.aes = F,family="Times",color = "firebrick")
ggsave(fig3e,filename = glue("{outdir}/Figure3e_Reverse_Manhattan_Plot_{Dat}.pdf"),width = 6,height=8)
#===== Source data output =====
write.xlsx(list("Panel a" = fig3a$data,
                "Panel b" = fig3b$data,
                "Panel d" = fig3d$data,
                "Panel e" = fig3e$data),
           glue("{outdir}/Figure3_sourceData_{Dat}.xlsx"),
           overwrite = T)