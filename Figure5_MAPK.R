setwd("~/ResubPrep/")
source("CodeOrg/00_Utilities.R")
#"The output directory you would like to put"
outdir <- "~/Desktop"
Dat <- Sys.Date()
#===== Manhattan Plot (Figure 5a) =====
# (shared 296/454 pairs, 46 negative control pairs)
minilibInt <- readRDS("Processed_Data/minilib_Init.rds")
saspInt <- readRDS("Processed_Data/SpSaLib_Init.rds")
thres <- -log10(1e-3)
subset_system <- c("IPC298_VCR1-WCR3","IPC298_WCR2-WCR3","IPC298_enCas12a","IPC298_spCas9-saCas9")
rawPair_Comb <- match_matrix(list(minilibInt$avg_Double,saspInt$avg_Double),subset_system)
Synergy_Comb <- match_matrix(list(-minilibInt$score,-saspInt$score),subset_system)
ncPairs_Shared <- intersect(minilibInt$negctrl_pair, saspInt$negctrl_pair)
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





#===== Synergy versus RAW gene pair LFC (Figure 5b) =====
ownInterests <- c("BRAF;RAF1","DUSP4;DUSP6","MAP2K1;MAP2K2",
                  "MAPK1;MAPK3","PPP2CA;PPP2CB",
                  "RPS6KA1;RPS6KB1","RPS6KA3;RPS6KB1","YWHAE;YWHAZ")
pS <- 4
lS <- 3.5
minilib <- readRDS("Processed_Data/minilib_Init.rds")
oneEss <- data.frame(gp=rownames(minilib$score))%>%
  filter(!grepl("nonTarget",gp))%>%
  separate(gp,c("g1",'g2'),remove=F)%>%
  filter(g1 %in% minilib$posctrl_gene|g2 %in% minilib$posctrl_gene)
cbind(minilib$avg_Double[rownames(minilib$score),"IPC298_VCR1-WCR3"],
      minilib$score[,"IPC298_VCR1-WCR3"])%>%
  as.data.frame()%>%
  set_colnames(c("AvgLFC","Synergy"))%>%
  rownames_to_column("GenePair")%>%
  filter(!grepl("nonTarget",GenePair))%>%
  mutate(highlights = case_when(GenePair %in% ownInterests~"OwnInterests",
                                GenePair %in% minilib$posctrl_pair~"Pos_Ctrl",
                                GenePair %in% oneEss$gp~"OneGeneEss",
                                TRUE~"Others")%>%factor(.,levels = c("OwnInterests","Pos_Ctrl","OneGeneEss","Others")))%>%
  ggplot(.,aes(x=AvgLFC,y=Synergy,color=highlights,fill=highlights,label=GenePair))+
  geom_rect(aes(xmin = -Inf, xmax = -1, ymin = -Inf, ymax = -1),
            fill = "#fbf7f7",color="grey",linetype="dashed")+
  geom_point(size=rel(pS), alpha=0.8)+
  scale_color_manual(values=c("firebrick","deepskyblue4","goldenrod2","grey50"))+
  scale_fill_manual(values=c("firebrick","deepskyblue4","goldenrod2","grey50"))+
  ggrepel::geom_text_repel(data=.%>%filter(highlights=="OwnInterests"),
                           size = rel(lS),
                           color="black",family="Helvetica",
                           point.padding = 1,arrow=T)+
  xlim(c(-5,1))+ylim(c(-5,1))+
  theme(plot.title = element_text(color='black', hjust = 0.5),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = rel(0.5), fill = NA),
        panel.grid = element_blank(),
        axis.text = element_text(color='black',family="Helvetica"),
        axis.title = element_text(color='black',family="Helvetica"),
        legend.position ="none")+
  geom_hline(yintercept  = 0,linetype = "longdash")+ 
  geom_vline(xintercept  = 0,linetype = "longdash")+
  labs(x="LFC",y="Sensitive Synergy Score",title="VCR1-WCR3: IPC298")
ggsave(filename = glue("{outdir}/Synergy_LFC_IPC298_3cols.pdf"),width = 6,height = 6)v