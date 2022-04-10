pacman::p_load(glue,eulerr,ggplot2,openxlsx)
source("CodeOrg/Utilities.R")
#"The output directory you would like to put"
outdir <- "~/Desktop/tmp"
Dat <- Sys.Date()
dirC(outdir)
#===== Mini-Lib sgRNA rank and pfam distribution (Extended Figure 1e-h) =====
Cas9_lib<- readRDS("Raw_Data/sgRNASource.rds")$Cas9_sgRNA
enCas12a_lib  <- readRDS("Raw_Data/sgRNASource.rds")$enCas12a_sgRNA
#S1e
pfam_cas9 <- Cas9_lib%>%filter(grepl("pfam",source))%>%nrow()
pdf(glue("{outdir}/Cas9_sgRNASourcePie_pFAM_{Dat}.pdf"),width = 6,height = 6)
pie(c(pfam_cas9,nrow(Cas9_lib)-pfam_cas9),labels = c("pFAM (62%)","non-pFAM"), 
    density=80 ,col=c("#0073C2FF", "#EFC000FF"))
dev.off()
Cas9_sources <- c("Avana","GPP")
pdf(glue("{outdir}/Cas9_sgRNASourceVenn_AvanaRuleSet2_{Dat}.pdf"))
lapply(Cas9_sources,function(s) Cas9_lib%>%filter(grepl(s,source))%>%pull(gRNA))%>%vennCount(.,c("Avana","Rule Set2"))%>%
  euler()%>%
  plot(., counts = TRUE, font=2, cex=2, alpha=0.6,quantities = TRUE,
       lty =1, fontfamily =4,fill=c("#0073C2FF", "#CD534CFF"))
dev.off()
#S1f
Cas9_lib%>%
  filter(grepl("GPP",source))%>%
  mutate(GPPrank = gsub(".*GPP: (.*)","\\1",eval)%>%as.numeric())%>%
  ggplot(., aes(x=GPPrank))+
  geom_histogram(position="identity", alpha=0.5,fill="#E69F00",color = "#E69F00",binwidth = 1)+
  geom_vline(xintercept = 7,linetype="longdash")+
  annotate("text", x = 75, y = 320,label = "Median Pick Order: 7")+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(color='black',family = "Helvetica"),
        axis.title = element_text(color='black',family = "Helvetica"),
        panel.border = element_rect(color = "black", size = 1, fill = NA))+
  labs(x="Pick Order (Rule Set2 + Off-target)",y="Frequency",title = "Cas9 sgRNAs")
ggsave(filename = glue("{outdir}/Cas9_RuleSet2_PickOrder_{Dat}.pdf"),width = 4,height = 4)
#S1g
pfam_cas12a <- enCas12a_lib%>%filter(pFAM!="")%>%nrow()
pdf(glue("{outdir}/enCas12a_sgRNASourcePie_pFAM_{Dat}.pdf"),width = 6,height = 6)
pie(c(pfam_cas12a,nrow(enCas12a_lib)-pfam_cas12a),labels = c("pFAM (65%)","non-pFAM") , density=80 ,col=c("#0073C2FF", "#EFC000FF"))
dev.off()
#S1h
ggplot(enCas12a_lib, aes(x=Order))+
  geom_histogram(position="identity", alpha=0.5,fill="#56B4E9",color = "#56B4E9",binwidth = 1)+
  geom_vline(xintercept = 5.5,linetype="longdash")+
  annotate("text", x = 30, y = 420,label = "Median Pick Order: 5.5")+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(color='black',family = "Courier"),
        axis.title = element_text(color='black',family = "Courier"),
        panel.border = element_rect(color = "black", size = 1, fill = NA))+
  labs(x="Pick Order (enPAM + GB + Off-target)",y="Frequency",title = "enCas12a sgRNAs")
ggsave(filename = glue("{outdir}/enCas12a_PickOrder_{Dat}.pdf"),width = 4,height = 4)

#===== Library Distribution (Extended Figure 1i) =====
pooled_counts <- readRDS("Processed_Data/Minilib_JointCounts.rds")
minis <- pooled_counts%>%
  filter(cell.line=="N/A")%>%
  pull(Collapsed_identifier)%>%unique()
count_list <- lapply(minis,function(s)
  pooled_counts%>%
    filter(cell.line=="N/A" & Collapsed_identifier==s)%>%
    pull(pooled_counts))
names(count_list) <- minis
minilib_libs<- ETP_diagnostic_plot(count_list)
mycolors = c(brewer.pal(name="Dark2", n = 8)%>%rev(),"#2baae0","#9e1010")
colorLabels <- data.frame(libs =c("Optimal screen",names(count_list)),auc = c(0.5,round(minilib_libs$auc,3)))%>%
  arrange(auc)%>%
  mutate(label  = paste0(libs," (",format(round(auc+1e-3,2), nsmall = 2),")"))%>%
  pull(label)
minilib_libs$plot+ 
  scale_color_manual(values=mycolors,labels=colorLabels,name="Library (AUC)")+
  theme(plot.title = element_text(color='black'),
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        plot.background = element_blank(),
        panel.background = element_blank(),
        text = element_text(color='black', size = 14),
        panel.grid = element_blank(),
        axis.text = element_text(color='black',family="Helvetica"),
        axis.title = element_text(color='black',family="Helvetica"),
        legend.key = element_rect(fill = "white"),
        legend.position = c(0.95, 0.05),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.title=element_text(family="Helvetica"),
        legend.text=element_text(family="Helvetica")) + 
  guides(color = guide_legend(order = 1,
                              family="Times",
                              override.aes = list(size=4,fill=NA)))+
  labs(x = "Relative rank", y = "Cumulative Fraction",title="Library distribution")
ggsave(glue("{outdir}/Library_distribution_minilib_{Dat}.pdf"),height = 6.8,width=6.5)


