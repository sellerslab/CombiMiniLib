pacman::p_load(glue,eulerr,ggplot2,openxlsx)
source("CodeOrg/Utilities.R")
#"The output directory you would like to put"
outdir <- "~/Desktop/tmp"
Dat <- Sys.Date()
dirC(outdir)
#===== Essential/Non-essential Gene/Gene-pair (Figure 1c) =====
minilibComp <- readRDS("Raw_Data/Minilib_Annotation.rds")
minilibComp%>%
  filter(essClass=="pan_essentialPairs")%>%
  mutate(pair = pmap_chr(list(symbol1,symbol2),~paste0(sort(c(...)),collapse = "_")))%>%
  pull(pair)%>%
  grep("AAVS1",.,invert=T,value=T)%>%
  unique()%>%gsub("\\_",";",.)%>%
  setdiff(.,symtem3_compare$minilib$posctrl_pair)

perc <- round(table(minilibComp$essClass)*100/nrow(minilibComp),0)
perc[1:2] <- perc[c(2,1)]
names(perc)[1:2] <-  names(perc)[c(2,1)]
labels <- c("Non-essentiall pair",
            "Non-essentiall gene",
            "Other",
            "Essential gene",
            "Essential pair")
pdf(glue("{outdir}/PieEssentialClass_Distribution.pdf"),width = 8,height = 8)
pie(perc,labels = paste0(labels," ",perc,"%") ,
    col=c("#547c85", "#99d1ac","#e2db79","#575e71","#b1d686"),
    init.angle=190)
dev.off()
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
