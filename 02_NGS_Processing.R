source("CodeOrg/Utilities.R")
fasta<- "{Your Folder Containing NGS data To Analyze}"
outdir <- "~/Desktop/tmp"
Dat <- Sys.Date()
dirC(outdir)
#===== Count Mapping =====
rlt_folder <- paste0(gsub(".*\\/(SN0.*)","\\1",fasta),"_count")
dir.create(rlt_folder)
oligo_mapping(path_to_fastq = fasta,output_path=rlt_folder)
#===== Count Pooling =====
mapped_info <- readRDS("Raw_Data/Barcode2Sample_Map.rds")
countFolders <- grep("SN0.*\\_count",
                     list.dirs("{Your Folder Storing All Count Folder}",
                               full.names = T,recursive = T),value=T)
counts_path <- data.frame(
  countdir =  countFolders)%>%
  mutate(countpattern = "pairs_count")%>%
  mutate(tracrFlag = F)
pooled_counts <- pbapply::pbapply(counts_path,FUN =counts_pooling_wrapper,
                                  barcodefile=mapped_info,MARGIN=1,cl=15)
pooled_counts%<>%bind_rows()%>%
  group_by(Collapsed_identifier,cell.line,pairs)%>%
  summarise(pooled_counts  = sum(pooled_counts))%>%
  as.data.frame()
saveRDS(pooled_counts,"Processed_Data/Minilib_JointCounts.rds")
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

