#===== Libraries =====
pacman::p_load(stringr,Biostrings,glue,mgsub,
               magrittr,dplyr,tidyr,purrr,
               readr,openxlsx,
               ggpubr,ggtext,eulerr,ggplot2,RColorBrewer,
               pbmcapply,pbapply,taigr)
import::from(data.table,"fread")
import::from(tibble,"rownames_to_column","column_to_rownames")
options(stringsAsFactors = F)
#===== General Unitilities ======
dirC <- function(fv){invisible(sapply(fv,function(s) if(!dir.exists(s)){dir.create(s)}))}
revComp <- function(x){paste(rev(strsplit(chartr("ATGC","TACG",x),"")[[1]]),collapse="")}
fixed_levels <- function(s,f){
  us <- unique(s)
  l <- c(us[us!=f],f)
  fs <- factor(s,levels=l)
  return(fs)
}
simple_auc <- function(y, x){
  # inputs already sorted, best scores first
  dy <- c(diff(y), 0)
  dx <- c(diff(x), 0)
  sum(y * dx) + sum(dy * dx)/2
}
#===== Specific Plots ======
vennCount <- function(l,n){
  vc <- c(lengths(l),length(Reduce(intersect,l)))
  vc[1:2] <- vc[1:2]-vc[3]
  names(vc) <- c(n,paste0(n,collapse ="&"))
  return(vc)
}
ETP_diagnostic_plot <- function(ETP_counts,lib=NULL){
  library(ggplot2)
  library(tidyverse)
  df.control = data.frame(cumulative_fraction = seq(0,1,0.0001),
                          relative_rank = seq(0,1,0.0001),
                          Screens = "Optimal screen")
  if(is.list(ETP_counts)){
    if(is.null(names(ETP_counts)))
      names(ETP_counts) <- 1:length(ETP_counts)
    
    df.lib <- lapply(1:length(ETP_counts), function(i){
      c <- ETP_counts[[i]]
      name <- names(ETP_counts)[i]
      etp <- sort(c, decreasing = T)
      cumfrac <- cumsum(etp) / sum(etp)
      ranketp <- (1:length(etp)) / length(etp)
      dat = data.frame(cumulative_fraction = cumfrac,
                       relative_rank = ranketp,
                       Screens = name)
      return(dat)
    }) 
    auc <- sapply(df.lib,function(s) simple_auc(sort(s$cumulative_fraction),sort(s$relative_rank)))
    df.lib%<>%bind_rows()
    df = rbind(df.control, df.lib)%>%
      mutate(Screens = factor(Screens,levels=c("Optimal screen",names(ETP_counts)[order(auc)])))
  }else if(is.vector(ETP_counts)){
    etp <- sort(ETP_counts, decreasing = T)
    cumfrac <- cumsum(etp) / sum(etp)
    ranketp <- (1:length(etp)) / length(etp)
    df.lib = data.frame(cumulative_fraction = cumfrac,
                        relative_rank = ranketp,
                        Screens = lib)
    auc <- simple_auc(sort(df.lib$cumulative_fraction),sort(df.lib$relative_rank))
    df = rbind(df.control, df.lib)
  }
  
  else{
    stop("Unable to process ETP_counts in current format.")
  }
  g = ggplot(data = df, aes(x = relative_rank, y = cumulative_fraction, color = Screens)) +
    geom_smooth(size=1)
  return(list(plot=g,auc=auc))
}
#===== fastq mapping ======
map_comparison <- function(s,ind1,ind2,GuideL_map,GuideR_map){
  a = ind1[[s]]
  b = ind2[[s]]
  output <- data.frame(l_map = length(a),
                       r_map = length(b))%>%
    mutate(pairs = case_when(l_map==1&r_map==1~paste0(GuideL_map[a[1],1],"_",GuideR_map[b[1],1]),
                             TRUE~"n/a"))
  return(output)
}
Cas12a_map <- function(x,DR1,backbone,lmap,rmap,tmpdir,output_path,oligos){
  FILE = x[[1]]
  sample_bc = paste0(x[[2]],"_l",x[[3]])
  fa <- paste0(tmpdir,"/",sample_bc,"_Cas12a_fastq.txt")
  system(paste0("gunzip -c ", FILE," | awk '{if(NR%4==2) print;}' > ",fa), intern = F)
  con = file(fa, open = 'r')
  fullfastq <- readLines(con)
  close(con)
  r <- DNAStringSet(fullfastq)
  n <- length(fullfastq)
  ind_backbone <- vwhichPDict(PDict(backbone, tb.start = 1, tb.end = -1), r) 
  AGATguide1_index <- vwhichPDict(
    PDict(lmap$guide1_substr, tb.start = 1, tb.end = -1),
    DNAStringSet(fullfastq)) 
  guide26T_index <- vwhichPDict(
    PDict(rmap$guide2_substr, tb.start = 1, tb.end = -1),
    DNAStringSet(fullfastq)) 
  ind_DR <- vwhichPDict(PDict(DR1, tb.start = 1, tb.end = -1), r) 
  map_rlts <- pblapply(1:n,FUN=map_comparison,ind1=AGATguide1_index,ind2 = guide26T_index,
                       GuideL_map = lmap, GuideR_map = rmap)%>%bind_rows()
  map_rlts%<>%mutate(pair_correct = case_when(pairs %in% oligos$Label~1,TRUE~0),DR = lengths(ind_DR))
  n_back = length(lengths(ind_backbone)[lengths(ind_backbone)==1])
  perc_calc <- map_rlts%>%
    mutate(final_note = case_when(l_map==1&r_map==1&pair_correct==1&DR==1~"both",
                                  l_map>=1&r_map==0&DR==1~"left_only",
                                  l_map==0&r_map>=1&DR==1~"right_only",
                                  l_map==1&r_map==1&pair_correct==0&DR==1~"mismatch",
                                  l_map==0&r_map==0&DR==1~"None",
                                  TRUE~"Others"))%>%
    group_by(final_note)%>%
    summarise(sum_n = n())%>%
    mutate(n_lines = n,barcode = sample_bc,n_backbone = n_back)
  rm(fullfastq)
  pair_counts <- map_rlts%>%filter(pair_correct==1&DR==1)%>%group_by(pairs)%>%summarise(count = n())%>%
    mutate(n_lines = n,barcode = sample_bc)
  saveRDS(pair_counts,paste0(output_path,"/",sample_bc,"_pairs_count.rds"))
  saveRDS(perc_calc,paste0(output_path,"/",sample_bc,"_mappability_calc.rds"))
}
Cas9_map <- function(x,lmap,rmap,tmpdir,output_path,oligos,tracrfasta){
  tracrflag <- tracrfasta!="" 
  if(tracrflag){
    tracrSeq <- read.delim(tracrfasta,header=F)
    tracr_mp <- data.frame(tracr = gsub(">","",tracrSeq$V1[seq(1,nrow(tracrSeq),2)]),
                           seq = tracrSeq$V1[seq(2,nrow(tracrSeq),2)])%>%
      mutate(len = nchar(seq),
             ind = case_when(grepl("Revcomp",tracr)~2,TRUE~1))
    tracrs <- unlist(x[6:7])%>%as.character()
    tracr_lens <- tracr_mp%>%
      filter(tracr %in% tracrs)%>%
      arrange(ind)%>%pull(len)
  }
  sample_bc = paste0(x[[1]],"_l",x[[2]])
  FILE_e1 = x[[3]]
  FILE_e2 = x[[4]]
  # Convert fastq to txt with each sequence on one line
  fa1 <- paste0(tmpdir,"/",sample_bc,"_fastq1.txt")
  fa2 <- paste0(tmpdir,"/",sample_bc,"_fastq2.txt")
  system(paste0("gunzip -c ", FILE_e1," | awk '{if(NR%4==2) print;}' > ",fa1), intern = F)
  system(paste0("gunzip -c ", FILE_e2," | awk '{if(NR%4==2) print;}' > ",fa2), intern = F)
  
  if(tracrflag){
    fasta1 <- paste0(tmpdir,"/",sample_bc,"_fastq1.fasta")
    fasta2 <- paste0(tmpdir,"/",sample_bc,"_fastq2.fasta")
    header_file <- paste0(tmpdir,"/",sample_bc,"_header.txt")
    tracrMap_1 <- paste0(tmpdir,"/",sample_bc,"_end1.csv")
    tracrMap_2 <- paste0(tmpdir,"/",sample_bc,"_end2.csv")
    system(paste0("gunzip -c ",FILE_e1," | gsed -n '1~4s/^@/>/p;2~4p' > ",fasta1),intern = F)
    system(paste0("gunzip -c ",FILE_e2," | gsed -n '1~4s/^@/>/p;2~4p' > ",fasta2),intern = F)
    system(paste0("gunzip -c ", FILE_e1," | awk '{if(NR%4==1) print;}' > ",header_file), intern = F)
    header_con = file(header_file,open='r')
    headers <- readLines(header_con)
    close(header_con)
  }
  
  # Read in full fasta files into memory - intensive.
  con1 = file(fa1, open = 'r')
  con2 = file(fa2, open = 'r')
  
  fullfasta1 <- readLines(con1)
  fullfasta2 <- readLines(con2)
  
  n1 <- length(fullfasta1)
  n2 <- length(fullfasta2)
  if(n1!=n2){message("File size not consistent!")}
  close(con1)
  close(con2)
  
  guideL_index_extended <- vwhichPDict(
    PDict(lmap$guide1_substr, tb.start = 1, tb.end = -1),
    DNAStringSet(fullfasta1)) 
  guideR_index_extended <- vwhichPDict(
    PDict(rmap$guide2_substr, tb.start = 1, tb.end = -1),
    DNAStringSet(fullfasta2)) 
  if(tracrflag){
    system(paste("blastn -query ",fasta1,
                 "-subject ",tracrfasta,
                 "-task megablast","-outfmt 6 -dust no -evalue 0.5 -gapopen 5 -gapextend 4 -strand plus",
                 "-out ",tracrMap_1),intern=F)
    system(paste("blastn -query ",fasta2,
                 "-subject ",tracrfasta,
                 "-task megablast","-outfmt 6 -dust no -evalue 0.5 -gapopen 5 -gapextend 4 -strand plus",
                 "-out ",tracrMap_2),intern=F)
    
    tracrRNA_map1 <- read_delim(tracrMap_1,col_names = F,delim = "\t")%>%
      set_colnames(c("query_id","subject_id","identity_perc",
                     "overlap_length","mismatch","no_gap_open",
                     "query_start", "query_end","subject_start","subject_end",
                     "evalue","bitscore"))%>%
      filter(subject_id==tracrs[1]&overlap_length>(tracr_lens[1]-3)&mismatch<2&no_gap_open<=2)
    
    tracrRNA_map2 <- read_delim(tracrMap_2,col_names = F,delim = "\t")%>%
      set_colnames(c("query_id","subject_id","identity_perc",
                     "overlap_length","mismatch","no_gap_open",
                     "query_start", "query_end","subject_start","subject_end",
                     "evalue","bitscore"))%>%
      filter(subject_id==tracrs[2]&overlap_length>(tracr_lens[2]-5)&mismatch<2&no_gap_open<=2)
    
    tracrRNA_f <- ifelse(gsub("@(.*) [12]:N:0:.*","\\1",headers) %in% intersect(tracrRNA_map1$query_id,tracrRNA_map2$query_id),1,0)
  }
  
  message("Guide mapping ... ...")
  map_rlts <- pblapply(1:n1,FUN=map_comparison,
                       ind1=guideL_index_extended,ind2 = guideR_index_extended,
                       GuideL_map = lmap,GuideR_map = rmap)%>%
    bind_rows()
  map_rlts%<>%mutate(pair_correct = case_when(pairs %in% oligos$Label~1,TRUE~0))
  perc_calc <- map_rlts%>%group_by(l_map,r_map,pair_correct)%>%summarise(n=n())%>%
    mutate(final_note = case_when(l_map==1&r_map==1&pair_correct==1~"both",
                                  l_map>=1&r_map==0~"left_only",
                                  l_map==0&r_map>=1~"right_only",
                                  l_map==1&r_map==1&pair_correct==0~"mismatch",
                                  l_map==0&r_map==0~"None",
                                  TRUE~"Others"))%>%
    group_by(final_note)%>%summarise(sum_n = sum(n))%>%
    mutate(n_lines = n1,barcode = sample_bc)
  pair_counts <- map_rlts%>%
    filter(pair_correct==1)%>%
    group_by(pairs)%>%
    summarise(pair_count = n())%>%
    mutate(n_lines = n1)%>%
    mutate(barcode = sample_bc)
  saveRDS(pair_counts,paste0(output_path,"/",sample_bc,"_pairs_count.rds"))
  saveRDS(perc_calc,paste0(output_path,"/",sample_bc,"_mappability_calc.rds"))
  file.remove(c(fa1,fa2))
  if(tracrflag){
    perc_calc_tracr <- map_rlts%>%cbind(tracrRNA_f)%>%
      group_by(l_map,r_map,pair_correct,tracrRNA_f)%>%summarise(n=n())%>%
      mutate(final_note = case_when(l_map==1&r_map==1&pair_correct==1&tracrRNA_f==1~"both",
                                    l_map>=1&r_map==0&tracrRNA_f==1~"left_only",
                                    l_map==0&r_map>=1&tracrRNA_f==1~"right_only",
                                    l_map==1&r_map==1&pair_correct==0&tracrRNA_f==1~"mismatch",
                                    l_map==0&r_map==0&tracrRNA_f==1~"None",
                                    TRUE~"Others"))%>%
      group_by(final_note)%>%summarise(sum_n = sum(n))%>%
      mutate(n_lines = n1,barcode = sample_bc)
    pair_counts_tracr <- map_rlts%>%
      filter(pair_correct==1&tracrRNA_f==1)%>%
      group_by(pairs)%>%
      summarise(pair_count = n())%>%
      mutate(n_lines = n1)%>%
      mutate(barcode = sample_bc)
    saveRDS(pair_counts_tracr,paste0(output_path,"/",sample_bc,"_pairs_count_tracrMAP.rds"))
    saveRDS(perc_calc_tracr,paste0(output_path,"/",sample_bc,"_mappability_calc_tracrMAP.rds"))
    file.remove(c(fa1,fa2,fasta1,fasta2,header_file,tracrMap_1,tracrMap_2))
  }
}
oligo_mapping <- function(path_to_fastq,output_path,tracrfasta,Cas="Cas12a",mcc=15,ended="paired",tmpdir="~/.tmp",
                          oligofile = "Raw_Data/MiniLib_OligoInfo.xlsx",Multiple.tracr=F,testmode=F,
                          tracrRNA_p = "[A-Z][0-9]{2}\\_.*\\_(.*)\\_PDNA"){
  tracrRNA <- tracrfasta!=""
  if(!ended %in% c("single","paired")){message("Please choose mode single or paired!")}
  if(ended=="paired"){
    file_list <- data.frame(file=dir(path = path_to_fastq,pattern = paste0('.*.unmapped.[12]{1}.fastq.gz'), full.names = T),
                            stringsAsFactors = F)%>%
      mutate(Barcode = gsub(".*\\_[A-Z0-9]{9}.[1,2,3,4]{1}.(.*).unmapped.*","\\1",file),
             Lane = gsub(".*/(.*)_[A-Z0-9]{9}.*","\\1",file),
             End = gsub(".*unmapped.(.*).fastq.gz","\\1",file))%>%
      spread(End,file)%>%
      set_colnames(c("Barcode","Lane","End1file","End2file"))
    
  }else{
    file_list <- data.frame(file=dir(path = path_to_fastq,pattern = paste0('.*.unmapped.[12]{1}.fastq.gz'), full.names = T),
                            stringsAsFactors = F)%>%
      mutate(Barcode = gsub(".*\\_[A-Z0-9]{9}.[1,2,3,4]{1}.(.*).unmapped.*","\\1",file),
             Lane = gsub(".*/(.*)\\_[A-Z0-9]{9}.*","\\1",file))
  }
  
  if(Cas=="Cas9"){
    oligos <-  openxlsx::read.xlsx(oligofile,sheet="Cas9")%>%
      select(Label,Sequence)%>%
      mutate(LabelL = gsub("(.*_L[0-9]{1,2})\\_.*\\_R[0-9]{1,2}","\\1",Label),
             LabelR = gsub(".*_L[0-9]{1,2}\\_(.*\\_R[0-9]{1,2})","\\1",Label),
             guide1_substr = str_sub(Sequence,33,52),
             guide2_substr = str_sub(Sequence,84,103)
      )
    GuideR_map <- oligos[,c(4,6)]%>%distinct()%>%
      mutate(guide2_substr = sapply(guide2_substr,revComp))
  }else{
    oligos <- openxlsx::read.xlsx(oligofile,sheet="enCas12a")%>%
      select(Label,Sequence)%>%
      mutate(LabelL = gsub("(.*_L[0-9]{1,2})\\_.*\\_R[0-9]{1,2}","\\1",Label),
             LabelR = gsub(".*_L[0-9]{1,2}\\_(.*\\_R[0-9]{1,2})","\\1",Label),
             guide1_substr = str_sub(Sequence,28,59),
             guide2_substr = str_sub(Sequence,70,103))
    GuideR_map <- oligos[,c(4,6)]%>%distinct()
  }
  GuideL_map <- oligos[,c(3,5)]%>%distinct()
  barcodeinfo <- readRDS("Raw_Data/Barcode_sample_mapinfo.rds")%>%
    filter(Seq==gsub(".*\\/(SN0.*)","\\1",path_to_fastq))%>%
    select(Collapsed_identifier,Barcode)%>%
    set_colnames(c("SampleName","Barcode"))
  file_list%<>%inner_join(barcodeinfo,by="Barcode")
  if(tracrRNA){
    if(Multiple.tracr & Cas=="Cas9"){
      file_list%<>%
        mutate(tracrRNA1 = gsub(tracrRNA_p,"\\1",SampleName)%>%gsub("\\_.*","",.)%>%toupper(),
               tracrRNA2 = paste0("Revcomp_",gsub(tracrRNA_p,"\\1",SampleName)%>%gsub(".*\\_","",.)%>%toupper()))
    }else{
      file_list%<>%
        mutate(tracrRNA1 = "VCR1",
               tracrRNA2 = "Revcomp_WCR3")
    }
  }
  file_list%<>%filter(!grepl("EMPTY",SampleName))
  if(testmode){file_list%<>%.[1:3,]}
  message("Starting mapping reads ......")
  mapFUN <- paste0(Cas,"_map")
  if(Cas=="Cas12a"){
    pbapply(file_list%>%arrange(SampleName),MARGIN = 1, 
            FUN=mapFUN,
            DR1 = "TAATTTCTACTGTCGTAGAT",
            backbone = "GAGACGCGTTTTTTCT",
            lmap = GuideL_map,
            rmap = GuideR_map,
            tmpdir = tmpdir,
            output_path = output_path,
            oligos = oligos,
            cl=mcc)
  }else{
    pbapply(file_list%>%arrange(SampleName),MARGIN = 1, 
            FUN=mapFUN,
            lmap = GuideL_map,
            rmap = GuideR_map,
            tmpdir = tmpdir,
            output_path = output_path,
            oligos = oligos,
            tracrfasta=tracrfasta,
            cl=mcc)
  }
}
#===== Count2LFC ======
counts_pooling <- function(countdir,removed_barcodes,outdir,barcodedir,countpattern,corlim=NULL,plot=T,tracrFlag=F){
  seqBatch <- gsub(".*(SN0[1-9]{6})\\_.*","\\1",countdir)
  barcodeinfo <- read.csv(list.files(barcodedir,pattern = ".csv",full.names = T),
                          stringsAsFactors = FALSE)[,1:2]%>%
    set_colnames(c("Sample.Name","Barcode"))
  counts_pooling <- pbapply::pblapply(dir(countdir,countpattern,full.names = T),
                                      function(s){readRDS(s)%>%
                                          mutate(barcode_lane =gsub(".*\\/.*([ATCG]{8}\\_l[1,2,3,4]{1}).*","\\1",s))%>%
                                          mutate(Barcode =  gsub("\\_.*","",barcode_lane))})%>%
    bind_rows()
  if(tracrFlag){counts_pooling%<>%filter(tracr==F)}
  count_name <- grep("count",colnames(counts_pooling),value=T)
  counts_pooling%<>%
    select(pairs,pair_count = count_name,barcode_lane,Barcode)%>%
    inner_join(barcodeinfo,by="Barcode")%>%
    filter(!Barcode %in% removed_barcodes)%>%
    mutate(Sample.Name = gsub("^[A-Z][0-9]{2}\\_(.*\\_.*\\_.*\\_REP.*)","\\1",Sample.Name))%>%
    group_by(Sample.Name,pairs)%>%
    summarise(pooled_counts  = sum(pair_count))%>%
    as.data.frame()
  if(plot){
    corCounts <- counts_pooling%>%
      spread(Sample.Name,pooled_counts)%>%
      mutate_all(~replace(., is.na(.), 0))%>%
      select(-pairs)%>%
      cor(.)
    if(is.null(corlim)){c(round(min(corCounts),1)-0.05, 1)}
    corCounts%>%heatmaply::heatmaply_cor(.,limits = corlim,
                                         file=paste0(outdir,"/",seqBatch,"_CountsCor.html"))
    
    counts_pooling%>%
      filter(pairs!="AAVS1_L2_AAVS1_R1")%>%
      ggplot(., aes(x=pooled_counts)) + 
      geom_density(alpha=.5,fill="dodgerblue") +
      facet_wrap(. ~ Sample.Name,nrow=3,scales = "free")+
      theme(
        plot.title = element_text(color='black', hjust = 0.5),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.grid = element_blank(),
        axis.text = element_text(color='black'))+
      labs(x="Pooled count of Cas9")+
      ggsave(paste0(outdir,"/",seqBatch,"_pooledCounts.pdf"),height = 14,width=18,dpi=2000)
  }
  counts_pooling%<>%mutate(Seq = seqBatch)
  return(counts_pooling)
}
counts_pooling_wrapper <- function(x,barcodefile){
  seqBatch <- gsub(".*(SN0[0-9]{6})\\_.*","\\1",x[[1]])
  counts_pooling <- pbapply::pblapply(dir(x[[1]],x[[2]],full.names = T),
                                      function(s){readRDS(s)%>%
                                          mutate(barcode_lane =gsub(".*\\/.*([ATCG]{8}\\_l[1,2,3,4]{1}).*","\\1",s))%>%
                                          mutate(Barcode =  gsub("\\_.*","",barcode_lane))})%>%
    bind_rows()%>%
    mutate(Seq = seqBatch)
  if(x[[3]]){counts_pooling%<>%filter(tracr==F)}
  count_name <- grep("count",colnames(counts_pooling),value=T)
  counts_pooling%<>%
    select(Seq,pairs,pair_count = count_name,barcode_lane,Barcode)%>%
    inner_join(barcodefile,by=c("Seq","Barcode"))%>%
    group_by(Collapsed_identifier,cell.line,pairs)%>%
    summarise(pooled_counts  = sum(pair_count))%>%
    as.data.frame()
  return(counts_pooling)
}
lfc_calc <- function(counts,replicate.map){
  message("Calculating LFC from counts and replicate.map ...")
  colnames(replicate.map) <- c("colname","samplename","TP")
  SCALE <- (sum(counts, na.rm = T) / ncol(counts))
  normed_counts <- apply(counts, 2, function(x){
    x = ((x/sum(x, na.rm = T)) * SCALE) + 32
  })
  median_normed_counts <- apply(normed_counts, 2, function(x){
    log2(x) - median(log2(x), na.rm = T)
  })
  ETP.map <- replicate.map%>%filter(TP=="ETP")
  LTP.map <- replicate.map%>%filter(TP=="LTP")
  n_ETP <- ETP.map%>%pull(colname)%>%unique()%>%length()
  if(n_ETP== 0){
    stop("No ETP samples identified in replicate.map")
  }else if(n_ETP == 1){
    ETP = median_normed_counts[,ETP.map$colname]
  }else if(n_ETP>1){
    ETP= median_normed_counts[,ETP.map$colname] %>%as.data.frame()%>%rowMeans(na.rm = T)
  }
  
  LFC <- LTP.map%>%
    group_by(samplename)%>%
    group_map(~ median_normed_counts[,.$colname]%>%
                as.data.frame(optional = T)%>%
                rowMeans(na.rm = T)%>%
                subtract(ETP))%>%
    set_names(sort(unique(LTP.map$samplename)))%>%
    bind_cols()%>%
    as.data.frame()%>%
    set_rownames(rownames(median_normed_counts))%>%
    .[,unique(LTP.map$samplename)]%>%
    as.matrix()
  colnames(LFC) <- unique(LTP.map$samplename)
  return(LFC)
}
#===== Data Fetch ======
nonexpr_Pub21q4 <- function(cl,CCLEpath="tmp/Public21q4_CCLE/",t1=0.1,p=0.9){
  CCLE_RNAseq <- data.table::fread(glue("{CCLEpath}/CCLE_expression_full.csv"))
  sample.info <- taigr::load.from.taiga(data.name='avana-public-21q4-e1b8', data.version=4, data.file='Achilles_sample_info')
  cl_depmap <- sample.info%>%
    filter(CCLE_name %in% cl)%>%
    select(CCLE_name,DepMap_ID)
  CCLE_RNAseq%<>%column_to_rownames("V1")
  colnames(CCLE_RNAseq) <- gsub(" \\(.*","",colnames(CCLE_RNAseq))
  notExpr_CCLE <- apply(CCLE_RNAseq,2,function(s) sum(s<t1)/nrow(CCLE_RNAseq))%>%.[. > p]%>%names()%>%as.character()
  notExpr_allcls <-  apply(CCLE_RNAseq[cl_depmap$DepMap_ID,notExpr_CCLE],2,function(s) sum(s<t1)/nrow(cl_depmap))%>%.[. ==1]%>%names()
  return(notExpr_allcls)
}
DepMap_DataFetch <- function(screenDatFolder=NULL,method="taigr"){
  AvanaScreen_DownloadLink <- list("Achilles_logfold_change" = "https://ndownloader.figshare.com/files/31315903",
                                   "Achilles_guide_map" = "https://ndownloader.figshare.com/files/31315819",
                                   "Achilles_replicate_map" = "https://ndownloader.figshare.com/files/31315876",
                                   "Achilles_sample_info"="https://ndownloader.figshare.com/files/31316011")
  if(method!="taigr"){
    message("Write out download file paths for wget ...")
    writeLines(as.character(AvanaScreen_DownloadLink),glue("{screenDatFolder}/downloadList.txt"))
    message("Using wget to download files ...")
    system(glue("wget -P {screenDatFolder} -i {screenDatFolder}/downloadList.txt"))
    fileNames <- gsub(".*files\\/","",as.character(AvanaScreen_DownloadLink))
    message("Load in data and assign names ...")
    for(s in seq_along(fileNames)){
      assign(names(AvanaScreen_DownloadLink)[s],
             data.table::fread(glue("{screenDatFolder}/",fileNames[s])),
             envir = .GlobalEnv)}
  }else{
    pacman::p_load(taigr)
    message("Use taigr to fetch data and assign names ...")
    # taigr github: https://github.com/broadinstitute/taigr
    # https://cds.team/taiga/dataset/avana-public-21q4-e1b8/4
    for(s in seq_along(AvanaScreen_DownloadLink)){
      assign(names(AvanaScreen_DownloadLink)[s],
             load.from.taiga(data.name='avana-public-21q4-e1b8', 
                             data.version=4,
                             data.file=names(AvanaScreen_DownloadLink)[s]),
             envir = .GlobalEnv)} 
  }
  
}
#===== Screen Data Initialization ======
DiagInit <- function(LFC,posctrl_gene,negctrl_gene,posctrl_pair,annot,cl,ignorePattern=NULL,gene_join=";",nc_gene="AAVS1",mcc=8){
  colnames(annot) <- c("rowname","symbol1","symbol2")
  annot%<>%mutate(pair = pmap_chr(list(symbol1,symbol2),~paste0(sort(c(...)),collapse = gene_join)))
  avg_Single <- pbmcmapply(function(g){
    nc.pairs <- sapply(nc_gene,function(s) paste(sort(c(g,s)),collapse = gene_join))%>%as.character()
    LFC[annot%>%filter(pair %in% nc.pairs)%>%pull(rowname),]%>%apply(.,2,mean,na.rm=T)},
    union(annot$symbol1,annot$symbol2),
    mc.cores = mcc)%>%
    t()%>%as.matrix()
  avg_Double <- pbmcmapply(function(p){
    LFC[annot%>%filter(pair==p)%>%pull(rowname),]%>%apply(.,2,mean,na.rm=T)},
    unique(annot$pair),
    mc.cores = mcc)%>%t()%>%as.matrix()
  pair_map <- data.frame(spairs = rownames(avg_Double),stringsAsFactors = F)%>%
    mutate(g1 = gsub(paste0(gene_join,".*"),"",spairs),
           g2 = gsub(paste0(".*",gene_join),"",spairs))
  Sensitive_score <- (avg_Double-pmin(avg_Single[pair_map$g1,],avg_Single[pair_map$g2,]))%>%
    as.data.frame()%>%
    set_rownames(rownames(avg_Double))%>%
    .[!grepl(nc_gene,rownames(avg_Double)),]
  pair_map_sgRNA <- annot%>%
    filter(symbol1!="AAVS1" & symbol2!="AAVS1")%>%
    mutate(sgRNA1 = gsub("\\_.*","",rowname),sgRNA2 = gsub(".*\\_","",rowname))%>%
    select(rowname,sgRNA1,sgRNA2)%>%
    gather(source,sgRNA,-rowname)%>%
    inner_join(
      annot%>%
        filter(symbol1=="AAVS1" |symbol2=="AAVS1")%>%
        mutate(sgRNA1 = gsub("\\_.*","",rowname),sgRNA2 = gsub(".*\\_","",rowname))%>%
        mutate(sgRNA = case_when(symbol2=="AAVS1"~sgRNA1,TRUE~sgRNA2))%>%
        select(rowname_aavs1 = rowname,sgRNA),by="sgRNA")%>%select(-sgRNA)
  kept_rowname <- pair_map_sgRNA%>%group_by(rowname)%>%summarise(n=n())%>%filter(n==2)%>%pull(rowname)
  pair_map_sgRNA%<>%
    filter(rowname %in% kept_rowname)%>%
    spread(source,rowname_aavs1)
  Sensitive_score_sgRNA <- LFC[pair_map_sgRNA$rowname,]- pmin(LFC[pair_map_sgRNA$sgRNA1,],LFC[pair_map_sgRNA$sgRNA2,])
  nonexpr_gs <- nonexpr_Pub21q4(cl)
  negctrl_pair <- annot%>%
    mutate(pair = pmap_chr(list(symbol1,symbol2),~paste0(sort(c(...)),collapse = ";")))%>%
    filter(symbol1 %in% nonexpr_gs & symbol2 %in% nonexpr_gs)%>%pull(pair)%>%unique()
  if(!is.null(ignorePattern)){negctrl_pair%<>%grep(pattern = ignorePattern,x=.,invert=T,value = T)}
  Diagmodel <- list(LFC=LFC,avg_Single = avg_Single,avg_Double = avg_Double,
                    posctrl_gene=posctrl_gene,negctrl_gene = negctrl_gene,
                    posctrl_pair = posctrl_pair,negctrl_pair = negctrl_pair,
                    annot =annot%>%select(-pair),score=Sensitive_score,
                    score_sgRNA = Sensitive_score_sgRNA)
  return(Diagmodel)
}
#===== Discriminant Analysis ======
Sep_ROC_AUC <- function(sepD){
  colnames(sepD) <- c("value","type","sample")
  pos_ctrl <- sepD%>%filter(type=="pos")%>%pull(sample)
  cut_seq <- seq(ceiling(min(sepD$value))-1,ceiling(max(sepD$value)),by=0.005)
  confusion_seq <- pbmclapply(cut_seq,function(t){
    our_hits <- sepD%>%filter(value <=t)%>%pull(sample)
    confusion_D <- data.frame(P = length(our_hits),
                              N = nrow(sepD) - length(our_hits),
                              TP = length(intersect(pos_ctrl, our_hits)),
                              FN = length(setdiff(pos_ctrl,our_hits)),
                              cutoff = t)%>%
      mutate(FP = P - TP,
             TN = N - FN)},mc.cores = 8)%>%
    bind_rows()%>%
    mutate(TPR =  TP/(TP+FN),
           TNR = TN/(TN+FP),
           FPR = FP/(FP+TN),
           PPV = TP/(TP+FP))
  AUC_ROC <- confusion_seq%>%
    filter(!is.nan(TPR)&!is.nan(FPR))%>%
    summarise(roc_auc =(simple_auc(TPR,FPR)%>%round(.,3)))
  return(list(AUC = AUC_ROC,confusion= confusion_seq))}
Sep_SSMD_NNMD <- function(sepD){
  colnames(sepD) <- c("value","type","sample")
  g1= sepD%>%filter(type=="pos")%>%pull(value)
  g2= sepD%>%filter(type=="neg")%>%pull(value)
  mu1 <- mean(g1, na.rm = T)
  s1 <- var(g1, na.rm = T)
  mu2 <- mean(g2, na.rm = T)
  s2 <- var(g2, na.rm = T)
  ssmd = (mu1 - mu2) / sqrt(s1 +s2)
  # Difference between the medians of the control set,normalized by the median absolute deviation of the negative controls
  nnmd = (median(g1, na.rm = T) - median(g2, na.rm = T)) / mad(g2)
  return(sepm = c(ssmd,nnmd))
}
SepEval <- function(scoreMat,ctrls){
  sep_calc <- lapply(colnames(scoreMat),function(s) {
    subD <- scoreMat%>%
      as.data.frame()%>%
      rownames_to_column("oligo")%>%
      gather(sample,score,-oligo)%>%
      mutate(cl=gsub("\\_.*","",sample))%>%
      inner_join(ctrls%>%mutate(cl=gsub("\\_.*","",cl)),by=c("cl","oligo"))%>%
      filter(sample==s)%>%
      select(score,type,oligo)
    ROC <- Sep_ROC_AUC(subD)
    metrics <- data.frame(metric=c("AUC-ROC","SSMD","NNMD"),
                          values=c(ROC$AUC,Sep_SSMD_NNMD(subD))%>%unlist(),
                          cl = s)%>%
      mutate(values = round(values,3))
    return(list(confusion = ROC$confusion%>%mutate(sample=s),metrics = metrics))})
  metrics_long <- lapply(sep_calc,function(s) s$metrics)%>%bind_rows()
  confusion_long <- lapply(sep_calc,function(s) s$confusion)%>%bind_rows()
  return(list(confusion = confusion_long,metrics = metrics_long))
}
EvalWrapper <- function(scoreList,gene=T){
  if(gene){
    list_sep <- lapply(scoreList,function(l) SepEval(scoreMat = l$LFC,
                                                     lapply(gsub("\\_.*","",colnames(l$LFC))%>%unique(),
                                                            function(s) l$annot%>%
                                                              filter(symbol1=="AAVS1"|symbol2=="AAVS1")%>%
                                                              mutate(Gene = case_when(symbol2=="AAVS1"~symbol1,TRUE~symbol2))%>%
                                                              mutate(type = case_when(Gene %in% l$posctrl_gene~"pos",
                                                                                      Gene %in% l$negctrl_gene~"neg",
                                                                                      TRUE~""))%>%
                                                              filter(type!="")%>%
                                                              mutate(cl=s))%>%
                                                       bind_rows()%>%
                                                       select(oligo=rowname,type,cl)))
  }else{
    list_sep <- lapply(scoreList,function(l){
      pair_sep <- SepEval(scoreMat = l$score_sgRNA,
                          lapply(gsub("\\_.*","",colnames(l$LFC))%>%unique(),
                                 function(s) l$annot%>%
                                   mutate(pair = pmap_chr(list(symbol1,symbol2),~paste0(sort(c(...)),collapse = ";")))%>%
                                   mutate(type = case_when(pair %in% l$posctrl_pair~"pos",
                                                           pair %in% l$negctrl_pair~"neg",
                                                           TRUE~""))%>%
                                   filter(type!="")%>%
                                   mutate(cl=s))%>%
                            bind_rows()%>%
                            select(oligo=rowname,type,cl))
      return(pair_sep)})}
  metrics_long <- lapply(list_sep,function(s) s$metrics)%>%bind_rows()
  confusion_long <- lapply(list_sep,function(s) s$confusion)%>%bind_rows()
  return(list(confusion = confusion_long,metrics = metrics_long))
}
SepBar <- function(sep_metrics,sortm = "AUC-ROC",hj=1,colorp="Dark2"){
  labelLevels <- sep_metrics%>%filter(metric ==sortm)%>%arrange(-values)%>%pull(cl)
  if(colorp=="Paired" & length(labelLevels)<=12){
    mycolors <- brewer.pal(length(labelLevels), colorp)
  }else{mycolors <- colorRampPalette(brewer.pal(8, colorp))(length(labelLevels))}
  mycolors <- mycolors[c(2,1,4,3,6,5,8,7,10,9)[1:length(labelLevels)]]
  barg <- ggplot(sep_metrics%>%mutate(cl =factor(cl,levels = labelLevels))%>%
                   mutate(lp = case_when(metric =="AUC-ROC"~0.03,
                                         metric=="NNMD"~-0.3,
                                         TRUE~-0.1))%>%
                   mutate(lp = values+lp),aes(x=cl,y=values,color=cl,fill=cl)) +
    geom_bar(stat="identity",alpha=0.8)+
    facet_wrap(. ~ metric, scales = 'free_y') + 
    scale_color_manual(values=mycolors,name="")+
    scale_fill_manual(values=mycolors,name="")+
    geom_text(aes(x=as.numeric(cl),y=lp,label=round(values,2)),show.legend = FALSE,family="Times")+
    theme(
      plot.title = element_text(color='black',family = "Helvetica"),
      plot.background = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(color = "black", size = 1, fill = NA),
      panel.grid = element_blank(),
      axis.text.x =  element_text(color='black',family = "Helvetica",angle = 30,hjust = hj,face="bold"),
      axis.text.y = element_text(color='black',family = "Helvetica"),
      axis.title = element_blank(),
      legend.key = element_rect(fill = "white"),
      axis.ticks.x = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(color='black',family = "Helvetica",size = rel(1)),
      legend.position = "none"
    )
  return(barg)
}
SepROC <- function(confusion,labelLevels,colorp="Paired",manucol=NULL,ptitle=""){
  confusion%<>%mutate(colorgrp = factor(sample,levels = labelLevels))
  if(is.null(manucol)){
    if(colorp=="Paired" & length(labelLevels)<=12){
      mycolors <- brewer.pal(length(labelLevels), colorp)
    }else{mycolors <- colorRampPalette(brewer.pal(8, colorp))(length(labelLevels))}
    mycolors <- mycolors[c(2,1,4,3,6,5,8,7,10,9)[1:length(labelLevels)]]
  }else{
    mycolors <- manucol
  }
  roc <- ggplot(data = confusion, aes(x = FPR, y = TPR,color = colorgrp)) +
    geom_segment(data = data.frame(), aes(x = 0, y = 0, xend = 1, yend = 1), inherit.aes = F, color = "grey", linetype = 3) +
    scale_color_manual(values=mycolors, name = "Library (AUC)") +
    geom_path() +
    xlim(0,1) +
    ylim(0,1) +
    geom_line(size=1.5) +
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.text = element_text(color='black',family = "Helvetica"),
          axis.title = element_text(color='black',family = "Helvetica"),
          plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(color = "black", size = 1, fill = NA),
          legend.key = element_rect(fill = "white"),
          legend.position = c(0.99, 0.01),
          legend.justification = c("right", "bottom"),
          legend.box.just = "right",
          legend.margin = margin(5, 5, 5, 5),
          legend.text=element_text(size=rel(1),family = "Helvetica",face="bold"))+ 
    guides(color = guide_legend(order = 1,override.aes = list(size=3,fill=NA)))+
    labs(title=ptitle)
  return(roc)
}
#===== Balance/Align Scatter Plot ======
left_right_comparison <- function(D,annot,posctrl,type="single",FUN="mean"){
  colnames(annot) <- c("rowname","symbol1","symbol2")
  if(type=="single"){
    mapping_12 <-annot%>%
      filter(symbol1=="AAVS1"|symbol2=="AAVS1")%>%
      mutate(Side = case_when(symbol1=="AAVS1"~"side2",TRUE~"side1"),
             symbol = case_when(symbol1=="AAVS1"~symbol2,TRUE~symbol1))%>%
      select(rowname,Side,symbol)%>%distinct()
    desireDesign <- mapping_12%>%group_by(symbol)%>%summarise(n=n())%>%filter(n==6)%>%pull(symbol)
    mapping_12%<>%filter(symbol %in% desireDesign)
  }else{
    mapping_12 <-annot%>%
      filter(!(symbol1=="AAVS1"|symbol2=="AAVS1"))%>%
      mutate(symbol = purrr::pmap_chr(list(symbol1,symbol2),~paste0(sort(c(...)),collapse = ";")))%>%
      mutate(Lead_g = gsub(";.*","",symbol))%>%
      mutate(Side = case_when(Lead_g==symbol1~"side1",TRUE~"side2"))%>%
      select(Side,symbol,rowname)
    desireDesign <- mapping_12%>%group_by(symbol)%>%summarise(n=n())%>%filter(n==18)%>%pull(symbol)
    mapping_12%<>%filter(symbol %in% desireDesign)
  }
  compare_data <- lapply(colnames(D),function(cl){
    D[,cl]%>%
      as.data.frame()%>%
      set_colnames(cl)%>%
      mutate(rowname = rownames(D))%>%
      inner_join(mapping_12,by="rowname")%>%
      group_by(symbol,Side)%>%
      summarise_at(vars(cl),FUN,na.rm=T)%>%
      spread(Side,cl)%>%
      na.omit()%>%
      mutate(Essentiality = case_when(symbol %in% posctrl~"yes",TRUE~"no"))%>%
      mutate(source = cl)})%>%bind_rows()
  return(compare_data)
}
left_right_scatter <- function(inputlist,type,FUN,fillColor="skyblue3",rowno=2,comb=F,cor.ordered=T,orders=NULL){
  if(type=="single"){
    lr_compare<- lapply(inputlist,function(s) left_right_comparison(D =s$LFC,annot = s$annot,posctrl = s$posctrl_gene,type=type,FUN=FUN)) 
  }else{
    lr_compare<- lapply(inputlist,function(s) left_right_comparison(D =s$LFC,annot = s$annot,posctrl = s$posctrl_pair,type=type,FUN=FUN))  
  }
  
  lr_compare%<>%bind_rows()%>%as.data.frame()%>%filter(source %in% orders)
  min_lim <- floor(apply(lr_compare%>%select(side1,side2),2,min,na.rm=T)%>%min())
  max_lim <- ceiling(apply(lr_compare%>%select(side1,side2),2,max,na.rm=T)%>%max())
  lims <- c(min_lim,max_lim)
  tailedX <-ifelse(type=="single",
                   "LFC( <span style='color:#009E73;'>Left</span> sgRNA - AAVS1)</span>",
                   "LFC( <span style='color:#009E73;'>sgGeneX</span> : <span style='color:#0072B2;'>sgGeneY</span>)</span>")
  tailedY <-ifelse(type=="single",
                   "LFC( AAVS1 - <span style='color:#0072B2;'>Right</span> sgRNA)</span>",
                   "LFC( <span style='color:#009E73;'>sgGeneY</span> : <span style='color:#0072B2;'>sgGeneX</span>)</span>")
  corsD <- lr_compare%>%
    group_by(source)%>%
    summarize(values = cor(side1, side2, use = "pairwise.complete",method="pearson"))%>%
    as.data.frame()%>%
    mutate(cors =  paste("r =", signif(values, 2)))%>%
    arrange(values)
  if(cor.ordered){
    lr_compare%<>%mutate(source=factor(source,levels = 
                                         (corsD%>%arrange(values)%>%pull(source)),
                                       ordered = TRUE))
    corsD%<>%mutate(
      source = factor(source,levels = 
                        (corsD%>%arrange(values)%>%pull(source)),
                      ordered = TRUE))
  }else{
    lr_compare%<>%mutate(source=factor(source,levels = 
                                         orders,
                                       ordered = TRUE))
    corsD%<>%mutate(
      source = factor(source,levels = orders,
                      ordered = TRUE))
  }
  
  
  p <- ggplot(lr_compare, aes(x=side1, y=side2,color=Essentiality)) +
    geom_rect(aes(xmin = -Inf, xmax = -1, ymin = -Inf, ymax = -1),
              fill = "#fbf7f7",linetype = "longdash",color = "tomato4")+
    geom_abline(intercept = 0)+ 
    geom_point(size=1.3,data = .%>%filter(Essentiality=="no"), alpha=0.8) +
    geom_point(size = 1.3,data = .%>%filter(Essentiality=="yes"), alpha=0.8)+
    #stat_smooth(color = "tomato4", method = 'lm', se = FALSE,linetype = "longdash") + 
    scale_color_manual(values=c("grey50","red3"),labels = c( "Others","Pan-essential genes"))+
    scale_x_continuous(limits=lims) + scale_y_continuous(limits=lims)+
    theme(plot.title = element_text(color='black', hjust = 0.5),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", size = 0.5, fill = NA),
          panel.grid = element_blank(),
          axis.text = element_text(color='black'),
          legend.key = element_rect(fill = NA),
          legend.position ="none",
          axis.title.x = element_markdown(),
          axis.title.y = element_markdown())+
    #xlim(lims)+ylim(lims)+
    guides(colour = guide_legend(override.aes = list(size=3)))+
    labs(x=tailedX,
         y = tailedY)+
    geom_hline(yintercept  = 0,linetype = "longdash")+ 
    geom_vline(xintercept  = 0,linetype = "longdash")
  if(comb){
    p <- p+
      geom_text(aes(x = min_lim+1.5, y = max_lim-0.5, color = NULL,label=
                      paste("r =", signif( cor(side1, side2, use = "pairwise.complete",method="pearson"), 2))),
                show.legend = FALSE,size=4.5,inherit.aes = F,family="Times")
  }else{
    p <- p+
      facet_wrap(~source,nrow=rowno,scales = "free")+
      geom_text(data = corsD,aes(x = min_lim+1.5, y = max_lim-0.5, color = NULL,label=cors),
                show.legend = FALSE,size=4.5,inherit.aes = F,family="Times")+
      theme(strip.text.x = element_text(
        size = 12, color = "white", family = "Times"),
        strip.background = element_rect(
          color="black", fill=fillColor, size=1, linetype="solid"))
  }
  return(p)}
scatterP_minilib <- function(d,xVAR,yVAR,cVAR,outdir,i="IPC298-VCR1-WCR3",
                             xlab="Average LFC in Achilles Avana",
                             ylab="Average LFC for Dual-Targetting sgRNAs",
                             Prefix="DualTarget_Avana",pt=1.8){
  Dat <- Sys.Date()
  plotD <- d%>%filter(condition==i)%>%select(xvar=all_of(xVAR),yvar=all_of(yVAR),cvar=all_of(cVAR))
  p <- ggplot(plotD, aes(x=xvar, y=yvar,color=cvar,shape=cvar)) +
    geom_rect(aes(xmin = -Inf, xmax = -1, 
                  ymin = -Inf, ymax = -1),
              fill = "#e2c4c4",color="black")+
    geom_point(size=pt,data = .%>%filter(cvar=="Non-essential"), alpha=0.8) +
    geom_point(size = pt,data = .%>%filter(cvar=="Pan-essential"), alpha=0.8)+
    geom_point(size = pt,data = .%>%filter(cvar=="Selective-essential"), alpha=0.8)+
    scale_color_manual(values=c("grey50","red3","deepskyblue4"),
                       labels = c( "Non-essential genes","Pan-essential genes",
                                   "Selective-essential genes"))+
    geom_abline(intercept = 0,linetype = "longdash")+ 
    geom_hline(yintercept  = 0,linetype = "longdash")+ 
    geom_vline(xintercept  = 0,linetype = "longdash")+ 
    theme(plot.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", size = 1, fill = NA),
          panel.grid = element_blank(),
          axis.text = element_text(color='black'),
          strip.background = element_blank(),
          legend.key = element_rect(fill = NA),
          legend.position ="none",
          axis.title = element_text(color='black'))+
    xlim(c(-6,0.8))+
    ylim(c(-6,0.8))+
    guides(colour = guide_legend(override.aes = list(size=3)))+
    labs(x=xlab,y = ylab,title = i)
  ggsave(glue("{outdir}/{Prefix}_Align_{i}_{Dat}.pdf"),width=5,height=5)
  return(p)
}
align_avana <- function(inputlist,rowno=2,cor.ordered=T,orders=NULL,clSingle=NULL,dataLoad=F){
  cols <- sapply(inputlist,function(s) colnames(s$LFC))%>%unlist()
  cls <- gsub("\\_.*","",cols)%>%unique()
  if(dataLoad){DepMap_DataFetch()}
  sample.sub <- Achilles_sample_info%>%
    mutate(stripped_cell_line_name = gsub("\\_.*","",CCLE_name))%>%
    filter(stripped_cell_line_name %in% cls)
  avana.sub <- Achilles_guide_map%>%
    mutate(gene =gsub(" \\(.*","",gene))%>%
    select(sgrna,gene)
  lfc.cls <- Achilles_replicate_map%>%filter(DepMap_ID %in% sample.sub$DepMap_ID)
  lfc.avana <- Achilles_logfold_change[,lfc.cls$replicate_ID]%>%
    as.data.frame()%>%
    rownames_to_column("sgrna")%>%
    gather(replicate_ID,lfc,-sgrna)%>%
    inner_join(lfc.cls%>%select(replicate_ID,DepMap_ID),by="replicate_ID")%>%
    inner_join(avana.sub,by="sgrna")%>%
    group_by(gene,DepMap_ID)%>%
    summarise(lfcAvana = mean(lfc))%>%
    as.data.frame()%>%
    inner_join(sample.sub%>%select(CCLE=stripped_cell_line_name,DepMap_ID),by="DepMap_ID")
  single_avgs <- lapply(inputlist,function(s) s$avg_Single%>%as.data.frame()%>%rownames_to_column("gene")%>%
                          gather(source,y,-gene)%>%
                          mutate(Essentiality = case_when(gene %in% s$posctrl_gene~"yes",TRUE~"no")))%>%bind_rows()
  single_avgs%<>%filter(source %in% orders)
  singleNull <- length(intersect(gsub("\\_.*","",single_avgs$source),cls))
  if(singleNull==0){
    single_avgs%<>%mutate(cl = clSingle)
  }else{
    single_avgs%<>%mutate(cl = gsub("\\_.*","",source)) 
  }
  avana_combined <- single_avgs%>%inner_join(lfc.avana,by="gene")%>%filter(cl==gsub("\\_.*","",CCLE))
  min_lim <- floor(apply(avana_combined%>%select(y,lfcAvana),2,min,na.rm=T)%>%min())
  max_lim <- ceiling(apply(avana_combined%>%select(y,lfcAvana),2,max,na.rm=T)%>%max())
  lims <- c(min_lim,max_lim)
  pearson_cors <- avana_combined%>%
    group_by(source)%>%
    summarize(cors =  paste("r =", signif(cor(y, lfcAvana, use = "pairwise.complete"), 3)))
  
  if(cor.ordered){
    pearson_cors%<>%mutate(
      source = factor(source,levels = 
                        (pearson_cors%>%arrange(cors)%>%pull(source)),
                      ordered = TRUE))
    avana_combined%<>%mutate(source=factor(source,levels = 
                                             (pearson_cors%>%arrange(cors)%>%pull(source)),
                                           ordered = TRUE))
  }else{
    pearson_cors%<>%mutate(
      source = factor(source,levels = orders,
                      ordered = TRUE))
    avana_combined%<>%mutate(source=factor(source,levels = orders,
                                           ordered = TRUE))
  }
  scatter_avana <- ggplot(avana_combined, aes(x=lfcAvana, y=y,color=Essentiality)) +
    geom_rect(aes(xmin = -Inf, xmax = -1, ymin = -Inf, ymax = -1),
              fill = "#fbf7f7",linetype = "longdash",color = "tomato4")+
    geom_abline(intercept = 0)+ 
    geom_point(size=1.3,data = .%>%filter(Essentiality=="no"), alpha=0.8) +
    geom_point(size = 1.3,data = .%>%filter(Essentiality=="yes"), alpha=0.8)+
    facet_wrap(~source,nrow=rowno,scales = "free")+
    scale_color_manual(values=c("grey50","red3"),labels = c( "Others","Pan-essential genes"))+
    geom_text(data = pearson_cors,aes(x = -4.8, y =1.8, color = NULL,label=cors),show.legend = FALSE,
              size=4.5,inherit.aes = FALSE,family="Times")+
    theme(plot.title = element_text(color='black', hjust = 0.5),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", size = 0.5, fill = NA),
          panel.grid = element_blank(),
          axis.text = element_text(color='black'),
          legend.key = element_rect(fill = NA),
          legend.position ="none"
    )+
    scale_x_continuous(limits=lims) + scale_y_continuous(limits=lims)+
    guides(colour = guide_legend(override.aes = list(size=3)))+
    labs(x="LFC (Avana)", y = "LFC (Single gene - AAVS1)")+
    theme(
      strip.background = element_blank(),
      strip.text.x = element_text(hjust = 0,colour = "black"))+
    geom_hline(yintercept  = 0,linetype = "longdash")+ 
    geom_vline(xintercept  = 0,linetype = "longdash")
  return(scatter_avana)
}
#===== Manhattan Plot =====
scores_SigTest <- function(pairLFC,score,nc_pairs){
  message("Calculating p values and FDR ...")
  ncp_flag <- !is.null(nc_pairs)
  if (ncp_flag) {
    if (!requireNamespace("mixtools")) {
      stop("Please install the mixtools package: install.packages('mixtools').")
    }
    # lethality
    nc_values <- unlist(score[nc_pairs, ])%>%na.omit()
    invisible(capture.output({nc_model_lethality <- mixtools::normalmixEM(nc_values)}))
    # calculate p-values
    pval_sensi_lethal <- nc_model_lethality$lambda[1] * pnorm(
      score,
      mean = nc_model_lethality$mu[1],
      sd = nc_model_lethality$sigma[1],
      lower.tail = FALSE
    ) +
      nc_model_lethality$lambda[2] * pnorm(
        score,
        mean = nc_model_lethality$mu[2],
        sd = nc_model_lethality$sigma[2],
        lower.tail = FALSE
      ) # null hypothesis: nc_model_lethality > score_sensi_lethal
  }
  if (ncp_flag) {fdr_sensi_lethal <- apply(pval_sensi_lethal, 2, function(x) p.adjust(x, method = "fdr"))}
  score_df = -log10(fdr_sensi_lethal) %>%
    as.data.frame() %>%
    rownames_to_column('genepair') %>%
    gather(sample, GEMINI.score, -genepair)
  lfc_df = pairLFC %>%
    as.data.frame() %>%
    rownames_to_column('genepair') %>%
    gather(sample, LFC, -genepair)
  plot_df <- score_df %>% 
    inner_join(lfc_df, by = c("genepair", "sample")) %>%
    filter(complete.cases(.)) %>%
    arrange(genepair, sample)
  return(plot_df)
}
match_matrix <- function(matrixList,cols){
  matrixList%<>%lapply(.,as.matrix)
  rowName_shared <- Reduce(intersect,lapply(matrixList,rownames))
  reduced_mat <- Reduce(cbind,lapply(matrixList,function(s) s[rowName_shared,]))
  oMat <- reduced_mat[,cols]
  return(oMat)
} 
#===== Boxplot Stack =====
boxplot_stack <- function(LFC,map,gene_pairs,cell_line,nc_gene="AAVS1",gene_join=":"){
  map%<>%mutate(pair = pmap_chr(list(symbol1,symbol2),~paste0(sort(c(...)),collapse = gene_join)))
  comb_D <- lapply(gene_pairs,function(p){
    ls <- c(strsplit(p,split=gene_join)[[1]],p)
    g <- gsub(paste0(gene_join,".*"),"",p)
    h <- gsub(paste0(".*",gene_join),"",p)
    nc_pairs <- expand.grid(nc_gene,c(g,h),stringsAsFactors = F)%>%
      mutate(spair = pmap_chr(list(Var1,Var2),~paste0(sort(c(...)),collapse = gene_join)))%>%
      pull(spair)
    lfc_boxplot <- data.frame(LFC[map%>%filter(pair %in% c(nc_pairs,p))%>%pull(rowname),cell_line])%>%
        rownames_to_column("rowname")%>%
        gather(sample,lfc_cl,-rowname)%>%
        inner_join(map%>%select(rowname,symbol1,symbol2,pair),by="rowname")%>%
        mutate(pair=gsub(paste0(paste0(gene_join,nc_gene),"|",paste0(nc_gene,gene_join)),"",pair))%>%
        mutate(spairs=p)%>%
        mutate(ind = factor(case_when(pair==ls[1]~1,pair==ls[2]~2,TRUE~3),levels=c(3,2,1)))
  })%>%bind_rows()
  all_levels <- lapply(gene_pairs,function(p) c(strsplit(p,split=gene_join)[[1]],p))%>%unlist()%>%unique()
  r_levels <- all_levels[!all_levels%in% c("VPS4A:VPS4B","VPS4A","VPS4B","KATNAL2","KATNAL2:VPS4A")]
  all_levels <- c(r_levels,"KATNAL2:VPS4A","VPS4A","VPS4B","KATNAL2","VPS4A:VPS4B")
  comb_D%<>%mutate(pair = factor(pair,levels=all_levels))
  comb_D%<>%mutate(sample = gsub(".*\\_","",sample))%>%
    mutate(sample =factor(sample,levels=c("enCas12a","WCR2.WCR3","VCR1.WCR3")))
  p <- ggplot(comb_D,aes(x = pair,y = lfc_cl,fill=ind))+
    facet_grid(sample~spairs, scales = "free")+
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0),
              fill = "#fbf7f7",linetype = "longdash",color = "tomato4")+
    geom_boxplot(outlier.shape = NA)+
    scale_fill_manual(values=c("#E69F00", "#666666","#999999"))+
    geom_jitter(size=1.5,width = 0.3,aes(colour=ind),shape=18)+
    scale_color_manual(values=c("#999999", "#E69F00","#E69F00"))+
    theme(plot.title = element_text(color='black', hjust = 0.5),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", size = 1, fill = NA),
          text = element_text(color='black'),
          panel.grid = element_blank(),
          axis.text.y = element_text(color='black',family = "Times"),
          axis.text.x = element_text(color='black',family = "Helvetica",angle=90,hjust=0.95,vjust=0.2),
          strip.background = element_blank(),
          axis.title.x  = element_blank(),
          legend.position = "none",
          strip.text.x  = element_blank()
    )+
    labs(y="LFC")+
    geom_hline(yintercept = 0,linetype="dashed",color="tomato3")+
    ylim(-6, 2)
  return(p)
}