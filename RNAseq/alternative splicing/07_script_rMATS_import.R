suppressPackageStartupMessages({
  library(readxl)
  library(tidyverse)
  library(openxlsx)
  library(biomaRt)
})

### see https://www.biostars.org/p/256949/ for rMATS headers explainations ####

### read in JCEC  reads that span splicing junctions and reads on target (striped regions on home page figure) ####
### read in JC only reads that span splicing junctions ####
parse_rmats <- function(z){
  filenames <- sort(system(paste0("ls"," results/",z,"/*JCEC.txt.gz"), intern=T))
  sstypes <- c("A3SS_JCEC","A5SS_JCEC","MXE_JCEC","RI_JCEC","SE_JCEC")
  read_rmats <- function(x,y){
    read.table(gzfile(x), sep="\t",header=T) %>% arrange(FDR) %>%
      separate(IJC_SAMPLE_1, c("IJC_SAMPLE1.1","IJC_SAMPLE1.2","IJC_SAMPLE1.3"), sep=",", convert=T) %>%
      separate(SJC_SAMPLE_1, c("SJC_SAMPLE1.1","SJC_SAMPLE1.2","SJC_SAMPLE1.3"), sep=",", convert=T) %>%
      separate(IJC_SAMPLE_2, c("IJC_SAMPLE2.1","IJC_SAMPLE2.2","IJC_SAMPLE2.3"), sep=",", convert=T) %>%
      separate(SJC_SAMPLE_2, c("SJC_SAMPLE2.1","SJC_SAMPLE2.2","SJC_SAMPLE2.3"), sep=",", convert=T) %>%
      separate(IncLevel1, c("IncLevel1.1","IncLevel1.2","IncLevel1.3"), sep=",", convert=T) %>%
      separate(IncLevel2, c("IncLevel2.1","IncLevel2.2","IncLevel2.3"), sep=",", convert=T) %>%
      mutate("startsite"=.[,6]-250, "endsite"=.[,7]+250) %>%
      mutate("ROI" = paste(chr,":",startsite,"-",endsite,":",strand,sep=""))
  }
  temp <- lapply(filenames, read_rmats)
  names(temp) <- paste0(z,sstypes)
  list2env(temp, envir = .GlobalEnv)
  rm(temp)
  
  ### make a DF for only sig events psi >0.05, fdr <0.05, avg count >10 in either group ####
  
  sig_data <- function(x){
    filter(x,FDR<0.05 & abs(IncLevelDifference)>0.05) %>% filter(rowMeans(cbind(IJC_SAMPLE1.1,IJC_SAMPLE1.2))>10|rowMeans(cbind(IJC_SAMPLE2.1,IJC_SAMPLE2.2))>10)
  }
  
  for (i in paste0(z,sstypes)){
    sstype <- get(i)
    tmp <- sig_data(sstype)
    assign(paste(i,"_sig", sep = "") , tmp, envir = globalenv())
  }
  
  ### SUMMARIZE ####
  
  ### write output ####  sample1=WT, sample2=MUT ####
  
  temp <- list("A3SS_JCEC"=get(paste0(z,"A3SS_JCEC")),
               "A5SS_JCEC"=get(paste0(z,"A5SS_JCEC")),
               "MXE_JCEC"=get(paste0(z,"MXE_JCEC")),
               "RI_JCEC"=get(paste0(z,"RI_JCEC")),
               "SE_JCEC"=get(paste0(z,"SE_JCEC")))
  openxlsx::write.xlsx(temp, file=paste0("20220304_",z,"_rMATS_JCEC.xlsx"))
  
  tempsig <- list("A3SS_JCEC_sig"=get(paste0(z,"A3SS_JCEC_sig")),
               "A5SS_JCEC_sig"=get(paste0(z,"A5SS_JCEC_sig")),
               "MXE_JCEC_sig"=get(paste0(z,"MXE_JCEC_sig")),
               "RI_JCEC_sig"=get(paste0(z,"RI_JCEC_sig")),
               "SE_JCEC_sig"=get(paste0(z,"SE_JCEC_sig")))
  openxlsx::write.xlsx(tempsig, file=paste0("20220304_",z,"_rMATS_JCEC_sig.xlsx"))
}




###########

