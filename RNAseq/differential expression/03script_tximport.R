suppressPackageStartupMessages({
  library(readxl)
  library(tidyverse)
  library(openxlsx)
  library(tximport)
  library(DESeq2)
  library(pheatmap)
  library(biomaRt)
})

if (!file.exists("biomart.rds")){
  mart <- useMart("ENSEMBL_MART_ENSEMBL")
  mart <- useDataset("mmusculus_gene_ensembl", mart)
  atributes <- c('ensembl_gene_id_version','ensembl_transcript_id_version', 'description','gene_biotype',"external_gene_name")

  ens_annot <- getBM(attributes=atributes, mart=mart)
  ens_annot$external_gene_name[ens_annot$external_gene_name==""] <- NA
  ens_annot$description[ens_annot$description==""] <- NA
  ens_annot[is.na(ens_annot$external_gene_name),]$external_gene_name <- ens_annot[is.na(ens_annot$external_gene_name),]$ensembl_gene_id_version
  saveRDS(ens_annot, "biomart.rds")
  
  colnames(ens_annot) <- c("GENEID","TXID","DESC","GENETYPE","GENENAME")
  
  gene_desc <- ens_annot[!grepl(c("Mt_|miRNA|rRNA|scRNA|snRNA|snoRNA|sRNA|scaRNA|vaultRNA"), ens_annot$GENETYPE),] %>%
    .[,c(3,5)] %>%
    separate(DESC, c("DESC","null"),sep="\\[") %>%
    dplyr::select(1,3) %>%
    unique()
  saveRDS(gene_desc, "biomart_desc.rds")
} else {
  ens_annot <- readRDS("biomart.rds")
  gene_desc <- readRDS("biomart_desc.rds")

}

if (!file.exists("txi.rds")){
  path <- "results/salmon" # replace with folder dir
  files <- list.files(path, pattern = "S")
  files_full <- paste0(path, "/", files, "/quant.sf.gz")
  names(files_full) <- files
  TX2Gene <- ens_annot[,c(2,5)]
  txi <- tximport(files_full, type = "salmon", tx2gene = ens_annot[,c(2,5)])
  saveRDS(txi, "txi.rds")
} else {
  txi <- readRDS("txi.rds")
}

if (!file.exists("txi_tx.rds")){
  path <- "results/salmon" # replace with folder dir
  files <- list.files(path, pattern = "S")
  files_full <- paste0(path, "/", files, "/quant.sf.gz")
  names(files_full) <- files
  TX2Gene <- ens_annot[,c(2,5)]
  txi <- tximport(files_full, type = "salmon", tx2gene = ens_annot[,c(2,2)])
  saveRDS(txi, "txi_tx.rds")
} else {
  txi <- readRDS("txi.rds")
}

counts <- data.frame(txi$counts) %>% rownames_to_column(var = "tx") %>%
  left_join(TX2Gene, by = c("tx" = "ensembl_transcript_id_version"))
