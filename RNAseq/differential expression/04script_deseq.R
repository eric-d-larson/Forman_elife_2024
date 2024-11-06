suppressPackageStartupMessages({
  library(readxl)
  library(tidyverse)
  library(openxlsx)
  library(tximport)
  library(DESeq2)
  library(pheatmap)
})

source("script_tximport.R")

meta <- read.xlsx("meta_data.xlsx")

if (!file.exists("dds.filtered.rds")){
  dds.raw <- DESeqDataSetFromTximport(txi,
                                colData = meta,
                                design = ~ genotype)

# filter based on average count of >= 5 in at least 1 group
keep <- rowSums(subset(counts(dds.raw), select=c(grep("\\bctrl\\b", meta$genotype)))) >= 10 | 
  rowSums(subset(counts(dds.raw), select=c(grep("\\btreated\\b", meta$genotype)))) >= 10
dds.filtered <- dds.raw[keep,]
saveRDS(dds.filtered, "dds.filtered.rds")
} else {
  dds.filtered <- readRDS("dds.filtered.rds")
}

#distance.m_rlog <- as.dist(1-cor(assay(rld), method="pearson")) 
#plot(hclust(distance.m_rlog), label=colnames(assay(rld)),main="rlog transformed read counts\ndistance: Pearson correlation") 
#df=as.data.frame(colData(dds)[,c("genotype")])
row.names(df) <- meta$lib
colnames(df) <- "genotype"
pheatmap(1-cor(assay(rld), method="pearson"),
         cluster_rows = T,
         cluster_cols = T,
         show_colnames = T,
         display_numbers = F,
         annotation_col = df)
### rlog transform the matrix

if (!file.exists("rld.rds")) {
  rld <- rlog(dds.filtered, blind=F)
  saveRDS(rld, "rld.rds")
} else {
  rld <- readRDS("rld.rds")
}

pcaData <- plotPCA(rld, intgroup=c("genotype"), returnData=TRUE, ntop = 500)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=group, label=name)) +
  geom_point(size=3) +
  geom_text(size=0) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme(legend.position = "top")


### make deseq

if (!file.exists("dds.rds")){
  dds <- DESeq(dds.filtered)
  saveRDS(dds, "dds.rds")
} else {
  dds <- readRDS("dds.rds")
}


### extact normalized count averages per genotype
if (!file.exists("norm_counts_avg.rds")){
  norm_counts_avg <- as.data.frame(counts(dds, normalized=T)) %>% 
    signif(3) %>%
    mutate("gene"=row.names(.)) %>% 
    dplyr::mutate("ctrl_avg"=rowMeans(subset(.,select=c(grep("\\bctrl\\b", meta$genotype)))),
                  "treated_avg"=rowMeans(subset(.,select=c(grep("\\btreated\\b", meta$genotype))))) %>% 
    remove_rownames() %>%
    column_to_rownames(var="gene") %>%
    signif(3) %>%
    rownames_to_column(var="gene") %>%
    left_join(readRDS("biomart_desc.rds"), by=c("gene"="GENENAME")) %>%
    dplyr::select(gene,DESC,grep('avg',colnames(.)))
  saveRDS(norm_counts_avg, "norm_counts_avg.rds")
} else {
  norm_counts_avg <- readRDS("norm_counts_avg.rds")
}


if (!file.exists("deseq.rds")){
  res <- results(dds,
                  alpha = 0.05, 
                  contrast = c("genotype","ctrl","treated")) %>%
    data.frame() %>%
    signif(3) %>%
    rownames_to_column(var="gene") %>% 
    as_tibble() %>% 
    arrange(padj) %>%
    left_join(norm_counts_avg) %>%
    dplyr::select(1,8:10, 3,7)
  write.xlsx(list("all_genes" = res,
                  "sig_up" = dplyr::filter(res, padj < 0.05 & log2FoldChange >0),
                  "sig_down" = dplyr::filter(res, padj < 0.05 & log2FoldChange <0)),
                  file="deseq.xlsx")
  saveRDS(res, "deseq.rds")
  res_shrink <- lfcShrink(dds, coef = "genotype_treated_vs_ctrl", type="apeglm")
  saveRDS(res_shrink, file="res_shrink.rds")
} else {
  res <- readRDS("deseq.rds")
  res_shrink <- readRDS("res_shrink.rds")
}
