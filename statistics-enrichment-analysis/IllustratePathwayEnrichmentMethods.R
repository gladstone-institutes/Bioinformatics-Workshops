rm(list = ls())
require(GSEABenchmarkeR)
require(clusterProfiler)
require(DESeq2)
require(tidyverse)
require(rSEA)
# tcga <- list()
tcga <- readRDS("bladder_cancer_tcga_summarized_experiment.rds")
#object of class summarized experiment
##phenotype data
colData(tcga)
tcga$GROUP <- as.factor(tcga$GROUP)
tcga$BLOCK <- as.factor(tcga$BLOCK)


##access count data
##FILL-IN

##create a DESeq data object
dds.bc <- DESeqDataSet(tcga, design = ~ GROUP + BLOCK)

##estimate normalization/size-factors and dispersions
dds.bc %<>% DESeq(.)


##variance stabilizing transformation to view the normalize data
vsd.bc <- dds.bc %>%
  vst(., blind=TRUE)

##generate the PCA plot
vsd.bc %>%
  plotPCA(., intgroup=c("GROUP"))

##differential expression association
diff.res <- dds.bc %>% 
  results(., contrast = c("GROUP", "1", "0"), pAdjustMethod="bonferroni")

##generate MA plot
diff.res %>%
  plotMA(.)

##output the results
diff.res %>%
  as.data.frame() %>%
  rownames_to_column('Gene') %>%
  write.csv(., "bladder_cancer_diff_exp_results.csv", row.names = FALSE)
##load the pathway gene set data-bases
database_lists <- load("databases.RData")#has wp, pfocr, go

##Check out WikiPathways annotation
head(wp_annotation)

##Check out WikiPathways list
head(wp_list)

##Run ORA analyses
##Choose set of differential expressed genes
##pick the differentially expressed genes using 0.05 threshold
diff_genes <- diff.res %>%
    as.data.frame() %>%
    rownames_to_column('gene') %>%
    filter(padj < 0.05) %>%
    .$gene 
##important to pick the universe of genes
universe_genes <- diff.res %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  .$gene

##run the ORA analyses
res_ora <- enricher(
  gene = diff_genes,
  universe = universe_genes,
  pAdjustMethod = "BH",
  pvalueCutoff = 1, #p.adjust cutoff
  qvalueCutoff = 1,
  minGSSize = 1,
  maxGSSize = 100000,
  TERM2GENE = wp_annotation[,c("set_id","gene")],
  TERM2NAME = wp_annotation[,c("set_id","name")])

res_ora <- res_ora@result

##Estimate the odds ratio
k <- sapply(res_ora$GeneRatio, function(x) as.numeric(strsplit(x, "/")[[1]][1]))
n <- sapply(res_ora$GeneRatio, function(x) as.numeric(strsplit(x, "/")[[1]][2]))
M <- sapply(res_ora$BgRatio, function(x) as.numeric(strsplit(x, "/")[[1]][1]))
N <- sapply(res_ora$BgRatio, function(x) as.numeric(strsplit(x, "/")[[1]][2]))
odds_ratio <- (k*(N-M-n+k))/((M-k)*(n-k))

res_ora %<>% mutate(odds_ratio=odds_ratio)
res_ora %>%
  write.csv(., "bladder_cancer_WikiPathways_ora.csv", row.names = FALSE)

##run GSEA: gene permutation
gene_list <- diff.res %>%
  as.data.frame() %>%
  rownames_to_column('Gene') %>%
  mutate(Score = sign(as.numeric(log2FoldChange)) * - log10(as.numeric(as.character(pvalue)))) %>%
  select(c("Score","Gene")) %>%
  arrange(desc(Score))

gene_list <- unlist(split(gene_list[, 1], gene_list[, 2]))
gene_list = sort(gene_list[unique(names(gene_list))], decreasing = TRUE)

res_gsea_gene_perm <- GSEA(
  gene_list,
  pAdjustMethod="BH",
  TERM2GENE = wp_annotation[,c("set_id","gene")],
  TERM2NAME = wp_annotation[,c("set_id","name")]    ,
  minGSSize = 1,
  maxGSSize = 100000,
  pvalueCutoff = 1,
  verbose=FALSE)

res_gsea_gene_perm <- res_gsea_gene_perm@result
View(res_gsea_gene_perm)
res_gsea_gene_perm %>%
  write.csv(., "bladder_cancer_WikiPathways_gsea_gene_perm.csv", row.names = FALSE)

##run GSEA sample label permutation
#library("devtools")
#install_github("GSEA-MSigDB/GSEA_R")
tcga.de <- readRDS("bladder_cancer_tcga_summarized_experiment_w_de_results.rds")
res_gsea_sample_perm <- runEA(tcga.de, method="gsea", gs=wp_list, perm=1000)
res_gsea <- res_gsea_sample_perm$gsea[[1]]$ranking %>% as.data.frame()
res_gsea %<>% mutate(set_id=GENE.SET)
res_gsea %<>% merge(wp_id_2_names,.) %>% slice(order(PVAL))
res_gsea %>%
  write.csv(., "bladder_cancer_WikiPathways_gsea_sample_perm.csv", row.names = FALSE)

res_safe_sample_perm <- runEA(tcga.de, method="safe", gs=wp_list, perm=1000)
res_safe <- res_safe_sample_perm$safe[[1]]$ranking %>% as.data.frame()
res_safe %<>% mutate(set_id=GENE.SET)
res_safe %<>% merge(wp_id_2_names,.) %>% slice(order(PVAL))
res_safe %>%
  write.csv(., "bladder_cancer_WikiPathways_safe_sample_perm.csv", row.names = FALSE)

res_padog_sample_perm <- runEA(tcga.de, method="padog", gs=wp_list, perm=1000)
res_padog <- res_padog_sample_perm$padog[[1]]$ranking %>% as.data.frame()
res_padog %<>% mutate(set_id=GENE.SET)
res_padog %<>% merge(wp_id_2_names,.) %>% slice(order(PVAL))
res_padog %>%
  write.csv(., "bladder_cancer_WikiPathways_padog_sample_perm.csv", row.names = FALSE)

##run rSEA method
res_rSEA <- SEA(diff.res$pvalue, universe_genes, pathlist = wp_list)

res_rSEA %<>% mutate(set_id=Name)
##get pathway names
wp_id_2_names <- wp_annotation %>%
  select(1,2) %>%
  unique()

res_rSEA %<>% merge(wp_id_2_names,.)
res_rSEA %>% slice(order(Comp.adjP)) %>%
  write.csv(., "bladder_cancer_WikiPathways_rSEA.csv", row.names = FALSE)

TDPestimate_full <- setTDP(diff.res$pvalue, universe_genes, alpha = 0.05)

##generate plot of counts of number of associated pathways
res_gsea <- read.csv("bladder_cancer_WikiPathways_gsea_sample_perm.csv", header = TRUE)
res_safe <- read.csv("bladder_cancer_WikiPathways_safe_sample_perm.csv", header = TRUE)
res_padog <- read.csv("bladder_cancer_WikiPathways_gsea_sample_perm.csv", header = TRUE)

res_gsea_01 <- res_gsea %>%
  mutate(gsea_sample_perm=(PVAL < 0.05)+0) %>%
  select(c(set_id, name, gsea_sample_perm))

res_gsea_gene_01 <- res_gsea_gene_perm %>%
  mutate(set_id=ID, name=Description) %>%
  mutate(gsea_gene_perm=(qvalues < 0.05)+0) %>%
  select(c(set_id, name, gsea_gene_perm))

res_ora_01 <- res_ora %>%
  mutate(set_id=ID, name=Description) %>%
  mutate(ora=(qvalue < 0.05)+0) %>%
  select(c(set_id, name, ora))

res_safe_01 <- res_safe %>%
  mutate(safe_sample_perm=(PVAL < 0.05)+0) %>%
  select(c(set_id, name, safe_sample_perm))

res_padog_01 <- res_padog %>%
  mutate(padog_sample_perm=(PVAL < 0.05)+0) %>%
  select(c(set_id, name, padog_sample_perm))

res_rsea_01 <- res_rSEA %>%
  mutate(rsea=(Comp.adjP < 0.05)+0) %>%
  select(c(set_id, name, rsea))

res_upset <- merge(res_gsea_01, res_gsea_gene_01)
res_upset %<>% merge(.,res_ora_01)
res_upset %<>% merge(.,res_safe_01)
res_upset %<>% merge(.,res_padog_01)
res_upset %<>% merge(.,res_rsea_01)
require(UpSetR)
upset(res_upset, nsets = 6, text.scale = 2) 

##Gastric cancer network 1 genes
gastric_genes <-row.names(diff.res) %in%  as.character(wp_list[names(wp_list)=="WP2361"]$WP2361) 
neural_genes <- row.names(diff.res) %in%  as.character(wp_list[names(wp_list)=="WP4565"]$WP4565)

diff.res.illustrate <- diff.res %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  mutate(gastric_genes=as.factor(gastric_genes), neural_genes=as.factor(neural_genes))

(ggplot(diff.res.illustrate, aes(x=log2FoldChange, col=gastric_genes)) + 
  geom_density() +
  geom_vline(xintercept = 0, lty=2) +
  theme(text = element_text(size=20))) %>%
  ggsave(filename = "log2_fold_change_gastric_genes.pdf",
         plot = .,
         width = 7, 
         height = 7)
(ggplot(diff.res.illustrate, aes(x=log2FoldChange, col=neural_genes)) + 
    geom_density() +
    geom_vline(xintercept = 0, lty=2) +
    theme(text = element_text(size=20))) %>%
  ggsave(filename = "log2_fold_change_neural_genes.pdf",
         plot = .,
         width = 7, 
         height = 7)

wp_gene_freq <- table(wp_annotation$gene) %>% as.integer()

wp_gene_freq <- data.frame(wp_gene_freq)

(ggplot(wp_gene_freq, aes(x=Freq)) +
    geom_histogram() +
    theme(text = element_text(size=20))) %>%
  ggsave(filename = "histogram_no_of_WikiPathways_Per_Gene.pdf",
         plot = .,
         width = 7, 
         height = 7)
