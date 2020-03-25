#setwd()
setwd("~/Dropbox (Gladstone)/Bioinformatics/Training_Workshops/Gladstone-internal/Intermediate_RNA-seq_Fall_2019")

library(magrittr)
library(edgeR)
library(org.Mm.eg.db)
library(ggplot2)
library(tidyverse)
library(vioplot)

phenotype_info_file <- "targets.txt"
raw_counts_file <- "GSE60450_Lactation-GenewiseCounts.txt.gz"

#12 samples and 6 categories
targets <- phenotype_info_file %>%
  read.delim(., stringsAsFactors=FALSE)

#This is equivalent to 
# targets <- read.delim(phenotype_info_file, stringsAsFactors = FALSE)

group <- targets %$%
  paste(CellType, Status, sep = ".") %>%
  factor()

#Length of a gene is the total number of bases in exons and UTRs for that gene.
GenewiseCounts <- raw_counts_file %>%
  read.delim(., row.names="EntrezGeneID")

colnames(GenewiseCounts) %<>% substring(.,1,7)

# colnames(GenewiseCounts)[-1] <- group

#------------------------
#Concept 1: MA plots
#------------------------

two_samples <- GenewiseCounts[, c(2, 3)] %>% #Replicate samples
  add(., 1) %>%
  log2()

plotData <- data.frame(M = two_samples[, 1] - two_samples[, 2],
                       A = (two_samples[, 1] + two_samples[, 2])/2)

ggplot(plotData, aes(x = A, y = M)) +
  geom_point() +
  geom_smooth()

#------------------------
#Create DGElist object and retrieve gene symbols.
#------------------------
y <- DGEList(counts = GenewiseCounts[,-1], 
             group=group,
             genes=GenewiseCounts[,1,drop=FALSE])

y$genes$Symbol <- mapIds(org.Mm.eg.db,
                         keys = rownames(y),
                         keytype="ENTREZID", 
                         column="SYMBOL")

#------------------------
#Independent filtering.
#------------------------

#Filter genes whose symbols are not found.
keep <- y$genes$Symbol %>%
  is.na() %>%
  not()

y <- y[keep, ]


#It is not possible to make reliable inference for genes ...
#... with very low counts.
#Can filter directly by counts but ...
#... better to account for difference in library sizes.

minimum_counts_reqd <- 10
cutoff <- y$samples$lib.size %>% 
  median() %>% 
  divide_by(., 10^6) %>% 
  divide_by(minimum_counts_reqd, .) %>% 
  round(., 1)

keep <- cpm(y) %>%
  is_greater_than(., cutoff) %>%
  rowSums() %>%
  is_weakly_greater_than(., 2)

#Alternatively, use counts to filter.
#Examples:
#keep <- rowSums(y$counts) > 50

y <- y[keep, , keep.lib.sizes=FALSE]

#-------------------------
#Normalizing the counts.
#-------------------------
y <- calcNormFactors(y)

#What's under the hood?

cnts <- y$counts %>% as.matrix()
lib.sizes <- apply(cnts, 2, sum)

cnts_adjst_libsize <- map_dfc(colnames(cnts), 
                              function(x) cnts[, x]/lib.sizes[x]) %>%
  set_colnames(., colnames(cnts))

vioplot(cnts_adjst_libsize)
#Intuitively, we should expect similar adjustments for similar samples.

f <- apply(cnts_adjst_libsize, 2, 
           function(x) quantile(x,p=0.75))

ref <- (f - mean(f)) %>% 
  abs() %>% 
  which.min()

TMM_norm_factors <- 
  map_dfc(colnames(cnts), 
          function(x, logratioTrim=.3, sumTrim=0.05, 
                   doWeighting=TRUE, Acutoff=-1e10) {
            
            #The following steps are excerpted from edgeR's source.
            nO <- lib.sizes[x]
            nR <- lib.sizes[ref]
            
            obs <- cnts[, x] %>% as.numeric()
            ref <- cnts[, ref] %>% as.numeric()
            
            logR <- log2(obs/nO) - log2(ref/nR)         # log ratio of expression, accounting for library size
            absE <- (log2(obs/nO) + log2(ref/nR))/2  # absolute expression
            v <- (nO-obs)/nO/obs + (nR-ref)/nR/ref   # estimated asymptotic variance
            
            #	remove infinite values, cutoff based on A
            fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)
            
            logR <- logR[fin]
            absE <- absE[fin]
            v <- v[fin]
            
            if(max(abs(logR)) < 1e-6) return(1)
            
            #	taken from the original mean() function
            n <- length(logR)
            loL <- floor(n * logratioTrim) + 1
            hiL <- n + 1 - loL
            loS <- floor(n * sumTrim) + 1
            hiS <- n + 1 - loS
            
            keep <- (rank(logR)>=loL & rank(logR)<=hiL) & 
              (rank(absE)>=loS & rank(absE)<=hiS)
            
            if(doWeighting)
              f <- sum(logR[keep]/v[keep], na.rm=TRUE) / sum(1/v[keep], na.rm=TRUE)
            else
              f <- mean(logR[keep], na.rm=TRUE)
            
            #	Results will be missing if the two libraries share no features with positive counts
            #	In this case, return unity
            if(is.na(f)) f <- 0
            2^f
          }) %>%
  data.frame() %>% 
  as.numeric() 

#Rescale norm factors for convenience of interpretation.
rescale <- TMM_norm_factors %>% 
  log() %>%
  mean() %>%
  exp()

TMM_norm_factors %<>% divide_by(., rescale)

#Plot the normalization factors by sample.
ggplot(y$samples %>% 
         cbind(., replicate = factor(1:2)), 
       aes(x = group, y = norm.factors, fill = replicate)) +
  geom_bar(stat= "identity", position = position_dodge())

#"A normalization factor below one indicates that a small number of high count genes 
#...are monopolizing the sequencing, causing the counts for other genes to be lower 
#...than would be usual given the library size." Chen et al., 2016
#Looks like L.lactating samples contain a number of very highly upregulated genes.

#------------------------
#Exploratory visualizations
#------------------------
#MDS plot
pch <- c(0,1,2,15,16,17)
colors <- rep(c("darkgreen", "red", "blue"), 2)
plotMDS(y, col=colors[group], pch=pch[group])
legend("top", legend=levels(group), pch=pch, col=colors, ncol=2, 
       text.width = 0.1)

#PCA plot
cpm <- cpm(y, log = T, prior.count = 0.01)
rv <- apply(cpm,1,var) 

#Select genes with highest variance.
keep <- order(rv, decreasing = TRUE)[1:500]
selected <- cpm[keep, ] %>% t()

#Transpose is needed to ensure that each row is a vector.
pca <- prcomp(selected, scale=T, center = T)

stddev <- pca$sdev
pc1_var <- round(100*stddev[1]^2/sum(stddev^2))
pc2_var <- round(100*stddev[2]^2/sum(stddev^2))
pc3_var <- round(100*stddev[3]^2/sum(stddev^2))
PlotData <- data.frame(cbind(PC1 = pca$x[,1], PC2 = pca$x[,2]))
PlotData <- targets[, c("CellType", "Status")] %>%
  cbind(PlotData, .)

ggplot(PlotData, aes(x=PC1, y=PC2, color=CellType, shape=Status)) + 
  geom_point(size=4.5) +  
  xlab(paste("PC1:", pc1_var, "% variance")) + 
  ylab(paste("PC2:", pc2_var, "% variance"))

#--------------------
#Fitting the model.
#--------------------
design <- model.matrix(~0+group) %>%   
  set_colnames(., levels(group))

y <- estimateDisp(y, design, robust=TRUE)
plotBCV(y)

fit <- glmQLFit(y, design, robust=TRUE) 
head(fit$coefficients)

plotQLDisp(fit) 

#--------------------
#Hypothesis testing
#--------------------

#Example 1:
B.LvsP <- makeContrasts(B.lactating-B.pregnant, levels=design)
res <- glmQLFTest(fit, contrast=B.LvsP)
topTags(res)
is.de <- decideTestsDGE(res)
summary(is.de)
plotMD(res, status=is.de, values=c(1,-1), col=c("red","blue"),
       legend="topright")

#Example 2:
B.LvsP <- makeContrasts(B.lactating-B.pregnant, levels=design)
res <- glmTreat(fit, contrast=B.LvsP, lfc=log2(1.5))
topTags(res)
is.de <- decideTestsDGE(res)
summary(is.de)
plotMD(res, status=is.de, values=c(1,-1), col=c("red","blue"),
       legend="topright")

#Example 3:
con <- makeContrasts(
  (L.lactating-L.pregnant)-(B.lactating-B.pregnant),
  levels=design)
res <- glmQLFTest(fit, contrast=con)
topTags(res)
is.de <- decideTestsDGE(res)
summary(is.de)
plotMD(res, status=is.de, values=c(1,-1), col=c("red","blue"),
       legend="topright")

#Example 4:
con <- makeContrasts(
  L.PvsL = L.pregnant - L.lactating,
  L.VvsL = L.virgin - L.lactating,
  L.VvsP = L.virgin - L.pregnant, levels=design)
res <- glmQLFTest(fit, contrast=con)
topTags(res)
is.de <- decideTestsDGE(res)
summary(is.de)

#----------------------------
#Save results in a table
result <- topTags(res, n = nrow(y$counts)) %>%
  data.frame()
write.table(result, 
            file = "DE_result.txt",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)



