#Read table.
library(gplots)

dat <- read.table("norm_counts.txt", sep = "\t", header = T)
dat <- as.matrix(dat)

#Save in pdf format.
pdf("Heatmap.pdf", width = 10, height = 10)
heatmap.2(dat, 
          dendrogram = "column") #Use labRow = "" to turn off row labels.
dev.off()

#Save in TIFF format.
tiff("Heatmap.tiff", width = 10, height = 10, units = 'in', res = 300)
heatmap.2(dat, 
          dendrogram = "column") #Use labRow = "" to turn off row labels.
dev.off()
