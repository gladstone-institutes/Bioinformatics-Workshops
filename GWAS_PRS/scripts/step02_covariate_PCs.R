#in Rstudio

########################################################################
#plot PCs and creating a second covariate file including PCS
########################################################################

setwd("~/Desktop/handson_plink/GWAS/")
cov <- read.csv("covariates/covariates.csv", header=T,comment.char="", check.names=F)
pc <- read.table("QC/unrelated_hg38_auto_maf_0.05_500k_filtered_indep_snp.eigenvec",header=T,  comment.char="", check.names=F)
colnames(pc)[1] <- "IID"
EUR <- read.csv("../PRS/EUR_input_for_PRS/EUR_indivs.csv", header=T,comment.char="", check.names=F)

updated_pc <- merge(cov, pc, by="IID")
updated_pc_EUR <- updated_pc[which(updated_pc$IID %in% EUR$IID),]

if (!dir.exists("R_plots")){
dir.create("R_plots")
} else {
    print("R_plots already exists!")
}


pdf("R_plots/Pop_stratification.pdf")
plot(updated_pc$PC1, updated_pc$PC2, xlab="PC1", ylab="PC2")
points(updated_pc_EUR$PC1, updated_pc_EUR$PC2, pch=19, col="forestgreen")
dev.off()

write.csv(updated_pc, "covariates/covariates_with_10PCs.csv",  quote=F, row.names=F)

########################################################################

#open and run step03_GWAS.sh

########################################################################
