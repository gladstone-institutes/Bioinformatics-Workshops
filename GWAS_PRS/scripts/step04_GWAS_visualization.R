
##################
#Visualization in R
##################
#manhattan plot
install.packages("qqman")
library(qqman)

setwd("~/Desktop/handson_plink/GWAS/")

if (!dir.exists("R_plots")){
dir.create("R_plots")
} else {
    print("R_plots already exists!")
}

gwas_obesity <-  read.table("plink_results/sex.age.PCs.adj.obesity.glm.logistic.hybrid", header=T, comment.char="", check.names=F)

pdf("R_plots/obesity.pdf")
manhattan(gwas_obesity, chr="#CHROM", bp="POS", snp="ID", p="P" , annotatePval = 0.0001, main="Manhattan plot - sex, age, PCs adjusted obesity")

#qqplot
qq(gwas_obesity$P, main="quantile-quantile plot - sex, age, PCs adjusted obesity")
dev.off()


gwas_BMI_no_cov <-  read.table("plink_results/no.covar.BMI.glm.linear", header=T,comment.char="", check.names=F)

gwas_BMI_sex_age_adj <-  read.table("plink_results/sex.age.adjusted.BMI.glm.linear", header=T,comment.char="", check.names=F)

gwas_BMI_sex_age_PCs_adj <-  read.table("plink_results/sex.age.PCs.adjusted.BMI.glm.linear", header=T,comment.char="", check.names=F)


pdf("R_plots/BMI.pdf")
manhattan(gwas_BMI_no_cov, chr="#CHROM", bp="POS", snp="ID", p="P" , annotatePval = 0.0001, main="Manhattan plot - no covar adjusted BMI")

manhattan(gwas_BMI_sex_age_adj, chr="#CHROM", bp="POS", snp="ID", p="P" , annotatePval = 0.0001, main="Manhattan plot - sex, age adjusted BMI")

manhattan(gwas_BMI_sex_age_PCs_adj, chr="#CHROM", bp="POS", snp="ID", p="P" , annotatePval = 0.0001, main="Manhattan plot - sex, age, PCs adjusted BMI")

#qqplot
qq(gwas_BMI_no_cov$P, main="quantile-quantile plot - no covar adjusted BMI")
qq(gwas_BMI_sex_age_adj$P, main="quantile-quantile plot - sex, age adjusted BMI")
qq(gwas_BMI_sex_age_PCs_adj$P, main="quantile-quantile plot - sex, age, PCs adjusted BMI")

dev.off()
#include
#pvals <- read.table("AMH_all_FVG_rs_pvalue", header=T)

#observed <- sort(pvals$pgc)
#lobs <- -(log10(observed))

#expected <- c(1:length(observed)) 
#lexp <- -(log10(expected / (length(expected)+1)))


########################################################################
#open and run step05_Plink_PRS.sh
