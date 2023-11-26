####################################################
#Best PRS model and Visualization in R
####################################################

setwd("~/Desktop/handson_plink/PRS/")

p.threshold <- c(0.001,0.05,0.1,0.2,0.3,0.4,0.5)
# Read in the phenotype file 
phenotype <- read.csv("EUR_input_for_PRS/EUR_BMI_phenofile.csv", header=T, check.names=F, comment.char="")
# Read in the PCs
#pcs <- read.table("EUR.eigenvec", header=F)
# The default output from plink does not include a header
# To make things simple, we will add the appropriate headers
# (1:6 because there are 6 PCs)
#colnames(pcs) <- c("FID", "IID", paste0("PC",1:6)) 
# Read in the covariates (here, it is sex and age)
covariate <- read.csv("EUR_input_for_PRS/EUR_covariates.csv", header=T,check.names=F, comment.char="")
# Now merge the files
pheno <- merge(phenotype, covariate, by=c("IID"))

#pheno <- merge(merge(phenotype, covariate, by=c("FID", "IID")), pcs, by=c("FID","IID"))
# We can then calculate the null model (model with PRS) using a linear regression 
# (as bmi is quantitative)
null.model <- lm(BMI~., data=pheno[,-1])
# And the R2 of the null model is 
null.r2 <- summary(null.model)$r.squared
prs.result <- NULL
for(i in p.threshold){
    # Go through each p-value threshold
    prs <- read.table(paste0("EUR_input_for_PRS/EUR.",i,".sscore"), header=T, check.names=F, comment.char="")
    # Merge the prs with the phenotype matrix
    # We only want the FID, IID and PRS from the PRS file, therefore we only select the 
    # relevant columns
    pheno.prs <- merge(pheno, prs[,c("IID", "SCORE1_AVG")], by= "IID")
    # Now perform a linear regression on BMI with PRS and the covariates
    # ignoring the FID and IID from our model
    model <- lm(BMI~., data=pheno.prs[,-1])
    # model R2 is obtained as 
    model.r2 <- summary(model)$r.squared
    # R2 of PRS is simply calculated as the model R2 minus the null R2
    prs.r2 <- model.r2-null.r2
    # We can also obtain the coeffcient and p-value of association of PRS as follow
    prs.coef <- summary(model)$coeff["SCORE1_AVG",]
    prs.beta <- as.numeric(prs.coef[1])
    prs.se <- as.numeric(prs.coef[2])
    prs.p <- as.numeric(prs.coef[4])
    # We can then store the results
    prs.result <- rbind(prs.result, data.frame(Threshold=i, R2=prs.r2, P=prs.p, BETA=prs.beta,SE=prs.se))
}
# Best result is:
prs.result[which.max(prs.result$R2),]
write.csv(prs.result, file="PRS_results/EUR_BMI_sex_age_adjusted_PRS_fit.csv", quote=F, row.names=F)

best_p <- prs.result[which.max(prs.result$R2),1]
best_prs <- read.table(paste0("EUR_input_for_PRS/EUR.",best_p,".sscore"), header=T,check.names=F, comment.char="")
# Rename the sex
covariate$sex <- as.factor(covariate$sex)
levels(covariate$sex) <- c("Male", "Female")
# Merge the files
dat <- merge(merge(best_prs, phenotype, by="IID"), covariate, by="IID")

# Start plotting
pdf("PRS_results/EUR_BMI_sex_age_adjusted_PRS_fit.pdf")
plot(x=dat$SCORE1_AVG, y=dat$BMI, col="white",
    xlab="Polygenic Score", ylab="BMI")
with(subset(dat, sex=="Male"), points(x=SCORE1_AVG, y=BMI, col="red", pch=19))
with(subset(dat, sex=="Female"), points(x=SCORE1_AVG, y=BMI, col="blue", pch=19))

plot(density(dat$SCORE1_AVG[dat$sex =="Male"]), col="red", lwd=3)
lines(density(dat$SCORE1_AVG[dat$sex =="Female"]), col="blue", lwd=3) 
dev.off()
####################################################
