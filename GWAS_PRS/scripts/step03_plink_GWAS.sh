########################################################################
#GWAS 
########################################################################
plink_dir=~/Desktop/handson_plink
cd ~/Desktop/handson_plink/GWAS
mkdir plink_results/
#linear model --glm or --linear
#no covariates
$plink_dir/plink2  --pfile genotypes/unrelated_hg38_auto_maf_0.05_500k_snp_filtered  \
--glm omit-ref  allow-no-covars  \
--pheno  phenotypes/phenofile.csv  --pheno-name BMI  \
--out plink_results/no.covar

#sex and age as covariate
$plink_dir/plink2  --pfile genotypes/unrelated_hg38_auto_maf_0.05_500k_snp_filtered  \
--covar    covariates/covariates.csv    --covar-variance-standardize   \
--glm omit-ref  hide-covar  \
--pheno  phenotypes/phenofile.csv  --pheno-name BMI  \
--out plink_results/sex.age.adjusted


#PCs, sex and age as covariate
$plink_dir/plink2 --pfile genotypes/unrelated_hg38_auto_maf_0.05_500k_snp_filtered  \
--covar    covariates/covariates_with_10PCs.csv    --covar-variance-standardize   \
--glm omit-ref  hide-covar  \
--pheno  phenotypes/phenofile.csv  --pheno-name BMI  \
--out plink_results/sex.age.PCs.adjusted




#logistic model --glm or --logistic
$plink_dir/plink2 --1 --pfile genotypes/unrelated_hg38_auto_maf_0.05_500k_snp_filtered  \
--covar    covariates/covariates_with_10PCs.csv    --covar-variance-standardize   \
--glm omit-ref cols=+orbeta,+ci hide-covar  \
--pheno  phenotypes/phenofile.csv  --pheno-name obesity --ci 0.95  \
--out plink_results/sex.age.PCs.adj

########################################################################
#open and run in Rstudio step04_GWAS_visualization.R