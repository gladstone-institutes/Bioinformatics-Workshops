#files
cd ~/Desktop/handson_plink/GWAS/

plink_dir=~/Desktop/handson_plink
#cd genotypes/
#unrelated_hg38_auto_maf_0.01_500k_snp.pvar
#unrelated_hg38_auto_maf_0.01_500k_snp.psam
#unrelated_hg38_auto_maf_0.01_500k_snp.pgen
#cd phenotypes/
#phenofile.csv
#cd covariates/
#covariate.csv



########################################################################
#QC and filtering
########################################################################
mkdir QC/

#frequency
$plink_dir/plink2 --pfile genotypes/unrelated_hg38_auto_maf_0.01_500k_snp \
--freq \
--nonfounders \
--out QC/unrelated_hg38_auto_maf_0.01_500k_snp


#MAF filtering  (--maf 0.05)
$plink_dir/plink2 --pfile genotypes/unrelated_hg38_auto_maf_0.01_500k_snp --min-af 0.05 \
--make-pgen \
--out genotypes/unrelated_hg38_auto_maf_0.05_500k_snp


#call rate and missingness
$plink_dir/plink2 --pfile genotypes/unrelated_hg38_auto_maf_0.05_500k_snp --geno 0.1 \
--mind 0.1 \
--make-pgen \
--out genotypes/unrelated_hg38_auto_maf_0.05_500k_snp_filtered


#LD pruning - cretae list of independent snps to include in PCs analysis 
$plink_dir/plink2 --pfile genotypes/unrelated_hg38_auto_maf_0.05_500k_snp_filtered --indep-pairwise 50 5 0.2  \
--out QC/unrelated_hg38_auto_maf_0.05_500k_snp_filtered

#recode genotypes  with indep SNPs only 
$plink_dir/plink2 --pfile genotypes/unrelated_hg38_auto_maf_0.05_500k_snp_filtered \
--extract QC/unrelated_hg38_auto_maf_0.05_500k_snp_filtered.prune.in \
--make-pgen \
--out genotypes/unrelated_hg38_auto_maf_0.05_500k_filtered_indep_snp

#pca with indep snps only
$plink_dir/plink2 --pfile genotypes/unrelated_hg38_auto_maf_0.05_500k_filtered_indep_snp  --nonfounders    --pca           \
  --read-freq  QC/unrelated_hg38_auto_maf_0.01_500k_snp.afreq \
  --out QC/unrelated_hg38_auto_maf_0.05_500k_filtered_indep_snp

########################################################################

#run step02_covariate_PCs.R

