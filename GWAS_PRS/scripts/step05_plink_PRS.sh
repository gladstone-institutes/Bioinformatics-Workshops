
########################################################################
#PRS
########################################################################
PRSdir=~/Desktop/handson_plink/PRS
GWASdir=~/Desktop/handson_plink/GWAS

cd ~/Desktop/handson_plink/PRS
mkdir PRS_results/
plink_dir=~/Desktop/handson_plink
#input files are in EUR_input_for_PRS/
#EUR summary stats from large GWAS meta analysis 
#Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt

#EUR list of indivs + covariate and pheno files
#EUR_input_for_PRS/EUR_BMI_phenofile.csv
#EUR_input_for_PRS/EUR_covariates.csv
#EUR_input_for_PRS/EUR_indivs.csv


#create geno files for european only as input for the PRS - we used maf 0.01 since we did not filter for call rate and missingness
$plink_dir/plink2 \
--pfile $GWASdir/genotypes/unrelated_hg38_auto_maf_0.01_500k_snp  \
--keep  EUR_input_for_PRS/EUR_indivs.csv  \
--make-bed \
--out EUR_input_for_PRS/EUR_unrelated_hg38_auto_maf_0.01_500k_snp_filtered

#clumping the summary stats from large GWAS meta analysis to get indep SNPs based on pvalue overlapping our genotypes
$plink_dir/plink2  \
    --bfile EUR_input_for_PRS/EUR_unrelated_hg38_auto_maf_0.01_500k_snp_filtered \
    --clump-p1 1 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --clump EUR_input_for_PRS/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt \
    --clump-snp-field SNP \
    --clump-field P \
    --out EUR_input_for_PRS/EUR_clumped
    
#valid independent snps    
awk 'NR!=1{print $3}' EUR_input_for_PRS/EUR_clumped.clumps >  EUR_input_for_PRS/EUR.valid.snp
    
    
#SNP and P from large GWAS meta analysis
awk '{print $3,$9}' EUR_input_for_PRS/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt > EUR_input_for_PRS/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.SNP.pvalue

#create threshold of SNPs to create scores for different subsets of SNPs
echo "0.001 0 0.001" > EUR_input_for_PRS/range_list 
echo "0.05 0 0.05" >> EUR_input_for_PRS/range_list
echo "0.1 0 0.1" >> EUR_input_for_PRS/range_list
echo "0.2 0 0.2" >> EUR_input_for_PRS/range_list
echo "0.3 0 0.3" >> EUR_input_for_PRS/range_list
echo "0.4 0 0.4" >> EUR_input_for_PRS/range_list
echo "0.5 0 0.5" >> EUR_input_for_PRS/range_list

#estimating the score for each range of P using the valid snps (for the score needs snp, allele and BETA)
$plink_dir/plink2  \
    --bfile EUR_input_for_PRS/EUR_unrelated_hg38_auto_maf_0.01_500k_snp_filtered \
    --score EUR_input_for_PRS/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt 3 4 7 header \
    --q-score-range EUR_input_for_PRS/range_list EUR_input_for_PRS/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.SNP.pvalue \
    --extract EUR_input_for_PRS/EUR.valid.snp \
    --out EUR_input_for_PRS/EUR
########################################################################
