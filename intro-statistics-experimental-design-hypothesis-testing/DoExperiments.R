rm(list = ls())
setwd("~/Dropbox (Gladstone)/scripts/Bioinformatics-Workshops/intro-experimental-design")

##load assay functions
source("SImulateData.R")

##do experiment comparing expression of gene between E9.5 and E11.5
PerformAssayGeneExpression(NsamplesPerTime = 100)

##Two models for expression of gene at E9_5: Gamma and Normal
AssayExpression_E9_5_two_models(NsamplesPerTime = 100)

##Repeat the experiment Nsim times to estimate mean gene expression at E9_5 under the two models of 
##data generating distribution
Repeat_AssayExpression_E9_5_two_models(NsamplesPerTime = 100, Nrepeat = 10000)

