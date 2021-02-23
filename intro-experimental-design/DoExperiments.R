rm(list = ls())
setwd("~/Dropbox (Gladstone)/Bioinformatics/Training_Workshops/Gladstone-internal/Introduction to statistics and experimental design - Feb 2021")

source("SImulateData.R")

##do experiment comparing expression of gene between E9.5 and E11.5
PerformAssayGeneExpression(NsamplesPerTime = 10)

##Two models for expression of gene at E9_5: Gamma and Normal
AssayExpression_E9_5_two_models(NsamplesPerTime = 4)

##Repeat the experiment Nsim times to estimate mean gene expression at E9_5 under the two models of 
##data generating distribution
Repeat_AssayExpression_E9_5_two_models(NsamplesPerTime = 10, Nrepeat = 10000)

