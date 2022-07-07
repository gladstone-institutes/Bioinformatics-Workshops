#Reading data file.
dat <- read.table("GSE60450_Lactation-GenewiseCounts.txt")

#Let us examine how our data looks.
View(dat)

#Seems like all data points are there. Can we improve appearance?
#Examine the details of read.table command.
?read.table

#Looks like we can inform read.table about separator type and presence of header.
dat <- read.table("GSE60450_Lactation-GenewiseCounts.txt",
                  header= TRUE, sep = "\t")

#Let us examine how data looks now.
View(dat)

#What are the observations represented in our data?
#colnames gives the names of columns of data.
colnames(dat)

#To check the number of rows and columns in table.
dim(dat)

#To check first few rows of table.
head(dat)

#To check last few rows of table.
tail(dat)

#To check basic stats for each column.
summary(dat)

#Let us extract a column of data.
#For example, EntrezGeneID.
geneIds <- dat$EntrezGeneID

#Check what kind of variable geneIds is.
class(geneIds)

#geneIds should be string type.
dat$EntrezGeneID <- as.character(dat$EntrezGeneID)

#Check the class of gene ids again
class(dat$EntrezGeneID)

#Information about samples is in another file.
phenotype_info <- read.table("targets.csv",
                             header = TRUE,
                             sep = ",")

#The column named GEO in the table represents sample id on GEO.
#Status and CellType are factor levels for statistical analysis.
phenotype_info$CellType <- as.factor(phenotype_info$CellType)

#Check class of CellType. It is a factor variable now.
class(phenotype_info$CellType)

#Currently Status and CellType has repetition of the same values.
#Get unique values.
celltypes <- unique(phenotype_info$CellType)

#Perhaps, no point in keeping spcs as factor in the above object.
#Convert celltypes to character variable.
celltypes <- as.character(celltypes)

#Checking if a text is present in a character variable?
"cardiomyocyte" %in% celltypes

#Let us say we want subset of data corresponding to B cells.
which_B <- phenotype_info$CellType == "B"
phenotype_info_B <- phenotype_info[which_B, ]

#Class of which_B? Logical
class(which_B)

#Alternative way to subset data.
phenotype_info_B <- subset(phenotype_info,
                           CellType == "B")

#Check mean counts for a sample.
mean(dat$MCL1.DG_BC2CTUACXX_ACTTGA_L002_R1)

#Check mean counts for all genes in the B cell samples.
clnames_dat <- colnames(dat)[-2:-1]
clnames_dat <- strsplit(clnames_dat, split = "_")
clnames_dat <- data.frame(clnames_dat)
clnames_dat <- t(clnames_dat)
rownames(clnames_dat) <- NULL
clnames_dat <- clnames_dat[, 1]
which_B <- which(clnames_dat %in% phenotype_info_B$X)
cnts_B <- dat[, which_B + 2]

#Estimate median counts for casein protein.
median(unlist(dat[dat$EntrezGeneID == "12992", -2:-1]))
median(unlist(dat[dat$EntrezGeneID == "12992", which_B + 2]))

#Estimate standard deviation.
sd(unlist(dat[dat$EntrezGeneID == "12992", -2:-1]))
sd(unlist(dat[dat$EntrezGeneID == "12992", which_B + 2]))

#Check histograms of data.
hist(log2(1 + dat$MCL1.DG_BC2CTUACXX_ACTTGA_L002_R1))
hist(log2(1 + dat$MCL1.LC_BC2CTUACXX_GCCAAT_L001_R1))

#These histograms are not easy to compare. 
#Perhaps, we can fix the axis limits.
hist(log2(1 + dat$MCL1.DG_BC2CTUACXX_ACTTGA_L002_R1),
     xlim = c(0, 20), ylim = c(0, 15000))
hist(log2(1 + dat$MCL1.LC_BC2CTUACXX_GCCAAT_L001_R1),
     xlim = c(0, 20), ylim = c(0, 15000))

#Let us look at boxplots.
boxplot(log2(1 + dat$MCL1.DG_BC2CTUACXX_ACTTGA_L002_R1),
        ylim = c(0, 25))
boxplot(log2(1 + dat$MCL1.LC_BC2CTUACXX_GCCAAT_L001_R1),
        ylim = c(0, 25))

#Scatter plots.
plot(x= log2(1 + dat$MCL1.DG_BC2CTUACXX_ACTTGA_L002_R1),
     y= log2(1 + dat$MCL1.LC_BC2CTUACXX_GCCAAT_L001_R1))
#Are higher counts associated with longer genes?
#Can  we color data points based on length?

#install.packages("ggplot2")
library(ggplot2)
qplot(x = log2(1 + MCL1.DG_BC2CTUACXX_ACTTGA_L002_R1), 
      y = log2(1 + MCL1.LC_BC2CTUACXX_GCCAAT_L001_R1), 
      data = dat[1:1000, ], color = log10(Length))

#Color scale does not give clear insight.
#Can we use discrete color scale of points instead?
dat_subset <- dat[1:1000, ]
dat_subset$genetype <- "smallGene"
dat_subset$genetype[dat_subset$Length > median(dat_subset$Length)] <- "longGene"
qplot(x = log2(1 + MCL1.DG_BC2CTUACXX_ACTTGA_L002_R1), 
      y = log2(1 + MCL1.LC_BC2CTUACXX_GCCAAT_L001_R1), 
      data = dat_subset, color = genetype)

#Check the boxplots for long and small genes simultaneously.
qplot(x = genetype, 
      y = log2(1 + MCL1.LC_BC2CTUACXX_GCCAAT_L001_R1), 
      data = dat_subset, geom = "boxplot")

#Additional stuff.
#Hypothesis test.
#Are counts for long genes significantly different from the small genes?
dat_small <- subset(dat_subset, genetype == "smallGene")
dat_long <- subset(dat_subset, genetype == "longGene")
test <- t.test(dat_small$MCL1.LC_BC2CTUACXX_GCCAAT_L001_R1, 
               dat_long$MCL1.LC_BC2CTUACXX_GCCAAT_L001_R1)
test$p.value

#Conditional statements take a condition and perform steps depending on validitiy of the statement.
if (test$p.value < 0.05) {
  print("Counts depend on gene length.")
} else {
  print("Counts don't depend on gene length.")
}

#Looping for repeating the same task but on different data.
for (item in c("smallGene", "longGene")) {
  y <- subset(dat_subset, genetype == item)
  smry <- summary(y)
  print(item)
  print(smry)
}





