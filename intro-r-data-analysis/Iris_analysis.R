#Reading data file.
dat <- read.table("iris.csv")

#Let us examine how our data looks.
View(dat)

#Seems like all data points are there. Can we improve appearance?
#Examine the details of read.table command.
?read.table

#Looks like we can inform read.table about separator type and presence of header.
dat <- read.table("iris.csv", header= TRUE, sep = ",")

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
#For example, sepal length.
spl_len <- dat$Sepal.Length

#Check what kind of variable spl_len is.
class(spl_len)

#Which species of Iris are represented in the data? 
spcs <- dat$Species

#Check class of spcs. It is a factor variable.
class(spcs)

#Current value of spcs has repetition of each spcs type.
#Get unique values.
spcs <- unique(spcs)

#Perhaps, no point in keeping spcs as factor now.
#Convert spcs to character variable.
spcs <- as.character(spcs)

#Checking if a text is present in a character variable?
"sapiens" %in% spcs


#Let us say we want subset of data corresponding to Iris setosa.
which_setosa <- dat$Species == "setosa"
dat_setosa <- dat[which_setosa, ]

#Class of which_setosa? Logical
class(which_setosa)

#Alternative way to subset data.
dat_setosa <- subset(dat, Species == "setosa")

#Check mean Sepal length for all observations.
mean(dat$Sepal.Length)

#Check mean Sepal length for Iris setosa only.
mean(dat_setosa$Sepal.Length)

#Estimate median.
median(dat$Petal.Width)
median(dat_setosa$Sepal.Length)

#Estimate standard deviation.
sd(dat$Petal.Length)
sd(dat_setosa$Sepal.Width)

#Check histograms of data.
hist(dat$Sepal.Length)
hist(dat_setosa$Sepal.Length)

#These histograms are not easy to compare. 
#Perhaps, we can fix the axis limits.
hist(dat$Sepal.Length, xlim = c(4, 8), ylim = c(0, 30))
hist(dat_setosa$Sepal.Length, xlim = c(4, 8), ylim = c(0, 30))

#Let us look at boxplots.
boxplot(dat$Sepal.Length, ylim = c(3, 9))
boxplot(dat_setosa$Sepal.Length, ylim = c(3, 9))

#Scatter plots.
plot(x= dat$Sepal.Length, y = dat$Petal.Length)
#There is a cluster of data points in the lower left corner. 
#Are these data points from one particular species?
#Can  we color data points based on species?

#install.packages("ggplot2")
library(ggplot2)
qplot(x = Sepal.Length, y = Petal.Length, data = dat, color = Species)

#But the journals charge extra for color figures.
#Can we use shapes to distinguish species?
qplot(x = Sepal.Length, y = Petal.Length, data = dat, shape = Species)

#Check the boxplots for all species simultaneously.
qplot(x = Species, y = Sepal.Length, data = dat, geom = "boxplot")

#Additional stuff.
#Hypothesis test.
#Is sepal length for Iris setosa significantly different from the other two species?
dat_setosa <- subset(dat, Species == "setosa")
dat_other <- subset(dat, Species != "setosa")
test <- t.test(dat_setosa$Sepal.Length, dat_other$Sepal.Length)
test$p.value

#Conditional statements take a condition and perform steps depending on validitiy of the statement.
if (test$p.value < 0.05) {
  print("Iris setosa can be classified based on sepal length.")
} else {
  print("Iris setosa may not be identified based on sepal length.")
}

#Looping for repeating the same task but on different data.
for (item in spcs) {
  y <- subset(dat, Species == item)
  smry <- summary(y)
  print(item)
  print(smry)
}
