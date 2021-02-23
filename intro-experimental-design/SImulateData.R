require(ggplot2)
require(gridExtra)

set.seed(1234)
##do experiment
PerformAssayGeneExpression <- function(NsamplesPerTime) {
  N=NsamplesPerTime
  Shape = 10
  Rate = 0.5
  delta = 0.2
  Y1 <- rgamma(N, shape = Shape, rate = Rate)
  Y2 <- rgamma(N, shape=(1+delta)*Shape, rate = Rate)
  Y <- append(Y1, Y2)
  Timepoint <- c(rep("E9_5", N), rep("E11_5", N))
  PlotData <- data.frame(Y, Timepoint)
  PlotData$Timepoint <- relevel(PlotData$Timepoint, ref="E9_5")
  print(paste0("Mean expression of Gene at E9_5 = ", round(mean(Y1), digits = 1)))
  print(paste0("Mean expression of Gene at E11_5 = ", round(mean(Y2), digits = 1)))
  
  print("=====")
  print("=====")

  if(NsamplesPerTime > 2) {
      print(paste0("SD expression of Gene at E9_5 = ", round(sd(Y1), digits = 1)))
      print(paste0("SD/sqrt(NsamplesPerTime) expression of Gene at E9_5 = ", round(sd(Y1)/sqrt(NsamplesPerTime), digits = 1)))
   }
  print(ggplot(PlotData, aes(x=Timepoint, y=Y)) + 
          geom_boxplot() + 
          geom_jitter(width=0.25) + 
          xlab("embryonic time") + 
          ylab("gene expression") +
          ylim(10,35))
}



##Two models for expression of gene at E9_5: Gamma and Normal
AssayExpression_E9_5_two_models <-  function(NsamplesPerTime = 100) {
  N = NsamplesPerTime
  Shape = 10
  Rate = 0.5
  ##Gamma distribution
  Yg <- rgamma(N, shape = Shape, rate=Rate)
  ##Normal distribution
  Yn <- rnorm(N, mean=Shape/Rate, sd=sqrt(Shape)/Rate)
  
  ##Organize the data
  Y <- append(Yg, Yn)
  Distribution <- c(rep("Gamma", N), rep("Normal", N))
  data <- data.frame(Y, Distribution, Condition=rep("E9_5", length(Y)))
  
  ##Visualize the data
  p1 <- (ggplot(data, aes(x=Condition, y=Y, color=Distribution))  + geom_boxplot() + geom_jitter(width=0.25))
  p2 <- (ggplot(data, aes(x=Condition, y=Y, color=Distribution))  + geom_violin() +  geom_jitter(width=0.25)) 
  p3 <- (ggplot(data, aes(x=Y, color=Distribution))  + geom_histogram()) 
  print(grid.arrange(p1,p2,p3,ncol=2))
  # ggplot(data, aes(x=Y, color=Distribution))  + geom_freqpoly() 
  
  print(paste0("Mean of gene expression under Gamma model = ",round(mean(Y[Distribution=="Gamma"]), digits = 1)))
  print(paste0("Mean of gene expression under Normal model = ", round(mean(Y[Distribution=="Normal"]), digits = 1)))

  print("=====")
  print("=====")
  
  print(paste0("SD of gene expression under Gamma model = ",round(sd(Y[Distribution=="Gamma"]), digits = 1)))
  print(paste0("SD of gene expression under Normal model = ", round(sd(Y[Distribution=="Normal"]), digits = 1)))

}


##Repeat the experiment Nsim times to estimate mean gene expression at E9_5 under the two models of 
##data generating distribution
Repeat_AssayExpression_E9_5_two_models <- function(NsamplesPerTime=4, Nrepeat=1000) {
  Shape = 10
  Rate = 0.5
  Nrepeat = Nrepeat
  Ygmean <- vector(mode = "numeric")
  Ygsd <- vector(mode = "numeric")
  Zg <- vector(mode = "numeric")
  
  Ynmean <- vector(mode = "numeric")
  Ynsd <- vector(mode = "numeric")
  Zn <- vector(mode = "numeric")
  
  N = NsamplesPerTime
  
  for(i in 1:Nrepeat) {
    Yg <- rgamma(N, shape = Shape, rate=Rate)
    Ygmean[i] <- mean(Yg)
    Ygsd[i] <- sd(Yg)
    Zg[i] <- sqrt(N)*(Ygmean[i] - (Shape/Rate))/(sqrt(Shape)/Rate)
    
    Yn <- rnorm(N, mean=Shape/Rate, sd=sqrt(Shape)/Rate)
    Ynmean[i] <- mean(Yn)
    Ynsd[i] <- sd(Yn)
    Zn[i] <- sqrt(N)*(Ynmean[i] - (Shape/Rate))/(sqrt(Shape)/Rate)
  
  }
  
  Ymean <- append(Ygmean, Ynmean)
  Ysd <- append(Ygsd, Ynsd)
  Z <- append(Zg, Zn)
  Distribution <- c(rep("Gamma", Nrepeat), rep("Normal", Nrepeat))
  Condition <- rep("E9_5", length(Ymean))
  CLT_data <- data.frame(Ymean, Ysd, Z, Distribution, Condition)
  
  p1 <-  (ggplot(CLT_data, aes(x=Condition, y=Ymean, color=Distribution))  + geom_violin() )
  p2 <- (ggplot(CLT_data, aes(x=Ymean, color=Distribution))  + geom_histogram() )
  print(grid.arrange(p1, ncol=1))
  
  print(paste0("SD or Precision of mean gene expression under Gamma model = ",  round(sd(Ymean[Distribution=="Gamma"]), digits = 1)))
  print(paste0("SD or Precision of gene expression under Normal model = ", round(sd(Ymean[Distribution=="Normal"]), digits = 1)))

}

