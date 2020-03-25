#making animations.

#Get a range of substrate concentration values. 
S <- 0:10000*0.01
      
#Half-maximal concentration constant.
K <- 50

#Define a container object for data.
#One column for substrate concentration and other 50 for Hill coefficients.
dat <- matrix(, 10001, 51) 

#Assign substrate concentration values to first column.
dat[, 1] <- S

#Loop through all Hill coefficients of interest.
for (i in 1:50) {
  
  #Get reaction velocity in terms of fraction of max. velocity.
  this_y <- (S^i)/((S^i)+(K^i)) #For formula refer Wikipedia page for Hill coefficient.
  
  #Assign values for this iteration to the appropriate column.
  dat[, i+1] <- this_y
}

#Give names to columns.
colnames(dat) <- c("S", 1:50)

#ggplot2 only accepts data frames. Convert matrix object to data.frame.
dat <- as.data.frame(dat)

#Open a pdf file for plotting.
pdf("Hill_equation.pdf")
#Loop through all Hill coefficients.
for (i in 1:50) {
  
  #Get the substrate and corresponding reaction velocity values for plotting.
  this_dat <- dat[, c("S", as.character(i))]
  
  #Rename columns.
  colnames(this_dat) <- c("Substrate", "Rate")
  
  #Generate figure.
  p <- ggplot(this_dat, aes(x = Substrate, y = Rate)) + 
    geom_line() +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 1))+
    theme(panel.border = element_rect(size = c(1,1,1,1), color = "black", 
                                      fill = NA))+
    ggtitle(paste("Hill coefficient =", i))
  
  #Display figure in pdf file.
  print(p)
}

#Close pdf file.
dev.off()
