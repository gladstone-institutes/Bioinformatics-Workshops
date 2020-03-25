library(ggplot2)

#Following library will be used to reformat data for plotting.
library(reshape2)

#Load the data for plotting.
load("nucleotide_resolution.RData")

#Open a pdf file for plotting.
pdf("Nucleotide_resolution.pdf")

#Loop through all genes.
for (i in 1:length(dat)) {
  
  #Get data for this gene.
  this_dat <- dat[[i]]
  
  #Reformat data for plotting.
  this_dat <- melt(this_dat, id.vars = "Position")
  
  #Store the plot in a variable.
  p <- ggplot(this_dat, aes(x = Position, y = value)) +
    geom_bar(stat = "identity") +
    facet_grid(variable ~ .)+
    scale_y_continuous(breaks = c(0, 0.5, 1))+
    ylab("Signal")+
    annotate(geom = "rect",
             xmin = to_annotate$Start[i],
             xmax = to_annotate$End[i], 
             ymin = -Inf,
             ymax = Inf, fill = "red",
             color = NA, 
             alpha = 0.3) +
    ggtitle(names(dat)[i])
  
  #Print plot to file.
  print(p)
}

#Close file.
dev.off()