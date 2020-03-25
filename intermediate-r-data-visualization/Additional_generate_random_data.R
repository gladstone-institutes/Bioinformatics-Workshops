dat <- list()


genes <- c("YPL1", "BRCA", "DBP2", "1p22", "RLF", "DCLRE1B",
           "ZNF268", "RNF220", "EPHA2", "SDC3")
to_annotate <- data.frame(Gene = genes, Start = NA, End = NA)

for (i in genes) {
  rows <- sample(50:100, 1)
  dat[[i]] <- matrix(, rows, 6)
  dat[[i]][, 1] <- runif(rows, 0, 1)
  dat[[i]][, 4] <- runif(rows, 0, 1)
  dat[[i]][, 2] <- dat[[i]][, 1] + runif(rows, 0, 0.1)
  dat[[i]][, 3] <- dat[[i]][, 1] + runif(rows, 0, 0.1)
  dat[[i]][, 5] <- dat[[i]][, 2] + runif(rows, 0, 0.1)
  dat[[i]][, 6] <- dat[[i]][, 2] + runif(rows, 0, 0.1)
  dat[[i]][dat[[i]] < 0 ] <- 0
  dat[[i]][dat[[i]] > 1 ] <- 1
  dat[[i]] <- cbind(sample(100:1000, 1) + 1:rows, dat[[i]])
  colnames(dat[[i]]) <- c("Position", "A1", "A2", "A3", "B1", "B2", "B3")
  dat[[i]] <- as.data.frame(dat[[i]])
  
  to_annotate[which(genes == i), "Start"] <- dat[[i]]$Position[1] + sample(1:20, 1)
  to_annotate[which(genes == i), "End"] <- to_annotate[which(genes == i), "Start"] + sample(1:20, 1)
}

save(list = c("dat", "to_annotate"), file = "Detailed_plotting_data.RData")
