# Converting loops to map with ChatGPT ------------------------------------
# Author: Natalie Elphick

# ** Important: **
# Before running this code, set your working directory to the same location as this file:
# Click Session -> Set Working Directory -> To Source File Location

# Run the code line by line with Ctrl+Enter or Cmd+Enter on Mac

# First we will create a simulated example and then use ChatGPT to convert the 
# loop into tidyverse style code that uses map()

library(tidyverse)

# Gene names
gene_names <- paste("Gene", 1:10)

# Simulate some data
data <- data.frame(
  Gene = rep(gene_names, each = 20),
  Treatment = rep(rep(c("Treatment", "Control"), each = 10), times = length(gene_names)),  
  Sample = rep(1:20, times = length(gene_names)),  
  Log_CPM = c(
    rnorm(10, mean = 5, sd = 1),
    rnorm(10, mean = 6, sd = 1)
  ) + rnorm(20 * length(gene_names), mean = 0, sd = 1)  # Add sample-specific variation
)

# Create a directory to save plots
dir.create("gene_boxplots", showWarnings = FALSE)

# Iterate through gene names and create boxplots
for (gene_name in gene_names) {
  # Subset data for the current gene
  gene_data <- data[data$Gene == gene_name, ]
  
  # Create a boxplot comparing treatment vs. control for Log CPM using ggplot2
  plot_filename <- paste("gene_boxplots/", gene_name, "_boxplot.png", sep = "")
  p <- ggplot(data = gene_data, aes(x = Treatment, y = Log_CPM, fill = Treatment)) +
    geom_boxplot() +
    labs(
      title = paste("Log CPM for", gene_name),
      x = "Treatment",
      y = "Log CPM"
    ) +
    scale_fill_manual(values = c("#009E73", "#D55E00")) +
    theme_minimal()
  ggsave(plot = p,plot_filename, width = 6, height = 4)
  
  # Print the plot filename
  cat("Boxplot saved as:", plot_filename, "\n")
}



# Converting the syntax ---------------------------------------------------

# Paste the following prompt into ChatGPT:
# Convert the following code into tidyverse style code that uses purr::map() :
# <paste in the code above>

# Put the result below and run it










# Here is what I got ------------------------------------------------------

library(purrr) # *This part was incorrect, it should still be library(tidyverse)*

# Gene names
gene_names <- paste("Gene", 1:10)

# Simulate some data
data <- tibble(
  Gene = rep(gene_names, each = 20),
  Treatment = rep(rep(c("Treatment", "Control"), each = 10), times = length(gene_names)),  
  Sample = rep(1:20, times = length(gene_names)),  
  Log_CPM = c(
    rnorm(10, mean = 5, sd = 1),
    rnorm(10, mean = 6, sd = 1)
  ) + rnorm(20 * length(gene_names), mean = 0, sd = 1)  # Add sample-specific variation
)

# Create a directory to save plots
dir.create("gene_boxplots", showWarnings = FALSE)

# Function to create and save a boxplot for a given gene
save_gene_plot <- function(gene_name) {
  gene_data <- data %>% filter(Gene == gene_name)
  plot_filename <- paste("gene_boxplots/", gene_name, "_boxplot.png", sep = "")
  p <- ggplot(data = gene_data, aes(x = Treatment, y = Log_CPM, fill = Treatment)) +
    geom_boxplot() +
    labs(
      title = paste("Log CPM for", gene_name),
      x = "Treatment",
      y = "Log CPM"
    ) +
    scale_fill_manual(values = c("#009E73", "#D55E00")) +
    theme_minimal()
  ggsave(plot = p, filename = plot_filename, width = 6, height = 4)
  
  # Print the plot filename
  cat("Boxplot saved as:", plot_filename, "\n")
}

# Iterate through gene names and create boxplots using map
map(gene_names, save_gene_plot)




