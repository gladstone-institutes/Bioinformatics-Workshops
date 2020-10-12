## Run this R script prior to the workshop. You can open it in RStudio and click "Source". 
## These pacakges require R version 3.6.0 or later. Check by typing: R.version.string
##
## Open this script in RStudio and step through it line-by-line, reading the comments as you go.
##
## Please confirm the "Success" message prior to attending the workshop.  If the script produced an "Error" 
## message, even after repeated attempts, please contact the instructor.

############
## PART 1 ##
############
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

load.libs <- c(
  "ggplot2",
  "EnhancedVolcano",
  "DOSE",
  "GO.db",
  "GSEABase",
  "org.Hs.eg.db",
  "clusterProfiler",
  "dplyr",
  "tidyr",
  "stringr",
  "RColorBrewer",
  "rWikiPathways",
  "RCy3")
options(install.packages.check.source = "no")
options(install.packages.compile.from.source = "never")
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(load.libs, update = TRUE, character.only = TRUE)
status <- sapply(load.libs,require,character.only = TRUE)
if(all(status)){
    print("SUCCESS: You have successfully installed and loaded all required libraries.")
} else{
    cat("ERROR: One or more libraries failed to install correctly. Check the following list for FALSE cases and try again...\n\n")
    status
}

## Troubleshooting Tips:
## * If you see errors related to 'data.table', 'units' or other missing packages, try these lines for each:
##      remove.packages("data.table")
##      install.packages("data.table")
##  Then run the entire script again.
## * If you continue to see problems. You may try restarting R and RStudio. Then trying again.
## * If you still have issues, please email the instructor as soon as possible before the workshop.

############
## PART 2 ##
############

## We will also be connecting to Cytoscape from R. You will also need to install and launch Cytoscape:
#
# * Download the latest Cytoscape from http://www.cytoscape.org/download.php
# * Complete installation wizard
# * Launch Cytoscape

# Run this line to confirm a successful connection with Cytoscape
cytoscapePing() 

# Then run these lines to install a few apps in Cytoscape
installApp('WikiPathways') 
installApp('CyTargetLinker') 
installApp('stringApp') 
