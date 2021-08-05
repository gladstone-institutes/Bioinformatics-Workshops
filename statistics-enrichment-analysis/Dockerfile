# Docker inheritance
FROM rocker/verse

RUN R -e 'install.packages("BiocManager")'
# Install required Bioconductor package
RUN R -e 'install.packages("batchtools")'
RUN R -e 'install.packages("tidyverse")'
RUN R -e 'install.packages("rSEA")'
RUN R -e 'install.packages("statmod")'
# Install required Bioconductor package
RUN R -e 'BiocManager::install("clusterProfiler")'
RUN R -e 'BiocManager::install("EnrichmentBrowser")'
RUN R -e 'BiocManager::install("GSEABenchmarkeR")'
RUN R -e 'BiocManager::install("org.Hs.eg.db")'
RUN R -e 'BiocManager::install("GEOquery")'
RUN R -e 'BiocManager::install("GO.db")'
RUN R -e 'BiocManager::install("edgeR")'
RUN R -e 'BiocManager::install("DESeq2")'
RUN R -e 'BiocManager::install("genefilter")'
RUN R -e 'BiocManager::install("geneplotter")'
RUN R -e 'BiocManager::install("ExperimentHub")'
RUN R -e 'BiocManager::install("GSEABase")'
RUN R -e 'BiocManager::install("PADOG")'
RUN R -e 'BiocManager::install("biomaRt")'


