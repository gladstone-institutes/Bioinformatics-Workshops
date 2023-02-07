FROM rocker/rstudio:4

RUN apt-get update && apt-get install -y \
    software-properties-common
RUN add-apt-repository universe

RUN apt-get update && apt-get install -y \
    apache2 \
    curl \
    git \
    libgsl-dev 

RUN apt-get update && apt-get install -y \
	zlib1g-dev





#RUN apt-get update && apt-get install -y \
#	python3-setuptools


RUN mkdir /etc/python3.8.16
WORKDIR /etc/python3.8.16
RUN curl https://www.python.org/ftp/python/3.8.16/Python-3.8.16.tgz --output Python-3.8.16.tgz
RUN tar xzf Python-3.8.16.tgz
RUN cd Python-3.8.16
RUN ./Python-3.8.16/configure


RUN make altinstall


ENV PATH "/etc/python3.8.16/:$PATH"

RUN rm -rf /usr/local/lib/python3.10
RUN rm /usr/bin/python3.10




RUN mkdir /etc/code
WORKDIR /etc/code
COPY MACS2-2.2.7.1 /etc/code/MACS2-2.2.7.1
WORKDIR /etc/code/MACS2-2.2.7.1
ENV PATH "/etc/code/MACS2-2.2.7.1/bin:$PATH"
ENV PYTHONPATH "/etc/code//local/lib/python3.8/dist-packages/:$PYTHONPATH"

RUN curl https://bootstrap.pypa.io/pip/get-pip.py --output get-pip.py
RUN python3.8 get-pip.py


RUN pip install numpy
RUN pip install --upgrade setuptools
RUN python setup.py install --prefix=/etc/code/



RUN apt-get update && apt-get install -y \
	liblzma-dev \
	libbz2-dev 


RUN apt-get update && apt-get install -y \
	libglpk-dev


#install specific versions of R packages
RUN R -e 'install.packages("remotes")'
RUN R -e 'install.packages("BiocManager")'
RUN R -e 'BiocManager::install(version = "3.16",ask = FALSE)'
RUN R -e 'BiocManager::install("glmGamPoi", version="3.16")'
RUN R -e 'library("remotes");install_version("ggplot2", "3.3.5")'
RUN R -e 'library("remotes");install_version("dplyr", "1.0.7")'
RUN R -e 'library("remotes");install_version("magrittr", "2.0.2")'
RUN R -e 'library("remotes");install_version("lme4", "1.1-28")'
RUN R -e 'library("remotes");install_version("gdata", "2.18.0")'
RUN R -e 'library("remotes");install_version("devtools", "2.4.3")'
RUN R -e 'devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())'
RUN R -e 'ArchR::installExtraPackages()'
RUN R -e 'BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", version="3.16")'
RUN R -e 'library("remotes");install_version("Seurat", "4.3.0")'

RUN cp /etc/code/MACS2-2.2.7.1/bin/macs2 /usr/local/bin/macs2
RUN cp -R /etc/code/lib/python3.8/site-packages/MACS2-2.2.7.1-py3.8-linux-x86_64.egg /usr/local/lib/python3.8/site-packages/

RUN cp -R /etc/code/lib/python3.8/site-packages/MACS2-2.2.7.1-py3.8-linux-x86_64.egg/MACS2 /usr/local/lib/python3.8/site-packages/

ENV PYTHONPATH "/usr/local/lib/python3.8/site-packages/:/etc/code/lib/python3.8/site-packages/:${PYTHONPATH}"

RUN R -e 'BiocManager::install("chromVAR", version="3.16")'
RUN R -e 'devtools::install_github("GreenleafLab/chromVARmotifs", ref="master")'
