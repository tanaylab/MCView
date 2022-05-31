FROM rocker/r-ver:4.2.0
RUN apt-get update && apt-get install -y  git-core libcurl4-openssl-dev libgit2-dev libicu-dev libpng-dev libssl-dev libxml2-dev make pandoc pandoc-citeproc python zlib1g-dev && rm -rf /var/lib/apt/lists/*
RUN echo "options(repos = c(CRAN = 'https://cran.rstudio.com/'), download.file.method = 'libcurl')" >> /usr/local/lib/R/etc/Rprofile.site

# install python
RUN apt-get update && \
    apt-get install -y software-properties-common && \
    add-apt-repository ppa:deadsnakes/ppa && \
    apt-get update && \
    apt-get -y install python3.9

# install pip
RUN apt-get install -y python3-pip

# install anndata
RUN pip3 install anndata

# Install MCView
RUN R -e 'install.packages("remotes")'
RUN R -e 'remotes::install_github("tanaylab/MCView", Ncpus = 40)'

RUN mkdir /MCView
ADD ./start_app.R /MCView/

WORKDIR /MCView

EXPOSE 3838
CMD ["./start_app.R", "/project", "3838"]