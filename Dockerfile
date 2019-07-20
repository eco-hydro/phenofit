FROM rocker/verse:latest

## geospatial
MAINTAINER "Dongdong Kong"
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    lbzip2 \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    libnetcdf-dev \
    libssl-dev \
    libudunits2-dev \
    sudo \
    gdebi-core \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    xtail \
    wget \
    tree \
    # libsqlite3-dev \
    # libprotobuf-dev \
    # libfftw3-dev \
    # libgsl0-dev \
    # libgl1-mesa-dev \
    # libglu1-mesa-dev \
    # libhdf4-alt-dev \
    # libhdf5-dev \
    # libjq-dev \
    # liblwgeom-dev \
    # netcdf-bin \
    # protobuf-compiler \
    # tk-dev \
    # unixodbc-dev \
  && install2.r --error \
    DT \
    plotly \
    RColorBrewer \
    mapdata \
    gstat \
    maptools \
    ncdf4 \
    proj4 \
    raster \
    rgdal \
    rgeos \
    sf \
    sp
    # RandomFields \
    # RNetCDF \
    # classInt \
    # deldir \
    # hdf5r \
    # lidR \
    # mapview \
    # rlas \
    # spacetime \
    # spatstat \
    # spdep \
    # tmap \
    # geoR \
    # geosphere \
    ## from bioconductor
    # && R -e "BiocInstaller::biocLite('rhdf5')"

## install shiny
# Download and install shiny server
RUN wget --no-verbose https://download3.rstudio.org/ubuntu-14.04/x86_64/VERSION -O "version.txt" && \
    VERSION=$(cat version.txt)  && \
    wget --no-verbose "https://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-$VERSION-amd64.deb" -O ss-latest.deb && \
    gdebi -n ss-latest.deb && \
    rm -f version.txt ss-latest.deb && \
    . /etc/environment && \
    R -e "install.packages(c('shiny', 'rmarkdown'), repos='$MRAN')" && \
    cp -R /usr/local/lib/R/site-library/shiny/examples/* /srv/shiny-server/
# COPY shiny-server.sh /usr/bin/shiny-server.sh
EXPOSE 3838

RUN installGithub.r kongdd/Ipaper \
    && R -e "devtools::install_github('kongdd/phenofit', ref='master')" \
    && rm -rf /tmp/downloaded_packages/
    # kongdd/phenofit \
