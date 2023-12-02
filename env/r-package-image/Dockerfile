FROM rocker/r-ver:4.2.1

MAINTAINER christopher.tastad@mssm.edu

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    cmake \
    pandoc \
    libcurl4-openssl-dev \
    libxml2-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libglpk-dev \
    libcairo2-dev \
    libxt-dev \
    libproj-dev \
    libgdal-dev \
    libgeos-dev \
    libhdf5-dev \
    xclip \
    && rm -rf /var/lib/apt/lists/* \
    && R -q -e 'install.packages("curl")'

VOLUME "${RENV_PATHS_CACHE_HOST}:${RENV_PATHS_CACHE_CONTAINER}"
ENV RENV_VERSION 0.15.5

RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

WORKDIR /dep
COPY renv.lock renv.lock

ENV RENV_PATHS_LIBRARY renv/library
RUN R -e "renv::restore()"

