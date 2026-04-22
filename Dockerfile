FROM rocker/shiny-verse:4.5.3
LABEL authors="Alex Lemenze" \
    description="Docker image for MaGIC Gene Set Enrichment Tool"

# ── System dependencies ────────────────────────────────────────────────────
RUN apt-get update && apt-get install -y \
    sudo \
    libhdf5-dev \
    build-essential \
    libcurl4-gnutls-dev \
    libxml2-dev \
    libssl-dev \
    libv8-dev \
    libsodium-dev \
    libglpk40 \
    libpng-dev \
    libjpeg-dev \
    libtiff-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# ── CRAN packages ──────────────────────────────────────────────────────────
RUN R -e "install.packages(c( \
    'BiocManager', \
    'shiny', \
    'shinyjs', \
    'shinythemes', \
    'shinycssloaders', \
    'shinyWidgets', \
    'DT', \
    'tidyverse', \
    'data.table', \
    'RColorBrewer', \
    'colourpicker', \
    'msigdbr', \
    'viridis' \
    ), repos='https://cran.rstudio.com/', dependencies=TRUE)"

# ── Bioconductor: annotation databases (large, cache separately) ──────────
RUN R -e "BiocManager::install(c( \
    'AnnotationDbi', \
    'org.Hs.eg.db', \
    'org.Mm.eg.db' \
    ), ask=FALSE, update=FALSE)"

# ── Bioconductor: analysis packages ───────────────────────────────────────
RUN R -e "BiocManager::install(c( \
    'DOSE', \
    'clusterProfiler', \
    'enrichplot', \
    'ReactomePA' \
    ), ask=FALSE, update=FALSE)"

# ── Copy application files ─────────────────────────────────────────────────
COPY ./app /srv/shiny-server/
COPY shiny-customized.config /etc/shiny-server/shiny-server.conf

# ── Permissions ────────────────────────────────────────────────────────────
RUN chown -R shiny:shiny /srv/shiny-server

EXPOSE 8080
USER shiny
CMD ["/usr/bin/shiny-server"]
