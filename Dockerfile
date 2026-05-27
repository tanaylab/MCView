# MCView Docker image.
#
# NOTE: build/run on a reasonably current Docker (>= 20.10.10). The image is
# Ubuntu-jammy based; older Docker seccomp profiles block syscalls its glibc/R
# need (symptoms: apt "NO_PUBKEY" at build time, "R_HOME not found" at runtime).
FROM rocker/r-ver:4.3.3

# MCView git ref to install (branch / tag / commit). Defaults to the latest master.
ARG MCVIEW_REF=master

# Posit Package Manager repo used to resolve R dependencies as binaries (fast,
# no compilation). Defaults to the latest snapshot; pass a dated snapshot URL
# (e.g. .../jammy/2026-04-10) via --build-arg for a fully reproducible build.
# NOTE: as of MCView 0.2.55 `qs` is an optional (Suggests) dependency, so the
# old constraint of staying before stringfish 0.19.0 (which broke qs's build) no
# longer applies.
ARG CRAN_SNAPSHOT=https://p3m.dev/cran/__linux__/jammy/latest

# System libraries needed by MCView's dependency tree
# (Shiny + ggplot/ragg/systemfonts + git2/gert + curl/openssl/xml2 + image libs).
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev libssl-dev libxml2-dev libgit2-dev libssh2-1-dev \
    libicu-dev libpng-dev libjpeg-dev libtiff-dev \
    libfontconfig1-dev libfreetype6-dev libharfbuzz-dev libfribidi-dev \
    zlib1g-dev libbz2-dev liblzma-dev \
    make cmake git pandoc \
    python3 python3-pip python3-venv python3-dev \
  && rm -rf /var/lib/apt/lists/*

# Python anndata in an isolated venv. The R 'anndata' package reads h5ad files
# through reticulate; point reticulate at this venv so it always finds anndata.
RUN python3 -m venv /opt/venv && \
    /opt/venv/bin/pip install --no-cache-dir --upgrade pip && \
    /opt/venv/bin/pip install --no-cache-dir anndata
ENV RETICULATE_PYTHON=/opt/venv/bin/python

# Resolve all R dependencies from the pinned snapshot.
RUN echo "options(repos = c(CRAN = '${CRAN_SNAPSHOT}'))" >> /usr/local/lib/R/etc/Rprofile.site

# Install MCView. pak resolves its GitHub Remotes (e.g. tgutil) and errors out on
# failure; the requireNamespace check is a hard backstop, because
# remotes::install_github silently returns success even when a dependency fails
# to build -- which would otherwise tag an image with no MCView in it.
RUN R -e 'install.packages("pak")'
RUN R -e "pak::pak('tanaylab/MCView@${MCVIEW_REF}')" \
 && R -e 'if (!requireNamespace("MCView", quietly = TRUE)) quit(status = 1); cat("MCView", as.character(packageVersion("MCView")), "installed OK\n")'

RUN mkdir /MCView
ADD ./start_app.R /MCView/

WORKDIR /MCView

EXPOSE 3838
CMD ["./start_app.R", "/project", "3838"]
