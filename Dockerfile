# syntax=docker/dockerfile:1.6
#
# QPX container image
# -------------------
# Provides the ``qpxc`` CLI built from this repository's source, including
# MuData export support.
#
# Build:
#   docker build -t ghcr.io/bigbio/qpx:dev .
#
# Run:
#   docker run --rm -v $(pwd):/data ghcr.io/bigbio/qpx:dev convert diann --help

FROM python:3.11-slim-bookworm

# Runtime system deps required by pyOpenMS (see environment.yml).
# hadolint ignore=DL3008
RUN apt-get update \
 && apt-get install -y --no-install-recommends libglib2.0-0 procps \
 && rm -rf /var/lib/apt/lists/*

# hatch-vcs derives the version from git history; when the build context lacks
# a .git directory (typical Docker builds), fall back to this placeholder.
ARG QPX_VERSION=0.0.0.dev0
ENV SETUPTOOLS_SCM_PRETEND_VERSION=${QPX_VERSION}

# Minimal build inputs required by ``pip install .`` under hatchling.
WORKDIR /src
COPY pyproject.toml README.md LICENSE ./
COPY qpx ./qpx

RUN pip install --no-cache-dir --upgrade pip==24.0 \
 && pip install --no-cache-dir . \
 && rm -rf /src

LABEL org.opencontainers.image.source="https://github.com/bigbio/qpx"
LABEL org.opencontainers.image.description="QPX: Quantitative Proteomics Parquet toolkit (qpxc CLI) with MuData export"
LABEL org.opencontainers.image.licenses="Apache-2.0"

WORKDIR /data
