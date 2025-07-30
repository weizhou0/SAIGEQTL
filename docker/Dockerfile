FROM ubuntu:22.04

LABEL maintainer="wzhou@broadinstitute.org" \
      description="SAIGE-QTL: Scalable and accurate eQTL mapping for single-cell studies" \
      version="0.3.2"

# Set environment variables
ENV LC_ALL=C.UTF-8 \
    LANG=C.UTF-8 \
    OMP_NUM_THREADS=1 \
    PIXI_NO_PATH_UPDATE=1 \
    DEBIAN_FRONTEND=noninteractive \
    TZ=UTC

# Install system dependencies in a single layer
RUN apt-get update -qq && \
    apt-get -y --no-install-recommends install \
        curl \
        time \
        ca-certificates \
        wget \
        build-essential \
        tzdata && \
    ln -fs /usr/share/zoneinfo/$TZ /etc/localtime && \
    dpkg-reconfigure -f noninteractive tzdata && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Set working directory
WORKDIR /app

# Copy dependency files first for better layer caching
COPY pixi.toml pixi.lock ./

# Install pixi and dependencies in a single layer
RUN curl -fsSL https://pixi.sh/install.sh | sh && \
    mv /root/.pixi/bin/pixi /usr/local/bin/ && \
    pixi install && \
    rm -rf /root/.cache /root/.pixi/cache

# Install GitHub packages not available in conda
RUN pixi run install-github-deps

# Copy source code
COPY . .

# Build and install the package
RUN pixi run build && \
    pixi run test

# Create wrapper scripts that work in both Docker and Singularity
RUN mkdir -p /usr/local/bin && \
    echo '#!/bin/bash' > /usr/local/bin/step1_fitNULLGLMM_qtl.R && \
    echo 'export PATH="/app/.pixi/envs/default/bin:$PATH"' >> /usr/local/bin/step1_fitNULLGLMM_qtl.R && \
    echo 'export LD_LIBRARY_PATH="/app/.pixi/envs/default/lib:$LD_LIBRARY_PATH"' >> /usr/local/bin/step1_fitNULLGLMM_qtl.R && \
    echo 'exec /app/.pixi/envs/default/bin/Rscript /app/extdata/step1_fitNULLGLMM_qtl.R "$@"' >> /usr/local/bin/step1_fitNULLGLMM_qtl.R && \
    chmod +x /usr/local/bin/step1_fitNULLGLMM_qtl.R && \
    \
    echo '#!/bin/bash' > /usr/local/bin/step2_tests_qtl.R && \
    echo 'export PATH="/app/.pixi/envs/default/bin:$PATH"' >> /usr/local/bin/step2_tests_qtl.R && \
    echo 'export LD_LIBRARY_PATH="/app/.pixi/envs/default/lib:$LD_LIBRARY_PATH"' >> /usr/local/bin/step2_tests_qtl.R && \
    echo 'exec /app/.pixi/envs/default/bin/Rscript /app/extdata/step2_tests_qtl.R "$@"' >> /usr/local/bin/step2_tests_qtl.R && \
    chmod +x /usr/local/bin/step2_tests_qtl.R && \
    \
    echo '#!/bin/bash' > /usr/local/bin/step3_gene_pvalue_qtl.R && \
    echo 'export PATH="/app/.pixi/envs/default/bin:$PATH"' >> /usr/local/bin/step3_gene_pvalue_qtl.R && \
    echo 'export LD_LIBRARY_PATH="/app/.pixi/envs/default/lib:$LD_LIBRARY_PATH"' >> /usr/local/bin/step3_gene_pvalue_qtl.R && \
    echo 'exec /app/.pixi/envs/default/bin/Rscript /app/extdata/step3_gene_pvalue_qtl.R "$@"' >> /usr/local/bin/step3_gene_pvalue_qtl.R && \
    chmod +x /usr/local/bin/step3_gene_pvalue_qtl.R && \
    \
    echo '#!/bin/bash' > /usr/local/bin/makeGroupFile.R && \
    echo 'export PATH="/app/.pixi/envs/default/bin:$PATH"' >> /usr/local/bin/makeGroupFile.R && \
    echo 'export LD_LIBRARY_PATH="/app/.pixi/envs/default/lib:$LD_LIBRARY_PATH"' >> /usr/local/bin/makeGroupFile.R && \
    echo 'exec /app/.pixi/envs/default/bin/Rscript /app/extdata/makeGroupFile.R "$@"' >> /usr/local/bin/makeGroupFile.R && \
    chmod +x /usr/local/bin/makeGroupFile.R

# Verify installation and create test
RUN step1_fitNULLGLMM_qtl.R --help > /dev/null && \
    step2_tests_qtl.R --help > /dev/null && \
    step3_gene_pvalue_qtl.R --help > /dev/null && \
    makeGroupFile.R --help > /dev/null && \
    echo "All executables verified successfully"

# Clean up source files to reduce image size
RUN rm -rf src/ R/ man/ && \
    rm -rf /root/.pixi/cache /root/.cache && \
    rm -rf /tmp/* /var/tmp/*

# Set default command
CMD ["/bin/bash"]
