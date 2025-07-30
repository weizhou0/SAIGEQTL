FROM ubuntu:20.04

LABEL maintainer="wzhou@broadinstitute.org"

ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8
ENV OMP_NUM_THREADS=1
# Force step_2 to use 1 single thread. More threads are ineffective
ENV PIXI_NO_PATH_UPDATE=1

RUN apt-get update -qq && \
    apt-get -y install curl time \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY . .

RUN curl -fsSL https://pixi.sh/install.sh | sh && \
    mv /root/.pixi/bin/pixi /bin && pixi install && \
    rm -rf /root/.cache && \
    pixi run Rscript -e 'install.packages("lintools", repos="https://cloud.r-project.org")' && \
    pixi run Rscript -e 'install.packages("remotes", repos="https://cloud.r-project.org"); remotes::install_github("barkasn/fastSave")'

#RUN pixi run R CMD INSTALL .
RUN pixi run Rscript -e 'install.packages(".", repos=NULL, type="source")'

# Create Singularity-compatible wrapper scripts
RUN for script in step1_fitNULLGLMM_qtl.R step2_tests_qtl.R step3_gene_pvalue_qtl.R makeGroupFile.R; do \
      echo '#!/bin/bash' > /usr/local/bin/$script && \
      echo 'export PATH="/root/.pixi/envs/default/bin:$PATH"' >> /usr/local/bin/$script && \
      echo 'export LD_LIBRARY_PATH="/root/.pixi/envs/default/lib:$LD_LIBRARY_PATH"' >> /usr/local/bin/$script && \
      echo 'exec /root/.pixi/envs/default/bin/Rscript /app/extdata/'$script' "$@"' >> /usr/local/bin/$script && \
      chmod +x /usr/local/bin/$script; \
    done && \
    step1_fitNULLGLMM_qtl.R --help && \
    step2_tests_qtl.R --help && \
    step3_gene_pvalue_qtl.R --help && \
    makeGroupFile.R --help

RUN rm -r src/
RUN rm -r R/
