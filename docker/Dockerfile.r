FROM bioconductor/bioconductor_docker:RELEASE_3_19

LABEL maintainer="OmicSage Contributors"
LABEL description="OmicSage R analysis environment"
LABEL version="0.1.0-dev"

RUN apt-get update && apt-get install -y --no-install-recommends \
        libhdf5-dev \
        libxml2-dev \
        libssl-dev \
        libcurl4-openssl-dev \
        libgit2-dev \
        libfftw3-dev \
        libgsl-dev \
        libglpk-dev \
        python3-pip \
    && rm -rf /var/lib/apt/lists/*

RUN R -e " \
    install.packages(c( \
        'Seurat', 'Signac', 'remotes', 'BiocManager', \
        'hdf5r', 'ggplot2', 'dplyr', 'patchwork', \
        'future', 'Matrix', 'data.table', \
        'jsonlite', 'yaml', 'argparse' \
    ), repos='https://cloud.r-project.org', Ncpus=4L) \
"

RUN R -e " \
    BiocManager::install(c( \
        'SingleR', 'celldex', 'DESeq2', 'edgeR', 'limma', \
        'clusterProfiler', 'org.Hs.eg.db', 'org.Mm.eg.db', \
        'ReactomePA', 'fgsea', 'BayesSpace', 'chromVAR', \
        'motifmatchr', 'JASPAR2024', 'TFBSTools', \
        'BSgenome.Hsapiens.UCSC.hg38', \
        'SummarizedExperiment', 'BiocParallel', \
        'ComplexHeatmap', 'scran', 'scuttle' \
    ), ask=FALSE, update=FALSE) \
"

RUN R -e " \
    remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', upgrade='never'); \
    remotes::install_github('GreenleafLab/ArchR', ref='master', repos=BiocManager::repositories(), upgrade='never'); \
    remotes::install_github('immunogenomics/harmony', upgrade='never'); \
"

RUN R -e "install.packages('SoupX', repos='https://cloud.r-project.org')"

WORKDIR /app
COPY . /app/

ENTRYPOINT ["Rscript"]
CMD ["--version"]
