

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

conda create --name transcriptomics fastp cutadapt trimmomatic hisat2 subread fastqc
