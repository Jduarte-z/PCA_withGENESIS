# PCA_withGENESIS
Pipeline to compute PCA with the R package GENESIS, with plotting and diagnostic logic implemented

For installing genesis and related needed packages using conda (assuming here that you have already downloaded miniconda3 in your machine):

```
#get things ready to run PC-Air and PC-Relate 
conda create -n genesis_pcair \
  -c conda-forge -c bioconda \
  r-base=4.5 \
  r-data.table \
  r-tidyverse \
  r-optparse \
  bioconductor-genesis \
  bioconductor-snprelate \
  bioconductor-gdsfmt \
  bioconductor-seqarray \
  bioconductor-seqvartools \
  bioconductor-gwastools \
  plink2 \
  bcftools \
  htslib \
  -y

#activate your environment
conda activate genesis 
```
