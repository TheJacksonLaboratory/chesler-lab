# transcriptome.imputation

Transcriptome Imputation Pipeline for Inbred Strains using Diversity Outbred Mice

## Installation

```r
# Install devtools if needed
install.packages("devtools")

# Install from GitHub
devtools::install_github("TheJacksonLaboratory/chesler-lab", 
                          subdir = "strain-recommender/transcriptome-imputation")
```

## Package Overview

Provides a complete pipeline for imputing inbred strain transcriptomes 
from Diversity Outbred (DO) mouse gene expression and genotype probability data. 
Implements preprocessing of metadata, expression, eQTL, and genotype probability 
inputs; ridge regression model training using phylogenetic distance weights; 
and prediction of target strain expression. Supports both local and HPC cluster 
execution via SLURM submission scripts.

## Quick Start

```r
# Load using the package name
library(transcriptome.imputation)

```

Brief example here.

## Dependencies

### Required
- [`data.table`](https://rdatatable.gitlab.io/data.table/)
- [`tidyverse`](https://www.tidyverse.org/)
- [`stringr`](https://stringr.tidyverse.org/)
- [`readr`](https://readr.tidyverse.org/)
- [`dplyr`](https://dplyr.tidyverse.org/)
- [`magrittr`](https://magrittr.tidyverse.org/)
- [`glmnet`](https://glmnet.stanford.edu/)
- [`ape`](https://cran.r-project.org/package=ape)
- [`R.utils`](https://cran.r-project.org/package=R.utils)

### Suggested
- `knitr` + `rmarkdown` — for building vignettes
- `withr` — for testing
- `testthat` — for running tests
