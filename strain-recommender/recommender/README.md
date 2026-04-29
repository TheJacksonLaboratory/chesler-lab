# recommender

Computational Tool for Identifying Biologically Relevant Inbred Mouse Strains

## Installation

```r
# Install devtools if needed
install.packages("devtools")

# Install from GitHub
devtools::install_github("TheJacksonLaboratory/chesler-lab", 
                          subdir = "strain-recommender/recommender")
```

## Package Overview

A computational tool that helps researchers identify which inbred 
mouse strains are most biologically relevant to a disease of interest. 
Rather than relying solely on traditional phenotypic similarity, this tool
uses transcriptomic data to match a disease-associated gene expression 
signature against predicted gene expression profiles for 657 inbred mouse
strains, enabling more informed strain selection for disease modeling.

## Quick Start

The main inputs are:

- **`gs`** — a data frame of gene scores representing a disease-associated gene 
  expression signature, with columns `gene_id` (mouse Ensembl IDs) and `score`
- **`meta`** — a data frame of sample metadata with columns `mouse_id`, `sex`, 
  and `strain`
- **`df`** — a data frame of gene expression data with `gene_id` as a column 
  and samples as remaining columns

```r
library(recommender)

# Gene scores: disease-associated signature
gs <- data.frame(
  gene_id = c("ENSMUSG00000000168", "ENSMUSG00000000544"),
  score   = c(1.5, -2.3)
)

# Sample metadata
meta <- data.frame(
  mouse_id = paste0("sample", 1:10),
  sex      = rep("both", 10),
  strain   = paste0("strain", LETTERS[1:10])
)

# Expression data
df <- data.frame(
  gene_id = gs$gene_id,
  matrix(rnorm(20, mean = 5, sd = 1.5), nrow = 2, ncol = 10,
         dimnames = list(NULL, paste0("sample", 1:10)))
)

# Run recommender
result <- recommend(gs = gs, meta = meta, df = df, B = 4000)

# Sort by vulnerability score
result[order(result$score, decreasing = TRUE), ]
```

## Learn More
See the [full vignette](vignettes/strain-recommender-workflow.Rmd) for advanced usage including multiple gene sets and reproducibility options.

## Dependencies

### Required
- [`dplyr`](https://dplyr.tidyverse.org/)
- [`magrittr`](https://magrittr.tidyverse.org/)

### Suggested
- `data.table` — for efficient data loading in vignettes
- `knitr` + `rmarkdown` — for building vignettes

## Citation

If you use this package, please cite the associated paper (DOI forthcoming)

To cite the package itself:

> Ball RL, Klein A, Chesler EJ (2026). *recommender*: Computational Tool for 
> Identifying Biologically Relevant Inbred Mouse Strains. R package version 1.0.0.
> https://github.com/TheJacksonLaboratory/chesler-lab

## License
This package is licensed under CC BY-NC 4.0. See [LICENSE.md](LICENSE.md) for details.
