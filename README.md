# Bisctree: Bayesian integrative somatic clonal tree

This R package contains methods to infer somatic clonal tree by integrating 
bulk DNA-seq and single-cell RNA- or DNA-seq.

## Installation

The `bisctree` package can be conveniently installed using the [`devtools`](https://www.rstudio.com/products/rpackages/devtools/) package thus:

```{R}
devtools::install_github("huangyh09/bisctree")
```

## Test and development

This package is still under fast development. Currently, we are implementing 
MCMC algorithms for clonal tree inference from bulk sample. Example tests 
including input data are located in folder [inst](https://github.com/huangyh09/bisctree/tree/master/inst).
