# let's install all the packages we will need for this training

# package list
packages <- c(
  "SNPRelate", "gdsfmt", "snpReady", "impute", "rgl", "stringr", "matrixcalc",
  "bestNormalize", "desplot", "car", "carData", "lme4", "Matrix", "PerformanceAnalytics",
  "xts", "zoo", "scales", "ASRgenomics", "superheat", "scatterplot3d", "pedigreemm",
  "factoextra", "cluster", "lubridate", "forcats", "dplyr", "purrr", "readr",
  "tidyr", "tibble", "ggplot2", "tidyverse", "reshape2", "breedR", "sp",
  "gaston", "RcppParallel", "Rcpp", "knitr", "CMplot", "EMMREML", "ape", "genetics",
  "mvtnorm", "MASS", "gtools", "gdata", "combinat", "gplots", "biganalytics",
  "biglm", "DBI", "foreach", "bigmemory", "GAPIT", "PCAtools", "ggrepel",
  "plotly", "ggpubr", "doMC", "doParallel", "iterators", "remotes", "metan",
  "statgenGxE", "devtools", "AlphaSimR", "R6"
)

# Remove duplicates
packages <- unique(packages)

# Install missing packages from CRAN
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Apply installation function
invisible(sapply(packages, install_if_missing))

# packages from other sources than CRAN
BiocManager::install("impute")
devtools::install_github(repo = 'italo-granato/snpReady', ref = 'dev')
remotes::install_github("Biometris/statgenGxE", ref = "develop", dependencies = TRUE)
devtools::install_github('famuvie/breedR')
