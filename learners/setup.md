---
title: Setup
---

## Software

This lesson uses R and RStudio, together with a few R/Bioconductor packages.

1. Install **R 4.6** or newer from [CRAN](https://cran.r-project.org/).
2. Install **RStudio Desktop** from
   [Posit](https://posit.co/download/rstudio-desktop/).

## Packages

We install the packages ourselves, step by step, so that you see where
each of them comes from. Run the commands below in the R console. The
installation can take several minutes. If R asks whether to update other
packages, you can answer `a` (all).

### Package installers

Bioconductor packages are installed with **BiocManager**, and packages that
are not yet in Bioconductor are installed from GitHub with **remotes**. Both
come from CRAN:

```r
install.packages(c("BiocManager", "remotes"))
```

The lesson needs the current **Bioconductor release**. Check your version
with:

```r
BiocManager::version()
```

### Bioconductor and CRAN packages

`BiocManager::install()` installs packages from both Bioconductor and CRAN:

```r
BiocManager::install(c(
    "mia",                  # microbiome analysis and the example data
    "miaViz",               # microbiome visualization
    "MultiAssayExperiment", # container for multi-omics data
    "scater",               # PCA
    "UpSetR",               # UpSet plots of the samples
    "MOFA2",                # multi-omics factor analysis
    "mixOmics",             # DIABLO
    "randomForest",         # random forest, used by IntegratedLearner
    "ggplot2"               # plots
))
```

MOFA2 uses a few more packages for some of its plots and tests:

```r
BiocManager::install(c("psych", "corrplot", "pheatmap"))
```

### GitHub packages

The [multiomics](https://github.com/microbiome/multiomics) and
[IntegratedLearner](https://github.com/himelmallick/IntegratedLearner)
packages are not yet in Bioconductor, so we install them from GitHub:

```r
remotes::install_github("microbiome/multiomics")
remotes::install_github("himelmallick/IntegratedLearner")
```

### Check the installation

Load each package. If a package was not installed, `library()` gives an
error that names it:

```r
library(mia)
library(miaViz)
library(MultiAssayExperiment)
library(scater)
library(UpSetR)
library(MOFA2)
library(mixOmics)
library(multiomics)
library(IntegratedLearner)
```

## Set up MOFA2

MOFA2 trains its models in Python. The first time you train a model, it
installs its own Python environment. You do not need to install Python
yourself, but this step can take **10 minutes or more**, so please do it
before the workshop. Run this code, which trains a small model on example
data from MOFA2:

```r
library(MOFA2)
model <- create_mofa(make_example_data()$data)
model <- prepare_mofa(model)
model <- run_mofa(model, outfile = tempfile(fileext = ".hdf5"),
                  use_basilisk = TRUE)
```

When it prints `Training finished`, you are ready.

::::::::::::::::::::::::::::::::::::: callout

## Problems with the installation?

If a package fails to install, check that your R version is recent enough with
`R.version.string`, and try the installation again. If the problem persists,
ask for help before the workshop, so that we can solve it together.

:::::::::::::::::::::::::::::::::::::::::::::::

## Project

Create a new RStudio project for the lesson (*File* → *New Project* → *New
Directory* → *New Project*). The data comes with the mia package, so there
is nothing to download.
