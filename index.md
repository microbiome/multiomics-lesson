---
site: sandpaper::sandpaper_site
---

This lesson teaches how to analyse microbiome data together with other data
types, such as metabolites and host biomarkers, measured from the same
samples. Using data from a diet experiment in rats, you will learn to store
several data tables in one object, to test whether they are related, to find
the shared patterns in them, and to predict an outcome from all of them.

The lesson assumes that you can already analyse a single microbiome table,
for example from the
[microbiome lesson](https://microbiome.github.io/microbiome-lesson/) (see the
prerequisites below). It uses the [OMA framework](https://microbiome.github.io/), a family of
R/Bioconductor packages for microbiome data science, together with other
Bioconductor packages for multi-omics. After the lesson, you can continue
learning with the online book
[Orchestrating Microbiome Analysis (OMA)](https://microbiome.github.io/OMA/).

## Learning objectives

After this lesson, learners will be able to:

| Learning objective | Episode |
|---|---|
| Explain what multi-omics data is, why the microbiome is studied together with other data types, where the data comes from, and which questions the main groups of integration methods answer. | [Introduction](episodes/01-introduction.Rmd) |
| Store, explore, subset and prepare several data types of the same samples in a `MultiAssayExperiment`. | [MultiAssayExperiment and the study data](episodes/02-multiassayexperiment.Rmd) |
| Test whether two data types are related, and find the features that link them, with the Mantel test, MiRKAT and pairwise associations. | [Exploratory methods](episodes/03-exploratory.Rmd) |
| Find the main patterns that several data types share with joint RPCA and MOFA. | [Unsupervised methods](episodes/04-unsupervised.Rmd) |
| Predict an outcome from several data types and find the features that predict it with IntegratedLearner and DIABLO. | [Supervised methods](episodes/05-supervised.Rmd) |

## Where does this lesson fit?

A microbiome study goes from collecting samples and sequencing them, through
processing the reads into tables, to analysing the tables. The analysis of
one table covers the data object, data processing and statistical analysis.
This lesson covers the next step: **data integration**, where the microbiome is
combined with other data types.

![Microbiome data science workflow. This lesson covers the data integration step. Source: [OMA](https://microbiome.github.io/OMA/).](episodes/fig/workflow-integration.png){alt="Microbiome data science workflow with the data integration step highlighted: taxonomic profiles and metabolic profiles are combined in a joint analysis." width="45%"}

## The data

We analyse samples from 40 rats that were fed one of four diets: a high-fat
or a low-fat diet, with or without a prebiotic supplement called
xylo-oligosaccharides (XOS)
([Hintikka et al. 2021](https://doi.org/10.3390/ijerph18084049)). For every
rat, the data contains three tables: the gut microbiota (16S rRNA
sequencing), metabolites (NMR spectroscopy) and host biomarkers. The data
comes with the mia package, so there is nothing to download.

:::::::::::::::::::::::::::::::::::::::::: prereq

## Prerequisites

- The [microbiome lesson](https://microbiome.github.io/microbiome-lesson/),
  or equivalent experience: the `TreeSummarizedExperiment` class,
  agglomeration and transformation, beta diversity and ordination, and the
  basics of machine learning.
- Basic R programming: data frames, factors, functions and simple `for`
  loops.
- Basic statistics: what a p-value and a correlation are.

::::::::::::::::::::::::::::::::::::::::::::::::::
