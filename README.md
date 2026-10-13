# Multi-omics microbiome data analysis with R/Bioconductor

A hands-on lesson on multi-omics microbiome data analysis with
R/Bioconductor, built with [The Carpentries Workbench][workbench]. The lesson
is in development (pre-alpha).

View the lesson at <https://microbiome.github.io/multiomics-lesson/>.

## About the lesson

Using data from a diet experiment in rats
([Hintikka et al. 2021](https://doi.org/10.3390/ijerph18084049)), learners
analyse the gut microbiota together with metabolites and host biomarkers
measured from the same samples. They store the data in a
`MultiAssayExperiment` and go through the main groups of integration methods:
exploratory, unsupervised and supervised. The lesson uses the packages of the
[OMA framework](https://microbiome.github.io/) together with other
Bioconductor packages for multi-omics. It takes about 7 hours and follows the
online book
[Orchestrating Microbiome Analysis (OMA)](https://bioconductor.org/books/release/OMA/).

| Episode | Topic |
|---|---|
| 01 | Introduction to multi-omics |
| 02 | `MultiAssayExperiment` and the study data |
| 03 | Exploratory methods: Mantel test, MiRKAT and pairwise associations |
| 04 | Unsupervised methods: joint RPCA and MOFA |
| 05 | Supervised methods: IntegratedLearner and DIABLO |
| 06 | Summary |

The lesson assumes that learners can analyse a single microbiome table, for
example from the
[microbiome lesson](https://microbiome.github.io/microbiome-lesson/).

## Building the lesson locally

The lesson is built with the [sandpaper][sandpaper] package. The R packages
it needs are recorded with renv in `renv/profiles/lesson-requirements/`. To
preview the lesson, run in R from the lesson folder:

```r
install.packages(
    c("sandpaper", "varnish", "pegboard", "tinkr"),
    repos = c("https://carpentries.r-universe.dev", getOption("repos"))
)
sandpaper::build_lesson()
```

The lesson data come with the mia package as the `HintikkaXOData` dataset,
so there is nothing to download.

## Contributing

Contributions are welcome. See [CONTRIBUTING.md](CONTRIBUTING.md) and our
[Code of Conduct](CODE_OF_CONDUCT.md).

## License

The lesson is licensed under [CC BY 4.0](LICENSE.md).

## Authors

Authors of this lesson are listed in [CITATION.cff](CITATION.cff).

## Current maintainers

- Tuomas Borman (<https://github.com/TuomasBorman>)
- Leo Lahti (<https://github.com/antagomir>)
- Himel Mallick (<https://github.com/himelmallick>)

[workbench]: https://carpentries.github.io/sandpaper-docs/
[sandpaper]: https://carpentries.github.io/sandpaper/
