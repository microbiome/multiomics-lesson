---
title: 'Reference'
---

## Glossary

16S amplicon sequencing
:   Sequencing one marker gene (16S rRNA) to find out which bacteria are in a
    sample. Gives counts of taxa.

AUC (area under the ROC curve)
:   A measure of how well a model separates two classes: 1 is perfect, 0.5 is
    no better than guessing.

Balanced accuracy
:   The average of the share of correctly predicted samples in each class.
    Unlike accuracy, it is not inflated when one class is much larger.

Batch effect
:   A technical difference between samples that were measured on different
    days, runs or plates.

Classification error rate
:   The share of samples that a model assigns to the wrong group. With two
    equally large groups, guessing gives 50 %.

Confounding
:   When two factors vary together (for example the batch and the study
    group), so that their effects cannot be told apart.

Cross-validation
:   Estimating how well a model predicts new samples by repeatedly training
    it on part of the data and testing it on the rest.

Detection limit
:   The smallest amount that an instrument can measure. Values below it are
    often reported as missing.

DIABLO
:   A supervised integration method in the mixOmics package
    (`block.splsda()`). It selects a few features from each data type that
    together separate the groups and are correlated between the data types.

Early fusion
:   Joining the tables of all data types side by side and training one model
    on them. Also called concatenation.

Experiment
:   One data type stored in a `MultiAssayExperiment`, for example the
    microbiota as a `TreeSummarizedExperiment`.

Factor (MOFA)
:   A hidden pattern that varies between the samples and is found by MOFA. A
    factor can be shared by several data types or be specific to one.

Factor model
:   A model in which a few hidden variables (factors), each with a weight per
    feature, explain most of the variation of the measured features.

Intermediate fusion
:   Training a model on each data type, joining what the models learned, and
    training one more model on top.

Holobiont / hologenomics
:   The host and its microbes seen as one biological system / the study of
    their genomes and molecules together.

Joint RPCA
:   Robust principal component analysis of two or more data types together.
    It draws one map of the samples, and fills in missing values, such as the
    zeros of the microbiota, itself.

Functional redundancy
:   Different taxa carrying out the same function, so that communities with
    different taxa can work in the same way.

Kernel
:   A matrix of similarities between samples, for example calculated from
    microbiome distances. Kernel methods such as MiRKAT use it to relate a
    variable to the whole microbiome.

Late fusion
:   Training a separate model on each data type and combining their
    predictions, for example as a weighted average. Also called stacking.

Mantel test
:   A test of whether two dissimilarity matrices, for example from two data
    types, are correlated.

Mass spectrometry (LC-MS, GC-MS)
:   A method that separates and weighs molecules. In metabolomics, it
    measures hundreds or thousands of metabolite features.

Metabolite
:   A small molecule produced or used in metabolism, such as a short-chain
    fatty acid.

MiRKAT
:   The microbiome regression-based kernel association test. Tests whether a
    single variable, such as a metabolite, is related to the microbiome as a
    whole.

MOFA (multi-omics factor analysis)
:   An unsupervised integration method that finds factors shared by several
    data types, like a PCA for several tables.

Multi-omics
:   Analysing several data types measured from the same samples, such as
    taxa, metabolites and host biomarkers.

MultiAssayExperiment (MAE)
:   The Bioconductor container for multi-omics data. It holds one experiment
    per data type, a shared table of sample information (`colData`) and a
    sample map that links them.

NMR spectroscopy
:   Nuclear magnetic resonance spectroscopy. In metabolomics, it measures
    the concentrations of tens of known metabolites.

Pairwise association
:   The correlation between every feature of one data type and every feature
    of another.

Permutation test
:   A test that shuffles the samples many times to see how often a result as
    large as the real one arises by chance.

Prebiotic
:   A fibre that the host cannot digest but some gut bacteria can, such as
    xylo-oligosaccharides (XOS).

primary / colname
:   The columns of `sampleMap()`: `primary` is the subject (a row of
    `colData(mae)`), and `colname` is the name of the column in an
    experiment.

rclr (robust centred log-ratio) transformation
:   A log-ratio transformation that ignores the zeros. By default mia fills
    them in by matrix completion; with `impute = FALSE` they stay missing.

Robust (in joint RPCA)
:   The method fits only the observed values and fills in zeros or missing
    values by matrix completion, so that a few extreme values do not
    dominate.

Sample map
:   The table in a `MultiAssayExperiment` that links each column of each
    experiment to a subject in the sample information.

Short-chain fatty acids (SCFA)
:   Metabolites such as butyrate, acetate and propionate that gut bacteria
    produce when they ferment fibre.

Shotgun metagenomics
:   Sequencing all DNA in a sample. Gives counts of taxa and of genes or
    pathways (functions).

Sparse partial least squares (sPLS)
:   A method that finds directions in two or more tables whose sample scores
    are strongly correlated (and, in DIABLO, related to the outcome), using
    only a few features.

Stacking
:   Late fusion in which a second model (the meta-learner) learns weights to
    combine the cross-validated predictions of the models of each data type.

Standardization
:   Scaling each feature to mean 0 and standard deviation 1, so that features
    measured on different scales can be compared.

Supervised method
:   A method that is given an outcome, such as the diet group, and looks for
    the features that predict it.

Targeted / untargeted metabolomics
:   Measuring a chosen list of known metabolites / everything that can be
    detected, including unknown molecules.

TreeSummarizedExperiment (TreeSE)
:   The Bioconductor container for microbiome data: abundance tables, taxonomy,
    sample information and a phylogenetic tree.

Unsupervised method
:   A method that looks for patterns in the data without knowing the study
    groups, such as PCA or MOFA.

UpSet plot
:   A plot that shows the sizes of the overlaps between several sets, here
    which samples have which data types.

Weight (MOFA)
:   How strongly a feature drives a factor. Like a loading in PCA.
