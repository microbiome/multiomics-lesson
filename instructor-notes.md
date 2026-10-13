---
title: 'Instructor Notes'
---

## General

- The lesson takes about 8 hours of teaching and exercises (see the front
  matter of each episode). It assumes that learners know
  `TreeSummarizedExperiment`, agglomeration, transformation, ordination and
  the basics of machine learning, for example from the
  [microbiome lesson](https://microbiome.github.io/microbiome-lesson/).
  Not every learner has taken that lesson, so do not refer to its examples
  or dataset; the recap callout in episode 1 lists what is assumed.
- Episode 1 has no code. Use it to introduce multi-omics, why the
  microbiome is studied with other data types, and the three groups of
  methods (exploratory, unsupervised, supervised), which structure the rest
  of the lesson. The section on how the data is generated is background:
  spend about 10 minutes on it, and stress the three take-homes: one
  transformation per data type, matched samples, and batches mixed up with
  the groups.
- Episode 2 introduces the study (HintikkaXO) and the
  `MultiAssayExperiment`. `primary` and `colname` of the `sampleMap()` are
  identical in HintikkaXO, so use the rat C1 example to show the concept.
- The summary page (episode 6) wraps up the whole lesson and how to choose
  a method. Leave about 10 minutes for it.
- The episodes build on each other through the `mae` object prepared in
  episode 2 and kept in the R session; nothing is saved to disk. If a
  learner restarts R or falls behind, they can catch up by running the code
  of episode 2 again (a few seconds).
- The setup page installs the packages with explicit commands instead of a
  script, so that learners see where each package comes from (Bioconductor,
  CRAN or GitHub). Go through it briefly at the start if needed.
- The lesson uses `plotJointRPCA()` from miaViz, which is in the
  Bioconductor release from 3.24 on. Learners need the current release.
- Live-code every chunk. The chunks are kept short on purpose, and the
  lesson uses only package functions, no custom helper functions.

## Before the workshop

- **MOFA2 needs Python.** The first `run_mofa()` installs a Python
  environment through basilisk, which can take 10 minutes or more. The setup
  page asks learners to train a test model before the workshop. Check this
  at the start of the workshop, at the latest before episode 4.
- multiomics and IntegratedLearner are installed from GitHub; check that
  they installed.

## The story of the lesson

All episodes use the HintikkaXO data: 40 rats on a high-fat or low-fat diet,
with or without the prebiotic XOS, with microbiota, metabolites and host
biomarkers. Two threads connect the episodes.

**The fat thread**

- Episode 2 (challenge 2): fat has a large effect on the microbiota (PCA)
  and on butyrate, acetate and propionate.
- Episode 3: the Mantel test is significant over all rats but weaker within
  the high-fat rats. Use challenge 1 to discuss how a shared factor (the
  diet) creates associations between data types. The butyrate challenge
  makes the same point for a single pair of features.
- Episode 4: joint RPCA roughly separates the fat groups, and MOFA factor 1,
  shared by all three data types, separates them perfectly.
- Episode 5: every data type predicts fat with AUC near 1.

**The XOS thread**

- XOS was the intervention of the study, but its effect is small. MOFA finds
  no XOS factor (episode 4, challenge 3), IntegratedLearner predicts it with
  cross-validated AUCs of only 0.5–0.6, and the DIABLO cross-validated error is about 50 %
  (episode 5). Use this to discuss that integration does not create a
  signal, and that cross-validation is essential.

**Consistent features**

- Glucose, butyrate and nicotinate come up in MiRKAT, MOFA factor 1 and
  DIABLO. Episode 5, challenge 2, asks learners to notice this.

**Episode-specific notes**

- Episode 3: within the diet groups, the butyrate association is very
  different (about 0 in the high-fat rats, about 0.4 in the low-fat rats).
  Use it to show that associations over all samples can hide or create
  links.
- Episode 4: `run_mofa()` warns that factor 3 follows the per-sample average
  of the metabolites; discuss it as a possible technical factor. The
  left/right direction of the joint RPCA map is computed inline, because the
  sign of an ordination axis is arbitrary.
- Episode 5: IntegratedLearner's `seed` argument only fixes the folds;
  `set.seed()` is needed for reproducible random forests. When every data
  type predicts the fat perfectly, the stacking weights change strongly with
  the seed. In challenge 3, the first DIABLO component of the XOS model
  separates the fat groups.

## Notes on the data

- The data comes with mia (`data("HintikkaXOData")`), so nothing is
  downloaded.
- The taxonomy has SILVA prefixes (`D_5__` for genus). The lesson keeps them
  to avoid custom code; mention that they do not change the analysis.
- Many genera are `uncultured`; mia makes their names unique with a number
  (`D_5__uncultured_2`).
- The metabolites have 27 missing values (below the detection limit), which
  episode 2 replaces with 0. The biomarkers have 95 missing values, which
  are kept: MOFA handles them, and IntegratedLearner drops those biomarkers
  (`na.rm = TRUE`).
- Joint RPCA and MOFA use rclr values with `impute = FALSE`, so that the
  methods fill in the zeros themselves. The default rclr of mia (and vegan)
  returns a low-rank reconstruction, which would make the microbiota look
  perfectly explained in MOFA.

## Common problems

- **Harmless warnings.**
  - Many multiomics functions and `getWithColData()` warn
    `'experiments' dropped` and print `harmonizing input`.
  - `upsetSamples()` prints deprecation warnings from UpSetR.
  - `correlate_factors_with_covariates()` warns that it converts the
    covariates to numbers.
  - MOFA warns about correlated factors; with 40 samples this is common.
  - `block.splsda()` prints that the design matrix has changed.
- **IntegratedLearner cannot find the outcome.** It reads the outcome and
  the sample IDs from the `colData` of each experiment. Run the three
  `getWithColData()` lines of episode 5 first.
- **Known package issues.** `getMantel(group.by = )` fails on a
  `MultiAssayExperiment` (multiomics), so the lesson subsets the samples
  instead. The IntegratedLearner ROC plot leaves an empty lower half when
  there is no validation set.
- **DIABLO `perf()` and the seed.** Use the `seed` argument of `perf()`,
  not `set.seed()`: with some versions of mixOmics and BiocParallel, a
  second `perf()` call after `set.seed()` fails with
  `invalid value of 'seed'`.
- **IntegratedLearner prints "training data".** The AUCs in its summary are
  cross-validated, despite the label.
