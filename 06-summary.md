---
title: "Summary"
teaching: 10
exercises: 0
---

:::::::::::::::::::::::::::::::::::::: questions

- What did we learn about our study?
- How do I choose an integration method?
- Where can I learn more?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Summarize the integration workflow of this lesson.
- Choose a group of integration methods for a research question.
- Find resources for further learning.

::::::::::::::::::::::::::::::::::::::::::::::::

## What we did

We integrated three data types of 40 rats, step by step:

| Episode | What we learned about the study |
|---|---|
| [MultiAssayExperiment](02-multiassayexperiment.Rmd) | The three data types fit into one object. Fat in the diet has a large effect on the microbiota and the metabolites; XOS a small one. |
| [Exploratory](03-exploratory.Rmd) | The microbiota and the metabolites are related, but weakly. Glucose, butyrate and nicotinate are most strongly linked to the microbiota. Much of the association comes from the diet. |
| [Unsupervised](04-unsupervised.Rmd) | The main pattern that all three data types share is the fat in the diet (joint RPCA, MOFA factor 1). No pattern matches XOS. |
| [Supervised](05-supervised.Rmd) | Every data type predicts the fat in the diet almost perfectly, but XOS is hard to predict. DIABLO selects a few genera and metabolites that tell the diets apart. |

The same features came up again and again with different methods. Results
that are consistent across methods are more convincing.

## Which method should I use?

To choose a method for your own study, start from the research question:

- *Are my data types related at all, and which features go together?*
  **Exploratory** methods: a Mantel test, then MiRKAT and pairwise
  associations for the details.
- *What are the main patterns in my data?* **Unsupervised** methods, such
  as joint RPCA and MOFA.
- *Can I predict an outcome, and which features predict it?*
  **Supervised** methods, such as IntegratedLearner and DIABLO. Always
  evaluate them with cross-validation or an independent dataset.

In every case, prepare each data type carefully, look at each of them alone
first, and remember that a factor shared by the data types, such as the
diet, creates associations between them.

## Where to go next

To continue learning, see the OMA book chapters on
[multi-omics](https://microbiome.github.io/OMA/docs/devel/pages/multiomics.html),
[cross-association](https://microbiome.github.io/OMA/docs/devel/pages/cross_correlation.html),
[multi-assay ordination](https://microbiome.github.io/OMA/docs/devel/pages/multiassay_ordination.html)
and
[machine learning](https://microbiome.github.io/OMA/docs/devel/pages/machine_learning.html),
the [ECCB 2026 hologenomics tutorial](https://microbiome.github.io/ECCB2026OHA/articles/OHA.html),
which applies the same methods to taxa and pathways of a colorectal cancer
study, and the [training slides](https://microbiome.github.io/outreach/) of
the OMA framework.

::::::::::::::::::::::::::::::::::::: keypoints

- Exploratory methods test whether data types are related, unsupervised methods find shared patterns, and supervised methods predict an outcome.
- Choose the group of methods by the research question, and look for results that are consistent across methods.
- Integrating data types does not create a signal that is not in the data.

::::::::::::::::::::::::::::::::::::::::::::::::
