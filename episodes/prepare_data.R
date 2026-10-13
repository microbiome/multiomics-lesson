# Helper used when the lesson is built.
#
# Each episode sources this file in a hidden chunk so that it can be built on
# its own. It recreates the `mae` object that learners prepare in episode 2
# and keep in their R session. Learners do not need this file.

suppressPackageStartupMessages({
    library(mia)
    library(MultiAssayExperiment)
})

data("HintikkaXOData", package = "mia")
mae <- HintikkaXOData
mae$XOS <- factor(mae$XOS, levels = c(0, 1), labels = c("control", "XOS"))

tse <- mae[["microbiota"]]
tse <- agglomerateByRank(tse, rank = "Genus")
tse <- transformAssay(tse, assay.type = "counts", method = "relabundance")
tse <- subsetByPrevalent(tse, assay.type = "relabundance",
                         detection = 0.001, prevalence = 0.2)
tse <- transformAssay(tse, assay.type = "counts", method = "rclr")
mae[["microbiota"]] <- tse

met <- mae[["metabolites"]]
assay(met, "nmr")[is.na(assay(met, "nmr"))] <- 0
met <- transformAssay(met, assay.type = "nmr", method = "log10",
                      pseudocount = TRUE)
mae[["metabolites"]] <- met
