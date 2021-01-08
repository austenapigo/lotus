context("quad.rarefied")
library(testthat)
library(lotus)

# Test whether the output is a data frame
test_that("structural.specificity() returns a data frame", {
  output_table <- structural.specificity(quad.rarefied, abundance.weighted = FALSE, trim = TRUE)
  expect_is(output_table, "data.frame")
})

# devtools::test()
# devtools::test_coverage()

# structural.object <- structural.specificity(quad.rarefied, abundance.weighted = TRUE, trim = TRUE)
# mean(structural.object$Structural.Specificity)
# matches source code: -0.7063895
###

## Internal Check ##
# mean(structural.dev[[1]]$Mean.Deviance)
# -0.7333748 

# network.object <- network.specificity(quad.rarefied, abundance.weighted = TRUE, trim = TRUE)
# mean(network.object$Network.Specificity)
# network.object <- network.specificity(quad.rarefied, abundance.weighted = FALSE, trim = TRUE)
# mean(network.object$Network.Specificity)
# network.object <- network.specificity(example.dat, abundance.weighted = TRUE, trim = TRUE)

# null.network.object <- null.network(quad.rarefied, iterations = 100, abundance.weighted = TRUE, randomization.method = "shuffle.web", trim = TRUE, notify = TRUE)
# null.network.object <- null.network(quad.rarefied, iterations = 100, abundance.weighted = FALSE, randomization.method = "shuffle.web", trim = TRUE, notify = TRUE)
# null.network(example.dat, iterations = 100, abundance.weighted = FALSE, randomization.method = "shuffle.web", trim = TRUE, notify = TRUE)

# network.dev <- deviance.network(quad.rarefied, randomized = null.network.object, model = "second", abundance.weighted = TRUE, trim = TRUE, notify = TRUE)
# head(network.dev[[1]])
# network.dev[[2]]
# network.dev[[81]]
# mean(network.dev[[1]]$Mean.Deviance)

# Internal Check
# phylogenetic.object <- phylogenetic.specificity(quad.rarefied, utree, abundance.weighted = TRUE, trim = TRUE)
# phylogenetic.object
# mean(phylogenetic.object$Phylogenetic.Specificity)
# -146.6859 

# phylogenetic.dev <- deviance.phylogenetic(quad.rarefied, utree, null.model = "taxa.labels", iterations = 100, model = "second", abundance.weighted = TRUE, trim = TRUE, notify = TRUE)
# head(phylogenetic.dev[[1]])
# phylogenetic.dev[[2]]
# mean(phylogenetic.dev[[1]]$Mean.Deviance)
# # -0.9141502 vs. -0.9135079

# beta.object <- beta.specificity(quad.rarefied, index = "morisita.horn", trim = TRUE, notify = TRUE)
# beta.object <- beta.specificity(example.dat, index = "morisita.horn", trim = TRUE, notify = TRUE)
# beta.object
# Symbiont A: 0.6666667; Symbiont B: 0.8902854; Symbiont C: 0.2496302; Symbiont D: 1.0000000
# mean(beta.object$Similarity.Index)
# mean morisita-horn 0.2008582 vs. 0.2008582

# null.beta.object <- null.beta(quad.rarefied, index = "morisita.horn", randomization.method = "shuffle.web", iterations = 100, trim = TRUE, notify = TRUE)
# null.beta.object <- null.beta(example.dat, index = "morisita.horn", randomization.method = "shuffle.web", iterations = 100, trim = TRUE, notify = TRUE)
# mean(null.beta.object$Similarity.Index)
# match 0.1069268 vs. 0.1066938



# Check to see if TRIM works 
# dim(subset(null.structural.object, null.structural.object$Randomization == 1))
# dim(subset(null.structural.object, null.structural.object$Randomization == 2))

# ## Internal check ##
# null.structural.object.ab <- null.structural(quad.rarefied, iterations = 100, abundance.weighted = TRUE, randomization.method = "shuffle.web", trim = TRUE, notify = TRUE)
# mean(null.structural.object.ab$Structural.Specificity)
# # mean shannon's h: -0.6223652
# null.structural.object.pa <- null.structural(quad.rarefied, iterations = 100, abundance.weighted = FALSE, randomization.method = "shuffle.web", trim = TRUE, notify = TRUE)
# mean(null.structural.object.pa$Structural.Specificity)
# # mean host richness: -4.939931