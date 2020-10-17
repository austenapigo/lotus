#############################
# Read in pacakges and data #
#############################
# library(tidyverse)
# library(ggpmisc)
# library(picante)
# library(vegan)
# library(bipartite)
# 
# # Read in data 
# dat <- read.csv("/Users/austenapigo/Desktop/github/lotus/otherdat/example_dat.csv", row.names = 1, header = TRUE)
# quad.rarefied <- read.csv("quad_rarefied.csv", row.names = 1, header = TRUE)
# quad.rarefied <- read.csv("/Users/austenapigo/Desktop/github/lotus/data/quad_rarefied.csv", row.names = 1, header = TRUE)
# save(quad.rarefied, file = "quad_rarefied.rda")
# utree <- read.tree("/Users/austenapigo/Desktop/github/lotus/otherdat/utree.txt")
# save(utree, file = "utree.rda")

##############################################################
# `structural.specificity`: calculate structural specificity #
##############################################################
#' Structural Specificity
#' 
#' Calculate structural specificity not corrected by null mdoels. 
#'
#' @param x Data frame of hosts populating rows and symbionts populating columns. 
#' 
#' @param abundance.weighted Logical. TRUE calculates Shannon's H per symbiont. FALSE calculates host richness per symbiont. 
#' 
#' @param trim Logical. TRUE removes symbionts that occupy one host sample. FALSE keeps all symbionts. 
#'
#' @return A data frame with symbiont identifiers and structural specificity values. 
#' @export
#' @examples
#' # Calculate host richness per symbiont
#' structural.specificity(quad.rarefied, abundance.weighted = TRUE, trim = TRUE)
#' 
#' # Calculate Shannon's H per symbiont
#' structural.specificity(quad.rarefied, abundance.weighted = TRUE, trim = TRUE)
structural.specificity <- function(x, abundance.weighted = TRUE, trim = TRUE) {
  # Calculate host richness or Shannon's H
  ifelse(abundance.weighted == TRUE, structural <- -1 * vegan::diversity(t(x)), structural <- -1 * vegan::specnumber(t(x)))
  # Make data frame
  structural.dat <- data.frame(Structural.Specificity = structural)
  # Trim noise
  ifelse(trim == TRUE, 
         # if calculating Shannon's H, 
         ifelse(abundance.weighted == TRUE, 
                # only consider symbionts with a Shannon's H less than 0
                structural.dat <- subset(structural.dat, Structural.Specificity < 0), 
                # otherwise, only consider symbionts with a host richness less than -1 
                structural.dat <- subset(structural.dat, Structural.Specificity < -1)), 
         structural.dat <- structural.dat) 
  # Return object
  return(structural.dat)
}

# structural.object <- structural.specificity(quad.rarefied, abundance.weighted = TRUE, trim = TRUE)
# mean(structural.object$Structural.Specificity)
# mean Shannon's H is -0.7925232 matches source code vs. -0.7925232
# structural.object <- structural.specificity(quad.rarefied, abundance.weighted = FALSE, trim = TRUE)
# mean(structural.object$Structural.Specificity)
# mean host richness is -5.842963 vs. -5.842963

#######################################################################
# `null.structural`: calculate null models for structural specificity #
#######################################################################
#' Structural Specificity Null Models 
#' 
#' Generate null models and calculate structural specificity per symbiont within each community randomization. 
#'
#' @param x Data frame of hosts populating rows and symbionts populating columns. 
#' 
#' @param iterations Integer. Indicate the number of randomized communities to generate. 
#' 
#' @param abundance.weighted Logical. TRUE calculates Shannon's H per symbiont. FALSE calculates host richness per symbiont. 
#' 
#' @param randomization.method Randomization method. Usage borrowed directly from bipartite::nullmodel. 
#' Specify as "r2dtable", "swap.web", "vaznull", "shuffle.web" or "mgen". 
#' 
#' @param trim Logical. TRUE removes symbionts that occupy one host sample. FALSE keeps all symbionts. 
#' 
#' @param notify Logical. TRUE prints the current iteration of the for loop. 
#'
#' @return A data frame with columns that refer to symbiont identifiers, absolute read abundance, structural specificities and randomization
#' identifiers. 
#' @export
#' @examples
#' # Generate randomized communities and calculate structural specificity per symbiont 
#' null.structural(quad.rarefied, iterations = 100, abundance.weighted = TRUE, randomization.method = "shuffle.web", trim = TRUE, notify = TRUE)
null.structural <- function(x, iterations = 100, abundance.weighted = TRUE, randomization.method = c("r2dtable", "swap.web", "vaznull", "shuffle.web", "mgen"), trim = TRUE, notify = TRUE) {
  # Set seed
  set.seed(123)
  # Match argument specified
  randomization.method <- match.arg(randomization.method)
  # Generate 100 randomized communities
  null.structural <- bipartite::nullmodel(x, N = iterations, method = randomization.method)
  # Make holding list
  null.dats <- list()
  # Make holding vectors 
  Symbiont <- rep()
  Abundance <- rep()
  # Calculate structural specificity for null models
  for (i in 1:length(null.structural)) {
    # Call a randomized community
    null <- null.structural[[i]]
    # Add row and column names
    rownames(null) <- rownames(x)
    colnames(null) <- colnames(x)
    # Make as a data frame
    null <- as.data.frame(null)
    # Subset symbiont name and calculate abundance per symbiont
    for (j in 1:ncol(null)) {
      # Pull symbiont name
      Symbiont[j] <- colnames(null)[j]
      # Calculate total read abundance per symbiont
      Abundance[j] <- sum(null[,j])
    }
    # Calculate structural specificity
    ifelse(abundance.weighted == TRUE, 
           Structural.Specificity <- -1 * vegan::diversity(t(null)), 
           Structural.Specificity <- -1 * vegan::specnumber(t(null)))
    # Make data frame
    null.temp <- data.frame(Symbiont, Abundance, Structural.Specificity)
    # Populate into holding list 
    null.dats[[i]] <- null.temp
    # Print iteration
    ifelse(notify == TRUE, print(i), NaN)
  }
  # Make into one data frame
  null.dats <- as.data.frame(do.call("rbind", null.dats))
  # Add randomization number as a new column 
  null.dats$Randomization <- as.factor(rep(1:iterations, each = ncol(x)))
  rownames(null.dats) <- NULL
  # Trim noise
  ifelse(trim == TRUE, 
         # If calculating Shannon's H, 
         ifelse(abundance.weighted == TRUE, 
                # only consider symbionts with a Shannon's H greater than 0
                null.dats <- subset(null.dats, null.dats$Structural.Specificity < 0), 
                # otherwise, only consider symbionts with a host richness greater than 1 
                null.dats <- subset(null.dats, null.dats$Structural.Specificity < -1)), 
         null.dats <- null.dats) 
  # Read out data frame 
  return(data.frame(null.dats))
}

# null.structural.object <- null.structural(quad.rarefied, iterations = 100, abundance.weighted = TRUE, randomization.method = "shuffle.web", trim = TRUE, notify = TRUE)
# mean(null.structural.object$Structural.Specificity)
# mean shannon's h: -0.6402127 vs. -0.6407739
# null.structural.object <- null.structural(quad.rarefied, iterations = 100, abundance.weighted = FALSE, randomization.method = "shuffle.web", trim = TRUE, notify = TRUE)
# mean(null.structural.object$Structural.Specificity)
# mean host richness: -4.086787 vs. -3.926589
# head(null.structural.object)

###########################################################################
# `deviance.structural`: calculate the deviance in structural specificity #
###########################################################################
#' Deviance in Structural Specificity to Null Models 
#' 
#' Calculate the deviance in observed structural specificity to a null model of structural specificity per symbiont. 
#' Deviance calculations are measured per symbiont and averaged per host sample. For example, all symbionts within a given host 
#' are evaluated for their host specificity across the entire host community. The host specificities of each symbiont are averaged
#' to calculate the mean host specificity for symbionts within a given host. 
#'
#' @param x Data frame of hosts populating rows and symbionts populating columns. 
#' 
#' @param randomized Data frame. Output from null.structural function. 
#' 
#' @param abundance.weighted Logical. TRUE calculates Shannon's H per symbiont. FALSE calculates host richness per symbiont. 
#' 
#' @param trim Logical. TRUE removes symbionts that occupy one host sample. FALSE keeps all symbionts. 
#' 
#' @param notify Logical. TRUE prints the current iteration of the for loop. 
#'
#' @return A list. First element in the list is a data frame with columns that refer to the host sample identifiers, mean deviance in 
#' structural specificity, standard error of structural specificities, number of symbionts per host sample and average symbiont read abundance. 
#' All subsequent elements of the list are plots of uncorrected host specificity as a function of log symbiont read abundance with the null mode
#' in blue relative to the symbionts host specificities in red. 
#' 
#' @export
#' @examples
#' # Calculate mean deviance per symbiont per host sample and visualize null vs. observed host specifities 
#' structural.dev <- deviance.structural(quad.rarefied, randomized = null.structural.object, abundance.weighted = TRUE, trim = TRUE, notify = TRUE)
#' structural.dev[[1]] # View data frame of output 
#' structural.dev[[2]] # View first graph 
#' structural.dev[[81]] # View last graph 
deviance.structural <- function(x, randomized = null.structural.object, abundance.weighted = TRUE, trim = TRUE, notify = TRUE) {
  # Make holding vectors 
  structural.plots <- list()
  mean.structural <- rep()
  se.structural <- rep()
  host.sample <- rep()
  num.symbionts <- rep()
  read.abund <- rep()
  # For every host sample
  for (i in 1:nrow(x)) {
    # Subset a host
    x.sub <- x[i, 1:dim(x)[2], drop = FALSE]
    # Remove symbionts with abundance of zero
    x.sub <- x.sub[ , colSums(x.sub) > 0, drop = FALSE]
    # Save column names 
    x.names <- colnames(x.sub)
    # Filter entire community
    x.input <- x[, colnames(x) %in% x.names, drop = FALSE]
    # Remove rows and columns that sum to zero
    x.input <- as.data.frame(x.input[rowSums(x.input) > 0, colSums(x.input) > 0])
    # Calculate structural specificity
    ifelse(abundance.weighted == TRUE, 
           Structural.Specificity <- -1 * diversity(t(x.input)), 
           Structural.Specificity <- -1 * specnumber(t(x.input)))
    # Make holding vectors
    Symbiont <- rep()
    Abundance <- rep()
    # For every symbiont
    for (j in 1:ncol(x.input)) {
      # Pull symbiont name
      Symbiont[j] <- colnames(x.input)[j]
      # Calculate total read abundance per symbiont
      Abundance[j] <- sum(x.input[,j])
    }
    # Make data frame
    structural.dat <- data.frame(Structural.Specificity, Symbiont, Abundance)
    # Trim noise
    ifelse(trim == TRUE, 
           # if calculating Shannon's H, 
           ifelse(abundance.weighted == TRUE, 
                  # only consider symbionts with a Shannon's H less than 0
                  structural.dat <- subset(structural.dat, structural.dat$Structural.Specificity < 0), 
                  # otherwise, only consider symbionts with a host richness less than -1 
                  structural.dat <- subset(structural.dat, structural.dat$Structural.Specificity < -1)), 
           structural.dat <- structural.dat) 
    # Plot null vs. empirical per sample
    structural.plots[[i+1]] <- 
      ggplot2::ggplot(structural.dat, aes(y = Structural.Specificity, x = log(Abundance))) +
      geom_point(data = randomized, aes(y = Structural.Specificity, x = log(Abundance)), color = "grey", alpha = 0.5, show.legend = TRUE, size = 3) +
      geom_smooth(data = randomized, aes(y = Structural.Specificity, x = log(Abundance)), color = "black", method = "lm", se = FALSE, lwd = 1, lty = "dashed", show.legend = FALSE, formula = y ~ x + I(x^2)) + 
      geom_point(color = "red", alpha = 1, show.legend = TRUE, size = 3) +
      stat_poly_eq(data = randomized, parse = TRUE, aes(label = ..eq.label..), formula = y ~ x + I(x^2), label.x = "left", label.y = "bottom", color = "black", size = 5) + 
      theme_bw() +
      ggtitle(rownames(x)[i]) + 
      theme_bw() +
      theme(
        axis.text.x = element_text(size = 13, color = "black"), 
        axis.text.y = element_text(size = 13, color = "black"),
        axis.title.x = element_text(size = 13, margin = margin(t = 5, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 13, margin = margin(t = 0, r = 5, b = 0, l = 0)),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 13, face = "bold.italic"),
        aspect.ratio = 0.85
      ) +
      labs(y = "Uncorrected Structural Specificity", x = "Log Absolute Symbiont Read Abundance")
    # Get model coefficients for null model
    null.eqn <- summary(lm(Structural.Specificity ~ log(Abundance) + I(log(Abundance)^2), data = randomized))
    null.eqn$coefficients[1, 1]
    null.eqn$coefficients[2, 1]
    null.eqn$coefficients[3, 1]
    # Calculate mean deviance
    mean.structural[i] <- mean(structural.dat$Structural.Specificity - (null.eqn$coefficients[1, 1] + null.eqn$coefficients[2, 1]*log(structural.dat$Abundance) + null.eqn$coefficients[3, 1]*log(structural.dat$Abundance)^2))
    # Calculate standard error of mean deviance
    se.structural[i] <- sd(structural.dat$Structural.Specificity - (null.eqn$coefficients[1, 1] + null.eqn$coefficients[2, 1]*log(structural.dat$Abundance) + null.eqn$coefficients[3, 1]*log(structural.dat$Abundance)^2)) / sqrt(length(structural.dat$Structural.Specificity - (null.eqn$coefficients[1, 1] + null.eqn$coefficients[2, 1]*log(structural.dat$Abundance) + null.eqn$coefficients[3, 1]*log(structural.dat$Abundance)^2)))
    # Host sample name
    host.sample[i] <- rownames(x)[i]
    # Number of symbionts column
    num.symbionts[i] <- length(structural.dat$Structural.Specificity)
    # Symbiont read abundance 
    read.abund[i] <- mean(structural.dat$Abundance)
    # Check current iteration
    ifelse(notify == TRUE, print(i), NaN)
    # Populate deviance data into data.frame
    structural.plots[[1]] <- data.frame(Host.Sample = host.sample, 
                                        Mean.Deviance = mean.structural, 
                                        Mean.Deviance.SE = se.structural,
                                        Number.of.Symbionts = num.symbionts,
                                        Avg.Symbiont.Abundance = read.abund)
  }
    return(structural.plots)
}

# structural.dev <- deviance.structural(quad.rarefied, randomized = null.structural.object, abundance.weighted = TRUE, trim = TRUE, notify = TRUE)
# head(structural.dev[[1]])
# structural.dev[[2]]
# structural.dev[[81]]
# mean(structural.dev[[1]]$Mean.Deviance)
# # -0.9141502 vs. -0.9135079

##################################################################
# `phylogenetic.specificity`: calculate phylogenetic specificity #
##################################################################
#' Phylogenetic Specificity
#' 
#' Calculate phylogenetic specificity not corrected by null mdoels. 
#'
#' @param x Data frame of hosts populating rows and symbionts populating columns. 
#' 
#' @param utree Newick formatted phylogenetic tree. 
#' 
#' @param abundance.weighted Logical. TRUE calculates abundance-weighted mean pairwise phylogenetic distance. 
#' FALSE calculates presence-absence mean pairwise phylogenetic distance. 
#' 
#' @param trim Logical. TRUE removes symbionts that occupy one host sample. FALSE keeps all symbionts. 
#'
#' @return A data frame with symbiont identifiers and phylogenetic specificity values. 
#' @export
#' @examples
#' # Calculate mean pairwise phylogenetic distance per symbiont
#' phylogenetic.specificity(quad.rarefied, utree, abundance.weighted = TRUE, trim = TRUE)
phylogenetic.specificity <- function(x, utree, abundance.weighted = TRUE, trim = TRUE) {
  # Calculate abundance-weighted or presence-absence mean pairwise phylogenetic distance
  ifelse(abundance.weighted == TRUE, 
         phylogenetic <- -1 * mpd(t(x), dis = cophenetic(utree), abundance.weighted = TRUE), 
         phylogenetic <- -1 * mpd(t(x), dis = cophenetic(utree), abundance.weighted = FALSE))
  # Make data frame
  phylogenetic.dat <- data.frame(Phylogenetic.Specificity = phylogenetic)
  rownames(phylogenetic.dat) <- colnames(x)
  # Trim noise
  ifelse(trim == TRUE, 
    # only consider symbionts with a mpd less than 0
    phylogenetic.dat <- subset(phylogenetic.dat, Phylogenetic.Specificity < 0), 
    phylogenetic.dat <- phylogenetic.dat) 
  # Return object
  return(phylogenetic.dat)
}

# phylogenetic.object <- phylogenetic.specificity(quad.rarefied, utree, abundance.weighted = TRUE, trim = TRUE)
# phylogenetic.object
# mean(phylogenetic.object$Phylogenetic.Specificity)
# -144.7806 vs. -144.7806

###############################################################################
# `deviance.phylogenetic`: calculate the deviance in phylogenetic specificity #
###############################################################################
#' Deviance in Phylogenetic Specificity to Null Models 
#' 
#' Calculate the deviance in observed phylogenetic specificity to a null model of structural specificity per symbiont. 
#' Deviance calculations are measured per symbiont and averaged per host sample. For example, all symbionts within a given host 
#' are evaluated for their host specificity across the entire host community. The host specificities of each symbiont are averaged
#' to calculate the mean host specificity for symbionts within a given host. 
#'
#' @param x Data frame of hosts populating rows and symbionts populating columns. 
#' 
#' @param utree Newick formatted phylogenetic tree. 
#' 
#' @param null.model Randomization method. Usage borrowed directly from picante::mpd 
#' Specify as "taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", "independentswap" or "trialswap". 
#' 
#' @param abundance.weighted Logical. TRUE calculates abundance-weighted mean pairwise phylogenetic distance. 
#' FALSE calculates presence-absence mean pairwise phylogenetic distance. 
#'  
#' @param trim Logical. TRUE removes symbionts that occupy one host sample. FALSE keeps all symbionts. 
#' 
#' @param notify Logical. TRUE prints the current iteration of the for loop.
#'
#' @return A data frame with columns that refer to the host sample identifiers, mean deviance in 
#' structural specificity, standard error of structural specificities, number of symbionts per host sample and average symbiont read abundance. 
#' 
#' @export
#' @examples
#' # Calculate mean deviance per symbiont per host sample and visualize null vs. observed host specifities 
#' deviance.phylogenetic(quad.rarefied, utree, null.model = "taxa.labels", iterations = 100, abundance.weighted = TRUE, trim = TRUE, notify = TRUE)
deviance.phylogenetic <- function(x, utree, null.model = c("taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", "independentswap", "trialswap"), iterations = 100, abundance.weighted = TRUE, trim = TRUE, notify = TRUE) {
  # Make holding vectors 
  phylogenetic.plots <- list()
  mean.phylogenetic <- rep()
  se.phylogenetic <- rep()
  num.symbionts <- rep()
  read.abund <- rep()
  # For every host sample
  for (i in 1:nrow(x)) {
    # Subset a host
    x.sub <- x[i, 1:dim(x)[2], drop = FALSE]
    # Remove symbionts with abundance of zero
    x.sub <- x.sub[ , colSums(x.sub) > 0, drop = FALSE]
    # Save column names 
    x.names <- colnames(x.sub)
    # Filter entire community
    x.input <- x[, colnames(x) %in% x.names, drop = FALSE]
    # Remove rows and columns that sum to zero
    x.input <- as.data.frame(x.input[rowSums(x.input) > 0, colSums(x.input) > 0])
    # Calculate phylogenetic specificity
    null.model <- match.arg(null.model)
    # Set a seed
    set.seed(123)
    # Calculate phylogenetic specificity
    ifelse(abundance.weighted == TRUE, 
           Phylogenetic.Specificity <- ses.mpd(t(x.input), dis = cophenetic(utree), null.model = null.model, abundance.weighted = TRUE, runs = iterations), 
           Phylogenetic.Specificity <- ses.mpd(t(x.input), dis = cophenetic(utree), null.model = null.model, abundance.weighted = FALSE, runs = iterations))
    # Make holding vectors
    Symbiont <- rep()
    Abundance <- rep()
    # For every symbiont
    for (j in 1:ncol(x.input)) {
      # Pull symbiont name
      Symbiont[j] <- colnames(x.input)[j]
      # Calculate total read abundance per ASV
      Abundance[j] <- sum(x.input[,j])
    }
    # Make data frame
    phylogenetic.dat <- data.frame(Symbiont, Abundance, Richness = Phylogenetic.Specificity$ntaxa, Phylogenetic.Specificity = -1* Phylogenetic.Specificity$mpd.obs.z)
    # Convert NA to 0
    phylogenetic.dat[is.na(phylogenetic.dat)] <- 0
    # Trim noise
    ifelse(trim == TRUE, 
           # otherwise, only consider symbionts with a host richness greater than 1 
           phylogenetic.dat <- subset(phylogenetic.dat, phylogenetic.dat$Richness > 1), 
           phylogenetic.dat <- phylogenetic.dat) 
    # Calculate mean deviance
    mean.phylogenetic[i] <- mean(phylogenetic.dat$Phylogenetic.Specificity)
    # Calculate the standard error of the mean
    se.phylogenetic[i] <- sd(phylogenetic.dat$Phylogenetic.Specificity) / sqrt(length(phylogenetic.dat$Phylogenetic.Specificity))
    # Number of symbionts column
    num.symbionts[i] <- length(phylogenetic.dat$Phylogenetic.Specificity)
    # Symbiont read abundance 
    read.abund[i] <- mean(phylogenetic.dat$Abundance)
    # Check current iteration
    ifelse(notify == TRUE, print(i), NaN)
  }
  return(data.frame(Host.Sample = rownames(x),
                    Mean.Deviance = mean.phylogenetic, 
                    Mean.Deviance.SE = se.phylogenetic,
                    Number.of.Symbionts = num.symbionts,
                    Avg.Symbiont.Abundance = read.abund))
}

# phylogenetic.dev <- deviance.phylogenetic(quad.rarefied, utree, null.model = "taxa.labels", iterations = 100, abundance.weighted = TRUE, trim = TRUE, notify = TRUE)
# head(phylogenetic.dev)
# mean(phylogenetic.dev$Mean.Deviance)
# # -0.9141502 vs. -0.9135079

##################################################
# `beta.specificity`: calculate beta-specificity #
##################################################
#' Beta-Specificity
#' 
#' Calculate structural specificity not corrected by null mdoels. 
#'
#' @param x Data frame of hosts populating rows and symbionts populating columns. 
#' 
#' @param index Character. Method for calculation with the Morisita-Horn, Horn or Sorensen Indices. 
#' 
#' @param trim Logical. TRUE removes symbionts that occupy one host sample. FALSE keeps all symbionts. 
#'
#' @return A data frame with symbiont identifiers and beta-specificity values. 
#' @export
#' @examples
#' # Calculate beta-specificity
#' beta.specificity(quad.rarefied, index = "morisita.horn", trim = TRUE, notify = TRUE)
beta.specificity <- function(x, index = c("morisita.horn", "horn", "sorensen"), trim = TRUE, notify = TRUE) {
  # Make holding vectors 
  output.vec <- rep()
  # For every column (symbionts)
  for (j in 1:ncol(x)) {
    # Match index to CqN value
    index = match.arg(index)
    ifelse(index == "sorensen", q <- 0.000000001, ifelse(index == "horn", q <- 0.99999, ifelse(index == "morisita.horn", q <- 2, stop("invalid distance method") )))
    # Subset each column
    col <- x[j]
    colnames(col)[1] <- "Symbiont.Abundance"
    # Make a new row of metadata from the sample names
    col$Sample <- row.names(col)
    # Separate the sample names
    col.sep <- col %>% separate(Sample, c("Host.Group", "Identifier"))
    # Aggregate by Identifier
    col.agg <- aggregate(col.sep[, 1] ~ Identifier, col.sep, sum)
    colnames(col.agg)[2] <- "Identifier.Abundance"
    # Aggregate by identifier
    col.agg.sep <- aggregate(col.sep[, 1] ~ Host.Group, col.sep, sum)
    # Remove Host.Group with zero
    col.agg.sep <- subset(col.agg.sep, col.agg.sep[2] > 0)
    # Filter by hosts that are present
    col.sep <- col.sep[col.sep$Host.Group %in% col.agg.sep$Host.Group, ]
    # Merge separated and aggregated by Identifier
    col.merge <- merge(col.sep, col.agg, by = "Identifier")
    # Calculate relative abundance per species per Identifier
    col.merge$Rel.Abund <- ifelse(is.na(col.merge[, 2] / col.merge[, 4]), 0, col.merge[, 2] / col.merge[, 4])
    # Set host species as factor
    col.merge$Host.Group <- as.factor(col.merge$Host.Group)
    # Calculate T = the number of sites (Identifiers)
    T <- as.numeric(length(unique(col.merge[["Identifier"]])))
    # Calculate CqN 
    # Set holding vector
    prod <- rep()
    sos <- rep()
      # For every symbiont
      for(k in 1:length(levels(col.merge$Host.Group))) {
        # Subset merged data frame by host group
        sp.k <- subset(col.merge, Host.Group == levels(col.merge$Host.Group)[k])
        # Calculate numerator
        prod[k] <- ( sum(sp.k$Rel.Abund)^q) - ( sum(sp.k$Rel.Abund^q) )
        # Calculate sum of squares
        sos[k] <- sum(sp.k$Rel.Abund^q)
      }
    # Calculate index
    output.vec[j] <- ( (1 / (T^q - T)) * (sum ( prod )) / ( (1/T) * (sum(sos)) ) )
    ifelse(notify == TRUE, print(j), NaN)
  }
  # Make data frame
  beta.dat <- data.frame(Symbiont = colnames(x), Similarity.Index = output.vec)
  # Trim noise
  ifelse(trim == TRUE, 
         beta.dat <- subset(beta.dat, Similarity.Index > 0), 
         beta.dat <- beta.dat)
  return(beta.dat)
}

# beta.object <- beta.specificity(quad.rarefied, index = "morisita.horn", trim = TRUE, notify = TRUE)
# beta.object
# mean(beta.object$Similarity.Index)
# mean morisita-horn 0.2008582 vs. 0.2008582

###########################################################
# `null.beta`: calculate null models for beta-specificity #
###########################################################
#' Beta-Specificity Null Models 
#' 
#' Generate null models and calculate beta-specificity per symbiont within each community randomization. 
#'
#' @param x Data frame of hosts populating rows and symbionts populating columns. 
#' 
#' @param index Character. Method for calculation with the Morisita-Horn, Horn or Sorensen Indices. 
#' 
#' @param iterations Integer. Indicate the number of randomized communities to generate. 
#' 
#' @param randomization.method Randomization method. Usage borrowed directly from bipartite::nullmodel. 
#' Specify as "r2dtable", "swap.web", "vaznull", "shuffle.web" or "mgen". 
#' 
#' @param trim Logical. TRUE removes symbionts that occupy one host sample. FALSE keeps all symbionts. 
#' 
#' @param notify Logical. TRUE prints the current iteration of the for loop. NOTE: This function can take some time. 
#'
#' @return A data frame with columns that refer to symbiont identifiers, absolute read abundance, beta-specificities and randomization
#' identifiers. 
#' @export
#' @examples
#' # Generate randomized communities and calculate beta-specificity per symbiont 
#' null.beta(quad.rarefied, index = "morisita.horn", randomization.method = "shuffle.web", iterations = 100, trim = TRUE, notify = TRUE)
null.beta <- function(x, index = c("morisita.horn", "horn", "sorensen"), randomization.method = c("r2dtable", "swap.web", "vaznull", "shuffle.web", "mgen"), iterations = 100, trim = TRUE, notify = TRUE) {
  # Set seed
  set.seed(123)
  # Make 100 randomized communities
  randomization.method <- match.arg(randomization.method)
  null.beta.specificity <- bipartite::nullmodel(x, N = iterations, method = randomization.method)
  # Make holding list
  null.dats <- list()
  # Make holding vectors for Symbiont identifer, read abundance and beta.specificity metric
  Symbiont <- rep()
  Abundance <- rep()
  # Calculate beta-specificity for null models
  for (i in 1:length(null.beta.specificity)) {
    # Call a randomized community
    null <- null.beta.specificity[[i]]
    # Add row and column names
    rownames(null) <- rownames(x)
    colnames(null) <- colnames(x)
    # Make as a data frame
    null <- as.data.frame(null)
    # Calculate beta.specificity per symbiont
    for (j in 1:ncol(null)) {
      # Pull Symbiont name
      Symbiont[j] <- colnames(null)[j]
      # Calculate total read abundance per symbiont
      Abundance[j] <- sum(null[,j])
    }
    # Calculate beta-specificity
    #index <- match.arg(index)
    Beta.Specificity <- beta.specificity(null, index = index, trim = FALSE, notify = FALSE)
    # Make data frame
    null.temp <- data.frame(Abundance, Beta.Specificity)
    # Populate into holding list
    null.dats[[i]] <- null.temp
    # Print iteration 
    ifelse(notify == TRUE, print(i), NaN)
  }
  # Total simulation data
  null.dats.beta <- as.data.frame(do.call("rbind", null.dats))
  # Add randomization number as a new column 
  null.dats.beta$Randomization <- as.factor(rep(1:iterations, each = ncol(x)))
  rownames(null.dats.beta) <- NULL
  # Trim noise
  ifelse(trim == TRUE,
         # only consider symbionts with a multiple-site overlap greater than 0
         null.dats.beta <- subset(null.dats.beta, null.dats.beta$Similarity.Index > 0),
         null.dats.beta <- null.dats.beta)
  # Read out data frame 
  return(data.frame(null.dats.beta))
}

# null.beta.object <- null.beta(quad.rarefied, index = "morisita.horn", randomization.method = "shuffle.web", iterations = 100, trim = TRUE, notify = TRUE)
# null.beta.object
# match 0.1069268 vs. 0.1066938

###############################################################
# `deviance.beta`: calculate the deviance in beta-specificity #
###############################################################
#' Deviance in Beta-Specificity to Null Models 
#' 
#' Calculate the deviance in observed beta-specificity to a null model of beta-specificity per symbiont. 
#' Deviance calculations are measured per symbiont and averaged per host sample. For example, all symbionts within a given host 
#' are evaluated for their host specificity across the entire host community. The host specificities of each symbiont are averaged
#' to calculate the mean host specificity for symbionts within a given host. 
#'
#' @param x Data frame of hosts populating rows and symbionts populating columns. 
#' 
#' @param randomized Data frame. Output from null.beta function. 
#' 
#' @param index Character. Method for calculation with the Morisita-Horn, Horn or Sorensen Indices. 
#' 
#' @param trim Logical. TRUE removes symbionts that occupy one host sample. FALSE keeps all symbionts. 
#' 
#' @param notify Logical. TRUE prints the current iteration of the for loop. 
#'
#' @return A list. First element in the list is a data frame with columns that refer to the host sample identifiers, mean deviance in 
#' beta-specificity, standard error of beta-specificities, number of symbionts per host sample and average symbiont read abundance. 
#' All subsequent elements of the list are plots of uncorrected host specificity as a function of log symbiont read abundance with the null mode
#' in blue relative to the symbionts host specificities in red. 
#' @export
#' @examples
#' # Calculate mean deviance per symbiont per host sample and visualize null vs. observed host specifities 
#' beta.dev <- deviance.beta(quad.rarefied, randomized = null.beta.object, index = "morisita.horn", trim = TRUE, notify = TRUE) 
#' beta.dev[[1]] # View data frame of output 
#' beta.dev[[2]] # View first graph 
#' beta.dev[[81]] # View last graph 
deviance.beta <- function(x, randomized = null.object, index = c("morisita.horn", "horn", "sorensen"), trim = TRUE, notify = TRUE) {
  # Make holding vectors 
  beta.plots <- list()
  mean.beta <- rep()
  se.beta <- rep()
  host.sample <- rep()
  num.symbionts <- rep()
  read.abund <- rep()
  # For every host sample
  for (i in 1:nrow(x)) {
    # Subset a host
    x.sub <- x[i, 1:dim(x)[2], drop = FALSE]
    # Remove symbionts with abundance of zero
    x.sub <- x.sub[ , colSums(x.sub) > 0, drop = FALSE]
    # Save column names 
    x.names <- colnames(x.sub)
    # Filter entire community
    x.input <- x[, colnames(x) %in% x.names, drop = FALSE]
    # Remove rows and columns that sum to zero
    x.input <- as.data.frame(x.input[rowSums(x.input) > 0, colSums(x.input) > 0])
    # Calculate beta-specificity
    index = match.arg(index)
    Similarity.Index <- beta.specificity(x.input, index = index, trim = FALSE, notify = FALSE)
    # Make holding vectors
    Symbiont <- rep()
    Abundance <- rep()
    # For every symbiont
    for (j in 1:ncol(x.input)) {
      # Pull Symbiont name
      Symbiont[j] <- colnames(x.input)[j]
      # Calculate total read abundance per Symbiont
      Abundance[j] <- sum(x.input[,j])
    }
    # Make data frame
    beta.specificity.dat <- data.frame(Symbiont, Abundance, Similarity.Index)
    # Remove noise
    ifelse(trim == TRUE, 
           beta.specificity.dat <- subset(beta.specificity.dat, Similarity.Index > 0), 
           beta.specificity.dat <- beta.specificity.dat) 
    # Plot null vs. empirical per sample
    beta.plots[[i+1]] <- 
      ggplot2::ggplot(beta.specificity.dat, aes(y = Similarity.Index, x = log(Abundance))) +
      geom_point(data = randomized, aes(y = Similarity.Index, x = log(Abundance)), color = "grey", alpha = 0.5, show.legend = TRUE, size = 3) +
      geom_smooth(data = randomized, aes(y = Similarity.Index, x = log(Abundance)), color = "black", method = "lm", se = FALSE, lwd = 1, lty = "dashed", show.legend = FALSE, formula = y ~ x + I(x^2)) + 
      geom_point(color = "red", alpha = 1, show.legend = TRUE, size = 3) +
      stat_poly_eq(data = randomized, parse = TRUE, aes(label = ..eq.label..), formula = y ~ x + I(x^2), label.x = "left", label.y = "bottom", color = "black", size = 5) + 
      theme_bw() +
      ggtitle(rownames(x)[i]) + 
      theme_bw() +
      theme(
        axis.text.x = element_text(size = 13, color = "black"), 
        axis.text.y = element_text(size = 13, color = "black"),
        axis.title.x = element_text(size = 13, margin = margin(t = 5, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 13, margin = margin(t = 0, r = 5, b = 0, l = 0)),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 13, face = "bold.italic"),
        aspect.ratio = 0.85
      ) +
      labs(y = "Uncorrected Structural Specificity", x = "Log Absolute Symbiont Read Abundance")
    # Get model coefficients for null model
    null.eqn <- summary(lm(Similarity.Index ~ log(Abundance) + I(log(Abundance)^2), data = randomized))
    null.eqn$coefficients[1, 1]
    null.eqn$coefficients[2, 1]
    null.eqn$coefficients[3, 1]
    # Calculate mean deviance
    mean.beta[i] <- mean(beta.specificity.dat$Similarity.Index - (null.eqn$coefficients[1, 1] + null.eqn$coefficients[2, 1]*log(beta.specificity.dat$Abundance) + null.eqn$coefficients[3, 1]*log(beta.specificity.dat$Abundance)^2))
    # Calculate standard error of mean deviance
    se.beta[i] <- sd(beta.specificity.dat$Similarity.Index - (null.eqn$coefficients[1, 1] + null.eqn$coefficients[2, 1]*log(beta.specificity.dat$Abundance) + null.eqn$coefficients[3, 1]*log(beta.specificity.dat$Abundance)^2)) / sqrt(length(beta.specificity.dat$Similarity.Index - (null.eqn$coefficients[1, 1] + null.eqn$coefficients[2, 1]*log(beta.specificity.dat$Abundance) + null.eqn$coefficients[3, 1]*log(beta.specificity.dat$Abundance)^2)))
    # Host sample name
    host.sample[i] <- rownames(x)[i]
    # Number of symbionts column
    num.symbionts[i] <- length(beta.specificity.dat$Similarity.Index)
    # Symbiont read abundance 
    read.abund[i] <- mean(beta.specificity.dat$Abundance)
    # Check current iteration
    ifelse(notify == TRUE, print(i), NaN)
  }
  # Populate deviance data into data.frame
  beta.plots[[1]] <- data.frame(Host.Sample = host.sample, 
                                Mean.Deviance = mean.beta, 
                                Mean.Deviance.SE = se.beta,
                                Number.of.Symbionts = num.symbionts,
                                Avg.Symbiont.Abundance = read.abund)
  return(beta.plots)
}



#########################
# Example set in ReadMe #
#########################

# Install lotus
devtools::install_github("austenapigo/lotus", auth_token = "aecbd6a15b658f307c23cbf296f6831b224b2e61")

# Load lotus
library(lotus)

# You can read more about each lotus function with the help function
help("structural.specificity")

# Calculate uncorrected host specificity (not relavitized to a null model)
hs.object <- structural.specificity(quad.rarefied, abundance.weighted = TRUE, trim = TRUE)
hs.object

# Explore data and identify whether negative or variance-decreasing relationships exist between host specificity and symbiont read abundance
plot(density(hs.object$Structural.Specificity)) # plot histogram

read.abund <- as.data.frame(colSums(phylocom$sample)) # get read abundances per symbiont
read.abund.trim <- read.abund[rownames(read.abund) %in% rownames(hs.object), ] # trim relative to hs.object

cor.test(hs.object$Structural.Specificity, read.abund.trim) # correlation test

plot(y = hs.object$Structural.Specificity, x = log(read.abund.trim), ylab = "Uncorrected Structural Specificity (HostRichness)", xlab = "Log Symbiont Read Abundance") # visualize host specificity - read abundance relationships
abline(lm(hs.object$Structural.Specificity~log(read.abund.trim)), col = "red")

# Randomize community matrix to generate a null model for deviance calculations
null.structural.object <- null.structural(quad.rarefied, iterations = 10, abundance.weighted = TRUE, randomization.method = "shuffle.web", trim = TRUE, notify = TRUE)
head(null.structural.object)
null.beta.object <- null.beta(quad.rarefied, index = "morisita.horn", randomization.method = "shuffle.web", iterations = 100, trim = TRUE, notify = TRUE)

# Calculate and plot the deviance of observed host specificity from the null boundary and get averages per host sample
structural.dev <- deviance.structural(quad.rarefied, randomized = null.structural.object, abundance.weighted = TRUE, trim = TRUE, notify = TRUE)
head(structural.dev[[1]]) # View data frame of output
structural.dev[[2]] # View occupancy-abundance model for the first sample
structural.dev[[81]] # View occupancy-abundance model for the last sample
beta.dev <- deviance.beta(quad.rarefied, randomized = null.beta.object, index = "morisita.horn", trim = TRUE, notify = TRUE)
head(beta.dev[[1]]) # View data frame of output
beta.dev[[2]] # View occupancy-abundance model for the first sample
beta.dev[[81]] # View occupancy-abundance model for the last sample
phylo.dev <- deviance.phylogenetic(as.data.frame(t(phylocom$sample)), phylocom$phylo, null.model = "taxa.labels", iterations = 100, abundance.weighted = TRUE, trim = TRUE, notify = TRUE)

# .rs.restartR()
