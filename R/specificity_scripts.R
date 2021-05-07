#######################################################################
# `structural.specificity`: calculate absolute structural specificity #
#######################################################################
#' Absolute Structural Specificity
#' 
#' Calculate absolute structural specificity not relativized by null mdoels. 
#'
#' @param x Data frame. Host by symbiont data frame with hosts populating rows and symbionts populating columns. 
#' 
#' @param abundance.weighted Logical. TRUE calculates Shannon's H per symbiont. FALSE calculates host richness per symbiont. 
#' 
#' @param trim Logical. TRUE removes symbionts that occupy one host species. FALSE keeps all symbionts. 
#'
#' @return A data frame with symbiont identifiers and absolute structural specificity values. 
#' 
#' @details Hosts are labeled by their species name with a period and number identifier (e.g., hostA.1) to differentiate host samples that are of the same species. This naming scheme is required because host specificity is quantified at the level of host species and not host samples.  If this naming scheme does not apply to your experimental design, you should still add  identifiers (e.g., .1, .2, .3, etc.) after each host species or sample identifer.
#' 
#' Host specificity for a symbiont is evaluated across the entire host community. More positive values indicate a narrower symbiont niche and higher host specificity. Structural specificity is negated (multipled by -1) to make this consistent across all metrics.
#' 
#' @references Austen Apigo and Ryoko Oono. 2021. Novel metrics reveal plant abundance, but not plant evolutionary history, shape host specificity in foliar fungal symbionts. In review. 
#' 
#' @export
#' 
#' @importFrom dplyr "%>%"
#' @importFrom tidyr "separate"
#' @importFrom vegan "diversity"
#' @importFrom vegan "specnumber"
#' @importFrom stats "aggregate"
#' 
#' @examples
#' # Calculate host species richness per symbiont
#' host.richness <- structural.specificity(comm.matrix, abundance.weighted = FALSE, trim = TRUE)
#' 
#' # Calculate host species richness per symbiont
#' shannons.h <- structural.specificity(comm.matrix, abundance.weighted = TRUE, trim = TRUE)
structural.specificity <- function(x, abundance.weighted = TRUE, trim = TRUE) {
  # Calculate host richness or Shannon's H
  sh.vector <- rep()
  for (j in 1:ncol(x)) {
    col <- x[j]
    colnames(col)[1] <- "Symbiont.Abundance"
    # Make a new row of metadata from the sample names
    col$Sample <- row.names(col)
    # Separate the sample names
    col.sep <- col %>% tidyr::separate(Sample, c("Host.Species", "Quadrat"))
    # Aggregate by quadrat
    col.agg <- stats::aggregate(col.sep[, 1] ~ Host.Species, col.sep, sum)
    colnames(col.agg)[2] <- "Abundance"
    rownames(col.agg) <- col.agg$Host.Species
    col.agg$Host.Species <- NULL
    # Make vector 
    ifelse(abundance.weighted == TRUE, sh.vector[j] <- vegan::diversity(t(col.agg)), sh.vector[j] <- vegan::specnumber(t(col.agg)))
  }
  # Make data frame
  structural.dat <- data.frame(Symbiont = colnames(x), Structural.Specificity = -1 * sh.vector)
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

################################################################################
# `null.structural`: calculate null models for absolute structural specificity #
################################################################################
#' Absolute Structural Specificity Null Models
#' 
#' Generate null models of structural specificity by calculating absolute structural specificity per symbiont within each community randomization. 
#'
#' @param x Data frame. Host by symbiont data frame with hosts populating rows and symbionts populating columns. 
#' 
#' @param iterations Integer. Indicate the number of randomized communities to generate. 
#' 
#' @param abundance.weighted Logical. TRUE calculates Shannon's H per symbiont. FALSE calculates host richness per symbiont. 
#' 
#' @param randomization.method Randomization method. Usage borrowed directly from bipartite::nullmodel. 
#' Specify as "r2dtable", "swap.web", "vaznull", "shuffle.web" or "mgen". 
#' 
#' @param trim Logical. TRUE removes symbionts that occupy one host species. FALSE keeps all symbionts. 
#' 
#' @param notify Logical. TRUE prints the current iteration of the for loop. NOTE: This function can take some time.
#' 
#' @param randomized.object List. If your preferred null model is not represented in bipartite::nullmodel you can use any other function to generate randomized communities (e.g., vegan::permatswap) outside of `lotus`. Make sure your output is formatted as a list of matrices (match output from bipartite::nullmodel) and you can supply this to the `randomized.object` argument. The `iterations` argument should match the number of iterations in the provided `randomized.object`.
#'
#' @return A data frame with columns that refer to symbiont identifiers, absolute read abundance, absolute structural specificities and randomization identifiers. 
#' 
#' @details Hosts are labeled by their species name with a period and number identifier (e.g., hostA.1) to differentiate host samples that are of the same species. This naming scheme is required because host specificity is quantified at the level of host species and not host samples. If this naming scheme does not apply to your experimental design, you should still add in identifiers (e.g., .1, .2, .3, etc.) after each host species or sample identifer.
#' 
#' Host specificity for a symbiont is evaluated across the entire host community. More positive values indicate a narrower symbiont niche and higher host specificity. Structural specificity is negated (multipled by -1) to make this consistent across all metrics.
#' 
#' @references Austen Apigo and Ryoko Oono. 2021. Novel metrics reveal plant abundance, but not plant evolutionary history, shape host specificity in foliar fungal symbionts. In review. 
#' 
#' @export
#' 
#' @importFrom dplyr "%>%"
#' @importFrom tidyr "separate"
#' @importFrom vegan "diversity"
#' @importFrom vegan "specnumber"
#' @importFrom bipartite "nullmodel"
#' @importFrom stats "aggregate"
#' 
#' @examples
#' # Generate randomized communities and calculate structural specificity per symbiont 
#' \donttest{null.str <- null.structural(comm.matrix, randomization.method = "shuffle.web")}
null.structural <- function(x, iterations = 100, abundance.weighted = TRUE, randomization.method = c("r2dtable", "swap.web", "vaznull", "shuffle.web", "mgen"), trim = TRUE, notify = TRUE, randomized.object = NULL) {
  # Match argument specified
  randomization.method <- match.arg(randomization.method)
  # Set seed
  set.seed(123)
  # Generate 100 randomized communities with option for user to supply their own list of randomized communities
  ifelse(is.null(randomized.object), 
         ifelse(abundance.weighted == TRUE, 
                null.object <- bipartite::nullmodel(x, N = iterations, method = randomization.method), 
                null.object <- bipartite::nullmodel((x > 0) + 0, N = iterations, method = randomization.method)), 
         null.object <- randomized.object)
  # Make holding list
  null.dats <- list()
  # Calculate structural specificity for null models
  for (i in 1:length(null.object)) {
    # Call a randomized community
    null <- null.object[[i]]
    # Add row and column names
    rownames(null) <- rownames(x)
    colnames(null) <- colnames(x)
    # Make as a data frame
    null <- as.data.frame(null)
    # Make holding vectors 
    Symbiont <- rep()
    Abundance <- rep()
    # Subset symbiont name and calculate abundance per symbiont
    for (j in 1:ncol(null)) {
      # Pull symbiont name
      Symbiont[j] <- colnames(null)[j]
      # Calculate total read abundance per symbiont
      Abundance[j] <- sum(null[,j])
    }
    # Calculate structural specificity
    Structural.Specificity <- rep()
    for (j in 1:ncol(null)) {
      col <- null[j]
      colnames(col)[1] <- "Symbiont.Abundance"
      # Make a new row of metadata from the sample names
      col$Sample <- row.names(col)
      # Separate the sample names
      col.sep <- col %>% tidyr::separate(Sample, c("Host.Species", "Quadrat"))
      # Aggregate by host species
      col.agg <- stats::aggregate(col.sep[, 1] ~ Host.Species, col.sep, sum)
      colnames(col.agg)[2] <- "Abundance"
      rownames(col.agg) <- col.agg$Host.Species
      col.agg$Host.Species <- NULL
      # Make vector 
      ifelse(abundance.weighted == TRUE, Structural.Specificity[j] <- -1 * vegan::diversity(t(col.agg)), Structural.Specificity[j] <- -1 * vegan::specnumber(t(col.agg)))
    }
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

####################################################################
# `relative.structural`: calculate relative structural specificity #
####################################################################
#' Relative Structural Specificity
#' 
#' Calculate the deviance in absolute structural specificity to a null model of structural specificity per symbiont and average output per host sample. 
#'
#' @param x Data frame. Host by symbiont data frame with hosts populating rows and symbionts populating columns. 
#' 
#' @param randomized Data frame. Output from null.structural function. 
#' 
#' @param abundance.weighted Logical. TRUE calculates Shannon's H per symbiont. FALSE calculates host richness per symbiont.
#' 
#' @param trim Logical. TRUE removes symbionts that occupy one host species from contributing to host specificity average per host sample. FALSE keeps all symbionts per host sample. 
#' 
#' @param notify Logical. TRUE prints the current iteration of the for loop. 
#'
#' @return A list. First element in the list is a data frame with columns that refer to the host sample identifiers, mean relative structural specificity, standard error of structural specificities, number of symbionts per host sample and average symbiont read abundance. All subsequent elements of the list are plots of absolute host specificity as a function of natural log symbiont read abundance with the null model in black (derived from host specificity of symbionts from randomized communities in grey) and absolute host specificities of the observed symbionts in red. Within graphs, there is an inset equation that refers to the null model expectation. 
#' 
#' @details Hosts are labeled by their species name with a period and number identifier (e.g., hostA.1) to differentiate host samples that are of the same species. This naming scheme is required because host specificity is quantified at the level of host species and not host samples. If this naming scheme does not apply to your experimental design, you should still add in identifiers (e.g., .1, .2, .3, etc.) after each host species or sample identifer.
#' 
#' Host specificity for a symbiont is evaluated across the entire host community. More positive values indicate a narrower symbiont niche and higher host specificity.
#' 
#' Deviance calculations are measured per symbiont and averaged per host sample. The host specificities of each symbiont are averaged to calculate the mean host specificity for symbionts within a given host.
#' 
#' A relative host specificity value greater than zero indicates that an symbiont was more host-specific relative to symbionts with the same read abundances within randomized communities. 
#' 
#' If your preferred null model is not represented in bipartite::nullmodel you can use any other function to generate randomized communities (e.g., vegan::permatswap) outside of `lotus`. Make sure your output is formatted as a list and you can supply this to the `randomized` argument in any deviance-related function.
#' 
#' @export
#' 
#' @import ggplot2
#' @importFrom ggpmisc "stat_poly_eq"
#' @importFrom dplyr "%>%"
#' @importFrom tidyr "separate"
#' @importFrom vegan "diversity"
#' @importFrom vegan "specnumber"
#' @importFrom stats "aggregate"
#' @importFrom stats "lm"
#' @importFrom stats "sd"
#' 
#' @examples
#' # Calculate mean relative structural specificity per symbiont per host sample 
#' \donttest{str.dev <- relative.structural(comm.matrix, randomized = null.str, abundance.weighted = TRUE)}
#' \donttest{str.dev} # View data frame of output 
relative.structural <- function(x, randomized = null.str, abundance.weighted = TRUE, trim = TRUE, notify = TRUE) {
######################### Calculate Absolute Structural Specificity #########################
  # Remove hosts or symbionts that sum to zero
  ifelse(trim == TRUE, x <- x[rowSums(x) > 0, specnumber(x) > 1], x <- x)
  # Make holding vectors 
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
    x.input <- as.data.frame(x.input[rowSums(x.input) > 0, colSums(x.input) > 0, drop = FALSE])
    # Calculate structural specificity
    Structural.Specificity <- rep()
    for (j in 1:ncol(x.input)) {
      col <- x.input[j]
      colnames(col)[1] <- "Symbiont.Abundance"
      # Make a new row of metadata from the sample names
      col$Sample <- row.names(col)
      # Separate the sample names
      col.sep <- col %>% tidyr::separate(Sample, c("Host.Species", "Quadrat"))
      # Aggregate by host species
      col.agg <- stats::aggregate(col.sep[, 1] ~ Host.Species, col.sep, sum)
      colnames(col.agg)[2] <- "Abundance"
      rownames(col.agg) <- col.agg$Host.Species
      col.agg$Host.Species <- NULL
      # Make vector 
      ifelse(abundance.weighted == TRUE, Structural.Specificity[j] <- -1* vegan::diversity(t(col.agg)), Structural.Specificity[j] <- -1 * vegan::specnumber(t(col.agg)))
    }
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
######################### Calculate Standardized Effect Sizes #########################
    # Make holding vectors
    mean.null <- rep()
    sd.null <- rep()
        # For each symbiont in the randomized communities calculate the mean and standard deviation of its host specificity
        for (j in 1:nrow(structural.dat)) {
          subset.null <- subset(randomized, Symbiont == as.character(structural.dat$Symbiont[j]))
          mean.null[j] <- mean(subset.null$Structural.Specificity)
          sd.null[j] <- sd(subset.null$Structural.Specificity)
        }
    # Make final data frame 
    total.df <- data.frame(structural.dat, mean.null, sd.null)
    # Filter if mean SES or SD is greater than 0
    total.df <- subset(total.df, mean.null < 0 & sd.null > 0)
    # Calculate standardized effect sizes per symbiont either compute the average or standard deviation
    mean.structural[i] <- mean((total.df$Structural.Specificity - total.df$mean.null) / total.df$sd.null)
    se.structural[i] <- sd((total.df$Structural.Specificity - total.df$mean.null) / total.df$sd.null) / sqrt(length(total.df$Structural.Specificity))
######################### Make Final Data Frame #########################
    # Host sample name
    host.sample[i] <- rownames(x)[i]
    # Number of symbionts column
    num.symbionts[i] <- length(total.df$Structural.Specificity)
    # Symbiont read abundance 
    read.abund[i] <- mean(total.df$Abundance)
    # Check current iteration
    ifelse(notify == TRUE, print(i), NaN)
    # Populate deviance data into data.frame
    structural.output <- data.frame(Host.Sample = host.sample, 
                                  Mean.Standardized.Effect.Size = mean.structural, 
                                  Standard.Error.of.Mean.SES = se.structural,
                                  Number.of.Symbionts = num.symbionts,
                                  Mean.Symbiont.Abundance = read.abund)
  }
    return(structural.output)
}

#################################################################
# `network.specificity`: calculate absolute network specificity #
#################################################################
#' Network Specificity
#' 
#' Calculate absolute network specificity not relativized by null models. 
#'
#' @param x Data frame. Host and symbiont data frame with hosts populating rows and symbionts populating columns. 
#' 
#' @param abundance.weighted Logical. TRUE calculates the Paired Difference Index per symbiont. FALSE calculates the Resource Range Index per symbiont. 
#' 
#' @param trim Logical. TRUE removes symbionts that occupy one host species. FALSE keeps all symbionts. 
#'
#' @return A data frame with symbiont identifiers and network specificity values. 
#' 
#' @details Hosts are labeled by their species name with a period and number identifier (e.g., hostA.1) to differentiate host samples that are of the same species. This naming scheme is required because host specificity is quantified at the level of host species and not host samples.  If this naming scheme does not apply to your experimental design, you should still add in identifiers (e.g., .1, .2, .3, etc.) after each host species or sample identifer.
#' 
#' Host specificity for a symbiont is evaluated across the entire host community. More positive values indicate a narrower symbiont niche and higher host specificity.
#' 
#' @references Austen Apigo and Ryoko Oono. 2021. Novel metrics reveal plant abundance, but not plant evolutionary history, shape host specificity in foliar fungal symbionts. In review. 
#' 
#' @export
#' 
#' @import ggplot2
#' @importFrom ggpmisc "stat_poly_eq"
#' @importFrom dplyr "%>%"
#' @importFrom tidyr "separate"
#' @importFrom bipartite "nullmodel"
#' @importFrom bipartite "PDI"
#' @importFrom stats "aggregate"
#' 
#' @examples
#' # Calculate Resource Range Index per symbiont
#' rri <- network.specificity(comm.matrix, abundance.weighted = FALSE, trim = TRUE)
#' 
#' # Calculate Paired Difference Index per symbiont
#' pdi <- network.specificity(comm.matrix, abundance.weighted = TRUE, trim = TRUE)
network.specificity <- function(x, abundance.weighted = TRUE, trim = TRUE) {
  # Calculate PDI or RRI
  net.vector <- rep()
  for (j in 1:ncol(x)) {
    col <- x[j]
    colnames(col)[1] <- "Symbiont.Abundance"
    # Make a new row of metadata from the sample names
    col$Sample <- row.names(col)
    # Separate the sample names
    col.sep <- col %>% tidyr::separate(Sample, c("Host.Species", "Quadrat"))
    # Aggregate by quadrat
    col.agg <- stats::aggregate(col.sep[, 1] ~ Host.Species, col.sep, sum)
    colnames(col.agg)[2] <- "Abundance"
    rownames(col.agg) <- col.agg$Host.Species
    col.agg$Host.Species <- NULL
    # Calculate PDI or RRI
    ifelse(abundance.weighted == TRUE, net.vector[j] <- bipartite::PDI(col.agg), net.vector[j] <- bipartite::PDI((col.agg > 0) + 0))
  }
  # Make data frame
  network.dat <- data.frame(Symbiont = colnames(x), Network.Specificity = net.vector)
  # Trim noise
  ifelse(trim == TRUE, 
         # only consider symbionts with a network specificity less than 1
         network.dat <- subset(network.dat, Network.Specificity < 1), 
         network.dat <- network.dat) 
  # Return object
  return(network.dat)
}

#################################################################
# `null.network`: calculate null models for network specificity #
#################################################################
#' Absolute Network Specificity Null Models 
#' 
#' Generate null models of network specificity by calculating absolute network specificity per symbiont within each community randomization. 
#'
#' @param x Data frame. Host by symbiont data frame with hosts populating rows and symbionts populating columns. 
#' 
#' @param iterations Integer. Indicate the number of randomized communities to generate. 
#' 
#' @param abundance.weighted Logical. TRUE calculates the Paired Difference Index per symbiont. FALSE calculates the Resource Range Index per symbiont. 
#' 
#' @param randomization.method Randomization method. Usage borrowed directly from bipartite::nullmodel. 
#' Specify as "r2dtable", "swap.web", "vaznull", "shuffle.web" or "mgen". 
#' 
#' @param trim Logical. TRUE removes symbionts that occupy one host species. FALSE keeps all symbionts. 
#' 
#' @param notify Logical. TRUE prints the current iteration of the for loop. NOTE: This function can take some time.
#' 
#' @param randomized.object List. If your preferred null model is not represented in bipartite::nullmodel you can use any other function to generate randomized communities (e.g., vegan::permatswap) outside of `lotus`. Make sure your output is formatted as a list of matrices (match output from bipartite::nullmodel) and you can supply this to the `randomized.object` argument. The `iterations` argument should match the number of iterations in the provided `randomized.object`.
#'
#' @return A data frame with columns that refer to symbiont identifiers, absolute read abundance, network specificities and randomization identifiers. 
#' 
#' @details Hosts are labeled by their species name with a period and number identifier (e.g., hostA.1) to differentiate host samples that are of the same species. This naming scheme is required because host specificity is quantified at the level of host species and not host samples. If this naming scheme does not apply to your experimental design, you should still add in identifiers (e.g., .1, .2, .3, etc.) after each host species or sample identifer.
#' 
#' Host specificity for a symbiont is evaluated across the entire host community. More positive values indicate a narrower symbiont niche and higher host specificity.
#' 
#' If your preferred null model is not represented in bipartite::nullmodel you can use any other function to generate randomized communities (e.g., vegan::permatswap) outside of `lotus`. Make sure your output is formatted as a list of matrices (match output from bipartite::nullmodel) and you can supply this to the `randomized.object` argument. The `iterations` argument should match the number of iterations in the provided `randomized.object`.
#' 
#' @references Austen Apigo and Ryoko Oono. 2021. Novel metrics reveal plant abundance, but not plant evolutionary history, shape host specificity in foliar fungal symbionts. In review. 
#' 
#' @export
#' 
#' @importFrom dplyr "%>%"
#' @importFrom tidyr "separate"
#' @importFrom bipartite "nullmodel"
#' @importFrom bipartite "PDI"
#' @importFrom stats "aggregate"
#' 
#' @examples
#' # Generate randomized communities and calculate network specificity per symbiont 
#' \donttest{null.net <- null.network(comm.matrix, randomization.method = "shuffle.web")}
null.network <- function(x, iterations = 100, abundance.weighted = TRUE, randomization.method = c("r2dtable", "swap.web", "vaznull", "shuffle.web", "mgen"), trim = TRUE, notify = TRUE, randomized.object = NULL) {
  # Match argument specified
  randomization.method <- match.arg(randomization.method)
  # Set seed
  set.seed(123)
  # Generate 100 randomized communities with option for user to supply their own list of randomized communities
  ifelse(is.null(randomized.object), 
         ifelse(abundance.weighted == TRUE, 
                null.object <- bipartite::nullmodel(x, N = iterations, method = randomization.method), 
                null.object <- bipartite::nullmodel((x > 0) + 0, N = iterations, method = randomization.method)), 
         null.object <- randomized.object)
  # Make holding list
  null.dats <- list()
  # For every randomized community
  for (i in 1:length(null.object)) {
    # Call a randomized community
    null <- null.object[[i]]
    # Add row and column names
    rownames(null) <- rownames(x)
    colnames(null) <- colnames(x)
    # Make as a data frame
    null <- as.data.frame(null)
    # Make holding vectors 
    Symbiont <- rep()
    Abundance <- rep()
    # Subset symbiont name and calculate abundance per symbiont
    for (j in 1:ncol(null)) {
      # Pull symbiont name
      Symbiont[j] <- colnames(null)[j]
      # Calculate total read abundance per symbiont
      Abundance[j] <- sum(null[,j])
    }
    # Calculate network specificity
    Network.Specificity <- rep()
    for (j in 1:ncol(null)) {
      col <- null[j]
      colnames(col)[1] <- "Symbiont.Abundance"
      # Make a new row of metadata from the sample names
      col$Sample <- row.names(col)
      # Separate the sample names
      col.sep <- col %>% separate(Sample, c("Host.Species", "Quadrat"))
      # Aggregate by quadrat
      col.agg <- aggregate(col.sep[, 1] ~ Host.Species, col.sep, sum)
      colnames(col.agg)[2] <- "Abundance"
      rownames(col.agg) <- col.agg$Host.Species
      col.agg$Host.Species <- NULL
      # Make vector 
      ifelse(abundance.weighted == TRUE, Network.Specificity[j] <- PDI(col.agg), Network.Specificity[j] <- PDI((col.agg > 0) + 0))
    }
    # Make data frame
    null.temp <- data.frame(Symbiont, Abundance, Network.Specificity)
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
         # only consider symbionts with a networks pecificity less than 1
         null.dats <- subset(null.dats, null.dats$Network.Specificity < 1),
         null.dats <- null.dats) 
  # Read out data frame 
  return(data.frame(null.dats))
}

##############################################################
# `relative.network`: calculate relative network specificity #
##############################################################
#' Relative Network Specificity 
#' 
#' Calculate the deviance in absolute network specificity to a null model of network specificity per symbiont and average output per host sample. 
#'
#' @param x Data frame. Host by symbiont data frame with hosts populating rows and symbionts populating columns. 
#' 
#' @param randomized Data frame. Output from null.network function. 
#' 
#' @param abundance.weighted Logical. TRUE calculates the Paired Difference Index per symbiont. FALSE calculates the Resource Range Index per symbiont. 
#' 
#' @param trim Logical. TRUE removes symbionts that occupy one host species from contributing to host specificity average per host sample. FALSE keeps all symbionts per host sample. 
#' 
#' @param notify Logical. TRUE prints the current iteration of the for loop. NOTE: This function can take some time. 
#'
#' @return A list. First element in the list is a data frame with columns that refer to the host sample identifiers, mean relative network specificity, standard error of network specificities, number of symbionts per host sample and average symbiont read abundance. All subsequent elements of the list are plots of absolute network specificity as a function of natural log symbiont read abundance with the null model in black (derived from host specificity of symbionts from randomized communities in grey) and the symbiont absolute network specificities in red. Within graphs, there is an inset equation that refers to the null model expectation.
#' 
#' @details Hosts are labeled by their species name with a period and number identifier (e.g., hostA.1) to differentiate host samples that are of the same species. This naming scheme is required because host specificity is quantified at the level of host species and not host samples. If this naming scheme does not apply to your experimental design, you should still add in identifiers (e.g., .1, .2, .3, etc.) after each host species or sample identifer.
#' 
#' Host specificity for a symbiont is evaluated across the entire host community. More positive values indicate a narrower symbiont niche and higher host specificity.
#' 
#' Deviance calculations are measured per symbiont and averaged per host sample. The host specificities of each symbiont are averaged to calculate the mean host specificity for symbionts within a given host. 
#' 
#' A relative host specificity value greater than zero indicates that an symbiont was more host-specific relative to symbionts with the same read abundances within randomized communities. 
#' 
#' If your preferred null model is not represented in bipartite::nullmodel you can use any other function to generate randomized communities (e.g., vegan::permatswap) outside of `lotus`. Make sure your output is formatted as a list and you can supply this to the `randomized` argument in any deviance-related function.
#' 
#' @export
#' 
#' @import ggplot2
#' @importFrom ggpmisc "stat_poly_eq"
#' @importFrom dplyr "%>%"
#' @importFrom tidyr "separate"
#' @importFrom bipartite "PDI"
#' @importFrom stats "aggregate"
#' @importFrom stats "lm"
#' @importFrom stats "sd"
#' 
#' @examples
#' # Calculate mean network specificity per symbiont per host sample
#' \donttest{net.dev <- relative.network(comm.matrix, randomized = null.net, abundance.weighted = TRUE)}
#' \donttest{net.dev} # View data frame of output 
relative.network <- function(x, randomized = null.net, abundance.weighted = TRUE, trim = TRUE, notify = TRUE) {
######################### Calculate Absolute Network Specificity #########################
  # Make holding vectors 
  mean.network <- rep()
  se.network <- rep()
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
    x.input <- as.data.frame(x.input[rowSums(x.input) > 0, colSums(x.input) > 0, drop = FALSE])
    # Calculate network specificity
    Network.Specificity <- rep()
    for (j in 1:ncol(x.input)) {
      col <- x.input[j]
      colnames(col)[1] <- "Symbiont.Abundance"
      # Make a new row of metadata from the sample names
      col$Sample <- row.names(col)
      # Separate the sample names
      col.sep <- col %>% tidyr::separate(Sample, c("Host.Species", "Quadrat"))
      # Aggregate by host species
      col.agg <- stats::aggregate(col.sep[, 1] ~ Host.Species, col.sep, sum)
      colnames(col.agg)[2] <- "Abundance"
      rownames(col.agg) <- col.agg$Host.Species
      col.agg$Host.Species <- NULL
      # Make vector 
      ifelse(abundance.weighted == TRUE, Network.Specificity[j] <- bipartite::PDI(col.agg), Network.Specificity[j] <- bipartite::PDI((col.agg > 0) + 0))
    }
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
    network.dat <- data.frame(Network.Specificity, Symbiont, Abundance)
    # Trim noise
    ifelse(trim == TRUE, 
           # only consider symbionts with a network specificity less than 0
           network.dat <- subset(network.dat, network.dat$Network.Specificity < 1), 
           network.dat <- network.dat) 
    ######################### Calculate Standardized Effect Sizes #########################
    # Make holding vectors
    mean.null <- rep()
    sd.null <- rep()
    # For each symbiont in the randomized communities calculate the mean and standard deviation of its host specificity
    for (j in 1:nrow(network.dat)) {
      subset.null <- subset(randomized, Symbiont == as.character(network.dat$Symbiont[j]))
      mean.null[j] <- mean(subset.null$Network.Specificity)
      sd.null[j] <- sd(subset.null$Network.Specificity)
    }
    # Make final data frame 
    total.df <- data.frame(network.dat, mean.null, sd.null)
    # Filter if mean SES or SD is greater than 0
    total.df <- subset(total.df, mean.null > 0 & sd.null > 0)
    # Calculate standardized effect sizes per symbiont either compute the average or standard deviation
    mean.network[i] <- mean((total.df$Network.Specificity - total.df$mean.null) / total.df$sd.null)
    se.network[i] <- sd((total.df$Network.Specificity - total.df$mean.null) / total.df$sd.null) / sqrt(length(total.df$Network.Specificity))
    ######################### Make Final Data Frame #########################
    # Host sample name
    host.sample[i] <- rownames(x)[i]
    # Number of symbionts column
    num.symbionts[i] <- length(network.dat$Network.Specificity)
    # Symbiont read abundance 
    read.abund[i] <- mean(network.dat$Abundance)
    # Check current iteration
    ifelse(notify == TRUE, print(i), NaN)
    # Populate deviance data into data.frame
    network.output <- data.frame(Host.Sample = host.sample, 
                                 Mean.Standardized.Effect.Size = mean.network, 
                                 Standard.Error.of.Mean.SES = se.network,
                                 Number.of.Symbionts = num.symbionts,
                                 Mean.Symbiont.Abundance = read.abund)
  }
  return(network.output)
}

##################################################################
# `phylogenetic.specificity`: calculate phylogenetic specificity #
##################################################################
#' Absolute Phylogenetic Specificity
#' 
#' Calculate phylogenetic specificity not relativized by null mdoels. 
#'
#' @param x Data frame. Host by symbiont data frame with hosts populating rows and symbionts populating columns. 
#' 
#' @param tree Newick formatted host phylogenetic tree. Cophenetic distances will be calculated the function. 
#' 
#' @param abundance.weighted Logical. TRUE calculates abundance-weighted mean pairwise phylogenetic distance. 
#' FALSE calculates presence-absence mean pairwise phylogenetic distance. 
#' 
#' @param trim Logical. TRUE removes symbionts that occupy one host species. FALSE keeps all symbionts. 
#'
#' @return A data frame with symbiont identifiers and phylogenetic specificity values.
#'
#' @details Hosts are labeled by their species name with a period and number identifier (e.g., hostA.1) to differentiate host samples that are of the same species. This naming scheme is required because host specificity is quantified at the level of host species and not host samples. If this naming scheme does not apply to your experimental design, you should still add in identifiers (e.g., .1, .2, .3, etc.) after each host species or sample identifer.
#' 
#' Host specificity for a symbiont is evaluated across the entire host community. More positive values indicate a narrower symbiont niche and higher host specificity.Phylogenetic specificity is negated (multipled by -1) to make this consistent across all metrics.
#' 
#' @references Austen Apigo and Ryoko Oono. 2021. Novel metrics reveal plant abundance, but not plant evolutionary history, shape host specificity in foliar fungal symbionts. In review. 
#'
#' @export
#' 
#' @importFrom picante "mpd"
#' @importFrom stats "cophenetic"
#' 
#' @examples
#' # Calculate mean pairwise phylogenetic distance per symbiont
#' mpd <- phylogenetic.specificity(comm.matrix, tree = phylo.tree, abundance.weighted = TRUE, trim = TRUE)
phylogenetic.specificity <- function(x, tree, abundance.weighted = TRUE, trim = TRUE) {
  # Calculate abundance-weighted or presence-absence mean pairwise phylogenetic distance
  ifelse(abundance.weighted == TRUE, 
         phylogenetic <- -1 * picante::mpd(t(x), dis = stats::cophenetic(tree), abundance.weighted = TRUE), 
         phylogenetic <- -1 * picante::mpd(t(x), dis = stats::cophenetic(tree), abundance.weighted = FALSE))
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

###############################################################################
# `relative.phylogenetic`: calculate the deviance in phylogenetic specificity #
###############################################################################
#' Deviance in Phylogenetic Specificity to Null Models 
#' 
#' Calculate the deviance in observed phylogenetic specificity to a null model of structural specificity per symbiont and average output per sample. 
#'
#' @param x Data frame. Host by symbiont data frame with hosts populating rows and symbionts populating columns. 
#' 
#' @param tree Newick formatted host phylogenetic tree. Cophenetic distances will be calculated the function. 
#' 
#' @param null.model Randomization method. Usage borrowed directly from picante::mpd 
#' Specify as "taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", "independentswap" or "trialswap". 
#' 
#' @param abundance.weighted Logical. TRUE calculates abundance-weighted mean pairwise phylogenetic distance. 
#' FALSE calculates presence-absence mean pairwise phylogenetic distance. 
#' 
#' @param iterations Integer. Indicate the number of randomized communities to generate. NOTE: This function can take some time. 
#'  
#' @param trim Logical. TRUE removes symbionts that occupy one host species. FALSE keeps all symbionts. 
#' 
#' @param notify Logical. TRUE prints the current iteration of the for loop.
#' 
#' @return A list. First element in the list is a data frame with columns that refer to the host sample identifiers, mean relative phylogenetic specificity, standard error of phylogenetic specificities, number of symbionts per host sample and average symbiont read abundance. All subsequent elements of the list are plots of absolute phylogenetic specificity as a function of natural log symbiont read abundance with the null model in black (derived from host specificity of symbionts from randomized communities in grey) and the symbiont absolute phylogenetic specificities in red. Within graphs, there is an inset equation that refers to the null model expectation.
#' 
#' @details Please make sure your host or sample identifiers match the labels on the host phylogenetic tree. 
#' 
#' Host specificity for a symbiont is evaluated across the entire host community. More positive values indicate a narrower symbiont niche and higher host specificity. 
#' 
#' Deviance calculations are measured per symbiont and averaged per host sample. The host specificities of each symbiont are averaged to calculate the mean host specificity for symbionts within a given host. Unlike relative structural, network and beta-specificity, relative phylogenetic specificity is not calculated as the vertical distance between absolute host specificity and the null model with respect to read abundance. Instead, relative phylogenetic specificity is calculated as the difference between absolute phylogenetic specificity and the mean value of phylogenetic specificity from randomized communities divided by the standard deviation of the null distribution of absolute phylogenetic specificity. 
#' 
#' For phylogenetic specificity, as the number of host species an symbiont occupies approaches the total number of host species in the community, the range of possible phylogenetic specificity values converges to a single value, or the MPD of the entire plant community. This produces a ‘funnel-shaped’ relationship between MPD and host species richness per symbiont and is often correlated with symbiont abundance. To account for this, relative phylogenetic specificity is calculated as the standardized effect size of MPD to quantify the difference between absolute phylogenetic specificity and the mean value of phylogenetic specificity from randomized communities.
#' 
#' A relative host specificity value greater than zero indicates that an symbiont was more host-specific relative to symbionts with the same read abundances within randomized communities. 
#' 
#' see documentation for `picante` for more documentation regarding the ses.mpd function 
#' 
#' @export relative.phylogenetic
#' 
#' @import ggplot2
#' @importFrom ggpmisc "stat_poly_eq"
#' @importFrom dplyr "%>%"
#' @importFrom picante "ses.mpd"
#' @importFrom stats "cophenetic"
#' @importFrom stats "lm"
#' @importFrom stats "sd"
#' 
#' @examples
#' # Calculate mean deviance per symbiont per host sample
#' \donttest{phy.dev <- relative.phylogenetic(comm.matrix, tree = phylo.tree, null.model = "taxa.labels", iterations = 100)}
#' \donttest{phy.dev} # View data frame of output 
relative.phylogenetic <- function(x, tree, null.model = c("taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", "independentswap", "trialswap"), iterations = 100, abundance.weighted = TRUE, trim = TRUE, notify = TRUE) {
######################### Calculate Absolute Phylogenetic Specificity #########################
  # Make holding vectors 
  mean.phylogenetic <- rep()
  se.phylogenetic <- rep()
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
    x.input <- as.data.frame(x.input[rowSums(x.input) > 0, colSums(x.input) > 0, drop = FALSE])
    # Calculate phylogenetic specificity
    null.model <- match.arg(null.model)
    # Set a seed
    set.seed(123)
    # Calculate phylogenetic specificity
    ifelse(abundance.weighted == TRUE, 
           Phylogenetic.Specificity <- picante::ses.mpd(t(x.input), dis = stats::cophenetic(tree), null.model = null.model, abundance.weighted = TRUE, runs = iterations), 
           Phylogenetic.Specificity <- picante::ses.mpd(t(x.input), dis = stats::cophenetic(tree), null.model = null.model, abundance.weighted = FALSE, runs = iterations))
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
    phylogenetic.dat <- data.frame(Symbiont, 
                                   Abundance, 
                                   Richness = Phylogenetic.Specificity$ntaxa, 
                                   Relative.Phylo.Specificity = -1 * Phylogenetic.Specificity$mpd.obs.z,
                                   Absolute.Phylo.Specificity = -1 * Phylogenetic.Specificity$mpd.obs,
                                   Randomized.Phylo.Specificity = -1 * Phylogenetic.Specificity$mpd.rand.mean)
    # Convert NA to 0
    phylogenetic.dat[is.na(phylogenetic.dat)] <- 0
    # Trim noise
    ifelse(trim == TRUE, 
           phylogenetic.dat <- subset(phylogenetic.dat, Absolute.Phylo.Specificity < 0), 
           phylogenetic.dat <- phylogenetic.dat) 
######################### Calculate Relative Phylogenetic Specificity #########################
    # Calculate mean deviance
    mean.phylogenetic[i] <- mean(phylogenetic.dat$Relative.Phylo.Specificity)
    # Calculate the standard error of the mean
    se.phylogenetic[i] <- sd(phylogenetic.dat$Relative.Phylo.Specificity) / sqrt(length(phylogenetic.dat$Relative.Phylo.Specificity))
    # Host sample name
    host.sample[i] <- rownames(x)[i]
    # Number of symbionts column
    num.symbionts[i] <- length(phylogenetic.dat$Relative.Phylo.Specificity)
    # Symbiont read abundance 
    read.abund[i] <- mean(phylogenetic.dat$Abundance)
    # Check current iteration
    ifelse(notify == TRUE, print(i), NaN)
    # Populate deviance data into data.frame
    phylogenetic.output <- data.frame(Host.Sample = host.sample, 
                                        Mean.Deviance = mean.phylogenetic, 
                                        Mean.Relative.SE = se.phylogenetic,
                                        Number.of.Symbionts = num.symbionts,
                                        Avg.Symbiont.Abundance = read.abund)
  }
  return(phylogenetic.output)
}

##################################################
# `beta.specificity`: calculate beta-specificity #
##################################################
#' Beta-Specificity
#' 
#' Calculate beta-specificity not corrected by null models. 
#'
#' @param x Data frame of hosts populating rows and symbionts populating columns. 
#' 
#' @param index Character. Method for calculation with the Morisita-Horn, Horn or Sorensen Indices. Sorensen is presence-absence. Horn and Morisita-Horn are abundance weighted. 
#' 
#' @param trim Logical. TRUE removes symbionts that occupy one host species. FALSE keeps all symbionts. 
#' 
#' @param notify Logical. TRUE prints the current iteration of the for loop. 
#'
#' @return A data frame with symbiont identifiers and beta-specificity values. 
#' 
#' @details Hosts are labeled by their species name with a period and number identifier (e.g., hostA.1) to differentiate host samples that are of the same species. This naming scheme is required because host specificity is quantified at the level of host species and not host samples. If this naming scheme does not apply to your experimental design, you should still add in identifiers (e.g., .1, .2, .3, etc.) after each host species or sample identifer.
#' 
#' Host specificity for a symbiont is evaluated across the entire host community. More positive values indicate a narrower symbiont niche and higher host specificity.
#' 
#' @references Austen Apigo and Ryoko Oono. 2021. Novel metrics reveal plant abundance, but not plant evolutionary history, shape host specificity in foliar fungal symbionts. In review. 
#' 
#' For further detail on the mathematics of beta-specificity, see Jost L, Chao A, Chazdon, L. R. Compositional similarity and β (beta) diversity. In: Magurran AE, McGill BJ, editors. Biological Diversity: Frontiers in Measurement and Assessment. Oxford University Press; 2011. p. 66–84.
#'
#' @export
#' 
#' @importFrom dplyr "%>%"
#' @importFrom tidyr "separate"
#' @importFrom stats "aggregate"
#' 
#' @examples
#' # Calculate Morisita-Horn beta-specificity
#' \donttest{beta.mor <- beta.specificity(comm.matrix, index = "morisita.horn", trim = TRUE, notify = TRUE)}
#' # Calculate Sorensen beta-specificity
#' \donttest{beta.sor <- beta.specificity(comm.matrix, index = "sorensen", trim = TRUE, notify = TRUE)}
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
    col.sep <- col %>% tidyr::separate(Sample, c("Host.Group", "Identifier"))
    # Aggregate by Identifier
    col.agg <- stats::aggregate(col.sep[, 1] ~ Identifier, col.sep, sum)
    colnames(col.agg)[2] <- "Identifier.Abundance"
    # Aggregate by Host.Group
    col.agg.sep <- stats::aggregate(col.sep[, 1] ~ Host.Group, col.sep, sum)
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
  beta.dat <- data.frame(Symbiont = colnames(x), Beta.Specificity = output.vec)
  # Trim noise
  ifelse(trim == TRUE, 
         beta.dat <- subset(beta.dat, Beta.Specificity < 1), 
         beta.dat <- beta.dat)
  return(beta.dat)
}

###########################################################
# `null.beta`: calculate null models for beta-specificity #
###########################################################
#' Beta-Specificity Null Models 
#' 
#' Generate null models of beta-specificity by calculating absolute beta-specificity per symbiont within each community randomization. 
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
#' @param trim Logical. TRUE removes symbionts that occupy one host species. FALSE keeps all symbionts. 
#' 
#' @param notify Logical. TRUE prints the current iteration of the for loop. NOTE: This function can take some time. 
#' 
#' @param randomized.object List. If your preferred null model is not represented in bipartite::nullmodel you can use any other function to generate randomized communities (e.g., vegan::permatswap) outside of `lotus`. Make sure your output is formatted as a list of matrices (match output from bipartite::nullmodel) and you can supply this to the `randomized.object` argument. The `iterations` argument should match the number of iterations in the provided `randomized.object`.
#'
#' @return A data frame with columns that refer to symbiont identifiers, absolute read abundance, beta-specificities and randomization identifiers. 
#' 
#' @details Hosts are labeled by their species name with a period and number identifier (e.g., hostA.1) to differentiate host samples that are of the same species. This naming scheme is required because host specificity is quantified at the level of host species and not host samples. If this naming scheme does not apply to your experimental design, you should still add in identifiers (e.g., .1, .2, .3, etc.) after each host species or sample identifer.
#' 
#' Host specificity for a symbiont is evaluated across the entire host community. More positive values indicate a narrower symbiont niche and higher host specificity.
#' 
#' If your preferred null model is not represented in bipartite::nullmodel you can use any other function to generate randomized communities (e.g., vegan::permatswap) outside of `lotus`. Make sure your output is formatted as a list of matrices (match output from bipartite::nullmodel) and you can supply this to the `randomized.object` argument. The `iterations` argument should match the number of iterations in the provided `randomized.object`.
#' 
#' @references Austen Apigo and Ryoko Oono. 2021. Novel metrics reveal plant abundance, but not plant evolutionary history, shape host specificity in foliar fungal symbionts. In review. 
#' 
#' @export
#' 
#' @importFrom dplyr "%>%"
#' @importFrom bipartite "nullmodel"
#' 
#' @examples
#' # Generate randomized communities and calculate beta-specificity per symbiont 
#' \donttest{null.beta <- null.beta(comm.matrix, index = "morisita.horn", randomization.method = "shuffle.web")}
null.beta <- function(x, iterations = 100, index = c("morisita.horn", "horn", "sorensen"), randomization.method = c("r2dtable", "swap.web", "vaznull", "shuffle.web", "mgen"), trim = TRUE, notify = TRUE, randomized.object = NULL) {
  # Match argument specified
  randomization.method <- match.arg(randomization.method)
  # Set seed
  set.seed(123)
  # Generate 100 randomized communities with option for user to supply their own list of randomized communities
  ifelse(is.null(randomized.object), 
         ifelse(index == "sorensen", 
                null.object <- bipartite::nullmodel((x > 0) + 0, N = iterations, method = randomization.method), 
                null.object <- bipartite::nullmodel(x, N = iterations, method = randomization.method)), 
         null.object <- randomized.object)
  # Make holding list
  null.dats <- list()
  # Calculate beta-specificity for null models
  for (i in 1:length(null.object)) {
    # Call a randomized community
    null <- null.object[[i]]
    # Add row and column names
    rownames(null) <- rownames(x)
    colnames(null) <- colnames(x)
    # Make as a data frame
    null <- as.data.frame(null)
    # Make holding vectors for Symbiont identifer, read abundance and beta.specificity metric
    Symbiont <- rep()
    Abundance <- rep()
    # Calculate beta.specificity per symbiont
    for (j in 1:ncol(null)) {
      # Pull Symbiont name
      Symbiont[j] <- colnames(null)[j]
      # Calculate total read abundance per symbiont
      Abundance[j] <- sum(null[,j])
    }
    # Calculate beta-specificity
    index <- match.arg(index)
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
         # only consider symbionts with a multiple-site overlap less than 1
         null.dats.beta <- subset(null.dats.beta, null.dats.beta$Beta.Specificity < 1),
         null.dats.beta <- null.dats.beta)
  # Read out data frame 
  return(data.frame(null.dats.beta))
}

######################################################################
# `relative.beta`: calculate the deviance in beta-specificity #
######################################################################
#' Deviance in Beta-Specificity to Null Models 
#' 
#' Calculate the deviance in observed beta-specificity to a null model of beta-specificity per symbiont and average output per host sample. 
#'
#' @param x Data frame. Host by symbiont data frame with hosts populating rows and symbionts populating columns. 
#' 
#' @param randomized Data frame. Output from null.beta function. 
#' 
#' @param index Character. Method for calculation with the Morisita-Horn, Horn or Sorensen Indices. 
#' 
#' @param trim Logical. TRUE removes symbionts that occupy one host species from contributing to host specificity average per host sample. FALSE keeps all symbionts per host sample. 
#' 
#' @param notify Logical. TRUE prints the current iteration of the for loop. 
#'
#' @return A list. First element in the list is a data frame with columns that refer to the host sample identifiers, mean relative beta-specificity, standard error of beta-specificities, number of symbionts per host sample and average symbiont read abundance. All subsequent elements of the list are plots of absolute beta-specificity as a function of natural log symbiont read abundance with the null model in black (derived from host specificity of symbionts from randomized communities in grey) and the symbiont absolute beta-specificities in red. Within graphs, there is an inset equation that refers to the null model expectation.
#'
#' @details Hosts are labeled by their species name with a period and number identifier (e.g., hostA.1) to differentiate host samples that are of the same species. This naming scheme is required because host specificity is quantified at the level of host species and not host samples. If this naming scheme does not apply to your experimental design, you should still add in identifiers (e.g., .1, .2, .3, etc.) after each host species or sample identifer.
#' 
#' Host specificity for a symbiont is evaluated across the entire host community. More positive values indicate a narrower symbiont niche and higher host specificity.
#' 
#' Deviance calculations are measured per symbiont and averaged per host sample. The host specificities of each symbiont are averaged to calculate the mean host specificity for symbionts within a given host. 
#' 
#' A relative host specificity value greater than zero indicates that an symbiont was more host-specific relative to symbionts with the same read abundances within randomized communities. 
#' 
#' @references Austen Apigo and Ryoko Oono. 2021. Novel metrics reveal plant abundance, but not plant evolutionary history, shape host specificity in foliar fungal symbionts. In review. 
#' 
#' @export
#' 
#' @import ggplot2
#' @importFrom ggpmisc "stat_poly_eq"
#' @importFrom stats "lm"
#' @importFrom stats "sd"
#' 
#' @examples
#' # Calculate relative beta specificity per symbiont per host sample 
#' \donttest{beta.dev <- relative.beta(comm.matrix, randomized = null.beta, index = "morisita.horn")}
#' \donttest{beta.dev} # View data frame of output 
relative.beta <- function(x, randomized = null.beta, index = c("morisita.horn", "horn", "sorensen"), trim = TRUE, notify = TRUE) {
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
    x.input <- as.data.frame(x.input[rowSums(x.input) > 0, colSums(x.input) > 0, drop = FALSE])
    # Calculate beta-specificity
    index = match.arg(index)
    Beta.Specificity <- beta.specificity(x.input, index = index, trim = FALSE, notify = FALSE)
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
    beta.dat <- data.frame(Symbiont, Abundance, Beta.Specificity)
    # Remove noise
    ifelse(trim == TRUE, 
           beta.dat <- subset(beta.dat, Beta.Specificity > 0), 
           beta.dat <- beta.dat) 
    
    ######################### Calculate Standardized Effect Sizes #########################
    # Make holding vectors
    mean.null <- rep()
    sd.null <- rep()
    # For each symbiont in the randomized communities calculate the mean and standard deviation of its host specificity
    for (j in 1:nrow(beta.dat)) {
      subset.null <- subset(randomized, Symbiont == as.character(beta.dat$Symbiont[j]))
      mean.null[j] <- mean(subset.null$Beta.Specificity)
      sd.null[j] <- sd(subset.null$Beta.Specificity)
    }
    # Make final data frame 
    total.df <- data.frame(beta.dat, mean.null, sd.null)
    # Filter if mean SES or SD is greater than 0
    total.df <- subset(total.df, mean.null > 0 & sd.null > 0)
    # Calculate standardized effect sizes per symbiont either compute the average or standard deviation
    mean.beta[i] <- mean((total.df$Beta.Specificity - total.df$mean.null) / total.df$sd.null)
    se.beta[i] <- sd((total.df$Beta.Specificity - total.df$mean.null) / total.df$sd.null) / sqrt(length(total.df$Beta.Specificity))
    ######################### Make Final Data Frame #########################
    # Host sample name
    host.sample[i] <- rownames(x)[i]
    # Number of symbionts column
    num.symbionts[i] <- length(beta.dat$Beta.Specificity)
    # Symbiont read abundance 
    read.abund[i] <- mean(beta.dat$Abundance)
    # Check current iteration
    ifelse(notify == TRUE, print(i), NaN)
  }
  # Populate deviance data into data.frame
  beta.plots[[1]] <- data.frame(Host.Sample = host.sample, 
                                Mean.Standardized.Effect.Size = mean.beta, 
                                Standard.Error.of.Mean.SES = se.beta,
                                Number.of.Symbionts = num.symbionts,
                                Mean.Symbiont.Abundance = read.abund)
  return(beta.plots)
}
