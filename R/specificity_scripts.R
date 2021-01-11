#######################################################################
# `structural.specificity`: calculate absolute structural specificity #
#######################################################################
#' Absolute Structural Specificity
#' 
#' Calculate absolute structural specificity not relativized by null mdoels. 
#'
#' @param x Data frame. Host x symbiont data frame with hosts populating rows and symbionts populating columns. 
#' 
#' @param abundance.weighted Logical. TRUE calculates Shannon's H per symbiont. FALSE calculates host richness per symbiont. 
#' 
#' @param trim Logical. TRUE removes symbionts that occupy one host species. FALSE keeps all symbionts. 
#'
#' @return A data frame with symbiont identifiers and structural specificity values. 
#' 
#' @details Hosts are labeled by their sampling origin with a period and number identifier after their species name (e.g., hostA.1). This naming scheme is required and helps differentiate host samples that are of the same species. If this naming scheme does not apply to your experimental design, you can add in 'dummy variables' (e.g., .1, .2, .3, etc. after each host species or sample identifer).
#' 
#' Host specificity is calculated for a symbiont across the entire host community. For example, a symbiont found in given host would be evaluated for its host specificity relative to all hosts that are present in a given dataset.
#' 
#' More positive values indicate a narrower symbiont niche and thus higher host specificity. Structural and phylogenetic specificity were negated (multipled by -1) to make this consistent across all metrics.
#' 
#' @export
#' 
#' @import ggplot2
#' @importFrom ggpmisc "stat_poly_eq"
#' @importFrom dplyr "%>%"
#' @importFrom tidyr "separate"
#' @importFrom vegan "diversity"
#' @importFrom vegan "specnumber"
#' @importFrom bipartite "nullmodel"
#' @importFrom bipartite "PDI"
#' @importFrom picante "mpd"
#' @importFrom picante "ses.mpd"
#' @importFrom stats "aggregate"
#' @importFrom stats "cophenetic"
#' @importFrom stats "lm"
#' @importFrom stats "sd"
#' 
#' @examples
#' # Calculate host species richness per symbiont
#' hr <- structural.specificity(quad.rarefied, abundance.weighted = FALSE, trim = TRUE)
structural.specificity <- function(x, abundance.weighted = TRUE, trim = TRUE) {
  # Calculate host richness or Shannon's H
  sh.vector <- rep()
  for (j in 1:ncol(x)) {
    col <- x[j]
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
#' Generate null models and calculate relative structural specificity per symbiont within each community randomization. 
#'
#' @param x Data frame. Host x symbiont data frame with hosts populating rows and symbionts populating columns. 
#' 
#' @param iterations Integer. Indicate the number of randomized communities to generate. 
#' 
#' @param abundance.weighted Logical. TRUE calculates Shannon's H per symbiont. FALSE calculates host richness per symbiont. 
#' 
#' @param randomization.method Randomization method. Usage borrowed directly from bipartite::nullmodel. 
#' Specify as "r2dtable", "swap.web", "vaznull", "shuffle.web" or "mgen". 
#' 
#' @param trim Logical. TRUE removes symbionts that occupy one host species FALSE keeps all symbionts. 
#' 
#' @param notify Logical. TRUE prints the current iteration of the for loop. 
#'
#' @return A data frame with columns that refer to symbiont identifiers, absolute read abundance, structural specificities and randomization
#' identifiers. 
#' 
#' @details Hosts are labeled by their sampling origin with a period and number identifier after their species name (e.g., hostA.1). This naming scheme is required and helps differentiate host samples that are of the same species. If this naming scheme does not apply to your experimental design, you can add in 'dummy variables' (e.g., .1, .2, .3, etc. after each host species or sample identifer).
#' 
#' Host specificity is calculated for a symbiont across the entire host community. For example, a symbiont found in given host would be evaluated for its host specificity relative to all hosts that are present in a given dataset.
#' 
#' More positive relative host specificity values indicate a narrower symbiont niche and thus higher host specificity. Structural and phylogenetic specificity were negated (multipled by -1) to make this consistent across all metrics.
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
#' @importFrom bipartite "nullmodel"
#' @importFrom bipartite "PDI"
#' @importFrom picante "mpd"
#' @importFrom picante "ses.mpd"
#' @importFrom stats "aggregate"
#' @importFrom stats "cophenetic"
#' @importFrom stats "lm"
#' @importFrom stats "sd"
#' 
#' @examples
#' # Generate randomized communities and calculate structural specificity per symbiont 
#' \donttest{null.str <- null.structural(quad.rarefied, randomization.method = "shuffle.web")}
null.structural <- function(x, iterations = 100, abundance.weighted = TRUE, randomization.method = c("r2dtable", "swap.web", "vaznull", "shuffle.web", "mgen"), trim = TRUE, notify = TRUE) {
  # Match argument specified
  randomization.method <- match.arg(randomization.method)
  # Set seed
  set.seed(123)
  # Generate 100 randomized communities
  ifelse(abundance.weighted == TRUE, 
         null.structural <- bipartite::nullmodel(x, N = iterations, method = randomization.method), 
         null.structural <- bipartite::nullmodel((x > 0) + 0, N = iterations, method = randomization.method))
  # Make holding list
  null.dats <- list()
  # Calculate structural specificity for null models
  for (i in 1:length(null.structural)) {
    # Call a randomized community
    null <- null.structural[[i]]
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
      col.sep <- col %>% separate(Sample, c("Host.Species", "Quadrat"))
      # Aggregate by quadrat
      col.agg <- aggregate(col.sep[, 1] ~ Host.Species, col.sep, sum)
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
#' Calculate the deviance in absolute structural specificity to a null model of structural specificity per symbiont. 
#'
#' @param x Data frame. Host x symbiont data frame with hosts populating rows and symbionts populating columns. 
#' 
#' @param randomized Data frame. Output from null.structural function. 
#' 
#' @param abundance.weighted Logical. TRUE calculates Shannon's H per symbiont. FALSE calculates host richness per symbiont.
#' 
#' @param model Character. Specify whether the null expectation should be approximated as a first-(linear) or second-(quadratic) order function. 
#' 
#' @param trim Logical. TRUE removes symbionts that occupy one host species FALSE keeps all symbionts. 
#' 
#' @param notify Logical. TRUE prints the current iteration of the for loop. 
#'
#' @return A list. First element in the list is a data frame with columns that refer to the host sample identifiers, mean deviance in structural specificity, standard error of structural specificities, number of symbionts per host sample and average symbiont read abundance. All subsequent elements of the list are plots of absolute host specificity as a function of natural log symbiont read abundance with the null mode in blue relative to the symbionts host specificities in red. Within graphs, there is an inset equation that refers to the null model expectation. 
#' 
#' @details Hosts are labeled by their sampling origin with a period and number identifier after their species name (e.g., hostA.1). This naming scheme is required and helps differentiate host samples that are of the same species. If this naming scheme does not apply to your experimental design, you can add in 'dummy variables' (e.g., .1, .2, .3, etc. after each host species or sample identifer).
#' 
#' Host specificity is calculated for a symbiont across the entire host community. For example, a symbiont found in given host would be evaluated for its host specificity relative to all hosts that are present in a given dataset. 
#' 
#' Deviance calculations are measured per symbiont and averaged per host sample. The host specificities of each symbiont are averaged to calculate the mean host specificity for symbionts within a given host.
#' 
#' A relative host specificity value greater than zero indicates that an endophyte was more host-specific relative to endophytes with the same read abundances within randomized communities. 
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
#' @importFrom bipartite "nullmodel"
#' @importFrom bipartite "PDI"
#' @importFrom picante "mpd"
#' @importFrom picante "ses.mpd"
#' @importFrom stats "aggregate"
#' @importFrom stats "cophenetic"
#' @importFrom stats "lm"
#' @importFrom stats "sd"
#' 
#' @examples
#' # Calculate mean relative structural specificity per symbiont per host sample 
#' \donttest{str.dev <- relative.structural(quad.rarefied, randomized = null.str, abundance.weighted = TRUE)}
#' \donttest{str.dev[[1]]} # View data frame of output 
#' # Visualize null vs. absolute host specifities 
#' \donttest{str.dev[[2]]} # View first graph 
#' \donttest{str.dev[[81]]} # View last graph 
relative.structural <- function(x, randomized = null.str, abundance.weighted = TRUE, model = c("first", "second"), trim = TRUE, notify = TRUE) {
######################### Calculate Absolute Structural Specificity #########################
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
    Structural.Specificity <- rep()
    for (j in 1:ncol(x.input)) {
      col <- x.input[j]
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
######################### Plot null models and absolute structural specificity #########################
    # Match model argument
    ifelse(model == "first", formula <- "y ~ x", formula <- "y ~ x + I(x^2)")
    # Plot null vs. empirical per sample
    structural.plots[[i+1]] <- 
      ggplot2::ggplot(structural.dat, aes(y = Structural.Specificity, x = log(Abundance))) +
      geom_point(data = randomized, aes(y = Structural.Specificity, x = log(Abundance)), color = "grey", alpha = 0.5, show.legend = TRUE, size = 3) +
      geom_smooth(data = randomized, aes(y = Structural.Specificity, x = log(Abundance)), color = "black", method = "lm", se = FALSE, lwd = 1, lty = "dashed", show.legend = FALSE, formula = formula) + 
      geom_point(color = "red", alpha = 1, show.legend = TRUE, size = 3) +
      stat_poly_eq(data = randomized, parse = TRUE, aes(label = ..eq.label..), formula = formula, label.x = "left", label.y = "bottom", color = "black", size = 5) + 
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
      labs(y = "Absolute Structural Specificity", x = "Log Absolute Symbiont Read Abundance")
######################### Calculate Relative Structural Specificity #########################
    # Get model coefficients for null model
    ifelse(model == "first",
           null.eqn <- summary(lm(Structural.Specificity ~ log(Abundance), data = randomized)), 
           null.eqn <- summary(lm(Structural.Specificity ~ log(Abundance) + I(log(Abundance)^2), data = randomized)))
    # Get null model vector
    ifelse(model == "first",
           null.vector <- null.eqn$coefficients[1, 1] + null.eqn$coefficients[2, 1]*log(structural.dat$Abundance), 
           null.vector <- null.eqn$coefficients[1, 1] + null.eqn$coefficients[2, 1]*log(structural.dat$Abundance) + null.eqn$coefficients[3, 1]*log(structural.dat$Abundance)^2)
    # Calculate mean deviance
    mean.structural[i] <- mean(structural.dat$Structural.Specificity - null.vector)
    # Calculate standard error of mean deviance
    se.structural[i] <- sd(structural.dat$Structural.Specificity - null.vector) / sqrt(length(structural.dat$Structural.Specificity - null.vector))
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
                                        Mean.relative.SE = se.structural,
                                        Number.of.Symbionts = num.symbionts,
                                        Avg.Symbiont.Abundance = read.abund)
  }
    return(structural.plots)
}

#################################################################
# `network.specificity`: calculate absolute network specificity #
#################################################################
#' Network Specificity
#' 
#' Calculate absolute network specificity not relativized by null models. 
#'
#' @param x Data frame. Host x symbiont data frame with hosts populating rows and symbionts populating columns. 
#' 
#' @param abundance.weighted Logical. TRUE calculates the Paired Difference Index per symbiont. FALSE calculates Resource Range Index per symbiont. 
#' 
#' @param abundance.weighted Logical. TRUE calculates the Paired Difference Index per symbiont. FALSE calculates the Resource Range Index per symbiont. 
#' 
#' @param trim Logical. TRUE removes symbionts that occupy one host species. FALSE keeps all symbionts. 
#'
#' @return A data frame with symbiont identifiers and network specificity values. 
#' 
#' @details Hosts are labeled by their sampling origin with a period and number identifier after their species name (e.g., hostA.1). This naming scheme is required and helps differentiate host samples that are of the same species. If this naming scheme does not apply to your experimental design, you can add in 'dummy variables' (e.g., .1, .2, .3, etc. after each host species or sample identifer).
#' 
#' Host specificity is calculated for a symbiont across the entire host community. For example, a symbiont found in given host would be evaluated for its host specificity relative to all hosts that are present in a given dataset. 
#' 
#' Deviance calculations are measured per symbiont and averaged per host sample. The host specificities of each symbiont are averaged to calculate the mean host specificity for symbionts within a given host. 
#' 
#' More positive values indicate a narrower symbiont niche and thus higher host specificity. Structural and phylogenetic specificity were negated (multipled by -1) to make this consistent across all metrics.
#' 
#' @export
#' 
#' @import ggplot2
#' @importFrom ggpmisc "stat_poly_eq"
#' @importFrom dplyr "%>%"
#' @importFrom tidyr "separate"
#' @importFrom vegan "diversity"
#' @importFrom vegan "specnumber"
#' @importFrom bipartite "nullmodel"
#' @importFrom bipartite "PDI"
#' @importFrom picante "mpd"
#' @importFrom picante "ses.mpd"
#' 
#' @examples
#' # Calculate Resource Range Index per symbiont
#' \donttest{rri <- network.specificity(quad.rarefied, abundance.weighted = FALSE, trim = TRUE)}
network.specificity <- function(x, abundance.weighted = TRUE, trim = TRUE) {
  # Calculate PDI or RRI
  net.vector <- rep()
  for (j in 1:ncol(x)) {
    col <- x[j]
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
    # Calculate PDI or RRI
    ifelse(abundance.weighted == TRUE, net.vector[j] <- PDI(col.agg), net.vector[j] <- PDI((col.agg > 0) + 0))
  }
  # Make data frame
  network.dat <- data.frame(Symbiont = colnames(x), Network.Specificity = net.vector)
  # Trim noise
  ifelse(trim == TRUE, 
         # only consider symbionts with a Shannon's H less than 0
         network.dat <- subset(network.dat, Network.Specificity < 1), 
         network.dat <- network.dat) 
  # Return object
  return(network.dat)
}

#################################################################
# `null.network`: calculate null models for network specificity #
#################################################################
#' Network Specificity Null Models 
#' 
#' Generate null models and calculate network specificity per symbiont within each community randomization. 
#'
#' @param x Data frame. Host x symbiont data frame with hosts populating rows and symbionts populating columns. 
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
#' @param notify Logical. TRUE prints the current iteration of the for loop. 
#'
#' @return A data frame with columns that refer to symbiont identifiers, absolute read abundance, network specificities and randomization
#' identifiers. 
#' 
#' @details Hosts are labeled by their sampling origin with a period and number identifier after their species name (e.g., hostA.1). This naming scheme is required and helps differentiate host samples that are of the same species. If this naming scheme does not apply to your experimental design, you can add in 'dummy variables' (e.g., .1, .2, .3, etc. after each host species or sample identifer).
#' 
#' Host specificity is calculated for a symbiont across the entire host community. For example, a symbiont found in given host would be evaluated for its host specificity relative to all hosts that are present in a given dataset. 
#' 
#' Deviance calculations are measured per symbiont and averaged per host sample. The host specificities of each symbiont are averaged to calculate the mean host specificity for symbionts within a given host. 
#' 
#' More positive values indicate a narrower symbiont niche and thus higher host specificity. 
#' 
#' #' If your preferred null model is not represented in bipartite::nullmodel you can use any other function to generate randomized communities (e.g., vegan::permatswap) outside of `lotus`. Make sure your output is formatted as a list and you can supply this to the `randomized` argument in any deviance-related function.
#' 
#' @export
#' 
#' @import ggplot2
#' @importFrom ggpmisc "stat_poly_eq"
#' @importFrom dplyr "%>%"
#' @importFrom tidyr "separate"
#' @importFrom vegan "diversity"
#' @importFrom vegan "specnumber"
#' @importFrom bipartite "nullmodel"
#' @importFrom bipartite "PDI"
#' @importFrom picante "mpd"
#' @importFrom picante "ses.mpd"
#' 
#' @examples
#' # Generate randomized communities and calculate network specificity per symbiont 
#' \donttest{null.net <- null.network(quad.rarefied, randomization.method = "shuffle.web")}
null.network <- function(x, iterations = 100, abundance.weighted = TRUE, randomization.method = c("r2dtable", "swap.web", "vaznull", "shuffle.web", "mgen"), trim = TRUE, notify = TRUE) {
  # Match argument specified
  randomization.method <- match.arg(randomization.method)
  # Set seed
  set.seed(123)
  # Generate 100 randomized communities
  ifelse(abundance.weighted == TRUE, 
         null.network <- bipartite::nullmodel(x, N = iterations, method = randomization.method), 
         null.network <- bipartite::nullmodel((x > 0) + 0, N = iterations, method = randomization.method))
  # Make holding list
  null.dats <- list()
  # For every randomized community
  for (i in 1:length(null.network)) {
    # Call a randomized community
    null <- null.network[[i]]
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
#' Calculate the deviance in absolute network specificity to a null model of network specificity per symbiont. 
#'
#' @param x Data frame. Host x symbiont data frame with hosts populating rows and symbionts populating columns. 
#' 
#' @param randomized Data frame. Output from null.network function. 
#' 
#' @param abundance.weighted Logical. TRUE calculates the Paired Difference Index per symbiont. FALSE calculates the Resource Range Index per symbiont. 
#' 
#' @param model Character. Specify whether the null expectation should be approximated as a first-(linear) or second-(quadratic) order function. 
#' 
#' @param trim Logical. TRUE removes symbionts that occupy one host species. FALSE keeps all symbionts. 
#' 
#' @param notify Logical. TRUE prints the current iteration of the for loop. 
#'
#' @return A list. First element in the list is a data frame with columns that refer to the host sample identifiers, mean deviance in 
#' network specificity, standard error of network specificities, number of symbionts per host sample and average symbiont read abundance. 
#' All subsequent elements of the list are plots of uncorrected host specificity as a function of log symbiont read abundance with the null mode
#' in blue relative to the symbionts host specificities in red. 
#' 
#' @details Hosts are labeled by their sampling origin with a period and number identifier after their species name (e.g., hostA.1). This naming scheme is required and helps differentiate host samples that are of the same species. If this naming scheme does not apply to your experimental design, you can add in 'dummy variables' (e.g., .1, .2, .3, etc. after each host species or sample identifer).
#' 
#' Host specificity is calculated for a symbiont across the entire host community. For example, a symbiont found in given host would be evaluated for its host specificity relative to all hosts that are present in a given dataset. 
#' 
#' Deviance calculations are measured per symbiont and averaged per host sample. The host specificities of each symbiont are averaged to calculate the mean host specificity for symbionts within a given host. 
#' 
#' A relative host specificity value greater than zero indicates that an endophyte was more host-specific relative to endophytes with the same read abundances within randomized communities. 
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
#' @importFrom bipartite "nullmodel"
#' @importFrom bipartite "PDI"
#' @importFrom picante "mpd"
#' @importFrom picante "ses.mpd"
#' @importFrom stats "aggregate"
#' @importFrom stats "cophenetic"
#' @importFrom stats "lm"
#' @importFrom stats "sd"
#' 
#' 
#' @examples
#' # Calculate mean network specificity per symbiont per host sample
#' \donttest{net.dev <- relative.network(quad.rarefied, randomized = null.net, abundance.weighted = TRUE)}
#' \donttest{net.dev[[1]]} # View data frame of output 
#' # Visualize null vs. absolute host specifities 
#' \donttest{net.dev[[2]]} # View first graph 
#' \donttest{net.dev[[81]]} # View last graph 
relative.network <- function(x, randomized = null.net, abundance.weighted = TRUE, model = c("first", "second"), trim = TRUE, notify = TRUE) {
######################### Calculate Absolute Network Specificity #########################
  # Make holding vectors 
  network.plots <- list()
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
    x.input <- as.data.frame(x.input[rowSums(x.input) > 0, colSums(x.input) > 0])
    # Calculate network specificity
    Network.Specificity <- rep()
    for (j in 1:ncol(x.input)) {
      col <- x.input[j]
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
######################### Plot Null Models and Absolute Network Specificity #########################
    # Match model argument
    ifelse(model == "first", formula <- "y ~ x", formula <- "y ~ x + I(x^2)")
    # Plot null vs. empirical per sample
    network.plots[[i+1]] <- 
      ggplot2::ggplot(network.dat, aes(y = Network.Specificity, x = log(Abundance))) +
      geom_point(data = randomized, aes(y = Network.Specificity, x = log(Abundance)), color = "grey", alpha = 0.5, show.legend = TRUE, size = 3) +
      geom_smooth(data = randomized, aes(y = Network.Specificity, x = log(Abundance)), color = "black", method = "lm", se = FALSE, lwd = 1, lty = "dashed", show.legend = FALSE, formula = formula) + 
      geom_point(color = "red", alpha = 1, show.legend = TRUE, size = 3) +
      stat_poly_eq(data = randomized, parse = TRUE, aes(label = ..eq.label..), formula = formula, label.x = "left", label.y = "bottom", color = "black", size = 5) + 
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
      labs(y = "Absolute Network Specificity", x = "Log Absolute Symbiont Read Abundance")
######################### Calculate Relative Network Specificity #########################
    # Get model coefficients for null model
    ifelse(model == "first",
           null.eqn <- summary(lm(Network.Specificity ~ log(Abundance), data = randomized)), 
           null.eqn <- summary(lm(Network.Specificity ~ log(Abundance) + I(log(Abundance)^2), data = randomized)))
    # Get null model vector
    ifelse(model == "first",
           null.vector <- null.eqn$coefficients[1, 1] + null.eqn$coefficients[2, 1]*log(network.dat$Abundance), 
           null.vector <- null.eqn$coefficients[1, 1] + null.eqn$coefficients[2, 1]*log(network.dat$Abundance) + null.eqn$coefficients[3, 1]*log(network.dat$Abundance)^2)
    # Calculate mean deviance
    mean.network[i] <- mean(network.dat$Network.Specificity - null.vector)
    # Calculate standard error of mean deviance
    se.network[i] <- sd(network.dat$Network.Specificity - null.vector) / sqrt(length(network.dat$Network.Specificity - null.vector))
    # Host sample name
    host.sample[i] <- rownames(x)[i]
    # Number of symbionts column
    num.symbionts[i] <- length(network.dat$Network.Specificity)
    # Symbiont read abundance 
    read.abund[i] <- mean(network.dat$Abundance)
    # Check current iteration
    ifelse(notify == TRUE, print(i), NaN)
    # Populate deviance data into data.frame
    network.plots[[1]] <- data.frame(Host.Sample = host.sample, 
                                        Mean.Deviance = mean.network, 
                                        Mean.relative.SE = se.network,
                                        Number.of.Symbionts = num.symbionts,
                                        Avg.Symbiont.Abundance = read.abund)
  }
  return(network.plots)
}

##################################################################
# `phylogenetic.specificity`: calculate phylogenetic specificity #
##################################################################
#' Absolute Phylogenetic Specificity
#' 
#' Calculate phylogenetic specificity not relativized by null mdoels. 
#'
#' @param x Data frame. Host x symbiont data frame with hosts populating rows and symbionts populating columns. 
#' 
#' @param utree Newick formatted host phylogenetic tree. Cophenetic distances will be calculated the function. 
#' 
#' @param abundance.weighted Logical. TRUE calculates abundance-weighted mean pairwise phylogenetic distance. 
#' FALSE calculates presence-absence mean pairwise phylogenetic distance. 
#' 
#' @param trim Logical. TRUE removes symbionts that occupy one host species. FALSE keeps all symbionts. 
#'
#' @return A data frame with symbiont identifiers and phylogenetic specificity values.
#'
#' @details Hosts are labeled by their sampling origin with a period and number identifier after their species name (e.g., hostA.1). This naming scheme is required and helps differentiate host samples that are of the same species. If this naming scheme does not apply to your experimental design, you can add in 'dummy variables' (e.g., .1, .2, .3, etc. after each host species or sample identifer).
#' 
#' Host specificity is calculated for a symbiont across the entire host community. For example, a symbiont found in given host would be evaluated for its host specificity relative to all hosts that are present in a given dataset.
#' 
#' More positive values indicate a narrower symbiont niche and thus higher host specificity. Structural and phylogenetic specificity were negated (multipled by -1) to make this consistent across all metrics.
#'
#' @export
#' 
#' @import ggplot2
#' @importFrom ggpmisc "stat_poly_eq"
#' @importFrom dplyr "%>%"
#' @importFrom tidyr "separate"
#' @importFrom vegan "diversity"
#' @importFrom vegan "specnumber"
#' @importFrom bipartite "nullmodel"
#' @importFrom bipartite "PDI"
#' @importFrom picante "mpd"
#' @importFrom picante "ses.mpd"
#' @importFrom stats "aggregate"
#' @importFrom stats "cophenetic"
#' @importFrom stats "lm"
#' @importFrom stats "sd"
#' 
#' @examples
#' # Calculate mean pairwise phylogenetic distance per symbiont
#' \donttest{mpd <- phylogenetic.specificity(quad.rarefied, utree, abundance.weighted = TRUE, trim = TRUE)}
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

###############################################################################
# `relative.phylogenetic`: calculate the deviance in phylogenetic specificity #
###############################################################################
#' Deviance in Phylogenetic Specificity to Null Models 
#' 
#' Calculate the deviance in observed phylogenetic specificity to a null model of structural specificity per symbiont.
#'
#' @param x Data frame. Host x symbiont data frame with hosts populating rows and symbionts populating columns. 
#' 
#' @param utree Newick formatted host phylogenetic tree. Cophenetic distances will be calculated the function. 
#' 
#' @param null.model Randomization method. Usage borrowed directly from picante::mpd 
#' Specify as "taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", "independentswap" or "trialswap". 
#' 
#' @param abundance.weighted Logical. TRUE calculates abundance-weighted mean pairwise phylogenetic distance. 
#' FALSE calculates presence-absence mean pairwise phylogenetic distance. 
#' 
#' @param iterations Integer. Indicate the number of randomized communities to generate. 
#' 
#' @param model Character. Specify whether the null expectation should be approximated as a first-(linear) or second-(quadratic) order function.
#'  
#' @param trim Logical. TRUE removes symbionts that occupy one host species. FALSE keeps all symbionts. 
#' 
#' @param notify Logical. TRUE prints the current iteration of the for loop.
#'
#' @return A data frame with columns that refer to the host sample identifiers, mean deviance in phylogenetic specificity, 
#' standard error of structural specificities, number of symbionts per host sample and average symbiont read abundance. 
#' 
#' @details Please make sure your host or sample identifiers match the labels on the host phylogenetic tree. 
#' 
#' Host specificity is calculated for a symbiont across the entire host community. For example, a symbiont found in given host would be evaluated for its host specificity relative to all hosts that are present in a given dataset. 
#' 
#' Deviance calculations are measured per symbiont and averaged per host sample. The host specificities of each symbiont are averaged to calculate the mean host specificity for symbionts within a given host. Unlike relative structural and network specificity, relative phylogenetic specificity is not calculated as the vertical distance between absolute host specificity and the null model with respect to read abundance. Instead, relative phylogenetic specificity is calculated as the difference between absolute phylogenetic specificity and the mean value of phylogenetic specificity from randomized communities divided by the standard deviation of the null distribution of absolute phylogenetic specificity. 
#' 
#' A relative host specificity value greater than zero indicates that an endophyte was more host-specific relative to endophytes with the same read abundances within randomized communities. 
#' 
#' see documentation for `picante` for more documentation regarding the ses.mpd function 
#' 
#' @export relative.phylogenetic
#' 
#' @import ggplot2
#' @importFrom ggpmisc "stat_poly_eq"
#' @importFrom dplyr "%>%"
#' @importFrom tidyr "separate"
#' @importFrom vegan "diversity"
#' @importFrom vegan "specnumber"
#' @importFrom bipartite "nullmodel"
#' @importFrom bipartite "PDI"
#' @importFrom picante "mpd"
#' @importFrom picante "ses.mpd"
#' @importFrom stats "aggregate"
#' @importFrom stats "cophenetic"
#' @importFrom stats "lm"
#' @importFrom stats "sd"
#' 
#' @examples
#' # Calculate mean deviance per symbiont per host sample
#' \donttest{phy.dev <- relative.phylogenetic(quad.rarefied, utree, null.model = "taxa.labels", iterations = 100)}
#' \donttest{phy.dev[[1]]} # View data frame of output 
#' # Visualize null vs. absolute host specifities 
#' \donttest{phy.dev[[2]]} # View first graph 
#' \donttest{phy.dev[[81]]} # View last graph 
relative.phylogenetic <- function(x, utree, null.model = c("taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", "independentswap", "trialswap"), iterations = 100, abundance.weighted = TRUE, model = c("first", "second"), trim = TRUE, notify = TRUE) {
######################### Calculate Absolute Phylogenetic Specificity #########################
  # Make holding vectors 
  phylogenetic.plots <- list()
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
    phylogenetic.dat <- data.frame(Symbiont, 
                                   Abundance, 
                                   Richness = Phylogenetic.Specificity$ntaxa, 
                                   Phylogenetic.Specificity = -1 * Phylogenetic.Specificity$mpd.obs.z,
                                   Observed.Phylo.Specificity = -1 * Phylogenetic.Specificity$mpd.obs,
                                   Randomized.Phylo.Specificity = -1 * Phylogenetic.Specificity$mpd.rand.mean)
    # Convert NA to 0
    phylogenetic.dat[is.na(phylogenetic.dat)] <- 0
    # Trim noise
    ifelse(trim == TRUE, 
           phylogenetic.dat <- subset(phylogenetic.dat, Observed.Phylo.Specificity < 0), 
           phylogenetic.dat <- phylogenetic.dat) 
######################### Plot Null Models and Absolute Phylogenetic Specificity #########################
    # Match model argument
    ifelse(model == "first", formula <- "y ~ x", formula <- "y ~ x + I(x^2)")
    # Plot null vs. empirical per sample
    phylogenetic.plots[[i+1]] <- 
      ggplot2::ggplot(phylogenetic.dat, aes(y = Observed.Phylo.Specificity, x = log(Abundance))) +
      geom_point(aes(y = Randomized.Phylo.Specificity, x = log(Abundance)), color = "grey", alpha = 0.5, show.legend = TRUE, size = 3) +
      geom_smooth(aes(y = Randomized.Phylo.Specificity, x = log(Abundance)), color = "black", method = "lm", se = FALSE, lwd = 1, lty = "dashed", show.legend = FALSE, formula = formula) + 
      geom_point(color = "red", alpha = 1, show.legend = TRUE, size = 3) +
      stat_poly_eq(parse = TRUE, aes(y = Randomized.Phylo.Specificity, x = log(Abundance), label = ..eq.label..), formula = formula, label.x = "left", label.y = "bottom", color = "black", size = 5) + 
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
      labs(y = "Absolute Phylogenetic Specificity", x = "Log Absolute Symbiont Read Abundance")
######################### Calculate Relative Network Specificity #########################
    # Calculate mean deviance
    mean.phylogenetic[i] <- mean(phylogenetic.dat$Phylogenetic.Specificity)
    # Calculate the standard error of the mean
    se.phylogenetic[i] <- sd(phylogenetic.dat$Phylogenetic.Specificity) / sqrt(length(phylogenetic.dat$Phylogenetic.Specificity))
    # Host sample name
    host.sample[i] <- rownames(x)[i]
    # Number of symbionts column
    num.symbionts[i] <- length(phylogenetic.dat$Phylogenetic.Specificity)
    # Symbiont read abundance 
    read.abund[i] <- mean(phylogenetic.dat$Abundance)
    # Check current iteration
    ifelse(notify == TRUE, print(i), NaN)
    # Populate deviance data into data.frame
    phylogenetic.plots[[1]] <- data.frame(Host.Sample = host.sample, 
                                        Mean.Deviance = mean.phylogenetic, 
                                        Mean.relative.SE = se.phylogenetic,
                                        Number.of.Symbionts = num.symbionts,
                                        Avg.Symbiont.Abundance = read.abund)
  }
  return(phylogenetic.plots)
}

##################################################
# `beta.specificity`: calculate beta-specificity #
##################################################
#' Beta-Specificity
#' 
#' Calculate structural specificity not corrected by null models. 
#'
#' @param x Data frame of hosts populating rows and symbionts populating columns. 
#' 
#' @param index Character. Method for calculation with the Morisita-Horn, Horn or Sorensen Indices. 
#' 
#' @param trim Logical. TRUE removes symbionts that occupy one host species. FALSE keeps all symbionts. 
#' 
#' @param notify Logical. TRUE prints the current iteration of the for loop. 
#'
#' @return A data frame with symbiont identifiers and beta-specificity values. 
#' 
#' @details Hosts are labeled by their sampling origin with a period and number identifier after their species name (e.g., hostA.1). This naming scheme is required and helps differentiate host samples that are of the same species. If this naming scheme does not apply to your experimental design, you can add in 'dummy variables' (e.g., .1, .2, .3, etc. after each host species or sample identifer).
#' 
#' Host specificity is calculated for a symbiont across the entire host community. For example, a symbiont found in given host would be evaluated for its host specificity relative to all hosts that are present in a given dataset.
#' 
#' More positive values indicate a narrower symbiont niche and thus higher host specificity. 
#'
#' @export
#' 
#' @import ggplot2
#' @importFrom ggpmisc "stat_poly_eq"
#' @importFrom dplyr "%>%"
#' @importFrom tidyr "separate"
#' @importFrom vegan "diversity"
#' @importFrom vegan "specnumber"
#' @importFrom bipartite "nullmodel"
#' @importFrom bipartite "PDI"
#' @importFrom picante "mpd"
#' @importFrom picante "ses.mpd"
#' @importFrom stats "aggregate"
#' @importFrom stats "cophenetic"
#' @importFrom stats "lm"
#' @importFrom stats "sd"
#' 
#' @examples
#' # Calculate beta-specificity
#' \donttest{beta.mh <- beta.specificity(quad.rarefied, index = "morisita.horn", trim = TRUE, notify = TRUE)}
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
  beta.dat <- data.frame(Symbiont = colnames(x), Beta.Specificity = output.vec)
  # Trim noise
  ifelse(trim == TRUE, 
         beta.dat <- subset(beta.dat, Beta.Specificity > 0), 
         beta.dat <- beta.dat)
  return(beta.dat)
}

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
#' @param trim Logical. TRUE removes symbionts that occupy one host species. FALSE keeps all symbionts. 
#' 
#' @param notify Logical. TRUE prints the current iteration of the for loop. NOTE: This function can take some time. 
#'
#' @return A data frame with columns that refer to symbiont identifiers, absolute read abundance, beta-specificities and randomization
#' identifiers. 
#' 
#' @export
#' 
#' @import ggplot2
#' @importFrom ggpmisc "stat_poly_eq"
#' @importFrom dplyr "%>%"
#' @importFrom tidyr "separate"
#' @importFrom vegan "diversity"
#' @importFrom vegan "specnumber"
#' @importFrom bipartite "nullmodel"
#' @importFrom bipartite "PDI"
#' @importFrom picante "mpd"
#' @importFrom picante "ses.mpd"
#' @importFrom stats "aggregate"
#' @importFrom stats "cophenetic"
#' @importFrom stats "lm"
#' @importFrom stats "sd"
#' 
#' @examples
#' # Generate randomized communities and calculate beta-specificity per symbiont 
#' \donttest{null.beta <- null.beta(quad.rarefied, index = "morisita.horn", randomization.method = "shuffle.web")}
null.beta <- function(x, iterations = 100, index = c("morisita.horn", "horn", "sorensen"), randomization.method = c("r2dtable", "swap.web", "vaznull", "shuffle.web", "mgen"), trim = TRUE, notify = TRUE) {
  # Match argument specified
  randomization.method <- match.arg(randomization.method)
  # Set seed
  set.seed(123)
  # Make 100 randomized communities
  ifelse(index == "sorensen", 
         null.beta.specificity <- bipartite::nullmodel((x > 0) + 0, N = iterations, method = randomization.method), 
         null.beta.specificity <- bipartite::nullmodel(x, N = iterations, method = randomization.method))
  # Make holding list
  null.dats <- list()
  # Calculate beta-specificity for null models
  for (i in 1:length(null.beta.specificity)) {
    # Call a randomized community
    null <- null.beta.specificity[[i]]
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
         # only consider symbionts with a multiple-site overlap greater than 0
         null.dats.beta <- subset(null.dats.beta, null.dats.beta$Beta.Specificity > 0),
         null.dats.beta <- null.dats.beta)
  # Read out data frame 
  return(data.frame(null.dats.beta))
}

###############################################################
# `relative.beta`: calculate the deviance in beta-specificity #
###############################################################
#' Deviance in Beta-Specificity to Null Models 
#' 
#' Calculate the deviance in observed beta-specificity to a null model of beta-specificity per symbiont. 
#' Deviance calculations are measured per symbiont and averaged per host sample. For example, all symbionts within a given host 
#' are evaluated for their host specificity across the entire host community. The host specificities of each symbiont are averaged
#' to calculate the mean host specificity for symbionts within a given host. 
#'
#' @param x Data frame. Host x symbiont data frame with hosts populating rows and symbionts populating columns. 
#' 
#' @param randomized Data frame. Output from null.beta function. 
#' 
#' @param index Character. Method for calculation with the Morisita-Horn, Horn or Sorensen Indices. 
#' 
#' @param model Character. Specify whether the null expectation should be approximated as a first-(linear) or second-(quadratic) order function. 
#' 
#' @param trim Logical. TRUE removes symbionts that occupy one host species. FALSE keeps all symbionts. 
#' 
#' @param notify Logical. TRUE prints the current iteration of the for loop. 
#'
#' @details Hosts are labeled by their sampling origin with a period and number identifier after their species name (e.g., hostA.1). This naming scheme is required and helps differentiate host samples that are of the same species. If this naming scheme does not apply to your experimental design, you can add in 'dummy variables' (e.g., .1, .2, .3, etc. after each host species or sample identifer).
#' 
#' Host specificity is calculated for a symbiont across the entire host community. For example, a symbiont found in given host would be evaluated for its host specificity relative to all hosts that are present in a given dataset. 
#' 
#' Deviance calculations are measured per symbiont and averaged per host sample. The host specificities of each symbiont are averaged to calculate the mean host specificity for symbionts within a given host. 
#' 
#' A relative host specificity value greater than zero indicates that an endophyte was more host-specific relative to endophytes with the same read abundances within randomized communities. 
#' 
#' @return A list. First element in the list is a data frame with columns that refer to the host sample identifiers, mean deviance in 
#' beta-specificity, standard error of beta-specificities, number of symbionts per host sample and average symbiont read abundance. 
#' All subsequent elements of the list are plots of uncorrected host specificity as a function of log symbiont read abundance with the null mode in blue relative to the symbionts host specificities in red. 
#' 
#' @export
#' 
#' @import ggplot2
#' @importFrom ggpmisc "stat_poly_eq"
#' @importFrom dplyr "%>%"
#' @importFrom tidyr "separate"
#' @importFrom vegan "diversity"
#' @importFrom vegan "specnumber"
#' @importFrom bipartite "nullmodel"
#' @importFrom bipartite "PDI"
#' @importFrom picante "mpd"
#' @importFrom picante "ses.mpd"
#' @importFrom stats "aggregate"
#' @importFrom stats "cophenetic"
#' @importFrom stats "lm"
#' @importFrom stats "sd"
#' 
#' @examples
#' # Calculate relative beta specificity per symbiont per host sample 
#' \donttest{beta.dev <- relative.beta(quad.rarefied, randomized = null.beta, index = "morisita.horn")}
#' \donttest{beta.dev[[1]]} # View data frame of output 
#' # Visualize null vs. absolute host specifities 
#' \donttest{beta.dev[[2]]} # View first graph 
#' \donttest{beta.dev[[81]]} # View last graph 
relative.beta <- function(x, randomized = null.beta, index = c("morisita.horn", "horn", "sorensen"), model = c("first", "second"), trim = TRUE, notify = TRUE) {
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
    # Match model argument
    ifelse(model == "first", formula <- "y ~ x", formula <- "y ~ x + I(x^2)")
    # Plot null vs. empirical per sample
    beta.plots[[i+1]] <- 
      ggplot2::ggplot(beta.dat, aes(y = Beta.Specificity, x = log(Abundance))) +
      geom_point(data = randomized, aes(y = Beta.Specificity, x = log(Abundance)), color = "grey", alpha = 0.5, show.legend = TRUE, size = 3) +
      geom_smooth(data = randomized, aes(y = Beta.Specificity, x = log(Abundance)), color = "black", method = "lm", se = FALSE, lwd = 1, lty = "dashed", show.legend = FALSE, formula = formula) + 
      geom_point(color = "red", alpha = 1, show.legend = TRUE, size = 3) +
      stat_poly_eq(data = randomized, parse = TRUE, aes(label = ..eq.label..), formula = formula, label.x = "left", label.y = "bottom", color = "black", size = 5) + 
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
      labs(y = "Absolute Structural Specificity", x = "Log Absolute Symbiont Read Abundance")
    # Get model coefficients for null model
    ifelse(model == "first",
           null.eqn <- summary(lm(Beta.Specificity ~ log(Abundance), data = randomized)), 
           null.eqn <- summary(lm(Beta.Specificity ~ log(Abundance) + I(log(Abundance)^2), data = randomized)))
    # Get null model vector
    ifelse(model == "first",
           null.vector <- null.eqn$coefficients[1, 1] + null.eqn$coefficients[2, 1]*log(beta.dat$Abundance), 
           null.vector <- null.eqn$coefficients[1, 1] + null.eqn$coefficients[2, 1]*log(beta.dat$Abundance) + null.eqn$coefficients[3, 1]*log(beta.dat$Abundance)^2)
    # Calculate mean deviance
    mean.beta[i] <- mean(beta.dat$Beta.Specificity - null.vector)
    # Calculate standard error of mean deviance
    se.beta[i] <- sd(beta.dat$Beta.Specificity - null.vector) / sqrt(length(beta.dat$Beta.Specificity - null.vector))
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
                                Mean.Deviance = mean.beta, 
                                Mean.relative.SE = se.beta,
                                Number.of.Symbionts = num.symbionts,
                                Avg.Symbiont.Abundance = read.abund)
  return(beta.plots)
}