#############################
# Read in pacakges and data #
#############################
library(tidyverse)
library(ggpmisc)
library(picante)
library(vegan)
library(bipartite)

# Read in data 
dat <- read.csv("example_dat.csv", row.names = 1, header = TRUE)
quad.rarefied <- read.csv("quad_rarefied.csv", row.names = 1, header = TRUE)
utree <- read.tree("utree.txt")

#############################################################
# structural.specificity`: Calculate Structural Specificity #
#############################################################
structural.specificity <- function(x, abundance.weighted = TRUE, trim = TRUE) {
  # Calculate host richness or Shannon's H
  ifelse(abundance.weighted == TRUE, structural <- -1 * diversity(as.data.frame(t(x))), structural <- -1 * specnumber(as.data.frame(t(x))))
  # Make data frame
  structural.dat <- data.frame(Structural.Specificity = structural)
  # Trim noise
  ifelse(trim == TRUE, 
         # if calculating Shannon's H, 
         ifelse(abundance.weighted == TRUE, 
                # only consider symbionts with a Shannon's H less than 0
                structural.dat <- subset(structural.dat, Structural.Specificity < 0), 
                # otherwise, only consider symbionts with a host richness less than -1 
                structural.dat <- subset(structural.dat, Structural.Specificity <-1)), 
         structural.dat <- structural.dat) 
  # Return object
  return(structural.dat)
}

structural.object <- structural.specificity(quad.rarefied, abundance.weighted = TRUE, trim = TRUE)
structural.object
mean(structural.object$Structural.Specificity)
# -0.7925232 matches source code

#######################################################################
# `null.structural`: calculate null models for structural specificity #
#######################################################################
null.structural <- function(x, iterations = 100, abundance.weighted = TRUE, randomization.method = "shuffle.web", notify = TRUE) {
  # Set seed
  set.seed(123)
  # Make 100 randomized communities
  null.structural<- bipartite::nullmodel(x, N = iterations, method = randomization.method)
  # Make holding list
  null.dats <- list()
  # Make holding vectors 
  Symbiont <- rep()
  Abundance <- rep()
  Structural.Specificity <- rep()
  # Calculate structural specificity for null models
  for (i in 1:length(null.structural)) {
    # Call a randomized community
    null <- null.structural[[i]]
    # Add row and column names
    rownames(null) <- rownames(x)
    colnames(null) <- colnames(x)
    # Make as a data frame
    null <- as.data.frame(null)
    # Subset symbiont name and calculate abundance per ASV
    for (j in 1:ncol(null)) {
      # Pull symbiont name
      Symbiont[j] <- colnames(null)[j]
      # Calculate total read abundance per ASV
      Abundance[j] <- sum(null[,j])
    }
    # Calculate structural specificity
    ifelse(abundance.weighted == TRUE, 
           Structural.Specificity <- -1 * diversity(t(null)), 
           Structural.Specificity <- -1 * specnumber(t(null)))
    null.temp <- data.frame(Symbiont, Abundance, Structural.Specificity)
    null.dats[[i]] <- null.temp
    ifelse(notify == TRUE, print(i), NaN)
  }
  # Make into one data frame
  null.dats <- as.data.frame(do.call("rbind", null.dats))
  null.dats$Randomization <- as.factor(rep(1:iterations, each = ncol(x)))
  rownames(null.dats) <- NULL
  # if calculating Shannon's H
  ifelse(abundance.weighted == TRUE, 
         # only consider symbionts with a Shannon's H greater than 0
         null.data.frame <- subset(null.dats, null.dats$Structural.Specificity < 0), 
         # otherwise, only consider symbionts with a host richness greater than 1 
         null.data.frame <- subset(null.dats, null.dats$Structural.Specificity < -1))
  # Read out data frame 
  return(data.frame(null.data.frame))
}

# Generate randomized communities for null model analysis 
null.structural.object <- null.structural(quad.rarefied, iterations = 100, abundance.weighted = TRUE, notify = TRUE)
null.structural.object
mean(null.structural.object$Structural.Specificity)
# I find that I get slightly different answers: -0.6402127 vs. -0.6407677 when run in original project even with set.seed... 

###########################################################################
# `deviance.structural`: calculate the deviance in structural-specificity #
###########################################################################
deviance.structural <- function(x, randomized = null.structural.object, abundance.weighted = TRUE, trim = TRUE, notify = TRUE) {
  # Make holding vectors 
  structural.plots <- list()
  mean.beta <- rep()
  se.beta <- rep()
  # For every host sample
  for (i in 1:nrow(x)) {
    # Subset a host
    x.sub <- x[i, 1:dim(x)[2]]
    # Remove ASVs with abundance of zero
    x.sub <- x.sub[ , colSums(x.sub) > 0]
    # Save column names 
    x.names <- colnames(x.sub)
    # Filter entire community
    x.input <- x[, colnames(x) %in% x.names]
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
      # Calculate total read abundance per ASV
      Abundance[j] <- sum(x.input[,j])
    }
    # Make data frame
    structural.dat <- data.frame(Structural.Specificity, Symbiont, Abundance)
    # Trim noise
    ifelse(trim == TRUE, 
           # if calculating Shannon's H, 
           ifelse(abundance.weighted == TRUE, 
                  # only consider symbionts with a Shannon's H greater than 0
                  structural.dat <- subset(structural.dat, structural.dat$Structural.Specificity < 0), 
                  # otherwise, only consider symbionts with a host richness greater than 1 
                  structural.dat <- subset(structural.dat, structural.dat$Structural.Specificity <-1)), 
           structural.dat <- structural.dat) 
    # Plot null vs. empirical per sample
    structural.plots[[i+1]] <- 
      ggplot(structural.dat, aes(y = Structural.Specificity, x = log(Abundance))) +
      geom_point(data = randomized, aes(y = Structural.Specificity, x = log(Abundance)), color = "blue", alpha = 0.1, show.legend = TRUE, size = 3) +
      geom_point(aes(y = Structural.Specificity, x = log(Abundance)), color = "red", alpha = 1, show.legend = TRUE, size = 3) +
      geom_smooth(aes(y = Structural.Specificity, x = log(Abundance)), color = "red", method = "lm", se = FALSE, lwd = 1, lty = "solid", show.legend = FALSE, formula = y ~ x + I(x^2)) + 
      geom_smooth(data = randomized, aes(y = Structural.Specificity, x = log(Abundance)), color = "blue", method = "lm", se = FALSE, lwd = 1, lty = "solid", show.legend = FALSE, formula = y ~ x + I(x^2)) + 
      stat_poly_eq(data = randomized, parse = TRUE, aes(label = ..eq.label..), formula=  y ~ x + I(x^2), label.x = "left", label.y = "top", color = "black", size = 5) + 
      theme_bw() +
      theme(
        axis.text.x = element_text(size = 15, color = "black"), 
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15, margin = margin(t = 5, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 5, b = 0, l = 0), vjust = 0.5),
        legend.title = element_text(size = 15), 
        legend.text = element_text(size = 15),  
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        panel.border = element_rect(linetype = "solid", size = 1),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold.italic"),
        text = element_text(),
        aspect.ratio = 0.85
      ) + 
      ggtitle(rownames(x)[i]) + 
      labs(y = "Uncorrected Host Specificity", x = "Log Absolute Endophyte Read Abundance")
    # Get model coefficients for null model
    null.eqn <- summary(lm(Structural.Specificity ~ log(Abundance) + I(log(Abundance)^2), data = randomized))
    null.eqn$coefficients[1, 1]
    null.eqn$coefficients[2, 1]
    null.eqn$coefficients[3, 1]
    # Calculate mean deviance
    mean.beta[i] <- mean(structural.dat$Structural.Specificity - (null.eqn$coefficients[1, 1] + null.eqn$coefficients[2, 1]*log(structural.dat$Abundance) + null.eqn$coefficients[3, 1]*log(structural.dat$Abundance)^2))
    # Calculate standard error of mean deviance
    se.beta[i] <- sd(structural.dat$Structural.Specificity - (null.eqn$coefficients[1, 1] + null.eqn$coefficients[2, 1]*log(structural.dat$Abundance) + null.eqn$coefficients[3, 1]*log(structural.dat$Abundance)^2)) / sqrt(length(structural.dat$Structural.Specificity - (null.eqn$coefficients[1, 1] + null.eqn$coefficients[2, 1]*log(structural.dat$Abundance) + null.eqn$coefficients[3, 1]*log(structural.dat$Abundance)^2)))
    # Check current iteration
    ifelse(notify == TRUE, print(i), NaN)
    structural.plots[[1]] <- data.frame(Mean.Deviance = mean.beta, Mean.Deviance.SE = se.beta)
  }
    return(structural.plots)
}

structural.dev <- deviance.structural(quad.rarefied, randomized = null.structural.object, trim = TRUE, notify = TRUE)
structural.dev
structural.dev[[1]]
structural.dev[[2]]
mean(structural.dev[[1]]$Mean.Deviance)
