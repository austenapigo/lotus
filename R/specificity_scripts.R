#############################
# Read in pacakges and data #
#############################
library(tidyverse)
library(ggpmisc)
library(picante)
library(vegan)

# Read in data 
dat <- read.csv("example_dat.csv", row.names = 1, header = TRUE)
quad.rarefied <- read.csv("quad_rarefied.csv", row.names = 1, header = TRUE)
utree <- read.tree("utree.txt")