#' Community matrix for host specificity analyses
#'
#' The study took place in the Mer Bleue peat bog of Ottawa, Canada in 1973. The paper is a preliminary evaluation of the pollination relationships of the major entomophilous plant species of the Mer Bleue.
#' 
#' The authors recorded their data by counting the number of individual flower visitors caught on each plant species. The total number of individuals collected on each plant species provide a rough estimate of the level of visitation that each species received. Data are presented as an interaction frequency matrix, in which cells with positive integers indicate the frequency of interaction between a pair of species, and cells with zeros indicate no interaction.
#' 
#' Description taken from (https://iwdb.nceas.ucsb.edu/html/small_1976.html). 
#' 
#' Note: This interaction matrix has been modified such that each host species is composed of two host samples. This change was made to make this dataset compatible for the measurement of beta-specificity.
#'
#' @docType data
#'
#' @usage data(comm.matrix)
#'
#' @format A data frame
#'
#' @keywords datasets
#'
#' @references Small, E. 1976. Insect pollinators of the Mer Bleue peat bog of Ottawa. Canadian Field Naturalist 90:22-28
#'
#' @source Small, E. 1976. Insect pollinators of the Mer Bleue peat bog of Ottawa. Canadian Field Naturalist 90:22-28
#'
#' @examples
#' data(comm.matrix)
#' 
#' structural.specificity(comm.matrix)
"comm.matrix"
