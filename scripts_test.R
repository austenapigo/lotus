structural.specificity <- function(data = x, abundance.weighted = TRUE, trim = TRUE) {
  
  ifelse(abundance.weighted == TRUE, structural <- -1 * diversity(t(x)), structural <- -1 * specnumber(t(x)))
  
  return(data.frame(Symbiont = colnames(x), Structural.Specificity = structural))
  
}

structural.specificity(data = x, abundance.weighted = TRUE)