# Host Specificity Scripts

High-throughput sequencing has accelerated the rate at which we can characterize host-associated microbial communities and sequencing-based datasets are ever accumulating. However, often times we have little to no information about the ecological or evolutionary context for many of these pairwise interactions. How should we try to make sense of these datasets rich with information about composition and interaction frequency, but (usually) without much else? 

We propose that host specificity can help fill this gap. Host specificity, or the degree to which symbionts are restricted in a host community, can provide fundamental information about the symbiont niche. Moreover, it can be leveraged to test hypotheses that seek to identify the underlying mechanisms that produce turnover in symbiont communities. 

In our most recent pre-print (Apigo and Oono 202X; _____), we quantified host specificity with a novel class of metrics to test hypotheses that integrate information about the endophyte lifestyle and features of the plant community that seek to understand general mechanisms that produce host-specific relationships. 

Host specificity analyses may be particularly vulnerable to biased inferences because read abundant is often 

Important notes:
1. hosts are rows; symbionts are columns
2. hostspecies.quadrat
3. singletons removed
4. symbionts that only appeared once in a sample were removed (noise parameter) 

Host Specificity Metrics 
1. Structural Specificity
 + structural.specificity
 + null.structural
 + deviance.structural
 + plot.structural 

2. Phylogenetic Specificity
 + phylogenetic.specificity
 + null.phylogenetic
 + deviance.phylogenetic
 + plot.phylogenetic 

3. Beta-Specificity
 + beta.specificity
 + null.beta
 + deviance.beta
 + plot.beta 



 
 

