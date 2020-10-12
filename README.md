# Lotus Package {R}

<i> This packaged is named `lotus` in honor and recognition of the first Black woman nominee for Vice President on a major party ticket, Kamala (meaning 'lotus' in Sanskrit) Devi Harris.

High-throughput sequencing has accelerated the rate at which we can characterize host-associated microbial communities and sequencing-based datasets are ever accumulating. However, often times we have little to no information about the ecological or evolutionary context for many of these pairwise interactions. How should we try to make sense of these datasets rich with information about composition and interaction frequency, but (usually) without much else? 

In our most recent pre-print (Apigo and Oono 202X; _____), we quantified host specificity with a novel class of metrics to test hypotheses that integrate information about the endophyte lifestyle and features of the plant community that seek to understand general mechanisms that produce host-specific relationships. 

We propose that host specificity can help fill this gap. Host specificity, or the degree to which symbionts are restricted in a host community, can provide fundamental information about the symbiont niche. We've modified a framework proposed from Poulin 2011 that describes three types of host specificity: structural specificity, phylogenetic specificity and beta-specificity in the context of host-parasite systems. We've added another metric called network specificity.

  + Structural specificity (FFE abundance and evenness)
  + Network specificity (interaction strength)
  + Phylogenetic specificity (evolutionary relationships)
  + Beta-specificity (spatial or temporal turnover)
  
For each of these metrics, a narrower host breadth indicates higher host specificity. Structural specificity quantifies the most fundamental ‘dimension’ of host specificity, the sum and evenness of abundance among hosts (Poulin et al. 2011). Network specificity quantifies the strength of plant-FFE interactions by accounting for all potential hosts a FFE could occupy in a plant community. Phylogenetic specificity quantifies host specificity relative to the phylogenetic scale of the plant hosts in a community, or the mean phylogenetic distance among occupied hosts (Webb et al. 2008). Structural, network and phylogenetic specificity quantify the degree of host specificity within a single locality, termed alpha-specificity (Fig. 1; Poulin et al. 2011). Analogous to alpha diversity (Whittaker 1972), these three host specificity metrics do not account for spatiotemporal variation of the interaction. Beta-specificity, however, quantifies the degree to which a given FFE displays consistent host specificity across a range of contexts. Each specificity metric has a presence-absence and abundance-weighted variant with higher values indicating higher host specificity. 

Host specificity analyses may be particularly vulnerable to biased inferences because read abundant is often 

Important notes:
1. hosts are rows; symbionts are columns
2. This dataset has host plants as rows labeled by their sampling origin (quadrat, transect, etc.; .1 = quadrat 1) and fungal symbionts as columns. If you do have a spatially-explicit sampling desing, you could simply label host conspecifics (or another taxonomic level) with different identifiers (.1, .2, .3, etc.). 
3. singletons removed
4. symbionts that only appeared once in a sample were removed (noise parameter) 
5. For all metrics, more positive values = higher host specificity. 

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

 Multiple-assemblage overlap measures 
Described in Chao et al. 2008 and Jost et al. 

this function is useful if you are interested in quantifying beta-diversity when samples are grouped by multiple assemblages or sites (quadrats, transects, etc.). In this case, species (fungi) reside within samples (host plants) and a grouped by sampling units (quadrats). 


per symbiont per sample
 

