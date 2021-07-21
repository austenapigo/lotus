# lotus 0.9.0 (beta version)
> I'm be appreciative of constructive feedback before we push this R package to CRAN. Please create an issue on GitHub or email me (aapigo@ucsb.edu) if you have questions or suggestions.

> This package is named `lotus` in recognition of Kamala (meaning 'lotus') Devi Harris, a longtime champion for gender equality and the first Black and Asian female Vice President of the United States.

<img align="right" width="350" height="350" src="https://github.com/austenapigo/lotus/blob/master/figures/lotus.png">

### Why use `lotus` and what does it do?

Host-symbiont relationships are ubiquitous throughout natural systems. Yet often times we have little to no information about the ecological or evolutionary context for many of the pairwise interactions between hosts and symbionts. We propose that host specificity, or the degree to which symbionts are restricted within a host community, provides a reasonable first step to understand the ecological and evolutionary limits to  the symbiont niche.  We developed this R package because we think these metrics are useful, but underutilized, and can be leveraged in many different systems to help answer questions about host specificity. See Poulin et al. 2011 (reference below) to read about the original coneceptualization of these metrics. 

**The purpose of `lotus` is to provide functions that:**
1. **Quantify host specificity metrics per symbiont**
2. **Relavitize these metrics to null host specificity models to account for variation in symbiont abundance**
3. **Visualize null expectations to observed host specificities within occupancy-abundance models**

### What are the metrics? 
Host specificity exists on a continuum of phylogenetic, spatial and temporal scales. For example, a given symbiont could be restricted in the number of hosts it associates with whereas another symbiont could be restricted in the phylogenetic or geographic breadth of available hosts. We've further developed a framework originally proposed in host-parasite systems by Poulin et al. (2011) that considers how host specificity can vary in a number of different ways: 
- **Structural specificity** quantifies the most fundamental perspective of host specificity, a symbiont's presence and evenness of abundance among hosts. Structural specificity asks how many hosts does a given symbiont occupy? 
- **Network specificity** quantifies the 'strength' of host-symbiont interactions by accounting for all potential hosts a symbiont could occupy in a host community. Network specificity asks how host-specific is a given symbiont when weighted by all potential interactions that could occur? 
- **Phylogenetic specificity** quantifies host specificity relative to the phylogenetic scale of the hosts in a community. Phylogenetic specificity asks what is the phylogenetic breadth of the hosts that a given symbiont occupies? 
- **Beta-specificity** quantifies host specificity relative to host spatial or temporal scales. Beta-specificity asks how consistent is a given symbiont over a host's geographical range (e.g., across quadrats or sites) or through time (e.g., sampling periods)? 

**For each of these metrics, a narrower symbiont niche indicates higher host specificity.**

![Host Specificity Figure](https://github.com/austenapigo/lotus/blob/master/figures/specificity.png)

<sup>**Conceptual diagram of host specificity metrics represented by relationships between plants and endophytes. However, these metrics are useful for any host-symbiont of choice!** Structural, network, phylogenetic and beta-specificity diagram depicting two endophytes (blue and grey shaded areas) occupying varying plant species (bryophyte, monilophyte, gymnosperm and angiosperm) across two sites (Sites A and B). **Structural specificity** (blue and grey shaded areas within solid lines) quantifies how endophytes vary in presence and evenness on hosts, the most fundamental feature of host-specific interactions. For example, Endophyte 1 (blue shade) occupies fewer plant species in Sites A and B (two plant species per site) relative to Endophyte 2 (four plant species in Site A and five hosts in Site B) and has higher structural specificity. **Phylogenetic specificity** (cladogram) quantifies host specificity in the context of plant evolutionary relatedness. Endophyte 1 occupies a phylogenetically narrower breadth of plant species in Sites A and B (gymnosperms only) relative to Endophyte 2 (angiosperm, bryophyte, gymnosperm and monilophyte) and has higher phylogenetic specificity. **Network specificity** (dashed lines) quantifies the strength of host-specific interactions by accounting for all potential hosts a symbiont could occupy in a host community. Endophyte 1 occupies a lower proportion of plants in Sites A and B (2/4 in Site A and 2/5 in Site B) relative to Endophyte 2 (4/4 in Site A and 5/5 in Site B) and has higher network specificity. In contrast to structural and phylogenetic specificity, **beta-specificity** (arrow between Sites A and B) quantifies how endophytes vary in the consistency of their specificity to plant species by geography (Krasnov et al. 2011) or time (not depicted). Endophyte 1 occupies the same set of hosts in Sites A and B and has higher beta-specificity, whereas Endophyte 2 occupies a more variable set of plant species across Sites A and B.</sup>

### Why do we need null host specificity models? 
Host specificity analyses could be particularly vulnerable to biased inferences by not accounting for relationships between host specificity and symbiont abundance. Rarer symbionts will always have a lower theoretical maximum of the number of hosts they could occupy and are frequently categorized as more host-specific relative to more abundant symbionts. If one didn't account for variation in symbiont read abundance, rare symbionts could be systematically biased to be more host-specific simply because they have a lower probablity of occurring in many hosts relative to more abundant symbionts. The relationship between a symbiont’s host specificity and its abundance is an extension of abundance-occupancy relationships, the widely documented biogeographic pattern in macroecology that predicts a positive relationship between a species abundance and its occurrence across sites. If more abundant symbionts tend to occupy more hosts, there is expected to be a positive relationship between an symbiont’s abundance and its host specificity, or occupancy across plants.

To account for bias that often occurs between symbiont read abundance and observed host specificity, we standarized the effect sizes of observed host specifities to null distributions. **Relativizing host specificity to null models allows one to ask, how does an observed symbiont's host specificity compare to its expected value under random community assembly?** We propose that this allowed us to compare the host specifities of symbionts that vary widely in read abundance instead of splitting rare vs. abundant symbiont communities at an arbitrary relative abundance threshold or not accounting for read abundance at all. 

We term observed host specificity as “absolute host specificity” because these values are not modified by comparison to a null model. Absolute host specificities are then relativized by the null expectation for a symbiont to account for the propensity of rarer symbionts to occupy fewer host species that we term “relative host specificity”. Relative host specificity is calculated as the standardized effect size of a symbiont’s absolute host specificity (obs - mean_null_host_specificity / standard deviation of the null model) and is averaged symbionts per host. Each symbiont is only compared to itself among the randomized communities. A relative host specificity value greater than zero indicates that an symbiont was more host-specific relative to itself within the randomized communities. We propose that relative host specificities are a better method to classify host-specific interactions because they do so depending on their expectation to occur by chance rather than being biased towards classifying symbionts as host-specific simply because they are rare.</sup>

### How to install `lotus`
```
# Install lotus
devtools::install_github("austenapigo/lotus")
```

### Getting Started
see an example workflow here: [lotus vignette](https://github.com/austenapigo/lotus/blob/master/vignettes/lotus_vignette.Rmd)

### Why not use multivariate techniques to evaluate host specificity? 

Contemporary microbial ecology studies almost exclusively infer host specificity by the degree of compositional similarity of symbiont communities per host in multivariate ordination space. Similarity among symbiont communities within each host is ultimately influenced by the degree of niche overlap across all hosts per symbiont. For example, specialist symbionts that are consistently found within the same host species, but not others, should increase compositional similarity within the same host species and increase dissimilarity among different host species. Few studies utilize multivariate approaches in tandem with univariate approaches measured per symbiont (e.g., the number of host species a symbiont occupies), perhaps because they are assumed to be correlated, but whether it's the **best and only tool** for every question related to host specificity is still an open question. Metrics calculated per symbiont can tell you more about *why* clustering in ordination space occurs (or does not occur) and can give a more transparent view of how symbionts vary in their host specificity. 

### References 
+ Pedro J. Aphalo (2020). ggpmisc: Miscellaneous Extensions to 'ggplot2'. R package version 0.3.6. https://CRAN.R-project.org/package=ggpmisc
+ Dormann CF (2011) How to be a specialist? Quantifying specialisation in pollination networks. Netw Biol 1:1–20
+ Dormann CF, Gruber B, Fruend J (2008) Introducing the bipartite package: analysing ecological networks. R News 8:8–11
+ Dormann CF, Fründ J, Blüthgen N, Gruber B (2009) Indices, graphs and null models: analyzing bipartite ecological networks. Open Ecol J 2:7–24
+ Jari Oksanen, F. Guillaume Blanchet, Michael Friendly, Roeland Kindt, Pierre Legendre, Dan McGlinn, Peter R. Minchin, R. B. O'Hara, Gavin L. Simpson, Peter Solymos, M. Henry H. Stevens, Eduard Szoecs and Helene Wagner (2019). vegan: Community Ecology Package. R package version 2.5-6. https://CRAN.R-project.org/package=vegan
+ Poisot T, Canard E, Mouquet N, Hochberg ME (2012) A comparative study of ecological specialization estimators. Methods Ecol Evol 3:537–544
+ Poulin R, Krasnov BR, Mouillot D (2011) Host specificity in phylogenetic and geographic space. Trends Parasitol 27:355–361
+ Shannon CE, Weaver W (1948) The mathematical theory of communication. University of Illinois Press, Urbana
+ S.W. Kembel, P.D. Cowan, M.R. Helmus, W.K. Cornwell, H. Morlon, D.D. Ackerly, S.P. Blomberg, and C.O. Webb. 2010. Picante: R tools for integrating phylogenies and ecology. Bioinformatics 26:1463-1464.
+ Webb CO, Ackerly DD, McPeek MA, Donoghue MJ (2002) Phylogenies and community ecology. Ann Rev Ecol Syst 33:475–505
+ Webb CO, Ackerly DD, Kembel SW (2008) Phylocom: software for the analysis of phylogenetic community structure and trait evolution. Bioinformatics 24:2098–2100
