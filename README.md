# lotus package {R}

>This packaged is named `lotus` in recognition of the first Black and Asian woman to be nominated for Vice President of the United States on a major party ticket, Kamala (meaning 'lotus') Devi Harris. 

### Introduction 

High-throughput sequencing has accelerated the rate at which we can characterize host-associated microbiomes and sequencing-based datasets are ever-accumulating. However, often times we have little to no information about the ecological or evolutionary context for many of these pairwise interactions. How should we try to make sense of these datasets rich with information about composition (i.e., taxonomy) and interaction frequency (i.e., read counts) but (usually) not much else? 

We propose that host specificity, or the degree to which microbial symbionts are restricted within a host community, provides a reasonable first step to generate ecological and evolutionary information about the symbiont niche. In our most recent pre-print [Apigo and Oono 202X](_____), we explore how host specificity can provide insight to the relationships between plants and their cryptic fungal symbionts and use these metrics to test hypotheses that seek to identify general mechanisms that can predict host-specific interactions. 

**The purpose of this package is to provide functions that:**
1. **Quantify host specificity metrics per symbiont**
2. **Relavitize these metrics to null host specificity models to account for variation in symbiont read abundance**
4. **Quantify the deviance in observed to null expectations in host specificity**
3. **Visualize null models to observed host specificities**

We've modified a framework proposed from Poulin 2011 to quantify host specificity that consider how host specificity can vary in a number of different ways: 
- ***Structural specificity*** quantifies the most fundamental ‘dimension’ of host speci- ficity, the sum and evenness of abundance among hosts. This metric asks how many hosts does a given symbiont occupy? 
- ***Network specificity*** quantifies the 'strength' of host-symbiont interactions by accounting for all potential hosts a symbiont could occupy in a host community. This metric asks how host-specific is a given symbiont when weighted by all potential interactions that could occur? 
- ***Phylogenetic specificity*** quantifies host specificity relative to the phylogenetic scale of the hosts in a community. This metric asks what is the phylogenetic breadth of the hosts that a given symbiont occupies? 
- ***Beta-specificity*** quantifies host specificity relative to host spatial or temporal scales. This metric asks: How consistent is a given symbiont to a host over that host's geographical range (e.g., across quadrats or transects) or through time (e.g., sampling periods)? More generally, beta-specificity cquantifies how host-specific a symbiont is to a given host species.  

***For each of these metrics, a narrower symbiont niche indicates higher host specificity.*** 

![Image of Yaktocat](https://github.com/austenapigo/lotus/blob/master/specificity.png)


Host specificity analyses could be particularly vulnerable to biased ecological inferences by not accounting for relationships between host specificity and endophyte abundance. For example, rarer symbionts will always have a lower theoretical maximum of the number of plants they could occupy and are frequently categorized as more host-specific relative to more abundant symbionts. The relationship between an symbionts’s host specificity and its abundance is an extension of abundance-occupancy relationships, the widely documented biogeographic pattern in macroecology that predicts a positive relationship between a species abundance and its occurrence across sites (Gaston et al. 2000, Shade et al. 2018). If more abundant symbionts tend to occupy more hosts, there is expected to be a negative relationship between an symbiont’s abundance and its host specificity, or occupancy across plants.

In our study, we found evidence for negative host specificity and symbiont abundance relationships. We compared empirical and null endophyte host specificities to account for bias that occurs between endophyte read abundance and observed host specificity (see Supplementary Figure 5 in our pre-print). Null models are often used in community ecology to test the null hypothesis that observed ecological structure is no different from a randomized structure (Gotelli and Graves 1996, Gotelli 2001, Ulrich et al. 2012). Here, the aim of null models is to shuffle endophyte read abundances among host samples to produce randomized distributions of host specificity. Observed host specificities were then relativized by the null expectation for an endophyte of that given read abundance to account for the propensity of rarer endophytes to occupy fewer plants (see Figure 1 and methods in our pre-print for more detailed explanations). 

```
# Install lotus
devtools::install_github("austenapigo/lotus")
```

#### Important notes:
1. The community matrix is organized with hosts populating the rows, while symbionts populate the columns.
2. The example dataset used here is the same as in our pre-print. Hosts are labeled by their sampling origin with a period and number identifier after their species name (e.g., host_species_name.1). If you do have a spatially-explicit sampling desing, you could simply label host conspecifics with these identifiers (.1, .2, .3, etc.) to calculate beta-specificity. Please make sure taxonomic labels before the periods match. These labels can vary depending on your questions and taxonomic-level of interest. For example, you might be interested in the beta-specificity of symbionts to a taxonomic family of hosts (see discussion in our pre-print that addresses this). 
3. These functions have a 'trim' parameter within them to remove symbionts that are either singletons (a symbiont with a read abundance of one) or appear in only one host sample. This is a logical statement where TRUE removes these symbionts from analysis. We think they should be removed because it's unclear how one could relativize either of those scenarios to a null model. A symbiont with a read abundance of one or that appears in only one sample will always have the highest host specificity value.  
4. More positive values always indicate higher host specificity and a narrower symbiont niche. We've negated (multipled by -1) structural and phylogenetic specificity to do this. 
6. Metrics are calculated per symbiont. In the case of 'deviance'-related functions, these per symbiont metrics are averaged per host sample. 
7. Host specificity is calculated for a symbiont across the entire community. For example, a symbiont found in given host would be evaluated for its host specificity relative to all hosts that are present in a given dataset. 

### Host Specificity Metrics 
+ Structural Specificity
    + `structural.specificity`: calculates 
    + `null.structural`
    + `deviance.structural`
    + `plot.structural` 

+ Phylogenetic Specificity
    + `phylogenetic.specificity`
    + `null.phylogenetic`
    + `deviance.phylogenetic`
    + `plot.phylogenetic`

+ Beta-Specificity
    + beta.specificity
    + null.beta
    + deviance.beta
    + plot.beta 
    
We think the metrics we use can be leveraged across many systems and have advantages over commonly-used multivariate techniques (e.g., ordination). An alternative approach to multivariate methods that quantify host specificity per endophyte community are univariate metrics that quantify host specificity per endophyte species or molecular taxonomic unit (Apigo and Oono 2018). Multivariate methods may unintentionally neglect different types of host specificity that have long been recognized to exist on a continuum of phylogenetic, spatial and temporal scales (Poulin and Mouillot 2003, Krasnov et al. 2011, Poulin et al. 2011). For example, an individual endophyte could be restricted in the number of plants they associate with whereas others could be restricted in the phylogenetic or geographic breadth of available plants (Supplementary Figure 3). 

#### Have Questions, Comments or Suggestions? 
Please create an issue or feel free to email me directly (aapigo@ucsb.edu). I would like this package to be useful to others, so please let me know what worked well for you or what didn't. 

#### References 
 Multiple-assemblage overlap measures 
Described in Chao et al. 2008 and Jost et al. 

 Structural specificity was calculated as host richness (the number of plants an endophyte occupies) or Shannon’s H diversity index (Shannon and Weaver 1948) to account for endophyte presence and evenness with the ‘diversity’ function in the vegan package (Oksanen et al. 2019). Phylogenetic specificity was calculated as the Mean Pairwise Phylogenetic Distance (Webb 2000) to quantify endophyte presence and evenness as a function of plant phylogenetic breadth with the ‘mpd’ function in the picante package (Kembel et al. 2010). An ultrametric phylogenetic tree was pruned from a backbone phylogeny of 74 533 vascular plant species representing all extant vascular plant families in North America with the ‘phylo.maker’ function in the V.PhyloMaker package (Qian and Jin 2016; Supplementary Figure 4). Beta-specificity was calculated as the Sørensen (Diserud and Odegaard 2007) or Morisita-Horn (Chao et al. 2008) Multiple-Assemblage Overlap Measure to quantify endophyte interaction consistency across five sampling quadrats with an R script that the authors wrote (see Data Accessibility section). Structural and phylogenetic specificity were negated (i.e., multiplied by negative one) such that a narrower endophyte niche among plants yielded a more positive host specificity value. 

We quantified structural and beta-specificity to test whether abundant plant species harbored more host-specific endophytes that either occupied fewer plants (structural specificity) or were more consistent in their interactions and occupied the same plant species across the landscape (beta-specificity). We quantified phylogenetic specificity to understand how endophytes varied in the evolutionary breadth of plants they occupied and performed analyses at varying plant phylogenetic scales (plant species, plant family, plant clade). 


 
Rare endophytes have been shown to display unique ecological signatures distinct from abundant taxa (Oono et al. 2017, Zhang et al. 2018b). To account for this, studies typically separate rare and abundant taxa at a variety of relative abundance thresholds (e.g., 0.1% - 1%; Logares et al. 2014, Zhang et al. 2018a, Xue et al. 2018) to analyze these groups separately. Instead of relative abundance cut-offs to partition endophyte communities, 

First, we regressed observed endophyte structural, phylogenetic and beta-specificity as a function of log-transformed endophyte read abundance and calculated Pearson r correlation coefficients (Supplementary Figure 5) to assess whether observed endophyte communities displayed negative endophyte abundance-host specificity relationships (i.e., analogous to positive abundance-occupancy relationships). Second, we randomized communities and calculated each host specificity metric per endophyte ASV within each community randomization. Host specificities per endophyte ASV were regressed as a function of log-transformed endophyte abundance to define a null model expectation to compare observed host specificities with. Lastly, we relativized host specificity values per endophyte ASV by calculating the difference between an endophyte’s observed host specificity to the null model expectation for an endophyte of that given read abundance within an endophyte abundance-host specificity model (further description in Figure 1). We also tested whether the averaged host specificity value of the observed community was statistically different from the mean host specificity of randomized communities with a one-sample t-test (Supplementary Figure 6) to confirm that the observed community deviated from random structure and host-specific patterns did not occur due to random chance. In addition to endophyte singletons, we excluded any endophyte that only appeared in one sample from our analyses, regardless of read abundance, because host specificity from an observation of one will always result in the highest host specificity value and cannot be relevatized in a meaningful way to a null expectation. 

For structural and beta-specificity, the plant-endophyte community was randomized (n = 100 randomizations) with the ‘shuffle. web’ algorithm using the ‘nullmodel’ function in the bipartite package (Dormann et al. 2009; Dormann 2011). This randomization method redistributed all abundance data among plant samples and endophyte ASVs, thereby changing marginal totals but maintained connectance such that the number of plants with and without endophytes remained constant. The average relationship between host specificity and log-transformed endophyte read abundance among the 100 randomizations was used to generate a null model expectation per host specificity metric. We chose a second-order polynomial function to define this relationship to better represent non-linear relationships that were observed between null host specificities and log-transformed endophyte read abundance. We then calculated the deviance along the y-axis from an endophyte’s observed host specificity to the null model expectation and averaged this value across endophytes as the mean deviance from the null model per plant sample. A mean deviance value greater than zero indicates that, on average, a community of endophytes was more host-specific relative to endophytes with the same read abundances within randomized communities (Figure 1). We repeated analyses with the ‘r2dtable’ randomization method that redistributed abundance data and kept marginal totals constant and the ‘mgen’ randomization method that kept neither marginal totals or connectance constant. 

For phylogenetic specificity, we used null models to account for a different type of known bias where variance in phylogenetic specificity decreases with the number of plants an endophyte occupies (Supplementary Figure 7), which is correlated with endophyte abundance (Supplementary Figure 5). As the number of plants an endophyte occupies approaches the total number of plants in the community, the range of possible phylogenetic specificity values converges to a single value, or the mean pairwise phylogenetic distance of the entire plant community, producing a ‘funnel-shaped’ relationship. To account for this, we calculated standardized effect sizes of mean pairwise phylogenetic distances that quantify the difference between observed phylogenetic specificity and the mean value of phylogenetic specificity from randomized communities. These differences were then divided by the standard deviation of the null distribution of phylogenetic specificity to account this variance-reduction relationship. We calculated standardized effect sizes of Mean Pairwise Phylogenetic Distances (MPD) per endophyte ASV with the ‘ses.mpd’ function in the picante package (Kembel et al. 2010). We used the ‘taxa.labels’ null model to randomize the labels of the plant species phylogenetic distance matrix. We repeated this analysis with the ‘phylogeny.pool’ null model that randomized the phylogenetic distance matrix labels by drawing plant species with equal probability. We specifically chose null models that randomized the phylogenetic structure of the distance matrix because endophytes that occupy plant lineages from more basal nodes (e.g., ferns; one species) within the plant phylogeny might be more likely to have inflated pairwise phylogenetic distances relative to endophytes that occur in plant lineages from more derived nodes (e.g., angiosperms; 33 species). For example, if an endophyte appeared in two plants, one of which was a fern, it might be more likely to appear in a second plant lineage more distantly-related to a fern because angiosperms have much higher species richness in this plant community. We repeated phylogenetic specificity analyses at the plant clade, plant family and plant species level because inference could vary depending on the phylogenetic scale of the analysis. We grouped angiosperms into clades based on nomenclature from Jansen et al. (2007) as asterids (14 plant species), rosids (12 plant species), commelinids (five plant species) and basal angiosperms (two plant species) along with conifers as Pinophyta (four plant species) and ferns as Polypodiophyta (one plant species). We define ‘basal angiosperms’ as extant lineages of angiosperms branching from nodes more basal within the angiosperm phylogeny that are sister-taxa to the asterid, rosid and commelinid clade.
