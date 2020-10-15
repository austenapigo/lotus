# lotus package

>This package is named `lotus` in recognition of Kamala (meaning 'lotus') Devi Harris, a longtime champion of LBGTQ+ rights and the first Black and Asian female nominee for Vice President of the United States on a major party ticket.

<img align="right" width="350" height="350" src="https://github.com/austenapigo/lotus/blob/master/figures/lotus.png">

### Why use `lotus`? 

High-throughput sequencing has accelerated the rate at which we can characterize host-associated microbiomes and sequencing-based datasets are ever-accumulating. Yet often times we have little to no information about the ecological or evolutionary context for many of these pairwise interactions. How should we try to make sense of these datasets rich with information about composition (taxonomy) and interaction frequency (read counts) but usually not much else? 

We propose that host specificity, or the degree to which microbial symbionts are restricted within a host community, provides a reasonable first step to generate ecological and evolutionary information about the symbiont niche. In our most recent pre-print [Apigo and Oono 2020](_____), (1) we explored how host specificity can provide insight to the cryptic relationships between plants and foliar fungal endophytes and (2) used these metrics to test hypotheses that sought to identify general mechanisms that explained host-specific interactions. We think these metrics can be leveraged in many different systems to help answer questions about host-symbiont ecology and evolution. 

**The purpose of `lotus` is to provide functions that:**
1. **Quantify host specificity metrics per symbiont**
2. **Relavitize these metrics to null host specificity models to account for variation in symbiont read abundance**
3. **Visualize null expectations to observed host specificities within occupancy-abundance models**

We've modified a framework proposed from Poulin 2011 to quantify host specificity that considers how host specificity can vary in a number of different ways: 
- **Structural specificity** quantifies the most fundamental perspective of host specificity, a symbiont's presence and evenness of abundance among hosts. For example, this metric asks how many hosts does a given symbiont occupy? 
- **Network specificity** quantifies the 'strength' of host-symbiont interactions by accounting for all potential hosts a symbiont could occupy in a host community. This metric asks how host-specific is a given symbiont when weighted by all potential interactions that could occur? 
- **Phylogenetic specificity** quantifies host specificity relative to the phylogenetic scale of the hosts in a community. This metric asks what is the phylogenetic breadth of the hosts that a given symbiont occupies? 
- **Beta-specificity** quantifies host specificity relative to host spatial or temporal scales. This metric asks how consistent is a given symbiont to a host over that host's geographical range (e.g., across quadrats or transects) or through time (e.g., sampling periods)? More generally, beta-specificity quantifies how host-specific a symbiont is to a given host group (species, genus, etc.).

**For each of these metrics, a narrower symbiont niche indicates higher host specificity.**

![Host Specificity Figure](https://github.com/austenapigo/lotus/blob/master/figures/specificity.png)

<sup>**Conceptual diagram of host specificity metrics represented by relationships between plants and endophytes.** Structural, phylogenetic and beta-specificity diagram depicting two endophytes (blue and grey shaded areas) occupying varying plant species (bryophyte, monilophyte, gymnosperm and angiosperm) across two sites (Sites A and B). Structural specificity (blue and grey shaded areas within solid lines) quantifies how endophytes vary in presence and evenness on hosts, the most fundamental feature of host-specific interactions. For example, Endophyte 1 (blue shade) occupies fewer plant species in Sites A and B (two plant species per site) relative to Endophyte 2 (four plant species in Site A and five hosts in Site B) and has higher structural specificity. Phylogenetic specificity (cladogram) quantifies host specificity in the context of plant evolutionary relatedness. Endophyte 1 occupies a phylogenetically narrower breadth of plant species in Sites A and B (gymnosperms only) relative to Endophyte 2 (angiosperm, bryophyte, gymnosperm and monilophyte) and has higher phylogenetic specificity. In contrast to structural and phylogenetic specificity, beta-specificity (arrow between Sites A and B) quantifies how endophytes vary in the consistency of their specificity to plant species by geography (Krasnov et al. 2011) or time (not depicted). Endophyte 1 occupies the same set of hosts in Sites A and B and has higher beta-specificity, whereas Endophyte 2 occupies a more variable set of plant species across Sites A and B.</sup>

### Null Models of Host Specificity 
Host specificity analyses could be particularly vulnerable to biased ecological inferences by not accounting for relationships between host specificity and endophyte abundance. For example, rarer symbionts will always have a lower theoretical maximum of the number of hosts they could occupy and are frequently categorized as more host-specific relative to more abundant symbionts. The relationship between a symbiont’s host specificity and its abundance is an extension of abundance-occupancy relationships, the widely documented biogeographic pattern in macroecology that predicts a positive relationship between a species abundance and its occurrence across sites (Gaston et al. 2000, Shade et al. 2018). By extending this to host-symbiont interactions, if more abundant symbionts tend to occupy more hosts, there is expected to be a negative relationship between an symbiont’s abundance and its host specificity, or occupancy across plants.

In our study, we found evidence of negative host specificity and symbiont abundance relationships. To account for bias that occurred between endophyte read abundance and observed host specificity, we relativized observed host specifities to null expectations (see Methods and Supplementary Figure 5 in our pre-print). **Doing so asked, how does a given symbiont's host specificity compare to the expected value for a symbiont of that read abundance under random community assembly?** Null models are often used in community ecology to test the null hypothesis that observed ecological structure is no different from a randomized structure (Gotelli and Graves 1996, Gotelli 2001, Ulrich et al. 2012). Here, the aim of null models is to shuffle symbiont read abundances among host samples to produce randomized distributions of host specificity. Observed host specificities are then relativized by the null expectation for a symbiont of that given read abundance to account for the propensity of rarer symbiont to occupy fewer hosts. There's a slight variation in this method for phylogenetic specificity. Basically, there can be convergence in this metric as the number of hosts a symbiont occupies approaches the total number of hosts in the community. We account for it using a slightly different model already implemented in the `picante` package (more details in the paper). 

![Null Model Figure](https://github.com/austenapigo/lotus/blob/master/figures/null_model.png)

<sup>**Conceptual diagram of host specificity metrics relativized by null models from Apigo and Oono 2020 represented by relationships between plants and endophytes.** (A-B) Each shape represents an individual endophyte. Shapes vertically positioned above each other have equal read abundance. Shapes positioned horizontal to each other occupy the same number of plants or have equal uncorrected host specificity values. The null expectation of the relationship between host specificity and read abundance is represented by a theoretical 1:1 relationship based on occupancy-abundance relationships (abundant endophytes likely occupy more plants). Asterisks represent the deviance in observed host specificity from a null model and are equal for all shapes. We describe three scenarios that highlight how null models can account for read abundance variation in host specificity analysis.</sup>

<sup>Scenario 1: Endophytes with different read abundances occupy the same number of plants and thus have the same uncorrected host specificities (green triangle vs. blue circle). The green triangle endophyte occurs in more plants than expected for a rarer endophyte, while the blue triangle occurs in fewer plants than expected for a more abundant endophyte. After relativizing their host specificities to the null model, the green triangle endophyte has lower host specificity and the blue circle has higher host specificity.</sup>

<sup>Scenario 2: Endophytes with the same read abundance occupy different numbers of plants and thus have different uncorrected host specificities (yellow square vs. green triangle or blue circle vs. orange diamond). The yellow square endophyte and the green triangle endophyte occupy either fewer or more plants than expected for endophytes with this read abundance. After relativizing their host specificities to the null model, differences in host specificity remain the same to uncorrected values because they exhibit the same magnitude of deviance from the null expectation.</sup>

<sup>Scenario 3: Endophytes with different read abundances occupy different numbers of plants and thus have different uncorrected host specificities (yellow square vs. orange diamond or yellow square vs. blue circle or green triangle vs. orange diamond). The yellow square endophyte occurs in fewer plants than expected for a rarer endophyte, while the orange diamond endophyte occurs in more plants than expected for a more abundant endophyte. After relativizing their host specificities to the null model, differences in their host specificities are less than their uncorrected host specificities but are equal in magnitude to Scenario 2. These endophytes display the same amount of deviance from a null expectation for endophytes of their given read abundance.</sup>

<sup>(B-D) Host specificity relativized by null models with empirical data. In this study, instead of a theoretical 1:1 null expectation, we calculated null models from randomized communities. We quantified how host specificity of endophytes within randomized communities varied as a function of read abundance. We then relavitized observed endophyte host specificity to null models by taking the difference between an observed endophyte’s host specificity to null expectations. The black line represents a null expectation generated from 100 community randomizations with an inset equation describing its relationship between uncorrected host specificity and endophyte abundance. Red points refer to individual endophytes that occur within a given plant sample, Galium angustifolium. The host specificity of its endophytes was calculated as Shannon’s H (structural specificity) across the entire plant community, plotted as a function of log-transformed absolute endophyte read abundance and relativized by the null model. For phylogenetic specificity, differences between observed and null host specificities were divided by the standard deviation of the null model to account for decreasing variance relationships (see Methods).</sup>

### How to install `lotus`
```
# Install lotus
devtools::install_github("austenapigo/lotus")
```

### R Functions
+ Structural Specificity
    + `structural.specificity`: calculates host richness (the number of hosts a symbionts occupies) or Shannon’s H diversity index (Shannon and Weaver 1948) to quantify symbiont presence and evenness among hosts uncorrected by symbiont read abundance.
    + `null.structural`: generates a null model by randomizing the community matrix and calculating structural specificity per symbiont within each randomized community.
    + `deviance.structural`: calculates the deviance in observed structural specificity from the null expectation per symbiont per sample.
    + `plot.structural`: plots null and observed structural specificities per symbiont as a function of symbiont read abundance.

+ Network Specificity
    + `network.specificity`: calculates the Resource Range Index or Paired Difference Index to quantify the 'strength' of host-symbiont interactions by accounting for all potential hosts a symbiont could occupy in a host community uncorrected by symbiont read abundance. 
    + `null.network`: generates a null model by randomizing the community matrix and calculating network specificity per symbiont within each randomized community
    + `deviance.network`: calculates the deviance in observed network specificity from the null expectation per symbiont per sample
    + `plot.network`: plots null and observed network specificities per symbiont as a function of symbiont read abundance.
    
+ Phylogenetic Specificity
    + `phylogenetic.specificity`: calculates as the Mean Pairwise Phylogenetic Distance (Webb 2000) to quantify symbiont presence and evenness as a function of host phylogenetic breadth uncorrected by symbiont read abundance. You must supply a phylogenetic distance matrix (see `cophenetic` in the `stats` package). 
    + `null.phylogenetic`: generates a null model by randomizing the community matrix and calculating phylogenetic specificity per symbiont within each randomized community
    + `deviance.phylogenetic`: calculates the deviance in observed phylogenetic specificity from the null expectation per symbiont per sample
    + `plot.phylogenetic`: plots null and observed phylogenetic specificities per symbiont as a function of symbiont read abundance.

+ Beta-Specificity
    + `beta.specificity`: calculates the Sørensen (Diserud and Odegaard 2007) or Morisita-Horn (Chao et al. 2008) Multiple-Assemblage Overlap Measure to quantify endophyte interaction consistency to a given host group (species, genus, family, etc.) across space or time uncorrected by symbiont read abundance.
    + `null.beta`: generates a null model by randomizing the community matrix and calculating beta-specificity per symbiont within each randomized community
    + `deviance.beta`: calculates the deviance in observed beta-specificity from the null expectation per symbiont per sample
    + `plot.beta`: plots null and observed beta-specificities per symbiont as a function of symbiont read abundance.
    
 ### Basic Workflow
 ```
 library(lotus)
 
 # Calculate uncorrected host specificity 
 hs.object <- structural.specificity(data = x, abundance.weighted = TRUE, trim = TRUE) 
 
 # Explore data and identify whether negative or variance-decreasing relationships exist between host specificity and symbiont read abundance
 # plot histogram
 pearson.cor <- cor.test(hs.object, colSums(quad.rarefied)) # correlation test
 # visualize host specificity - read abundance relationships 
 
 # Randomize your community matrix to generate a null model boundary
 null.structural.object <- null.structural(x, iterations = 10, abundance.weighted = TRUE, randomization.method = "shuffle.web", notify = TRUE)
 
 # Calculate the deviance of observed host specificity from the null boundary and get averages per host sample 
deviance.structural(data = x, randomized = null.structural.object, abundance.weighted = TRUE, trim = TRUE, notify = TRUE)
 
 # Plot null or observed host specificity as a function of log-symbiont read abundance
 plot.structural(data = x, randomized = null.structural.object, abundance.weighted = TRUE, trim = FALSE, notify = TRUE) 
```
#### Expected Output:


#### Important Notes:
1. The community matrix for analysis should be organized with hosts populating the rows, while symbionts populate the columns.
2. The example dataset used here is the same as in our pre-print. Hosts are labeled by their sampling origin with a period and number identifier after their species name (e.g., host_species_name.1). If you do not have a spatially-explicit sampling desing, you could simply label host conspecifics with similar identifiers (.1, .2, .3, etc.) to calculate beta-specificity. Please make sure taxonomic labels before the periods match. These labels can vary depending on your questions and taxonomic-level of interest. For example, you might be interested in the beta-specificity of symbionts to a taxonomic family of hosts (see The Discussion in our pre-print that addresses this) and could label hosts as Pinaceae.1, Pinaceae.2, Pinaceae.3, etc., for example. 
3. Host specificity is calculated for a symbiont across the entire host community. For example, a symbiont found in given host would be evaluated for its host specificity relative to all hosts that are present in a given dataset. 
4. More positive values always indicate a narrower symbiont niche and thus higher host specificity. We've negated (multipled by -1) structural and phylogenetic specificity to do this. 
In the case of 'deviance'-related functions, these per symbiont metrics are averaged per host sample and report standard error of the mean in host specificity. 
5. `null`- and `deviance`-related functions are calculated with log base 2 transformed read abundance. 
6. `null`-related functions accept `randomization.methods` from the bipartite::nullmodel() function (Dormann et al. 2009; Dormann 2011). If you wish to use other randomization functions, like vegan::permatswap(), save your output as a list and you can supply that it to any `deviance`-related function. 
7. `deviance`- and `plot`- related function use a second-order polynomial null boundary generated from the host specifities of randomized communities to calculate deviance. 
5. `lotus` functions have a 'trim' parameter within them to remove symbionts that are either singletons (a symbiont with a read abundance of one) or appear in only one host sample. This is a logical statement where TRUE removes these symbionts from analysis. We think they should be removed because it's unclear how one could relativize either of those scenarios to a null model. A symbiont with a read abundance of one or that appears in only one sample will always have the highest host specificity value. 
6. `lotus` functions have a 'notify' parameter within them to tell you which iteration the for() loop is on. 

### Why not use multivariate techniques to answer this question? 

We'd be appreciative of constructive feedback before we push this R package to CRAN. Please create an issue or feel free to email me directly (aapigo@ucsb.edu) if you have questions or suggestions.

#### References 
