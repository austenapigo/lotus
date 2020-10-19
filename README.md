# lotus package (still under development as of 10/18/2020)

> This package is named `lotus` in recognition of Kamala (meaning 'lotus') Devi Harris, a longtime champion for gender equality and the first Black and Asian female nominee for Vice President of the United States on a major party ticket.

<img align="right" width="350" height="350" src="https://github.com/austenapigo/lotus/blob/master/figures/lotus.png">

### Why use `lotus` and what does it do?

High-throughput sequencing has accelerated the rate at which we can characterize host-associated microbiomes and sequencing-based datasets are ever-accumulating. Yet often times we have little to no information about the ecological or evolutionary context for many of the pairwise interactions between hosts and symbionts. How should we try to make sense of these datasets rich with information about composition (taxonomy) and interaction frequency (read counts) but usually not much else? 

We propose that host specificity, or the degree to which microbial symbionts are restricted within a host community, provides a reasonable first step to generate ecological and evolutionary information about the symbiont niche. In our most recent pre-print [Apigo and Oono 2020](_____), we explored how host specificity can provide insight to the cryptic relationships between plants and foliar fungal endophytes (fungi that live asymptomatically in plant photosynthetic tissues, primarily within leaves). We used novel metrics to test hypotheses that sought to identify general mechanisms that explained why host-specific interactions occur. We developed this package because we think these metrics are useful, but underutilized, and can be leveraged in many different systems to help answer questions about host specificity. 

**The purpose of `lotus` is to provide functions that:**
1. **Quantify host specificity metrics per symbiont**
2. **Relavitize these metrics to null host specificity models to account for variation in symbiont read abundance**
3. **Visualize null expectations to observed host specificities within occupancy-abundance models**

### What are the metrics? 
Host specificity exists on a continuum of phylogenetic, spatial and temporal scales. For example, a given symbiont could be restricted in the number of hosts they associate with whereas another symbiont could be restricted in the phylogenetic or geographic breadth of available hosts. We've further developed a framework originally proposed in host-parasite systems by Poulin et al. (2011) to quantify host specificity that considers how host specificity can vary in a number of different ways: 
- **Structural specificity** quantifies the most fundamental perspective of host specificity, a symbiont's presence and evenness of abundance among hosts. Structural specificity asks how many hosts does a given symbiont occupy? 
- **Network specificity** quantifies the 'strength' of host-symbiont interactions by accounting for all potential hosts a symbiont could occupy in a host community. Network specificity asks how host-specific is a given symbiont when weighted by all potential interactions that could occur? 
- **Phylogenetic specificity** quantifies host specificity relative to the phylogenetic scale of the hosts in a community. Phylogenetic specificity asks what is the phylogenetic breadth of the hosts that a given symbiont occupies? 
- **Beta-specificity** quantifies host specificity relative to host spatial or temporal scales. Beta-specificity asks how consistent is a given symbiont over a host's geographical range (e.g., across quadrats or sites) or through time (e.g., sampling periods)? 

**For each of these metrics, a narrower symbiont niche indicates higher host specificity.**

+ We originally described these metrics here: [Apigo and Oono 2018](https://link.springer.com/chapter/10.1007/978-3-319-89833-9_2).
+ We used these metrics to test hypotheses in a plant-endophyte system in our recent pre-print here: [Apigo and Oono 2020](_____). 

![Host Specificity Figure](https://github.com/austenapigo/lotus/blob/master/figures/specificity.png)

<sup>**Conceptual diagram of host specificity metrics represented by relationships between plants and endophytes.** Structural, network, phylogenetic and beta-specificity diagram depicting two endophytes (blue and grey shaded areas) occupying varying plant species (bryophyte, monilophyte, gymnosperm and angiosperm) across two sites (Sites A and B). **Structural specificity** (blue and grey shaded areas within solid lines) quantifies how endophytes vary in presence and evenness on hosts, the most fundamental feature of host-specific interactions. For example, Endophyte 1 (blue shade) occupies fewer plant species in Sites A and B (two plant species per site) relative to Endophyte 2 (four plant species in Site A and five hosts in Site B) and has higher structural specificity. **Phylogenetic specificity** (cladogram) quantifies host specificity in the context of plant evolutionary relatedness. Endophyte 1 occupies a phylogenetically narrower breadth of plant species in Sites A and B (gymnosperms only) relative to Endophyte 2 (angiosperm, bryophyte, gymnosperm and monilophyte) and has higher phylogenetic specificity. **Network specificity** (dashed lines) quantifies the strength of host-speciifc interactions by accounting for all potential hosts a symbiont could occupy in a host community. Endophyte 1 occupies a lower proportion of plants in Sites A and B (2/4 in Site A and 2/5 in Site B) relative to Endophyte 2 (4/4 in Site A and 5/5 in Site B) and has higher network specificity. In contrast to structural and phylogenetic specificity, **beta-specificity** (arrow between Sites A and B) quantifies how endophytes vary in the consistency of their specificity to plant species by geography (Krasnov et al. 2011) or time (not depicted). Endophyte 1 occupies the same set of hosts in Sites A and B and has higher beta-specificity, whereas Endophyte 2 occupies a more variable set of plant species across Sites A and B.</sup>

### Why do we need null host specificity models? 
Host specificity analyses could be particularly vulnerable to biased ecological inferences by not accounting for relationships between host specificity and symbiont abundance. Rarer symbionts will always have a lower theoretical maximum of the number of hosts they could occupy and are frequently categorized as more host-specific relative to more abundant symbionts. For example, is a symbiont with a read abundance of 10 that appears in one host ***really*** equally as host-specific to a symbiont with a read abundance of 10,000 that appears in only one host? We'd argue that the more abundant symbiont in this hypothetical scenario is more host-specific because it appears in far fewer hosts than expected for an symbiont with such high read abundance.  

If one didn't account for variation in symbiont read abundance, rare symbionts could be systematically biased to be more host-specific simply because they have a lower probablity of occurring in many hosts relative to more abundant endophytes. The relationship between a symbiont’s host specificity and its abundance is an extension of abundance-occupancy relationships, the widely documented biogeographic pattern in macroecology that predicts a positive relationship between a species abundance and its occurrence across sites. If more abundant symbionts tend to occupy more hosts, there is expected to be a negative relationship between an symbiont’s abundance and its host specificity, or occupancy across plants.

In our study, we found evidence of negative relationships between host specificity and symbiont abundance. To account for bias that occurred between endophyte read abundance and observed host specificity, we relativized observed host specifities to null expectations (see Methods and Supplementary Figure 5 in our pre-print). **Normalizing host specificity to null models allowed us to ask, how does an observed symbiont's host specificity compare to its expected value under random community assembly?** We propose that this allowed us to compare the host specifities of endophytes that varied widely in read abundance instead of splitting rare vs. abundant endophyte communities at an arbitrary relative abundance threshold (common in other studies) or not accounting for read abundance at all. There's a slight modification of this null model approach for phylogenetic specificity because relationships between phylogenetic specificity and read abundance are usually decreasing-variance (rather than negative) and you can read more about the specifics in our paper. 

![Null Model Figure](https://github.com/austenapigo/lotus/blob/master/figures/null_model.png)

<sup>**Conceptual diagram of host specificity metrics relativized by null models from Apigo and Oono 2020 represented by relationships between plants and endophytes.** (A-B) Each shape represents an individual endophyte. Shapes vertically positioned above each other have equal read abundance. Shapes positioned horizontal to each other occupy the same number of plants or have equal uncorrected host specificity values. The null expectation of the relationship between host specificity and read abundance is represented by a theoretical 1:1 relationship based on occupancy-abundance relationships (abundant endophytes likely occupy more plants). Asterisks represent the deviance in observed host specificity from a null model and are equal for all shapes. We describe three scenarios that highlight how null models can account for read abundance variation in host specificity analysis.</sup>

<sup>Scenario 1: Endophytes with different read abundances occupy the same number of plants and thus have the same uncorrected host specificities (green triangle vs. blue circle). The green triangle endophyte occurs in more plants than expected for a rarer endophyte, while the blue triangle occurs in fewer plants than expected for a more abundant endophyte. After relativizing their host specificities to the null model, the green triangle endophyte has lower host specificity and the blue circle has higher host specificity.</sup>

<sup>Scenario 2: Endophytes with the same read abundance occupy different numbers of plants and thus have different uncorrected host specificities (yellow square vs. green triangle or blue circle vs. orange diamond). The yellow square endophyte and the green triangle endophyte occupy either fewer or more plants than expected for endophytes with this read abundance. After relativizing their host specificities to the null model, differences in host specificity remain the same to uncorrected values because they exhibit the same magnitude of deviance from the null expectation.</sup>

<sup>Scenario 3: Endophytes with different read abundances occupy different numbers of plants and thus have different uncorrected host specificities (yellow square vs. orange diamond or yellow square vs. blue circle or green triangle vs. orange diamond). The yellow square endophyte occurs in fewer plants than expected for a rarer endophyte, while the orange diamond endophyte occurs in more plants than expected for a more abundant endophyte. After relativizing their host specificities to the null model, differences in their host specificities are less than their uncorrected host specificities but are equal in magnitude to Scenario 2. These endophytes display the same amount of deviance from a null expectation for endophytes of their given read abundance.</sup>

<sup>(B-D) Host specificity relativized by null models with empirical data. In this study, instead of a theoretical 1:1 null expectation, we calculated null models from randomized communities. We quantified how host specificity of endophytes within randomized communities varied as a function of read abundance. We then relavitized observed endophyte host specificity to null models by taking the difference between an observed endophyte’s host specificity to null expectations. The black line represents a null expectation generated from 100 community randomizations with an inset equation describing its relationship between uncorrected host specificity and endophyte abundance. Red points refer to individual endophytes that occur within a given plant sample, *Galium angustifolium*. The host specificity of its endophytes was calculated as Shannon’s H (structural specificity) across the entire plant community, plotted as a function of log-transformed absolute endophyte read abundance and relativized by the null model. For phylogenetic specificity, differences between observed and null host specificities were divided by the standard deviation of the null model to account for decreasing-variance relationships.</sup>

### How to install `lotus`
```
# Install lotus
devtools::install_github("austenapigo/lotus")
```
Download process should look something like this:
![Download](https://github.com/austenapigo/lotus/blob/master/figures/download.png)

### R Functions
+ Structural Specificity
    + `structural.specificity`: calculates host richness (the number of hosts a symbionts occupies) or Shannon’s H diversity index (Shannon and Weaver 1948) to quantify symbiont presence and evenness among hosts uncorrected by symbiont read abundance.
    + `null.structural`: generates a null occupancy-abundance model by randomizing the community matrix and calculating structural specificity per symbiont within each randomized community.
    + `deviance.structural`: calculates and plots the deviance in observed structural specificity from the null expectation per symbiont per sample.

+ Network Specificity
    + `network.specificity`: calculates the Resource Range Index or Paired Difference Index (Poisot et al. 2012) to quantify the 'strength' of host-symbiont interactions by accounting for all potential hosts a symbiont could occupy in a host community uncorrected by symbiont read abundance. 
    + `null.network`: generates a null occupancy-abundance model by randomizing the community matrix and calculating network specificity per symbiont within each randomized community
    + `deviance.network`: calculates and plots the deviance in observed network specificity from the null expectation per symbiont per sample
    
+ Phylogenetic Specificity
    + `phylogenetic.specificity`: calculates as the Mean Pairwise Phylogenetic Distance (Webb 2000) to quantify symbiont presence and evenness as a function of host phylogenetic breadth uncorrected by symbiont read abundance. You must supply a phylogenetic distance matrix (see `cophenetic` in the `stats` package). 
    + `deviance.phylogenetic`: calculates and plots the deviance in observed phylogenetic specificity from the null expectation per symbiont per sample
        + Note: `null.phylogenetic` is combined within deviance.phylogenetic. The `picante` package conveniently calculates null phylogenetic models and calculates deviance (i.e., standardized effect sizes).

+ Beta-Specificity
    + `beta.specificity`: calculates the Sørensen (Diserud and Odegaard 2007) or Morisita-Horn (Chao et al. 2008) Multiple-Assemblage Overlap Measure to quantify endophyte interaction consistency to a given host group (species, genus, family, etc.) across space or time uncorrected by symbiont read abundance.
    + `null.beta`: generates a null occupancy-abundance model by randomizing the community matrix and calculating beta-specificity per symbiont within each randomized community
    + `deviance.beta`: calculates and plots the deviance in observed beta-specificity from the null expectation per symbiont per sample
    
 ### Example Workflow
 ```
# Load lotus
library(lotus)

# You can read more about each lotus function with the help function
help("structural.specificity")

# Calculate uncorrected host specificity (not relavitized to a null model)
hs.object <- structural.specificity(quad.rarefied, abundance.weighted = TRUE, trim = TRUE)
hs.object

# Explore data and evaluate relationships between host specificity and symbiont read abundance
plot(density(hs.object$Structural.Specificity)) # plot histogram

read.abund <- as.data.frame(colSums(quad.rarefied)) # get read abundances per symbiont
read.abund.trim <- read.abund[rownames(read.abund) %in% rownames(hs.object), ] # trim relative to hs.object

cor.test(hs.object$Structural.Specificity, read.abund.trim) # correlation test

plot(y = hs.object$Structural.Specificity, x = log(read.abund.trim), ylab = "Uncorrected Structural Specificity (HostRichness)", xlab = "Log Symbiont Read Abundance") # visualize host specificity - read abundance relationships
abline(lm(hs.object$Structural.Specificity~log(read.abund.trim)), col = "red")

# Randomize community matrix to generate a null model for deviance calculations
null.structural.object <- null.structural(quad.rarefied, iterations = 100, abundance.weighted = TRUE, randomization.method = "shuffle.web", trim = TRUE, notify = TRUE)

# Calculate and plot the deviance of observed host specificity from the null boundary and get averages per host sample
structural.dev <- deviance.structural(quad.rarefied, randomized = null.structural.object, abundance.weighted = TRUE, trim = TRUE, notify = TRUE)
structural.dev[[1]] # View data frame of output
structural.dev[[2]] # View occupancy-abundance model for the first sample
structural.dev[[81]] # View occupancy-abundance model for the last sample
 ```
#### Expected Output:
+ Output from structural.specificity()
![Output1](https://github.com/austenapigo/lotus/blob/master/figures/structural.output.png)

+ Output from cor.test()
![Output2](https://github.com/austenapigo/lotus/blob/master/figures/correlation_graph.png)

+ Output from null.structural()
![Output3](https://github.com/austenapigo/lotus/blob/master/figures/null.output.png)

+ Output from structural.dev()
![Output4](https://github.com/austenapigo/lotus/blob/master/figures/deviance.output.png)

+ Output from structural.dev()
![Output5](https://github.com/austenapigo/lotus/blob/master/figures/deviance.png)

#### Important Notes:
+ The community matrix for analysis should be organized with hosts populating the rows, while symbionts populate the columns. Please make sure your input data is a data frame with `is.data.frame()`. 

+ The example dataset used here is the same as in our pre-print. Hosts are labeled by their sampling origin with a period and number identifier after their species name (e.g., host_species_name.1). This format is required to calculated beta-specificity. 

+ Host specificity is calculated for a symbiont across the entire host community. For example, a symbiont found in given host would be evaluated for its host specificity relative to all hosts that are present in a given dataset. 

+ More positive values always indicate a narrower symbiont niche and thus higher host specificity. We've negated (multipled by -1) structural and phylogenetic specificity to make this consistent across all metrics. 

+ For beta-specificity, labels can vary depending on your questions and taxonomic-level of interest. For example, you might be interested in the beta-specificity of symbionts to a taxonomic family of hosts, rather than host species (see the Discussion in our pre-print that addresses this), and could label hosts as Pinaceae.1, Pinaceae.2, Pinaceae.3, etc., for example. 

+ Common `lotus` arguments: 
    + abundance.weighted Logical. TRUE calculates abundance-weighted metrics per symbiont. FALSE calculates presence-absence metrics per symbiont. 
    + model Character. Specify whether the null expectation should be approximated as a first-(linear) or second-(quadratic) order function. 
    + trim Logical. TRUE removes symbionts that occupy one host sample. FALSE keeps all symbionts. Note: We think they should be removed because host specificity from an observation of one will always result in the highest host specificity value and cannot be relevatized in a meaningful way to a null expectation. 
    + notify Logical. TRUE prints the current iteration of the for loop. 

### Why not use multivariate techniques to evaluate host specificity? 

You might find that it's common to see host specificity evaluted with multivariate methods that quantify and visualize differences in compositional dissimilarity (e.g., Bray-Curtis) among symbiont communities in multivariate space. Clustering among host samples from the same species or host taxonomic group in ordination space suggests they share similar communities of symbionts and this is usually interpreted as host specificity. However, we pose that there are many reasons host samples could be clustered.  

For example, if all symbionts in a given group of host samples were only found in those hosts, but not any others, we would expect these samples to cluster in ordination space and have high structural specificity (i.e., a low number of occupied hosts per symbiont across the entire host community). As another scenario, these same host samples could harbor many rare symbiont generalists found in other samples, but one very abundant specialist symbiont only found in this particular group of hosts. In this case, we might still see clustering in ordination space (if using an abundance-weighted metric like Bray-Curtis) but symbionts within these hosts would have lower structural specificity (relative to the previous scenario) because of all the rare generalists they harbor occupy more hosts. 

These two (of many possible) hypothetical scenarios highlight how metrics calculated per symbiont can tell you more about *why* clustering in ordination space occurs (or does not occur) and can give a more transparent view of how symbiont communities vary in their host specificity beyond differences in a given dissimilarity metric. We still think multivariate methods are extremely useful tools to identify host specificity, but whether it's the **best** tool for every question related to host specificity is still an open question. We advocate for future studies to consider using the alternative metrics we propose in conjunction with multivariate techniques that each provide useful perspectives to host specificity in host-associated microbiomes.

### Contact 

We'd be appreciative of constructive feedback before we push this R package to CRAN. Please create an issue on GitHub or email me (aapigo@ucsb.edu) if you have questions or suggestions.

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
