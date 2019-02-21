#  _doseR_ 
## A software package for dosage compensation analysis and much more...

*doseR* is bioinformatics software primarily designed for analysis of sex chromosome dosage compensation using RNA-seq data. However, the software package and underlying statistical framework can be applied broadly to detect shifts in gene expression among an arbitray number of pre-defined groups of loci.

The *doseR* package, implemented in R as part of the bioconductor project, is being developed in the [Walters Lab](http://walterslab.org/) at the University of Kansas. The project is funded by the National Science Foundation (ABI-1661454, beginning September, 2017). The work will be done in close collaboration with [Thomas Hardcastle](http://people.ds.cam.ac.uk/tjh48/) at the University of Cambridge. 

The software is currently in early stages of development and not yet ready for public release. Examples of the underlying principle can be found in the following publications:
* [Walters, J.R., Hardcastle, T.J. & Jiggins, C.D. (2015). Sex Chromosome Dosage Compensation in Heliconius Butterflies: Global yet Still Incomplete? Genome Biology and Evolution, 7, 2545–2559.](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4607515/)
* [Hardcastle, T.J. & Lewsey, M.G. (2016). Mobile small RNAs and their role in regulating cytosine methylation of DNA. RNA Biology, Vol 13, Issue 11.](http://www.tandfonline.com/doi/abs/10.1080/15476286.2016.1218591)


### NSF project summary

#### Overview
This research aims to develop and make broadly available a novel linear-modeling statistical methodology for analyzing sex chromosome dosage compensation using genome-wide RNA-seq expression data. Performance of the new statistical model and its software implementation relative to previous methods of assessing dosage compensation will be evaluated through extensive simulations of RNA-sequencing data. Application to specific empirical data sets relevant to dosage compensation will also be examined and evaluated. The software implementation, named doseR, will be written in the R statistical programming language and distributed as part of the Bioconductor suite of bioinformatic software tools. 

#### Intellectual merit
Genome-wide assessment of gene expression via RNA-sequencing is increasingly used to assess patterns of sex chromosome dosage compensation.  However, the statistical approaches currently employed for such analyses are far from ideal given the nature of the data and the desired set of inferences. Currently, biological replicates are averaged into a single measurement per gene and heavily normalized. Then particular effects of gene expression on the sex chromosome relative to autosomes are evaluated using absolute expression while gene dosage effects are assessed using expression ratios, in both cases using non-parametric statistical tests. These approaches have distinct shortcomings. With two independent tests, one using ratios and the other expression values, it is impossible to assess the relative magnitudes of a dosage versus sex chromosome effect on gene expression. Also, because replicates are averaged before analysis, no accounting is made for variation between replicates in assessing dosage compensation. These analyses thus neglect the natural variation in gene expression between individuals that may cause excessive confidence in the presence and estimation of dosage effect sizes.
A more statistically robust approach is to employ linear mixed-effects modeling of gene expression. This provides a unified statistical framework to assess magnitude and significance of both chromosome-specific and dosage effects on gene expression. Moreover, it is applied directly to the sequencing read counts for each gene and incorporates scaling factors such as sequencing depth and transcript length into the models describing the data, as is done in most analyses of differential expression. As such, extensive normalization is avoided and statistical replicates are readily incorporated into the analysis.  

Implementing this linear modeling in the doseR software will make this sophisticated statistical framework widely available to biologists researching dosage compensation. This is the first effort to create a bioinformatic tool specifically for dosage analysis and should substantially increase repeatability and comparability across analyses that are currently done primarily with custom scripting. While the immediate motivation for software development is dosage compensation analysis, the proposed methodology can be employed in any analytical scenario requiring the detection of directional shifts in expression for multiple, specific subsets of genes. It therefore provides a tool with broad utility in systems biology research.


#### Broader impacts
The grant will support the development of an eight-hour workshop providing an introductory overview of genome biology and bioinformatics methods for University of Kansas undergraduates. The aim of this workshop is to raise student awareness concerning relevant resources and research methods representative of genome biology and bioinformatics.  
