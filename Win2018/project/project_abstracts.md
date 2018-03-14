---
layout: page
mathjax: true
permalink: /Win2018/project/abstracts
---

## EE 372: Project Abstracts

#### Next-gen Basecallers
_Andre Cornman and Logan Spear_  

Our project was on basecallers for ONT data. When ONT reads were first released, basecalling was done based off of stochastic models involving multiple bases being inside of the pore during each measurement. As discussing in class, for this model, the raw data is first turned into levels, which are supposed to correspond to the presence of a particular context in the pore, and then bases are called using Viterbi. Modern basecallers have made the switch to directly using the raw signal and from using HMM models to Neural Networks (RNN's). Basecallers that employ RNN's typically see a greater accuracy than HMM based basecallers, although, as our results and results of others suggest, the RNN basecallers are very sensitive to the species of DNA that they are trained on.

#### Improved Long-Read Assembly from Hi-C information
_Mark Nishimura_  

Hi-C-seq is a useful tool for determining the chromosome's three-dimensional conformation in the nucleus. In this work, we attempt to use Hi-C information, in the form of a contact map, to achieve improved long-read assembly. By solving an LP relaxation of the Traveling Salesman Problem, HINGE + Hi-C provides better repeat resolution on simulated data.

#### Parallelizing Medoid Calculations through Spark
_Daniel Hsu and Brijesh Patel_  

We seek to find a representative cell from a single cell gene expression dataset from mice brain
cells with high dimensions. The Med-dit algorithm uses the concept of multi-armed bandits to find
medoids in almost linear time. We want to see the effect on runtime of this algorithm put in a Spark
framework and in conjunction with Google Cloud. Spark is an open source MapReduce framework for
distributed systems that is used in most modern systems for distributed computing. It provides a clean
API for managing datasets and running parallelizable algorithms, which is what we are hoping to utilize
for speedup in our algorithm.

#### Single Cell Co-Expression Subnetwork Detection
_Laura Shen and Po-Hsuan Wei_  

The transcriptome reflects the genes being expressed at any given time since it can vary depending on external environmental conditions. Many genes do not function in isolation, but work in a coordinated fashion with other genes to create complex proteins and other organic molecules required for biological function. Recent work has been done to assess different detection methods to find these gene networks with single-cell data. This project aims to find possible gene networks with single-cell data using different subnetwork detection methods, with corresponding metrics to evaluate the performances.

#### Limitations of Single Cell RNA-Seq in Neuronal Tissues: The Effect of Dissociation on Gene Expression

_Michelle Drews_

With the increasing feasibility and versatility of high throughput sequencing technologies,
single cell RNA sequencing has arisen as a very popular method in neuroscience for
interrogating the cellular compositions of different brain regions. However, a huge
potential limitation of single cell sequencing in tissues like the brain, is that these tissues do
not exist naturally as single cells, but instead exist as a highly interconnected and anchored
tissue. To make single cells, a fairly lengthy dissociation process must be undertaken,
where a digestive enzyme breaks down the solid tissue into single cells. This then raises
the question what effect does this dissociation process have on gene expression, and could
it be altering downstream clustering analyses. To address this question, we have taken
brains from five adult mice and immediately harvested RNA from the left half, while the
right half was dissociated into single cells like one would do for single cell RNAseq, and
RNA harvested from the single cells in bulk. These RNA were then sequenced by the Quake
lab, and raw read counts mapped to the _mus musculus_ RefSeq transcriptome were returned.
Differential expression analysis in DESeq2 revealed that, even with Bonferroni correction
at a threshold of 0.001, approximately 15% of the transcriptome was perturbed following
dissociation. Included among the upregulated genes were inflammatory and damage
response genes such as HSPa1a and KLF2, as well as neuronal immediate early genes such
as Fos and Jun. Included among the downregulated genes were synaptic organizing
proteins and synaptic transmission proteins such as CAMK2a and KCNAB2. When these
genes were removed from single cell clustering analysis, different clusters were obtained,
indicating that these genes may be driving single cell clustering – not due to inherent
biological variation, but rather as an artifact of the preparation to single cell sequencing.
