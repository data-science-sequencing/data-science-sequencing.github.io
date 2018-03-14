---
layout: page
mathjax: true
permalink: /Win2018/project/guidelines
---

## EE 372 : Project Guidelines and Ideas

Projects for the course should be broadly related to the overarching theme of the class: sequencing. Projects can have a creative (or innovative component) or be expository. We expect expository projects also to have some component involving real data. Up to two students can collaborate on a project. Some datasets may be larger and therefore take longer to process, so start early! Please keep the following dates in mind:

1. Office hours to discuss projects  (Friday, 2 February 2018)
2. Proposal deadline (Friday, 9 February 2018)  
3. Milestone report and meeting (Friday, 2 March 2018)  
4. Poster presentations (Thursday, 15 March 2018)  
5. Project write-up deadline (Wednesday, 21 March 2018)

Students may also find [GEO](http://www.ncbi.nlm.nih.gov/geo/) useful in obtaining datasets.

### Deliverables

#### Office hours

Office hours will be held in Packard 264 (David's office) from 3:00-5:00pm. Please sign up for a 20-minute slot using [this Google Doc](https://docs.google.com/spreadsheets/d/1ySZWc7qNUYyR0Ayj3adITypp-Ea2u8L9vF9aRAEKbv8/edit?usp=sharing). The purpose of this office hours is for us to help shape your ideas into a good project. Before coming to the meeting, make sure you have done some initial research and have an idea of what problem you want to work on.

#### Proposal

The project proposal should be about a page in length (not including figures and references) and briefly explain the following:

1. Objective: Why is your problem interesting? What do you hope to accomplish?
2. Data: What is the dataset you will be working with?  
3. Methodology: What is your approach?  
4. Potential problems: What major obstacles do you foresee? What will you need to overcome these obstacles?

#### Milestone report and meeting

The purpose of the milestone meeting is to ensure that you have made progress on your project. During the meeting, we will discuss your milestone report, which should be about 2 pages in length and hit the following points:

1. Introduction: Brief summary of your proposal. What has changed?  
2. Initial results: What have you tried? What has not worked? What relevant figures or equations have you produced? What do they mean?
3. Next steps: How will you overcome initial obstacles? What else do you want to accomplish for this project?

Like the initial project office hours, the milestone meeting will also take place from 3:00-5:00pm in Packard 264 (David's office). Please sign up for a 20-minute slot using [this Google doc](https://docs.google.com/spreadsheets/d/1ySZWc7qNUYyR0Ayj3adITypp-Ea2u8L9vF9aRAEKbv8/edit?usp=sharing).

#### Poster Presentation

The poster presentation will take place in class on the last day of lecture. The exact location is to be determined and will be announced later. The poster should summarize your findings for the project and be 3' by 2' in size.

#### Writeup

The write-up should be about 4 pages in length (not including figures and references) and include the following sections:

1. Introduction: Background information about your project. What was the objective of this project? Why is your problem interesting?  
2. Data: What datasets did you look at for this project?  
3. Methodology: How did you go about accomplishing your objective?  
4. Results: What relevant figures or equations have you produced? What do they mean? What went wrong?  
5. Future Work: How do you build upon the work you have done?

## Project ideas

A project broadly related to any of the following topics will be appropriate for the course. Students can pick projects outside the ones listed below but should check suitability with the course staff.

1. **Cancer genomics:** Cancer has been labeled "a disease of the genome,"
and several unusual genomic events have been tied to cancer. These events include
 [chimeric reads](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3107329/) and
 [copy number aberrations](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2703871/).
 [The Cancer Genome Atlas](https://cancergenome.nih.gov/) offers a plethora of
data for exploring cancer genomics.
1. **Single cell clusteroids:** With single-cell data scaling up quickly, almost linear-time algorithms
for clustering the datasets become more valuable. Finding clusteroids
(see [this](http://infolab.stanford.edu/~ullman/mmds/ch7.pdf) write up or
[this](https://www.youtube.com/watch?v=YdqTScQFKQs) and
[this](https://www.youtube.com/watch?v=U3sdWVqMWEc) video for an overview),
or representatives of clusters, is one widely used approach for handling these large datasets.
Doing so if we have want to cluster cells as points in large dimensions can be done in
almost linear time as shown
 [here](https://arxiv.org/pdf/1711.00817.pdf) recently.
However if one wants reduce the dimension of the problem by say computing
a k-nearest neighbour graph and then find the clusteroid almost-linear
algorithms are not known. The best algorithm known is a $$O(n^{\frac{3}{2}})$$ algorithm
called [Trimed](https://arxiv.org/abs/1605.06950). Can one do better
in this case?
1. **Single cell cancer genomics:** One advantage of the resolution
 granted by single-cell assays is that we can observe how a cellular population
 differentiates. Looking at single-cell datasets for cancer cells may provide
 insight on how cancer evolves and which genes are important for this evolution.
 [This work](https://www.nature.com/articles/ncomms15081), for example, attempts
 to profile breast cancer using single-cell sequencing.
1. **Single cell assays:** Recent technologies have allowed researchers
 to extract cell signatures at the resolution of individual cells (e.g. single-cell RNA-seq).
  What can we learn about cell populations from these datasets? Prominent works in this area
  include [DropSeq](http://www.cell.com/abstract/S0092-8674(15)00549-8) and
  [10x](https://www.nature.com/articles/ncomms14049).
1. **Single-cell simulator:** A single-cell data simulator with known ground truth would be useful for evaluating single-cell methods (and understanding the mechanism underlying data generation). [Splatter](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1305-0) is a tool that simulates single-cell gene count data using a parametric model. Extending this to read count data or other modeling techniques could be useful.
1. **Doublet/multiplet detection:** Single-cell technology relies on isolating individual cells in droplets, but occasionally more than one cell will end up in a droplet. Can we detect these multiplet events, and can we further identify the types of cells in the multiplet? See [DropSeq](http://www.cell.com/abstract/S0092-8674(15)00549-8) and
[10x](https://www.nature.com/articles/ncomms14049) papers. When evaluating a new technology, gauging the doublet rate is important.
1. **Barcode indel correction:** In single-cell RNA-Seq experiments, cells are tagged with unique barcodes. Grouping reads based on their barcodes is an important computational step but is hindered by insertion and deletion events. [Here](https://www.biorxiv.org/content/biorxiv/early/2017/05/09/136242.full.pdf) is an example strategy for handling these barcode errors.
1. **Long read diploid assembly on read overlap graphs:** While long reads are a powerful
alternative to short-read sequencing, assembling genomes from these reads
come with their own challenges.
See [HINGE](https://www.ncbi.nlm.nih.gov/pubmed/28320918)
for an approach to long-read assembly to output a graph for the haploid case.
Can one improve it for the diploid case using the graph information.
1. **Gene networks:** Genes often act in conjunction with several other genes. With the advent of single-cell technology, we have enough data to estimate the pairwise interactions of different genes. A recent work in [Nature](https://www.nature.com/articles/s41598-017-15525-z) explores this. Additionally, comparing these networks for different datasets could yield insight into the underlying biology.
1. **Variant calling:** Diploid organisms such as humans have two chromosomes. Variant calling is the problem of figuring out the differences between the two chromosomes. A great tool for comparing variants is [PLINK](http://zzz.bwh.harvard.edu/plink/).
1. **Peak calling for epigenetics:**  While cells have roughly identical genomes, they express different parts of the genome due to epigenetic factors. Assays such as ChIP-seq and ATAC-seq help scientists study which parts of the genome are accessible for transcription. This problem involves determining which parts of the genome are accessible based on obtained reads. Two popular peak calling tools are [MACS](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3120977/) and [ZINBA](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2011-12-7-r67).  
1. **Short Tandem Repeat (STR) calling:**  STRs are short strings of DNA which appear one behind the other. These are known to be important for biological function, but calling them is challenging. A tool for performing STR calling is [LoBSTR](http://www.ncbi.nlm.nih.gov/pubmed/22522390).
1. **Dealing with batch effects:** While scientists are careful when executing experimental protocols, there will always be differences between repeated experiments. Batch effects are observed when the experiment ID (e.g. day the experiment was done) overcomes signal from the biological variable of interest. Tools like [ComBat](http://biostatistics.oxfordjournals.org/content/8/1/118.abstract) have been designed to mitigate batch effects.  
1. **Meta-genomics:** When one has the genomes of multiple organisms (or many genomes from the same type of organism), a natural question to ask is: how do these organisms relate to each other on a phylogenetic level? Some interesting work done here include [MetaGene](http://www.ncbi.nlm.nih.gov/pubmed/17028096) and [MetaBat](https://peerj.com/articles/1165/). [This](http://www.nature.com/nbt/journal/v32/n8/full/nbt.2939.html) is another interesting paper on this. This problem has become more relevant recently due to growing interest in microbiome research.
1. **Multi-omics:** With RNA-Seq, ChIP-Seq, Hi-C-Seq, etc., we are often inundated with many different categories of data for the same type of biological sample. Could we combine the information from these different types of data to discover something more? Check out [this review paper](https://www.nature.com/articles/nrg3868) to get started.
1. **Hi-C-Seq:** Hi-C-Seq is an assay used to study the 2D and 3D structure of the chromatin, giving insights on how the genome folds and how distant regions of the genome can interact with each other. This technology has been used for [assembling the genomes of polyploid organisms](https://www.nature.com/articles/nature22043), for example.
1. **Long range information:** In class we have covered short reads and long reads. There are other read types including read clouds ([10x Genomics](http://www.10xgenomics.com/)) and reads with long range information ([BioNano Genomics](http://www.bionanogenomics.com/)). From a theoretical perspective, what is possible and what is not with such reads?
1. **Pseudotime:** Single-cell sequencing can shed light on how cells evolve. The problem of rearranging individual cells into some sort of evolutionary trajectory (sometimes involving branching) has been tackled by tools such as [Monocle](http://cole-trapnell-lab.github.io/monocle-release/papers/) and [Wanderlust](https://www.c2b2.columbia.edu/danapeerlab/html/wanderlust.html).
1. **Base-calling for fourth generation sequencing technologies:** This is an area of active research with interesting recent papers involving methods such as [deep learning](http://arxiv.org/pdf/1603.09195.pdf).
1. **Reanalysis of published data:** Several datasets related to the above topics already exist. What can you discover by reanalyzing one of these datasets using a novel method?
