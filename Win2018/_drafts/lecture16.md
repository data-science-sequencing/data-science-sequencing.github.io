---
layout: page
mathjax: true
permalink: /Win2018/lectures/lecture16/
---
## Lecture 16: Single-cell RNA-seq - Introduction

Thursday 1 March 2018

_scribed by Po-Hsuan Wei and edited by the course staff_

-----------------

### Topics



1.	<a href='#bulk'>Limitations of bulk RNA-Seq</a>
1.	<a href='#early'>Early single-cell RNA-Seq</a>
1.	<a href='#modern'>Modern single-cell RNA-Seq</a>
    - <a href='#umi'>Unique molecular identifiers</a>
1.	<a href='#error'>Sources of error</a>
    - <a href='#errbar'>Barcoding errors</a>
    - <a href='#errencap'>Encapsulation errors</a>
1.	<a href='#conclusion'>Conclusion</a>

### <a id='bulk'></a>Limitations of bulk RNA-Seq

From a high level, for bulk RNA-Seq we extract the RNA content for several cells in a tissue by pooling together the cells together. After preparing the genomic content of the cells, we feed them through a sequencer. We obtain reads, perform alignment and quantification (e.g. via expectation-maximization), and ultimately obtain a vector of gene or transcript expressions describing the sample. Because we are sequencing all the cells in the sample, we lose the cell-level information and obtain some sort of average of all the cells.

[IMAGE of population density v. transcript_1 abundance with bulk measurement]

We give two examples for when bulk measurements fail:

1. Bulk measurements cannot distinguish between a change in regulation (i.e. difference in gene expression for a certain cell type in the tissue) and a change in composition (i.e. difference in proportions of cell types in the tissue).

2. If we are interested in how a cell differentiates from some base type (e.g. stem cell) to a specialized type (e.g. myoblast), bulk measurements will lose the composition of cells in different _stages_ of differentiation. Some cells differentiate faster than others, for example. The axis of differentiation has also been referred to as _pseudotime_. Waddington in 1957 visually drew an analogy between the cell differentiation process and a ball falling down an epigenetic "landscape:"

[IMAGE: epigenetic landscape]

### <a id='early'></a>Early single-cell RNA-Seq

The first study that succeeded in performing RNA-Seq for an individual cell was published in 2009 by [Tang et al.](https://www.ncbi.nlm.nih.gov/pubmed/19349980) A fundamental chemistry challenge here is the microscopic amount of starting material. One single cell, however, does not tell us much about the underlying biology.

[IMAGE of population density v. transcript_1 abundance with bulk measurement]

A key challenge here is that the cost associated with sequencing multiple cells is quite high in terms of both time and money if we need to address each cell in isolation. This problem is alleviated by _microfluidics_, which parallelize the sequencing of multiple cells. A microfluidics chip can isolate individual cells and compartamentalize the chemical reactions. One of the earlier chips was Fluidigm C1, which was designed by the company [Fluidigm](https://www.fluidigm.com/) and can sequence 96 cells in parallel:

[VIDEO: microfluidics chip: Fluidigm C1]

From a high level level, the single-cell RNA-Seq method can be described as follows:
1. Lysis: breaking the cells to release the genomic content
2. Reverse transcription: transforming RNA transcripts back to DNA, known as _cDNA_
3. Preamplification
4. Fragmentation: chopping up longer fragments into shorter ones more appropriate for sequnecing
5. Polymerase chain reaction ([PCR](https://en.wikipedia.org/wiki/Polymerase_chain_reaction)): to amplify the fragments
6. Sequencing

### <a id='modern'></a>Modern single-cell RNA_Seq

In 2015, researchers created methods for highly parallelizing the single-cell RNA-Seq process by isolating individual cells in droplets using microfluidics. The primary two methods are [DropSeq](http://mccarrolllab.com/dropseq/) and [inDrop](http://www.cell.com/cell/abstract/S0092-8674(15)00500-0). We will discuss DropSeq in more detail.

[IMAGE: DropSeq pipeline]

In each droplet, we also have a _bead_, and several short fragments of DNA are attached to the bead by design. All fragments attached the bead contain a short region that is identical, and this region is known as a _barcode_. Note that during the droplet generation step, we can end up with several droplets with either no cells or no beads, or droplets with more than one cells.

[VIDEO: Droplet generation]

After isolating cells within droplets, the fragments within cells will hybridize to the fragments attached to the bead. We generate hundreds of millions of _paired-end reads_ using modern sequencing technology.

[IMAGE: barcode-UMI-cDNA]

Note that the first read in each pair is attached the barcode, allowing us to identify which cell the read came from. The second read in each pair can be aligned to a transcriptome, allowing us to identify which transcript the read came from. We will discuss the _UMIs_ in more detail in the next section. More specifically, the computational workflow looks something like the following:

1. Group reads by the cells they come from using the barcodes
2. Align cDNA reads to identify transcripts
3. Count reads per cell per gene

[IMAGE: Computational workflow]

#### <a id='umi'></a>Unique molecular identifiers

Recall that each bead is designed such that several DNA fragments are attached to it, and these fragments are used for capturing the fragments from the lysed cell. Each of these bead fragments share the same cell barcode, but they also each contain a different random short sequence called a _unique molecular identifier_ (UMI). During the PCR step, we amplify each transcript several folds, but due to biases in the PCR process, some transcripts may certain PCR cycles, resulting in exponentially less copies of that transcripts than expected. More importantly, comparing the counts of different transcripts will be difficult now due to the PCR noise. UMIs allow us to know exactly how many pre-PCR molecules corresponding to each transcript there are. Using UMIs instead of reads, the computational workflow can be adjusted to:

1. Group reads by the cells they come from using the barcodes
2. Align cDNA reads to identify transcripts
3. Count _unique molecules_ per cell per gene

For DropSeq, the UMIs were all length 8, resulting in $$B = 2^{16} \approx 64000$$ possible barcodes. Each cell has roughly $$N = 100000$$ transcripts. It seems that we do not have enough combinations to distinguish all of the transcripts, since the probability of two transcripts having the same UMI is upper-bounded by $$N^2/2B$$. Luckily, we can leverage information from the second read. The barcodes are only necessary for distinguishing molecules _of the same gene type_. If we only look at molecules aligning to the same gene, then our $$N$$ becomes significantly smaller.

### <a id='error'></a>Sources of error

While the single-cell RNA-Seq process as presented so far seems straightforward, there are a few sources of noise than make the computational problem more interesting:

1. Barcodes can have errors
2. Errors in cell encapsulation

#### <a id='errbar'></a>Barcoding errors

One method for generating barcodes is known as _splitting and pooling_. In this process, we start with a collection of beads and four "buckets" of nucleotides, one for each nucleotide. We split the beads randomly into four parts, dip them into the corresponding buckets, re-mix them, and repeat. This will result in several random barcodes.

A fundamental problem here is identifying which barcodes are due to sequencing errors and which ones are actual barcodes. To solve this problem, we can simply count how many of each barcode we have. Let's say we have 1000 cells and 50000 reads per cell with a sequencing error rate of 20%. By plotting the frequency of each barcode, we can visually separate true barcodes (which appear about 80% of the time) from erroneous barcodes (which appear more less frequently). We can then fix the erroneous barcodes.

#### <a id='errencap'></a>Encapsulation errors

As the cells move through the microfluidic chip, there are a few factors we want to control. We can model the cell capture success rate using a Poisson process with some Poisson rate $$\lambda$$, which we tune using the concentration of cells in the solution we feed to the microfluidics. With a high Poisson rate (and therefore a high cell concentration), we would obtain more multiplets, less singlets, and less empty.

[IMAGE: cell encapsulation plot]

For example, one way DropSeq was evaluated involved mixing mouse cells with human cells to get a sense of the technology's doublet rate with respect to cell concentration.

[IMAGE: DropSeq doublet plots]

### <a id='conclusion'></a>Concluding remarks

While sequencing cost is going down quickly (faster than Moore's law), there is also a different curve of interest that's been growing exponentially: the size of single-cell RNA-Seq experiments (in terms of number of cells) has increased exponentially. Most recently, [1.3 million cells](https://community.10xgenomics.com/t5/10x-Blog/Our-1-3-million-single-cell-dataset-is-ready-to-download/ba-p/276) have been sequenced by 10X genomics. At approximately $0.00001 per read, this would cost on the order of half a million dollars. Sequencing cost is becoming the bottleneck again.

[IMAGE: Perspective curve]
