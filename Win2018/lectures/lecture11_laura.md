---
layout: page
mathjax: true
---
## Lecture 11: Haplotype Phasing & RNA-Seq

Tuesday 13 February 2018

-----------------

## Topics
1. Spectral algorithm for phasing
2. RNA-seq

### Spectral algorithm for phasing
The problem of phasing involing determining the pattern of 0's and 1's. This can be reformulated as a community recovery problem where the community identity is detmined by [1-0] or [0-1] sequence. The data relationship comes from the linking reads where having the same parity indicates same community while having a different parity means the read are not likely in the same community. For example, the middle [0-1] linking on Chromosome 1 indicates that these SNPs are not in the same community. However, noisy data can cause errors in correctly identifying communities.

<div class="fig figcenter fighighlight">
  <img src="/Win2018/assets/lecture11/commrec.png" width="80%">
	<div class="figcaption">Haplotype phasing as a community recovery problem.</div>
</div>

Computating all possible partitions and finding the links is exponential in the number of nodes. Thus, we would like to take advantage of global information, but we need a more efficient method of computation. The adjacency matrix $$A$$ is an equivalent representation of the data where each column and row represents a SNP. This matrix is usually sparse (i.e. many entries are zero). The above example can be represented by the following adjacency matrix:

$$
A =
\begin{bmatrix}
\star & -1 & 0 & 0 \\
-1 & \star & 1 & -1 \\
0 & 1 & \star & 1 \\
0 & -1 & 1 & \star
\end{bmatrix}.
$$

For each entry of the adjacency matrix:
<ul>
  <li>$$A_{ij}=-1$$: different parity</li>
  <li>$$A_{ij}=0$$: no read</li>
  <li>$$A_{ij}=1$$: same parity</li>
</ul>

In the uniform linking model, $$q$$ is the probability of $$A_{ij}$$ being non-zero and $$p$$ is the probability of error:
<ul>
  <li>$$Pr(A_{ij}=0)=(1-q)</li>
  <li>$$Pr(A_{ij}\neq0)=q</li>
  <ul>
    <li>Pr(no error)=(1-p)</li>
    <ul>
      <li>$$A_{ij}=1$$ if $$(i,j)$$ have same parity</li>
      <li>$$A_{ij}=-1$$ if $$(i,j)$$ have different parity</li>
    </ul>
    <li>Pr(error)=p</li>
    <ul>
      <li>$$A_{ij}=1$$ if $$(i,j)$$ have different parity</li>
      <li>$$A_{ij}=-1$$ if $$(i,j)$$ have same parity</li>
    </ul>
  </ul>
</ul>

#### Efficient approximate recovery
In the ideal situation of no noise and every link is sampled, we can get a rank1 matrix with a unique eigen vector. The sampling process and errors are similar to noise and add randomness which allows us to compute the expectation of the adjacency matrix.

<div class="fig figcenter fighighlight">
  <img src="Win2018/assets/lecture11/expectation_adj_matrix.PNG" width="80%">
	<div class="figcaption">Expectation adjacency matrix with no errors</div>
</div>

This becomes a random matrix fluctuating about the mean matrix. For each random matrix, the principal eigen vector (inner product) is computed. We can then take the principal eigen vector (usually a non-integer value) and round to either +1 or -1). This is the basis of the spectral algorithm where we take the adjacency matrix, compute the principal eigen vector and threshold each element to determine the community. With more samples, the expectation of the eigen vector of the random matrix becomes closer and closer to the community vector. Even with a sparse matrix, spectral is quite accurate.

#### Genie-aided lower bound
We want to compute a lower bound on the number of reads required to produce a correct result. Suppose a genie told us the correct community of all SNPs except one (the blue node):

<div class="fig figcenter fighighlight">
  <img src="Win2018/assets/lecture11/genie_phasing.PNG" width="80%">
	<div class="figcaption">Genie tells us communities of all SNPs expect one.</div>
</div>

It is still possible to make a mistake if the majority of the linking reads are wrong. A lower bound such that this event happens with very lower probability:
$$\text{# of linking reads} > \frac{n \log n}{2[1-e^{-D(0.5\|p)}]}.$$

#### Two-step algorithm
<b>Step 1:</b> approximate recovery using spectral algorithm
<ul>
  <li>A good percentage of nodes are correct, but there could still be error, thus the need for step 2.</li>
</ul>
<b>Step 2:</b> refinement using majority vote for each SNP
<ul>
  <li>Assuming node neighbours are correct, we can do a majority vote and switch current community identity value if different from the majority. This can be done for all elements (nodes) in linear time.</li>
</ul>

The lower bound can be achieved by this two-step algorithm. This two-step algorithm was applied to real linking data obtained from 10X company where the reads are partitioned into clouds of reads (same side of chromosomes) of roughly few hundred SNPs and the random expectation matrices (also known as contact maps) were computed. Shown below are the contact maps for chromosome 16 where white = link and black = no link. The contact map shows the distributioni of "contact"/linking information.

<div class="fig figcenter fighighlight">
  <img src="/Win2018/assets/lecture11/contactmap.png" width="90%">
	<div class="figcaption">The contact map obtained on real 10X data in comparison to the contact map obtained from data sampled from the uniform linking model.</div>
</div>

The left picture is the actual contact map produced by the two-step algorithm from real linking data while the right picture is the theoretical uniform linking model. Thus, theory does not line up with practice. There are several reasons for this discrepancy:
<ol>
  <li>3k SNPs is ~3million genomic locations apart (3million segments)</li>
  <li>10X reads are 100k apart (read cloud)</li>
  <li>All linking information is diagonal centric while uniform linking model is more spread out => community recovery</li>
  <li>There is locality in communities (i.e. social graphs with friends)</li>
</ol>

If we zoom in closer to diagonal of the left picture contact map for the real linking data, we see that the zoomed contact map more closely resembled local uniform linking model. Computing the principal eigen vector by sliding the window along the diagonal and performing thresholding leads to much better results. However, we want community of all nodes, not just in the diagonal. How do we bridge this gap?

<div class="fig figcenter fighighlight">
  <img src="/Win2018/assets/lecture11/stitch.png" width="80%">
	<div class="figcaption">Motivation for a sliding window/stitching spectral algorithm (top right: recovered eigenvector, bottom right: thresholded recovered eigenvector).</div>
</div>

#### Spectral stitching algorithm
We now introduce a new algorithm called "spectral stitching":
<b>Step 1:</b> break adjacency up into many sub-matrices and perform approximate recovery on overlapping blocks via the spectral algorithm
	<ul>
    <li>Slide window over by half the window size to interleave the blocks</li>
  </ul>
<b>Step 2:</b> stitching across the blocks
	<ul>
	  <li>two community assignments (one from top and one from bottom) => might need to switch communities</li>
  </ul>
<b>Step 3:</b> perform refinement by taking majority vote on each SNP


This algorithm was evaulated on 10X real linking data using 4 metrics:
<ul>
  <li># of unphased SNPs: Some areas of chromosome have very few SNPs, so the reads are not long enough to pass the deserts meaning some SNPs are impossible to phase because they are so isolated.</li>
  <li>N50 of phased blocks</li>
  <li>short switch error rate</li>
  <li>long switch error rate</li>
</ul>

##### N50 of phased blocks
The output of phasing algorithm is typically in terms of the phase blocks (organized from longest to shorted). N50 is length of block such that that half of chromosome is contained.
<div class="fig figcenter fighighlight">
  <img src="/Win2018/assets/lecture11/n50.png" width="70%">
	<div class="figcaption">Visual definition of the N50 phasing metric.</div>
</div>

##### Switch errors within phased block
There are two flavours of switch errors: short and long. Short switch errors are the occasional error which a pair of SNPs have the incorrect community assignments. Long switch error are long sequences of errors due to incorrect community (reverses everything) which is very different from 5 independently isolated switches.
<div class="fig figcenter fighighlight">
  <img src="/Win2018/assets/lecture11/switching.png" width="70%">
	<div class="figcaption">Switch errors that can happen within a phased block.</div>
</div>

##### Comparison
For unphased SNPs, the spectral stitching algorithm designed to use data efficiently, so the 10X data is well covered. What happens with lower coverage? As coverage reduces, the unphased SNPs increases but spectral stitching is still better. Spectral stitching is more stable with long switch errors.
<div class="fig figcenter fighighlight">
  <img src="/Win2018/assets/lecture11/stitchresult.png" width="90%">
	<div class="figcaption">Performance of spectral-stitching algorithm on 10X data in comparison to 10X's algorithm.</div>
</div>

The time to perform spectral stitching for a few million SNPs is ~5500s on entire dataset (22 chromosomes) due to linear time.

### RNA-Seq
Unlike DNA which is fixed for a cell, RNA varies due to external environmental factors and thus is a dynamic indicator of how the cell is doing. RNA can be analyzed similar to the genome, but the focus on determining the abundance of the transcripts.
<div class="fig figcenter fighighlight">
  <img src="/Win2018/assets/lecture11/Figure1.png" width="80%">
	<div class="figcaption">The central dogma of molecular biology.</div>
</div>

Given reads of RNA, there are three possible scenarios:
<ol>
  <li>The transcriptome is known, but the abundances is not known.</li>
  <li>The genome is known, but the transcripts and their abundanes are not known.</li>
  <li>Neither the genome, transcripts, nor their abundances are known.</li>
</ol>

Reads could come from overlapping transcripts (multiple mappings). If we used unique reads only then we would waste a lot of potentially useful data and transcripts of consistings of entire reads would not be unique and thus be discarded. Better statistical methods for inference will be examined in the next lecture.
