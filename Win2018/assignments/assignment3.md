---
layout: page
mathjax: true
permalink: /Win2018/assignments/assignment3/
---
## Assignment 3

posted Tuesday 27 February 2018  
due Tuesday 12 March 2018 at 11:59pm

**Submission policy**: Report all plots and your code in [this Jupyter notebook](/Win2018/assets/assignment3/ee372_assignment3.ipynb). Print your notebook as a PDF and attach it to the rest of your assignment. Turn in your assignment through [Gradescope](https://gradescope.com/).

### Question I: Minhashing

In class we discussed briefly how minhashing can be used for aligning reads with large error rates. In this question, we will explore the minhashing concept.

1. We can describe a read as a set of unique overlapping $$k$$-mers, and we would expect similar reads to have similar sets. Write a function that takes $$k$$ and a read as inputs and outputs a dictionary indicating which of the $$4^k$$ $$k$$-mers are in the read.

2. If we think of each read as a set of $$k$$-mers, a natural metric of similarity between two reads $$R_1$$ and $$R_2$$ is the Jaccard similarity, which is defined as
\\[
J(R_1,R_2) = \frac{|R_1 \cap R_2|}{|R_1|+|R_2|-|R_1 \cap R_2|}.
\\]
Explain how this metric captures similarity between two sets and how you might use this metric to align reads (2-3 sentences). Compute the minimum and maximum possible Jaccard similarity between any two sets.

3. Write a function to compute the Jaccard similarity between two dictionaries outputted by your function from part 1. Using this function and the one you wrote for part 1, compute the Jaccard similarity between the reads <tt>CATGGACCGACCAG</tt> and <tt>GCAGTACCGATCGT</tt> for $$k = 3$$. What is the runtime complexity of your function? If you have $$N$$ reads of length $$L$$ each, what is the worst-case runtime complexity for computing the Jaccard similarity between every pair of reads?

4. Suppose you have a function that hashes a $$k$$-mer to a value between $$1$$ and $$4^k$$. For minhashing, you would use this hash function to map each unique $$k$$-mer in a read to an index, ultimately returning the smallest index. Prove that the probability that two sets will generate the same minhash index is equal to their Jaccard similarity.

5. In practice, we would use multiple $$H$$ hash functions to compute multiple minhash indices for $$R_1$$ and $$R_2$$. Write down an estimator that uses the $$H$$ minhash indices to estimate $$J(R_1,R_2)$$. What is the runtime complexity of obtaining this estimation? How does this compare to the runtime you obtained for part 3?

### Question II: Haplotype phasing coverage

In this problem we examine a simplified version of the haplotype
assembly problem. We assume that a genome has $$n$$ SNPs, and each
SNP is heterozygous. Every mate-pair read covers a pair of
adjacent SNPs.
Download the Adjacency matrix of the SNPs from a 5000 base pair
region
(SNP 10000-1500) of chromosome 16 of NA12878
[here](/Win2018/assets/assignment3/matrix_sparse.pkl).
Use [pickle](https://wiki.python.org/moin/UsingPickle)
to load the matrix. The matrix is a sparse format
([scipy.sparce.csc_matrix](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csc_matrix.html)).
One can
convert a sparse matrix $$\texttt{A}$$ to dense by using
$$\texttt{A.to_dense()}\ $$.  


1. How many measurements are in the matrix?

2. What is the maximum number of measurements between two SNPs?

2. Is the graph connected?

3. What would
you expect the maximum number of measurments  to be if
you had  an Erdos-Renyi
graph with the same number of edges?

3. Download the ground truth of SNPs from
[here](/Win2018/assets/assignment3/ground_truth.pkl). Get the primary
eigenvector of the window of the first 200 SNPs in the
adjacency matrix and find its correlation coefficient with the ground
truth.


### Question III: RNA-seq Quantification

In class we discussed the basic EM algorithm for RNA-seq quantification in the simple case when the transcript lengths are all the same and there are no errors in the read. In this question, we will consider extensions to unequal transcript lengths and read errors. We start with the same RNA-seq model as discussed in class.

1. Implement the EM algorithm for the error-free case where all transcripts have the same length. Fill in the code in the last cell of the [Jupyter notebook](/Win2018/assets/assignment2/ee372_assignment3.ipynb).

2. Instead of equal transcript length $$\ell$$, let us now consider the case when the transcript lengths are $$\ell_1, \ell_2, \dots, \ell_K$$. The reads are still error-free.
- Develop the log likelihood model
- Derive the EM iterative algorithm for this model, specializing from the general EM discussed in the lecture.

3. Now suppose the reads have errors: each base is read correctly with probability $$1-\delta$$ and incorrectly with probability $$\delta$$, and if incorrect the base can be read equally likely as any of the other three possibilities.
- Generalize the log likelihood model in Part 1 to this case.
- Derive the EM iterative algorithm for this model, again specializing from the general EM algorithm.
- Suppose the alignment tool at your disposal can compute all exact alignments and all approximate alignments up to one base different. If the error rate $$\delta$$ is small such that the chance of a read having two or more errors is negligible, explain how you would use your alignment tool in implementing the EM algorithm derived in the previous part.

### Question IV: Single-cell RNA-seq

In this problem, we will analyze the single-cell RNA-Seq dataset published by Zeisel et al. in 2015. Download the dataset [here](/Win2018/assets/assignment3/Zeisel_data.zip).

1. You can find the paper describing their experiment [here](http://science.sciencemag.org/content/347/6226/1138). Answer the following questions:
- What kind of cells were they looking at? (and what organism did they come from?)
- How many cells were there?
- What is a unique molecular identifier and how was it used?
- How did they validate their findings?
- What was their motivation for performing the experiment?
- What did they discover?

2. Plot two histograms: log(total counts) across cells and log(total counts) across genes. Why do we need to take a log here? What does this say about the data? Remove the columns of the matrix corresponding to genes that are expressed less than 25 times across all cells. How many genes are removed?

3. Transform the data by adding 1 to all entries and taking the log (this is common practice). Two common dimensionality reduction strategies are principal component analysis (PCA) and t-distributed stochastic neighbor embedding (tSNE). You can implement PCA using the Python [sklearn package](http://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html). We recommend [this](https://github.com/danielfrg/tsne) tSNE package for generating a tSNE embedding. List a pro and a con of each of the two methods. Create the following three plots and report what you see (for all plots, color points by their label):
- Run PCA and plot the first two principal components
- Run tSNE and plot the first two tSNE components (Note: this may take a few minutes)
- Run tSNE on the top 50 principal components and plot the first two tSNE components

4. Using sklearn, we will test four different [clustering algorithms](http://scikit-learn.org/stable/modules/clustering.html) on the dataset and evaluate their performance using the labels provided by the authors. The four clustering methods are K-means clustering, spectral clustering, affinity propagation, and density-based spatial clustering of applications with noise (DBSCAN). Using the top 50 principal components for the log-transformed data, for each clustering method:
- Describe the hyperparameters one needs to tune
- Tune the hyperparameters until you obtain reasonable labels (i.e. you obtain about 9 clusters, just like the authors did)
- Report your final hyperparameters
- Compare your labels to the true labels using sklearn's [adjusted Rand Index](http://scikit-learn.org/stable/modules/generated/sklearn.metrics.adjusted_rand_score.html).

5. In a few sentences, describe how you might go about determining which genes most effectively distinguish a cluster from the other clusters (this is also known as _differential expression analysis_).
