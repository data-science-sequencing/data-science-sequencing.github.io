---
layout: page
mathjax: true
permalink: /Win2018/lectures/lecture18_jz/
---
## Lecture 18: Empirical Bayes and Principal Component Analysis

Thursday 8 March 2018

*scribed by Brijesh Patel and edited by the course staff*

-----------------

### Topics

1. <a href='#eb'>Empirical Bayes</a>
2. <a href='#pca'>Principal component analysis</a>

### <a id='eb'></a>Empirical Bayes

Recall that for differential expression analysis, we need more than one replicate in each condition. Otherwise, we would not be able to estimate variance or standard error (i.e. we would not be able to draw error bars). We obtain estimates for the mean and variances under the two different conditions A and B:

$$
\begin{align*}
\hat{\mu}^A & = \frac{1}{n}\sum_{i=1}^n X_i^A \\
\hat{\mu}^B & = \frac{1}{n}\sum_{i=1}^n X_i^B \\
(\hat{\sigma}^A)^2 & = \frac{1}{n-1} \sum_{i=1}^n (X_i^A - \hat{\mu}^A )^2 \\
(\hat{\sigma}^B)^2 & = \frac{1}{n-1} \sum_{i=1}^n (X_i^B - \hat{\mu}^B )^2.
\end{align*}
$$

Because getting replicates is expensive, we are interested in understanding the minimum number of replicates required such that we can draw reasonable conclusions. From a statistical testing point of view, our null hypothesis is that the two means and variances are equal:

$$
\begin{align*}
H_0: & X_i^A \sim N(\mu, \sigma^2) \\
& X_i^B \sim N(\mu, \sigma^2)
\end{align*}
$$

Please see [lecture 15](/Win2018/lectures/lecture15/) for a more detailed discussion on using the Student's $$t$$-test for differential expression analysis. We start by comparing the $$t$$-distribution to the standard normal distribution $$N(0, 1)$$:

[IMAGE: Gaussian v. t-distribution]

Note that the standard normal has smaller tails. This means that controlling the false positive rate for the $$t$$-distribution is more difficult, and we must use a larger cutoff in order to control the false positive rate. If we knew the variance, then the distribution of the $$T$$ statistic would be Gaussian (a linear combination of Gaussians is Gaussian):

$$ \tilde{T} = \frac{\hat{\mu}^A - \hat{\mu}^B}{\sigma}. $$

This highlights how the estimate of the variance dictates how good our test would be. Estimating this quantity introduces a lot of uncertainty in our distribution (hence fattening the tails). As shown in the figure above, having only 3 replicates does not give us a reliable estimate of the variance, resulting in a low-power test.

Empirical Bayes tells us that we can gain power by estimating some statistics associated all the data (i.e. by estimating some global prior using the data). [DESeq](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-10-r106) leverages the information across genes to power a statistical test for differential expression. With 20000 genes, for example, we can plot the variance with respect to the mean. While the estimate for an individual gene may be poor, we can fit some (often parametric) model for all genes. This process is also known as _shrinkage_ because we shrink the points towards our model, resulting in a tightening of the distribution. As shown in the figure, a Poisson model (blue) underestimates the variance in real data while a negative binomial distribution (orange) tends to fit the data better.

[IMAGE: DESeq figure]

### <a id='pca'></a>Principal component analysis

For the last part of the course, we will discuss various computational problems that arise from single-cell RNA-Seq analysis. A [recent dataset](https://community.10xgenomics.com/t5/10x-Blog/Our-1-3-million-single-cell-dataset-is-ready-to-download/ba-p/276) published by 10x genomics consists of 1.3 millions cells and 28000 genes. Each cell is represented by a 28000-dimensional expression vectors.

[IMAGE: 1.3 million dataset tsne]

A typical pipeline would generate the above image as an output. After sequencing all these cells, we would like gain a sense of how many cell types there are. In contrast, bulk RNA-Seq gives us an average of all the cells in the sample. For the above experiment, the result is a clustering of individual cells into 16 cell types.

We note that several single-cell analysis pipelines exist because the technology is still relatively young and hence evolving. An example single-cell analysis pipeline is:

1. Keep high-variance genes, reducing the dimension (e.g. going from 28000 genes to 2000 genes)
2. PCA to further reduce the dimension
3. Clustering the principal components

Principal component analysis (PCA) consists of the following steps:

1. Estimate the mean vector and covariance matrix ($$K$$ will be size 2000-by-2000 for the above example)
$$
\begin{align*}
\hat{\mu} & = \frac{1}{n}\sum_{i=1}^n x_i \\
\hat{K} & = \frac{1}{n} \sum_{i=1}^n (x_i-\hat{\mu}) (x_i-\hat{\mu})^T
\end{align*}
$$
2. Find $$v_1, \dots, v_d$$, the principal eigenvectors of $$\hat{K}$$, giving us the matrix $$V = [v_1 \ \dots \ v_d]$$. As a point of reference, $$d = 50$$ for the 10x dataset, and $$d = 32$$ for the DropSeq dataset.
3. Project the original data onto the principal eigenvectors: $$x_i \rightarrow V^T x_i = \tilde{x}_i$$

In the next lecture, we will discuss the clustering problem in more detail.
