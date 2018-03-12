---
layout: page
mathjax: true
permalink: /Win2018/lectures/lecture18/
---
## Lecture 18: Empirical Bayes and Principal Component Analysis

Thursday 8 March 2018

*scribed by Brijesh Patel and edited by the course staff*

-----------------

### Topics

1. <a href='#fdr'>Recap of False Discovery Rate</a>
1. <a href='#eb'>Empirical Bayes</a>
1. <a href='#pca'>Principal component analysis</a>

### <a id='fdr'></a>Recap of False Discovery Rate

We have a $$m$$-hypothesis testing problem where we test whether each transcript is differentially expressed or not in condition $$A$$ vs condition $$B$$. We obtain $$m$$ $$p$$-values from $$m$$ statistical tests (one for each transcript), and we want to know which transcripts are actually differentially expressed.

One approach for justifying the Benjamini-Hochberg (BH) procedure is by viewing the $$p$$-values as drawings of some random variable $$P$$ from an underlying distribution $$F$$. Recall that under the null hypothesis, the $$p$$-value is distributed $$U[0, 1]$$. Under the alternate hypothesis, the $$p$$-value will come from the some non-uniform (but unknown) distribution. The BH procedure attempts to estimate an empirical distribution $$\hat{F}$$ from the data.

Once $$\hat{F}$$ is computed, we can perform a transcript-by-transcript test where we evaluate if $$p_i \leq \theta$$ for some decision threshold $$\theta$$. From the discussion [last lecture](/Win2018/lectures/lecture17/), we showed how $$\theta$$ can be obtained using the expression

$$ \frac{\pi_0 \theta^*}{\hat{F}(\theta^*)} = \alpha$$

where $$\alpha$$ is our lever of FDR control. If you know $$F$$, then you can set your threshold to satisfy this equation. With our estimated $$\hat{F}$$, we can rearrange the above expression to obtain

$$ \frac{\pi_0 \theta^*}{\alpha} = \hat{F}(\theta^*). $$

At any point $$\theta$$, $$\hat{F}$$ is (roughly speaking) the number of $$p_i$$'s which is less than $$\theta$$:

$$ \hat{F}(\theta) = \frac{\text{# of } p_i \leq \theta}{m}. $$

Sorting the $$p_i$$'s, we obtain $$ p_{(1)} < p_{(2)} < \dots < p_{(m)}$$ where the inequalities are strict because the $$p$$-values are continuous. This results in a step-like function for $$\hat{F}$$ as shown in the figure below.

<div class="fig figcenter fighighlight">
  <img src="/Win2018/assets/lecture18/fdr.png" width="40%">
	<div class="figcaption">Justifying the decision boundary for the BH procedure.</div>
</div>

In conclusion, the BH procedure was a significant development in statistics because it says underlines the importance of multiple testing. Multiple testing can be viewed as an advantage of estimating, allowing us to estimate the distribution $$F$$ so that instead of bounding a $$p$$-value, which is the probability we get a wrong discovery given that the null distribution is true, we can instead bound the probability we've made a mistake given that we've committed to a discovery. The latter requires modeling $$F$$, which is a combination of the null and the alternate distributions. We cannot obtain a good estimate of $$F$$ without multiple testing.

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

Here, $$X_i^A$$ represents the expression level in the $$i$$th replicate under condition $$A$$. Because getting replicates is expensive, we are interested in understanding the minimum number of replicates required such that we can draw reasonable conclusions. From a statistical testing point of view, our null hypothesis is that the two means and variances are equal:

$$
\begin{align*}
H_0: & X_i^A \sim N(\mu, \sigma^2) \\
& X_i^B \sim N(\mu, \sigma^2)
\end{align*}
$$

Under the null, the variation stems from both technical and biological variation for both conditions. We use the mean and variance estimates above to compute a $$t$$-statistic:

$$ \hat{T}_n = \frac{\hat{\mu}^A-\hat{\mu}^B}{\sqrt{\frac{(\hat{\sigma}^A)^2 + (\hat{\sigma}^B)^2}{2}}}. $$

The bigger the difference in means (i.e. the bigger the $$T$$-statistic), the more likely the null is invalid. In other words, whenever $$\hat{T}_n$$ is large in magnitude, we are "suspicious of the null." We start by comparing the $$t$$-distribution to the standard normal distribution $$N(0, 1)$$:

<div class="fig figcenter fighighlight">
  <img src="/Win2018/assets/lecture18/tails.png" width="40%">
	<div class="figcaption">Comparing the t-distribution with 2 degrees of freedom (red) to the standard normal distribution (blue).</div>
</div>

Note that the standard normal has smaller tails. This means that controlling the false positive rate for the $$t$$-distribution is more difficult, and we must use a larger cutoff in order to control the false positive rate. If we knew the variance $$\sigma^2$$, then the distribution of the $$T$$ statistic would be Gaussian (a linear combination of Gaussians is Gaussian):

$$ \tilde{T} = \frac{\hat{\mu}^A - \hat{\mu}^B}{\sigma}. $$

This highlights how the estimate of the variance dictates how good our test would be. Estimating this quantity introduces a lot of uncertainty in our distribution (hence fattening the tails). As shown in the figure above, having only 3 replicates does not give us a reliable estimate of the variance, resulting in a low-power test. A larger variance estimate results in a smaller $$t$$-statistic. We can take a closer look at how the variance changes as a function of the data for $$n = 2$$:

$$ (\hat{\sigma}^A)^2 = (x_1 - \hat{\mu}^A)^2 + (x_1 - \hat{\mu}^B)^2. $$

<div class="fig figcenter fighighlight">
  <img src="/Win2018/assets/lecture18/sigma.png" width="40%">
	<div class="figcaption">Variance estimate with respect to the data.</div>
</div>

Because there is a large probability that $$x_1$$ is small (red area under the curve), our estimate of the variance and therefore the $$T$$ statistic will be poor. As we increase the number of replicates, our performance gets better quickly.

<div class="fig figcenter fighighlight">
  <img src="/Win2018/assets/lecture18/tailsthin.png" width="80%">
	<div class="figcaption">Comparing the t-distribution with varying degrees of freedom (red) to the standard normal distribution (blue) as the number of replicates increases.</div>
</div>

Empirical Bayes tells us that we can gain power by estimating some statistics associated all the data (i.e. by estimating some global prior using the data). [DESeq](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-10-r106) leverages the information across genes to power a statistical test for differential expression. With 20000 genes, for example, we can plot the variance with respect to the mean. While the estimate for an individual gene may be poor, we can fit some (often parametric) model for all genes. This process is also known as _shrinkage_ because we shrink the points towards our model, resulting in a tightening of the distribution. As shown in the figure, a Poisson model (blue) underestimates the variance in real data while a negative binomial distribution (orange) tends to fit the data better.

<div class="fig figcenter fighighlight">
  <img src="/Win2018/assets/lecture18/deseq.png" width="40%">
	<div class="figcaption">Modeling the expression across all genes in an RNA-Seq experiment using a Poisson model (purple) and a negative binomial model (orange). DESeq uses the latter model. </div>
</div>

### <a id='pca'></a>Principal component analysis

For the last part of the course, we will discuss various computational problems that arise from single-cell RNA-Seq analysis. A [recent dataset](https://community.10xgenomics.com/t5/10x-Blog/Our-1-3-million-single-cell-dataset-is-ready-to-download/ba-p/276) published by 10x genomics consists of 1.3 millions cells and 28000 genes. Each cell is represented by a 28000-dimensional expression vectors.

<div class="fig figcenter fighighlight">
  <img src="/Win2018/assets/lecture18/10x.png" width="80%">
	<div class="figcaption">The 10x Genomics 1.3 million cell dataset visualized in 2 dimensions using t-stochastic neighbor embedding.</div>
</div>

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
