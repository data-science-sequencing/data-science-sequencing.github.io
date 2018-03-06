---
layout: page
mathjax: true
permalink: /Win2018/lectures/lecture15/
---
## Lecture 15: Differential Analysis and Multiple Testing

Tuesday 27 February 2018

_scribed by Michelle Drews and edited by the course staff_

-----------------

## Topics

1.	<a href='#diff'>Differential analysis</a>
    - <a href='#classical'>Classical statistics approach and the _p_-value</a>
1.	<a href='#mt'>Multiple testing</a>


### <a id='diff'></a> Differential analysis

Our course so far has followed the high-throughput sequencing pipeline:

<div class="fig figcenter fighighlight">
  <img src="/Win2018/assets/lecture15/Fig1_EENotes.jpg" width="50%">
	<div class="figcaption">RNA-seq pipeline.</div>
</div>

Especially in the last few lectures, we have discussed the pipeline up to the transcript abundances in great detail. We now have several numbers describing the transcript or gene abundances of a biological sample. We can now use these numbers in downstream analysis to extract some information governing the biology of our sample.

We will first talk about the problem of _differential analysis_ where we attempt to identify the genes or transcripts that distinguish the samples from two populations (e.g. two conditions before/after drug treatment). Visually, differential analysis looks like the following:

<div class="fig figcenter fighighlight">
  <img src="/Win2018/assets/lecture15/Fig2_EENotes.jpg" width="50%">
	<div class="figcaption">Differntial analysis.</div>
</div>

Examples of software that handle differential analysis are [DESeq](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-10-r106) and [DESeq2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8). A key idea here is the idea of _significance_. It's not enough for two expression levels to be different; we need some notion of significance to quantify whether the difference observed is actually meaningful.

As an aside, while we are talking about gene/transcript abundances here, we will note that other units may be used in practice. For example, raw integer counts (e.g. molecule or read counts aligning to a particular gene) are also commonly used. The discussions here will apply to those units as well.

Suppose we sequence samples from two different conditions: sample A is a control group and sample B is a group treated with a drug. Suppose we sequence the samples from each sample, quantify transcript abundances with some EM-based algorithm, and obtain the following plots:

<div class="fig figcenter fighighlight">
  <img src="/Win2018/assets/lecture15/Fig3_EENotes.jpg" width="50%">
	<div class="figcaption">Differntial analysis by eyeballing.</div>
</div>

Just by eyeballing these plots, we might immediately conclude that transcript 1 is differentially expressed and transcripts 2 and 3 are not. Note that there are no error bars on this plot, and therefore we cannot really make any strong conclusion. Perhaps the observed difference, while large, is due to some natural variation and not due to the drug. With only one sample in each population, we cannot get error bars, which capture a sense of the variation in the estimate.

We will assume that we have more than one sample in the form of experimental _replicates_ (repeated sampling from the same condition). Assume we have $$n$$ replicates. Let $$X_i^A$$ represent the abundance of transcript $$i$$ under condition $$A$$. We can estimate the mean abundance for each population as

$$
\begin{align*}
\hat{\mu}^A & = \frac{1}{n}\sum_{i=1}^n X_i^A \\
\hat{\mu}^B & = \frac{1}{n}\sum_{i=1}^n X_i^B.
\end{align*}
$$

We can also estimate sample variances using:

$$
\begin{align*}
(\hat{\sigma}^A)^2 & = \frac{1}{n-1} \sum_{i=1}^n (X_i^A - \hat{\mu}^A )^2 \\
(\hat{\sigma}^B)^2 & = \frac{1}{n-1} \sum_{i=1}^n (X_i^B - \hat{\mu}^B )^2
\end{align*}
$$

The $$n-1$$ here (rather than $$n$$) ensures that our variance estimator is unbiased. Note that in reality, $$n$$ scales proportionally with the cost of the experiment (sequencing is expensive!). Therefore we would like to know how large we should make our $$n$$ so that we can minimize our cost while getting a reasonable estimate of our sample variance. We will discuss this in more detail later.

#### <a id='classical'></a> Classical statistics approach and the _p_-value

Assuming we have more than one replicate in each condition, we can draw error bars for each transcript. Using the estimated mean and standard deviation, we would like to test if the two populations are different. This problem was originally solved by [RA Fisher](https://en.wikipedia.org/wiki/Ronald_Fisher), who introduced the idea of a [_p_-value](https://en.wikipedia.org/wiki/P-value):

[IMAGE: muA & sigmaA, muB & sigmaB -> [testing] -> p-value]

Critically, the notion of a $$p$$-value only makes sense in the context of a null hypothesis. Consider the following procedure

1. Calculate a statistic $$T$$
2. Set a null hypothesis: abundances in A & B have the same distribution
3. Under the null, suppose $$T \sim G$$ where $$G$$ is a cumulative distribution function (CDF). Then $$p = 1-G(T)$$

In other words, the _p_-value is the probability that statistic described by the null distribution exceeds the observed statistic. Typically, $$G$$ is known (e.g. normal with the same mean and variance across the two populations).

**Example**

$$T \triangleq \frac{\sqrt{n}(\hat{\mu}^A-\hat{\mu}^A)}{\sqrt{(\hat{\sigma}^A)^2 + (\hat{\sigma}^A)^2}} $$

The statistic here is known as the _t_-statistic and is used for a _two-sample t-test_. Under the null hypothesis, we assume

$$ X_i^A, X_i^B \sim N(\mu, \sigma^2), $$

in which case $$p \sim U[0, 1]$$, a uniform distribution between 0 and 1.

**Fact**: $$Y$$ with cdf $$F(y) \triangleq \text{Pr}(Y \leq y) $$, then $$F(Y) = Z U[0, 1]$$.

We can use such a test to determine if the mean value of a transcript is the same between conditions A and B based on the _p_-value obtained. While we only described one statistic here, there are several types of statistics one could use depending on the assumptions one would like to make on the null hypothesis (e.g. non-equal variance).

### <a id='mt'></a> Multiple testing

For each of the, say, 30000 genes, we can perform a two-sample _t_-test to determine if that gene is different between the two populations. We would therefore obtain 30000 _p_-values between 0 and 1. Intuitively, the smaller the _p_-value, the more differentially expressed the gene (the less likely we are to observe a value at least as extreme as the actual observed value under the null hypothesis). Fisher gave also gave us a _p_-value threshold 0.05 where a test is only deemed significant if it produces a _p_-value under this threshold.

Notice that with 30000 tests running in parallel, we would expect certain genes to appear significant due to sheer randomness. A _p_-value of 0.05, after all, means that 1 in 20 samples drawn from the null distribution would appear to be significant despite being drawn from the null distribution.

Let $$m = $$ the number of transcripts. Let $$H_i$$ represent the null hypothesis that transcript $$i$$ is not expressed differentially, and let $$p_i$$ be the associated _p_-value. A _discovery_ indicates a transcript found to be significant. We would like to introduce some sort of correction to reduce the amount of _false discoveries_ (i.e. incorrect rejections of the null hypothesis) without significantly compromising the amount of _true discoveries_.

Let $$V = $$ total number of false discoveries. One way of controlling $$V$$ is by controlling the probability that we get at least one false discovery, also known as the _family-wise error rate_ (FWER):

$$ \text{Pr}[V \geq 1] \leq 0.05.$$

Note that this event is a union of the events of making false discover in the first test, or the second test, or the $$n$$-th test. If we say that

$$ P_i \leq \theta$$

then the probability of falsely rejecting the null is $$\theta$$. Taking a union, we getting

$$
\begin{align*}
\text{Pr}[V \geq 1] & \leq m \theta \leq 0.05 \\
\implies \theta & = \frac{0.05}{m} \\
& = 5 * 10^{-6}
\end{align*}
$$

The procedure of dividing $$\theta$$ by $$m$$ is known as the _bonferroni correction_. The resulting threshold is perhaps too conservative in practice, but it's the price we must pay in order to guarantee that in 10000s of tests, we have almost no false discoveries. A fundamental question here is: would making no discoveries under this stringent criteria be better than making a few false discoveries in addition to some true discoveries?

In 1995, [Benjamini and Hochberg](https://www.jstor.org/stable/2346101) proposed looking the _false discovery rate_ (FDR) instead of the FWER. Letting $$T = $$ the total number of discoveries,

$$ \text{FDR} \triangleq E \left[ \frac{V}{R} \right]. $$

With FDR $$ \leq 0.05 $$, we would be saying that of the discoveries we have reported, 5% of them are false. Of course, we do not know which ones are false. Note that both $$V$$ and $$R$$ are random variables, and therefore we need to take an expectation. Importantly, Benjamini and Hochberg does not only give us a new metric to control, but also a procedure for controlling the metric. We will discuss this in further detail next time.
