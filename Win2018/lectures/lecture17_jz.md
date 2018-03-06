---
layout: page
mathjax: true
permalink: /Win2018/lectures/lecture17_jz/
---
## Lecture 17: Statistics for High-Dimensional Data

Tuesday 6 March 2018

*scribed by Andre Cornman and edited by the course staff*

-----------------

### Topics

1. <a href='#mt'>Multiple testing</a>
  - <a href='#bh'>Benjamini-Hochberg test</a>
  - <a href='#math'>Justifying the BH test</a>
3. <a href='#bayes'>Empirical Bayes theory</a>

### <a id='mt'></a>Multiple testing

For today's lecture, we will continue from where we left off at the end of last Tuesday's lecture. We will go back to talking about bulk RNA-Seq. Recall that this whole procedure can be summarized using the following block diagram:

<div class="fig figcenter fighighlight"> <img src="/Win2018/assets/lecture17/blocks.png" width="80%"> <div class="figcaption">The statistical testing framework for RNA-Seq.</div> </div>

Last time, we discussed the concept of a $$p$$-value where we let $$P_i$$ represent the $$p$$-value for the $$i$$th transcript. We have two important facts about $$p$$-values:

1. A $$p$$-value is a random variable
2. Under the null hypothesis, a $$p$$-value has uniform distribution

We discussed two methods for tackling the _multiple testing_ problem. Let $$V$$ represent the number of false discoveries out of $$m$$ tests (e.g. $$m$$ can represent the number of transcripts).:

1. **Control FWER** (Bonferroni procedure): The family-wise error rate FWER is defined as
$$Pr( V > 0)$$.
To control the FWER at 5% (i.e.
$$ FWER \leq 0.005$$
), we simply change our rejection threshold from 0.05 to 0.05/$$m$$. This approach is rather conservative, however. Perhaps we are willing to allow a few false positives if we can make several true discoveries.

2. **Control FDR** (Benjamini-Hochberg procedure): The false detection rate is defined as
$$ E\left[\frac{V}{\max(R, 1)} \right]$$
where $$R$$ represents the number of discoveries. We will discuss a procedure for controlling the FDR in detail below.

#### <a id='bh'></a>Benjamini-Hochberg test

Consider a list of $$p$$-values obtained from a large set of hypothesis tests (e.g. one for each transcript): $$p_1, \dots, P_m$$. The [Benjamini-Hochberg test](https://www.jstor.org/stable/2346101?seq=1#page_scan_tab_contents) was created in 1995 and consists of the following two steps:

1. Sort the $$p$$-values to obtain $$P_{(1)} \leq P_{(2)} \leq \dots \leq P_{(m)} $$

2. Looking at the $$\alpha/m$$ line, we reject the null for all tests (ordered by $$p$$-value) with $$p$$-values below this line

<div class="fig figcenter fighighlight"> <img src="/Win2018/assets/lecture17/bh.png" width="80%"> <div class="figcaption">Benjamini-Hochberg procedure figure borrowed from the PH525x series notes.</div> </div>

Intuitively, if we decrease $$\alpha$$, the slope of the line decreases, resulting in less $$p$$-values lying under the $$\alpha/m$$ curve and therefore less rejections. This makes sense since we expect fewer rejections with a more stringent $$\alpha$$.

For the Bonferroni procedure, we consider the hypotheses $$H_1, \dots, H_m$$. Importantly, the Bonferroni procedure compares each of the $$m$$ $$p$$-values _independently_ to a fixed threshold $$\alpha/m$$. For the Benjamini-Hochberg procedure, the decision we make on a test $$H_i$$ depends both on $$P_i$$ and the other $$p$$-values (captured by the sorting step). Therefore even though a $$p$$-value can be small, rejecting that $$p$$-value becomes harder if several other tests produce even smaller $$p$$-values.

**Example**

Consider the following set of $$p$$-values obtained from a differential expression test: 0.02, 0.03, 0.035, 0.006, 0.055, 0.047, 0.01, 0.04, 0.015, 0.025.

1. Under the Bonferroni procedure with $$\alpha = 0.05$$, we make discoveries on none of the tests.
2. Under the BH procedure with $$\alpha = 0.05$$, we make discoveries on the 1st, 2nd, 3rd, 4th, 7th, 8th, 9th, and 10th tests.

#### <a id='math'></a>Justifying the BH test

The BH-test, while straightforward, seems rather arbitrary at first glance. We will attempt to justify why the BH procedure is indeed valid and controls FDR at $$\alpha$$. For simplicity, we start with just one test. We know that under the null $$H_0$$, $$P \sim U[0, 1]$$ (i.e. the $$p$$-value has a uniform distribution). Under the alternate $$H_1$$, then $$P \sim f_1$$.

[IMAGE of f_1 with theta on horizontal axis, f_1 on vertical axis]

The false discovery rate captures the fraction of times the null is actually true despite being rejected:

$$
\begin{align*}
Pr(\text{null} | P \leq \theta) & = \frac{Pr(\text{null & } P \leq \theta)}{Pr (P \leq \theta)} \\
& = \frac{Pr(\text{null}) Pr(P \leq \theta | \text{null})}{Pr(\text{null}) Pr(P \leq \theta | \text{null}) + Pr(\text{alternate}) Pr(P \leq \theta | \text{alternate}) } \\
& = \frac{\pi_0 \theta}{\pi_0 \theta + (1-\pi_0) F_1(\theta)} \\
& = \frac{\pi_0 \theta}{F(\theta)}
\end{align*}
$$

where $$F_1$$ represents the CDF under the alternate, and $$F$$ represents the mixture of the null and alternative distributions. To control FDR to be $$\alpha$$, set $$\theta$$ such that

$$
\begin{align*}
\frac{\pi_0 \theta}{F(\theta)} = \alpha
\end{align*}
$$

where we have one equation and one unknown. Taking a step back, recall that Fisher was interested in

$$Pr(P \leq \theta | \text{null}) = \alpha.$$

Compare this to the probability expression we started with and note that they are not the same! In Fisher's case, we can just see $$\theta = \alpha$$, which is how he came up his rejection procedure. Importantly, to solve this equation, we do not need knowledge of $$F_1$$. For our probability of interest, however, we need both $$F_1$$ and $$\pi_0$$. Therefore Fisher's approach cannot solve the single-hypothesis-testing problem without making strong assumptions on what the alternative looks like.

For the multiple testing scenario, let's assume that $$P_1 \sim U[0, 1]$$ under the null, and $$P_1 \sim f_1$$ under the alternate. With, say, $$m = 20000$$ tests, we can now solve the above equation because we have 20000 observations drawn from $$F(\theta)$$, the mixture distribution. We note that

$$ \frac{\pi_0 \theta}{\alpha} = F(\theta) $$

has a CDF on the right-hand side and a line with slope $$\pi_0/\alpha$$ on the left-hand side.

[IMAGE of F, a cdf, intersecting with the pi/alpha slope, and an empirical approximation of $$F$$. Indicate the width and height of the stairs]

While we do not have $$F$$, we have enough samples to compute an empirical approximation of $$F$$. Note that for the empirical approximation, the height of each step is $$1/m$$. The width of the first step is exactly the smallest $$p$$-value. Intuitively, the sorting step of the BH procedure is captured by computing the approximate CDF.

We flip the two axes and stretch the (new) horizontal axis by $$m$$, resulting in a line with slope $$\alpha/(\pi_0 m)$$. How do we handle the $$\pi_0$$? We can replace $$\pi_0$$ by 1, resulting in a more conservative threshold. Since the less conservative threshold already controlled FDR at rate $$\alpha$$, increasing the stringency will still be valid.

A key assumption here is that the alternate has the same distribution across the tests. Second, we are assuming that the tests are independent. We will continue this discussion next time.
