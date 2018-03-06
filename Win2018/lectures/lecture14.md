---
layout: page
mathjax: true
permalink: /Win2018/lectures/lecture14/
---

## Lecture 14 - RNA-seq - Quantification and the EM Algorithm Part 2

Tuesday 20 February 2018

_scribed by Mira Moufarrej and edited by the course staff_

-----------------


## Topics

1. <a href='#int'>Estimating RNA Abundance</a>
1. <a href='#hs'>Hard vs Soft EM</a>
1. <a href='#mm'>EM as alternative maximization</a>
1. <a href='#pseudo'> Pseudo-alignment</a>

### <a id='int'></a>Estimating RNA Abundance
When we sequence RNA, we would like to quantify the abundances of distinct transcripts.
As we have seen previously, NGS yields short reads that can map back to multiple transcripts.
For instance, if we have a read composed of an exon $$B$$, and we have two transcripts where
transcript 1 is $$ABC$$ and transcript 2 is $$BCD$$, the read maps to both transcripts
1 and 2. This is known as the _multiple mapping problem_.  

Here, we want to infer each transcript's abundance based on the read data.
Let $$Y$$ be the data (reads). Let $$\rho$$ be a vector of unknown abundances for
each transcript, which we want to estimate. Using maximum likelihood (ML), we can
infer the most likely estimate of $$\rho$$ based on the data $$Y$$. Using the ML
formula, we would get the following,

$$\hat{\rho}_{ML} = \text{argmax}_\rho \log P(Y;\rho)$$


However, this
formula is computationally intractable. We would like to infer $$\rho$$ based
on the data, but $$Y$$ inherently depends on $$\rho$$.
Therefore, we try to solve an easier problem.

Let $$Z = (Z_{ik})$$ be a hidden variable such that $$Z_{ik} = 1$$ if read $$i$$ comes from transcript $$k$$. If we knew $$Z$$, we will see that
ML inference becomes a lot easier, and we could write the ML problem as:

$$\hat{\rho}_{ML} = \text{argmax}_\rho \log P(Y,Z;\rho)$$

Using the
chain rule to expand this joint probability, we get that:

$$\hat{\rho}_{ML} = \text{argmax}_\rho \log P(Y|Z;\rho)P(Z;\rho)$$

Since $$P(Y|Z)$$ has no
dependence on $$\rho$$, we can treat it as a constant and rewrite the expression
above as maximizing only the second term:

$$\hat{\rho}_{ML} = \text{argmax}_\rho \log P(Z;\rho)$$

But what is $$P(Z;\rho)$$ and is it tractable? Assuming all transcripts are the
same length, we can rewrite $$P(Z;\rho)$$ as follows:

$$P(Z;\rho) = \prod_{i=1}^N \prod_{k=1}^K \rho_k^{Z_{ik}}$$

where $$Z_{ik}$$ is an indicator function that equals 1 when read $$i$$
comes from transcript $$k$$. If transcripts were of different lengths,
we would scale by length. Taking the log of the expression for $$P(Z;\rho)$$,
we get the following:

$$P(Z;\rho) = \sum_{k=1}^K\sum_{i=1}^N Z_{ik} \log(\rho_k)$$

Notice that in this expression, each transcript contributes one log term and
$$\sum_{i=1}^N Z_{ik}$$ is the number of reads sampled from transcript $$k$$.
This expression is computationally tractable. Therefore if we knew $$Z$$,
we would be able to solve for $$\hat{\rho}_{ML}$$. This is the same as the estimation problem for when all reads are uniquely mapped. By revealing $$Z$$, one converts the multiple mapping problem to a unique mapping problem.

But because we don't know $$Z$$, we can use EM to maximize both $$\rho$$ and $$Z$$ iteratively. Why does this work?
Notice that in the original ML formula, we had that:

$$\rho_{ML} = \max_\rho \log P(Y;\rho)$$

This is essentially equivalent to marginalizing $$Z$$ out of the sum
(computationally intractable), which can be approximated as maximizing over the
joint probability of $$Y$$ and $$Z_{MAP}$$ where $$Z_{MAP}$$ is the _maximum a
posteriori_ estimate of $$Z$$ given $$Y$$.

$$ \hat{\rho}_{ML} = \text{argmax}_\rho \log P(Y;\rho) = \text{argmax}_\rho \log \sum_Z P(Y,Z;\rho) \approx \text{argmax}_\rho \log P(Y,Z_{MAP};\rho)$$

However, this is a chicken and egg problem since computing $$Z_{MAP}$$ requires knowledge of $$\rho$$.


### <a id='hs'></a>Hard vs Soft EM

We start with some initial value for $$\rho$$ called $$\rho^{(m)}$$. Given these values of $$\rho$$, we find the values of $$Z$$ that maximizes the probability in the objective (_maximum a posteriori_).

$$ \max_\rho \log P(Y, Z_\text{map}(Y); \rho). $$

This type of EM is called _hard_ EM. This inference process is attempting to estimate the hidden variable (the transcript a read comes from) for each read. If a read could potentially come from, say, 3 different transcripts, we pick the transcript with maximum probability. In _hard_ EM, we make a hard decision on the hidden variable at each iteration. The drawback with this approach is that we do know for sure that a read comes from the max-probability transcript, and it may feel a bit aggressive.

Alternatively, instead of making a hard decision, we can instead average over the posterior distribution before continuing the inference step, resulting in a _soft_ version. This results in the EM algorithm we discussed in the last lecture.

### <a id='mm'></a>EM as alternative maximization

So far, we have attempted to justify EM from an intuitive point of view. We still have a problem, however, as the current method still seems somewhat heuristic. How do we know if this algorithm actually gives us the ML solution? To understand the connection of EM with maximum likelihood, we can define

$$F(Q_Z, \rho) \triangleq \log P(Y; \rho) - D(Q_Z \| P(Z|Y; \rho)) $$

where $$D(P_1 \| P_2)$$ is the Kullback-Leibler divergence ([KL-divergence](https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence)):

$$ D(P_1 \| P_2)  = \sum_z P_1(z) \log \frac{P_1(z)}{P_2(z)} \geq 0 $$

and intuitively measures the similarity between the two distributions $$P_1$$ and $$P_2$$. Note that if we maximize over $$Q_Z$$, we will get
$$D(Q_Z \| P(Z|Y; \rho)) = 0 $$
, resulting in $$F(Q_Z, \rho) = \log P(Y; \rho)$$. Therefore

$$ \max_\rho \log P(Y; \rho) \iff \max_\rho \max_{Q_Z} F(Q_Z; \rho). $$

This gives us an alternate maximization procedure:

1. Fix $$\rho$$, max over $$Q_Z$$
2. Fix $$Q_Z$$, max over $$\rho$$

and the problem now becomes: is this procedure easy to do? We see that

$$
\begin{align*}
F(Q_Z, \rho) & = \log P(Y; \rho) - \sum_Z Q(Z) \log \frac{Q(Z)}{P(Z|Y; \rho)} \\
& = \log P(Y; \rho) + \sum_Z Q(Z) \log \frac{P(Z|Y; \rho)}{Q(Z)} \\
& = \sum_Z Q(Z) \log P(Y; \rho) + \sum_Z Q(Z) \log \frac{P(Z|Y; \rho)}{Q(Z)} \\
& = \sum_Z Q(Z) \left[\log \frac{P(Y; \rho) P(Z|Y; \rho)}{Q(Z)} \right] \\
& = \sum_Z Q(Z) \left[\log \frac{P(Y, Z; \rho)}{Q(Z)} \right]
\end{align*}
$$

For step 2, we have $$Q_Z$$ fixed, and therefore maximizing $$F(Q_Z, \rho)$$ is just maximizing this final expectation over $$\rho$$. In other words, step 2 is equivalent to

$$ \max_\rho E_Q [\log P(Y, Z; \rho)] $$

for fixed $$Q$$. We see that the two steps of the alternate maximization procedure correspond exactly to the E and M steps in the EM algorithm. Technically, the E step here also requires maximization, and therefore perhaps MM is a more suitable name (see [here](https://en.wikipedia.org/wiki/MM_algorithm)).

As a note about existing software, [Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/) and [RSEM](https://deweylab.github.io/RSEM/) are both early tools that implement this EM approach.


<!-- ## Hard EM

**EM** is an iterative method to find the maximum likelihood (ML) estimate for a parameter of interest. Here, we want to quantify RNA abundance based on reads.  

EM is composed of two steps: (1) The expectation step or E-step (2) The maximization step or M-step. In hard EM, for a given iteration $m$, we infer $Z$ using MAP and our current estimate for $\rho$ as seen below:

$$\rho^{(m)} = max_Z log P(Y,Z;\rho^{(m)}) = max_Z log P(Z|Y;\rho^{(m)})P(Y;\rho^{(m)})$$

Notice that here, we are using a hard decision to estimate $\rho$ each time based on our current estimate. One drawback of this is that we do not necessarily know that a given read comes from a given transcript and therefore a hard decision as seen above might be too strong of a decision. Instead, we could think of using soft EM.

## Soft EM

Soft EM does not make a hard cutoff for $\rho$ on each iteration. Instead, it takes the expectation of $p(Y,Z;\rho)$ and then averages over the posterior.  

Before we go into the nuts and bolts of the algorithm, let's first take a step back and ask the following. How do we know that this algorithm will converge on the ML solution?  

To answer this question, we can think of EM as an alternative maximization to ML. In ML, we want to maximize the following $max_\rho log P(Y;\rho)$. To define our alternative maximization, let us define a new function, $F(Q_Z, \rho)$, which is defined as follows: $$F(Q_Z,\rho)  = logP(Y;\rho) - D(Q|| P(Z|Y;\rho)$$

The term $D(Q_Z|| P(Z|Y;\rho) = \sum_Z Q_Z(z) log \frac{Q_Z(Z=z)}{P(Z=z|Y;\rho}$ is the Kullback-Leibler (KL) divergence between two distributions, $Q_Z$ and $P(Z|Y;\rho)$. Importantly the KL divergence is always non-negative and it equals 0 if two distributions are equivalent.  

To maximize $F(Q_Z, \rho)$, we would have to maximize over $Q_Z$ such that $Q_Z = P(Z|Y,\rho)$ and $D=0$. Notice that in the case where $D=0$, $F(Q_Z,\rho) = logP(Y;\rho)$, which is exactly the expression we wish to maximize over in original ML. Therefore, we can rewrite our original ML expression in terms of $F$ as follows: $$max_\rho log P(Y;\rho) \Longleftrightarrow max_\rho max_{Q_Z} F(Q_Z; \rho)$$

Strangely, this means we have converted a single maximization problem into a double maximization. As Prof. Tse says, "Life is strange. Simplification sometimes only comes with complication." To compute this double maximization, we iterate over two steps (E-Step and M-Step) as follows:  

**E-Step**: Fix $\rho$ and maximize $Q_Z$  
**M-Step**: Fix $Q_Z$ and maximize $\rho$  

Are both of these steps computationally tractable? Let's expand and simplify the expression we wish to maximize.  

$$F(Q_Z, \rho) = logP(Y;\rho) - \sum_Z Q(z)log \frac{Q(z)}{p(Z=z|Y;\rho} = logP(Y;\rho) + \sum_Z Q(z)log \frac{p(Z=z|Y;\rho}{Q(z)}$$  

In the M-Step when $Q(z)$ is fixed, we have that $Q(z)$ is a constant, and we can rewrite the expression above as follows using the chain rule.  

$$F(Q_Z, \rho) = \sum_Z Q(z)log\frac{P(Y;\rho)P(Z|Y;\rho)}{Q(z)} = \sum_ZQ(z)logP(Y,Z)$$  

Notice that the expression on the right is precisely the expectation of a logarithm, hence where the algorithm gets its name. Maximizing for $\rho$ in the M-step, we get the following:  

$$max_\rho \sum_Z Q(z)log P(Y,Z;\rho) = max_\rho E_Q log P(Y,Z;\rho)$$  

Similarly, we could fix $\rho$ and solve for $Q_z$ as done in the E-step. Popular software that implements EM includes Cufflinks and RSEM. -->

<!-- ## Psuedo-alignment

When scientists started quantifying abundances, they noticed it took a very long time. It turns out that alignment is a time-consuming algorithm (quadratic in read length and number of transcripts) . For a transcriptome, each read is mapped to each transcript, which creates a bottleneck. To alleviate this bottleneck, they came up with pseudo-alignment, which performs only what is needed to do EM after alignment.

For alignment, we start with ML again, $max_\rho \sum_{i=1}^N log(Y\rho)_i$ where Y is a matrix where each row is a read and each column is a transcript. $Y_{ij} = 1$ if a read $i$ maps to a transcript $j$. The idea behind pseudo-alignment is as follows: With a reliable estimate of $Y$, we can estimate $\rho$.  

So, can we compute Y without full alignment? What we really want is the following: For each read $r$, compute the subset $S$ of transcripts from which $r$ can come.  To do so, we can get rough information using fingerprinting with k-mers derived from the read. For each k-mer, we can ask where does it appear on the genome, which is indexed, and then use a hashtable to look up the k-mer later. For a given read's k-mers, we could find the set of transcripts that intersect those k-mers. Since we are not keeping the position of successive k-mer indices, the true transcripts that a read comes from are a strict subset of this intersection.

So how does this cut down run-time? One k-mer could be contained in many transcripts and therefore would yield redundant information. We could build a DeBruijn graph using the k-mers and overlay transcripts to find transcripts that correspond to certain reads. There is software that implements this. -->

### <a id='pseudo'></a>Pseudo-alignment

Taking a step back, recall that the RNA-Seq problem involves determining which
of several transcripts in a transcriptome a read comes from. The EM step is
actually computationally cheap compared to aligning 100M reads, and therefore
_alignment_ is the computational bottleneck.
But for the purposes of EM, do we really need alignment?

The fundamental problem here is the ML problem

$$ \max_\rho \sum_{i=1}^N \log (Y \rho)_i $$

where we have a data matrix

$$ Y =
\begin{bmatrix}
0 & 0 & 1 & 0 & 1 & 0 \\
\vdots
\end{bmatrix}
$$

with $$Y_{ij} = 1$$ if read $$i$$ aligns to transcript $$j$$.

Starting in about 2013, people started asking if we could recover $$Y$$ without
performing full alignment. For each read $$r$$, we need to compute $$S$$ of
transcripts from which $$r$$ can come. Recall that for aligning reads to a long
genome, we broke each read up into $$k$$-mers. For each $$k$$-mer, we could
quickly find where it maps to using a hash table. In other words, we indexed
the genome first by building a hash table where $$k$$-mers are keys. This is
much faster than doing full-scale alignment. We exploit the fact that even
in the event of errors, a shorter $$k$$-mer sequence is less likely to have
errors, and therefore we have a reasonable chance of getting exact mappings
of the $$k$$-mers.

Similarly, we can also index the transcriptome by building a hash tables of
all the $$k$$-mers. As an aside, we cannot choose even values for $$k$$
(see assignment 3). Each 31-mer, for instance, will map to set of transcripts
that $$k$$-mer can come from. Note that we do not concern ourselves too much
with the computational cost associated with building this table as we only need
to do it once (and thus the cost is amortised).
The procedure will look something like

1. Build index (hash table)
2. Quantification

For each read, we will attempt to find the transcripts that the read belongs to.
For each $$k$$-mer in the read, we obtain a set of transcripts (a few $$k$$-mers
will be erroneous, but they will likely not map to any entry in the table).
To find the unique set of transcripts that all the $$k$$-mers can come from,
we can just take an intersection of the sets. Therefore the "alignment"
step here consists of

1. Break read into $$k$$-mers
2. For each $$k$$-mer, get the set of transcripts that $$k$$-mer maps to
3. Take intersection of all transcript sets

Note that since we are not tracking positional information (yet), these
transcript sets are a bit bigger than they could be. To incorporate the
positional information, we can use a data structure similar to a
$$k$$-mer overlap graph:

<div class="fig figcenter fighighlight">
  <img src="/Win2018/assets/lecture14/pseudo.png" width="80%">
	<div class="figcaption">Pseudo-alignment figure from kallisto paper.</div>
</div>

This strategy, known as _pseudo-alignment_, is used by
[kallisto](https://www.nature.com/articles/nbt.3519) and
[Salmon](https://www.nature.com/articles/nmeth.4197).
