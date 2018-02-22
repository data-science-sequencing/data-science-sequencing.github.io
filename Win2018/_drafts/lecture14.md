---
layout: page
mathjax: true
permalink: /Win2018/lectures/lecture14/
---
## Lecture 14: RNA-seq - Wrapping up EM and Pseudo-alignment

Tuesday 22 February 2018

-----------------

## Topics

1. <a href='#hs'>Hard vs Soft EM</a>
1. <a href='#mm'>EM as alternative maximization</a>
1. <a href='#pseudo'> Pseudo-alignment</a>

### <a id='hs'></a>Hard vs Soft EM

For the expectation-maximization (EM) framework, we have data $$Y$$, unknowns $$\rho$$, and some hidden variables $$Z$$. In the context of RNA-Seq, we continue to assume that all transcript lengths are the same length. Recall that the EM algorithm consists of two steps:

1. E-step: compute
$$P(Z|Y; \rho^{(m)})$$
2. M-step:
$$\max_{\rho} E_{P(Z|Y; \rho^{(m)})} [\log P(Y, Z; \rho)] \rightarrow \rho^{(m+1)}$$

When viewing this problem from a maximum likelihood (ML) perspective, we need to solve the optimization problem
$$\max_{\rho} \log P(Y; \rho).$$

We can simplify this problem by introducing hidden variables $$Z$$, resulting in

$$ \max_\rho \log(P(Y, Z; \rho)) = \max_\rho \log(P(Z; \rho) P(Y|Z; \rho)) \iff \max_\rho \log P(Z; \rho) $$

where we exploited the fact that
$$P(Y | Z; \rho) = 1$$. Now,

$$P(Z; \rho) = \prod_{i=1}^N \prod_{k=1}^K \rho_k^{Z_{ik}} $$

where $$Z_{ik} = 1$$ if read $$i$$ comes from transcript $$k$$. Note that this product will consist of exactly one $$\rho$$ term since read $$i$$ can only equal come from one transcript. Therefore

$$\log P(Z; \rho) = \log \prod_{i=1}^N \prod_{k=1}^K \rho_k{Z_{ik}} = \sum_{k=1}^K \left(\sum_{i=1}^N z_{ik}\right) \log \rho_k$$

We introduce $$Z$$ into the problem to remove the original ambiguity; however we do not know the value of $$Z$$. Going back to the original ML problem

$$
\begin{align*}
& \max_\rho \log P(Y; \rho) \\
= & \max_\rho \log \sum_Z P(Y, Z; \rho)
\end{align*}
$$

which gives us

$$\rho^{(m)} \rightarrow \max_Z P(Y, Z; \rho^{(m)}) = \max_Z P(Z| Y; \rho^{(m)}) P(Y; \rho^{(m)}) $$

We start with some initial value for $$\rho$$ called $$\rho^{(m)}$$. Given these values of $$\rho$$, we find the values of $$Z$$ that maximizes the probability in the objective (_maximum a posteriori_).

$$ \max_\rho \log P(Y, Z_\text{map}(Y); \rho). $$

This type of EM is called _hard_ EM. This inference process is attempting to estimate the hidden variable (the transcript a read comes from) for each read. If a read could potentially come from, say, 3 different transcripts, we pick the transcript with maximum probability. In _hard_ EM, we make a hard decision on the hidden variable at each iteration. The drawback with this approach is that we do know for sure that a read comes from the max-probability transcript, and it may feel a bit aggressive.

Alternatively, instead of making a hard decision, we can instead average over the posterior distribution before continuing the inference step, resulting in a _soft_ version.

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

### <a id='pseudo'></a>Pseudo-alignment

Taking a step back, recall that the RNA-Seq problem involves determining which of several transcripts in a transcriptome a read comes from. The EM step is actually computationally cheap compared to aligning 100M reads, and therefore _alignment_ is the computational bottleneck. But for the purposes of EM, do we really need alignment?

The fundamental problem here is the ML problem

<!-- $$ \max_\rho \sum_{i=1}^N \log (Y \rho)_i $$ -->

where we have a data matrix

$$ Y =
\begin{bmatrix}
0 & 0 & 1 & 0 & 1 & 0 \\
\vdots
\end{bmatrix}
$$

with $$Y_{ij} = 1$$ if read $$i$$ aligns to transcript $$j$$.

Starting in about 2013, people started asking if we could recover $$Y$$ without performing full alignment. For each read $$r$$, we need to compute $$S$$ of transcripts from which $$r$$ can come. Recall that for aligning reads to a long genome, we broke each read up into $$k$$-mers. For each $$k$$-mer, we could quickly find where it maps to using a hash table. In other words, we indexed the genome first by building a hash table where $$k$$-mers are keys. This is much faster than doing full-scale alignment. We exploit the fact that even in the event of errors, a shorter $$k$$-mer sequence is less likely to have errors, and therefore we have a reasonable chance of getting exact mappings of the $$k$$-mers.

Similarly, we can also index the transcriptome by building a hash tables of all the $$k$$-mers. As an aside, we cannot choose even values for $$k$$ (see assignment 3). Each 31-mer, for instance, will map to set of transcripts that $$k$$-mer can come from. Note that we do not concern outselves too much with the computational cost associated with building this table as we only need to do it once (or at least, very infrequently). The procedure will look something like

1. Build index (hash table)
2. Quantification

For each read, we will attempt to find the transcripts that the read belongs to. For each $$k$$-mer in the read, we obtain a set of transcripts (a few $$k$$-mers will be erroneous, but they will likely not map to any entry in the table). To find the unique set of transcripts that all the $$k$$-mers can come from, we can just take an intersection of the sets. Therefore the "alignment" step here consists of

1. Break read into $$k$$-mers
2. For each $$k$$-mer, get the set of transcripts that $$k$$-mer maps to
3. Take intersection of all transcript sets

Note that since we are not tracking positional information (yet), these transcript sets are a bit bigger than they could be. To incorporate the positional information, we can use a data structure similar to a $$k$$-mer overlap graph:

<div class="fig figcenter fighighlight">
  <img src="/Win2018/assets/lecture14/pseudo.png" width="80%">
	<div class="figcaption">Pseudo-alignment figure from kallisto paper.</div>
</div>

This strategy, known as _pseudo-alignment_ is used by the tool [kallisto](https://www.nature.com/articles/nbt.3519).
