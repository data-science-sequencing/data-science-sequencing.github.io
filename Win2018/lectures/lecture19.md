---
layout: page
mathjax: true
permalink: /Win2018/lectures/lecture19/
---
## Lecture 19: Clustering

Tuesday 13 March 2018

*scribed by Daniel Hsu and edited by the course staff*

-----------------

### Topics
1. <a href='#clust'>Introduction</a>
1. <a href='#kmeans'>K-means</a>
  - <a href='#kmeansformulation'>K-means formulation</a>
  - <a href='#kmeansconv'>K-means convergence</a>
  - <a href='#kmeanscomplex'>K-means time complexity</a>
1. <a href='#em'>EM</a>
  - <a href='#model'>Probabilistic Model</a>
  - <a href='#E'>E-step</a>
  - <a href='#M'>M-step</a>
  - <a href='#comp'>Final comparisons between k-means and EM</a>


### <a id='clust'></a>Introduction

  The clustering problem can be posed as follows: we start with several data points of dimension $$p$$. For the single-cell dataset mentioned last lecture, we had $$n = 1,300,000$$ cells and $$p = 2000$$ after some preprocessing. The goal of clustering is to partition the $$n$$ cells based on some notion of distance.

  <div class="fig figcenter fighighlight">
    <img src="/Win2018/assets/lecture18/10x.png" width="80%">
  	<div class="figcaption">The 10x Genomics 1.3 million cell dataset visualized in 2 dimensions using t-stochastic neighbor embedding and colored based on computed clusters.</div>
  </div>

### <a id='kmeans'></a>K-means

We will use K-means to to do clustering of cells. For reference, a typically used dataset of cells we want to cluster based on gene-expression has $n = 1.3$ million data-points and $$p = 28000\ $$ gene expressions.

Before we apply a clustering algorithm, we need to define a distance metric. Since our goal is to cluster points that are close together into the same cluster, this depends on how we represent "closeness" with our distance metric.

Given a set of points $$G = \{y_1, y_2, y_3, ... , y_n\}$ with $y_i \in \mathbb{R} \  \ $$ one typical distance metric is the squared Euclidean distance ($$\ell_2$$ distance):

$$\text{dist}(y_1, y_2) = ||y_1 - y_2||^2$$

And with this distance metric, we will try to group our set of points $$G$$ into $$k$$ clusters (each denoted as $$G_i$$) such that $$G = G_1 \cup G_2 \cup G_3 \cup ... \cup G_k\ $$. We should note that the number of clusters $$k$$ is not always known, especially if we don't have any information about the structure of the data beforehand. Therefore you may want to run an our loop outside of you clustering algorithm that tries to find the optimal number of clusters $$k$$. We won't be discussing this approach in this lecture, and will instead just assume that we are already given a $k$ from our outer loop.

#### <a id='kmeansformulation'></a>K-means formulation

Each cluster (represented by $$G_j$$ with $$j \in \{1,2,3,...,k\}$$) has a
representative point
$$r_j \in \mathbb{R}$$. This can be thought as a center of $$G_j$$. The goal of
clustering is to assign each point $y_i$ to the cluster label $$Z_i$$
with $$Z_i \in \{1,2,3,...,k\}$$. In order to find the most optimal label
assignment for our set of points, we will seek to minimize the error:

$$J = \min_{\{Z_i\}, \{r_j\}}{\sum_{i=1}^{n}||y_i - r_{z_i}||^2}$$

Since we see that there is a double minimization in our error expression, we will do iterations of two steps each to minimize this expression.

Given $$\{r_j\}$$ (centers are fixed), for every $$i$$ we will minimize
$$||y_i - r_j||^2$$ over $$j$$. That is, we will assign $$y_i$$ to the closet cluster.

Given $$\{z_i\}$$ (assignments are fixed), we need to compute center of each
cluster. Let's focus on computing the center of one cluster $$G_j$$,
since we can apply this same procedure to computing the centers of all other
clusters as well. So given cluster $$G_j$$, we will need to minimize
$$\sum_{i\in G_j}||y_i - r||^2$$ over $$r$$ (space of potential
centers/representative points for our cluster $$j$$). We can express the
optimal $$r$$ as

$$r_j^* = \frac{1}{|G_j|}\sum_{i \in G_j}y_i$$

This can be thought of as the center of mass of cluster $$G_j$$.

So in every step we'll do two calculations (one over
$$r$$ and one on $$j$$). The iterative algorithm can be thought of as

$$r_j^{(0)} \rightarrow z_i^{(0)} \rightarrow r_j^{(1)} \rightarrow z_i^{(1)} \rightarrow r_j^{(2)} \rightarrow z_i^{(2)} \rightarrow ...$$

#### <a id='kmeansconv'></a>K-means convergence

One more thing to consider is whether this algorithm will converge. We see that
in the first step, the total cost will go down. In the second step, we will
improve the objective found as well. We also know that the error term is always
positive. So taking all this into account we see that because we are always
reducing the nonzero error, we must converge to a minimum. Another thing to
consider is that the initialization of the centers can greatly affect the rate
of convergence of k-means. There is a good amount of research into how to pick
initial centers more optimally. We should note that just because k-means
converges to a minimum, it doesn't mean that it converges to the global minimum.
There are structures based on the data and initial clusters that might
force k-means to converge to a local minimum.

#### <a id='kmeanscomplex'></a>K-means time complexity

When calculating the time complexity of the k-means computation, we
only really care about the the complexity in terms of the number of points
$$n$$ because $n$ is typically so much larger than all other quantities.

During the first step we see that we need to compute k distances per point, and
this is linear in the number of points. Since we choose to express complexity
in terms of $$n$$, we see that this time complexity is $$O(n)$$. During the
second step we compute the center of mass, which again requires one distance
computation per point. Thus this time complexity is also $$O(n)$$.
The total time complexity for a single step of k-means if $$O(n)$$.

If we change the distance measure from the squared Euclidean distance to some
other measure, does the runtime change from $$O(n)$$ per iteration?
We recall that the formulate we used for calculating
$$r^*$$ is tied to using $$\ell_2$$ as the distance measure.
If you use a different distance measure that requires a more complex
computation of $$r^*$$, then it maybe not just require one distance
computation per point. This might cause the time complexity to different
from $$O(n)$$.

### <a id='em'></a>EM

We can see that k-means might look similar to the EM algorithm we have talked
about in previous lectures. K-means is purely an optimization problem,
while EM requires a probabilistic model.

In applying EM to clustering, we will first need a generative model.
This will provide use with a probability of a connection between clusters and
data-points. More specifically, we can use a Gaussian mixture model. So for
every $$i$$-th datapoint associated with cluster $$Z_i$$ with
$$Z_i \in \{1,2,3,...,k\}\ $$ have define

$$P(Z_i = j) = p_j$$

as the probability that the $$i$$-th point belongs to cluster $$j$$.
Now conditioning on $$Z_i = j$$ we say

$$Y_i \sim \mathcal{N}(\mu, I)$$ where $$\mu$$ is the mean and $$I$$ is the covariance.

So for example if a point comes from $$\mu_1$$, we can generate some random
point centered around $$\mu_1$$ based on the above distribution. Therefore
$$Y_i$$ comes from a mixture of $$k$$ different Gaussians. We should not
that $$Y_i$$ itself is not a pure Gaussian.

#### <a id='model'></a>Probabilistic Model

We define the model of $$y$$ as

$$f(y;\theta)$$

The means are unknown and the $$p_j$$'s are also unknown
(how "big" the cluster $$j$$ is). We can encapsulate the unknown parameters in
$$\theta$$ as

$$\theta = (\mu_1, \mu_2, ... , \mu_k, p_1, p_2, ... , p_k)$$

Now we express a model for a single point $$y_i$$ as

$$f(y;\theta) = \sum_{j}p_j \frac{1}{\left(\sqrt{2\pi}\right)^{p_j}}e^{-\frac{1}{2}||y_i-\mu_j||^2}$$

The summation over $$j$$ with $$p_j$$ comes from the fact that $$y_i$$ comes
from each cluster $$j$$ with probability $$p_j$$. Expressing the log likelihood we get

$$\log(f(y_i;\theta)) = \log\left(\sum_{j}p_j \frac{1}{\left(\sqrt{2\pi}\right)^{p_j}}e^{-\frac{1}{2}||y_i-\mu_j||^2}\right)$$

Now to consider all the datapoints we will need to combine all the above models
for points independently. We want to express $$\log(f(Y, \theta))$$ with
$$Y = \{y_1, y_2, ... , y_n\}$$. To do this we simply sum the log likelihood
expressions given above across all $n$ points. Doing so gives us

$$\log(f(Y, \theta)) = \sum_{i=1}^{n}\log\left(\sum_{j}p_j \frac{1}{\left(\sqrt{2\pi}\right)^{p_j}}e^{-\frac{1}{2}||y_i-\mu_j||^2}\right)$$

We see that this is a non-convex function, so if we try to maximize this
expression over $$\theta$$ we may only reach at local maximum.

#### <a id='E'></a>E-step

For the E-step of the EM algorithm we fix $$\theta$$ to get $$P(Z|Y;\theta)$$
where $$Z$$ represents all cluster labels of all points and $$Y$$
represents all the points, assumed to be independently generated.

$$P(Z|Y;\theta) = \prod_{i=1}^{n}P(Z_i|Y_i;\theta)$$

Now with the posterior (soft EM) we have

$$P(Z_i = j|Y_i; \theta) \propto p_j \frac{e^{-||y_i - \mu_j||^2}}{\sum_{j}p_j e^{-||y_i - \mu_j||^2}}$$

The $$p_i$$ term represents how likely cluster $$j$$ is. The fractional with
the exponentials represented the normalized likelihood of being in a cluster
$$j$$.

#### <a id='M'></a>M-step

In the M-step we want to commuterise (will chek later..)
$$\max_{\theta}E_{Z|Y}[\log(P(Y,Z;\theta))]\ \ $$
when fixing $$P(Z|Y)$$. Now the hidden variable $$Z$$ is exposed and so we no
longer have a summation. Now we write out the full expression.

$$
\begin{align*}
    \max_{\theta}E_{Z|Y}[\log(P(Y,Z;\theta))] &= \max_{\theta}E_{Z|Y}\left[\sum_{i=1}^{n}\log\left(\frac{P_{Z_i}}{(\sqrt{2\pi})^{p_{Z_i}}} e^{-\frac{1}{2}||y_i - u_{Z_i}||^2} \right)]\right] \\
    &= \max_{\theta} E_{Z|Y} \sum_{i=1}^{n}\left( \log(P_{Z_i}) - \frac{1}{2}||y_i - \mu_{Z_i}||^2 - \log\left((\sqrt{2\pi})^{p_{Z_i}}\right) \right)
\end{align*}
$$

Now we can estimate the $$p$$ probabilities. For example if we focus on cluster
1, each point has a probabilistic assignment to it. We can add all them up in
each cluster, and normalize by the total to get our estimates. Additionally
given the probability (soft) of clustering assignments we can say

$$p_s^* \propto \text{total fractional assignments to cluster j summed over all data points}$$

and

$$\mu_j^* = \text{weighted average of $y_i$'s}$$

where the weight is proportional the probability that $$y_i$$ is assigned to
cluster $$j$$.

#### <a id='comp'></a>Final comparisons between k-means and EM

Though k-means and EM are similar, they are not quite equivalent. For the E step
the analogy is assigning points to clusters in k-means, but in EM we assign
points probabilistically to cluster (soft assignment). Another difference
is that we are not making hard assignments in EM, since we are doing soft EM.
