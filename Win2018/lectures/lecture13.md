---
layout: page
mathjax: true
permalink: /Win2018/lectures/lecture13/
---
## Lecture 13: RNA-seq - Quantification and the EM Algorithm

Tuesday 20 February 2018

_scribed by Jacob Blum, adapted from Paraskevas Deligiannis' [scribed notes](/Spr2016/lectures/lecture15/), and edited by the course staff_

-----------------

## Topics

1.	<a href='#equal'>Algorithm for equal length transcripts</a>
1.	<a href='#ml'>Maximum Likelihood estimation</a>
    - <a href='#convex'>Convexity</a>
    - <a href='#unique'>Special case: Uniquely-mapped reads</a>
1.	<a href='#em'>Expectation-Maximization (EM) algorithm</a>
    - <a href='#optimality'>Convergence & Optimality</a>
    - <a href='#steps'>Iterative steps</a>
    - <a href='#application'>Application to our problem</a>


### <a id='equal'></a> Algorithm for equal length transcripts

For simplicity, we first we will assume that every transcript is of equal length. Let $$S_i$$ represent the set of transcripts that contain read $$r_i$$, and let $$f_{ik}$$ represent the fraction of read $$i$$ that we allocate to some $$k \in S_i$$. Our quantification algorithm will be performed in two steps:

1. Initialize abundance estimates for each transcript:
  $$\rho_1^0, \dots , \rho_K^0$$

2. For step $$m=1,2,...$$ repeat:  
- For $$i=1,2,..,N$$, split $$R_i$$ into fractional counts for each transcript $$k \in S_i$$, equal to the relative abundances of the transcripts in $$S_i$$. This splitting is done as follows:  
$$
f_{ik}^{(m)}=\begin{cases}
\frac{\rho_k^{(m)}}{\sum_{j \in S_i}{\rho_j^{(m)}}} &\text{if }k \in S_i \\
0 & \text{otherwise}  \end{cases}
$$  
- The updated estimate is  
$$
\hat{\rho}_k^{(m+1)}=\frac{1}{N}\sum_{i=1}^{N}{f_{ik}^{(m)}}
$$

Please refer to [lecture 12](/Win2018/lectures/lecture12/#eg2) for an example where transcript lengths are not equal and one transcript is a substring of the other.

### <a id='ml'></a>Maximum Likelihood estimation

Let $$r_1, \dots, r_N$$ be our collection of reads (data).
The distribution of the reads depends on $$\mathbf{\rho}$$.
Since the reads are sampled independently, we have that

$$
Pr(R_1=r_1,R_2=r_2,...,R_N=r_N;\mathbf{\rho}) = \prod_{i=1}^N{Pr(R_i=r_i;\mathbf{\rho})}.
$$

We note that $$R_i=r_i$$ if and only if the read $$R_i$$ maps to a transcript
$$t_k$$ that contains the sequence $$r_i$$ and $$R_i$$ starts from the same position as
$$r_i$$. (Recall our assumption that
  _each read can come from at most one location on each transcript_.)
Therefore, the probability

$$
\begin{equation}
\begin{split}
Pr(R_i=r_i;\mathbf{\rho}) = \sum_{k=1}^K&[{Pr(R_i\text{ maps to } t_k;\mathbf{\rho}) \times Pr(r_i\text{ belongs to } t_k) \times ...}  \\
&...Pr(R_i\text{ starts at the position where } r_i\text{ starts in } t_k)].
\end{split}
\end{equation}
$$

The reads are sampled uniformly, so $$Pr(R_i\text{ maps to } t_k;\mathbf{\rho})=\rho_k$$.
Also, each transcript $$t_k$$ has length $$l$$ and each read of length L coming from $$t_k$$
has an equal chance of starting at each of the possible $$l-L+1$$
positions, hence

$$
Pr(R_i\text{ starts at the position where } r_i\text{ starts in } t_k)=\frac{1}{l-L+1}.
$$

Finally, let

$$
y_{ik}=\left\{ \begin{array}[cc]\\
1 &\text{if }r_i \in t_k\\
0 & \text{otherwise.}  \end{array}\right.
$$

As a result,

$$
Pr(R_i=r_i;\mathbf{\rho}) = \sum_{k=1}^K{\frac{y_{ik}}{l-L+1}\rho_k},
$$

which means that

$$
\begin{equation}
\begin{split}
\log{Pr(R_1=r_1,R_2=r_2,...,R_N=r_N;\mathbf{\rho})} &= \sum_{i=1}^N{\log{\left(\sum_{k=1}^K{\frac{y_{ik}}{l-L+1}\rho_k}\right)}} \\
  &= \sum_{i=1}^N{\log{\left(\sum_{k=1}^K{y_{ik}\rho_k}\right)}}-N\log{(l-L+1)}.
\end{split}
\end{equation}
$$

Maximizing the likelihood function is equivalent to maximizing the log-likelihood,
and since $$-N\log{(l-L+1)}$$ is a constant, the ML estimation for $$k=1,...,K$$ is

$$
\DeclareMathOperator*{\argmax}{arg\!\max}
\hat{\rho}_{k}^{ML}=\argmax_{\rho_k \ge 0\ :\ \sum_{k=1}^K{\rho_k}=1}{\left(\sum_{i=1}^N{\log \left( \sum_{k=1}^K{y_{ik}\rho_k} \right)}\right)}.
$$

If we define the matrix $$Y \in \{0,1\}^{N \times K}$$ with elements $$y_{ik}$$, then

$$
\mathbf{\hat{\rho}^{ML}}=\argmax_{\mathbf{\rho} \succeq 0\ :\ \|\mathbf{\rho}\|_1=1} {\sum_{i=1}^N{\log{(Y\mathbf{\rho})}}}.
$$

Along each row $$i$$ of matrix $$Y$$ the positions with value 1 represent the transcripts to which $$R_i$$ maps. Along each column $$k$$ of $$Y$$, the positions with value 1 represent the reads that map to $$t_k$$.

As we can see, the effect of the data (reads) on the ML estimate can be summarized by matrix $$Y$$, or equivalently, its entries $$y_{ik}$$. Therefore, $$Y$$ is a [_sufficient statistic_](https://en.wikipedia.org/wiki/Sufficient_statistic) for this ML inference problem.

#### <a id='convex'></a>Convexity

The above maximization problem is a convex optimization problem, which is desired, since convex problems are well studied and there exist tools to solve them efficiently.

The convexity arises from the combination of the following factors:  

 - it involves the _maximization of a concave objective_ function, as  

     - $$\sum_{k=1}^K{y_{ik}\rho_k}$$ is an affine function of $$\mathbf{\rho}$$  

     - $$\log(\cdot)$$ is concave  

     - the composition of a concave with an affine function is concave  

     - the sum of concave functions is concave  

 - both $$\rho_k \ge 0$$ and $$\sum_{k=1}^K{\rho_k}=1$$ are _linear constraints_

#### <a id='unique'></a>Special case: Uniquely-mapped reads

Let's consider the special case when there are no multiply-mapped reads, i.e. every read maps to at  most one transcript, or equivalently, each row of matrix $$Y$$ has at most one element with value 1. If we denote by $$N_k$$ the number of reads that map to transcript $$t_k$$, then the above maximization problem becomes

$$
\hat{\rho}_{k}^{ML\ unique}=\argmax_{\rho_k \ge 0\ :\ \sum_{k=1}^K{\rho_k}=1}{(\sum_{k=1}^K{N_k\log{\rho_k}})}.
$$

This problem is easy to solve by equating to 0 the gradient of the objective function w.r.t $$\mathbf{\rho}$$ and using [Lagrange multipliers](https://en.wikipedia.org/wiki/Lagrange_multiplier). The result is obvious and reassuring:

$$
\hat{\rho}_{k}^{ML\ unique}=\frac{N_k}{\sum_{k=1}^K{N_k}}
$$

After all, we knew that the naive approach with split reads works when all reads are uniquely mapped. Note that when all reads map to exactly one transcript, $$\sum_{k=1}^K{N_k}=N$$.

### <a id='em'></a>Expectation-Maximization (EM) algorithm

The [EM algorithm](https://en.wikipedia.org/wiki/Expectation%E2%80%93maximization_algorithm) was thoroughly described in a 1977 [article](https://www.jstor.org/stable/2984875) written by  [Arthur Dempster](https://en.wikipedia.org/wiki/Arthur_P._Dempster), [Nan Laird](https://en.wikipedia.org/wiki/Nan_Laird), and [Donald Rubin](https://en.wikipedia.org/wiki/Donald_Rubin). It provides an iterative method to compute -- locally optimal -- ML parameters of a statistical model. If $$\mathbf{\theta}$$ denotes the unknown parameters to be estimated (in our case, $$\mathbf{\rho}$$) and $$\mathbf{x}$$ are the observed variables (data, in our case, the reads), then the EM algorithm computes

$$
\argmax_{\mathbf{\theta}}{\ Pr(\mathbf{x};\mathbf{\theta})}.
$$

The basic observation is that in many cases of statistical inference, the problem would become significantly easier, or even trivial, if we knew some unobserved _hidden_ variable(s) $$Z$$. For instance, in our case, if we somehow knew for every read which transcript it came from, then our problem would reduce to the simple problem with uniquely mapped reads. This information can be represented by a matrix $$Z$$ such that

$$
Z_{ik}=\left\{ \begin{array}[cc]\\
1 &\text{if }R_i \text{ comes from } t_k\\
0 & \text{otherwise.}  \end{array}\right.
$$

The idea behind EM is to iteratively try to use $$\mathbf{x}$$ and the current estimate $$\mathbf{\hat{\theta}}$$ to
compute the posterior distribution: $$p(Z \mid \mathbf{x} ,\mathbf{\hat{\theta}}) $$.
Then we  update our estimate $$\mathbf{\hat{\theta}}$$ to maximize the expected
value of the likelihood $$p(\mathbf{x},Z \mid \mathbf{\theta})$$ w.r.t. $$p(Z\mid\mathbf{x};\mathbf{\hat{\theta}})$$ and so on:

$$
\mathbf{\hat{\theta}^{(0)}} \to p(Z|\mathbf{x};\mathbf{\hat{\theta}}^{(0)}),\\
\max_{\mathbf{\theta}}\mathbb{E}[\log{p(\mathbf{x},Z|\mathbf{\theta})}] \to \mathbf{\hat{\theta}^{(1)}},\\ ...
$$

#### <a id='optimality'></a>Convergence & Optimality

It is proven that the EM algorithm  _always_ improves the likelihood $$Pr(\mathbf{x};\mathbf{\hat{\theta}})$$ at every iteration. Since $$Pr(\mathbf{x};\mathbf{\theta})$$ is bounded above by 1, this means that the EM algorithm converges to a local maximum.

However, we saw from our ML analysis in a previous section that the ML optimization problem is concave. This means that there is only one local maximum, the global one. Thus, the EM algorithm is guaranteed to find the optimal ML estimate.

#### <a id='steps'></a>Iterative steps

Let's formally describe the EM algorithm. The algorithm consists of two successive steps, the expectation step (E step) and the maximization step (M step), after which it is named. These steps are repeated until the algorithm converges.

For iteration $$m = 0,1,2,...$$  

1. __E step:__ Compute the expected value of the log-likelihood function w.r.t the conditional distribution of $$Z|\mathbf{x}$$ using the current estimate of the unknown parameters $$\mathbf{\hat{\theta}^{(m)}}$$:  
$$
J(\mathbf{\theta};\mathbf{\hat{\theta}^{(m)}}) = \mathbb{E}_{Z|\mathbf{x};\mathbf{\hat{\theta}^{(m)}}}[\log{p(\mathbf{x},Z|\mathbf{\theta})}]
$$  

2. __M step:__ Find the updated estimate $$\mathbf{\hat{\theta}}$$ that maximizes $$J$$:  
$$
\mathbf{\hat{\theta}^{(m+1)}} = \argmax_{\mathbf{\theta}}{J(\mathbf{\theta};\mathbf{\hat{\theta}^{(m)}})}
$$

#### <a id='application'></a>Application to our problem

In our case, it is easy to see that the probability that a read $$R_i \to S_i$$ comes from the transcript $$t_k$$ is proportional to the relative abundance of $$t_k$$ in $$S_i$$. Namely, given our reads $$\mathbf{R}$$ and our current estimate $$\mathbf{\rho}^{(m)}$$,

$$
Pr(Z_{ik}=1|\mathbf{R},\mathbf{\rho}^{(m)}) = f_{ik}^{(m)}.
$$

We have that  

$$
\begin{split}
p(\mathbf{R},Z|\mathbf{\rho}) &= \prod_{k=1}^K \prod_{i=1}^N \left( \rho_k^{Z_{ik}} \right) \mathbb{I} (\mathbf{R} \text{ is consistant with } Z) \\
\implies  \log p(\mathbf{R},Z|\mathbf{\rho}) &= \sum_{k=1}^K \left( \sum_{i=1}^NZ_{ik} \right)\log \rho_k \\
  &= \sum_{k=1}^K{N_k\log{\rho_k}},
\end{split}
$$  

where in the first step, we note that since $$Z$$ represents the ground truth, there is at most one non-zero term in the sum $$\sum_{k=1}^K{Z_{ik}}$$ (similar to the case with uniquely mapped reads). We drop the indicator function in the second step assuming that $$\mathbf{R}$$ and Z are consistant.

Consequently,  

$$
\begin{align}
\mathbf{\hat{\rho}^{(m+1)}}
&= \argmax_{\mathbf{\rho}}{J(\mathbf{\rho};\mathbf{\hat{\rho}^{(m)}})}, \\
&= \argmax_{\mathbf{\rho}}{\mathbb{E}_{Z|\mathbf{R};\mathbf{\hat{\rho}^{(m)}}}[\log{p(\mathbf{R},Z|\mathbf{\rho})}]},\\
&= \argmax_{\mathbf{\rho}}{\sum_{k=1}^K{(\sum_{i=1}^N{f_{ik}^{(m)})\log{\rho_k}}}} \\
\implies \mathbf{\hat{\rho}^{(m+1)}_k} &= \frac{1}{N}\sum_{i=1}^N{f_{ik}^{(m)}},
\end{align}
$$  


where we have omitted the linear constraints on $$\mathbf{\rho}$$ for simplicity. One can verify that if we define $$\hat{N}_k = \sum_{i=1}^N{f_{ik}^{(m)}}$$, the above results follow readily from our analysis for the case of uniquely mapped reads. This is, after all, the reason why we use the _hidden_ variable $$Z$$.
