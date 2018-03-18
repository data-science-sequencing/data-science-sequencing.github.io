---
layout: page
mathjax: true
permalink: /Win2018/assignments/assignment2/
---
## Assignment 2

posted Tuesday 6 February 2018  
due Tuesday 20 February 2018 at 11:59pm

[Solutions](/Win2018/assets/assignment2/ee372_assignment2_solutions/)

**Submission policy**: Report all plots and your code in [this Jupyter notebook](/Win2018/assets/assignment2/ee372_assignment2.ipynb). Print your notebook as a PDF and attach it to the rest of your assignment. Turn in your assignment through [Gradescope](https://gradescope.com/).

### Question I: Lander-Waterman and repeat statistics

In Question III of the last problem set, you aligned reads from an _E. coli_ sequencing experiment to the _E. coli_ genome using Bowtie2, obtaining 627862 aligned reads (out of 644022 total reads). Download the set of aligned reads [here](/Win2018/assets/assignment2/ecoli_aligned_reads.bam) (the alignments are saved as a \*.bam file). You will be using SAMtools for this problem. You can install SAMtools on the Stanford Rice cluster using the following commands:

```
wget http://data-science-sequencing.github.io/Win2018/assets/assignment2/install_samtools.sh
bash install_samtools.sh <SUNETID>
```

1. You can view the first 15 lines of the file using the command ```samtools view ecoli_aligned_reads.bam | head -15```. Refer to the documentation for the \*.sam file format [here](https://samtools.github.io/hts-specs/SAMv1.pdf). Convert the \*.bam file to a \*.sam file using the ```samtools view``` command. What position does the first read align to?
2. The length of the _E. coli_ genome is 4639675. Write a function to compute the proportion of the genome by covered by the reads in a given \*.sam file. What proportion of the genome is covered by the reads in ecoli_aligned_reads.bam?
3. Subsample the \*.sam file using the ```samtools view -s``` command. You should generate a new \*.sam file for each subsampling. Using the function you wrote above, plot the proportion of the genome covered as a function of $$p$$, the proportion of reads kept. Interpret what you see (1-2 sentences).
4. Compute the triple repeat and interleaved statistics of the _E. coli_ genome (available [here](http://portal.nersc.gov/dna/microbial/assembly/uploads/dtse/Mock-Community/E.coli_K12_ATCC_700926.fasta)). Report the number of triple and interleaved repeats of length more than 200. Consider only the forward strand of the reference. _Hint_: Use the software [MUMmer](http://mummer.sourceforge.net/manual/). The ```repeat-match``` command might be helpful. You can install MUMmer on the Stanford Rice cluster using the following commands:

```
wget http://data-science-sequencing.github.io/Win2018/assets/assignment2/install_mummer.sh
bash install_mummer.sh
```


You may also use [pysam](http://pysam.readthedocs.io/en/latest/api.html) for this problem.

### Question II: de Bruijn graph assembly

1. Give the 5-mer spectrum of TAAAAACCCCAAAAAG. How many different assemblies are consistent with the 5-mer spectrum?

2. The support of a $$k$$-mer spectrum is the $$k$$-mer spectrum with the value of all non-zero $$k$$-mers set to 1. Give the assembly of TAAAAACCCCAAAAAG from the support of its 5-mer spectrum. How many different assemblies are consistent with the support of this 5-mer spectrum?

3. Study the implementation of the de Bruijn graph assembler by Ben Langmead [here](http://nbviewer.jupyter.org/github/BenLangmead/comp-genomics-class/blob/master/notebooks/CG_deBruijn.ipynb). You should copy and paste the code from the top cell into your notebook as you will use a tweak of this  class to perform assembly. (You will need to make sure that the assembly is done using the k-mers of non-zero support. That is, even if a k-mer is seen twice one would add only an edge between the k-1 mers.) Note that you will need to pass a list of reads (strings) as the ```strIter``` argument when initializing an instance of the class (see ```__iter___```). You can use the ```eulerianWalkOrCycle``` method to obtain a list of $$k-1$$-mers corresponding to an Eulerian walk through the graph. Write a function that obtains the assembly from this list of $$k-1$$-mers (i.e. if the list is ['ABCD','BCDE'] with $$k=5$$, then your function should return 'ABCDE').

4. Write a function to generate random reads. The input should be the number of reads generated $$N$$ and the length $$L$$ of each read generated. The output should be $$N$$ random length-$$L$$ sequences of nucleotides. Generate a random length-100 genome.

5. Write a function to sample reads from a genome. The input should be the genome, the number of reads generated $$N$$, and the length $$L$$ of each read generated. Assuming that $$L = 10$$, how many reads do you need to achieve a coverage depth of 30? Generate this number of reads and give the assembly using your code from part 3.

6. Write a modified version of the previous function for sampling reads from a genome with error. Generate random length-10 reads with 5% error rate and a coverage of 30. Give the assembly using your code for part 3. What do you observe? You may want to rerun your code for 4-6 to make sure your observations are consistent.

<!-- ### Question III: Alignment of random strings

1. Write a function that computes the edit distance between two strings using dynamic programming. You can use standard libraries such as  [editdistance](https://pypi.python.org/pypi/editdistance) to check the correctness of your function. The input should be two strings each of at least length 2. The output should be a float representing the edit distance between the two input strings.

2. How can you use the above function to perform alignment between two reads? If you have $$N$$ total reads each of length $$L$$, what is the runtime of this approach in terms of $$N$$ and $$L$$?

3. For $$L=10$$, generate two reads of length $$L$$ and compute their edit distance. Average your result over 100 runs to obtain an estimate of the average edit distance $$\hat{d}_L$$ of two randomly generated reads of length $$L$$. Repeat for $$L = 1, 100, 1000, 10000$$. Plot $$\hat{d}_L$$ as a function of $$L$$ with error bars. Do you observed any trends? What can you say about how well random strings align? (2-3 sentences). -->

### Question III: Viterbi Algorithm

Consider the following state diagram for a simplified Nanopore sequencer only sequencing two bases 0 and 1.

<div class="fig figcenter fighighlight">
  <img src="/Win2018/assets/lecture4/lecture4-figure1.png" width="40%">
</div>

1. The state space here is  $$\{00, 01, 10, 11\} \$$. What can you say about the
inter-symbol interference? In other words, how many nucleotides are in the pore at each sampling?

2. Consider the output sequence $$1.34, 0.23, 1.45, 0.5, 0.11 $$. Using the above state diagram, compute the sequence obtained by Viterbi decoding. Assume that each observed $$y_i$$ is a weighted linear combination of the symbols in the context plus some noise: $$y_i = 0.9 s_i + 0.25 s_{i-1} + n_i$$. Assume that the added noise is Gaussian with zero mean.

3. Simulate several length-20 $$s$$. Use these to simulate $$y$$ where $$n$$ is sampled from a $$\mathcal{N}(0, \sigma^2)$$ distribution. Decode $$y$$ using the your Viterbi algorithm. Report the average edit distance between your recovered $$s$$ and the actual $$s$$ as a function of $$\sigma$$. Test for $$\sigma = $$ 1e-3, 1e-2, 1e-1, 1e0, 1e1. You can use the [editdistance](https://github.com/aflc/editdistance) Python package to compute edit distance.

4. Repeat 3 except randomly delete a base in $$s$$ before generating your $$y$$. What do you observe? From a high level, what might you need to change to fix this?
