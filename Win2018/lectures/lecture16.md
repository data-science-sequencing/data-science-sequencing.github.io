# Lecture 15: Introduction to Single Cell RNA-Seq

Thursday 1 March 2018

*scribed by Po-Hsuan Wei*

###Topics
1. <a href='#motivation'>Motivation</a>
2. <a href='#platforms'>Sequencing Platforms</a>
  - <a href='#fluidigm'>Fluidigm C1</a>
  - <a href='#dropseq'>Drop-Seq</a>
3. <a href='#workflow'>Computational Workflow of Drop-Seq</a>
4. <a href='#future'>Outlook</a>


###<a id='type'></a>Motivation
In previous lectures, we discussed bulk RNA-seq, which is done by **sample &rarr; cell lysis &rarr; reverse transcription &rarr; amplification &rarr; sequencing &rarr; read-alignment &rarr; quantification &rarr; transcript population **. Since multiple cells were lysed together, the transcript abundance we get at the final step is a population average of all cells in the sample. This prevents us from differentiating between changes in regulation within cells and compostion of the sample of cells. Neither can we infer the cell differentiation stages, as the cells in a sample have various differentiation rates. 

###<a id='platforms'></a>Single Cell Sequencing Procedure and Platforms

Tang et al. (2009) obtained a single-cell transcript abundance by isolating the cell before lysis. The abundance data from one cell, however, is not sufficient to represent the transcriptome distribution. Since multiple single-cell sequencing is necessary, the time-consuming and expensive procedure carried out by Tang et al. were replaced by platforms that parallelize the preparation process up to the sequencing step. Fluidigm C1, speeds up the isolation process by trapping and releasing single cells into individual microfluidic channels, whereas Drop-Seq used droplet and "tags" to allow for joint downstream processing.

<div class="fig figcenter fighighlight"> <img src="/Win2018/assets/lecture16/abundance.png" width="40%"> <div class="figcaption"> Transcript Abudnace Distribution</div> </div>
</br>
####<a id='fluidigm'></a>Fluidigm C1

The microfluidic chip traps a cell per capture site, and the cells are subsequently flushed to obtain single cell samples per pipe. The cells then undergo cell lysis, reverses transcription and amplification steps within its own chamber. The resulting products are then ready to be sequenced.

<div class="fig figcenter fighighlight"> <img src="/Win2018/assets/lecture16/fluidigm.png" width="80%"> <div class="figcaption"> Fluidigm C1</div> </div>
</br>
####<a id='dropseq'></a>Drop-Seq


<div class="fig figcenter fighighlight"> <img src="/Win2018/assets/lecture16/dropseq2.png" width="90%"> <div class="figcaption"> Drop-Seq</div> </div>
</br>
Drop-Seq uses microfluidics to compartmentalize a barcoded bead, a cell and lysis buffer. Each bead has millions of primer attached, and each primer consists of a PCR handle, a cell barcode, a UMI (unique molecular identifiers) and a long T-tail. The cells are lysed within each droplet, and their mRNAs hybrizes to the long T-tail on the primers. Subsequently, the droplets are broken, and all mRNAs are reverse transcribed, amplified with PCR and are ready for sequencing. The sequencing step is done using paired-end reads, where the first read of the pair always coincides with the cell (barcode + UMI) part of the primer. 

Drop-Seq aims to trace molecules back to each cell, and this is achieved by using a 12bp cell barcode to distinguish different cells, and an 8bp UMI to differentiate molecules within the same cell. UMI is particularly useful for removing the preferential amplification of PCR.

Since expression information is encoded in the cDNA, we don't need to worry about collision of UMIs between different transcripts in the same cell. Still, we can calculate the expected value of collision:


N = # of molecules in the same cell ( approx. 100k)
B = # of distinct UMIS ($$4^8$$, approx. 64k) 
$$X_i$$ = The UMI of the ith cell
$$P(X_i = k, X_j=k)  = 1/B^2$$
$$P(X_i = X_j)  = 1/B$$
$$E[# of collisions] = N(N-1)/2B = 77.3k $$
</br>

<div class="fig figcenter fighighlight"> <img src="/Win2018/assets/lecture16/splitnpool.jpg" width="70%"> <div class="figcaption"> Split-and-Pool for Generating Random Cell Barcodes</div> </div>
</br>

<div class="fig figcenter fighighlight"> <img src="/Win2018/assets/lecture16/biasedreads.png" width="60%"> <div class="figcaption">Biased Paired Reads</div> </div>
</br>

###<a id='workflow'></a>Computational Workflow of Drop-Seq

Now we have obtained reads consisting of cell barcode, UMI and cDNA, we can estimante the transcript abundances. We first group reads by cell barcode, align cDNA reads, and then count unique molecules per cell per gene using UMI.

<div class="fig figcenter fighighlight"> <img src="/Win2018/assets/lecture16/biasedreads2.png" width="40%"> <div class="figcaption">Recovering Cells by Barcode Grouping</div> </div>
</br>

$$\text{The final output is a gene expression table, whose columns are cells, and rows are genes. There are 2 major sources of errors. An error can occur during sequencing and synthesis of barcodes, or each droplet can accidentally capture more than one cell. If the error occurs during sequencing, since the barcode is short, we can simply filter out the erroneous barcodes by choosing the ones with many read counts. For the case where a droplet captures more than one cell, we can model droplet arrivals with Poisson distribution. The parameter }\lambda\text{ is controlled by the concentration of the cells, all else being the same. There is a tradeoff between singlet encapsulation rate and mutiplet encapsulation rate. Different techniques can be adopted to detect doublets. If 2 cells with similar mRNA counts trapped together, then the read count effectively doubles for reads with the same cell barcode. If the 2 cells in a doublet have very different mRNA counts, then the support would be longer than either of the cells.}$$


<div class="fig figcenter fighighlight"> <img src="/Win2018/assets/lecture16/tradeoff.png" width="40%"> <div class="figcaption"> </div>Choosing a Desgin Point </div>
</br>
###<a id='future'></a>Outlook

While the sequencing cost is decreasing, the number of single-cell studies increases exponentially. To put the cost into perspective , suppose we sequence 1M cells with 50k reads per cell, 0.00001 USD per read, then the overall cost is 500k.
<div class="fig figcenter fighighlight"> <img src="/Win2018/assets/lecture16/future.png" width="80%"> <div class="figcaption"> </div> </div>
</br>