# scBatch
Correct scRNA-seq count matrix subject to batch effects by sample distance matrix correction

scBatch utilizes previous correction on sample distance matrices, such as [QuantNorm](github.com/tengfei-emory/QuantNorm), to further correct the count matrix. For the results generated for the manuscript, the relevant scripts are available at [this repository](github.com/tengfei-emory/scBatch-paper-scripts). We implemented the method with RcppArmadillo for higher efficiency.

# Installation
The package requires R version 3.3.0 with prerequisite package [Rcpp](https://CRAN.R-project.org/package=Rcpp), stats and utils. The package can be installed using the following code. The installation will typically complete within a minute.
```{r}
devtools::install_github('tengfei-emory/scBatch')
library(scBatch)
```

# Toy example
The following script utilizes scBatch to conduct batch effect correction on a simulated scRNA-seq data by [splatter](https://bioconductor.org/packages/release/bioc/html/splatter.html). It took around two minutes to correct a relative small data set with 200 samples and 1000 genes.  

```{r}
library(splatter)
library(scBatch)

#sample size n, number of genes p and the probability of being differentially expressed de
n=200
p=1000
de=0.1

#Simulate a scRNA-seq data with four batches and four biological groups by splatter.
sim.groups <- splatSimulate(nGenes=p, batchCells = c(n/4,n/4,n/4,n/4), seed = 1234, group.prob = c(0.4,0.3,0.2,0.1),
                            de.prob = c(de,de,de,de), method = "groups", verbose = F)
sim.groups <- normalise(sim.groups)

#Extract batch, biological information
batch = sim.groups@colData$Batch
cell.type = sim.groups@colData$Group

#normalized count matrix
exp <- exprs(sim.groups)

#Distance matrix correction by QuantNorm   
correctedD <- QuantNorm(exp,as.numeric(as.factor(batch)),logdat=F,method='row/column',cor_method='pearson',max=5)

#Corrected count based on the corrected distance matrix.
correctedmatrix <-scBatchCpp(c=exp,d=correctedD,w=diag(n),m=5,max=20,tol=1e-10,step=0.0001,derif=scBatch::derif,verbose=T)

```

# References
Fei, Teng, et al. "Mitigating the adverse impact of batch effects in sample pattern detection", Bioinformatics 34(15):2634–2641. (2018).

Zappia, Luke, Belinda Phipson, and Alicia Oshlack. "Splatter: simulation of single-cell RNA sequencing data." Genome biology 18.1 (2017): 174.

Eddelbuettel, Dirk, and Conrad Sanderson. "RcppArmadillo: Accelerating R with high-performance C++ linear algebra." Computational Statistics & Data Analysis 71 (2014): 1054-1063.
