---
title: "Mapping Many Gene Expression Traits"
teaching: 30
exercises: 30
---

:::::::::::::::::::::::::::::::::::::: questions 

- How do I map many genes?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- To map several genes at the same time

::::::::::::::::::::::::::::::::::::::::::::::::

<!-- DMG: Notes from sumner run. 11 hours with 20 cores & 100GB of memory.
11.4 G output file. 69005 rows x 21771 columns. -->

### Load Libraries  

Load the libraries and source two other R scripts.


``` r
library(tidyverse)
library(knitr)
library(broom)
library(qtl2)
library(qtl2ggplot)
library(RColorBrewer)
library(AnnotationHub)
library(rtracklayer)
#source("../code/gg_transcriptome_map.R")
#source("../code/qtl_heatmap.R")
```

### Load Data


``` r
# expression data
load("../data/attie_DO500_expr.datasets.RData")

# data from paper
load("../data/dataset.islet.rnaseq.RData")

# phenotypes
load("../data/attie_DO500_clinical.phenotypes.RData")

# mapping data
load("../data/attie_DO500_mapping.data.RData")

# genotype probabilities
probs = readRDS("../data/attie_DO500_genoprobs_v5.rds")
```

### Data Selection

For this lesson, lets choose a random set of 50 gene expression phenotypes.




``` r
genes = colnames(norm)

sams <- sample(length(genes), 50, replace = FALSE, prob = NULL)
genes <- genes[sams]

gene.info <- dataset.islet.rnaseq$annots[genes,]
rownames(gene.info) = NULL
kable(gene.info[1:10,])
```

### Expression Data

Lets check the distribution for the first 20 gene expression phenotypes. If you 
would like to check the distribution of all 50 genes, change 
`for(gene in genes[1:20])` in the code below to `for(gene in genes)`.


``` r
par(mfrow=c(3,4))
for(gene in genes[1:20]){
  hist(norm[,gene], main = gene)
  }
```

Check the distributions.  Do they all have a normal distribution?

You will notice that the distribution of some genes are skewed to the left. This 
means that that only a small amount of samples have data and therefore, will 
need to be removed.  A suitable qc would be keeping expression data that have at 
least 5% of the samples with more than 10 reads.


``` r
genes_qc <- which(as.numeric(colSums(counts[ , genes] > 10)) >= 0.05 * nrow(counts[,genes]))
genes <- genes[genes_qc]
```

### Genotype probabilities  

We have explored this earlier in th previous [lesson](https://smcclatchy.github.io/eqtl-mapping/review-mapping-steps/index.html#genotype-probabilities).  But, as a reminder, we have already calculated genotype 
probabilities which we loaded above called `probs`.  This contains the 8 state 
genotype probabilities using the 69k grid  map of the same 500 DO mice that also 
have clinical phenotypes. 

### Covariates    

Now let's add the necessary covariates. For these 50 gene expression data, we 
will correct for `DOwave`,`sex` and `diet_days`.





``` r
# convert sex and DO wave (batch) to factors
pheno_clin$sex = factor(pheno_clin$sex)
pheno_clin$DOwave = factor(pheno_clin$DOwave)
pheno_clin$diet_days = factor(pheno_clin$DOwave)

covar = model.matrix(~sex + DOwave + diet_days, data = pheno_clin)[,-1]
```

### QTL Scans


``` r
qtl.file = "../results/gene.norm_qtl_cis.trans.Rdata"

if(file.exists(qtl.file)) {
  load(qtl.file)
  } else {
    qtl = scan1(genoprobs = probs, 
                pheno = norm[,genes, drop = FALSE],
                kinship = K, 
                addcovar = covar, 
                cores = 2)
    save(qtl, file = qtl.file)
    }
```

### QTL plots

Let's plot the first 20 gene expression phenotypes.  If you would like to plot 
all 50, change `for(i in 1:20)` in the code below to `for(i in 1:ncol(qtl))`.




``` r
par(mfrow=c(3,4))
for(i in 1:20) {
  plot_scan1(x = qtl, 
             map = map, 
             lodcolumn = i, 
             main = colnames(qtl)[i])
  abline(h = 6, col = 2, lwd = 2)
  }
```

### QTL Peaks

We are also going to save our peak results so we can use these again else where.  
First, lets get out peaks with a LOD score greater than 6. 


``` r
lod_threshold = 6
peaks = find_peaks(scan1_output = qtl, 
                   map = map, 
                   threshold = lod_threshold, 
                   peakdrop = 4, 
                   prob = 0.95)
```

We will save these peaks into a csv file. 


``` r
kable(peaks[1:10,] %>% 
        dplyr::select(-lodindex) %>% 
        arrange(chr, pos), caption = "Expression QTL (eQTL) Peaks with LOD >= 6")

# write_csv(peaks, "../results/gene.norm_qtl_peaks_cis.trans.csv")
```



### QTL Peaks Figure


``` r
qtl_heatmap(qtl = qtl, map = map, low.thr = 3.5)
```

::::::::::::::::::::::::::::::::::::: challenge 

## Challenge 1: 


What do the qtl scans for all gene exression traits look like? *Note:* Don't 
worry, we've done the qtl scans for you!!!
You can read in this file, `../data/gene.norm_qtl_all.genes.Rdata`, which are 
the `scan1` results for all gene expression traits. 

:::::::::::::::::::::::: solution 


``` r
load("../data/gene.norm_qtl_all.genes.Rdata")

lod_threshold = 6
peaks = find_peaks(scan1_output = qtl.all, 
                map = map, 
                threshold = lod_threshold, 
                peakdrop = 4, 
                prob = 0.95)
write_csv(peaks, "../results/gene.norm_qtl_all.genes_peaks.csv")

## Heat Map
qtl_heatmap(qtl = qtl, map = map, low.thr = 3.5)
```

:::::::::::::::::::::::::::::::::


::::::::::::::::::::::::::::::::::::::::::::::::


::::::::::::::::::::::::::::::::::::: keypoints 

- 

::::::::::::::::::::::::::::::::::::::::::::::::
