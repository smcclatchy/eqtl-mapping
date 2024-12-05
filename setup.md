---
title: Setup
---

## Software Setup

R is a programming language that is especially powerful for data exploration, 
visualization, and statistical analysis. To interact with R, we use RStudio. 

1. Install the latest version of R from [CRAN](https://cran.r-project.org/).

2. Install the latest version of [RStudio](https://www.rstudio.com/products/rstudio/download/). 
Choose the free RStudio Desktop version for Windows, Mac, or Linux. 

3. Start RStudio. 

4. Install R and Bioconductor packages. 

```r
install.packages(c("tidyverse", "ggbeeswarm", "knitr", "qtl2", "remotes"))

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("AnnotationHub")
BiocManager::install("DeSeq2")
BiocManager::install("rtracklayer")

remotes::install_github("StoreyLab/qvalue")
remotes::install_github("churchill-lab/intermediate")
```

Once the installation is complete, load the libraries to make sure that they 
installed correctly. 

```r
library(tidyverse)
library(ggbeeswarm)
library(knitr)
library(intermediate)
library(qvalue)
library(qtl2)
library(AnnotationHub)
library(DESeq2)
library(rtracklayer)
```

If the libraries don't load and you received errors during the installation,
please contact the workshop instructors before the workshop to help you.

## Project organization

1. Create a new project in your Desktop called `eqtl_mapping`. 
- Click the `File` menu button, then `New Project`.
- Click `New Directory`. 
- Click `New Project`.
- Type `eqtl_mapping` as the directory name. Browse to your Desktop to create the project there.
- Click the `Create Project` button.

2. Use the `Files` tab to create  a `data` folder to hold the data, a `scripts` folder to 
house your scripts, and a `results` folder to hold results. Alternatively, you can use the 
R console to run the following commands for step 2 only. You still need to create a 
project with step 1.

```r
dir.create("./data")
dir.create("./scripts")
dir.create("./results")
```

## Data Sets

For this course, we will have several data files which you will need to 
download to the `data` directory in the project folder on your Desktop.

1. Download files from the 
[Github lesson repository](https://github.com/smcclatchy/eqtl-mapping/tree/main/episodes/data).
You will need to download them one by one using the direct links below. For each
file, locate the download button at upper right.

![Select the download button](../episodes/fig/download-button.png){alt="Graphic showing the download button at right on the Github data file page"}

Repeat this process for each file. Then move the files from wherever your 
downloads go (*e.g.* `Downloads`) to the `data` directory in the `eqtl_mapping` 
project. You can use a graphical user interface (*e.g.* Windows File Explorer, 
Mac Finder) to move the files.

1. [physiological phenotypes](https://github.com/smcclatchy/eqtl-mapping/blob/main/episodes/data/attie_do_pheno.rds)
1. [phenotype dictionary](https://github.com/smcclatchy/eqtl-mapping/blob/main/episodes/data/attie_do_pheno_dict.rds)
1. [covariates](https://github.com/smcclatchy/eqtl-mapping/blob/main/episodes/data/attie_do_covar.rds)
1. [gene annotations](https://github.com/smcclatchy/eqtl-mapping/blob/main/episodes/data/attie_do_expr_annot.rds)
1. [raw gene expression](https://github.com/smcclatchy/eqtl-mapping/blob/main/episodes/data/attie_do_expr_raw.rds)
1. [map](https://github.com/smcclatchy/eqtl-mapping/blob/main/episodes/data/attie_do_map.rds)
1. [kinship](https://github.com/smcclatchy/eqtl-mapping/blob/main/episodes/data/attie_do_kinship.rds)
1. [chromosome 11 insulin blups](https://github.com/smcclatchy/eqtl-mapping/blob/main/episodes/data/ins_tauc_blup_chr11.rds)
1. [chromosome 11 Hnf1b blups](https://github.com/smcclatchy/eqtl-mapping/blob/main/episodes/data/hnf1b_blup_chr11.rds)

2. Copy, paste, and run the following code in the RStudio console to download
the genotype probabilities for the gene expression study we will explore in this
lesson.

```r
download.file(url      = "https://thejacksonlaboratory.box.com/shared/static/4hy4hbjyrxjbrzh570i4g02r62bx3lgk.rds",
              destfile = "data/attie_DO500_genoprobs_v5.rds",
              mode     = "wb")
```
