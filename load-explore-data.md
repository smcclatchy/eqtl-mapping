---
title: "Load and explore the data"
teaching: 15
exercises: 30
---

:::::::::::::::::::::::::::::::::::::: questions 

- What data are required for eqtl mapping?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- To provide an example and exploration of data used for eqtl mapping.

::::::::::::::::::::::::::::::::::::::::::::::::

Load the libraries.

``` r
library(ggbeeswarm)
library(tidyverse)
library(knitr)
library(corrplot)
# the following analysis is derived from supplementary 
# File S1 Attie_eQTL_paper_physiology.Rmd 
# by Daniel Gatti. See Data Dryad entry for more information.
```

## Physiological Phenotypes

The complete data used in these analyses are available from 
[Data Dryad](https://doi.org/10.5061/dryad.pj105). 

Load in the clinical phenotypes.


``` r
# load the data
load("../data/attie_DO500_clinical.phenotypes.RData")
```

``` warning
Warning in readChar(con, 5L, useBytes = TRUE): cannot open compressed file
'../data/attie_DO500_clinical.phenotypes.RData', probable reason 'No such file
or directory'
```

``` error
Error in readChar(con, 5L, useBytes = TRUE): cannot open the connection
```

See the [data dictionary](../data/Attie-232_Attie_DO_Islets-dictionary.csv) to 
see a description of each of these phenotypes. You can also view a table of
the data dictionary.


``` r
pheno_clin_dict %>% 
  select(description, formula) %>% 
  kable()
```

``` error
Error: object 'pheno_clin_dict' not found
```

### Phenotype Distributions

Boxplots are a great way to view the distribution of the data and to identify 
any outliers. We will be using the total area under the curve of insulin from 
the glucose tolerance test (Ins_tAUC). We will also log-transform the data 
using the [scale_y_log10()][scale_y_log10] function.


``` r
# plot Insulin on a log 10 scale
ggplot(pheno_clin, aes(sex, Ins_tAUC)) +
  geom_boxplot() +
  scale_y_log10() +
  labs(title = "Insulin tAUC", y = "Insulin tAUC")
```

``` error
Error: object 'pheno_clin' not found
```

Another visualization that has become popular is the 
[Violin Plot][https://en.wikipedia.org/wiki/Violin_plot]. We can create one 
using ggplot's 
[geom_violin][https://ggplot2.tidyverse.org/reference/geom_violin.html].
Whereas the boxplot automatically adds the median, we must tell `geom_violin()`
which quantiles that we want to draw using the argument 
`draw_quantiles = c(0.25, 0.5, 0.75)`. We have also overlaid the data points 
using ggbeeswarm's
[geom_beeswarm][https://www.rdocumentation.org/packages/ggbeeswarm/versions/0.5.3/topics/geom_beeswarm].
We have told `geom_beeswarm()` to plot the points using the argument 
`alpha = 0.1`. The `alpha` argument ranges between 0 (completely transparent) to
1 (completely opaque). A value of 0.1 means mostly transparent.


``` r
# plot Insulin on a log 10 scale
ggplot(pheno_clin, aes(sex, Ins_tAUC)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_beeswarm(alpha = 0.1) +
  scale_y_log10() +
  labs(title = "Insulin tAUC", y = "Insulin tAUC")
```

``` error
Error: object 'pheno_clin' not found
```

::::::::::::::::::::::::::::::::::::: challenge 

## Challenge 1: 

How many orders of magnitude (powers of 10) does Insulin tAUC span?


:::::::::::::::::::::::: solution 

Insulin tAUC spans three orders of magnitude, from near 10 to over 1000.

:::::::::::::::::::::::::::::::::


## Challenge 2: 

Which sex has higher median Insulin tAUC values?

:::::::::::::::::::::::: solution 

Males have higher Insulin tAUC than females.  

:::::::::::::::::::::::::::::::::
::::::::::::::::::::::::::::::::::::::::::::::::

Both of the boxplot and the violin plot are useful visualizations which you can
use to get some sense of the distribution of your data.

### Quality Control of Data

Many statistical tests rely upon the data having a "normal" (or Gaussian) 
distribution. Many biological phenotypes do not follow this distribution and
must be transformed before analysis. This is why we log-transformed the data
in the plots above. 

While we can "eyeball" the distributions in the violin plot, it would be 
better to use a "quantile-quantile" plot. 


``` r
pheno_clin %>% 
  ggplot(aes(sample = Ins_tAUC)) +
    stat_qq() +
    geom_qq_line() +
    facet_wrap(~sex)
```

``` error
Error: object 'pheno_clin' not found
```

In these plots, the "quantiles" of the normal distribution are plotted on the
X-axis and the data are plotted on the Y-axis. The line indicates the 
quantiles that would be followed by a normal distribution. The untransformed
data do **not** follow a normal distribution because the points are far from
the line.  

Next, we will loag-transform the data and then create a quantile-quantile plot.


``` r
pheno_clin %>% 
  mutate(Ins_tAUC = log(Ins_tAUC)) %>% 
  ggplot(aes(sample = Ins_tAUC)) +
    stat_qq() +
    geom_qq_line() +
    facet_wrap(~sex)
```

``` error
Error: object 'pheno_clin' not found
```

::::::::::::::::::::::::::::::::::::: challenge 

## Challenge 3: 

Does the log transformation make the data more normally distributed? Explain 
your answer.

:::::::::::::::::::::::: solution 

Yes. The log transformation makes the data more normally distributed because the 
data points follow the normality line more closely. 

:::::::::::::::::::::::::::::::::


## Challenge 4: 

Do any data points look suspicious to you? Explain your answer.

:::::::::::::::::::::::: solution 

The data points that deviate from the normality line would be worth 
investigating. All data deviates somewhat from normality, but the three lowest 
points in the male data plot would be worth investigating. They may be real, but 
there may also have been mishap in the assay.

:::::::::::::::::::::::::::::::::
::::::::::::::::::::::::::::::::::::::::::::::::

Another way to identify outliers is to standardize the data and look for data 
points that are more than four standard deviations from the mean.

To do this, we will log transform and standardize Insulin tAUC. 


``` r
ins_tauc = pheno_clin %>% 
             select(mouse, sex, Ins_tAUC) %>%
             group_by(sex) %>% 
             mutate(Ins_tAUC = log(Ins_tAUC),
                    Ins_tAUC = scale(Ins_tAUC))
```

``` error
Error: object 'pheno_clin' not found
```

``` r
ins_tauc %>% 
  ggplot(aes(x = sex, y = Ins_tAUC)) +
    geom_boxplot() +
    geom_hline(aes(yintercept = -4), color = 'red') +
    geom_hline(aes(yintercept =  4), color = 'red') +
    labs(title = "Distribution of Standardized Ins_tAUC")
```

``` error
Error: object 'ins_tauc' not found
```

There are no data points outside of the four standard deviation limits.

## Gene Expression Phenotypes


``` r
# load the expression data along with annotations and metadata
load("../data/dataset.islet.rnaseq.RData")
```

``` warning
Warning in readChar(con, 5L, useBytes = TRUE): cannot open compressed file
'../data/dataset.islet.rnaseq.RData', probable reason 'No such file or
directory'
```

``` error
Error in readChar(con, 5L, useBytes = TRUE): cannot open the connection
```

``` r
names(dataset.islet.rnaseq)
```

``` error
Error: object 'dataset.islet.rnaseq' not found
```


``` r
# look at gene annotations
dataset.islet.rnaseq$annots[1:6,]
```

``` error
Error: object 'dataset.islet.rnaseq' not found
```

``` r
# look at raw counts
dataset.islet.rnaseq$raw[1:6,1:6]
```

``` error
Error: object 'dataset.islet.rnaseq' not found
```

``` r
# look at sample metadata
# summarize mouse sex, birth dates and DO waves
table(dataset.islet.rnaseq$samples[, c("sex", "birthdate")])
```

``` error
Error: object 'dataset.islet.rnaseq' not found
```

``` r
table(dataset.islet.rnaseq$samples[, c("sex", "DOwave")])
```

``` error
Error: object 'dataset.islet.rnaseq' not found
```

In order to make reasonable gene comparisons between samples, the count data 
need to be normalized. In the quantile-quantile (Q-Q) plot below, count data for 
the first gene are plotted over a diagonal line tracing a normal distribution 
for those counts. Notice that most of the count data values lie off of this 
line, indicating that these gene counts are not normally distributed. 


``` error
Error: object 'dataset.islet.rnaseq' not found
```

Q-Q plots for the first six genes show that count data for these genes are not
normally distributed. They are also not on the same scale. The y-axis values for
each subplot range to 20,000 counts in the first subplot, 250 in the second, 90
in the third, and so on. 


``` r
dataset.islet.rnaseq$raw %>% 
  as.data.frame() %>%
  select(ENSMUSG00000000001:ENSMUSG00000000058) %>% 
  pivot_longer(cols = everything(), names_to = 'gene', values_to = 'value') %>% 
  ggplot(aes(sample = value)) +
    stat_qq() +
    geom_qq_line() +
    facet_wrap(~gene, scales = 'free') +
    labs(title = 'Count distribution for six genes',
         xlab = 'Normal percentiles', y = 'Count percentiles')
```

``` error
Error: object 'dataset.islet.rnaseq' not found
```

Q-Q plots of the normalized expression data for the first six genes show that 
the data values match the diagonal line well, meaning that they are now normally
distributed. They are also all on the same scale now as well.


``` r
dataset.islet.rnaseq$expr %>% 
  as.data.frame() %>%
  select(ENSMUSG00000000001:ENSMUSG00000000058) %>% 
  pivot_longer(cols = everything(), names_to = 'gene', values_to = 'value') %>% 
  ggplot(aes(sample = value)) +
    stat_qq() +
    geom_qq_line() +
    facet_wrap(~gene, scales = 'free') +
    labs(title = 'Normalized count distribution for six genes',
         xlab = 'Normal percentiles', y = 'Count percentiles')
```

``` error
Error: object 'dataset.islet.rnaseq' not found
```

Boxplots of raw counts for six example genes are shown at left below. Notice 
that the median count values (horizontal black bar in each boxplot) are not 
comparable between the genes because the counts are not on the same scale. At
right, boxplots for the same genes show normalized count data on the same 
scale.


``` r
raw = dataset.islet.rnaseq$raw %>% 
        as.data.frame() %>% 
        select(ENSMUSG00000000001:ENSMUSG00000000058) %>% 
        pivot_longer(cols = everything(), names_to = 'gene', values_to = 'value') %>% 
        mutate(type = 'raw')
```

``` error
Error: object 'dataset.islet.rnaseq' not found
```

``` r
norm = dataset.islet.rnaseq$expr %>% 
         as.data.frame() %>% 
         select(ENSMUSG00000000001:ENSMUSG00000000058) %>% 
         pivot_longer(cols = everything(), names_to = 'gene', values_to = 'value') %>% 
         mutate(type = 'normalized')
```

``` error
Error: object 'dataset.islet.rnaseq' not found
```

``` r
bind_rows(raw, norm) %>%
  mutate(type = factor(type, levels = c('raw', 'normalized'))) %>% 
  ggplot(aes(gene, value)) +
    geom_boxplot() +
    facet_wrap(~type, scales = 'free') +
    labs(title = 'Count distributions for example genes') +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1))
```

``` error
Error in `bind_rows()` at dplyr/R/mutate.R:146:3:
! Argument 1 must be a data frame or a named atomic vector.
```

``` r
rm(raw, norm)
```

``` warning
Warning in rm(raw, norm): object 'raw' not found
```

``` warning
Warning in rm(raw, norm): object 'norm' not found
```

Have a look at the first several rows of normalized count data.


``` r
# look at normalized counts
dataset.islet.rnaseq$expr[1:6,1:6]
```

``` error
Error: object 'dataset.islet.rnaseq' not found
```

The expression data loaded provides LOD peaks for the eQTL analyses performed in
this study. As a preview of what you will be doing next, look at the first 
several rows of LOD peak values and extract the LOD peaks for chromosome 11.


``` r
# look at LOD peaks
dataset.islet.rnaseq$lod.peaks[1:6,]
```

``` error
Error: object 'dataset.islet.rnaseq' not found
```

``` r
# look at chromosome 11 LOD peaks
chr11_peaks <- dataset.islet.rnaseq$annots %>% 
   select(gene_id, chr) %>% 
   filter(chr=="11") %>%
   left_join(dataset.islet.rnaseq$lod.peaks, 
             by = c("chr" = "chrom", "gene_id" = "annot.id")) 
```

``` error
Error: object 'dataset.islet.rnaseq' not found
```

``` r
# look at the first several rows of chromosome 11 peaks
head(chr11_peaks)
```

``` error
Error: object 'chr11_peaks' not found
```

``` r
# how many rows?
dim(chr11_peaks)
```

``` error
Error: object 'chr11_peaks' not found
```

``` r
# how many rows have LOD scores?
chr11_peaks %>% filter(!is.na(lod)) %>% dim()
```

``` error
Error: object 'chr11_peaks' not found
```

``` r
# sort chromosome 11 peaks by LOD score
chr11_peaks %>% arrange(desc(lod)) %>% head()
```

``` error
Error: object 'chr11_peaks' not found
```

``` r
# range of LOD scores and positions
range(chr11_peaks$lod, na.rm = TRUE)
```

``` error
Error: object 'chr11_peaks' not found
```

``` r
range(chr11_peaks$pos, na.rm = TRUE)
```

``` error
Error: object 'chr11_peaks' not found
```

``` r
# view LOD scores by position
chr11_peaks %>% arrange(desc(lod)) %>% 
  ggplot(aes(pos, lod)) + geom_point()
```

``` error
Error: object 'chr11_peaks' not found
```

::::::::::::::::::::::::::::::::::::: keypoints 

- Use `.md` files for episodes when you want static content
- Use `.Rmd` files for episodes when you need to generate output
- Run `sandpaper::check_lesson()` to identify any issues with your lesson
- Run `sandpaper::build_lesson()` to preview your lesson locally

::::::::::::::::::::::::::::::::::::::::::::::::
