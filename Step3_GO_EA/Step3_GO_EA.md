---
title: "GO Enrichment Analysis"
author: "Andrea Nieto-Aliseda Sutton"
date: "07/04/2020"
output:
  html_document:
    keep_md: true
---




### __Summary__

A GO Enrichment Analysis sets out to identify GO terms that are more related to a subset of genes than expected. The expectation of how frequent a GO term should appear is defined by the __universe__, being a sort of a background model.  

The list of genes we are going to use as query has been obtained by the integration between B-ALL associated DMPs falling at enhancer regions as defined by the ensembl regulatory build and PCHi-C data (for now we are using Naive B).  
The genes therefore are those which are possibly being dysregulated by long-range regulatory elements (enhancers).  
<br>  
 
### __Clusterprofiler__

To start, the subset of genes has to be read in. 
<br>  

The subset of genes we want to study consist of __``1607`` genes__
<br>  

__``29139``__ genes will be used as background  
<br>  

The only two GO terms which pass the FDR cutoff < __0.05__ for BP ontology are the following:

```
## # A tibble: 2 x 6
##   ID      Description                        GeneRatio   pvalue p.adjust  qvalue
##   <chr>   <chr>                              <chr>        <dbl>    <dbl>   <dbl>
## 1 GO:000… antigen processing and presentati… 7/920      1.13e-6  0.00284 0.00280
## 2 GO:001… antigen processing and presentati… 7/920      1.13e-6  0.00284 0.00280
```
<br>
With a p-value cutoff of __0.1__ the following 3 enriched terms are found:

```
## # A tibble: 3 x 6
##   ID      Description                        GeneRatio   pvalue p.adjust  qvalue
##   <chr>   <chr>                              <chr>        <dbl>    <dbl>   <dbl>
## 1 GO:000… antigen processing and presentati… 7/920      1.13e-6  0.00284 0.00280
## 2 GO:001… antigen processing and presentati… 7/920      1.13e-6  0.00284 0.00280
## 3 GO:001… antigen processing and presentati… 7/920      5.04e-5  0.0847  0.0836
```
<br>
With a p-value and q-value cutoff of __0.25__ the following enriched terms are found; 0.25 has been decided as a cutoff given the explanation given here: https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/FAQ#Why_does_GSEA_use_a_false_discovery_rate_.28FDR.29_of_0.25_rather_than_the_more_classic_0.05.3F

```
## # A tibble: 6 x 6
##   ID      Description                        GeneRatio   pvalue p.adjust  qvalue
##   <chr>   <chr>                              <chr>        <dbl>    <dbl>   <dbl>
## 1 GO:000… antigen processing and presentati… 7/920      1.13e-6  0.00284 0.00280
## 2 GO:001… antigen processing and presentati… 7/920      1.13e-6  0.00284 0.00280
## 3 GO:001… antigen processing and presentati… 7/920      5.04e-5  0.0847  0.0836 
## 4 GO:002… central nervous system neuron axo… 9/920      1.07e-4  0.135   0.133  
## 5 GO:000… antigen processing and presentati… 14/920     1.89e-4  0.191   0.188  
## 6 GO:004… antigen processing and presentati… 12/920     2.93e-4  0.246   0.243
```
<br>  
In essence, the following table shows the GO terms enriched for different cutoff values:  

ID                  | Description                 | p-adjust    | qvalue
--------------------|-----------------------------|-----------------|-----------
GO:0002483      | antigen processing and presentation of endogenous peptide antigen | 0.002837516 | 0.002798439
GO:0019885      | antigen processing and presentation of endogenous peptide antigen via MHC class I | 0.002837516 | 0.002798439
--------------------|-----------------------------|-----------------|-----------
GO:0019883      |	antigen processing and presentation of endogenous antigen                       | 0.084738535 | 0.083571577
--------------------|-----------------------------|-----------------|-----------
GO:0021955      | central nervous system neuron axonogenesis                                      | 0.134774977 | 0.132918952	
GO:0002474      | antigen processing and presentation of peptide antigen via MHC class I          | 0.190618215 | 0.187993157	
GO:0042590	    | antigen processing and presentation of exogenous peptide antigen via MHC class I | 0.246106157 | 0.242716959	

link to MHC wiki: https://en.wikipedia.org/wiki/Major_histocompatibility_complex  
   
<br>
Lets try with a different ontology like Molecular Function (MF) (pvalue < 0.5, qvalue < 0.5)

```
## # A tibble: 10 x 6
##    ID      Description                         GeneRatio  pvalue p.adjust qvalue
##    <chr>   <chr>                               <chr>       <dbl>    <dbl>  <dbl>
##  1 GO:001… sodium channel regulator activity   8/897     4.40e-4    0.308  0.305
##  2 GO:003… heat shock protein binding          16/897    7.16e-4    0.308  0.305
##  3 GO:001… hydrolase activity, acting on carb… 4/897     1.72e-3    0.472  0.468
##  4 GO:000… ribonuclease activity               14/897    2.68e-3    0.472  0.468
##  5 GO:000… RNA polymerase I activity           4/897     3.69e-3    0.472  0.468
##  6 GO:001… hydrolase activity, acting on carb… 14/897    3.79e-3    0.472  0.468
##  7 GO:004… peptide antigen binding             5/897     4.64e-3    0.472  0.468
##  8 GO:004… protein binding involved in protei… 5/897     4.64e-3    0.472  0.468
##  9 GO:014… catalytic activity, acting on RNA   33/897    4.99e-3    0.472  0.468
## 10 GO:000… endoribonuclease activity           9/897     5.49e-3    0.472  0.468
```
These terms have very low statistical significance  

<br>
Lets try with a different ontology like Cellular Component (CC) (pvalue < 0.5, qvalue < 0.5)

```
## # A tibble: 4 x 6
##   ID       Description                         GeneRatio  pvalue p.adjust qvalue
##   <chr>    <chr>                               <chr>       <dbl>    <dbl>  <dbl>
## 1 GO:0042… MHC protein complex                 6/980     5.14e-4    0.326  0.320
## 2 GO:0030… integral component of endoplasmic … 18/980    1.16e-3    0.347  0.341
## 3 GO:0045… keratin filament                    12/980    2.03e-3    0.347  0.341
## 4 GO:0031… intrinsic component of endoplasmic… 18/980    2.19e-3    0.347  0.341
```

<br>  

**What about using transcripts instead of genes to gain statistical significance?**  





```r
library(GOplot)
```

```
## Loading required package: ggdendro
```

```
## Loading required package: gridExtra
```

```
## 
## Attaching package: 'gridExtra'
```

```
## The following object is masked from 'package:dplyr':
## 
##     combine
```

```
## The following object is masked from 'package:Biobase':
## 
##     combine
```

```
## The following object is masked from 'package:BiocGenerics':
## 
##     combine
```

```
## Loading required package: RColorBrewer
```

```r
r <- results10[,c(1,2,8,6,7,9)]
colnames(r) <- c("ID","Term","Genes","adj_pval","zscore","count")

GOBar(data = r[1:150,])
```

![](Step3_GO_EA_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

```r
goplot(bpt, showCategory = 2)
```

![](Step3_GO_EA_files/figure-html/unnamed-chunk-10-2.png)<!-- -->
<br>
<br>


__What about the Network of Cancer Gene?__  
Network of Cancer Gene (NCG)(A. et al. 2016) is a manually curated repository of cancer genes. NCG release 5.0 (Aug. 2015) collects 1,571 cancer genes from 175 published studies

```r
bitr <- bitr(g, "ENSEMBL", "ENTREZID", OrgDb = org.Hs.eg.db)
```

```
## 'select()' returned 1:many mapping between keys and columns
```

```
## Warning in bitr(g, "ENSEMBL", "ENTREZID", OrgDb = org.Hs.eg.db): 26.92% of input
## gene IDs are fail to map...
```

```r
ncg <- enrichNCG(bitr$ENTREZID, pvalueCutoff = 0.3, minGSSize = 100)
ncg@result[,c(4:7)]
```

```
##                BgRatio      pvalue   p.adjust     qvalue
## head_and_neck 102/1571 0.003416837 0.02391786 0.02158003
## melanoma      134/1571 0.061839006 0.21643652 0.19528107
## leukemia      212/1571 0.155671395 0.36323325 0.32772925
## lymphoma      162/1571 0.413442643 0.54181370 0.48885447
## glioblastoma  104/1571 0.414072107 0.54181370 0.48885447
## pancreas      149/1571 0.464411745 0.54181370 0.48885447
## lung          143/1571 0.576795415 0.57679542 0.52041692
```
Here, count of genes related to LEUKEMIA is the highest (14/78)  

<br>  
<br>  
<br>  
<br> 
















