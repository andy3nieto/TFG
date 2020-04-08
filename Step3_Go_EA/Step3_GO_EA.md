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
With a p-value cutoff of __0.1__ the following enriched terms are found:

```
## # A tibble: 3 x 6
##   ID      Description                        GeneRatio   pvalue p.adjust  qvalue
##   <chr>   <chr>                              <chr>        <dbl>    <dbl>   <dbl>
## 1 GO:000… antigen processing and presentati… 7/920      1.13e-6  0.00284 0.00280
## 2 GO:001… antigen processing and presentati… 7/920      1.13e-6  0.00284 0.00280
## 3 GO:001… antigen processing and presentati… 7/920      5.04e-5  0.0847  0.0836
```
<br>
With a p-value and q-value cutoff of __0.5__ the following enriched terms are found:

```
## # A tibble: 13 x 6
##    ID      Description                       GeneRatio   pvalue p.adjust  qvalue
##    <chr>   <chr>                             <chr>        <dbl>    <dbl>   <dbl>
##  1 GO:000… antigen processing and presentat… 7/920      1.13e-6  0.00284 0.00280
##  2 GO:001… antigen processing and presentat… 7/920      1.13e-6  0.00284 0.00280
##  3 GO:001… antigen processing and presentat… 7/920      5.04e-5  0.0847  0.0836 
##  4 GO:002… central nervous system neuron ax… 9/920      1.07e-4  0.135   0.133  
##  5 GO:000… antigen processing and presentat… 14/920     1.89e-4  0.191   0.188  
##  6 GO:004… antigen processing and presentat… 12/920     2.93e-4  0.246   0.243  
##  7 GO:002… central nervous system neuron de… 13/920     3.56e-4  0.256   0.253  
##  8 GO:000… transcription elongation from RN… 7/920      9.50e-4  0.495   0.489  
##  9 GO:003… Arp2/3 complex-mediated actin nu… 7/920      9.50e-4  0.495   0.489  
## 10 GO:000… cell cycle arrest                 25/920     9.90e-4  0.495   0.489  
## 11 GO:001… regulation of epithelial cell mi… 28/920     1.12e-3  0.495   0.489  
## 12 GO:000… termination of RNA polymerase I … 7/920      1.18e-3  0.495   0.489  
## 13 GO:000… response to hypoxia               33/920     1.28e-3  0.496   0.489
```
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
<br>  
<br>  
<br>  
<br> 
__What about a Disease Ontology?__
For this it is necessary to map our ENSEMBL gene IDs to ENTREZ gene IDs, which fails to map 26.92% in the process, leaving us with 1146 genes

#### NOTHING
<br>  
<br>  
<br>  
<br> 



























