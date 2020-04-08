---
title: "Enhancer - CpG statistics"
author: "Andrea Nieto-Aliseda Sutton"
date: "4/6/2020"
output: 
  html_document:
    keep_md: true
---



## STEP 1.1 
### Data  

We have a regulatory build from ensembl (homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20180925.gff). From this, we have to retreive all the enhancer regions. To do this we can use bash command line:



The *enhancer_regions.txt* file will look like this:

```
##                ID chr  start    end
## 1 ENSR00000424863   1  22801  23000
## 2 ENSR00000424874   1  89201  89800
## 3 ENSR00000424875   1  90401  91800
## 4 ENSR00000424880   1 132001 133000
## 5 ENSR00000424886   1 143601 144400
## 6 ENSR00000424895   1 246801 247800
```
In this table we have <span style="color: red;">140,648 enhancer regions</span>
  
------------------ 
  
## STEP 1.2  
### Exploratory Analysis of Regulatory Build (specifically of enhancer regions)  

The following list of statistics are of our interest to put our integration analysis in perspective:  

* Number of Enhancers  (140,648)
* Total number of base pairs composed by enhancers
* Percentage of genome covered by enhancers
* Number of CpGs (of genome) located at these enhancers
* Number of CpGs in enhancers overlapping a captured HindIII fragments (These CpGs will be analysed on the Captured-Captured interactions!)
* Number of CpGs located at enhancers covered by 450K array
* Number of DMPs located at enhancers
* Number of target genes potentially altered by these DMPs -- this is obtained from integration 

> Total number of base pairs composed by enhancers


<span style="color: red;">``result: 107309123bp``</span>  
  
  
  
> Percentage of genome covered by enhancers
  
The reference for the number of base pairs composing the genome as of version GRCh37.p13 is https://grch37.ensembl.org/Homo_sapiens/Info/Annotation  
The length of the genome is 3,326,743,047  


<span style="color: red;">``result: 3.23%``</span>  
  
  
  
> Number of CpGs (of genome) located at these enhancers  
  
To calculate this the CpG positions in the genome must be obtained and overlapped with the enhancer regions we have from ENSEMBL's regulatory build.


<span style="color: red;">`28217009 CpG sites identified in the genome`</span>

This number of corroborated by the following link: https://emea.illumina.com/content/dam/illumina-marketing/documents/products/datasheets/truseq-methyl-capture-epic-sequencing-panel-data-sheet-470-2016-004.pdf"




<span style="color: red;">``932247 genome CpGs fall at enhancer regions``</span>  
<span style="color: red;">``This is 3.3% of the CpGs in the genome``</span>  
  
  
> Number of CpGs in enhancers overlapping a captured HindIII fragments (These CpGs will be analysed on the Captured-Captured interactions!)
  

<span style="color: red;">``56672 CpGs related to enhancers fall in captured fragments``</span>  
<span style="color: red;">``This is 6.08% of the CpGs falling at enhancer regions``</span>  
<span style="color: red;">``This is 0.2% of the CpGs in the genome``</span>  
  
  
> Number of CpGs located at enhancers covered by 450K array 
  
For this, we need the coordinates of all the CpGs in the 450K array, which I have already saved in a file. This information has originally been obtained from the manifesto (annotation R package)

```
## # A tibble: 6 x 3
##   Name       chr     pos
##   <chr>      <chr> <int>
## 1 cg13869341 1     15865
## 2 cg14008030 1     18827
## 3 cg12045430 1     29407
## 4 cg20826792 1     29425
## 5 cg00381604 1     29435
## 6 cg20253340 1     68849
```
The CpG positions can then be overlapped with the enhancer coordinates, obtaining the number of CpGs at enhancers regions in the process

<span style="color: red;"> ``15973 CpGs from 450K array fall in enhancer regions``</span>  

  
> Number of DMPs located at enhancers   

This number will be obtained in the next step (step1.3)
  
  
  
  
### SUMMARY TABLE

Characteristic                   | Value
---------------------------------|---------------------------
Length of genome (bp)            | `3326743047`
mean width of enhancers (bp)     | ``763``
N^o^ bp covered by enhancers     | ``107309123``
% of genome covered by enhancers | ``3.23``
---------------------------------|---------------------------
N^o^ of CpG sites in genome      | ``28217009``
N^o^ of CpG sites at enhancers   | ``932247``
% of CpG sites at enhancers      | ``3.3%``
---------------------------------|---------------------------
N^o^ of CpGs covered by 450K array    | ``485512``
N^o^ of enhancer CpGs covered by 450K | ``15973``
% genome CpGs covered by 450K array   | ``0.01%``
% 450K CpGs falling in enhancers      | ``3.29%``
% enhancer CpGs covered by 450K       | ``1.71%``
---------------------------------|---------------------------
N^o^ of enhancer CpGs at captured HindIII fragments      | ``56672``
% of covered enhancer CpGs at captured HindIII fragments | ``11.67%``
---------------------------------|---------------------------
N^o^ of DMPs                        | `225785`
N^o^ of DMPs at enhancer regions    | `5837`
% DMPs of covered enhancer CpGs     | ``1.2%``
% DMPs of genome enhancer CpGs      | ``0.63%``
