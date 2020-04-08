## STEP 1.3
### Integration  

We have <span style="color: red;">225,785 DMPs</span> from out differential methylation analysis. These must be integrated with our enhancer regions in order to obtain a list of DMPs that fall within enhancers.  


```
## GRanges object with 5837 ranges and 2 metadata columns:
##          seqnames    ranges strand |        name dm_status
##             <Rle> <IRanges>  <Rle> | <character> <integer>
##      [1]        1   1856073      * |  cg19356311        -1
##      [2]        1   1856512      * |  cg08440544         1
##      [3]        1   2023482      * |  cg19906777         1
##      [4]        1   2032058      * |  cg07912723        -1
##      [5]        1   2098929      * |  cg09144608        -1
##      ...      ...       ...    ... .         ...       ...
##   [5833]        9 139040286      * |  cg14384803        -1
##   [5834]        9 139040738      * |  cg13830329         1
##   [5835]        9 139081755      * |  cg21229482        -1
##   [5836]        9 139383171      * |  cg13428009        -1
##   [5837]        9 140401455      * |  cg25812043         1
##   -------
##   seqinfo: 22 sequences from an unspecified genome; no seqlengths
```

Above we see an excerpt of the DMPs falling in enhancer regions. This is a total of 5837.  
The number of enhancer regions involved is 5,234. This means that some enhancers have more than one DMP.  


  
  
## STEP 2.1  
### Data  

This next step involves the integration of enhancer DMPs with PCHi-C data (be that of Naive B, HSC, pre-B or pro-B). For now the only PCHi-C we have available is that of Naive B, which we have to load into the session using the novel JavierreLabR package :)  


  
  
A total of 192,104 interactions are annotated from the Naive-B cell PCHi-C data.  
However, for now, we will exclude from the analysis any interaction involving fragments enominated as UCEs or non-annotated. This leaves us with a <span style="color: red;">total of 189284 interactions to work with</span>

  
## STEP 2.2 
### Integration  

The PCHi-C interaction data and the enhancer DMPs are overlapped at this point, also using JavierreLabR package.

<span style="color: red;">Integration of Naive-B interactions with enhancer DMPs results in 4284 interactions.</span>  
  
  
  
## STEP 3  
### Exploration  

**What types of interactions to we have in this subset?**

```
## 
##  OE_P_I  P_OE_I P_OE_II   P_P_I  P_P_II 
##     948    2697      22     608       9
```
  
Knowing that the DMPs fall in enhancers, we can assume that those interactions with peaks in P (captured regions) dont actually fall within promoter regions.  

The only types of interactions that will not be taken into account are the P_OE_I.  

We have two subgroup of interactions:  
1. OE_P_I + P_OE_II (970 interactions)  
2. P_P_I + P_P_II (617 interactions)  
  
  
  
## STEP 4
### Listing possible transcripts/genes being dysregulated

> SUGROUP 1  


Number of transcripts possibly being dysregulated by DMPs at long-range enhancers: 4632.  
This translates to <span style="color: red;">1107 genes</span>  
  
  
> SUGROUP 2  


Number of transcripts possibly being dysregulated by DMPs at long-range enhancers: 2015    
This translates to <span style="color: red;">500 genes</span>  
  
  
Subgroup   | Interactions    | N^o^ transcripts | N^o^ genes
-----------|-----------------|------------------|-------------
1          | OE_P_I, P_OE_II | 4632             | 1107
2          | P_P_I, P_P_II   | 2015             | 500
both       | all above       | 6647             | 1607  
  
  
