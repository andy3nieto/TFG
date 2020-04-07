################################################################
#### METHYLATION ANALYSIS PIPELINE FUNCTIONAL NORMALISATION ####         
################################################################

library(knitr)
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(Gviz)
library(DMRcate) #dmr's 
library(stringr)
library(RColorBrewer)
library(plyr)
library(dplyr)

#########################################
########## Reading in the data ##########
#########################################

setwd("/media/andrea/SAMSUNG/IJC/Aim2/methylation")

# folder arrays_0102 specifically created for clustering and identification of B-ALL samples. This means that a special sample sheet is created for the cause
# also, another column is added to know which samples come from which array (Unkown/T)
both_all <- read.metharray.sheet("arrays_0102/", pattern = "sample_sheet_forclustering.csv")

rgSet <- read.metharray.exp(targets = both_all)

#####################################
########## Quality Control ##########
#####################################

# quality of the signal of every CpG is defined by the p-values
detP <- detectionP(rgSet)

keep <- colMeans(detP) < 0.05
rgSet <- rgSet[,keep]

both_all <- both_all[keep,]

detP <- detP[,keep]

##############################################
########## Functional Normalisation ##########
##############################################

mSetSq <-preprocessFunnorm(rgSet)

mSetRaw <- preprocessRaw(rgSet)


################################
########## Clustering ##########
################################

# get methlyation data matrix 
# rows are CpGs and columns are samples
M <- getM(mSetSq) 

library(RColorBrewer)
pal <- brewer.pal(8,"Dark2")

par(mfrow=c(1,1))
pdf("/media/andrea/SAMSUNG/IJC/Aim2/methylation/FUNCTIONALNORMALISATION/plotMDS_separationT.B-ALL.pdf", width = 12, height = 12)
a <- plotMDS(M, top=1000, gene.selection="common", 
             col=pal[factor(both_all$Cell_Type)])
dev.off()


select_data <- names(a$y[a$y > 0])
length(select_data) #39 samples of B-ALL to use in analysis

#save(select_data, file = "/media/andrea/SAMSUNG/IJC/Aim2/methylation/b_all_samples_funnorm.RData")

############################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################

#########################################
########## Reading in the data ##########
#########################################

# load selected b-all samples
x <- load(file = "/media/andrea/SAMSUNG/IJC/Aim2/methylation/b_all_samples_funnorm.RData")

# (1) reading b-all cell samples
B_all <- read.metharray.sheet("arrays_01", pattern = "sample_sheet_modified.csv")

B_all$ID <- paste("Unknown",row.names(B_all),sep=".")
keep <- B_all$ID %in% select_data
B_all <- B_all[keep,]

B_all$ID <- paste("B.all",row.names(B_all),sep=".")

# (2) reading precursor cell samples
precursor_cells <- read.metharray.sheet("idat_files", pattern = "sample_sheet.csv")
keep <- grepl('^S1|S2', precursor_cells$Sample_Name)
precursor_cells <- precursor_cells[keep,]

precursor_cells$ID <- precursor_cells$Sample_Name

#joining tables
targets <- rbind(precursor_cells,B_all[,c(1,2,8,9,10,11,12)])


#adding cell type column
s1 <- grepl('^S1', targets$Sample_Name)
s2 <- grepl('^S2', targets$Sample_Name)
#pbmc <- grepl('^PBMC', targets$Sample_Name)

targets$cell_type <- rep("B_all", nrow(targets))
targets$cell_type[s1] <- "S1"
targets$cell_type[s2] <- "S2"
#targets$cell_type[pbmc] <- "PBMC"


#creating rgset of all samples selected
rgset <- read.metharray.exp(targets = targets)  # 622399 probes, 51 samples
sampleNames(rgset) <- targets$ID


#####################################
########## Quality Control ##########
#####################################


# examine mean detection p-values across all samples to identify any failed samples
detP <- detectionP(rgset)
head(detP)

# plots for proper visualisation 
# (pvalues_barplots_qc.pdf)
pal <- brewer.pal(8,"Dark2")
par(mfrow=c(1,2))
barplot(colMeans(detP), col=pal[factor(targets$cell_type)], las=2, 
        cex.names=0.8, ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$cell_type)), fill=pal,
       bg="white")

barplot(colMeans(detP), col=pal[factor(targets$cell_type)], las=2, 
        cex.names=0.8, ylim=c(0,0.002), ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$cell_type)), fill=pal, 
       bg="white")

# quality control report (minfi package) generates plots for quality control
# read minfi vignette for a proper understanding of the plots
qcReport(rgset, sampNames=targets$Sample_Name, sampGroups=targets$cell_type, 
         pdf="/media/andrea/SAMSUNG/IJC/Aim2/methylation/FUNCTIONALNORMALISATION/qcReport.pdf")


# remove poor quality samples --> nothing to remove 
keep <- colMeans(detP) < 0.05
table(keep) # NO FALSE


###################################
########## Normalisation ##########
###################################

# normalising data, genomicratioset it output 
msetsq <- preprocessFunnorm(rgset) 

# create a MethylSet object from the raw data for plotting
msetraw <- preprocessRaw(rgset)

# visualise what the data looks like before and after normalisation
par(mfrow=c(1,2))

densityPlot(rgset, sampGroups=targets$cell_type,main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targets$cell_type)), 
       text.col=brewer.pal(8,"Dark2"))

densityPlot(getBeta(msetsq), sampGroups=targets$cell_type,
            main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(targets$cell_type)), 
       text.col=brewer.pal(8,"Dark2"))


## DATA EXPLORATION ##
# multi-dimensional scaling plots used to visualise data: (as before)

par(mfrow=c(1,1))
plotMDS(getM(msetsq), top=1000, gene.selection="common", 
        col=pal[factor(targets$cell_type)])
legend("top", legend=levels(factor(targets$cell_type)), text.col=pal,
       bg="white", cex=0.7)


###############################
########## Filtering ##########
###############################

# ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(msetsq),rownames(detP)),] 
head(detP)


# remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.01) == ncol(msetsq) 
table(keep)

msetsqflt <- msetsq[keep,]
# FILTERED FROM 485512 -> 472368 #



# if your data includes males and females, remove probes on the sex chromosomes
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

keep <- !(featureNames(msetsqflt) %in% ann450k$Name[ann450k$chr %in% 
                                                      c("chrX","chrY")])
table(keep)
msetsqflt <- msetsqflt[keep,]
# FILTERED FROM 472368 -> 461260 #



# remove probes with SNPs at CpG site
msetsqflt <- dropLociWithSnps(msetsqflt)
# FILTERED FROM  461260 -> 445372 



# re-inspection of data
par(mfrow=c(1,1))
plotMDS(getM(msetsqflt), top=1000, gene.selection="common", 
        col=pal[factor(targets$cell_type)], cex=0.8)
legend("top", legend=levels(factor(targets$cell_type)), text.col=pal,
       bg="white", cex=0.7)
# PC1 creates quite a separation between b-all and precursors clusters,


#### CALCULATION OF VALUES ####

# calculate M-values for statistical analysis
mVals <- getM(msetsqflt)
head(mVals[,1:5])

#                S2_11     S1_12     S2_15     S1_16      S2_19
# cg13869341  2.880938  2.831305  2.860820  2.646219  2.7490868
# cg14008030  1.488088  1.591628  1.427523  1.204298  1.5126614
# cg12045430 -3.259536 -2.749680 -2.867267 -3.447351 -4.0054758
# cg20826792 -1.468585 -1.721705 -1.898119 -2.362716 -1.3189292
# cg00381604 -4.389305 -5.042657 -3.811603 -4.746788 -4.9247094
# cg20253340  1.121894  1.317122  2.377750  3.052654 -0.1858459

bVals <- getBeta(msetsqflt)
head(bVals[,1:5])

#                 S2_11      S1_12      S2_15      S1_16      S2_19
# cg13869341 0.88047267 0.87680447 0.87899737 0.86226362 0.87051520
# cg14008030 0.73719957 0.75086531 0.72898592 0.69735962 0.74048617
# cg12045430 0.09454702 0.12943848 0.12052817 0.08397522 0.05861375
# cg20826792 0.26542784 0.23265227 0.21153859 0.16277688 0.28613870
# cg00381604 0.04554522 0.02944616 0.06648369 0.03590814 0.03187471
# cg20253340 0.68517076 0.71360535 0.83863821 0.89244266 0.46783982

par(mfrow=c(1,2))
densityPlot(bVals, sampGroups=targets$cell_type, main="Beta values", 
            legend=FALSE, xlab="Beta values")
legend("top", legend = levels(factor(targets$cell_type)), 
       text.col=brewer.pal(8,"Dark2"))
densityPlot(mVals, sampGroups=targets$cell_type, main="M-values", 
            legend=FALSE, xlab="M values")
legend("topleft", legend = levels(factor(targets$cell_type)), 
       text.col=brewer.pal(8,"Dark2"))

write.table(mVals ,"/media/andrea/SAMSUNG/IJC/Aim2/methylation/FUNCTIONALNORMALISATION/mvals.txt", row.names = T, col.names = T, sep = "\t", quote = F)
write.table(bVals, "/media/andrea/SAMSUNG/IJC/Aim2/methylation/FUNCTIONALNORMALISATION/betavals.txt", row.names = T, col.names = T, sep = "\t", quote = F)

###########################################################################
########## Probe-wise differential methylation analysis // DMP's ##########
###########################################################################

# this is the factor of interest
cellType <- factor(targets$cell_type)

design <- model.matrix(~0+cellType, data=targets) #defining in matrix which sample is from which cell type
colnames(design) <- c(levels(cellType))
rownames(design) <- c(targets$ID)

#fit linear model for methylation values obtained above
fit <- lmFit(mVals, design)

#create contrast matrix for specific comparisons between cell types
contMatrix <- makeContrasts("B_all-(S1+S2)/2", levels = design)

#fitting the contrasts (comparisons between cell types)
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit2))

# A total of 116122+109663 = 225785 CpGs are found as significantly differentially methylated 
# 
# 
#        B_all-(S1+S2)/2
# Down            116122
# NotSig          219587
# Up              109663

#obtaining DMPs
ann450kSub <- ann450k[match(rownames(mVals),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]

cpgs <- topTable(fit2, num=Inf, coef=1, genelist=ann450kSub)
dmps <- cpgs[cpgs$adj.P.Val <= 0.05,] #225785x29x

dmps$cpg_name <- rownames(dmps)
rownames(dmps) <- c(1:length(rownames(dmps)))
dmps[dmps == ""] <- NA

write.table(dmps, file="/media/andrea/SAMSUNG/IJC/Aim2/methylation/FUNCTIONALNORMALISATION/DMPs.txt", sep="\t", quote = F, col.names = T, row.names = F)
# there are + or - 14000 dmps less than with pre-process quantile normalisation approach

dmps <- dplyr::as_tibble(dmps)

DMstatus <- decideTests(fit2)
DMstatus <- as.data.frame(DMstatus)

DMstatus$cpg_name <- rownames(DMstatus)
rownames(DMstatus) <- c(1:length(rownames(DMstatus)))
colnames(DMstatus)[1] <- "DMstat"
DMstatus <- dplyr::as_tibble(DMstatus)

joined <- inner_join(dmps, DMstatus)

write.table(joined, file="/media/andrea/SAMSUNG/IJC/Aim2/methylation/FUNCTIONALNORMALISATION/DMPs_with_DMstatus.txt", sep="\t", quote = F, col.names = T, row.names = F)


table(joined$Enhancer) #44675
























































