############################################################################################
############################## METHYLATION ANALYSIS PIPELINE  ##############################         
############################################################################################

library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(stringr)
library(RColorBrewer)
library(plyr)
library(dplyr)


#########################################
########## Reading in the data ##########
#########################################

# load b-all samples
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


######################################
########## Data Exploration ##########
######################################

# multi-dimensional scaling plots used to visualise data:

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
# PC1 creates quite a separation between b-all and precursors clusters


#### CALCULATION OF VALUES ####

# calculate M-vals and B-vals for statistical analysis
mVals <- getM(msetsqflt)
bVals <- getBeta(msetsqflt)

# write.table(mVals ,"/media/andrea/SAMSUNG/IJC/Aim2/methylation/FUNCTIONALNORMALISATION/mvals.txt", row.names = T, col.names = T, sep = "\t", quote = F)
# write.table(bVals, "/media/andrea/SAMSUNG/IJC/Aim2/methylation/FUNCTIONALNORMALISATION/betavals.txt", row.names = T, col.names = T, sep = "\t", quote = F)


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


# obtaining DMPs
ann450kSub <- ann450k[match(rownames(mVals),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]

cpgs <- topTable(fit2, num=Inf, coef=1, genelist=ann450kSub)
dmps <- cpgs[cpgs$adj.P.Val <= 0.05,]

dmps$cpg_name <- rownames(dmps)
rownames(dmps) <- c(1:length(rownames(dmps)))
dmps[dmps == ""] <- NA

# write.table(dmps, file="/media/andrea/SAMSUNG/IJC/Aim2/methylation/FUNCTIONALNORMALISATION/DMPs.txt", sep="\t", quote = F, col.names = T, row.names = F)


# Adding DM status (hyper/hypo) + formatting table to export
dmps <- dplyr::as_tibble(dmps)

DMstatus <- as.data.frame(decideTests(fit2))
DMstatus$cpg_name <- rownames(DMstatus)

rownames(DMstatus) <- c(1:length(rownames(DMstatus)))
colnames(DMstatus)[1] <- "DMstat"

DMstatus <- dplyr::as_tibble(DMstatus)
joined <- inner_join(dmps, DMstatus)

# write.table(joined, file="/media/andrea/SAMSUNG/IJC/Aim2/methylation/FUNCTIONALNORMALISATION/DMPs_with_DMstatus.txt", sep="\t", quote = F, col.names = T, row.names = F)

























































