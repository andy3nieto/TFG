####################################################################################################
########################### INTEGRATION: DMP + ENHANCER REGIONS + PCHiC   ##########################
####################################################################################################

library(dplyr)
library(GenomicInteractions)
library(JavierreLabR)


# Integration of DMPs with enhancer regions are done separately for hypermethylated/hypomethylated DMPs.


#########################################
########## Reading in the data ##########
#########################################

# read in enhancer coordinates
enhancers <- read.table("/media/andrea/SAMSUNG/IJC/Regulatory_Builds/general/rmarkdown/enhancer_regions.txt", stringsAsFactors = F)
colnames(enhancers) <- c("ID","chr","start","end")
enhancers <- arrange(enhancers, chr, start)
e <- makeGRangesFromDataFrame(enhancers, keep.extra.columns = T)


#read in DMP coordinates
dmps <- as_tibble(read.table("/media/andrea/SAMSUNG/IJC/Aim2/methylation/FUNCTIONALNORMALISATION/DMPs_with_DMstatus.txt", header = T, stringsAsFactors = F)) #reading dmps
dmps <- dmps[,c(4,30,1,2)]
dmps <- cbind(dmps,dmps[,4])
colnames(dmps) <- c("name","dm_status","chr","start","end")
d <- makeGRangesFromDataFrame(dmps, seqnames.field = "chr", start.field = "start", end.field = "end", keep.extra.columns = T) 

#reading in enhancer and dmp annotation table 
annot <- read.table("/media/andrea/SAMSUNG/IJC/integration_annotations21.04.txt", header = T, stringsAsFactors = F)
annot <- annot[!is.na(annot$DMP),c(1,3)]

#reading in pchic interactions
pchic <- JavierreLabR::load_interactions("/media/andrea/SAMSUNG/IJC/NB/NB_Merged_Step2.txt_seqmonk.txt")
pchic <- JavierreLabR::annotate_interactions(interactions = pchic, annotation = "/media/andrea/SAMSUNG/IJC/superquery/digest_and_probes_homo_sapiens_hg19_updated_16_02_2020.txt")

#reading in captured transcripts annotation table (ENSEMBL annotation)
superquery <- read.table("/media/andrea/SAMSUNG/IJC/superquery/realfinal_superquery.txt", header = T, stringsAsFactors = F)


######################################
########## DMPs + Enhancers ##########
######################################

#overlap of enhancer and dmp genomic coordinates
overlap <- findOverlaps(e, d)
overlap <- d[subjectHits(overlap),] #5837 dmps fall in known enhancers
dmps_e <- as.data.frame(overlap@elementMetadata)


ex <- annot[1,]
df <- cbind(ex$enhancer_ID, unlist(strsplit(ex$DMP,",")))

annot <- annot[-1,]
for(i in 1:nrow(annot)){
  
  x <- cbind(annot$enhancer_ID[i], unlist(strsplit(annot$DMP[i],",")))
  df <- rbind(df, x)
}

df <- as.data.frame(df)
colnames(df) <- c("enhancer_ID","DMP")

df <- dplyr::left_join(df,dmps,by=c("DMP"="name"))

#separate hypermethylated and hypomethylated DMPs 
hyper_enhancers <- df$enhancer_ID[df$dm_status == 1]
hypo_enhancers <- df$enhancer_ID[df$dm_status == -1]

#obtaining enhancer regions that are known to be differentially methylated, separated by hyper/hypo
hyper_enhancerregions <- enhancers[enhancers$ID %in% hyper_enhancers,-1]
hypo_enhancerregions <- enhancers[enhancers$ID %in% hypo_enhancers,-1]



###########################################
########## DM Enhancers + PCHi-C ##########
###########################################

#removing interactions involving UCEs or non-annotated fragments
dat <- pchic@elementMetadata
pchic <- pchic[dat$int == "P_OE"|dat$int == "P_P",]
dat <- dat[dat$int == "P_OE"|dat$int == "P_P",]
pchic <- pchic[!(dat$gene_I == "non-annotated"),]
dat <- dat[!(dat$gene_I == "non-annotated"),]
pchic <- pchic[dat$gene_II != "non-annotated",]
dat <- dat[!(dat$gene_II == "non-annotated"),]

###############################################################################################################################################################################################################
################################################################################### (1) HYPERMETHYLATED #######################################################################################################
###############################################################################################################################################################################################################


#################################
########## Integration ##########
#################################

write.table(hyper_enhancerregions, "/media/andrea/SAMSUNG/IJC/Aim3/WORKFLOW-31-03/auxtable_hyper.txt", col.names = F, row.names = F, quote = F, sep = "\t")
integration <- JavierreLabR::integrate_regions(interactions = pchic, regions = "/media/andrea/SAMSUNG/IJC/Aim3/WORKFLOW-31-03/auxtable_hyper.txt")

####################################################################
########## Obtaining list of target genes to DM enhancers ##########
####################################################################

# a) first handle OE-P-I and P-OE-II interactions
subset <- integration[integration@elementMetadata$int == "OE_P_I" | integration@elementMetadata$int == "P_OE_II",]
d <- subset@elementMetadata

t <- c(unlist(strsplit(d$gene_I[d$overlap_II == TRUE], ",")))
file_t <- unique(t)

# map ensembl transcript IDs to ensembl gene IDs 
file_g <- unique(superquery$ensembl_gene_id[superquery$ensembl_transcript_id %in% t])





# b) second handle P-P-I and P-P-II interactions
subset <- integration[integration@elementMetadata$int == "P_P_I" | integration@elementMetadata$int == "P_P_II",]
d <- subset@elementMetadata

t <- c(unlist(strsplit(d$gene_I[d$overlap_II == TRUE], ",")))
t2 <- c(unlist(strsplit(d$gene_II[d$overlap_I == TRUE], ",")))
t <- unique(c(t,t2))

# join a + b subsets to obtain list of target transcripts
file_t <- unique(c(file_t,t))
write(file_t, file = "/media/andrea/SAMSUNG/IJC/Aim3/WORKFLOW-31-03/transcripts_all_hyper.txt")

# map ensembl transcript IDs to ensembl gene IDs 
g <- unique(superquery$ensembl_gene_id[superquery$ensembl_transcript_id %in% t])

# join a + b subsets to obtain list of target genes
file_g <- unique(c(file_g,g))
write(file_g, file = "/media/andrea/SAMSUNG/IJC/Aim3/WORKFLOW-31-03/genes_all_hyper.txt")



##############################################################################################################################################################################################################
################################################################################### (2) HYPOMETHYLATED #######################################################################################################
##############################################################################################################################################################################################################


#################################
########## Integration ##########
#################################

write.table(hyper_enhancerregions, "/media/andrea/SAMSUNG/IJC/Aim3/WORKFLOW-31-03/auxtable_hypo.txt", col.names = F, row.names = F, quote = F, sep = "\t")
integration <- JavierreLabR::integrate_regions(interactions = pchic, regions = "/media/andrea/SAMSUNG/IJC/Aim3/WORKFLOW-31-03/auxtable_hypo.txt")

####################################################################
########## Obtaining list of target genes to DM enhancers ##########
####################################################################

# a) first handle OE-P-I and P-OE-II interactions
subset <- integration[integration@elementMetadata$int == "OE_P_I" | integration@elementMetadata$int == "P_OE_II",]
d <- subset@elementMetadata

t <- c(unlist(strsplit(d$gene_I[d$overlap_II == TRUE], ",")))
file_t <- unique(t)

# map ensembl transcript IDs to ensembl gene IDs 
file_g <- unique(superquery$ensembl_gene_id[superquery$ensembl_transcript_id %in% t])





# b) second handle P-P-I and P-P-II interactions
subset <- integration[integration@elementMetadata$int == "P_P_I" | integration@elementMetadata$int == "P_P_II",]
d <- subset@elementMetadata

t <- c(unlist(strsplit(d$gene_I[d$overlap_II == TRUE], ",")))
t2 <- c(unlist(strsplit(d$gene_II[d$overlap_I == TRUE], ",")))
t <- unique(c(t,t2))

# join a + b subsets to obtain list of target transcripts
file_t <- unique(c(file_t,t))
write(file_t, file = "/media/andrea/SAMSUNG/IJC/Aim3/WORKFLOW-31-03/transcripts_all_hypo.txt")

# map ensembl transcript IDs to ensembl gene IDs 
g <- unique(superquery$ensembl_gene_id[superquery$ensembl_transcript_id %in% t])

# join a + b subsets to obtain list of target genes
file_g <- unique(c(file_g,g))
write(file_g, file = "/media/andrea/SAMSUNG/IJC/Aim3/WORKFLOW-31-03/genes_all_hypo.txt")