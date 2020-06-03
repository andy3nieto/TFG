################################################################################
######################### MOTIF ENRICHMENT ANALYSIS   ##########################
################################################################################

#devtools::install_github('robertamezquita/marge', ref = 'master')
options('homer_path' = "/media/andrea/SAMSUNG/IJC/Aim4/HOMER/")
options(stringsAsFactors = F)
library(marge)
library(purrr)
library(dplyr)

# The aim of this task is to identify enriched Transcription Factor Binding Motifs at regions we are currently analysing. That is to say, for now we are defining the regions to study as the following:  
#   
# 1) DMPs falling in __enhancer__ regions, regardless of being involved in interactions
# + Hyper DMPs  
# + Hypo DMPs  
# 
# 2) DMPs falling in __enhancer__ regions found in HindIII fragments involved in interactions 
# + Hyper DMPs  
# + Hypo DMPs  
# 
# 3) DMPs falling in __promoter__ regions, regardless of being involved in interactions
# + Hyper DMPs  
# + Hypo DMPs  
# 
# 4) DMPs falling in __promoter__ regions found in HindIII fragments involved in interactions 
# + Hyper DMPs  
# + Hypo DMPs 

#########################################
########## Reading in the data ##########
#########################################

# read in background regions
cpgs <- read.table("/media/andrea/SAMSUNG/IJC/Aim2/methylation/characterised_450Kcpgs.txt", header = T, stringsAsFactors = F)
bg <- cpgs[,c(2:3)]
bg$start <- bg$pos-250
bg$end <- bg$pos+250
bg <- bg[,c(1,3:4)]
bg$chr <- paste0("chr",bg$chr)


#preparing regions of interest 
enhancers <- read.table(file = "/media/andrea/SAMSUNG/IJC/integration_annotations21.04.txt", header = T, stringsAsFactors = F)
e <- unique(unlist(strsplit(enhancers$DMP[!is.na(enhancers$genes)],",")))

promoters <- read.table(file = "/media/andrea/SAMSUNG/IJC/integration_promoter_annotations22.04.txt", header = T, stringsAsFactors = F)
p <- unique(unlist(strsplit(promoters$dmp_in_promoter[promoters$interactions==TRUE],",")))


ed <- cpgs[cpgs$in_enhancer_region==TRUE & cpgs$dmp==TRUE,1:3]
ed$start <- ed$pos-250
ed$end <- ed$pos+250
edi <- ed[ed$Name %in% e,]

pd <- cpgs[cpgs$in_HindIII_promoter_region==TRUE & cpgs$dmp==TRUE,1:3]
pd$start <- pd$pos-250
pd$end <- pd$pos+250
pdi <- pd[pd$Name %in% p,]

pd <- pd[,c(1:2,4:5)] #dmps in promoters
pdi <- pdi[,c(1:2,4:5)] #dmps in promoters in interactions 
ed <- ed[,c(1:2,4:5)] #dmps in enhancers
edi <- edi[,c(1:2,4:5)] #dmps in enhancers in interactions

pd$chr <- paste0("chr",pd$chr) #dmps in promoters
pdi$chr <- paste0("chr",pdi$chr) #dmps in promoters in interactions 
ed$chr <- paste0("chr",ed$chr) #dmps in enhancers
edi$chr <- paste0("chr",edi$chr)  #dmps in enhancers in interactions



######################################################
######### Running Motif Enrichment Analysis ##########
######################################################

#selecting those DMPs which are hypermethylated
keep <- cpgs$Name[cpgs$dmp==TRUE & cpgs$DMstat=="hyper"]

#regions of interest based on those selected, i.e keep the hypermethylated dmps which are in enhancers 
region <- ed[ed$Name %in% keep ,-1]

#motif enrichment 
find_motifs_genome(region,path = "/media/andrea/SAMSUNG/IJC/Aim4/motif_ea_5_output/1.2",genome ='hg19r',
                   scan_size ="given",background =bg,
                   local_background =FALSE,only_known =TRUE,only_denovo =FALSE,
                   fdr_num =5, cores = 1,cache =100,overwrite =TRUE,keep_minimal = FALSE)


#read results
result <- read_known_results("/media/andrea/SAMSUNG/IJC/Aim4/motif_ea_5_output/1.2")

my_motifs_file <- paste0("/media/andrea/SAMSUNG/IJC/Aim4/motif_ea_5_output/1.2/","my-motifs.motif")

#write PWM motif file based on enriched motifs in given subset of regions
pwalk(.l = list(motif_pwm  = r$motif_pwm, motif_name = r$motif_name,log_odds_detection = r$log_odds_detection,consensus = r$consensus),
      .f = write_homer_motif,
      file = my_motifs_file,
      append = TRUE)

keep <- cpgs$Name[cpgs$dmp==TRUE & cpgs$DMstat=="hyper"]
region <- ed[ed$Name %in% keep,]

#identify motif instances (coordinates of motifs) based on all above
motif_instances_file <- paste0("/media/andrea/SAMSUNG/IJC/Aim4/motif_ea_5_output/1.2/","motif-instances.txt")
find_motifs_instances(x = region[,c(2:4,1)],path = motif_instances_file, genome = "hg19r", my_motifs_file)





