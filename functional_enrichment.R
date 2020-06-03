#################################################################################################
######################### FUNCTIONAL ENRICHMENT ANALYSIS (GO + KEGG)   ##########################
#################################################################################################
# SAME PIPELINE USED FOR ALL SUBSETS 
# 1 - TARGET GENES TO HYPERMETHYLATED ENHANCERS 
# 2 - TARGET GENES TO HYPOMETHYLATED ENHANCERS 
# 3 - TARGET GENES TO HYPERMETHYLATED ENHANCERS WITH P73 BINDING MOTIF 


######### Code used based on ViSEAGO vignette


library(ViSEAGO)
library(clusterProfiler)
library(org.Hs.eg.db)

# genes of interest
selection <- read.table("/media/andrea/SAMSUNG/IJC/Aim3/WORKFLOW-31-03/genes.txt", stringsAsFactors = F)
selection <- selection$V1

#background genes to use
capt <- read.table("/media/andrea/SAMSUNG/IJC/superquery/realfinal_superquery.txt", header = T, stringsAsFactors = F)
background <- unique(capt$ensembl_gene_id[!is.na(capt$ensembl_gene_id)]) 

# ensembl database used as reference
Ensembl <- ViSEAGO::Ensembl2GO(host = "grch37.ensembl.org")
myGENE2GO <- ViSEAGO::annotate("hsapiens_gene_ensembl", Ensembl)


########################################################
######### Gene Ontology: Biological Process   ##########
########################################################

# enrichment analysis
BP <- ViSEAGO::create_topGOdata(geneSel=selection, allGenes=background, gene2GO=myGENE2GO, ont="BP", nodeSize=5)
classic <- topGO::runTest(BP, algorithm ="classic", statistic = "fisher")
BP_sResults<-ViSEAGO::merge_enrich_terms(Input=list(classic=c("BP","classic")))

# view results table
ViSEAGO::show_table(BP_sResults)

# calculating semanitic similarity
myGOs<-ViSEAGO::build_GO_SS(gene2GO=myGENE2GO, enrich_GO_terms=BP_sResults)
myGOs<-ViSEAGO::compute_SS_distances(myGOs, distance="Wang")

# multidimensional scaling of enriched go terms
ViSEAGO::MDSplot(myGOs)

# heatmap of enriched go terms
Wang_clusters_wardD2<-ViSEAGO::GOterms_heatmap(myGOs, showIC=TRUE,showGOlabels=TRUE,
                                               GO.tree=list(tree=list(distance="Wang",aggreg.method="ward.D2"),
                                                            cut=list(dynamic=list(PamStage=TRUE,PamRespectsDendro=TRUE,deepSplit=2,minClusterSize =2))),
                                               samples.tree=NULL)
ViSEAGO::show_heatmap(Wang_clusters_wardD2,"GOterms")


# display colored MDSplot
ViSEAGO::MDSplot(Wang_clusters_wardD2,"GOterms")

#compute clusters of go terms based on semantic similarity
Wang_clusters_wardD2<-ViSEAGO::compute_SS_distances(Wang_clusters_wardD2, distance=c("max", "avg","rcmax", "BMA"))

# heatmap of clusters with hclust
Wang_clusters_wardD2<-ViSEAGO::GOclusters_heatmap(Wang_clusters_wardD2,tree=list(distance="BMA",aggreg.method="ward.D2"))
ViSEAGO::show_heatmap(Wang_clusters_wardD2, "GOclusters")


#######################################################
######### Gene Ontology: Molecular Function  ##########
#######################################################

# enrichment analysis
MF <- ViSEAGO::create_topGOdata(geneSel=selection, allGenes=background, gene2GO=myGENE2GO, ont="MF", nodeSize=5)
mfclassic <- topGO::runTest(MF, algorithm ="classic", statistic = "fisher")
MR_sResults<-ViSEAGO::merge_enrich_terms( Input=list(classic=c("MF","mfclassic")))

# view results table
ViSEAGO::show_table(MR_sResults)

 


#######################################################
######### Gene Ontology: Cellular Component  ##########
#######################################################

# enrichment analysis 
CC <- ViSEAGO::create_topGOdata(geneSel=selection, allGenes=background, gene2GO=myGENE2GO, ont="CC", nodeSize=5)
CCclassic <- topGO::runTest(CC, algorithm ="classic", statistic = "fisher")
CC_sResults<-ViSEAGO::merge_enrich_terms( Input=list(classic=c("CC","CCclassic")))

# view results table
ViSEAGO::show_table(CC_sResults)




###################################
######### KEGG pathways  ##########
###################################

#map ensembl ID to ENTREZ ID - package KEGG enrichment doesn't accept ENSEMBL
e <- bitr(selection, fromType = "ENSEMBL", toType = "ENTREZID",  OrgDb = org.Hs.eg.db)

#kegg pathway enrichment  
k <- enrichKEGG(e$ENTREZID, organism = "hsa", keyType = "ncbi-geneid")





