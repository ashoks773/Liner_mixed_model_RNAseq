setwd("~/Box/David_Casero_Lab/Laxmi_Yeruva/RNA_seq/Complex_ana_Final/Clustering_on_Sel_Genes")

#-- Genes for LI and SI were selected using Variance Partition analysis
#- Genes which shows maximum variance for Condition + SAC were considred
#- Furthermore Genes were filtered based on maximum variance of 0.20
#- This lead in total 2731 genes for LI and 2469 genes for SI
#- Counts of these genes were used for the Clustering analysis

#-- Load the SI gene counts file
SI_sel_counts <- read.csv(file="~/Box/David_Casero_Lab/Laxmi_Yeruva/Linear_mixed_models_RNAseq_data/SI_top_genes_conditionSAC_counts.txt", sep="\t", row.names = 1, header = T)
SI_sel_counts <- data.frame(t(SI_sel_counts))

SI_meta <- read.csv(file="~/Box/David_Casero_Lab/Laxmi_Yeruva/Linear_mixed_models_RNAseq_data/SI_meta.txt", sep = "\t", row.names = 1, header = T)

#---- Merge Them and Sort Them based on Experimental Group and then divide again
SI_sel_counts_meta <- merge(SI_meta, SI_sel_counts, by=0, all=F)
rownames(SI_sel_counts_meta) <- SI_sel_counts_meta$Row.names; SI_sel_counts_meta$Row.names <- NULL

#--- Sort based on Experimental Groups
SI_sel_counts_meta <- SI_sel_counts_meta[with(SI_sel_counts_meta, order(ExperimentalGroups)),]


#---- Now get the counts Data
SI_sel_counts <- data.frame (t(SI_sel_counts_meta[,13:2490]))
GeneidCluster <- row.names(SI_sel_counts)
Normalizer = rep(1,47)

Treatment <- SI_sel_counts_meta$ExperimentalGroups
Treatment <- c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,8,8,8,8) 

library ("MBCluster.Seq")
DataCluster <- RNASeq.Data(SI_sel_counts, Normalize =NULL, Treatment, GeneidCluster)
save (DataCluster, file = "DataCluster_SI.Rdata")

#Start clustering, the problem is that the number of clusters needs to be pre-selected. For 2,000 genes we can start with 20
nclusters=20
c0TargetsALL = KmeansPlus.RNASeq(DataCluster,nK=nclusters,model = "nbinom")
clsTargetsALL = Cluster.RNASeq(data= DataCluster,model="nbinom",centers= c0TargetsALL$centers,method="SA")

#plot heatmap to visualize clusters
trTargetsALL = Hybrid.Tree(data= DataCluster,cluster= clsTargetsALL$cluster,model="nbinom")

jpeg("SI_clusters.jpg", height = 7, width = 15, units = 'in', res = 600)
plotHybrid.Tree(merge=trTargetsALL,cluster=clsTargetsALL$cluster,logFC=DataCluster$logFC,tree.title=NULL)
dev.off ()

# cluster membership for each gene is stored in clsTargetsALL$cluster
# the logFC for each gene is stored in DataCluster$logFC

#--- Export Clustering Information of Genes along with thier names-
#-- For each cluster export the VSD values and lavel them with Cluster ID to make Individual Cluster Boxplot

Fc <- data.frame (DataCluster$logFC)
clusterInfo <- data.frame(clsTargetsALL$cluster)
Gname_cluter_FC <- cbind (clusterInfo, Fc)
write.table (Gname_cluter_FC, file = "SI_Gname_cluter_FC.txt", sep = "\t")
#-- Now we have Gname - clustering Info and FC values Across coditions
#-- Now take only Ganme and Cluster ID from this file and Get VSD from this file


#---- Load Genes VSD file along with Clusters - whihc includes all samples - filter only SI ones
#perl Get_cluster_info.pl SI_Gname_cluter_FC.txt ../AllSamples_Filtered_VST_NormaSIzedCoutns.txt
# Add header and C in front of Cluster ID
VSD_counts <- read.csv(file="~/Box/David_Casero_Lab/Laxmi_Yeruva/RNA_seq/Complex_ana_Final/Clustering_on_Sel_Genes/SI_cluster_genes_VSD.txt", sep="\t", header = T)
SI_VSD_counts <- VSD_counts[, grepl("_SI", colnames(VSD_counts))]
SI_VSD_counts_clusters <- cbind (VSD_counts[,1:2], SI_VSD_counts)
write.table(SI_VSD_counts_clusters, file="SI_cluster_genes_VSD_Final.txt", sep="\t")

#############################################
######## Merging Modules ####################
#############################################
SI_module_genes_vsd <- SI_VSD_counts_clusters[,c(1,3:49)]
SI_module_genes_vsd_aggre <- aggregate(SI_module_genes_vsd, by = list(SI_module_genes_vsd$Cluster), FUN = mean)
rownames(SI_module_genes_vsd_aggre) <- SI_module_genes_vsd_aggre$Group.1; SI_module_genes_vsd_aggre$Group.1 <- NULL
SI_module_genes_vsd_aggre <- data.frame(t(SI_module_genes_vsd_aggre[,2:48]))
SI_module_genes_vsd_aggre <- SI_module_genes_vsd_aggre[order(rownames(SI_module_genes_vsd_aggre)),]

SI_module_genes_vsd_aggre_meta <- cbind(SI_module_genes_vsd_aggre, SI_meta)

jpeg("SI_clusters/C1_Expression.jpg", height = 3, width = 4, units = 'in', res = 600)
ggplot(aes(y = C1, x = SAC, fill = Condition), data = SI_module_genes_vsd_aggre_meta) + geom_boxplot() + ggtitle("") + labs(x="",y="C1 Genes")
dev.off ()
jpeg("SI_clusters/C2_Expression.jpg", height = 3, width = 4, units = 'in', res = 600)
ggplot(aes(y = C2, x = SAC, fill = Condition), data = SI_module_genes_vsd_aggre_meta) + geom_boxplot() + ggtitle("") + labs(x="",y="C2 Genes")
dev.off ()
jpeg("SI_clusters/C3_Expression.jpg", height = 3, width = 4, units = 'in', res = 600)
ggplot(aes(y = C3, x = SAC, fill = Condition), data = SI_module_genes_vsd_aggre_meta) + geom_boxplot() + ggtitle("") + labs(x="",y="C3 Genes")
dev.off ()
jpeg("SI_clusters/C4_Expression.jpg", height = 3, width = 4, units = 'in', res = 600)
ggplot(aes(y = C4, x = SAC, fill = Condition), data = SI_module_genes_vsd_aggre_meta) + geom_boxplot() + ggtitle("") + labs(x="",y="C4 Genes")
dev.off ()
jpeg("SI_clusters/C5_Expression.jpg", height = 3, width = 4, units = 'in', res = 600)
ggplot(aes(y = C5, x = SAC, fill = Condition), data = SI_module_genes_vsd_aggre_meta) + geom_boxplot() + ggtitle("") + labs(x="",y="C5 Genes")
dev.off ()
jpeg("SI_clusters/C6_Expression.jpg", height = 3, width = 4, units = 'in', res = 600)
ggplot(aes(y = C6, x = SAC, fill = Condition), data = SI_module_genes_vsd_aggre_meta) + geom_boxplot() + ggtitle("") + labs(x="",y="C6 Genes")
dev.off ()
jpeg("SI_clusters/C7_Expression.jpg", height = 3, width = 4, units = 'in', res = 600)
ggplot(aes(y = C7, x = SAC, fill = Condition), data = SI_module_genes_vsd_aggre_meta) + geom_boxplot() + ggtitle("") + labs(x="",y="C7 Genes")
dev.off ()
jpeg("SI_clusters/C8_Expression.jpg", height = 3, width = 4, units = 'in', res = 600)
ggplot(aes(y = C8, x = SAC, fill = Condition), data = SI_module_genes_vsd_aggre_meta) + geom_boxplot() + ggtitle("") + labs(x="",y="C8 Genes")
dev.off ()
jpeg("SI_clusters/C9_Expression.jpg", height = 3, width = 4, units = 'in', res = 600)
ggplot(aes(y = C9, x = SAC, fill = Condition), data = SI_module_genes_vsd_aggre_meta) + geom_boxplot() + ggtitle("") + labs(x="",y="C9 Genes")
dev.off ()
jpeg("SI_clusters/C10_Expression.jpg", height = 3, width = 4, units = 'in', res = 600)
ggplot(aes(y = C10, x = SAC, fill = Condition), data = SI_module_genes_vsd_aggre_meta) + geom_boxplot() + ggtitle("") + labs(x="",y="C10 Genes")
dev.off ()
jpeg("SI_clusters/C11_Expression.jpg", height = 3, width = 4, units = 'in', res = 600)
ggplot(aes(y = C11, x = SAC, fill = Condition), data = SI_module_genes_vsd_aggre_meta) + geom_boxplot() + ggtitle("") + labs(x="",y="C11 Genes")
dev.off ()
jpeg("SI_clusters/C12_Expression.jpg", height = 3, width = 4, units = 'in', res = 600)
ggplot(aes(y = C12, x = SAC, fill = Condition), data = SI_module_genes_vsd_aggre_meta) + geom_boxplot() + ggtitle("") + labs(x="",y="C12 Genes")
dev.off ()
jpeg("SI_clusters/C13_Expression.jpg", height = 3, width = 4, units = 'in', res = 600)
ggplot(aes(y = C13, x = SAC, fill = Condition), data = SI_module_genes_vsd_aggre_meta) + geom_boxplot() + ggtitle("") + labs(x="",y="C13 Genes")
dev.off ()
jpeg("SI_clusters/C14_Expression.jpg", height = 3, width = 4, units = 'in', res = 600)
ggplot(aes(y = C14, x = SAC, fill = Condition), data = SI_module_genes_vsd_aggre_meta) + geom_boxplot() + ggtitle("") + labs(x="",y="C14 Genes")
dev.off ()
jpeg("SI_clusters/C15_Expression.jpg", height = 3, width = 4, units = 'in', res = 600)
ggplot(aes(y = C15, x = SAC, fill = Condition), data = SI_module_genes_vsd_aggre_meta) + geom_boxplot() + ggtitle("") + labs(x="",y="C15 Genes")
dev.off ()
jpeg("SI_clusters/C16_Expression.jpg", height = 3, width = 4, units = 'in', res = 600)
ggplot(aes(y = C16, x = SAC, fill = Condition), data = SI_module_genes_vsd_aggre_meta) + geom_boxplot() + ggtitle("") + labs(x="",y="C16 Genes")
dev.off ()
jpeg("SI_clusters/C17_Expression.jpg", height = 3, width = 4, units = 'in', res = 600)
ggplot(aes(y = C17, x = SAC, fill = Condition), data = SI_module_genes_vsd_aggre_meta) + geom_boxplot() + ggtitle("") + labs(x="",y="C17 Genes")
dev.off ()
jpeg("SI_clusters/C18_Expression.jpg", height = 3, width = 4, units = 'in', res = 600)
ggplot(aes(y = C18, x = SAC, fill = Condition), data = SI_module_genes_vsd_aggre_meta) + geom_boxplot() + ggtitle("") + labs(x="",y="C18 Genes")
dev.off ()
jpeg("SI_clusters/C19_Expression.jpg", height = 3, width = 4, units = 'in', res = 600)
ggplot(aes(y = C19, x = SAC, fill = Condition), data = SI_module_genes_vsd_aggre_meta) + geom_boxplot() + ggtitle("") + labs(x="",y="C19 Genes")
dev.off ()
jpeg("SI_clusters/C20_Expression.jpg", height = 3, width = 4, units = 'in', res = 600)
ggplot(aes(y = C20, x = SAC, fill = Condition), data = SI_module_genes_vsd_aggre_meta) + geom_boxplot() + ggtitle("") + labs(x="",y="C20 Genes")
dev.off ()

