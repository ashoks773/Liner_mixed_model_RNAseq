#---- Useful to Make 
#1. GSEA PLOT
#2. Pathway PLOT
#3. Rank Distribution Scatter plot using Two different factors

setwd("~/Box/David_Casero_Lab/Laxmi_Yeruva/RNA_seq/Complex_ana_Final/Clustering_on_Sel_Genes/Make_Final_Figures")

#######**************************
###### Figure1
#######**************************
#1A -> Flow Chart

#1B -> Variance partition analysis /Users/sharmaa4/Box/David_Casero_Lab/Laxmi_Yeruva/RNA_seq/Complex_ana_Final/Step1_LinearModels.R (take violin plot)
#-- If I need to make my violin plot then take data from here.
# Result stored at "/Users/sharmaa4/Box/David_Casero_Lab/Laxmi_Yeruva/Linear_mixed_model_New/Final_files/MasterFile2_Variance_partition_Results.xlsx"

#1C -> HeatMap
# -> take top 300 genes Based on Location from "MasterFile2_Variance_partition_Results.xlsx" and get VSD counts to make the heatMap
# grep -w -f top300genes_basedon_TissueType.txt ../../AllSamples_Filtered_VST_NormalizedCoutns.txt > top300genes_basedon_TissueType_VSD.txt

top300_vsd <- read.csv(file = "~/Box/David_Casero_Lab/Laxmi_Yeruva/RNA_seq/Complex_ana_Final/Clustering_on_Sel_Genes/Make_Final_Figures/top300genes_basedon_TissueType_VSD.txt", sep = "\t", row.names = 1, header = T)
top300_vsd <- top300_vsd[,-1]
top300_vsd <- data.frame(t(top300_vsd))
top300_vsd <- top300_vsd[order(rownames(top300_vsd)),]


meta <- read.csv(file="~/Box/David_Casero_Lab/Laxmi_Yeruva/RNA_seq/Complex_ana_Final/Metadata_final_edited.txt", sep = "\t", row.names = 1, header = T)
meta <- data.frame(t(meta))
meta = select(meta, -X62_SI_MY008)
meta <- data.frame(t(meta))
meta <- meta[order(rownames(meta)),]

top300_vsd_meta <- cbind(meta, top300_vsd)
top300_vsd_meta <- top300_vsd_meta[order(top300_vsd_meta$Tissue.type),]

#colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
colors <- colorRampPalette(c('red', 'black', 'yellow'))
colors <- colors(20)
Group <- data.frame (top300_vsd_meta$Tissue.type)
rownames(Group) <- rownames(top300_vsd_meta)

#--- Make a HeatMap of Selected Taxa
jpeg("Top300Genes_VSD_HeatmapBlack.png", height = 12, width = 8, units = 'in', res = 600)
aheatmap(top300_vsd_meta[,13:312], col=colors,   main = "Association", 
         distfun = "spearman", hclustfun = "complete", fontsize=10, 
         scale = "column", Rowv = F, Colv = TRUE, border_color = NA, legend = TRUE, cexRow = .7, cexCol = .9, 
         annRow = Group, annLegend = TRUE)
dev.off ()

#---AlternativeHeatMap
z <- scale(top300_vsd_meta[,13:312],center=TRUE,scale=TRUE)
z <- data.frame(z)
#---- To cover the Z-score range only from -2 to +2 (what ever gene has greater then +2 or lesser than -2 make it +2 or -2)
z[z > 2 ] <- 2
z[z < -2 ] <- -2

jpeg("Colored_Top300Genes_VSD_Heatmap.png", height = 7, width = 9, units = 'in', res = 600)
#jpeg("BlackWhite_Top300Genes_VSD_Heatmap.png", height = 8, width = 8, units = 'in', res = 600)
aheatmap(z, col=colors, annRow = Group, annLegend = TRUE)
#aheatmap(z, color = 1L, annRow = Group, annLegend = TRUE)
dev.off ()

library(dplyr)
z_filtered <- z%>% filter(z %in% (-2:2))
aheatmap(z_filtered)

#######**************************
#1D -> GSEA analyis
#-- For LI and SI seperately
#--- Gene rank file was made based on the Variance explained by SAC - .rnk file (Rank File)
#--- Top 300 genes were considered based on Variacne expained by Condition (Treatment) - .gmx file

#- For GSEA
#1. Gene set database : gmx file
#2. Ranked list: rnk file
#3. Chip platform: Mouse_Gene_Symbol_Remapping_MSigDB.v7.0.chip


#---- Make your own GSEA Figure
LI_genes <- read.csv(file="GSEA_analsyis_again/LI_GSEA.GseaPreranked.1636149142687/LI_GENES.tsv", sep = "\t", row.names = 1, header = T)
SI_genes <- read.csv(file="GSEA_analsyis_again/SI_GSEA.GseaPreranked.1636148787715/SI_GENES.tsv", sep = "\t", row.names = 1, header = T)
SI_genes <- SI_genes[1:299,]
df <- data.frame(LI_genes, SI_genes)

jpeg("GSEA_analsyis_again/ES_Score_LI_SI_New.png", height = 4, width = 6, units = 'in', res = 600)
#plot (LI_genes$RUNNING.ES ~ LI_genes$RANK.IN.GENE.LIST, type="l", col="#6495ED", ylim=range(c(0,0.85)), lwd = 3, xlab="Gene Rank", ylab="Enrichment Score (ES)")
plot (LI_genes$RUNNING.ES ~ LI_genes$RANK.IN.GENE.LIST, type="l", col="blue", ylim=range(c(0,0.85)), lwd = 3, xlab="Gene Rank", ylab="Enrichment Score (ES)")
#par(new=TRUE)
#plot(SI_genes$RUNNING.ES ~ SI_genes$RANK.IN.GENE.LIST, type="l", col="darkgreen", add=TRUE)
points(SI_genes$RUNNING.ES ~ SI_genes$RANK.IN.GENE.LIST,col='darkred',type="l", lwd = 3)
dev.off ()

#--- Rank Bars
jpeg("GSEA_analsyis_again/ES_Score_RankBars.png", height = 3, width = 6, units = 'in', res = 600)
plot(LI_genes$RANK.IN.GENE.LIST,rep(1,length(LI_genes$RANK.IN.GENE.LIST)),col="white")
segments(LI_genes$RANK.IN.GENE.LIST,rep(0,length(LI_genes$RANK.IN.GENE.LIST)),LI_genes$RANK.IN.GENE.LIST,rep(1,length(LI_genes$RANK.IN.GENE.LIST)),col="blue")
segments(SI_genes$RANK.IN.GENE.LIST,rep(1,length(SI_genes$RANK.IN.GENE.LIST)),SI_genes$RANK.IN.GENE.LIST,rep(2,length(SI_genes$RANK.IN.GENE.LIST)),col="darkred")
dev.off ()

#######**************************
#### Figure2
#######**************************
#------- Plot 2A and 2D - Selected Pathway BarPLOTS for SI and LI
#head -9 ~/Box/CASERO_LAB/UAMS_PAPER/SuppTable1_Sheet5_SI_FINAL_GO_Condition_Age.csv > SI_pathway.txt
SI_pathways <- read.csv(file="SI_pathway.txt", header = T)
SI_pathways <- SI_pathways[,c(2:20,23:24)]
#--- Take cluster3,5,6 and 8
SI_pathways <- SI_pathways[,c(1,13,15:16,18:21)]
colnames(SI_pathways) <- c("Age","C3","C5","C6","C8","Condition","Description","LogP")

# library
library(ggplot2)

# create a dataset
SI_pathways_Melted <- melt(SI_pathways, id.vars = "Description")

# Grouped
jpeg("SI_pathways.png", height = 5, width = 5.7, units = 'in', res = 600)
ggplot(SI_pathways_Melted, aes(fill=variable, y=value, x=Description)) + 
  geom_bar(position="dodge", stat="identity") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + coord_flip() + ylim(0,-18)
dev.off ()


#head -15 ~/Box/CASERO_LAB/UAMS_PAPER/SuppTable1_Sheet6_LI_FINAL_GO_Condition_Age.csv > LI_pathway.txt
LI_pathways <- read.csv(file="LI_pathway.txt", header = T)
LI_pathways <- LI_pathways[,c(2:22,23:26)]
#--- Take cluster 13,15,7 and 8
LI_pathways <- LI_pathways[,c(1,6,8,18,19,21,24:25)]
colnames(LI_pathways) <- c("Age","C13","C15","C7","C8","Condition","Description","LogP")

# create a dataset
LI_pathways_Melted <- melt(LI_pathways, id.vars = "Description")

# Grouped
jpeg("LI_pathways.png", height = 7, width = 6.3, units = 'in', res = 600)
ggplot(LI_pathways_Melted, aes(fill=variable, y=value, x=Description)) + 
  geom_bar(position="dodge", stat="identity") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + coord_flip() + ylim(0,-27)
dev.off ()

#------- Plot 2C and 2E - Selected rank bi-plot for genes in the clusters on Figure 2B and 3E, if not too crowed

#---- To Check the status of Genes involved in Pathways of Interest
#-- For LI and SI seperately

#1. Make two extra columns SAC_rank (based on variance explained by SAC) and Cond_rank (based on variance explained by Condition/Treatment)
#2. Select the cluster of Interest and take all Gene ranks belonging to that cluster make scatter plot based on their ranks in SAC and Condition
#---------- Plot SI Genes
SI_ranks <- read.csv(file="Rank_bi_Plots/SI_SAC_Condtion_Rank_reverse.txt", sep ="\t", row.names = 1, header = T)
C3_SI <- subset(SI_ranks, MB_Cluster == "C3")
C3_SI$SAC_Rev_rank <- as.numeric(as.character(C3_SI$SAC_Rev_rank))
C3_SI$Cond_Rev_rank <- as.numeric(as.character(C3_SI$Cond_Rev_rank))
jpeg("Rank_bi_Plots/C3_SI_rev.png", height = 3.5, width = 4, units = 'in', res = 600)
p <- ggplot(C3_SI, aes(x=SAC_Rev_rank, y=Cond_Rev_rank)) +geom_point(colour = "orchid1",size = 3,show.legend=FALSE)+geom_vline(xintercept=9000, linetype="dotted")+geom_hline(yintercept=9000, linetype="dotted")
p + geom_label_repel(aes(label=rownames(C3_SI)),max.overlaps=20,force = 200,size=5)+xlim(0,18000)+ylim(0,18000)+theme_bw()+guides(x = "none", y = "none")+labs(x="",y="")
dev.off ()


C5_SI <- subset(SI_ranks, MB_Cluster == "C5")
C5_SI$SAC_Rev_rank <- as.numeric(as.character(C5_SI$SAC_Rev_rank))
C5_SI$Cond_Rev_rank <- as.numeric(as.character(C5_SI$Cond_Rev_rank))
jpeg("Rank_bi_Plots/C5_SI_rev.png", height = 3.5, width = 4, units = 'in', res = 600)
p <- ggplot(C5_SI, aes(x=SAC_Rev_rank, y=Cond_Rev_rank)) +geom_point(colour = "orchid1",size = 3,show.legend=FALSE)+geom_vline(xintercept=9000, linetype="dotted")+geom_hline(yintercept=9000, linetype="dotted")
p + geom_label_repel(aes(label=rownames(C5_SI)),max.overlaps=20,force = 200,size=5)+xlim(0,18000)+ylim(0,18000)+theme_bw()+guides(x = "none", y = "none")+labs(x="",y="")
dev.off ()

C6_SI <- subset(SI_ranks, MB_Cluster == "C6")
C6_SI$SAC_Rev_rank <- as.numeric(as.character(C6_SI$SAC_Rev_rank))
C6_SI$Cond_Rev_rank <- as.numeric(as.character(C6_SI$Cond_Rev_rank))
jpeg("Rank_bi_Plots/C6_SI_rev.png", height = 3.5, width = 4, units = 'in', res = 600)
p <- ggplot(C6_SI, aes(x=SAC_Rev_rank, y=Cond_Rev_rank)) +geom_point(colour = "orchid1",size = 3,show.legend=FALSE)+geom_vline(xintercept=9000, linetype="dotted")+geom_hline(yintercept=9000, linetype="dotted")
p + geom_label_repel(aes(label=rownames(C6_SI)),max.overlaps=20,force = 200,size=5)+xlim(0,18000)+ylim(0,18000)+theme_bw()+guides(x = "none", y = "none")+labs(x="",y="")
dev.off ()

C8_SI <- subset(SI_ranks, MB_Cluster == "C8")
C8_SI$SAC_Rev_rank <- as.numeric(as.character(C8_SI$SAC_Rev_rank))
C8_SI$Cond_Rev_rank <- as.numeric(as.character(C8_SI$Cond_Rev_rank))
jpeg("Rank_bi_Plots/C8_SI_rev.png", height = 3.5, width = 4, units = 'in', res = 600)
p <- ggplot(C8_SI, aes(x=SAC_Rev_rank, y=Cond_Rev_rank)) +geom_point(colour = "orchid1",size = 3,show.legend=FALSE)+geom_vline(xintercept=9000, linetype="dotted")+geom_hline(yintercept=9000, linetype="dotted")
p + geom_label_repel(aes(label=rownames(C8_SI)),max.overlaps=20,force = 200,size=5)+xlim(0,18000)+ylim(0,18000)+theme_bw()+guides(x = "none", y = "none")+labs(x="",y="")
dev.off ()

#---------- Plot LI Genes
LI_ranks <- read.csv(file="Rank_bi_Plots/LI_SAC_Condition_Rank_reverse.txt", sep ="\t", row.names = 1, header = T)
C7_LI <- subset(LI_ranks, MB_Cluster == "C7")
C7_LI$SAC_Rev_rank <- as.numeric(as.character(C7_LI$SAC_Rev_rank))
C7_LI$Cond_Rev_rank <- as.numeric(as.character(C7_LI$Cond_Rev_rank))

jpeg("Rank_bi_Plots/C7_LI_rev.png", height = 3.5, width = 4, units = 'in', res = 600)
p <- ggplot(C7_LI, aes(x=SAC_Rev_rank, y=Cond_Rev_rank)) +geom_point(colour = "#6495ED",size = 3,show.legend=FALSE)+geom_vline(xintercept=9000, linetype="dotted")+geom_hline(yintercept=9000, linetype="dotted")
p + geom_label_repel(aes(label=rownames(C7_LI)),max.overlaps=20,force = 200,size=5)+xlim(0,18000)+ylim(0,18000)+theme_bw()+guides(x = "none", y = "none")+labs(x="",y="")
dev.off ()

C8_LI <- subset(LI_ranks, MB_Cluster == "C8")
C8_LI$SAC_Rev_rank <- as.numeric(as.character(C8_LI$SAC_Rev_rank))
C8_LI$Cond_Rev_rank <- as.numeric(as.character(C8_LI$Cond_Rev_rank))
jpeg("Rank_bi_Plots/C8_LI_rev.png", height = 3.5, width = 4, units = 'in', res = 600)
p <- ggplot(C8_LI, aes(x=SAC_Rev_rank, y=Cond_Rev_rank)) +geom_point(colour = "#6495ED",size = 3,show.legend=FALSE)+geom_vline(xintercept=9000, linetype="dotted")+geom_hline(yintercept=9000, linetype="dotted")
p + geom_label_repel(aes(label=rownames(C8_LI)),max.overlaps=20,force = 200,size=5)+xlim(0,18000)+ylim(0,18000)+theme_bw()+guides(x = "none", y = "none")+labs(x="",y="")
dev.off ()

C13_LI <- subset(LI_ranks, MB_Cluster == "C13")
C13_LI$SAC_Rev_rank <- as.numeric(as.character(C13_LI$SAC_Rev_rank))
C13_LI$Cond_Rev_rank <- as.numeric(as.character(C13_LI$Cond_Rev_rank))
jpeg("Rank_bi_Plots/C13_LI_rev.png", height = 4, width = 4, units = 'in', res = 600)
p <- ggplot(C13_LI, aes(x=SAC_Rev_rank, y=Cond_Rev_rank)) +geom_point(colour = "#6495ED",size = 3,show.legend=FALSE)+geom_vline(xintercept=9000, linetype="dotted")+geom_hline(yintercept=9000, linetype="dotted")
p + geom_label_repel(aes(label=rownames(C13_LI)),max.overlaps=20,force = 200,size=5)+xlim(0,18000)+ylim(0,18000)+theme_bw()+guides(x = "none", y = "none")+labs(x="",y="")
dev.off ()

C15_LI <- subset(LI_ranks, MB_Cluster == "C15")
C15_LI$SAC_Rev_rank <- as.numeric(as.character(C15_LI$SAC_Rev_rank))
C15_LI$Cond_Rev_rank <- as.numeric(as.character(C15_LI$Cond_Rev_rank))
jpeg("Rank_bi_Plots/C15_LI_rev.png", height = 3.5, width = 4, units = 'in', res = 600)
p <- ggplot(C15_LI, aes(x=SAC_Rev_rank, y=Cond_Rev_rank)) +geom_point(colour = "#6495ED",size = 3,show.legend=FALSE)+geom_vline(xintercept=9000, linetype="dotted")+geom_hline(yintercept=9000, linetype="dotted")
p + geom_label_repel(aes(label=rownames(C15_LI)),max.overlaps=20,force = 200,size=5)+xlim(0,18000)+ylim(0,18000)+theme_bw()+guides(x = "none", y = "none")+labs(x="",y="")
dev.off ()
