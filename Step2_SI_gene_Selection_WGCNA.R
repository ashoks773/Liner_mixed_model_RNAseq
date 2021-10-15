#############################################################################
#--- Approach2: Highest variance for Condition + SAC then Any other Variables
#perl Sel_genes_condition_SAC.pl SI_model2_Symbols.txt
#-- Select Top genes on the basis of variance explained by SAC + CONDITION (var=0.2)- File:"TopGenes_var02_condition_SAC.txt"
#-- Get counts data for these genes for WGCNA analysis
# grep -w -f TopGenes_var02_condition_SAC.txt ../../ProCoding_Gene_counts_mitoremove_Riboremove_Filtered.txt > Top_genes_conditionSAC_Counts.txt
# Add header
# Cut Gid, counts and only SI SAMPLES: Final file: "SI_Top_genes_conditionSAC_Counts"
#-- Prepare metadata file for only SI samples
#--- Remove X62_SI_MY008 Sample- Outlier No sequencing Reads
#- From 41_S1_, 62_SI_, 63_SI_ ; _ was removed
#--- Run WGCNA analysis

LI_counts <- read.csv(file="LI_Top_genes_conditionSAC_Counts.txt", sep = "\t", row.names = 1, header = T)
LI_counts <- data.frame(t(LI_counts))
LI_counts[1:1000] = lapply(LI_counts[1:1000], as.numeric)

#########################################
# Choose a set of soft-thresholding powers
#########################################
powers = c(c(1:10), seq(from = 2, to=80, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(LI_counts, powerVector = powers, verbose = 5)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.8

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#########################################
# Calcualte similarity and plot modules
#########################################
softPower = 20;
adjacency = adjacency(LI_counts, power = softPower)

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.08)


# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 4, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
png("LI_Modules.png", height = 7, width = 7, units = 'in', res = 600)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.04,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off ()

#########################################
# Quantifying module–trait associations
#########################################
# Define numbers of genes and samples 
nGenes = ncol(LI_counts)
nSamples = nrow(LI_counts)

# Recalculate MEs with color labels
MEs0alldata = moduleEigengenes(LI_counts, dynamicColors)
MEs0 = moduleEigengenes(LI_counts, dynamicColors)$eigengenes 
MEs = orderMEs(MEs0) 
barplot(MEs$MEred,col="red",main="Red",cex.main=2)
