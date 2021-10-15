setwd("~/Box/David_Casero_Lab/Laxmi_Yeruva/RNA_seq/Complex_ana_Final")

##---- Packages Required
library(dplyr)
library(tidyr)
#library(DataCombine)
library(DESeq2)
library(ggplot2)
library("apeglm")
library("pheatmap")
library('ashr')
library(vegan)
#library("factoextra")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("variancePartition")
library('variancePartition')
#
#---- This file was generated in Step1_RNAseq_ana.R script
ProCoding_Gene_counts_mitoremove <-  read.csv(file ="../ProCoding_Gene_counts_mitoremove_Riboremove.txt", sep="\t", header = T, row.names = 1)#Only non mitochondrial protien coding genes were considered
Cumulative_counts <- rowSums(ProCoding_Gene_counts_mitoremove[,4:97])
ProCoding_Gene_counts_mitoremove_Cum <- cbind(Cumulative_counts, ProCoding_Gene_counts_mitoremove)
ProCoding_Gene_counts_mitoremove_Cum_filtered <- subset(ProCoding_Gene_counts_mitoremove_Cum, Cumulative_counts > 100) #To consider at least one count per sample
#write.table(ProCoding_Gene_counts_mitoremove_Cum_filtered, file="ProCoding_Gene_counts_mitoremove_Riboremove_Filtered.txt", sep="\t")
ProCoding_Gene_counts_mitoremove_Cum_filtered <- read.csv(file = "../ProCoding_Gene_counts_mitoremove_Riboremove_Filtered.txt", sep="\t", row.names = 1, header = TRUE)

##################################
#---- Final DataFrames to be Used
##################################
GenesPC_final_filtered <- ProCoding_Gene_counts_mitoremove_Cum_filtered[,5:98]
maplength <- ProCoding_Gene_counts_mitoremove_Cum_filtered$MappabilityLength
Gsymbols <- ProCoding_Gene_counts_mitoremove_Cum_filtered$Gname

#################################################
#----- Run DESeq- before that remove sample and reorder
#################################################
GenesPC_final_filtered_new = select(GenesPC_final_filtered, -X62_SI_MY008)
GenesPC_final_filtered_new <- data.frame(t(GenesPC_final_filtered_new))
GenesPC_final_filtered_new <- GenesPC_final_filtered_new[order(rownames(GenesPC_final_filtered_new)),]
GenesPC_final_filtered_new <- data.frame(t(GenesPC_final_filtered_new))

meta <- read.csv(file="Metadata_final_edited.txt", sep = "\t", row.names = 1, header = T)
meta <- data.frame(t(meta))
meta_new = select(meta, -X62_SI_MY008)
meta_new <- data.frame(t(meta_new))
meta_new <- meta_new[order(rownames(meta_new)),]

library (DESeq2)
#dds <- DESeqDataSetFromMatrix(countData = GenesPC_final_filtered_new, colData = meta_new, design = ~ Condition + Tissue.type + Condition:Tissue.type)
# Compute variance-stabilized data (vsd) using DESeq2. We first create a DESeq object with a "flat" design:
rnaSeq_dds <- DESeqDataSetFromMatrix(GenesPC_final_filtered_new, colData= meta_new,design= ~1)
rnaSeq_dds <- DESeq(rnaSeq_dds)
rnaSeq_dds_vsd <- varianceStabilizingTransformation(rnaSeq_dds,blind=FALSE)
#write.table(assay(rnaSeq_dds_vsd), file="VST_normalized_coutns.txt", sep="\t")

#Run a quick check into how the data looks like samples clustering:
sampleDists <- dist(t(assay(rnaSeq_dds_vsd)))
library("RColorBrewer")
library("pheatmap")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- rownames(meta_new)
colnames(sampleDistMatrix) <- rownames(meta_new)
colors <- colorRampPalette( rev(brewer.pal(9,"Blues")))(255)
jpeg("heatMap.jpg", height = 30, width = 20, units = 'in', res = 600)
thisheat <- pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,clustering_distance_cols=sampleDists,col=colors, annotation_row= meta_new, show_rownames =FALSE,show_colnames = TRUE,legend = TRUE,legend.cex = .05)##
dev.off ()

#################################################
#----- Variance Partition analysis
#################################################
#  https://www.bioconductor.org/packages/release/bioc/vignettes/variancePartition/inst/doc/variancePartition.pdf
VSD_counts <- read.csv(file="VST_normalized_coutns.txt", sep = "\t", row.names = 1, header = T)
model2 <- ~ (1|SAC) + (1|Mouse) + (1|Condition) + (1|Tissue.type) + (1|Sex) + (1|Isolator) + (1|Immunization) + (1|Epithelial_cell)
#-As there are both random/categorical and fixed effects (SAC, and Proportion of Epithelial_cell), this is called a mixed-effect linear model:
#model2 <- ~ (1|Condition) + (1+SAC) + (1|Mouse)  + (1|Tissue.type) + (1|Sex) + (1|Isolator) + (1|Immunization) + Epithelial_cell
#Suggestion: rescale fixed effect variables.
#varPart_model2 <- fitExtractVarPartModel( assay(rnaSeq_dds_vsd), model2, meta_new )
varPart_model2 <- fitExtractVarPartModel(VSD_counts, model2, meta_new)
dim(varPart_model2)
#varPart_model2@row.names <- ProCoding_Gene_counts_mitoremove_Cum_filtered$Gname
rownames(varPart_model2) = make.names(ProCoding_Gene_counts_mitoremove_Cum_filtered$Gname, unique=TRUE)
varPart_model2_df <- data.frame(varPart_model2)
write.table(varPart_model2_df, file="varPart_model2_df.txt", sep = "\t")
#sort the terms of the model by average variance explained across all genes, so when we plot they will be sorted by overall importance:
vp2 <- sortCols( varPart_model2 )
# Violin plot
jpeg("variance_partision.jpg", height = 7, width = 7, units = 'in', res = 600)
plotVarPart( vp2  ,label.angle = 90)
dev.off ()
#sort genes based on variance explained by Condition (#Here it is a treatment)
head(varPart_model2[order(varPart_model2$Condition, decreasing=TRUE),])


# Correlation between the factors of the model using Canonical Correlation analysis:
form_cor2= ~  SAC + Mouse + Condition + Tissue.type + Sex + Isolator + Immunization + Epithelial_cell
C2 = canCorPairs(form_cor2, meta_new)
jpeg("Var_Correlations.jpg", height = 6, width = 6, units = 'in', res = 600)
plotCorrMatrix( C2 )
dev.off ()

#Screen for potentially relevant genes is to retrieve those more strongly associated (highest % of variance) with a given factor (e.g. Conidtion). 
#For instance, lets plot the top 60 genes with the highest variance explained by the severity score (pucainum):
jpeg("Top60genes_condition.jpg", height = 12, width = 6, units = 'in', res = 600)
plotPercentBars( vp2[order(varPart_model2$Condition, decreasing=TRUE)[1:60],] )
dev.off ()


#----- Use Interaction Terms specifically for Condition and Tissue Type - Remove other variables which are not having much effect ex: immunization and sex
model_sel <- ~ (Tissue.type+0|Condition) + (1|SAC) + (1|Mouse) + (1|Tissue.type) + (1|Isolator) + (1|Epithelial_cell)
varPart_model_sel <- fitExtractVarPartModel( VSD_counts, model_sel, meta_new, showWarnings=FALSE )
dim(varPart_model_sel)
rownames(varPart_model_sel) = make.names(ProCoding_Gene_counts_mitoremove_Cum_filtered$Gname, unique=TRUE)
varPart_model_sel_df <- data.frame(varPart_model_sel)
write.table(varPart_model_sel_df, file="varPart_model_sel_df.txt", sep = "\t")
#sort the terms of the model by average variance explained across all genes, so when we plot they will be sorted by overall importance:
vp2_sel <- sortCols( varPart_model_sel )
jpeg("variance_partision_InteractionModel.jpg", height = 7, width = 7, units = 'in', res = 600)
plotVarPart( vp2_sel  ,label.angle = 90)
dev.off ()
#sort genes based on variance explained by Condition (#Here it is a treatment)
head(varPart_model_sel[order(varPart_model_sel$`Tissue.typeLI/Condition`, decreasing=TRUE),])
head(varPart_model_sel[order(varPart_model_sel$`Tissue.typeSI/Condition`, decreasing=TRUE),])

#Screen for potentially relevant genes is to retrieve those more strongly associated (highest % of variance) with a given factor (e.g. Conidtion). 
#For instance, lets plot the top 60 genes with the highest variance explained by the severity score (pucainum):
jpeg("Top60genes_conditionLI.jpg", height = 12, width = 6, units = 'in', res = 600)
plotPercentBars( vp2_sel[order(varPart_model_sel$`Tissue.typeLI/Condition`, decreasing=TRUE)[1:60],] )
dev.off ()
jpeg("Top60genes_conditionSI.jpg", height = 12, width = 6, units = 'in', res = 600)
plotPercentBars( vp2_sel[order(varPart_model_sel$`Tissue.typeSI/Condition`, decreasing=TRUE)[1:60],] )
dev.off ()

#Finally, there is an easy way to fit a model and retrieve the actual residuals of the model, that could be used for downstream analysis:
modelA <- ~  (1|Sex) + Age +IG+MT
varPart_modelA <- fitVarPartModel( assay(ProtectCntGMaskFactor_vsd), modelA, sampleTableGSVA_PC, fxn=residuals )
residMatrixA = do.call(rbind, varPart_modelA)


###########################
#--Seperate Model for LI
###########################
VSD_counts_transposed <- data.frame(t(VSD_counts))
VSD_counts_transposed_meta <- cbind(VSD_counts_transposed, meta_new)
VSD_counts_transposed_meta_LI <- subset(VSD_counts_transposed_meta, Tissue.type=="LI")

LI_VSD_counts <- data.frame(t(VSD_counts_transposed_meta_LI[,1:16506]))
LI_meta_new <- VSD_counts_transposed_meta_LI[,16507:16518]

#----
LI_VSD_counts <- data.frame (t(LI_VSD_counts))
nzv_cols <- nearZeroVar(LI_VSD_counts)
if(length(nzv_cols) > 0) LI_VSD_counts_new <- LI_VSD_counts[, -nzv_cols]
LI_VSD_counts <- data.frame(t(LI_VSD_counts_new))
#LI_VSD_counts <- na.omit(LI_VSD_counts)

LI_model2 <- ~ (1|SAC) + (1|Condition) + (1|Sex) + (1|Isolator) + (1|Immunization)
varPart_LI_model2 <- fitExtractVarPartModel(LI_VSD_counts, LI_model2, LI_meta_new)
dim(varPart_LI_model2)
varPart_LI_model2_df <- data.frame(varPart_LI_model2)
write.table(varPart_LI_model2_df, file="LI_model2.txt", sep = "\t")

vp2_LI <- sortCols(varPart_LI_model2 )
jpeg("LI.jpg", height = 7, width = 5, units = 'in', res = 600)
plotVarPart( vp2_LI  ,label.angle = 90)
dev.off ()

###########################
#--Seperate Model for SI
###########################
VSD_counts_transposed_meta_SI <- subset(VSD_counts_transposed_meta, Tissue.type=="SI")

SI_VSD_counts <- data.frame(t(VSD_counts_transposed_meta_SI[,1:16506]))
SI_meta_new <- VSD_counts_transposed_meta_SI[,16507:16518]

#----
SI_VSD_counts <- data.frame (t(SI_VSD_counts))
nzv_cols <- nearZeroVar(SI_VSD_counts)
if(length(nzv_cols) > 0) SI_VSD_counts_new <- SI_VSD_counts[, -nzv_cols]
SI_VSD_counts <- data.frame(t(SI_VSD_counts_new))
#LI_VSD_counts <- na.omit(LI_VSD_counts)

SI_model2 <- ~ (1|SAC) + (1|Condition) + (1|Sex) + (1|Isolator) + (1|Immunization)

varPart_SI_model2 <- fitExtractVarPartModel(SI_VSD_counts, SI_model2, SI_meta_new)
dim(varPart_SI_model2)
varPart_SI_model2_df <- data.frame(varPart_SI_model2)
write.table(varPart_SI_model2_df, file="SI_model2.txt", sep = "\t")

vp2_SI <- sortCols(varPart_SI_model2 )
jpeg("SI.jpg", height = 7, width = 5, units = 'in', res = 600)
plotVarPart( vp2_SI  ,label.angle = 90)
dev.off ()
