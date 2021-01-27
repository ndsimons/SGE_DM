#### code for calculating DM from SGE animals ####

## Load packages ##
.libPaths( c("/data/tunglab/crc/software/rLibs_3.6.1", .libPaths() ) )
workDir <- "/data/tunglab/crc/macaqueRNA/sgeDM/"

#load(paste0(workDir,"generatePseudobulkMeans.R"))

library(data.table)
library(car)
library(Seurat)
library(dplyr)
library(limma)
library(edgeR)
library(EMMREML)
library(SingleCellExperiment)
library(ggplot2)
library(scater)
library(scran)

## read main SGE seurat object and seperate NC and LPS data ##
sge_pbmc <- readRDS("/data/tunglab/nds/2021/sge_dm/pbmc_integrated_sctransform_tsne_umap_run_censored_res_0.2_new.RDS")
sge_pbmc_NC <- subset(sge_pbmc, subset = trt == 'NC')
sge_pbmc_LPS <- subset(sge_pbmc, subset = trt == 'LPS')
rm(sge_pbmc)

#4) pseudobulk mean expression per gene

## extract and assign the raw single-cell count data for each individual ## 
ind_list <- unique(sge_pbmc_NC$ID)

for (i in ind_list){
  print(i)
  tmp <- subset(sge_pbmc_NC, subset = ID == i)
  #pseudobulk straight from sge object
  obj <- SubsetData(object = tmp, do.clean = TRUE, do.scale = TRUE, do.center = TRUE)
  sce <- as.SingleCellExperiment(obj, assay="RNA")
  ## QC metrics ##
  names(assays(sce))
  sce <- addPerCellQC(sce)
  sce <- addPerFeatureQC(sce)
  sce_filt <- sce
  sce_filt <- computeSumFactors(sce_filt)
  sce_filt$size_factor <- sizeFactors(sce_filt)
  colData(sce_filt)$group <- as.factor(colData(sce_filt)$group)
  sce_filt <- logNormCounts(sce_filt, pseudo_count = 1)
  normalized_data <- exprs(sce_filt)
  meta_data <- as.data.frame(colData(sce_filt))
  meta_data$sample_condition <- paste(meta_data$ID, meta_data$trt, sep="")
  sample_colname <- "sample_condition"
  IDs <- as.data.frame(meta_data)[, sample_colname]
  unique_ID_list <- as.list(unique(IDs))
  pseudobulk <- as.data.frame(lapply(unique_ID_list, FUN = function(x){rowMeans(normalized_data[,IDs == x, drop = FALSE])}))
  colnames(pseudobulk) <- unique(IDs)
  rownames(pseudobulk) <- rownames(normalized_data)
  pseudobulk <- as.data.frame(pseudobulk)
  object_pseudobulkMeans <- pseudobulk
  assign(value = object_pseudobulkMeans, x = paste(i,"_NC_counts_pbMean", sep=""))
  
  #move the counts to another object for mean, var, DM calcs
  tmp2 <- as.matrix(GetAssayData(tmp, assay = "RNA", slot = "counts"))    
  assign(value = tmp2, x = paste(i,"_NC_counts", sep=""))  
}

mo_NC_counts_pbMean <- `77_06_NC_counts_pbMean`
rm(`77_06_NC_counts_pbMean`)


sge_NC_pbMean_matrix <- cbind(Ai10_NC_counts_pbMean,Az15_NC_counts_pbMean,Bi10_NC_counts_pbMean,Bm7_NC_counts_pbMean,Bw12_NC_counts_pbMean,Cg16_NC_counts_pbMean,Cl13_NC_counts_pbMean,Dj15_NC_counts_pbMean,DV2J_NC_counts_pbMean,Ed8_NC_counts_pbMean,Et12_NC_counts_pbMean,Fw12_NC_counts_pbMean,Fy11_NC_counts_pbMean,Ga13_NC_counts_pbMean,Gl11_NC_counts_pbMean,Gq10_NC_counts_pbMean,Gu10_NC_counts_pbMean,Ht8_NC_counts_pbMean,Ie6_NC_counts_pbMean,Ik6_NC_counts_pbMean,JE11_NC_counts_pbMean,JVA_NC_counts_pbMean,Kk13_NC_counts_pbMean,Lm8_NC_counts_pbMean,Lo9_NC_counts_pbMean,Mg12_NC_counts_pbMean,mo_NC_counts_pbMean,Nn10_NC_counts_pbMean,Ph7_NC_counts_pbMean,Pt8_NC_counts_pbMean,Pz13_NC_counts_pbMean,Qt8_NC_counts_pbMean,Qv5_NC_counts_pbMean,Rn9_NC_counts_pbMean,Rz6_NC_counts_pbMean,Sd4_NC_counts_pbMean,Sm8_NC_counts_pbMean,Ss10_NC_counts_pbMean,Ta11_NC_counts_pbMean,Tm13_NC_counts_pbMean,Tr13_NC_counts_pbMean,Vt12_NC_counts_pbMean,Wm14_NC_counts_pbMean,Yr14_NC_counts_pbMean,Yz6_NC_counts_pbMean)
sge_NC_pbMean_matrix <- na.omit(sge_NC_pbMean_matrix)
colnames(sge_NC_pbMean_matrix) <- c("Ai10", "Az15", "Bi10", "Bm7", "Bw12", "Cg16", "Cl13", "Dj15", "DV2J", "Ed8",  "Et12", "Fw12", "Fy11", "Ga13", "Gl11", "Gq10", "Gu10", "Ht8",  "Ie6",  "Ik6", "JE11", "JVA",  "Kk13", "Lm8", "Lo9",  "Mg12", "mo",   "Nn10", "Ph7",  "Pt8",  "Pz13", "Qt8", "Qv5",  "Rn9",  "Rz6",  "Sd4", "Sm8",  "Ss10", "Ta11", "Tm13", "Tr13", "Vt12", "Wm14", "Yr14", "Yz6")

write.table(sge_NC_pbMean_matrix,file=paste0(workDir,'sge_NC_pbMean_matrix'),row.names=T,col.names=T,quote=F,sep='\t')

### filter low expression genes (if wanted) ###
#for (i in ls(pattern=".raw.data")){
#  tmp <- get(i)
#  tmp2 <- tmp[rowSums(tmp>0)>ncol(tmp)*0.005,]
#  assign(value = tmp2, x = paste(i,".filtered",sep=""))
#}

#### Now do the same for LPS samples ####
## extract and assign the raw single-cell count data for each individual ## 
ind_list <- unique(sge_pbmc_LPS$ID)

for (i in ind_list){
  print(i)
  tmp <- subset(sge_pbmc_LPS, subset = ID == i)
  #pseudobulk straight from sge object
  obj <- SubsetData(object = tmp, do.clean = TRUE, do.scale = TRUE, do.center = TRUE)
  sce <- as.SingleCellExperiment(obj, assay="RNA")
  ## QC metrics ##
  names(assays(sce))
  sce <- addPerCellQC(sce)
  sce <- addPerFeatureQC(sce)
  sce_filt <- sce
  sce_filt <- computeSumFactors(sce_filt)
  sce_filt$size_factor <- sizeFactors(sce_filt)
  colData(sce_filt)$group <- as.factor(colData(sce_filt)$group)
  sce_filt <- logNormCounts(sce_filt, pseudo_count = 1)
  normalized_data <- exprs(sce_filt)
  meta_data <- as.data.frame(colData(sce_filt))
  meta_data$sample_condition <- paste(meta_data$ID, meta_data$trt, sep="")
  sample_colname <- "sample_condition"
  IDs <- as.data.frame(meta_data)[, sample_colname]
  unique_ID_list <- as.list(unique(IDs))
  pseudobulk <- as.data.frame(lapply(unique_ID_list, FUN = function(x){rowMeans(normalized_data[,IDs == x, drop = FALSE])}))
  colnames(pseudobulk) <- unique(IDs)
  rownames(pseudobulk) <- rownames(normalized_data)
  pseudobulk <- as.data.frame(pseudobulk)
  object_pseudobulkMeans <- pseudobulk
  assign(value = object_pseudobulkMeans, x = paste(i,"_LPS_counts_pbMean", sep=""))
  
  #move the counts to another object for mean, var, DM calcs
  tmp2 <- as.matrix(GetAssayData(tmp, assay = "RNA", slot = "counts"))    
  assign(value = tmp2, x = paste(i,"_LPS_counts", sep=""))  
}

mo_LPS_counts_pbMean <- `77_06_LPS_counts_pbMean`
rm(`77_06_LPS_counts_pbMean`)

sge_LPS_pbMean_matrix <- cbind(Ai10_LPS_counts_pbMean,Az15_LPS_counts_pbMean,Bi10_LPS_counts_pbMean,Bm7_LPS_counts_pbMean,Bw12_LPS_counts_pbMean,Cg16_LPS_counts_pbMean,Cl13_LPS_counts_pbMean,Dj15_LPS_counts_pbMean,DV2J_LPS_counts_pbMean,Ed8_LPS_counts_pbMean,Et12_LPS_counts_pbMean,Fw12_LPS_counts_pbMean,Fy11_LPS_counts_pbMean,Ga13_LPS_counts_pbMean,Gl11_LPS_counts_pbMean,Gq10_LPS_counts_pbMean,Gu10_LPS_counts_pbMean,Ht8_LPS_counts_pbMean,Ie6_LPS_counts_pbMean,Ik6_LPS_counts_pbMean,JE11_LPS_counts_pbMean,JVA_LPS_counts_pbMean,Kk13_LPS_counts_pbMean,Lm8_LPS_counts_pbMean,Lo9_LPS_counts_pbMean,Mg12_LPS_counts_pbMean,mo_LPS_counts_pbMean,Nn10_LPS_counts_pbMean,Ph7_LPS_counts_pbMean,Pt8_LPS_counts_pbMean,Pz13_LPS_counts_pbMean,Qt8_LPS_counts_pbMean,Qv5_LPS_counts_pbMean,Rn9_LPS_counts_pbMean,Rz6_LPS_counts_pbMean,Sd4_LPS_counts_pbMean,Sm8_LPS_counts_pbMean,Ss10_LPS_counts_pbMean,Ta11_LPS_counts_pbMean,Tm13_LPS_counts_pbMean,Tr13_LPS_counts_pbMean,Vt12_LPS_counts_pbMean,Wm14_LPS_counts_pbMean,Yr14_LPS_counts_pbMean,Yz6_LPS_counts_pbMean)
sge_LPS_pbMean_matrix <- na.omit(sge_LPS_pbMean_matrix)
colnames(sge_LPS_pbMean_matrix) <- c("Ai10", "Az15", "Bi10", "Bm7", "Bw12", "Cg16", "Cl13", "Dj15", "DV2J", "Ed8",  "Et12", "Fw12", "Fy11", "Ga13", "Gl11", "Gq10", "Gu10", "Ht8",  "Ie6",  "Ik6", "JE11", "JVA",  "Kk13", "Lm8", "Lo9",  "Mg12", "mo",   "Nn10", "Ph7",  "Pt8",  "Pz13", "Qt8", "Qv5",  "Rn9",  "Rz6",  "Sd4", "Sm8",  "Ss10", "Ta11", "Tm13", "Tr13", "Vt12", "Wm14", "Yr14", "Yz6")

write.table(sge_LPS_pbMean_matrix,file=paste0(workDir,'sge_LPS_pbMean_matrix'),row.names=T,col.names=T,quote=F,sep='\t')

## interate through each individual calculating 
#1) mean expression per gene
#2) variance across all cells per gene,
#3) DM per gene. Assign these objects ##


countsFiles <- c(ls(pattern="NC_counts"))

for (i in countsFiles){
  tmp<-get(i)
  meanGenes <- rowMeans(tmp)
  CV2Genes <- apply(tmp , 1, var) / meanGenes^2
  varGenes <- rowVars(as.matrix(tmp))
  names(meanGenes) <- rownames(tmp)
  names(CV2Genes) <- rownames(tmp)
  names(varGenes) <- rownames(tmp)
  DMlevels <- DM(meanGenes, CV2Genes) 
  assign(value = DMlevels, x = paste(i,"_DM", sep=""))
  assign(value = meanGenes, x = paste(i,"_mean", sep=""))
  assign(value = varGenes, x = paste(i,"_var", sep=""))
}

## R doesn't like objects starting with numbers so just rename 77_06 to mo for the three objects ##
mo_NC_counts_DM <- `77_06_NC_counts_DM`
rm(`77_06_NC_counts_DM`)

mo_NC_counts_mean <- `77_06_NC_counts_mean`
rm(`77_06_NC_counts_mean`)

mo_NC_counts_var <- `77_06_NC_counts_var`
rm(`77_06_NC_counts_var`)


## for DM, mean, and variance - generate matrices of values containing all individuals, and remove genes with NAs ##
sge_NC_DM_matrix <- t(bind_rows(Ai10_NC_counts_DM,Az15_NC_counts_DM,Bi10_NC_counts_DM,Bm7_NC_counts_DM,Bw12_NC_counts_DM,Cg16_NC_counts_DM,Cl13_NC_counts_DM,Dj15_NC_counts_DM,DV2J_NC_counts_DM,Ed8_NC_counts_DM,Et12_NC_counts_DM,Fw12_NC_counts_DM,Fy11_NC_counts_DM,Ga13_NC_counts_DM,Gl11_NC_counts_DM,Gq10_NC_counts_DM,Gu10_NC_counts_DM,Ht8_NC_counts_DM,Ie6_NC_counts_DM,Ik6_NC_counts_DM,JE11_NC_counts_DM,JVA_NC_counts_DM,Kk13_NC_counts_DM,Lm8_NC_counts_DM,Lo9_NC_counts_DM,Mg12_NC_counts_DM,mo_NC_counts_DM,Nn10_NC_counts_DM,Ph7_NC_counts_DM,Pt8_NC_counts_DM,Pz13_NC_counts_DM,Qt8_NC_counts_DM,Qv5_NC_counts_DM,Rn9_NC_counts_DM,Rz6_NC_counts_DM,Sd4_NC_counts_DM,Sm8_NC_counts_DM,Ss10_NC_counts_DM,Ta11_NC_counts_DM,Tm13_NC_counts_DM,Tr13_NC_counts_DM,Vt12_NC_counts_DM,Wm14_NC_counts_DM,Yr14_NC_counts_DM,Yz6_NC_counts_DM))
sge_NC_DM_matrix <- na.omit(sge_NC_DM_matrix)
colnames(sge_NC_DM_matrix) <- c("Ai10", "Az15", "Bi10", "Bm7", "Bw12", "Cg16", "Cl13", "Dj15", "DV2J", "Ed8",  "Et12", "Fw12", "Fy11", "Ga13", "Gl11", "Gq10", "Gu10", "Ht8",  "Ie6",  "Ik6", "JE11", "JVA",  "Kk13", "Lm8", "Lo9",  "Mg12", "mo",   "Nn10", "Ph7",  "Pt8",  "Pz13", "Qt8", "Qv5",  "Rn9",  "Rz6",  "Sd4", "Sm8",  "Ss10", "Ta11", "Tm13", "Tr13", "Vt12", "Wm14", "Yr14", "Yz6")

sge_NC_mean_matrix <- t(bind_rows(Ai10_NC_counts_mean,Az15_NC_counts_mean,Bi10_NC_counts_mean,Bm7_NC_counts_mean,Bw12_NC_counts_mean,Cg16_NC_counts_mean,Cl13_NC_counts_mean,Dj15_NC_counts_mean,DV2J_NC_counts_mean,Ed8_NC_counts_mean,Et12_NC_counts_mean,Fw12_NC_counts_mean,Fy11_NC_counts_mean,Ga13_NC_counts_mean,Gl11_NC_counts_mean,Gq10_NC_counts_mean,Gu10_NC_counts_mean,Ht8_NC_counts_mean,Ie6_NC_counts_mean,Ik6_NC_counts_mean,JE11_NC_counts_mean,JVA_NC_counts_mean,Kk13_NC_counts_mean,Lm8_NC_counts_mean,Lo9_NC_counts_mean,Mg12_NC_counts_mean,mo_NC_counts_mean,Nn10_NC_counts_mean,Ph7_NC_counts_mean,Pt8_NC_counts_mean,Pz13_NC_counts_mean,Qt8_NC_counts_mean,Qv5_NC_counts_mean,Rn9_NC_counts_mean,Rz6_NC_counts_mean,Sd4_NC_counts_mean,Sm8_NC_counts_mean,Ss10_NC_counts_mean,Ta11_NC_counts_mean,Tm13_NC_counts_mean,Tr13_NC_counts_mean,Vt12_NC_counts_mean,Wm14_NC_counts_mean,Yr14_NC_counts_mean,Yz6_NC_counts_mean))
sge_NC_mean_matrix <- na.omit(sge_NC_mean_matrix)
colnames(sge_NC_mean_matrix) <- c("Ai10", "Az15", "Bi10", "Bm7", "Bw12", "Cg16", "Cl13", "Dj15", "DV2J", "Ed8",  "Et12", "Fw12", "Fy11", "Ga13", "Gl11", "Gq10", "Gu10", "Ht8",  "Ie6",  "Ik6", "JE11", "JVA",  "Kk13", "Lm8", "Lo9",  "Mg12", "mo",   "Nn10", "Ph7",  "Pt8",  "Pz13", "Qt8", "Qv5",  "Rn9",  "Rz6",  "Sd4", "Sm8",  "Ss10", "Ta11", "Tm13", "Tr13", "Vt12", "Wm14", "Yr14", "Yz6")

sge_NC_var_matrix <- t(bind_rows(Ai10_NC_counts_var,Az15_NC_counts_var,Bi10_NC_counts_var,Bm7_NC_counts_var,Bw12_NC_counts_var,Cg16_NC_counts_var,Cl13_NC_counts_var,Dj15_NC_counts_var,DV2J_NC_counts_var,Ed8_NC_counts_var,Et12_NC_counts_var,Fw12_NC_counts_var,Fy11_NC_counts_var,Ga13_NC_counts_var,Gl11_NC_counts_var,Gq10_NC_counts_var,Gu10_NC_counts_var,Ht8_NC_counts_var,Ie6_NC_counts_var,Ik6_NC_counts_var,JE11_NC_counts_var,JVA_NC_counts_var,Kk13_NC_counts_var,Lm8_NC_counts_var,Lo9_NC_counts_var,Mg12_NC_counts_var,mo_NC_counts_var,Nn10_NC_counts_var,Ph7_NC_counts_var,Pt8_NC_counts_var,Pz13_NC_counts_var,Qt8_NC_counts_var,Qv5_NC_counts_var,Rn9_NC_counts_var,Rz6_NC_counts_var,Sd4_NC_counts_var,Sm8_NC_counts_var,Ss10_NC_counts_var,Ta11_NC_counts_var,Tm13_NC_counts_var,Tr13_NC_counts_var,Vt12_NC_counts_var,Wm14_NC_counts_var,Yr14_NC_counts_var,Yz6_NC_counts_var))
sge_NC_var_matrix <- na.omit(sge_NC_var_matrix)
colnames(sge_NC_var_matrix) <- c("Ai10", "Az15", "Bi10", "Bm7", "Bw12", "Cg16", "Cl13", "Dj15", "DV2J", "Ed8",  "Et12", "Fw12", "Fy11", "Ga13", "Gl11", "Gq10", "Gu10", "Ht8",  "Ie6",  "Ik6", "JE11", "JVA",  "Kk13", "Lm8", "Lo9",  "Mg12", "mo",   "Nn10", "Ph7",  "Pt8",  "Pz13", "Qt8", "Qv5",  "Rn9",  "Rz6",  "Sd4", "Sm8",  "Ss10", "Ta11", "Tm13", "Tr13", "Vt12", "Wm14", "Yr14", "Yz6")



## save objects for posterity ##
write.table(sge_NC_DM_matrix,file=paste0(workDir,'sge_NC_DM_matrix'),row.names=T,col.names=T,quote=F,sep='\t')
write.table(sge_NC_mean_matrix,file=paste0(workDir,'sge_NC_mean_matrix'),row.names=T,col.names=T,quote=F,sep='\t')
write.table(sge_NC_var_matrix,file=paste0(workDir,'sge_NC_var_matrix'),row.names=T,col.names=T,quote=F,sep='\t')




### filter low expression genes (if wanted) ###
#for (i in ls(pattern=".raw.data")){
#  tmp <- get(i)
#  tmp2 <- tmp[rowSums(tmp>0)>ncol(tmp)*0.005,]
#  assign(value = tmp2, x = paste(i,".filtered",sep=""))
#}


## interate through each individual calculating 
#1) mean expression per gene
#2) variance across all cells per gene,
#3) DM per gene. Assign these objects ##
#4) pseudobulk mean expression per gene

countsFiles <- c(ls(pattern="LPS_counts"))

for (i in countsFiles){
  tmp<-get(i)
  meanGenes <- rowMeans(tmp)
  CV2Genes <- apply(tmp , 1, var) / meanGenes^2
  varGenes <- rowVars(as.matrix(tmp))
  names(meanGenes) <- rownames(tmp)
  names(CV2Genes) <- rownames(tmp)
  names(varGenes) <- rownames(tmp)
  DMlevels <- DM(meanGenes, CV2Genes) 
  assign(value = DMlevels, x = paste(i,"_DM", sep=""))
  assign(value = meanGenes, x = paste(i,"_mean", sep=""))
  assign(value = varGenes, x = paste(i,"_var", sep=""))
}

## R doesn't like objects starting with numbers so just rename 77_06 to mo for the three objects ##
mo_LPS_counts_DM <- `77_06_LPS_counts_DM`
rm(`77_06_LPS_counts_DM`)

mo_LPS_counts_mean <- `77_06_LPS_counts_mean`
rm(`77_06_LPS_counts_mean`)

mo_LPS_counts_var <- `77_06_LPS_counts_var`
rm(`77_06_LPS_counts_var`)


## for DM, mean, and variance - generate matrices of values containing all individuals, and remove genes with NAs ##
sge_LPS_DM_matrix <- t(bind_rows(Ai10_LPS_counts_DM,Az15_LPS_counts_DM,Bi10_LPS_counts_DM,Bm7_LPS_counts_DM,Bw12_LPS_counts_DM,Cg16_LPS_counts_DM,Cl13_LPS_counts_DM,Dj15_LPS_counts_DM,DV2J_LPS_counts_DM,Ed8_LPS_counts_DM,Et12_LPS_counts_DM,Fw12_LPS_counts_DM,Fy11_LPS_counts_DM,Ga13_LPS_counts_DM,Gl11_LPS_counts_DM,Gq10_LPS_counts_DM,Gu10_LPS_counts_DM,Ht8_LPS_counts_DM,Ie6_LPS_counts_DM,Ik6_LPS_counts_DM,JE11_LPS_counts_DM,JVA_LPS_counts_DM,Kk13_LPS_counts_DM,Lm8_LPS_counts_DM,Lo9_LPS_counts_DM,Mg12_LPS_counts_DM,mo_LPS_counts_DM,Nn10_LPS_counts_DM,Ph7_LPS_counts_DM,Pt8_LPS_counts_DM,Pz13_LPS_counts_DM,Qt8_LPS_counts_DM,Qv5_LPS_counts_DM,Rn9_LPS_counts_DM,Rz6_LPS_counts_DM,Sd4_LPS_counts_DM,Sm8_LPS_counts_DM,Ss10_LPS_counts_DM,Ta11_LPS_counts_DM,Tm13_LPS_counts_DM,Tr13_LPS_counts_DM,Vt12_LPS_counts_DM,Wm14_LPS_counts_DM,Yr14_LPS_counts_DM,Yz6_LPS_counts_DM))
sge_LPS_DM_matrix <- na.omit(sge_LPS_DM_matrix)
colnames(sge_LPS_DM_matrix) <- c("Ai10", "Az15", "Bi10", "Bm7", "Bw12", "Cg16", "Cl13", "Dj15", "DV2J", "Ed8",  "Et12", "Fw12", "Fy11", "Ga13", "Gl11", "Gq10", "Gu10", "Ht8",  "Ie6",  "Ik6", "JE11", "JVA",  "Kk13", "Lm8", "Lo9",  "Mg12", "mo",   "Nn10", "Ph7",  "Pt8",  "Pz13", "Qt8", "Qv5",  "Rn9",  "Rz6",  "Sd4", "Sm8",  "Ss10", "Ta11", "Tm13", "Tr13", "Vt12", "Wm14", "Yr14", "Yz6")

sge_LPS_mean_matrix <- t(bind_rows(Ai10_LPS_counts_mean,Az15_LPS_counts_mean,Bi10_LPS_counts_mean,Bm7_LPS_counts_mean,Bw12_LPS_counts_mean,Cg16_LPS_counts_mean,Cl13_LPS_counts_mean,Dj15_LPS_counts_mean,DV2J_LPS_counts_mean,Ed8_LPS_counts_mean,Et12_LPS_counts_mean,Fw12_LPS_counts_mean,Fy11_LPS_counts_mean,Ga13_LPS_counts_mean,Gl11_LPS_counts_mean,Gq10_LPS_counts_mean,Gu10_LPS_counts_mean,Ht8_LPS_counts_mean,Ie6_LPS_counts_mean,Ik6_LPS_counts_mean,JE11_LPS_counts_mean,JVA_LPS_counts_mean,Kk13_LPS_counts_mean,Lm8_LPS_counts_mean,Lo9_LPS_counts_mean,Mg12_LPS_counts_mean,mo_LPS_counts_mean,Nn10_LPS_counts_mean,Ph7_LPS_counts_mean,Pt8_LPS_counts_mean,Pz13_LPS_counts_mean,Qt8_LPS_counts_mean,Qv5_LPS_counts_mean,Rn9_LPS_counts_mean,Rz6_LPS_counts_mean,Sd4_LPS_counts_mean,Sm8_LPS_counts_mean,Ss10_LPS_counts_mean,Ta11_LPS_counts_mean,Tm13_LPS_counts_mean,Tr13_LPS_counts_mean,Vt12_LPS_counts_mean,Wm14_LPS_counts_mean,Yr14_LPS_counts_mean,Yz6_LPS_counts_mean))
sge_LPS_mean_matrix <- na.omit(sge_LPS_mean_matrix)
colnames(sge_LPS_mean_matrix) <- c("Ai10", "Az15", "Bi10", "Bm7", "Bw12", "Cg16", "Cl13", "Dj15", "DV2J", "Ed8",  "Et12", "Fw12", "Fy11", "Ga13", "Gl11", "Gq10", "Gu10", "Ht8",  "Ie6",  "Ik6", "JE11", "JVA",  "Kk13", "Lm8", "Lo9",  "Mg12", "mo",   "Nn10", "Ph7",  "Pt8",  "Pz13", "Qt8", "Qv5",  "Rn9",  "Rz6",  "Sd4", "Sm8",  "Ss10", "Ta11", "Tm13", "Tr13", "Vt12", "Wm14", "Yr14", "Yz6")

sge_LPS_var_matrix <- t(bind_rows(Ai10_LPS_counts_var,Az15_LPS_counts_var,Bi10_LPS_counts_var,Bm7_LPS_counts_var,Bw12_LPS_counts_var,Cg16_LPS_counts_var,Cl13_LPS_counts_var,Dj15_LPS_counts_var,DV2J_LPS_counts_var,Ed8_LPS_counts_var,Et12_LPS_counts_var,Fw12_LPS_counts_var,Fy11_LPS_counts_var,Ga13_LPS_counts_var,Gl11_LPS_counts_var,Gq10_LPS_counts_var,Gu10_LPS_counts_var,Ht8_LPS_counts_var,Ie6_LPS_counts_var,Ik6_LPS_counts_var,JE11_LPS_counts_var,JVA_LPS_counts_var,Kk13_LPS_counts_var,Lm8_LPS_counts_var,Lo9_LPS_counts_var,Mg12_LPS_counts_var,mo_LPS_counts_var,Nn10_LPS_counts_var,Ph7_LPS_counts_var,Pt8_LPS_counts_var,Pz13_LPS_counts_var,Qt8_LPS_counts_var,Qv5_LPS_counts_var,Rn9_LPS_counts_var,Rz6_LPS_counts_var,Sd4_LPS_counts_var,Sm8_LPS_counts_var,Ss10_LPS_counts_var,Ta11_LPS_counts_var,Tm13_LPS_counts_var,Tr13_LPS_counts_var,Vt12_LPS_counts_var,Wm14_LPS_counts_var,Yr14_LPS_counts_var,Yz6_LPS_counts_var))
sge_LPS_var_matrix <- na.omit(sge_LPS_var_matrix)
colnames(sge_LPS_var_matrix) <- c("Ai10", "Az15", "Bi10", "Bm7", "Bw12", "Cg16", "Cl13", "Dj15", "DV2J", "Ed8",  "Et12", "Fw12", "Fy11", "Ga13", "Gl11", "Gq10", "Gu10", "Ht8",  "Ie6",  "Ik6", "JE11", "JVA",  "Kk13", "Lm8", "Lo9",  "Mg12", "mo",   "Nn10", "Ph7",  "Pt8",  "Pz13", "Qt8", "Qv5",  "Rn9",  "Rz6",  "Sd4", "Sm8",  "Ss10", "Ta11", "Tm13", "Tr13", "Vt12", "Wm14", "Yr14", "Yz6")



## save objects for posterity ##
write.table(sge_LPS_DM_matrix,file=paste0(workDir,'sge_LPS_DM_matrix'),row.names=T,col.names=T,quote=F,sep='\t')
write.table(sge_LPS_mean_matrix,file=paste0(workDir,'sge_LPS_mean_matrix'),row.names=T,col.names=T,quote=F,sep='\t')
write.table(sge_LPS_var_matrix,file=paste0(workDir,'sge_LPS_var_matrix'),row.names=T,col.names=T,quote=F,sep='\t')




#### At this stage I move to my local machine so I can make some plotzz ####
