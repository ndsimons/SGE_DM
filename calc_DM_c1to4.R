#### code for calculating DM from SGE animals ####
#### clusters 1-4 ####


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

#for NC and LPS separately:

#NC first
treatment <- "NC"
for (clusterNum in 1:2) {

  ## extract and assign the raw single-cell count data for each individual ## 
  ind_list <- unique(sge_pbmc_NC$ID)
  
  #clusterNum <- 1
  cNumDouble <- paste0(floor(clusterNum/10),clusterNum %% 10)
  
  for (i in ind_list){
    print(i)
    tmp <- subset(sge_pbmc_NC, subset = ID == i)
    #pseudobulk straight from sge object
    obj <- SubsetData(object = tmp, do.clean = TRUE, do.scale = TRUE, do.center = TRUE, ident.use = clusterNum)
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
    assign(value = object_pseudobulkMeans, x = paste(i,"_NC_counts_c",cNumDouble,"_pbMean", sep=""))
    
    #move the counts for JUST THIS CLUSTER to another object for mean, var, DM calcs
    tmp2 <- as.matrix(GetAssayData(obj, assay = "RNA", slot = "counts"))    
    assign(value = tmp2, x = paste(i,"_NC_counts_c",cNumDouble, sep=""))
    
    # now calculate mean, cv2, var, and DM
    meanGenes <- rowMeans(tmp2)
    CV2Genes <- apply(tmp2 , 1, var) / meanGenes^2
    varGenes <- rowVars(as.matrix(tmp2))
    names(meanGenes) <- rownames(tmp2)
    names(CV2Genes) <- rownames(tmp2)
    names(varGenes) <- rownames(tmp2)
    DMlevels <- DM(meanGenes, CV2Genes) 
    assign(value = DMlevels, x = paste(i,"_NC_counts_c",cNumDouble,"_DM", sep=""))
    assign(value = meanGenes, x = paste(i,"_NC_counts_c",cNumDouble,"_mean", sep=""))
    assign(value = varGenes, x = paste(i,"_NC_counts_c",cNumDouble,"_var", sep=""))
  }
  
  ## R doesn't like objects starting with numbers so just rename 77_06 to mo for the three objects ##
  
  i <- "77_06"
  #mo_NC_counts_DM <- `77_06_NC_counts_DM`
  assign(value = get(paste(i,"_NC_counts_c",cNumDouble,"_DM", sep="")), x = paste("mo","_NC_counts_c",cNumDouble,"_DM", sep=""))
  #rm(`77_06_NC_counts_DM`)
  
  #mo_NC_counts_mean <- `77_06_NC_counts_mean`
  assign(value = get(paste(i,"_NC_counts_c",cNumDouble,"_mean", sep="")), x = paste("mo","_NC_counts_c",cNumDouble,"_mean", sep=""))
  #rm(`77_06_NC_counts_mean`)
  
  #mo_NC_counts_var <- `77_06_NC_counts_var`
  assign(value = get(paste(i,"_NC_counts_c",cNumDouble,"_var", sep="")), x = paste("mo","_NC_counts_c",cNumDouble,"_var", sep=""))
  #rm(`77_06_NC_counts_var`)
  
  #mo_NC_counts_pbMean <- `77_06_NC_counts_pbMean`
  assign(value = get(paste(i,"_NC_counts_c",cNumDouble,"_pbMean", sep="")), x = paste("mo","_NC_counts_c",cNumDouble,"_pbMean", sep=""))
  #rm(`77_06_NC_counts_pbMean`)
  
  #mo_NC_counts <- `77_06_NC_counts`
  assign(value = get(paste(i,"_NC_counts_c",cNumDouble, sep="")), x = paste("mo","_NC_counts_c",cNumDouble, sep=""))
  #rm(`77_06_NC_counts_DM`)
  
  #replace 77_06 in id_list
  ind_listBAK <- ind_list
  ind_list[which(ind_list == "77_06")] <- "mo"
  
  
  ## for DM, mean, variance, and pbMeans - generate matrices of values containing all individuals, and remove genes with NAs ##
  treatment <- "NC"
  #leaving pbMean out for now... cbind v. bind_rows needs to be switched)
  
  for (metric in c("DM","mean","var")) {
    #metric <- "DM"
    
    tmpOutput <- t(bind_rows(get(paste0("Ai10","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("Az15","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("Bi10","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("Bm7","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("Bw12","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("Cg16","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("Cl13","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("Dj15","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("DV2J","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("Ed8","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("Et12","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("Fw12","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("Fy11","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("Ga13","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("Gl11","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("Gq10","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("Gu10","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("Ht8","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("Ie6","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("Ik6","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("JE11","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("JVA","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("Kk13","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("Lm8","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("Lo9","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("Mg12","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("mo","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("Nn10","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("Ph7","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("Pt8","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("Pz13","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("Qt8","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("Qv5","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("Rn9","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("Rz6","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("Sd4","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("Sm8","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("Ss10","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("Ta11","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("Tm13","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("Tr13","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("Vt12","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("Wm14","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("Yr14","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                              get(paste0("Yz6","_",treatment,"_counts_c",cNumDouble,"_",metric))))
    tmpOutput <- na.omit(tmpOutput)
    
    colnames(tmpOutput) <- c("Ai10", "Az15", "Bi10", "Bm7", "Bw12", "Cg16", "Cl13", "Dj15", "DV2J", "Ed8",  "Et12", "Fw12", "Fy11", "Ga13", "Gl11", "Gq10", "Gu10", "Ht8",  "Ie6",  "Ik6", "JE11", "JVA",  "Kk13", "Lm8", "Lo9",  "Mg12", "mo",   "Nn10", "Ph7",  "Pt8",  "Pz13", "Qt8", "Qv5",  "Rn9",  "Rz6",  "Sd4", "Sm8",  "Ss10", "Ta11", "Tm13", "Tr13", "Vt12", "Wm14", "Yr14", "Yz6")
    
    assign(value = tmpOutput, x = paste0("sge_",treatment,"_",cNumDouble,"_",metric,"_matrix"))
    
    write.table(tmpOutput,file=paste0(workDir,paste0("sge_NC_c",cNumDouble,"_",metric,"_matrix")),row.names=T,col.names=T,quote=F,sep='\t')
  }
}

###Now LPS

treatment <- "LPS"
for (clusterNum in 1:4) {
  
  ## extract and assign the raw single-cell count data for each individual ## 
  ind_list <- unique(sge_pbmc_LPS$ID)
  
  #clusterNum <- 1
  cNumDouble <- paste0(floor(clusterNum/10),clusterNum %% 10)
  
  for (i in ind_list){
    print(i)
    tmp <- subset(sge_pbmc_LPS, subset = ID == i)
    #pseudobulk straight from sge object
    obj <- SubsetData(object = tmp, do.clean = TRUE, do.scale = TRUE, do.center = TRUE, ident.use = clusterNum)
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
    assign(value = object_pseudobulkMeans, x = paste(i,"_LPS_counts_c",cNumDouble,"_pbMean", sep=""))
    
    #move the counts for JUST THIS CLUSTER to another object for mean, var, DM calcs
    tmp2 <- as.matrix(GetAssayData(obj, assay = "RNA", slot = "counts"))    
    assign(value = tmp2, x = paste(i,"_LPS_counts_c",cNumDouble, sep=""))
    
    # now calculate mean, cv2, var, and DM
    meanGenes <- rowMeans(tmp2)
    CV2Genes <- apply(tmp2 , 1, var) / meanGenes^2
    varGenes <- rowVars(as.matrix(tmp2))
    names(meanGenes) <- rownames(tmp2)
    names(CV2Genes) <- rownames(tmp2)
    names(varGenes) <- rownames(tmp2)
    DMlevels <- DM(meanGenes, CV2Genes) 
    assign(value = DMlevels, x = paste(i,"_LPS_counts_c",cNumDouble,"_DM", sep=""))
    assign(value = meanGenes, x = paste(i,"_LPS_counts_c",cNumDouble,"_mean", sep=""))
    assign(value = varGenes, x = paste(i,"_LPS_counts_c",cNumDouble,"_var", sep=""))
  }
  
  ## R doesn't like objects starting with numbers so just rename 77_06 to mo for the three objects ##
  
  i <- "77_06"
  #mo_LPS_counts_DM <- `77_06_LPS_counts_DM`
  assign(value = get(paste(i,"_LPS_counts_c",cNumDouble,"_DM", sep="")), x = paste("mo","_LPS_counts_c",cNumDouble,"_DM", sep=""))
  #rm(`77_06_LPS_counts_DM`)
  
  #mo_LPS_counts_mean <- `77_06_LPS_counts_mean`
  assign(value = get(paste(i,"_LPS_counts_c",cNumDouble,"_mean", sep="")), x = paste("mo","_LPS_counts_c",cNumDouble,"_mean", sep=""))
  #rm(`77_06_LPS_counts_mean`)
  
  #mo_LPS_counts_var <- `77_06_LPS_counts_var`
  assign(value = get(paste(i,"_LPS_counts_c",cNumDouble,"_var", sep="")), x = paste("mo","_LPS_counts_c",cNumDouble,"_var", sep=""))
  #rm(`77_06_LPS_counts_var`)
  
  #mo_LPS_counts_pbMean <- `77_06_LPS_counts_pbMean`
  assign(value = get(paste(i,"_LPS_counts_c",cNumDouble,"_pbMean", sep="")), x = paste("mo","_LPS_counts_c",cNumDouble,"_pbMean", sep=""))
  #rm(`77_06_LPS_counts_pbMean`)
  
  #mo_LPS_counts <- `77_06_LPS_counts`
  assign(value = get(paste(i,"_LPS_counts_c",cNumDouble, sep="")), x = paste("mo","_LPS_counts_c",cNumDouble, sep=""))
  #rm(`77_06_LPS_counts_DM`)
  
  #replace 77_06 in id_list
  ind_listBAK <- ind_list
  ind_list[which(ind_list == "77_06")] <- "mo"
  
  
  ## for DM, mean, variance, and pbMeans - generate matrices of values containing all individuals, and remove genes with NAs ##
  treatment <- "LPS"
  #leaving pbMean out for now... cbind v. bind_rows needs to be switched)
  
  for (metric in c("DM","mean","var")) {
    #metric <- "DM"
    
    tmpOutput <- t(bind_rows(get(paste0("Ai10","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("Az15","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("Bi10","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("Bm7","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("Bw12","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("Cg16","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("Cl13","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("Dj15","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("DV2J","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("Ed8","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("Et12","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("Fw12","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("Fy11","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("Ga13","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("Gl11","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("Gq10","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("Gu10","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("Ht8","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("Ie6","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("Ik6","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("JE11","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("JVA","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("Kk13","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("Lm8","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("Lo9","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("Mg12","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("mo","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("Nn10","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("Ph7","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("Pt8","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("Pz13","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("Qt8","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("Qv5","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("Rn9","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("Rz6","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("Sd4","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("Sm8","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("Ss10","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("Ta11","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("Tm13","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("Tr13","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("Vt12","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("Wm14","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("Yr14","_",treatment,"_counts_c",cNumDouble,"_",metric)),
                             get(paste0("Yz6","_",treatment,"_counts_c",cNumDouble,"_",metric))))
    tmpOutput <- na.omit(tmpOutput)
    
    colnames(tmpOutput) <- c("Ai10", "Az15", "Bi10", "Bm7", "Bw12", "Cg16", "Cl13", "Dj15", "DV2J", "Ed8",  "Et12", "Fw12", "Fy11", "Ga13", "Gl11", "Gq10", "Gu10", "Ht8",  "Ie6",  "Ik6", "JE11", "JVA",  "Kk13", "Lm8", "Lo9",  "Mg12", "mo",   "Nn10", "Ph7",  "Pt8",  "Pz13", "Qt8", "Qv5",  "Rn9",  "Rz6",  "Sd4", "Sm8",  "Ss10", "Ta11", "Tm13", "Tr13", "Vt12", "Wm14", "Yr14", "Yz6")
    
    assign(value = tmpOutput, x = paste0("sge_",treatment,"_",cNumDouble,"_",metric,"_matrix"))
    
    write.table(tmpOutput,file=paste0(workDir,paste0("sge_LPS_c",cNumDouble,"_",metric,"_matrix")),row.names=T,col.names=T,quote=F,sep='\t')
  }
}