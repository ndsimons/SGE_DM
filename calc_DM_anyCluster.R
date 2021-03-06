#### code for calculating DM from SGE animals ####
#### any clusters ####


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
#treatment <- "NC"
for (clusterNum in c(0:7,14)) {
#for (clusterNum in 14) {
  print(clusterNum)
  cNumDouble <- paste0(floor(clusterNum/10),clusterNum %% 10)
  
  #creating new cluster nomenclature "14" which is c1 + c4 = all CD4 cells
  if (clusterNum == 14) {
    clusterSubset <- c(1,4)
  } else {
  clusterSubset <- clusterNum
  }
  
  #first, get a list of individuals for use in both treatments
  
  for(treatment in c("NC","LPS")) {
    
    cObj <- SubsetData(object = get(paste0("sge_pbmc_",treatment)), do.clean = TRUE, do.scale = TRUE, do.center = TRUE, ident.use = clusterSubset)
    
    cMeta <- cObj@meta.data %>% as.data.table
    
    ind_list <- vector()
    
    #need more than 4 cells to make the cutoff
    cellThresh <- 5
    
    for (c in unique(cMeta$ID)) {
      if ( dim(cMeta[ID == c])[1] > (cellThresh - 1) ) {
        ind_list <- c(ind_list, c)
      }
    }
    
    assign(value = ind_list, x = paste("ind_list_",treatment, sep=""))
  }
  
  final_ind_list <- ind_list_NC[ind_list_NC %in% ind_list_LPS]
  
  
  for(treatment in c("NC","LPS")) {
    
    print(paste0(clusterNum,"_",treatment))
    
    for (i in final_ind_list){
      #print(i)
      
      if ( i == "mo") {
        i <- "77_06"
      }
      
      tmp <- subset(get(paste0("sge_pbmc_",treatment)), subset = ID == i)
      #pseudobulk straight from sge object
      obj <- SubsetData(object = tmp, do.clean = TRUE, do.scale = TRUE, do.center = TRUE, ident.use = clusterSubset)
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
      assign(value = object_pseudobulkMeans, x = paste(i,"_",treatment,"_counts_c",cNumDouble,"_pbMean", sep=""))
      
      #move the counts for JUST THIS CLUSTER to another object for mean, var, DM calcs
      tmp2 <- as.matrix(GetAssayData(obj, assay = "RNA", slot = "counts"))    
      assign(value = tmp2, x = paste(i,"_",treatment,"_counts_c",cNumDouble, sep=""))
      
      # now calculate mean, cv2, var, and DM
      meanGenes <- rowMeans(tmp2)
      CV2Genes <- apply(tmp2 , 1, var) / meanGenes^2
      varGenes <- rowVars(as.matrix(tmp2))
      names(meanGenes) <- rownames(tmp2)
      names(CV2Genes) <- rownames(tmp2)
      names(varGenes) <- rownames(tmp2)
      DMlevels <- DM(meanGenes, CV2Genes)
      
      CV2pseudoGenes <- apply(tmp2 , 1, var) / pseudobulk^2
      DMpb <- DM(pseudobulk, CV2pseudoGenes)
      
      names(DMpb) <- rownames(tmp2)
      
      assign(value = DMlevels, x = paste(i,"_",treatment,"_counts_c",cNumDouble,"_DM", sep=""))
      assign(value = meanGenes, x = paste(i,"_",treatment,"_counts_c",cNumDouble,"_mean", sep=""))
      assign(value = varGenes, x = paste(i,"_",treatment,"_counts_c",cNumDouble,"_var", sep=""))
      assign(value = DMpb, x = paste(i,"_",treatment,"_counts_c",cNumDouble,"_DMpb", sep=""))
    }
    
    ## R doesn't like objects starting with numbers so just rename 77_06 to mo for the three objects ##
    
    # CHECK FOR mo FIRST
    
    if ("mo" %in% final_ind_list) {
        
      i <- "77_06"
      #mo_NC_counts_DM <- `77_06_NC_counts_DM`
      assign(value = get(paste(i,"_",treatment,"_counts_c",cNumDouble,"_DM", sep="")), 
             x = paste("mo","_",treatment,"_counts_c",cNumDouble,"_DM", sep=""))
      #rm(`77_06_NC_counts_DM`)
      
      #mo_NC_counts_mean <- `77_06_NC_counts_mean`
      assign(value = get(paste(i,"_",treatment,"_counts_c",cNumDouble,"_mean", sep="")), 
             x = paste("mo","_",treatment,"_counts_c",cNumDouble,"_mean", sep=""))
      #rm(`77_06_NC_counts_mean`)
      
      #mo_NC_counts_var <- `77_06_NC_counts_var`
      assign(value = get(paste(i,"_",treatment,"_counts_c",cNumDouble,"_var", sep="")), 
             x = paste("mo","_",treatment,"_counts_c",cNumDouble,"_var", sep=""))
      #rm(`77_06_NC_counts_var`)
      
      #mo_NC_counts_pbMean <- `77_06_NC_counts_pbMean`
      assign(value = get(paste(i,"_",treatment,"_counts_c",cNumDouble,"_DMpb", sep="")), 
             x = paste("mo","_",treatment,"_counts_c",cNumDouble,"_DMpb", sep=""))
      #rm(`77_06_NC_counts_pbMean`)
      
      #mo_NC_counts_pbMean <- `77_06_NC_counts_pbMean`
      assign(value = get(paste(i,"_",treatment,"_counts_c",cNumDouble,"_pbMean", sep="")), 
             x = paste("mo","_",treatment,"_counts_c",cNumDouble,"_pbMean", sep=""))
      #rm(`77_06_NC_counts_pbMean`)
      
      #mo_NC_counts <- `77_06_NC_counts`
      assign(value = get(paste(i,"_",treatment,"_counts_c",cNumDouble, sep="")), 
             x = paste("mo","_",treatment,"_counts_c",cNumDouble, sep=""))
      #rm(`77_06_NC_counts_DM`)
      
      #replace 77_06 in id_list
      final_ind_listBAK <- final_ind_list
      final_ind_list[which(final_ind_list == "77_06")] <- "mo"
    }
    
    ## for DM, mean, variance, and pbMeans - generate matrices of values containing all individuals, and remove genes with NAs ##
    #treatment <- "NC"
    #leaving pbMean out for now... cbind v. bind_rows needs to be switched)
    
    for (metric in c("DM","mean","var","pbMean","DMpb")) {
      print(paste0(clusterNum," ",metric))
      
      tmpOutput <- numeric()
      
      for (j in final_ind_list) {
        tmpOutput <- cbind(tmpOutput,as.matrix(get(paste0(j,"_",treatment,"_counts_c",cNumDouble,"_",metric))))
      }
      
      naCount <- vector()
      
      for (r in 1:dim(tmpOutput)[2]) {
        naCount[r] <- sum(is.na(tmpOutput[r,]))
      }
      
      tmpOutput <- tmpOutput[naCount != dim(tmpOutput)[2],]
      
      #colnames(tmpOutput) <- c("Ai10", "Az15", "Bi10", "Bm7", "Bw12", "Cg16", "Cl13", "Dj15", "DV2J", "Ed8",  "Et12", "Fw12", "Fy11", "Ga13", "Gl11", "Gq10", "Gu10", "Ht8",  "Ie6",  "Ik6", "JE11", "JVA",  "Kk13", "Lm8", "Lo9",  "Mg12", "mo",   "Nn10", "Ph7",  "Pt8",  "Pz13", "Qt8", "Qv5",  "Rn9",  "Rz6",  "Sd4", "Sm8",  "Ss10", "Ta11", "Tm13", "Tr13", "Vt12", "Wm14", "Yr14", "Yz6")
      
      colnames(tmpOutput) <- final_ind_list
      
      
      assign(value = tmpOutput, x = paste0("sge_",treatment,"_",cNumDouble,"_",metric,"_matrix"))
      
      write.table(tmpOutput,file=paste0(workDir,paste0("sge_",treatment,"_c",cNumDouble,"_",metric,"_matrix")),row.names=T,col.names=T,quote=F,sep='\t')
    }
    
  }
  
}

