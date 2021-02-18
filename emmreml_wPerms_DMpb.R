
## Load packages ##
.libPaths( c("/data/tunglab/crc/software/rLibs_3.6.1", .libPaths() ) )
workDir <- "/data/tunglab/crc/macaqueRNA/sgeDM/"

#bring in two arguments, cluster "number" and number of iterations
args <- commandArgs(trailingOnly = TRUE)
print(args)

library(ggplot2)
library(limma)
library(edgeR)
library(EMMREML)
library(cobs)
library(cowplot)
library(snakecase)
library(stringr)
#library(qvalue)
library(scales)
library(reshape2)
#library(topGO)
#library(biomaRt)
#library(parallel)
#library(mashr)
#library(ggcorrplot)

nCluster <- as.numeric(args[1])
iters <- as.numeric(args[2])

print(nCluster)

figDir <- "/data/tunglab/crc/macaqueRNA/sgeDM/figures/"
workDir <- "/data/tunglab/crc/macaqueRNA/sgeDM/"


perm.fdr=function(input_df,perm_df,Pvals_col_name,name){
  
  pvals_index=which(colnames(input_df)==Pvals_col_name)
  ro<-input_df[order(input_df[,pvals_index]),]
  p_obs <- data.frame(pvalue=ro[,pvals_index])
  p_vector<-matrix(as.matrix(perm_df),ncol=1)
  p_vector=data.frame(p_vector[order(p_vector)])
  
  F<-p_obs[,1]
  F_o<-p_obs[,1]
  pi_hat<-p_obs[,1]
  
  j=1
  observed<-length(p_obs[,1])
  randoms<-length(p_vector[,1])
  
  for(i in 1:observed)
  {
    repeat
    {
      if((p_vector[j,1]<p_obs[i,1])&j<randoms){j<-j+1}else{break}
    }
    F[i]=i/observed
    F_o[i]=(j-1)/randoms
    if(F_o[i]<1){pi_hat[i]=(1-F[i])/(1-F_o[i])}else{pi_hat[i]=1}
  }
  tabla <-data.frame(pi_hat,pval=p_obs[,1])
  
  tabla[1,]=c(1,0)
  last_percentile_average=mean(tabla$pi_hat[as.integer(min((length(tabla[,1])*0.99),(nrow(tabla)-1)):length(tabla[,1]))])
  tabla[nrow(tabla),]=c(last_percentile_average,1)
  constraint_matrix=as.matrix(data.frame(c(0,2),c(0,1),c(1,0)))
  f_hat<-suppressWarnings(cobs(tabla$pval,tabla$pi_hat,constraint="convex",pointwise=constraint_matrix,maxiter=1000,print.warn=FALSE,print.mesg=FALSE))
  
  f_hat_serie=f_hat$fitted
  pi_o=f_hat_serie[length(f_hat_serie)]
  pi_o=min(pi_o,1)
  pi_o=max(pi_o,0)
  
  Fdr_ST_perm=pi_o*F_o/F
  
  for(i in 1:length(p_obs[,1]))
  {
    Fdr_ST_perm[i]=pi_o*F_o[i]/F[i]
    if(i>1)
    {
      for(j in 1:(i-1))
      {
        if(Fdr_ST_perm[i-j]>Fdr_ST_perm[i]){Fdr_ST_perm[i-j]=Fdr_ST_perm[i]}else{break}
      }
    }
    if(Fdr_ST_perm[i]>1)  Fdr_ST_perm[i]=1
  }
  
  fdrs_df <-data.frame(ro,q_ST_perm=Fdr_ST_perm)
  rownames(fdrs_df)=rownames(ro)
  colnames(fdrs_df)[ncol(fdrs_df)]=paste0("fdr_",name)
  
  return(fdrs_df)
}


plotModelHists <- function(emmremlDF, folderName = "~/Downloads/", modelName = "model") {
  #modelName <- "DMpb_model_AgAsym_c02"
  #emmremlDF <- res_full_aggAsym_DMpb_c02
  #folderName <- figDir
  
  nVars <- dim(emmremlDF)[2] / 3
  
  colNums <- (2 * nVars + 1) : (3 * nVars)
  
  for (c in colNums) {
    ggplot() +
      geom_histogram(aes(x = emmremlDF[,c]), bins = 100) +
      xlab(colnames(emmremlDF)[c]) + ggtitle(paste0("emmreml run ",modelName))
    ggsave(paste0(folderName,"emmreml_hist_",modelName,"_",str_replace(colnames(emmremlDF)[c],":","_"),".png"), width = 6, height = 4)
  }
}

#number of permutations
#iters <- 100

#only run cluster 7 as a test
for (cluster in c("00","14","02","05","06","07")[nCluster]) {
  #cluster <- "07"
  
  sge_LPS_pbMean_matrix <- read.table(paste0(workDir,"data/sge_LPS_c",cluster,"_pbMean_matrix"))
  sge_LPS_DMpb_matrix <- read.table(paste0(workDir,"data/sge_LPS_c",cluster,"_DMpb_matrix"))
  
  sge_NC_pbMean_matrix <- read.table(paste0(workDir,"data/sge_NC_c",cluster,"_pbMean_matrix"))
  sge_NC_DMpb_matrix <- read.table(paste0(workDir,"data/sge_NC_c",cluster,"_DMpb_matrix"))
  
  ## remove genes w/NAas DM values for any individuals, from NC & LPS DM
  numberNAs <- vector()
  
  for (r in 1:dim(sge_NC_DMpb_matrix)[1]) {
    numberNAs[r] <- sum(is.na(sge_NC_DMpb_matrix[r,]))
  }
  sge_NC_DMpb_matrix <- sge_NC_DMpb_matrix[numberNAs < 1,]
  
  numberNAs <- vector()
  for (r in 1:dim(sge_LPS_DMpb_matrix)[1]) {
    numberNAs[r] <- sum(is.na(sge_LPS_DMpb_matrix[r,]))
  }
  sge_LPS_DMpb_matrix <- sge_LPS_DMpb_matrix[numberNAs < 1,]
  
  
  ## only use common set of genes - must exist in all matrices (8)
  geneList <- rownames(table(c(rownames(sge_NC_pbMean_matrix),
                               rownames(sge_NC_DMpb_matrix),
                               rownames(sge_LPS_pbMean_matrix),
                               rownames(sge_LPS_DMpb_matrix))))[table(c(rownames(sge_NC_pbMean_matrix),
                                                                        rownames(sge_NC_DMpb_matrix),
                                                                        rownames(sge_LPS_pbMean_matrix),
                                                                        rownames(sge_LPS_DMpb_matrix))) == 4]
  
  geneList <- geneList[geneList != "chrMT"]
  
  #for any objects that start sge_NC or sge_LP, cut down their gene list and replace 77_06 with Mo
  
  for (fObj in apropos("^sge_[NL][CP]")) {
    tmpObj <- get(fObj)
    ## subset mean and variance data frames to include only the genes present in the DM data frame ##
    tmpObj <- subset(tmpObj, rownames(tmpObj) %in% geneList)
    
    ## fix MO - 77_06 ID
    if (sum(colnames(tmpObj) == "X77_06") > 0) {
      colnames(tmpObj)[which(colnames(tmpObj) == "X77_06")] <- "Mo"
    }
    #sort the columns by ID to match metadata
    tmpObj <- tmpObj[,order(colnames(tmpObj))]
    
    #reassign the temp object to the sge_NC, etc
    assign(x = fObj, value = tmpObj)
  }
  
  
  ## generate a dataframe of NC and LPS combined ##
  sge_pbMean_matrix_nested <- merge(sge_NC_pbMean_matrix,sge_LPS_pbMean_matrix, by=0)
  rownames(sge_pbMean_matrix_nested) <- sge_pbMean_matrix_nested$Row.names
  sge_pbMean_matrix_nested$Row.names <- NULL
  
  #and raw means and variances to use in the variance model
  sge_DMpb_matrix_nested <- merge(sge_NC_DMpb_matrix,sge_LPS_DMpb_matrix, by=0)
  rownames(sge_DMpb_matrix_nested) <- sge_DMpb_matrix_nested$Row.names
  sge_DMpb_matrix_nested$Row.names <- NULL
  
  indIDs <- colnames(sge_NC_DMpb_matrix)
  
  #pick metadata based on samples in the matrix
  
  ## Here I'm doing some modeling to ask what are the effects of treatment and rank on DM ##
  ## get some metadata ##
  metadata <- read.table(paste0(workDir,"data/sge_metadata_matrix"), header = T) 
  # this file is also in SGE_DM/dataFiles
  metadata$ID <- as.character(metadata$ID)
  metadata[1,1] <- "Mo"
  
  metadata <- metadata[order(metadata$ID),]
  
  metadata <- metadata[metadata$ID %in% indIDs,]
  
  ## make sure metadata order matches data frames ##
  all.equal(metadata$ID,colnames(sge_NC_DMpb_matrix))
  
  ## double the metadata so it includes both NC and LPS ##
  metadata <- rbind(metadata,metadata)
  metadata$trt <- rep(c('NC','LPS'), times=c(length(indIDs),length(indIDs)))
  
  ## run a nested model with trt, and age and elo nested within trt ##
  ## elo and age have to be mean centered for nested model ##
  metadata$elo_centered <- scale(metadata$elo, center = T, scale = F)
  #metadata$elo_centered <- metadata$elo - mean(metadata$elo)
  metadata$age_centered <- scale(metadata$age, center = T, scale = F)
  #metadata$age_centered <- metadata$age - mean(metadata$age)
  
  behaveData <- read.table(paste0(workDir,"data/SGEII_METADATA_MASTER.txt"), header = T)
  
  #phase 1 only
  behaveData <- subset(behaveData, phase == 1)
  
  #animals we need data for
  # behave data is 2 character upper, genomic data is 4 char upper/lower
  needIDs <- to_any_case(unique(substr(metadata$ID, start = 1, stop  = 2)), case = "screaming_snake")
  
  #limit behavioral data to correct 45 individs
  behaveData <- subset(behaveData, ID %in% needIDs)
  
  #sort by ID
  behaveData <- behaveData[order(behaveData$ID),]
  
  #match size of metadata object, repeats for NC & LPS
  behaveData <- rbind(behaveData,behaveData)
  
  #get the important variables, centered
  metadata$aggAsymm_cent <- behaveData$gc.agg.Diff
  metadata$aggRec_cent <- behaveData$gc.overall.agg.rec
  metadata$groom_cent <- behaveData$gc.overall.groom
  
  #flip factors in metadata
  metadata$trt <- factor(metadata$trt, levels = c("NC","LPS"))
  
  ## regress out group effects from DMpb ##
  design <- model.matrix(~metadata$group)
  fit <-lmFit(sge_DMpb_matrix_nested,design)
  intercepts=data.frame(eBayes(fit))[,1]
  resid_DMpb=apply(residuals.MArrayLM(object=fit, sge_DMpb_matrix_nested),2,function(x){x+intercepts})
  
  ## regress out group effects from pbMean ##
  design <- model.matrix(~metadata$group)
  fit <-lmFit(sge_pbMean_matrix_nested,design)
  intercepts=data.frame(eBayes(fit))[,1]
  resid_pbMean=apply(residuals.MArrayLM(object=fit, sge_pbMean_matrix_nested),2,function(x){x+intercepts})
  
  #DMpb v pbMean instead of raw mean
  p1 <- ggplot(sge_NC_pbMean_matrix, aes(sge_NC_pbMean_matrix[,1],sge_NC_DMpb_matrix[,1])) +
    geom_point(alpha=0.5) +
    #ylim(0,250) +
    #xlim(0,50) +
    stat_smooth(method='lm', color='black', fullrange = T) +
    xlab('mean ge - individual 1') +
    ylab('variance in gene expression - individual 1') +
    theme_minimal()
  p2 <- ggplot(sge_NC_pbMean_matrix, aes(sge_NC_pbMean_matrix[,2],sge_NC_DMpb_matrix[,2])) +
    geom_point(alpha=0.5) +
    #ylim(0,250) +
    #xlim(0,50) +
    stat_smooth(method='lm', color='black', fullrange = T) +
    xlab('mean ge - individual 2') +
    ylab('variance in gene expression - individual 2') +
    theme_minimal()
  p3 <- ggplot(sge_NC_pbMean_matrix, aes(sge_NC_pbMean_matrix[,3],sge_NC_DMpb_matrix[,3])) +
    geom_point(alpha=0.5) +
    #ylim(0,250) +
    #xlim(0,50) +
    stat_smooth(method='lm', color='black', fullrange = T) +
    xlab('mean ge - individual 3') +
    ylab('variance in gene expression - individual 3') +
    theme_minimal()
  
  plot_grid(p1,p2,p3, ncol = 3)
  #ggsave(paste0(figDir,"c",cluster,"/","DMpb_v_pbMean_ge__by_individual_c",cluster,".png"))
  
  
  ###### Fit an elo model for each gene using emmreml
  
  ### modeling DMpb
  ###### Fit an elo model for each gene using emmreml to DM pseudobulk
  
  design = model.matrix(~trt+trt:elo_centered+trt:age_centered,data=metadata)
  
  #Declare object res_full to store in it the model results: beta coefficients, standard deviations and p values
  res_full=resid_DMpb[,1:(3*ncol(design))]
  colnames(res_full)[1:ncol(design)]=paste0("beta_",colnames(design))
  colnames(res_full)[(ncol(design)+1):(2*ncol(design))]=paste0("sdev_",colnames(design))
  colnames(res_full)[((2*ncol(design))+1):(3*ncol(design))]=paste0("p_value_",colnames(design))
  
  #Declare object random_effects to store in it the individual-wise u random effects
  random_effects=resid_DMpb[,1:length(indIDs)]
  colnames(random_effects)=metadata[1:length(indIDs),'ID']
  
  #Define matrix Z describing the sample-to-individual mapping
  K <- read.table(paste0(workDir,"data/SGE2kinshipmatrix.txt"))
  # correcting capitalization so it matches
  rownames(K)[34] <- 'JE11'
  colnames(K)[34] <- 'JE11'
  K=K[unique(metadata$ID),unique(metadata$ID)]
  
  Z=matrix(rep(0,nrow(metadata)*ncol(K)),nrow=nrow(metadata),ncol=ncol(K))
  rownames(Z)=metadata$ID
  colnames(Z)=colnames(K)
  for(i in 1:ncol(Z)){
    set=which(metadata$ID == colnames(Z)[i])
    Z[set,i]=1
  }
  
  for(i in 1:nrow(resid_DMpb))
  {
    emma=emmreml(y=resid_DMpb[i,],X=design,Z=as.matrix(Z),K=as.matrix(K),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
    random_effects[i,]=t(emma$uhat)
    res_full[i,]=t(c(emma$betahat,emma$varbetahat,emma$pvalbeta[,"none"]))
  }
  
  assign(value = data.frame(res_full), x = paste0("res_full_DMpb_c",cluster))
  
  print(date())
  print("finished elo model")
  
  #plotModelHists(get(paste0("res_full_DMpb_c",cluster)), folderName = paste0(figDir,"c",cluster,"/"), modelName = paste0("DMpb_model_c",cluster))
  
  ######Now the other behavioral metrics
  ######agAsymm | AgRec | Groom
  
  
  #Agonism Received
  design = model.matrix(~trt+trt:aggRec_cent + trt:age_centered, data=metadata)
  
  #Declare object res_full to store in it the model results: beta coefficients, standard deviations and p values
  res_full=resid_DMpb[,1:(3*ncol(design))]
  colnames(res_full)[1:ncol(design)]=paste0("beta_",colnames(design))
  colnames(res_full)[(ncol(design)+1):(2*ncol(design))]=paste0("sdev_",colnames(design))
  colnames(res_full)[((2*ncol(design))+1):(3*ncol(design))]=paste0("p_value_",colnames(design))
  
  #Declare object random_effects to store in it the individual-wise u random effects
  random_effects=resid_DMpb[,1:length(indIDs)]
  colnames(random_effects)=metadata[1:length(indIDs),'ID']
  
  #Define matrix Z describing the sample-to-individual mapping
  K <- read.table(paste0(workDir,"data/SGE2kinshipmatrix.txt"))
  # correcting capitalization so it matches
  rownames(K)[34] <- 'JE11'
  colnames(K)[34] <- 'JE11'
  K=K[unique(metadata$ID),unique(metadata$ID)]
  
  Z=matrix(rep(0,nrow(metadata)*ncol(K)),nrow=nrow(metadata),ncol=ncol(K))
  rownames(Z)=metadata$ID
  colnames(Z)=colnames(K)
  for(i in 1:ncol(Z)){
    set=which(metadata$ID == colnames(Z)[i])
    Z[set,i]=1
  }
  
  for(i in 1:nrow(resid_DMpb))
  {
    emma=emmreml(y=resid_DMpb[i,],X=design,Z=as.matrix(Z),K=as.matrix(K),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
    random_effects[i,]=t(emma$uhat)
    res_full[i,]=t(c(emma$betahat,emma$varbetahat,emma$pvalbeta[,"none"]))
  }
  
  assign(value = data.frame(res_full), x = paste0("res_full_aggRec_DMpb_c",cluster))
  #plotModelHists(get(paste0("res_full_aggRec_DMpb_c",cluster)), folderName = paste0(figDir,"c",cluster,"/"), modelName = paste0("DMpb_model_AgRec_c",cluster))
  print(date())
  print("finished agRec model")
  
  #Overall Grooming
  
  design = model.matrix(~trt+trt:groom_cent + trt:age_centered, data=metadata)
  
  #Declare object res_full to store in it the model results: beta coefficients, standard deviations and p values
  res_full=resid_DMpb[,1:(3*ncol(design))]
  colnames(res_full)[1:ncol(design)]=paste0("beta_",colnames(design))
  colnames(res_full)[(ncol(design)+1):(2*ncol(design))]=paste0("sdev_",colnames(design))
  colnames(res_full)[((2*ncol(design))+1):(3*ncol(design))]=paste0("p_value_",colnames(design))
  
  #Declare object random_effects to store in it the individual-wise u random effects
  random_effects=resid_DMpb[,1:length(indIDs)]
  colnames(random_effects)=metadata[1:length(indIDs),'ID']
  
  #Define matrix Z describing the sample-to-individual mapping
  K <- read.table(paste0(workDir,"data/SGE2kinshipmatrix.txt"))
  # correcting capitalization so it matches
  rownames(K)[34] <- 'JE11'
  colnames(K)[34] <- 'JE11'
  K=K[unique(metadata$ID),unique(metadata$ID)]
  
  Z=matrix(rep(0,nrow(metadata)*ncol(K)),nrow=nrow(metadata),ncol=ncol(K))
  rownames(Z)=metadata$ID
  colnames(Z)=colnames(K)
  for(i in 1:ncol(Z)){
    set=which(metadata$ID == colnames(Z)[i])
    Z[set,i]=1
  }
  
  for(i in 1:nrow(resid_DMpb))
  {
    emma=emmreml(y=resid_DMpb[i,],X=design,Z=as.matrix(Z),K=as.matrix(K),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
    random_effects[i,]=t(emma$uhat)
    res_full[i,]=t(c(emma$betahat,emma$varbetahat,emma$pvalbeta[,"none"]))
  }
  
  assign(value = data.frame(res_full), x = paste0("res_full_groom_DMpb_c",cluster))
  #plotModelHists(get(paste0("res_full_groom_DMpb_c",cluster)), folderName = paste0(figDir,"c",cluster,"/"), modelName = paste0("DMpb_model_groom_c",cluster))
  
  print(date())
  print("finished groom model")
  
  ## Permutations for emprirical null ##
  ################################################################################################################################
  ## Run iters iterations of the model after permutting Elo ratings to retrieve an empiric distribution of p-values for rank effects
  ################################################################################################################################
  
  #go through each of the 3 models (elo, agRec, groom) and generate a null
  for (model in c("elo","agRec","groom")) {
    
    #generate N sets of pvalues
    for(iter in 1:iters) {
      
      #randomize treatment
      #first for the NC samples
      metaRandTrt <- subset(metadata, trt == "NC")
      #pick from a list of 50:50 NC:LPS
      metaRandTrt$trt <- sample(metadata$trt)[1:length(metaRandTrt$trt)]
      
      #duplicate the data, but flip the NC/LPS
      metaRandLPS <- metaRandTrt
      metaRandLPS$trt <- ifelse(metaRandLPS$trt == "NC", "LPS", "NC")
      
      metaRandTrt <- rbind(metaRandTrt,metaRandLPS)
      
      #randomize behavioral variable
      if ( model == "elo" ) {
        varCol <- 6
      } else if ( model == "agRec" ) {
        varCol <- 9
      } else {
        varCol <- 10
      }
      
      #Step 1 - pull out just the NC condition
      metaRandVar <- subset(metadata, trt == "NC")
      
      #Step 2 - randomize behavoral variable
      metaRandVar[,varCol] <- sample(metaRandVar[,varCol])
      
      #Step 3 - create the "LPS" metadata, e.g. the second half of the metadata table
      #by repeating the table, but assigning the opposite treatment
      metaRandLPS <- metaRandVar
      metaRandLPS$trt <- "LPS"
      
      #Step 4 - bring them back together for a single metaRand object
      metaRandVar <- rbind(metaRandVar,metaRandLPS)
      
      ###### Fit an elo model for each gene using emmreml to DM pseudobulk
      if ( model == "elo") {
        designVar <- model.matrix(~trt + trt:elo_centered + trt:age_centered, data=metaRandVar)
        designTrt <- model.matrix(~trt + trt:elo_centered + trt:age_centered, data=metaRandTrt)
      } else if ( model == "agRec" ) {
        designVar <- model.matrix(~trt + trt:aggRec_cent + trt:age_centered, data=metaRandVar)
        designTrt <- model.matrix(~trt + trt:aggRec_cent + trt:age_centered, data=metaRandTrt)
      } else {
        designVar <- model.matrix(~trt + trt:groom_cent + trt:age_centered, data=metaRandVar)
        designTrt <- model.matrix(~trt + trt:groom_cent + trt:age_centered, data=metaRandTrt)
      }
      
      ### Treatment Model
      
      #Declare object res_null to store in it the model results: beta coefficients, standard deviations and p values
      res_null=resid_DMpb[,1:(3*ncol(designTrt))]
      colnames(res_null)[1:ncol(designTrt)]=paste0("beta_",colnames(designTrt))
      colnames(res_null)[(ncol(designTrt)+1):(2*ncol(designTrt))]=paste0("sdev_",colnames(designTrt))
      colnames(res_null)[((2*ncol(designTrt))+1):(3*ncol(designTrt))]=paste0("p_value_",colnames(designTrt))
      
      #Declare object random_effects to store in it the individual-wise u random effects
      random_effects=resid_DMpb[,1:length(indIDs)]
      colnames(random_effects)=metadata[1:length(indIDs),'ID']
      
      #Define matrix Z describing the sample-to-individual mapping
      K <- read.table(paste0(workDir,"data/SGE2kinshipmatrix.txt"))
      # correcting capitalization so it matches
      rownames(K)[34] <- 'JE11'
      colnames(K)[34] <- 'JE11'
      K=K[unique(metadata$ID),unique(metadata$ID)]
      
      Z=matrix(rep(0,nrow(metadata)*ncol(K)),nrow=nrow(metadata),ncol=ncol(K))
      rownames(Z)=metadata$ID
      colnames(Z)=colnames(K)
      for(i in 1:ncol(Z)){
        set=which(metadata$ID == colnames(Z)[i])
        Z[set,i]=1
      }
      
      #Fit a model for each gene using emmreml
      for(i in 1:nrow(resid_DMpb)) {
        emma=emmreml(y=resid_DMpb[i,],X=designTrt,Z=as.matrix(Z),K=as.matrix(K),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
        random_effects[i,]=t(emma$uhat)
        res_null[i,]=t(c(emma$betahat,emma$varbetahat,emma$pvalbeta[,"none"]))
      }
      
      #we register p-values of the associations to treatment
      if( iter == 1 ) {
        shuffled_elo_trt_pvals <- data.frame(x=res_null[,"p_value_trtLPS"])
        rownames(shuffled_elo_trt_pvals) <- rownames(res_null)
      } else {
        shuffled_elo_trt_pvals <- cbind(shuffled_elo_trt_pvals,x=res_null[,"p_value_trtLPS"])
      }
      
      ### Behavioral Variable Model
      
      #Declare object res_null to store in it the model results: beta coefficients, standard deviations and p values
      res_null=resid_DMpb[,1:(3*ncol(designVar))]
      colnames(res_null)[1:ncol(designVar)]=paste0("beta_",colnames(designVar))
      colnames(res_null)[(ncol(designVar)+1):(2*ncol(designVar))]=paste0("sdev_",colnames(designVar))
      colnames(res_null)[((2*ncol(designVar))+1):(3*ncol(designVar))]=paste0("p_value_",colnames(designVar))
      
      #Declare object random_effects to store in it the individual-wise u random effects
      random_effects=resid_DMpb[,1:length(indIDs)]
      colnames(random_effects)=metadata[1:length(indIDs),'ID']
      
      #Define matrix Z describing the sample-to-individual mapping
      K <- read.table(paste0(workDir,"data/SGE2kinshipmatrix.txt"))
      # correcting capitalization so it matches
      rownames(K)[34] <- 'JE11'
      colnames(K)[34] <- 'JE11'
      K=K[unique(metadata$ID),unique(metadata$ID)]
      
      Z=matrix(rep(0,nrow(metadata)*ncol(K)),nrow=nrow(metadata),ncol=ncol(K))
      rownames(Z)=metadata$ID
      colnames(Z)=colnames(K)
      for(i in 1:ncol(Z)){
        set=which(metadata$ID == colnames(Z)[i])
        Z[set,i]=1
      }
      
      #Fit a model for each gene using emmreml
      for(i in 1:nrow(resid_DMpb)) {
        emma=emmreml(y=resid_DMpb[i,],X=designVar,Z=as.matrix(Z),K=as.matrix(K),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
        random_effects[i,]=t(emma$uhat)
        res_null[i,]=t(c(emma$betahat,emma$varbetahat,emma$pvalbeta[,"none"]))
      }
      
      #we register p-values of the associations of the behavior at NC and at LPS alone.
      if( iter == 1 ) {
        shuffled_elos_pvals_NC <- data.frame(x=res_null[,15])
        shuffled_elos_pvals_LPS <- data.frame(x=res_null[,16])
        rownames(shuffled_elos_pvals_NC) <- rownames(res_null)
        rownames(shuffled_elos_pvals_LPS) <- rownames(res_null)
      } else {
        shuffled_elos_pvals_NC <- cbind(shuffled_elos_pvals_NC,x=res_null[,15])
        shuffled_elos_pvals_LPS <- cbind(shuffled_elos_pvals_LPS,x=res_null[,16])
      }
      
    }
    
    #after the permutations are finished, use the function "perm.fdr" to add permutation based qvalues to the output
    if ( model == "elo") {
      emmremlResults <- get(paste0("res_full_DMpb_c",cluster))
      emmremlResults <- perm.fdr(emmremlResults,shuffled_elo_trt_pvals,"p_value_trtLPS","eloTrt")
      emmremlResults <- perm.fdr(emmremlResults,shuffled_elos_pvals_NC,"p_value_trtNC.elo_centered","eloNestNC")
      emmremlResults <- perm.fdr(emmremlResults,shuffled_elos_pvals_LPS,"p_value_trtLPS.elo_centered","eloNestLPS")
      assign(value = emmremlResults, x = paste0("res_full_DMpb_c",cluster))

    } else if ( model == "agRec" ) {
      emmremlResults <- get(paste0("res_full_aggRec_DMpb_c",cluster))
      emmremlResults <- perm.fdr(emmremlResults,shuffled_elo_trt_pvals,"p_value_trtLPS","agRecTrt")
      emmremlResults <- perm.fdr(emmremlResults,shuffled_elos_pvals_NC,"p_value_trtNC.aggRec_cent","agRecNestNC")
      emmremlResults <- perm.fdr(emmremlResults,shuffled_elos_pvals_LPS,"p_value_trtLPS.aggRec_cent","agRecNestLPS")
      assign(value = emmremlResults, x = paste0("res_full_aggRec_DMpb_c",cluster))
      
    } else {
      emmremlResults <- get(paste0("res_full_groom_DMpb_c",cluster))
      emmremlResults <- perm.fdr(emmremlResults,shuffled_elo_trt_pvals,"p_value_trtLPS","groomTrt")
      emmremlResults <- perm.fdr(emmremlResults,shuffled_elos_pvals_NC,"p_value_trtNC.groom_cent","groomNestNC")
      emmremlResults <- perm.fdr(emmremlResults,shuffled_elos_pvals_LPS,"p_value_trtLPS.groom_cent","groomNestLPS")
      assign(value = emmremlResults, x = paste0("res_full_groom_DMpb_c",cluster))
      
    }

    write.table(shuffled_elo_trt_pvals, file = paste0(workDir,"permutations/shuffled_",model,"_trt_pval_c",cluster,".txt"), quote = F, row.names = F)
    write.table(shuffled_elos_pvals_NC, file = paste0(workDir,"permutations/shuffled_",model,"_LPS_pval_c",cluster,".txt"), quote = F, row.names = F)
    write.table(shuffled_elos_pvals_LPS, file = paste0(workDir,"permutations/shuffled_",model,"_NC_pval_c",cluster,".txt"), quote = F, row.names = F)
    
    write.table(emmremlResults, file=paste0(workDir,paste0("results/results_c",cluster,"_",model,"_emmreml")),row.names=T,col.names=T,quote=F,sep='\t')
    
    print(date())
    print(paste0("finished ",model," model permutations"))
    
  }
  
  
}
