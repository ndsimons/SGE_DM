library(ggplot2)
library(limma)
library(edgeR)
library(EMMREML)
library(cobs)
library(cowplot)

## read in data frames with DM, mean, and variance for each gene ; these files are in SGE_DM/dataFiles ## 
sge_LPS_DM_matrix <- read.table('~/Desktop/SGE_DM/sge_LPS_DM_matrix')
sge_LPS_mean_matrix <- read.table('~/Desktop/SGE_DM/sge_LPS_mean_matrix')
sge_LPS_var_matrix <- read.table('~/Desktop/SGE_DM/sge_LPS_var_matrix')
sge_NC_DM_matrix <- read.table('~/Desktop/SGE_DM/sge_NC_DM_matrix')
sge_NC_mean_matrix <- read.table('~/Desktop/SGE_DM/sge_NC_mean_matrix')
sge_NC_var_matrix <- read.table('~/Desktop/SGE_DM/sge_NC_var_matrix')

## subset mean and variance data frames to include only the genes present in the DM data frame ##
sge_NC_mean_matrix <- subset(sge_NC_mean_matrix, rownames(sge_NC_mean_matrix) %in% rownames(sge_NC_DM_matrix))
sge_NC_var_matrix <- subset(sge_NC_var_matrix, rownames(sge_NC_var_matrix) %in% rownames(sge_NC_DM_matrix))
sge_LPS_mean_matrix <- subset(sge_LPS_mean_matrix, rownames(sge_LPS_mean_matrix) %in% rownames(sge_LPS_DM_matrix))
sge_LPS_var_matrix <- subset(sge_LPS_var_matrix, rownames(sge_LPS_var_matrix) %in% rownames(sge_LPS_DM_matrix))


## Does the DM metric adequately remove the mean/variance relationship? ##
## Plot mean vs variance for a few individuals to observe initial pattern ##
p1 <- ggplot(sge_NC_mean_matrix, aes(sge_NC_mean_matrix[,1],sge_NC_var_matrix[,1])) +
  geom_point(alpha=0.5) +
  ylim(0,1000) +
  xlim(0,50) +
  stat_smooth(method='lm', color='black', fullrange = T) +
  xlab('mean gene expression - individual 1') +
  ylab('variance in gene expression - individual 1') +
  theme_minimal()
p2 <- ggplot(sge_NC_mean_matrix, aes(sge_NC_mean_matrix[,2],sge_NC_var_matrix[,2])) +
  geom_point(alpha=0.5) +
  ylim(0,1000) +
  xlim(0,50) +
  stat_smooth(method='lm', color='black', fullrange = T) +
  xlab('mean gene expression - individual 2') +
  ylab('variance in gene expression - individual 2') +
  theme_minimal()
p3 <- ggplot(sge_NC_mean_matrix, aes(sge_NC_mean_matrix[,3],sge_NC_var_matrix[,3])) +
  geom_point(alpha=0.5) +
  ylim(0,1000) +
  xlim(0,50) +
  stat_smooth(method='lm', color='black', fullrange = T) +
  xlab('mean gene expression - individual 3') +
  ylab('variance in gene expression - individual 3') +
  theme_minimal()

plot_grid(p1,p2,p3, ncol = 3)

## You can pretty clearly see the mean variance relationship, now plot mean vs DM for same individuals ##
p1 <- ggplot(sge_NC_mean_matrix, aes(sge_NC_mean_matrix[,1],sge_NC_DM_matrix[,1])) +
  geom_point(alpha=0.5) +
  xlim(0,50) +
  stat_smooth(method='lm', color='black') +
  xlab('mean gene expression - individual 1') +
  ylab('DM - individual 1') +
  theme_minimal()
p2 <- ggplot(sge_NC_mean_matrix, aes(sge_NC_mean_matrix[,2],sge_NC_DM_matrix[,2])) +
  geom_point(alpha=0.5) +
  xlim(0,50) +
  stat_smooth(method='lm', color='black') +
  xlab('mean gene expression - individual 2') +
  ylab('DM - individual 2') +
  theme_minimal()
p3 <- ggplot(sge_NC_mean_matrix, aes(sge_NC_mean_matrix[,3],sge_NC_DM_matrix[,3])) +
  geom_point(alpha=0.5) +
  xlim(0,50) +
  stat_smooth(method='lm', color='black') +
  xlab('mean gene expression - individual 3') +
  ylab('DM - individual 3') +
  theme_minimal()

plot_grid(p1,p2,p3, ncol = 3)

## 
meanVarRs <- vector()
for (i in 1:45){
  meanVarRs[[i]] <- cor.test(sge_NC_mean_matrix[,i],sge_NC_var_matrix[,i])$estimate 
}
meanDMRs <- vector()
for (i in 1:45){
  meanDMRs[[i]] <- cor.test(sge_NC_mean_matrix[,i],sge_NC_DM_matrix[,i])$estimate 
}


ggplot() +
  geom_density(aes(meanVarRs, fill = 'mean/variance')) +
  geom_density(aes(meanDMRs, fill = 'mean/DM')) +
  xlim(-0.01,1) +
  xlab('R') +
  scale_fill_discrete(name = "")


ggplot() +
  geom_density(aes(abs(sge_NC_DM_matrix[,1]), fill='NC'), alpha=0.5) +
  geom_density(aes(abs(sge_LPS_DM_matrix[,1]), fill='LPS'), alpha=0.5) +
  xlim(0,0.4)


## Here I'm doing some modeling to ask what are the effects of treatment and rank on DM ##
## get some metadata ##
metadata <- read.table('~/Desktop/SGE_DM/sge_metadata_matrix', header = T) # this files is also in SGE_DM/dataFiles
metadata[1,1] <- 'Mo'
metadata <- metadata[order(metadata$ID),]

## make sure metadata order matches data frames ##
all.equal(metadata$ID,colnames(sge_NC_DM_matrix))

## double the metadata so it includes both NC and LPS ##
metadata <- rbind(metadata,metadata)
metadata$trt <- rep(c('NC','LPS'), times=c(45,45))

## generate a dataframe of NC and LPS combined ##
sge_DM_matrix_nested <- merge(sge_NC_DM_matrix,sge_LPS_DM_matrix, by=0)
rownames(sge_DM_matrix_nested) <- sge_DM_matrix_nested$Row.names
sge_DM_matrix_nested$Row.names <- NULL

## regress out group effects from DM ##
design <- model.matrix(~metadata$group)
fit <-lmFit(sge_DM_matrix_nested,design)
intercepts=data.frame(eBayes(fit))[,1]
resid_DM=apply(residuals.MArrayLM(object=fit, sge_DM_matrix_nested),2,function(x){x+intercepts})

## run a nested model with trt, and age and elo nested within trt ##
## elo and age have to be mean centered for nested model ##
metadata$elo_centered <- scale(metadata$elo, center = T, scale = F)
metadata$age_centered <- scale(metadata$age, center = T, scale = F)

design = model.matrix(~trt+trt:elo_centered+trt:age_centered,data=metadata)

#Declare object res_full to store in it the model results: beta coefficients, standard deviations and p values
res_full=resid_DM[,1:(3*ncol(design))]
colnames(res_full)[1:ncol(design)]=paste0("beta_",colnames(design))
colnames(res_full)[(ncol(design)+1):(2*ncol(design))]=paste0("sdev_",colnames(design))
colnames(res_full)[((2*ncol(design))+1):(3*ncol(design))]=paste0("p_value_",colnames(design))

#Declare object random_effects to store in it the individual-wise u random effects
random_effects=resid_DM[,1:45]
colnames(random_effects)=metadata[1:45,'ID']

#Define matrix Z describing the sample-to-individual mapping
K <- read.table('~/Downloads/SGE2kinshipmatrix.txt')
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
for(i in 1:nrow(resid_DM))
{
  emma=emmreml(y=resid_DM[i,],X=design,Z=as.matrix(Z),K=as.matrix(K),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
  random_effects[i,]=t(emma$uhat)
  res_full[i,]=t(c(emma$betahat,emma$varbetahat,emma$pvalbeta[,"none"]))
}


## Permutations for emprirical null ##
################################################################################################################################
## Run iters iterations of the model after permutting Elo ratings to retrieve an empiric distribution of p-values for rank effects
################################################################################################################################
metadata$flag=0
for(i in 1:45)
{
  metadata$flag[which(metadata$ID == metadata[1:45,'ID'][i])[1]]=1
}

iters = 3
cols_random<-metadata

for(iter in 1:iters)
{
  print(iter)
  
  #Permute Elo among the set of samples flagged as 1 (one sample per individual)
  cols_random$elo[which(cols_random$flag==1)]=sample(cols_random$elo[which(cols_random$flag==1)])
  
  #Complete the permuted Elos of the other samples (this way the mapping individual-to-rank is conserved in the permutations)
  for(i in 1:45)
  {
    set=which(cols_random$ID==metadata[1:45,'ID'][i] & cols_random$flag==0)
    set_ref=which(cols_random$ID==metadata[1:45,'ID'][i] & cols_random$flag==1)
    if(length(set)==1)
      cols_random$elo[set]=cols_random$elo[set_ref]
  }	
  
  #Declare null (i.e. based on permutations) nested design for fixed effects.
  design = model.matrix(~trt+trt:elo_centered+trt:age_centered,data=cols_random)
  
  #Declare object res_null to store in it the permutations' p-values:
  res_null=resid_DM[,1:(ncol(design))]
  colnames(res_null)[1:ncol(design)]=paste0("p_value_",colnames(design))
  
  #Fit a model for each gene using emmreml
  for(i in 1:nrow(resid_DM))
  {
    emma=emmreml(y=resid_DM[i,],X=design,Z=as.matrix(Z),K=as.matrix(K),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
    res_null[i,]=t(c(emma$pvalbeta[,"none"]))
  }
  
  #we register p-values of the associations to Elo at NC and LPS alone.
  if(iter==1)
  {
    shuffled_elos_pvals_NC <-data.frame(x=res_null[,"p_value_trtNC:elo_centered"])
    shuffled_elos_pvals_LPS <-data.frame(x=res_null[,"p_value_trtLPS:elo_centered"])
    
    rownames(shuffled_elos_pvals_NC)=rownames(res_null)
    rownames(shuffled_elos_pvals_LPS)=rownames(res_null)
  } else {
    shuffled_elos_pvals_NC <- cbind(shuffled_elos_pvals_NC,x=res_null[,"p_value_trtNC:elo_centered"])
    shuffled_elos_pvals_LPS <- cbind(shuffled_elos_pvals_LPS,x=res_null[,"p_value_trtLPS:elo_centered"])
  }
}

###################################################################################################
## Run iters iterations of the model after permuting Condition to retrieve an empirical distribution
## of p-values for LPS stimulation effect (at average rank, age and tissue composition)
###################################################################################################

cols_random<-metadata

for(iter in 1:iters)
{
  print(iter)
  
  #Permute Elo among the set of samples flagged as 1 (one sample per individual)
  cols_random$trt  <- sample(cols_random$trt)
  
  #Declare null (i.e. based on permutations) nested design for fixed effects.
  design = model.matrix(~trt+trt:elo_centered+trt:age_centered,data=cols_random)
  
  #Declare object res_null to store in it the permutations' p-values:
  res_null=resid_DM[,1:(ncol(design))]
  colnames(res_null)[1:ncol(design)]=paste0("p_value_",colnames(design))
  
  #Fit a model for each gene using emmreml
  for(i in 1:nrow(resid_DM))
  {
    emma=emmreml(y=resid_DM[i,],X=design,Z=as.matrix(Z),K=as.matrix(K),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
    res_null[i,]=t(c(emma$pvalbeta[,"none"]))
  }
  
  #we register p-values of the associations to Elo at NC and LPS alone.
  if(iter==1)
  {
    shuffled_elos_pvals_LPS_effect <-data.frame(x=res_null[,"p_value_trtNC"])
    rownames(shuffled_elos_pvals_LPS_effect)=rownames(res_null)
  } else {
    shuffled_elos_pvals_LPS_effect <- cbind(shuffled_elos_pvals_LPS_effect,x=res_null[,"p_value_trtNC"])
  }
}


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

### Ryan, this is the part that is breaking, let me know if it works for you ###

res_full=perm.fdr(data.frame(res_full),shuffled_elos_pvals_NC,'p_value_trtNC:elo_centered',"eloNC")
res_full=perm.fdr(data.frame(res_full),shuffled_elos_pvals_LPS,"p_value_trtLPS:elo_centered","eloLPS")
res_full=perm.fdr(data.frame(res_full),shuffled_elos_pvals_LPS_effect,"p_value_trtNC","test")
