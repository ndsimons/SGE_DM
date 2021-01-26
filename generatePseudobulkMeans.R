generatePseudobulkMeans <- function(clusterNum, print=FALSE){
  object <- SubsetData(object = sge_seurat_object, ident.use = clusterNum, do.clean = TRUE, do.scale = TRUE, do.center = TRUE)
  sce <- as.SingleCellExperiment(object, assay="RNA")
  ## QC metrics ##
  names(assays(sce))
  sce <- addPerCellQC(sce)
  sce <- addPerFeatureQC(sce)
  logCPM_filter <- 0.2
  keep_gene <- rowMeans(log2(calculateCPM(sce) + 1)) > logCPM_filter
  sce_filt <- sce[keep_gene,]
  sce_filt <- computeSumFactors(sce_filt)
  sce_filt$size_factor <- sizeFactors(sce_filt)
  colData(sce_filt)$group <- as.factor(colData(sce_filt)$group)
  sce_filt <- sce_filt[, (sizeFactors(sce_filt) < 8 & sizeFactors(sce_filt) > 0.125)]
  sce_filt <- logNormCounts(sce_filt, pseudo_count = 1)
  normalized_data <- exprs(sce_filt)
  meta_data <- as.data.frame(colData(sce_filt))
  meta_data$sample_condition <- paste(meta_data$ID, meta_data$trt, sep="")
  sample_colname <- "sample_condition"
  IDs <- as.data.frame(meta_data)[, sample_colname]
  unique_ID_list <- as.list(unique(IDs))
  pseudobulk <- as.data.frame(pblapply(unique_ID_list, FUN = function(x){rowMeans(normalized_data[,IDs == x, drop = FALSE])}))
  colnames(pseudobulk) <- unique(IDs)
  rownames(pseudobulk) <- rownames(normalized_data)
  pseudobulk <- as.data.frame(pseudobulk)
  object_pseudobulkMeans <- pseudobulk
  return(object_pseudobulkMeans)
}