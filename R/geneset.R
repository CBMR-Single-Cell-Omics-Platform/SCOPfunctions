#' Generate activity scores (embeddings) for a feature set in a sample matrix
#'
#' @param mat_datExpr a m feature * n sample matrix
#' @param vec_geneWeights a numeric vector of weights with feature names. Weights should be pre-normalized
#' @param min_feats_present a minimum number of features shared between mat_datExpr and vec_geneWeights
#'
#' @return returns a n vector of gene weights, scaled to 1 standard deviation, but not centred.
#' @export
#'
geneset_embed <- function(
  mat_datExpr,
  vec_geneWeights,
  min_feats_present = 5) {
#
#   require("magrittr")
#   require("data.table")

  vec_logic_all_0 = apply(mat_datExpr, MARGIN=1, FUN=function(eachRow) all(eachRow==0))

  mat_datExpr = mat_datExpr[rownames(mat_datExpr) %in% names(vec_geneWeights) & !vec_logic_all_0,]

  mat_datExpr %>% t %>% scale(., center=FALSE, scale=apply(.,2,sd,na.rm=T)) %>% t  -> mat_datExpr_scaled

  # find corresponding rows of the expression matrix
  vec_idxRow <- match(names(vec_geneWeights), rownames(mat_datExpr_scaled))

  # remove genes that aren't in datExpr
  vec_geneWeights <- vec_geneWeights[!is.na(vec_idxRow)]
  vec_idxRow <- vec_idxRow[!is.na(vec_idxRow)]

  if (length(vec_idxRow) < min_feats_present) {
    stop(paste0("Only ", length(vec_idxRow), " of the features are present in datExpr"))
  }

  # compute weighted sum of normalized and scaled expression
  return(as.numeric(vec_geneWeights %*% mat_datExpr_scaled[vec_idxRow,]))
}


#' Generate activity scores (embeddings) for a list of feature sets and add to a Seurat object
#'
#' A variant of geneset_embed that works with a list of feature sets and takes and returns a Seurat object
#'
#' @param seurat_obj A Seurat object
#' @param list_vec_geneWeights a list of numeric vectors of weights with feature names. Weights should be pre-normalized
#' @param slot Seurat object slot to use
#' @param assay Seurat object assay to use
#' @param min_feats_present a minimum number of features shared between mat_datExpr and vec_geneWeights
#' @param n_cores_max max number of cores to use in multicore processing. To use sequential processing, set to 1
#'
#' @return a Seurat object with feature embeddings added to metadata
#' @export
#'
geneset_embed_list_seurat <- function(
  seurat_obj,
  list_vec_geneWeights,
  slot="scale.data",
  assay="SCT",
  min_feats_present = 5,
  n_cores_max = Inf) {

  # require("magrittr")
  # require("data.table")
  # require("furrr")
  # require("Seurat")

  op <- options("mc.cores"=min(length(list_vec_geneWeights), n_cores_max))
  on.exit(options(op), add = T)

  #vec_genes_union = Reduce(f=function(a,b) unique(c(names(a),names(b))), x=list_vec_geneWeights)
  vec_genes_union = c()
  for (vec in list_vec_geneWeights) {
    vec_genes_union = unique(c(names(vec),vec_genes_union))
  }
  # scale matrix to unit standard deviation
  mat_datExpr = Seurat::GetAssayData(seurat_obj, assay=assay, slot=slot)
  # only keep relevant genes
  idx_match = match(vec_genes_union, rownames(mat_datExpr))
  idx_match = idx_match[!is.na(idx_match)]
  mat_datExpr = mat_datExpr[idx_match,]

  if (!grepl("scale", slot)) mat_datExpr %>% t %>% scale(., center=FALSE, scale=apply(.,2,sd,na.rm=T)) %>% t -> mat_datExpr

  strategy = if (getOption("mc.cores")>1) "multicore" else "sequential"
  future::plan(strategy=strategy)

  list_vec_embed = furrr::future_map(
    .x=list_vec_geneWeights,
    .f=function(vec_geneWeights) {
      # find corresponding rows of the expression matrix
      vec_idxRow <- match(names(vec_geneWeights), rownames(mat_datExpr))

      # remove genes that aren't in datExpr
      vec_geneWeights <- vec_geneWeights[!is.na(vec_idxRow)]
      vec_idxRow <- vec_idxRow[!is.na(vec_idxRow)]


      if (length(vec_idxRow) < min_feats_present) {
        stop(paste0("only ", length(vec_idxRow), " of the features are present in datExpr"))
      }

      # compute weighted sum of normalized and scaled expression
      vec_embed = vec_geneWeights %*% mat_datExpr[vec_idxRow,] %>% as.numeric
      vec_embed = scale(vec_embed, center=F, scale=sd(vec_embed))
      return(vec_embed)
    })

  names(list_vec_embed) = names(list_vec_geneWeights)
  df_embed = data.frame(list_vec_embed,
                        row.names = colnames(seurat_obj))

  seurat_obj <- AddMetaData(seurat_obj, metadata=df_embed)

  return(seurat_obj)
}


#' Compute Variance Inflation Factors for gene sets
#'
#' Utility function to compute Variance Inflation Factors for gene sets using an expression dataset
#'        in order to account for inter-gene correlation in gene sets derived from expression data
#'
#' @param datExpr gene * cell expression matrix with row and column names
#' @param list_genesets list of genesets, using same gene names as datExpr
#' @param min_feats_present a minimum number of features shared between mat_datExpr and vec_geneWeights
#'
#' @return list of 2-vectors, one per geneset. First vector element is VIF, second is mean Pearson's rho correlation.
#' @export
#'
#' @references modified from https://rdrr.io/bioc/qusage/src/R/qusage.R
#' . drawing on https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3458527/
#'
geneset_VIF <- function(
  datExpr,
  list_genesets,
  min_feats_present = 5) {

  list_vec_vif = lapply(names(list_genesets), function(genesetname) {
    vec_logicalgenes <-rownames(datExpr) %in% list_genesets[[genesetname]]
    if (sum(vec_logicalgenes) < min_feats_present) {
      warning(paste0("GeneSet '", genesetname, "' contains fewer than", min_feats_present, " overlapping genes. NAs produced."))
      return(c("vif"=NA, "mean.cor"=NA))
    }
    cor.mat <- cor(t(datExpr[vec_logicalgenes, ]), use = "pairwise.complete.obs")
    cor.mat[is.na(cor.mat)] <- 0
    (cor.mat - Diagonal(x = diag(cor.mat))) %>% mean -> mean.cor
    vif <- 1+(sum(vec_logicalgenes)-1)*mean.cor
    return(c("vif"=vif, "mean.cor"=mean.cor))
  })
  names(list_vec_vif) <- names(list_genesets)
  return(list_vec_vif)
}
