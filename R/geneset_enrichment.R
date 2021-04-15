#' Compute Variance Inflation Factors for gene sets
#'
#' Utility function to compute Variance Inflation Factors for gene sets using an expression dataset
#'        in order to account for inter-gene correlation in gene sets derived from expression data
#'
#' @param datExpr gene * cell expression matrix with row and column names
#' @param list_genesets list of genesets, using same gene names as datExpr
#'
#' @return list of 2-vectors, one per geneset. First vector element is VIF, second is mean Pearson's rho correlation.
#' @export
#'
#' @references modified from https://rdrr.io/bioc/qusage/src/R/qusage.R
#' . drawing on https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3458527/
#' @examples list_vec_VIF <- geneset_VIF(datExpr=mat_counts, list_genesets=mylist)
geneset_VIF <- function(datExpr, list_genesets) {
  list_vec_vif = lapply(names(list_genesets), function(genesetname) {
    vec_logicalgenes <-rownames(datExpr) %in% list_genesets[[genesetname]]
    if (sum(vec_logicalgenes) < 2) {
      warning("GeneSet '", genesetname, "' contains one or zero overlapping genes. NAs produced.")
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
