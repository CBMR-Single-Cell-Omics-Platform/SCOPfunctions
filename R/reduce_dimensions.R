
#' Estimate dataset dimensionality
#'
#' @param x : Seurat or SingleCellExperiment object
#' @param max_dim : max number of dimensions to try
#' @param seed.use : random seed
#'
#' @return integer number of dimensions, <= max_dim
#' @export
#'
#' @examples ##estimate_dimensions(seurat_obj)
estimate_dimensions <- function(x,
                       max_dim=75,
                       seed.use=12345) {

 if(is(x, "SingleCellExperiment")) {
   x <- as.Seurat(x)
 }

 x %>% RunPCA(., npcs=max_dim, seed.use=seed.use, verbose = F) -> x

 dims <- round(as.numeric(intrinsicDimension::maxLikGlobalDimEst(data = x@reductions$pca[, 1:max_dim], k = 20)))

 return(dims)
}


