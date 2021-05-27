#' Make an area plot of proportions of cells assigned to singlets, doublets, and negatives
#' as a function of the positive.quantile provided to Seurat::HTODemux
#'
#' @param seurat_obj A seurat object
#' @param vec_range_quantile a numeric vector of quantile values (0-1)
#' @param n_cores_max max number of cores to use in multicore processing
#'
#' @return a ggplot2 object
#' @export
#'
prep_HTO_q_area_plot <- function(
  seurat_obj,
  vec_range_quantile=seq(0.8,0.99,0.01),
  n_cores_max=Inf
) {

  # require("furrr")
  # require("Seurat")
  # require("ggplot2")
  # require("RColorBrewer")

  # https://www.tidyverse.org/blog/2020/04/self-cleaning-test-fixtures/#the-onexit-pattern
  # set new value and capture old in op
  op <- options(
    "mc.cores"=min(n_cores_max, length(vec_range_quantile)),
    "future.globals.maxSize" = object.size(seurat_obj)*1.1)
  on.exit(options(op), add = T)

  current_assay = Seurat::DefaultAssay(seurat_obj)
  on.exit(expr= Seurat::DefaultAssay(seurat_obj) <- current_assay)
  Seurat::DefaultAssay(seurat_obj)  <- "HTO"

  singlet_labels = unique(rownames(seurat_obj))
  doublet_labels = purrr::cross2(unique(rownames(seurat_obj)),unique(rownames(seurat_obj))) %>%
    lapply(X = ., FUN=function(vec) paste(vec, collapse="_")) %>% unlist(use.names=F)
  doublet_labels = doublet_labels[!doublet_labels %in% stringr::str_c(singlet_labels,singlet_labels, sep="_")]

  strategy = if (getOption("mc.cores")>1) "multicore" else "sequential"
  future::plan(strategy=strategy)

  # iterate over quantile cutoffs
  # compute the ratio of singlets to doublets
  list_vec_classif = furrr::future_map(
    .x = vec_range_quantile,
    .f = function(q) {
      seurat_obj_tmp = Seurat::HTODemux(
        seurat_obj,
        assay = "HTO",
        positive.quantile = q,
        seed=TRUE)

    vec_out = sapply(c(singlet_labels, doublet_labels, "Negative"), function(cl) {
      sum(seurat_obj_tmp$HTO_classification == cl)
    })
    # vec_out["Singlet_doublet_ratio"] = sum(vec_out[names(vec_out) %in% singlet_labels])/
    #                                sum(vec_out[names(vec_out) %in% doublet_labels])
    vec_out
    })

  # prepare results for plotting
  tbbl = sapply(list_vec_classif, FUN= I ) %>% t %>% data.frame %>% tibble
  colnames(tbbl) = c(singlet_labels, doublet_labels, "Negative")#, "Singlet_doublet_ratio")
  tbbl$quantile = vec_range_quantile

  # prepare plot
  cols = c(singlet_labels, doublet_labels, "Negative")[c(singlet_labels, doublet_labels, "Negative") %in% colnames(tbbl)]
  tbbl_long <- tidyr::pivot_longer(tbbl,cols = cols, names_to = "classification")
  tbbl_long$classification = factor(tbbl_long$classification, levels = c(singlet_labels, doublet_labels, "Negative"), ordered = T)
  vec_colors_singlet = colorRampPalette(colors = RColorBrewer::brewer.pal(n=min(11,length(singlet_labels)), name = "Spectral"))(length(singlet_labels))
  vec_colors_doublet = colorRampPalette(colors = RColorBrewer::brewer.pal(n=min(11,length(doublet_labels)), name = "BrBG"))(length(doublet_labels))
  names(vec_colors_singlet) = singlet_labels
  names(vec_colors_doublet) = doublet_labels
  vec_colors = c(vec_colors_singlet,vec_colors_doublet)
  vec_colors["Negative"] = "grey90"

  p <- ggplot(tbbl_long, mapping=aes(x=quantile, y=value, fill=classification)) +
    geom_area(stat = "identity",
              position = "stack") +
    scale_fill_manual(values=vec_colors)

  # find the best ratio of singlets to doublets
  # idx_max = which.max(tbbl$Singlet_doublet_ratio)
  # q_max = tbbl$quantile[idx_max]

  # p + geom_vline(xintercept=q_max)

  return(p)
}


#' A subroutine of prep_intrahash_doub to compute number of doublet neighbours
#'
#' @param x an SingleCellExperiment object
#'
#' @return a logical vector of prediction scores
#'
.dub_cutoff <- function(x) {

  # calculate number of neighbors at each proportion that are doublets
  data.frame("prop"=x$proportion_dub_neighbors) %>%
    group_by(prop) %>%
    summarize(n=n()) %>%
    mutate(pct = n/sum(n)) -> data
  # find point at which we gain very few doublets as proportion increases
  cut <- data$prop[PCAtools::findElbowPoint(variance = sort(data$n, decreasing = T))+1]
  vec <- if_else(x$proportion_dub_neighbors <= cut, F, T)
  return(vec)
}


#' Estimate intra-hashtag doublets by finding cells whose nearest neighbours are doublets
#'
#' @param seurat_obj a Seurat object
#' @param assay Seurat object assay to use
#' @param npcs  number of principal components to use
#' @param randomSeed random seed for random number generator
#'
#' @return Seurat object with a predicted_dub_intra logical column added to metadata
#' @export
#'
prep_intrahash_doub = function(
  seurat_obj,
  assay = "RNA",
  npcs=20,
  randomSeed = 12345
) {

  # require("Seurat")
  # require("scuttle")
  # require("scran")
  # require("scater")
  # require("scDblFinder")

  seurat_obj_cp = seurat_obj
  # set assay
  # current_assay = Seurat::DefaultAssay(seurat_cp)
  # on.exit(expr= Seurat::DefaultAssay(seurat_obj) <- current_assay)
  Seurat::DefaultAssay(seurat_obj_cp)  <- assay

  # this apparently only keeps the current assay
  sce.full = Seurat::as.SingleCellExperiment(seurat_obj_cp)
  rm(seurat_obj_cp)
  sce.full$doublet <- if_else(sce.full$HTO_classification.global == "Doublet", true = T, false = F)
  sce.full <- scuttle::logNormCounts(sce.full)
  dec.hash <- scran::modelGeneVar(sce.full)
  top.hash <- scran::getTopHVGs(dec.hash, n=1000)
  set.seed(randomSeed)
  sce.full <- scater::runPCA(
    sce.full,
    subset_row=top.hash,
    ncomponents=npcs)

  # Recovering the intra-sample doublets:
  hashed.doublets <- scDblFinder::recoverDoublets(
    sce.full,
    use.dimred="PCA",
    doublets=sce.full$doublet,
    samples=table(sce.full$HTO_maxID))


  sce.full$proportion_dub_neighbors <- hashed.doublets$proportion
  sce.full$predicted_dub_std <- hashed.doublets$predicted
  sce.full$predicted_dub_intra <- .dub_cutoff(sce.full)

  # add metadata to original seurat object
  seurat_obj$predicted_dub_intra <- sce.full$predicted_dub_intra

  return(seurat_obj)

}



#' Perform QC on RNA, classifying cells as keepers or outliers to be filtered
#'
#' @param seurat_obj Seurat object
#' @param assay assay to use
#'
#' @return Seurat object with logical qc_discard column added to metadata
#' @export
#'
prep_qc_rna <- function(
  seurat_obj,
  assay = "RNA") {

  # require("Seurat")
  # require("scuttle")
  # require("dplyr")

  # set assay
  current_assay = Seurat::DefaultAssay(seurat_obj)
  on.exit(expr= Seurat::DefaultAssay(seurat_obj) <- current_assay)
  Seurat::DefaultAssay(seurat_obj)  <- assay

  sce.full = Seurat::as.SingleCellExperiment(seurat_obj)

  # qc filtering
  sce.full <- scuttle::addPerCellQC(sce.full)
  reason <- data.frame(
    "qc.umi.low" = scuttle::isOutlier(sce.full$sum, log=TRUE, type="lower"),
    "qc.umi.high" = scuttle::isOutlier(sce.full$sum, type="higher"),
    "qc.nexprs.low" = scuttle::isOutlier(sce.full$detected, log=TRUE, type="lower"),
    "qc.nexprs.high" = scuttle::isOutlier(sce.full$detected, type="higher")
  ) %>%
    mutate(reason = pmap_int(.,~any(.==T)),
           reason = if_else(reason == 1, true=T, false=F)) %>%
    dplyr::pull(reason)

  # annotate object to discard
  sce.full$qc_discard <- reason

  return(Seurat::as.Seurat(sce.full))
}
