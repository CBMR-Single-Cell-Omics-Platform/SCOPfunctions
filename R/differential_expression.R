#' @param cells.1 Vector of cell names belonging to group 1
#' @param cells.2 Vector of cell names belonging to group 2
#' @param features Features to calculate fold change for.
#' If NULL, use all features
#' @importFrom Matrix rowSums
#' @rdname FoldChange
#' @concept differential_expression
#' @export
#' @method FoldChange default
#' . Stuart and Butler et al. Comprehensive Integration of Single-Cell Data. Cell (2019)
#'
#' . Copied from Seurat 3 source at https://github.com/satijalab/seurat/blob/master/R/differential_expression.R
FoldChange.default <- function(
  object,
  cells.1,
  cells.2,
  mean.fxn,
  fc.name,
  features = NULL,
  ...
) {
  features <- features %||% rownames(x = object)
  # Calculate percent expressed
  thresh.min <- 0
  pct.1 <- round(
    x = rowSums(x = object[rownames(object) %in% features, colnames(object) %in% cells.1, drop = FALSE] > thresh.min) /
      length(x = cells.1),
    digits = 3
  )
  pct.2 <- round(
    x = rowSums(x = object[rownames(object) %in% features, colnames(object) %in% cells.2, drop = FALSE] > thresh.min) /
      length(x = cells.2),
    digits = 3
  )
  # Calculate fold change
  data.1 <- mean.fxn(object[rownames(object) %in% features, colnames(object) %in% cells.1, drop = FALSE])
  data.2 <- mean.fxn(object[rownames(object) %in% features, colnames(object) %in% cells.2, drop = FALSE])
  fc <- (data.1 - data.2)
  fc.results <- as.data.frame(x = cbind(fc, pct.1, pct.2))
  colnames(fc.results) <- c(fc.name, "pct.1", "pct.2")
  return(fc.results)
}

#' MAST differential expression test with random effect for Seurat object
#'
#' Based on Seurat::FindMarkers MAST test but with a single random effects variable
#'
#' @param object Seurat >=3 object
#' @param random_effect.var Random effects to add in MAST model i.e. (1|random_effect.var) .
#' @param ident.1 Identity class to define markers for
#' @param ident.2 A second identity class for comparison. Leave as NULL to compare with all other cells
#' @param logfc.threshold Only return results with a DE exceeding threshold
#' @param base base for log when computing mean and for output, default exp(1). NB: Seurat has changed to log2 in V4.
#' @param group.by Regroup cells into a different identity class prior to performing differential expression
#' @param assay Assay to pull data from; defaults to default assay
#' @param slot Slot to pull data from; defaults to "data" as MAST expects log-normalized data
#' @param features Genes to test. Default is to use all genes
#' @param min.pct Only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expressed. Default is 0.1
#' @param max.cells.per.ident Down sample each identity class to a max number. Default is no downsampling. Not activated by default (set to Inf)
#' @param random.seed Random seed to use for down sampling
#' @param latent.vars Variables to test
#' @param n_cores How many cores to use (temporarily sets options(mc.cores=n_cores))
#' @param verbose whether to print stdout from MAST functions
#' @param p.adjust.method passed to p.adjust. Note that n is all the genes in the object
#' @param \dots Additional parameters (other than formula, sca, method, ebayes, strictConvergence) to pass to MAST::zlm
#' @return data.frame with column names p_val, avg_log[base]FC, pct.1, pct.2, p_val_adj (identical format to Seurat::FindMarkers)
#' @export
#'
#' @examples
#'
#' @references Zimmerman, K.D., Espeland, M.A. & Langefeld, C.D.
#' . A practical solution to pseudoreplication bias in single-cell studies.
#' . Nat Commun 12, 738 (2021). https://doi.org/10.1038/s41467-021-21038-1
#'
#' . McDavid A, Finak G, Yajima M (2020). MAST: Model-based Analysis of Single Cell Transcriptomics. R package version 1.16.0, https://github.com/RGLab/MAST/.
#'
#' . Stuart and Butler et al. Comprehensive Integration of Single-Cell Data. Cell (2019)
#'
#' . Seurat 3 source at https://github.com/satijalab/seurat/blob/master/R/differential_expression.R
MAST_DE_rand_effect_seurat = function(
  object,
  random_effect.var,
  ident.1,
  ident.2 = NULL,
  group.by = NULL,
  logfc.threshold = 0.25,
  base = exp(1),
  assay=NULL,
  slot="data",
  features = NULL,
  min.pct = 0.1,
  max.cells.per.ident = NULL,
  random.seed = 1,
  latent.vars = NULL,
  n_cores=NULL,
  verbose=TRUE,
  p.adjust.method="fdr",
  ...
) {

  if (!is.null(n_cores)) old_opt <- options("mc.cores"=n_cores)

  fc.name  <- if (base == exp(1)) "avg_logFC" else paste0("avg_log", base, "FC")
  #======== check inputs ========================================

  stopifnot(!is.null(random_effect.var))

  stopifnot(random_effect.var %in% colnames(object@meta.data))

  if (slot != "data") warning(paste0("MAST uses the logNormalised counts which are usually in the 'data' slot, but you are using ",slot))

  if (verbose) message("preparing data")
  if (verbose) message(paste0("using ", round(base,2), " as log base. Make sure that the data slot has been log-transformed using this base"))
  #======== resolve idents and get cells ===================

  logical.cells.1 = if (!is.null(group.by)) {object@meta.data[[group.by]] == ident.1} else {Idents(object) == ident.1}
  if (sum(logical.cells.1)==0) {
    stop("no cells found matching ident.1. Did you forget to set Idents(object) or use group.by?")
  }
  logical.cells.2 = if (is.null(ident.2)) {!logical.cells.1} else {if (!is.null(group.by)) {object@meta.data[[group.by]] == ident.2} else {Idents(object) == ident.2}}
  if (sum(logical.cells.2)==0) {
    stop("no cells found matching ident.2. Did you forget to set Idents(object) or use group.by?")
  }

  cells.1 = colnames(object)[logical.cells.1]
  cells.2 = colnames(object)[logical.cells.2]

  #======== features =================================

  features <- features %||% rownames(object)

  if (any(!features %in% rownames(object))) {
    warning(paste0(paste(features[!features %in% rownames(object)], collapse = ", "), " not found in data!"))
  }

  features = features[features %in% rownames(object)]
  vec_logical_features = rep(T, length(features))

  data = Seurat::GetAssayData(object = object, assay=assay, slot=slot)

  # compute average log fold change
  pseudocount.use = 1
  mean.fxn = function(x) {
    return(log(x = rowMeans(x = expm1(x = x)) + pseudocount.use, base = base))
  }
  # outputs a data.frame with columns fc.name, "pct.1", "pct.2"
  fc.results = FoldChange.default(
      object=data,
      cells.1=cells.1,
      cells.2=cells.2,
      mean.fxn=mean.fxn,
      fc.name=fc.name,
      features = features)

  if (!is.null(logfc.threshold)) {
    vec_logical_features = vec_logical_features & fc.results[[fc.name]]>logfc.threshold  # need to convert back to non-log space to take mean
  }
  if (sum(vec_logical_features)<2) stop(paste0("Fewer than two features have a log fold change above ", logfc.threshold))
  if (!is.null(min.pct)) {
    vec_logical_features = vec_logical_features & (fc.results$pct.1 > min.pct | fc.results$pct.2 > min.pct)
  }
  if (sum(vec_logical_features)<2) stop(paste0("Fewer than two features are expressed in ", min.pct, " of cells in any condition"))

  features <- features[vec_logical_features]

  data = data[rownames(data) %in% features, , drop=F]

  #======== down sample cells ==============================

  idx.cells.1 = which(logical.cells.1)
  idx.cells.2 = which(logical.cells.2)

  if (!is.null(max.cells.per.ident)) {
    if (sum(logical.cells.1) > max.cells.per.ident) {
      set.seed(randomSeed)
      idx.cells.1 = sample(x = idx.cells.1, size = max.cells.per.ident, replace = F)
    }
    if (sum(!logical.cells.1) > max.cells.per.ident) {
      set.seed(randomSeed)
      idx.cells.2 = sample(x = idx.cells.2, size = max.cells.per.ident, replace = F)
    }
  }

  cells.1 = colnames(data)[idx.cells.1]
  cells.2 = colnames(data)[idx.cells.2]
  data = data[,c(idx.cells.1,idx.cells.2), drop=F]

  #======== filter out all zero features (if not already removed) ==

  data = data[apply(X = data, MARGIN=1, FUN = sum)>0,]

  #======== prep column (cell/sample) data =================

  df_coldata = data.frame(
    row.names = c(cells.1, cells.2),
    "wellKey" = c(cells.1, cells.2)
    ) # wellKey is hardcoded feature name in MAST

  # DE group factor
  df_coldata[cells.1, "group"] <- "Group1"
  df_coldata[cells.2, "group"] <- "Group2"
  df_coldata[, "group"] <- factor(x = df_coldata[, "group"])

  # latent vars
  if (!is.null(latent.vars)) {
    for (latent.var in latent.vars) {
      df_coldata[[latent.var]] <- object@meta.data[c(idx.cells.1,idx.cells.2), latent.var]
    }
  }

  # random var
  df_coldata[[random_effect.var]] <- as.factor(object@meta.data[c(idx.cells.1,idx.cells.2),
                                                                random_effect.var])

  # check whether random var levels overlap entirely with test condition
  for (re_lvl in levels(df_coldata[[random_effect.var]])) {
    if (all(df_coldata[[random_effect.var]]==re_lvl & 1:ncol(data) %in% idx.cells.1)) {
      stop(paste0("random_effect.var level ", re_lvl, " overlaps entirely with ident.1"))
    }
    else if (all(df_coldata[[random_effect.var]]==re_lvl & 1:ncol(data) %in% idx.cells.2)) {
      stop(paste0("random_effect.var level ", re_lvl, " overlaps entirely with ident.2"))
    }
  }

  #======== prep  feature data =========================
  df_feature = data.frame("primerid"=rownames(data)) # primerid is hardcoded feature name in MAST

  #======== prep SingleCellAssay object ============================

  expr = quote(MAST::FromMatrix(exprsArray=as.matrix(data),
                                 check_sanity = TRUE,
                                 cData=df_coldata,
                                 fData=df_feature))
  sca <- if (verbose) eval(expr) else suppressMessages(eval(expr))

  # set group as a factor and make group1 the reference
  #cond <- factor(x = SummarizedExperiment::colData(sca)$group)
  # in Seurat source this is Group1!?
  # see https://github.com/satijalab/seurat/blob/master/R/differential_expression.R
  SummarizedExperiment::colData(sca)$group <- relevel(SummarizedExperiment::colData(sca)$group, ref="Group2")
  # cond <- relevel(x = cond, ref = "Group2")
  # SummarizedExperiment::colData(sca)$condition <- cond

  # add the centered number of expressed genes to the model
  # cdr2 <- colSums(SummarizedExperiment::assay(sca)>0)
  # SummarizedExperiment::colData(sca)$cngeneson <- scale(cdr2)

  #======== run MAST test ======================

  str_plus = if (!is.null(latent.vars)) " + " else ""
  fmla <- as.formula(
    object = paste0(" ~ group", str_plus, paste(latent.vars, collapse = "+"), " + (1|", random_effect.var, ")")
  )

  # fit model parameters
  expr = quote(MAST::zlm(formula = fmla,
                   sca = sca,
                   method = 'glmer',
                   ebayes = F,
                   strictConvergence = FALSE,
                   ...))

  zlmCond <- if (verbose) eval(expr) else suppressMessages(eval(expr))

  if (verbose) {
    print(zlmCond)
    message("Compute likelihoods and p-values")
  }
  # call a likelihood ratio test on the fitted object
  expr = quote(summary(object = zlmCond, doLRT = 'groupGroup1'))
  summaryCond <- if (verbose) eval(expr) else suppressMessages(eval(expr))

  # print(summaryCond,n=4)
  # Fitted zlm with top 4 genes per contrast:
  #   ( log fold change Z-score )
  # primerid           conditionGroup1 cngeneson
  # Akr1c20              -4.1            13.2*
  # Bclaf3               18.5*            4.8
  # Cmss1                -0.3            -9.9*
  # Cux2                -16.5*            4.6
  # ENSMUSG00000115022  -19.5*            7.6
  # Pot1b                13.5*            5.3
  # Siah1a                4.6            17.6*
  # Sipa1l2               1.8            10.5*
  #
  #
  # summaryCond$datatable has colnames
  # "primerid" feature name
  # "component" if no latent.vars, "C" (continuous) "D" (discrete)  "H" (hurdle) "S" (RE?) "logFC"
  # "contrast" if no latent.vars,  conditionGroup1 (Intercept) cngeneson
  # "ci.hi" confidence interval (high)
  # "ci.lo" confidence interval (low)
  # "coef" coefficient estimate
  # "z" z-scores

  # uses
  # https://www.bioconductor.org/packages/release/bioc/vignettes/MAST/inst/doc/MAITAnalysis.html#differential-expression-using-a-hurdle-model
  # and
  # https://github.com/kdzimm/PseudoreplicationPaper/blob/master/Type_1_Error/Type%201%20-%20MAST%20RE.Rmd

  de.results <- data.frame(
    "p_val" = summaryCond$datatable[contrast=='groupGroup1' & component=='H', `Pr(>Chisq)`],
    fc.results[vec_logical_features,])  #setDF(summaryCond$datatable[contrast=='groupGroup1' & component=='logFC', .(coef)])
  #colnames(df_out)[1] = fc.name
  de.results$p_val_adj = p.adjust(de.results$p_val, method=p.adjust.method, n=length(vec_logical_features))
  #rownames(df_out) = summaryCond$datatable[contrast=='groupGroup1' & component=='H', primerid]

  de.results = de.results[order(de.results$p_val, -de.results[[fc.name]]),]
  # reset mc.cores option
  if (!is.null(n_cores)) options(mc.cores=old_opt$mc.cores)

  return(de.results)
}
