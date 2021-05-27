#' Convert a large sparse matrix into a dense matrix without errors
#'
#' Avoid the following error
#' . Error in asMethod(object) :
#' . Cholmod error 'problem too large' at file ../Core/cholmod_dense.c,
#' . by slicing the matrix into submatrices, converting and cbinding them
#' . Increases number of slices until they succeed.
#'
#' @param sparseMat a big sparse matrix of a type coercible to dense Matrix::Matrix
#' @param n_slices_init initial number of slices. Default value 1, i.e. whole matrix
#' @param verbose print progress
#'
#' @return a dense matrix
#' @export
#'
utils_big_as.matrix <- function(
  sparseMat,
  n_slices_init=1,
  verbose=T
  ) {

  n_slices <- n_slices_init-1
  while (TRUE) {
    list_densemat = list()
    n_slices = n_slices+1
    if (verbose) message(paste0("n_slices=",n_slices))
    idx_to = 0
    for (slice in 1:n_slices) {
      if (verbose) message(paste0("converting slice ",slice,"/",n_slices))
      idx_from <- idx_to+1
      idx_to <- if (slice<n_slices) as.integer(ncol(sparseMat)*slice/n_slices) else ncol(sparseMat)
      if (verbose) message(paste0("columns ", idx_from,":", idx_to))
      densemat_sub = try(
        expr = {
          as.matrix(sparseMat[,idx_from:idx_to])
        }, silent = if (verbose) FALSE else TRUE)
      if ("try-error" %in% class(densemat_sub)) {
        break # exit to while loop
      } else {
        list_densemat[[slice]] = densemat_sub
      }
    }
    if (length(list_densemat)==n_slices) break # exit while loop
  }
  if (verbose) message("cbind dense submatrices")
  densemat <- Reduce(f=cbind, x=list_densemat)
  return(densemat)
}




#' Map  dataframe column, vector, or list of vectors
#'
#' Use a user-provided look-up table, df_mapping, to map a vector. Handle NAs and collapse duplicates.
#'
#' @param dataIn data.frame or data.table with a column or rownames containing genes to remap,
#' @param colRemap if dataIn is a dataframe, the name of the column to remap; in this case NULL implies rownames
#' @param df_mapping data.frame or coercible with columns corresponding to 'from' and 'to' arguments (using grep partial matching)
#' @param from df_mapping colname to map from
#' @param to df_mapping colname to map to
#' @param replace boolean; if dataIn is a data.frame, TRUE replaces original, FALSE adds a new column to the data.frame
#' @param na.rm boolean; remove rows that fail to map or leave them as NAs?; defaults to TRUE
#'
#' @return an object of the same format as dataIn with new gene names
#' @export
#'
remap <- function(dataIn,
                  colRemap = NULL,
                  df_mapping,
                  from="hgnc",
                  to="ensembl",
                  replace = F,
                  na.rm = T) {


  stopifnot(any(class(df_mapping)%in%c("data.frame", "data.table")))
  stopifnot(length(from)>0 & length(to)>0)

  fromMapCol <- if (from %in% colnames(df_mapping)) from else grep(pattern=from, x = colnames(df_mapping), ignore.case=T, value = T)
  toMapCol <- if (to %in% colnames(df_mapping)) to else grep(pattern=to, x = colnames(df_mapping), ignore.case=T, value = T)

  if (length(fromMapCol)==0) stop(paste0(from, " not found in df_mapping column names"))
  if (length(toMapCol)==0) stop(paste0(to, " not found in df_mapping column names"))

  if (!is.null(dim(dataIn))) {

    if(is.null(colRemap) & replace==T) message("Duplicate genes will be averaged and merged to keep row.names unique")

    genes_from <- if (is.null(colRemap)) rownames(dataIn) else dataIn[[colRemap]]
    idx_match <- match(toupper(gsub("-|_", ".", genes_from)), toupper(gsub("-|_", ".", df_mapping[[fromMapCol]])))
    genes_to <- df_mapping[[toMapCol]][idx_match]

    if (replace) { # remove NAs
      if (is.null(colRemap)) {
        # average identical gene names to ensure unique row names
        dataIn_aggr <- aggregate(dataIn, by= list(genes_to), FUN=mean, na.rm=T)
        rownames(dataIn_aggr) <- dataIn_aggr[["Group.1"]]
        dataIn <- within(dataIn_aggr, rm("Group.1"))
      } else {
        if (na.rm) {
          dataIn <- dataIn[!is.na(genes_to),]
          genes_to <- genes_to[!is.na(genes_to)]
        }
        dataIn[[colRemap]] <- genes_to
        colnames(dataIn)[which(colnames(dataIn)==colRemap)] <- to
      }
    }  else {
      if (na.rm) {
        vec_logicalNA <- is.na(genes_to)
        dataIn <- dataIn[!vec_logicalNA,,drop=F]
        dataIn[[to]] <- genes_to[!vec_logicalNA]
      } else {
        dataIn[[to]] <- genes_to
      }
    }
  } else if (class(dataIn)=="list") {
    dataIn <- lapply(dataIn, function(eachVec) {
      oldNames <- if (class(eachVec)== "numeric") {
        names(eachVec)
      } else if (class(eachVec)=="character") {
        eachVec
      }
      newNames <- df_mapping[[toMapCol]][match(toupper(gsub("-|_", ".", oldNames)), toupper(gsub("-|_", ".", df_mapping[[fromMapCol]])))]
      if (na.rm) {
        vec_logicalNA <- is.na(newNames)
        eachVec <- eachVec[!vec_logicalNA]
        newNames <- newNames[!vec_logicalNA]
      }
      if (class(eachVec)=="numeric") names(eachVec) <- newNames else eachVec <- newNames
      return(eachVec)
    })
  } else if (class(dataIn)=="numeric") {
    newNames <- df_mapping[[toMapCol]][match(toupper(gsub("-|_", ".", names(dataIn))), toupper(gsub("-|_", ".", df_mapping[[fromMapCol]])))]
    if (na.rm) {
      vec_logicalNA <- is.na(newNames)
      newNames <- newNames[!vec_logicalNA]
      dataIn <- dataIn[!vec_logicalNA]
    }
    names(dataIn) <- newNames
  } else if (class(dataIn)=="character") {
    dataIn <- df_mapping[[toMapCol]][match(toupper(gsub("-|_", ".", dataIn)), toupper(gsub("-|_", ".", df_mapping[[fromMapCol]])))]
    if (na.rm) dataIn <- dataIn[!is.na(dataIn)]
  }
  return(dataIn)
}
