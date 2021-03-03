
f_gene_map <- function(dataIn,
                       colGene = NULL,
                       df_mapping,
                       from="hgnc",
                       to="ensembl",
                       replace = F,
                       na.rm = T) {
  #' @usage map genes in a dataframe column, vector, or list of vectors between naming schemes in df_mapping
  #' @param dataIn data.frame or data.table with a column or rownames containing genes to remap,
  #' matrix with gene rownames, a list of vectors, or a vector, either named numeric or character
  #' If a list, if the vectors are numeric, the vector names are assumed to be genes, if the vectors are character
  #' the vector values are assumed to be genes
  #' @param colGene if dataIn is a dataframe, the name of the gene column; in this case NULL implies rownames
  #' @param df_mapping data.frame or matrix with columns corresponding to 'from' and 'to' arguments (using grep partial matching)
  #' @param from df_mapping colnames
  #' @param to df_mapping colnames
  #' @param replace boolean; if dataIn is a data.frame, TRUE replaces original gene names, FALSE adds a new column to the data.frame
  #' @param na.rm boolean; remove genes that fail to map or leave them as NAs?; defaults to TRUE
  #' @value an object of the same format as dataIn with new gene names
  
  #' @example :
  # df_test <- gene_map(dataIn=load_obj("/projects/jonatan/tmp-mousebrain/tables/mousebrain_Vascular_ClusterName_1_PER3_kIM.csv"),
  #                     colGene ="genes",
  #                     df_mapping=fread("/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/Mus_musculus.GRCm38.90.gene_name_version2ensembl.txt.gz"),
  #                     from="ensembl",
  #                     to="gene_name_optimal",
  #                     replace=T,
  #                     na.rm=T)
  # ortholog mapping:
  # /projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/gene_annotation.hsapiens_mmusculus_unique_orthologs.GRCh37.ens_v91.txt.gz
  
  
  stopifnot(any(class(df_mapping)%in%c("data.frame", "data.table")))
  stopifnot(length(from)>0 & length(to)>0)
  
  fromMapCol <- if (from %in% colnames(df_mapping)) from else grep(pattern=from, x = colnames(df_mapping), ignore.case=T, value = T)
  toMapCol <- if (to %in% colnames(df_mapping)) to else grep(pattern=to, x = colnames(df_mapping), ignore.case=T, value = T)
  
  if (length(fromMapCol)==0) stop(paste0(from, " not found in df_mapping column names"))
  if (length(toMapCol)==0) stop(paste0(to, " not found in df_mapping column names"))
  
  if (!is.null(dim(dataIn))) {
    
    if(is.null(colGene) & replace==T) message("Duplicate genes will be averaged and merged to keep row.names unique")
    
    genes_from <- if (is.null(colGene)) rownames(dataIn) else dataIn[[colGene]]
    idx_match <- match(toupper(gsub("-|_", ".", genes_from)), toupper(gsub("-|_", ".", df_mapping[[fromMapCol]])))
    genes_to <- df_mapping[[toMapCol]][idx_match]
    
    if (replace) { # remove NAs
      if (is.null(colGene)) {
        # average identical gene names to ensure unique row names
        dataIn_aggr <- aggregate(dataIn, by= list(genes_to), FUN=mean, na.rm=T)
        rownames(dataIn_aggr) <- dataIn_aggr[["Group.1"]]
        dataIn <- within(dataIn_aggr, rm("Group.1"))
      } else {
        if (na.rm) {
          dataIn <- dataIn[!is.na(genes_to),]
          genes_to <- genes_to[!is.na(genes_to)]
        }
        dataIn[[colGene]] <- genes_to
        colnames(dataIn)[which(colnames(dataIn)==colGene)] <- to
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
