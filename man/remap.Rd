% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{remap}
\alias{remap}
\title{Map  dataframe column, vector, or list of vectors}
\usage{
remap(
  dataIn,
  colRemap = NULL,
  df_mapping,
  from = "hgnc",
  to = "ensembl",
  replace = F,
  na.rm = T
)
}
\arguments{
\item{dataIn}{data.frame or data.table with a column or rownames containing genes to remap,}

\item{colRemap}{if dataIn is a dataframe, the name of the column to remap; in this case NULL implies rownames}

\item{df_mapping}{data.frame or coercible with columns corresponding to 'from' and 'to' arguments (using grep partial matching)}

\item{from}{df_mapping colname to map from}

\item{to}{df_mapping colname to map to}

\item{replace}{boolean; if dataIn is a data.frame, TRUE replaces original, FALSE adds a new column to the data.frame}

\item{na.rm}{boolean; remove rows that fail to map or leave them as NAs?; defaults to TRUE}
}
\value{
an object of the same format as dataIn with new gene names
}
\description{
Use a user-provided look-up table, df_mapping, to map a vector. Handle NAs and collapse duplicates.
}
