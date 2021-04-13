#' Detach all packages from current session
#'
#' @return None
#' @export
#' @references https://stackoverflow.com/questions/7505547/detach-all-packages-while-working-in-r
#' @examples utils_detachAllPkgs()
utils_detachAllPkgs <- function() {
  while (!is.null(names(sessionInfo()$otherPkgs))) lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
}


#' Check if a file matches given hash, else download it and check again
#'
#' @param url passed to download.file()
#' @param destfile passed to download.file().
#' @param correct_checksum file md5 checksum
#'
#' @return
#' @export
#'
#' @examples utils_dlFile(url=https:/.mysite.com/files/test.csv, destfile="./data", correct_checksum="2c3294272a968ffc91c0c2af2c296988")
utils_dlFile <- function(url, destfile, correct_checksum) {
  if(!file.exists(destfile)){
    print("Downloading file")
    download.file(url, destfile, quiet = FALSE)
    print("Verifying file integrity...")
    checksum = try(md5sum(destfile))
    if(checksum != correct_checksum){
      warning("downloaded file checksum does not match correct_checksum!")
    }
  } else{
    print("Verifying file integrity...")
    checksum = md5sum(destfile)
    if(checksum != correct_checksum){
      print("Existing file looks corrupted or is out of date, downloading again.")
      try(download.file(url, destfile, quiet = FALSE))
      print("Verifying file integrity...")
      checksum = try(md5sum(destfile))
      if(checksum != correct_checksum){
        warning("downloaded file checksum does not match correct_checksum!")
      }
    } else{
      print("Latest file already exists.")
    }
  }
}


#' Show the memory size of the top n objects
#'
#' @param n how many objects to list
#'
#' @return
#' @export
#'
#' @examples utils_printObjectSizes(10)
utils_printObjectSizes <- function(n) {
  sapply(ls(), function(x) object.size(eval(parse(text=x)))) %>% sort(., decreasing=T) %>% head(., n)
}

#' Detect cores with RAM constraints
#'
#' @description Given a ceiling on RAM to use and predicted extra RAM needed in
#' . each worker in addition to current usage, return safe number of cores
#' @param Gb_max ceiling on session memory usage in Gb, assuming that each worker duplicates the session memory
#' @param additional_Gb max additional memory requirement for new (temporary) objects created within a parallel session
#'
#' @return max number of cores (integer)
#' @export
#'
#' @examples n_cores_use = utils_detectCoresRAM(Gb_max=250,additional_Gb=1)
utils_detectCoresRAM <- function(Gb_max=250,
                                additional_Gb=1) {
  obj_size_Gb <- as.numeric(sum(sapply(ls(envir = .GlobalEnv), function(x) object.size(x=eval(parse(text=x))))) / 1024^3)
  max(1, min(parallel::detectCores(), Gb_max %/% (obj_size_Gb + additional_Gb))-1)
}


#' Wrapper to make a function more verbose
#'
#' @param fnc function to execute
#' @param \dots named arguments to function
#'
#' @return value returned by fnc
#' @export
#'
#' @examples utils_verboseFnc(rowSums, x=Matrix::Matrix(data=rnorm(9), nrow = r), na.rm = TRUE)
utils_verboseFnc <- function(fnc,...) {
  # TODO: fix output to outfile. Currently NULL
  # TODO: deparse(substritute(fnc)) just prints fnc

  message(paste0("FUNCTION CALL: ", deparse(substitute(fnc))))

  args = list(...)
  #if (!is.null(path_outFile)) cat(text = paste0("\nFUNCTION CALL: ", deparse(substitute(fnc))), file=path_outFile, append=T, sep="\n")

  message("ARGUMENTS:")
  #if (!is.null(path_outFile)) cat(text = "\nARGUMENTS:", file=path_outFile, append=T, sep="\n")

  # first print head / str / first cols and rows of args
  for (i in 1:length(args)) {
    anArg = args[[i]]
    print(paste0(names(args)[i]))
    #if (!is.null(path_outFile)) cat(text = paste0("\n",names(args)[i]) , file =  path_outFile, append=T, sep="\n")
    if (any(class(anArg) %in% c("data.frame", "data.table", "matrix", "Matrix", "data.frame", "dgCMatrix"))) {
      #if(ncol(anArg)>8) print(anArg[1:5,1:(min(nrow(anArg), 8))]) else head(anArg, n=5)
      print(anArg[1:3,1:(min(nrow(anArg), 5))])
      #if (!is.null(path_outFile)) cat(text = anArg[1:3,1:(min(nrow(anArg), 5))] , file =  path_outFile, append=T, sep="\n")

    } else if (any(c("standardGeneric", "function") %in% class(anArg))) {
      print(deparse(anArg, nlines=3))
      # TODO: can we concatenate something printed?
      #if (!is.null(path_outFile)) {
      #  cat(text = format(str(anArg, vec.len=4, nchar.max = 40, list.len=3)), file =  path_outFile, append=T, sep="\n")
    } else {
      str(anArg, vec.len=4, nchar.max = 40, list.len=3)
    }
  }


  # call fnc
  out <- do.call(what=fnc, args = args)

  # print value
  message("VALUE:")

  #if (!is.null(path_outFile)) cat(text = "\nVALUE: ", file=path_outFile, append=T, sep="\n")

  if (any(class(out) %in% c("data.frame", "data.table", "matrix", "Matrix", "data.frame", "dgCMatrix"))) {
    #if(ncol(out)>8) print(out[1:5,1:(min(nrow(out), 8))]) else head(out, n=5)
    #if (!is.null(path_outFile)) cat(text = out[1:3,1:(min(nrow(out), 5))], file =  path_outFile, append=T, sep="\n")
    print(out[1:3,1:(min(nrow(out), 5))])
  } else if (any(c("standardGeneric", "function") %in% class(out))) {
    print(deparse(out, nlines=3))
    # TODO: can we concatenate something printed?
    #if (!is.null(path_outFile)) {
    #  cat(text = format(str(anArg, vec.len=4, nchar.max = 40, list.len=3)), file =  path_outFile, append=T, sep="\n")
  } else {
    #if (!is.null(path_outFile)) { cat(text = format(str(out, vec.len=4, nchar.max = 40, list.len=4)), file =  path_outFile, append=T, sep="\n")
    #  } else {
    str(out, vec.len=4, nchar.max = 40, list.len=4)
    #   }
  }

  # print any warnings
  if (!is.null(warnings())) {
    #if (doPrint) {
    message("\nWARNINGS:")
    print(format(warnings()))
    #}
    #if (!is.null(path_outFile)) {
    #  cat(text="WARNINGS: ", file=path_outFile, append=T, sep="\n")
    #  cat(text = format(warnings()), file=path_outFile, append=T, sep="\n")
    #}
  }
  # return
  return(out)

}



#' Save file with metadata
#'
#' @description saves the following metadata
#' . filename (if savefnc is not NULL) or else
#' . current date and time
#' . sessionInfo() #sessioninfo::session_info()
#' . github log, if current or parent dir is a github dir
#' @param savefnc a function to write some file to disk. If given as object is converted to character, default NULL
#' @param doPrint print output to screen? useful if directin Rscript stdout to a log file (&>), default F
#' @param path_log specify log file path; defaults to creating a file in current working dir, default NULL
#' @param \dots named arguments to savefnc
#'
#' @return None
#' @export
#'
#' @examples utils_saveMeta(savefnc=ggplot2::ggsave, doPrint=F, path_log = "~/log.txt", width=10)
utils_saveMeta <- function(savefnc=NULL, doPrint=F, path_log=NULL, ...) {
  #' TODO add pander pandoc markdown output option?
  #' TODO add MD5SUM save
  #' TODO add overwrite argument

  # check args
  if (is.null(savefnc) & is.null(path_log) & !doPrint) stop("saveMeta: savefnc and path_log are NULL and doPrint is FALSE, no output")

  args = list(...)

  gsub("\\:", ".", gsub("\\ ", "_",as.character(Sys.time())))-> flag_date

  if (!is.null(savefnc)) {
    # convert savefnc arg to character
    if (!"character" %in% class(savefnc)) fnc <- as.character(quote(savefnc))

    # call savefnc to save the file
    do.call(what=savefnc, args=args, quote = F)

    # infer the path or file argument given to savefnc and make a modified path for the log file
    name_path_log <- sapply(c("con", "file", "path"), function(str) grep(str, names(args), value=T))
    name_path_log <- Filter(x=name_path_log,f=length)[[1]]
    path_log <- if ("character"%in%class(args[[name_path_log]])) {
      args[[name_path_log]]
    } else { #if not given as a character capture it otherwise
      args[[name_path_log]] %>% capture.output %>% '['(2) %>% gsub('description|\\ |\\"',"",.)
    }
    # replace the file extention suffix with datestamp and .txt file extention
    path_log <- gsub("\\.[a-z]+$", paste0("_meta_", substr(flag_date,1,10),".txt"),path_log, ignore.case = T)
  }

  # If savefnc is NULL and no path_log provided, save to a file in the working dir
  if (is.null(path_log)) path_log <- paste0(getwd(), "log_meta_", substr(flag_date,1,10),".txt")

  # Write to log file
  # session_info()
  utils::capture.output(sessioninfo::session_info()) %>% writeLines(text=.,con = path_log) #this has to come first since writeLines doesn't append
  #utils::capture.output(devtools::session_info()) %>% pander %>% cat(., file=path_log)
  cat("\n", file=path_log, append=T)
  cat(text="-Environment-----------------------------------", file=path_log, sep = "\n", append=T)
  capture.output(ls.str(envir = .GlobalEnv)) %>% cat(... = ., file=path_log, sep="\n", append=T)

  # if (!is.null(msg)) {
  #   path_msg = "~/msg_tmp.txt"
  #   write.table(x = msg, file = path_msg, quote = F, sep="\t")
  #   system2(command="cat", args=c(path_msg, path_log, paste0(" > ", path_log)), stdout=F, stderr=F)
  #   system2(command="rm", args=c(path_msg), stdout=F, stderr=F)
  # }

  # git commit: check parent directory for .git file
  path_parent <- gsub("[^/]*$","", path_log)
  path_parent <- substr(x=path_parent,1,nchar(path_parent)-1)
  path_parent <- gsub("[^/]*$","", path_parent)

  if (length(dir(path = path_parent, pattern="\\.git", all.files = T, recursive = F))>1) {
    try({
      gitCommitEntry <- try(system2(command="git", args=c("log", "-n 1 --oneline"), stdout=F, stderr = F))
      if (!"try-error" %in% class(gitCommitEntry)) {
        cat("\n", file=path_log, append=T)
        cat(text = "-Git commit------------------------------------", file =  path_log,  sep = "\n", append=T)
        cat(text = gitCommitEntry, file =  path_log, append=T, sep = "\n")
      }
    })
  } else {
    gitCommitEntry <- NULL
  }

  #date
  # cat("\n", file=path_log, append=T)
  # cat(text="-Date, time------------------------------------", file=path_log, sep = "\n", append=T)
  # cat(text = flag_date, file = path_log, append=T, sep = "\n")
  # cat(text="-----------------------------------------------", file=path_log, sep = "\n", append=T)

  if (doPrint) {
    # same as above, but print

    #date
    #message(x = flag_date, appendLF = T)
    #if (!is.null(msg)) print(msg)

    sessioninfo::session_info() %>% print

    message("")
    print(x="-Environment-----------------------------------")
    ls.str(envir = .GlobalEnv) %>% print

    if (!"try-error" %in% class(gitCommitEntry) & !is.null(gitCommitEntry)) {
      message("")
      print(x = "-Git commit------------------------------------")
      print(x = gitCommitEntry)
    }
  }
}

#' Convert a large sparse matrix into a dense matrix
#'
#' Avoid the error
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
#' @examples
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
