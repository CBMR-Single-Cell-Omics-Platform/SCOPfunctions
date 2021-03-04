# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


f_rmPkgs <- function() {
  #' @usage detach and unload all packages
  #' @return none
  #' @references https://stackoverflow.com/questions/7505547/detach-all-packages-while-working-in-r
  while (!is.null(names(sessionInfo()$otherPkgs))) lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
}

f_dlfile <- function(url, destfile, correct_checksum) {
  #' @usage: checks if destfile is present; if not, downloads and checks md5sum against correct_md5sum
  #' @param url: passed to download.file().
  #' @param destfile: passed to download.file().
  #' @param correct_checksum: to check file integrity
  #' @value: none

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


f_top_objectsizes <- function(n) {
  #' @usage show the memory size of the top n objects
  #' @param n how many objects to list
  #' @return NULL
  sapply(ls(), function(x) object.size(eval(parse(text=x)))) %>% sort(., decreasing=T) %>% head(., n)
}

f_detectCores_plus <- function(Gb_max=250,
                                additional_Gb=1) {
  #' @param Gb_max ceiling on session memory usage in Gb, assuming that each worker duplicates the session memory
  #' @param additional_Gb: max additional memory requirement for new (temporary) objects created within a parallel session
  #' @returns: max number of cores (integer)
  #' @depends parallel package

  require("parallel")
  obj_size_Gb <- as.numeric(sum(sapply(ls(envir = .GlobalEnv), function(x) object.size(x=eval(parse(text=x))))) / 1024^3)
  max(1, min(parallel::detectCores(), Gb_max %/% (obj_size_Gb + additional_Gb))-1)
}



f_verboseFnc <- function(fnc,
                       args) {
  #' @usage Function wrapper to show head of inputs and outputs
  #' @param fnc: a function object
  #' @param args: a list of named arguments
  #' @value value returned by fnc, if any

  # TODO: fix output to outfile. Currently NULL
  # TODO: deparse(substritute(fnc)) just prints fnc

  message(paste0("FUNCTION CALL: ", deparse(substitute(fnc))))

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

#' @TODO add pander pandoc markdown output option?
#' @TODO add MD5SUM save
#' @TODO add overwrite argument

saveMeta <- function(savefnc=NULL, doPrint=F, path_log=NULL, ...) {
  #' @usage If savefnc provided, write a timestamped metadata file along with file; else write log anyway
  #'        If path_log is not given, i.e. NULL, makes a file in the working directory
  #'        Optionally print to screen
  #'        Includes
  #'         filename (if savefnc is not NULL) or else
  #'         current date and time
  #'         sessionInfo() #devtools::session_info()
  #'         github log, if current or parent dir is a github dir
  #' @param savefnc a function to write some file to disk. If given as object is converted to character, default NULL
  #' @param doPrint print output to screen? useful if directin Rscript stdout to a log file (&>), default F
  #' @param path_log specify log file path; defaults to creating a file in current working dir, default NULL
  #  #' @param msg string,  message to add at the top of the log file
  #' @param ... arguments to pass on to savefnc
  #' @value NULL
  #' @examples
  #' saveMeta(savefnc= ggsave, plot=p, filename= "/projects/myname/RNAplot.pdf", height=10,width=12)
  #' saveMeta(savefnc= save.image,  file= "/projects/myname/session.image.Rdata.gz", compress="gzip")
  #' saveMeta(doPrint=T)

  # packages
  require(magrittr)
  require(utils)
  require(devtools) # for devtools::session_info()
  #require(pander)

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
  utils::capture.output(devtools::session_info()) %>% writeLines(text=.,con = path_log) #this has to come first since writeLines doesn't append
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

    devtools::session_info() %>% print

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
