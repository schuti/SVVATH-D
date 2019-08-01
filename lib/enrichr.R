# Adopted from enrichr package with a bug fixed 
library(httr)

enrichr <- function(genes, databases = NULL) {
  cat("Uploading data to Enrichr... ")
  if (is.vector(genes) & ! all(genes == "") & length(genes) != 0) {
    temp <- POST(url="http://amp.pharm.mssm.edu/Enrichr3/enrich", # replacing Enrichr with Enrichr3
                 body=list(list=paste(genes, collapse="\n")))
  } else if (is.data.frame(genes)) {
    temp <- POST(url="http://amp.pharm.mssm.edu/Enrichr3/enrich", # replacing Enrichr with Enrichr3
                 body=list(list=paste(paste(genes[,1], genes[,2], sep=","),
                                      collapse="\n")))
  } else {
    warning("genes must be a non-empty vector of gene names or a dataframe with genes and score.")
  }
  GET(url="http://amp.pharm.mssm.edu/Enrichr3/share") # replacing Enrichr with Enrichr3
  cat("Done.\n")
  dbs <- as.list(databases)
  dfSAF <- options()$stringsAsFactors
  options(stringsAsFactors = FALSE)
  result <- lapply(dbs, function(x) {
    cat("  Querying ", x, "... ", sep="")
    r <- GET(url="http://amp.pharm.mssm.edu/Enrichr3/export", # replacing Enrichr with Enrichr3
             query=list(file="API", backgroundType=x))
    r <- gsub("&#39;", "'", intToUtf8(r$content))
    tc <- textConnection(r)
    r <- read.table(tc, sep = "\t", header = TRUE, quote = "", comment.char="")
    close(tc)
    cat("Done.\n")
    return(r)
  })
  options(stringsAsFactors = dfSAF)
  cat("Parsing results... ")
  names(result) <- dbs
  cat("Done.\n")
  return(result)
}
