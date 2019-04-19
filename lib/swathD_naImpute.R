library(impute)

swathD_naImpute <- function(dat = NULL){

dat <- norm_all[[choose_norm]]
colnames(dat) <- sample_label  
ind <- which(is.na(dat), arr.ind=TRUE)
  
if(naImpute == "Zero replacement"){
  # Zero replacement
  dat[ind] <- 0
  df3 <<- dat
}else if(naImpute == "Minimum"){
  # Minimum (MinDet, deterministic minimum imputation)
  #dat[ind] <- apply(dat, 1, min, na.rm = TRUE)[ind[ ,1]]
  dat[ind] <- apply(dat, 2, min, na.rm = TRUE)[ind[ ,2]] # the left censored focusing
  df3 <<- dat
}else if(naImpute == "Median"){
  # Median      
  #dat[ind] <- apply(dat, 1, median, na.rm = TRUE)[ind[ ,1]]
  dat[ind] <- rowMedians(dat, na.rm = TRUE)[ind[ ,1]]
  df3 <<- dat
}else if(naImpute == "Average"){
  # Average      
  dat[ind] <- rowMeans(dat, na.rm = TRUE)[ind[ ,1]]
  df3 <<- dat
}else if(naImpute == "kNN"){
  # kNN impute
  dat[ind] <- impute.knn(dat, k = 10, rowmax = 0.5, colmax = 0.9, maxp = 1500, rng.seed = 1234)$data[ind[ ,1]]
  df3 <<- dat
}}