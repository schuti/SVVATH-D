

swathD_imputeFirst <- function(dat = NULL){
  
  impBeforeNorm <- as.matrix(expr_all)
  impBeforeNorm[impBeforeNorm == 0] <- NA
  colnames(impBeforeNorm) <- sample_label
  ind <- which(is.na(impBeforeNorm), arr.ind=TRUE)
  
  if(naImpute == "Minimum"){
    # Minimum 
    impBeforeNorm[ind] <- apply(impBeforeNorm, 1, min, na.rm = TRUE)[ind[ ,1]]
  }else if(naImpute == "Median"){
    # Median      
    #impBeforeNorm[ind] <- apply(impBeforeNorm, 1, median, na.rm = TRUE)[ind[ ,1]]
    impBeforeNorm[ind] <- rowMedians(impBeforeNorm, na.rm = TRUE)[ind[ ,1]]
  }else if(naImpute == "Average"){
    # Average      
    impBeforeNorm[ind] <- rowMeans(impBeforeNorm, na.rm = TRUE)[ind[ ,1]]
  }else if(naImpute == "kNN"){
    # kNN impute
    impBeforeNorm[ind] <- impute.knn(impBeforeNorm, k = 10, rowmax = 0.5, colmax = 0.9, maxp = 1500, rng.seed = 1234)$data[ind[ ,1]]
  }
  expr_all <<- impBeforeNorm %>% as.data.frame()
}
