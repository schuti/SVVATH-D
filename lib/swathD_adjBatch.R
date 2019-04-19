

# for testing
    #df <- clean_ds[ , 5:length(clean_ds)] 
    #m_batch <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 1, 1, 1, 2, 2, 2, 3, 3, 3) 
    #m_batch <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)

# Adjust batch effect by ComBat 
swathD_adjBatch <- function(dat = 'NULL'){
  
 # require(sva)
  require(limma)
  
# get batch number from the metadata file  
  m_batch <- as.factor(group_label$batch)
  assign("m_batch", m_batch, envir = .GlobalEnv)
  dat <- dat %>% as.matrix()
# Adjust batch effect  
#  df2 <- ComBat(dat = df, batch = m_batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
  df3 <- removeBatchEffect(x = dat, batch = m_batch)
  
  # when ComBat produces values below zero (minus values), replaced each minus value by 1 (which will be 1 after log2 transform)
  for (i in seq_along(df3)){
    if(df3[i] < 0){
      df3[i] <- 1
    } else {
      df3[i] = df3[i]
    }}

# overwrite df for subsequent processes
  expr_processed <- as.data.frame(df3)
  assign("expr_processed", expr_processed, envir = .GlobalEnv)
  
}
