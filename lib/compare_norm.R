# Compare_norm function

require(SwathXtend)
require(preprocessCore)
require(vsn)

compare_norm <- function(x){

expr_all <- expr_all

none <- expr_all %>% as.matrix() %>% log2() %>% as.data.frame() 
assign("none", none, envir = .GlobalEnv)

tmp <- sapply(expr_all, sum, na.rm=TRUE) %>% log10() %>% round() 
TAS <- expr_all %>% 
  lapply(function(x){x/sum(x, na.rm = TRUE)}) %>% 
  data.frame()
TAS <- (10^tmp * TAS) %>% as.matrix() %>% log2() %>% as.data.frame()
assign("TAS", TAS, envir = .GlobalEnv)

MLR <- expr_all %>% as.matrix() %>% log2() %>% mlrGroup(Group = group) 
MLR <- MLR[["mat.norm"]] %>% mlrrep()
MLR <- MLR[["mat.norm"]] %>% as.data.frame()
assign("MLR", MLR, envir = .GlobalEnv)

Quantile <- expr_all %>% as.matrix() %>% log2() %>% normalize.quantiles() %>% as.data.frame() 
colnames(Quantile) <- sample_label
assign("Quantile", Quantile, envir = .GlobalEnv)

VSN <- expr_all %>% as.matrix() %>% justvsn() %>% as.data.frame()  
assign("VSN", VSN, envir = .GlobalEnv)

norm_all <- list(none = as.matrix(none), 
                 TAS = as.matrix(TAS), 
                 Quantile = as.matrix(Quantile), 
                 MLR = as.matrix(MLR), 
                 VSN = as.matrix(VSN)
                 ) 

assign("norm_all", norm_all, envir = .GlobalEnv)

#plot_meanSd <- lapply(norm_all, meanSdPlot)
#source("./lib/qc_box.R", local = TRUE)
qc_box()
#source("./lib/qc_density.R", local = TRUE)
qc_density()

pmad <- list()

#for (i in seq_along(norm_all)){
 ##pmad[[i]] <- lapply(as.data.frame(norm_all[[i]]), mad, na.rm = TRUE) %>% unlist()
 #pmad[[i]] <- apply(as.data.frame(norm_all[[i]]), 1, function(x){mad(x, na.rm=TRUE)})

 
# names(pmad)[i] <- names(norm_all[i])
 
# pmad <<- pmad
#}

group_mad <- function(x){      
  tmp <- data.frame(group = group, t(x))
  tmp %>% gather(tmp_id, expression, -group) %>%
    dplyr::group_by(group, tmp_id) %>%
    dplyr::summarize(group_mad = mad(expression, na.rm = TRUE)) %>%
    .$group_mad
  }

pmad <- list()
for (i in seq_along(norm_all)){
  pmad[[i]] <-  group_mad(norm_all[[i]])
  names(pmad)[i] <- names(norm_all[i])
  pmad <<- pmad
} 



#pmad_plot <- graphics::boxplot(pmad, main = 'Pooled median absolute deviation (PMAD)', outline = FALSE)
#pmad_plot <<- pmad_plot

group_cv <- function(x){      
  tmp <- data.frame(group = group, t(x))
  tmp %>% gather(tmp_id, expression, -group) %>%
    dplyr::group_by(group, tmp_id) %>%
    dplyr::summarize(group_cv = 100*sd(expression, na.rm = TRUE)/mean(expression, na.rm = TRUE) ) %>%
    .$group_cv
}

pcv <- list()
for (i in seq_along(norm_all)){
  pcv[[i]] <-  group_cv(norm_all[[i]])
  names(pcv)[i] <- names(norm_all[i])
  pcv <<- pcv
} 

group_variance <- function(x){      
  tmp <- data.frame(group = group, t(x))
  tmp %>% gather(tmp_id, expression, -group) %>%
    dplyr::group_by(group, tmp_id) %>%
    dplyr::summarize(group_variance = sd(expression, na.rm = TRUE)^2 ) %>%
    .$group_variance
}

pev <- list()
for (i in seq_along(norm_all)){
  pev[[i]] <-  group_variance(norm_all[[i]])
  names(pev)[i] <- names(norm_all[i])
  pev <<- pev
} 

}

