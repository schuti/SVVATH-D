

qc_box <- function(x){
  require(graphics)
  box_all <- list()
  for (i in 1:length(norm_all)){
    box_all[i] <- lapply(norm_all[i], graphics::boxplot, main = names(norm_all)[i], outline = FALSE, las = 2,
                    # ylim = c(0, (max(norm_all[[i]][,1])+min(norm_all[[i]][,1]))),
                     col = colorRampPalette(c('purple', 'lightblue', 'green', 'orange', 'red'))(length(unique(group)))[group]) %>% 
                  lapply(recordPlot)
  }
  names(box_all) <- names(norm_all)
  box_all <<- box_all
  }
  