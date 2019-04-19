

qc_density <- function(x){

require(grDevices)
colramp <- colorRampPalette(c("purple", "blue", "green", "orange", "red"), alpha = TRUE)(length(expr_all))
range <- c( (min(none[,1], na.rm = TRUE)-5), (max(none[,1], na.rm = TRUE)+5) )
yht <- (density(none[,1], na.rm = TRUE) %>% .$y %>% max(na.rm = TRUE))*1.5

plot(stats::density(none[, 1], na.rm = TRUE), xlab = "", lty = 3, col = colramp[1], main = "None", lwd = 2, xlim = range, ylim = c(0, yht)) %>% suppressWarnings()  #, xlim = c(0, 30), ylim = c(0, 0.25)
#ind <- which(is.na(none), arr.ind = TRUE)
#none[ind] <- 0
for(i in 2:ncol(none))(lines(stats::density(none[,i], na.rm = TRUE), lty = 3, col = colramp[i], lwd = 2))
density_none <- recordPlot()

plot(density(TAS[, 1], na.rm = TRUE), xlab = "", lty = 3, col = colramp[1], main = "TAS", lwd = 2, xlim = range, ylim = c(0, yht)) %>% suppressWarnings()   # , xlim = c(0, 30), ylim = c(0, 0.25)   c( (min(TAS[,1])-5), (max(TAS[,1])+5) )
for(i in 2:ncol(TAS))(lines(density(TAS[,i], na.rm = TRUE), lty = 3, col = colramp[i], lwd = 2))
density_TAS <- recordPlot()

plot(density(MLR[, 1], na.rm = TRUE), xlab = "", lty = 3, col = colramp[1], main = "MLR", lwd = 2, xlim = range, ylim = c(0, yht)) %>% suppressWarnings()   # , xlim = c(0, 30), ylim = c(0, 0.25)
for(i in 2:ncol(MLR))(lines(density(MLR[,i], na.rm = TRUE), lty = 3, col = colramp[i], lwd = 2))
density_MLR <- recordPlot()

plot(density(Quantile[, 1], na.rm = TRUE), xlab = "", lty = 3, col = colramp[1], main = "Quantile", lwd = 2, xlim = range, ylim = c(0, yht)) %>% suppressWarnings()  # , xlim = c(0, 30), ylim = c(0, 0.25)
for(i in 2:ncol(Quantile))(lines(density(Quantile[,i], na.rm = TRUE), lty = 4, col = colramp[i], lwd = 2))
density_Quantile <- recordPlot()

plot(density(VSN[, 1], na.rm = TRUE), xlab = "", lty = 3, col = colramp[1], main = "VSN", lwd = 2, xlim = range, ylim = c(0, yht)) %>% suppressWarnings()  # , xlim = c(0, 30), ylim = c(0, 0.25)
for(i in 2:ncol(VSN))(lines(density(VSN[,i], na.rm = TRUE), lty = 3, col = colramp[i], lwd = 2))
density_VSN <- recordPlot()

density_all <- list(none = density_none, 
                    TAS = density_TAS, 
                    MLR = density_MLR, 
                    Quantile = density_Quantile, 
                    VSN = density_VSN)
assign("density_all", density_all, envir = .GlobalEnv)
}