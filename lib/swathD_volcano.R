library(ggplot2)
library(ggrepel)


swathD_volcano <- function(data = 'NULL',
                           re_cutFC = 2,
                           re_cutP = 0.05,
                           #cutFC_label = 4,
                           cutFC_label = 2,
                           cutP_label = 4,
                           sig_col = "red",
                           non_sig_col = "grey",
                           range_x = c(-5, 5),
                           range_y = c(0, 5)){

#cutFC <- 2
#cutP <- 0.05

#cutFC_label <- 4
#cutP_label <- 4

#sig_col <- "red"
#non_sig_col <- "grey"

#data <- read.csv("~/Desktop/long_pfc.csv") 

#long_p.fc <- data

#long_p.fc[ ,2] <- -log10(long_p.fc[ , 2])
#long_p.fc[ , 3] <- log2(long_p.fc[ ,3])
#colnames(as.data.frame(long_p.fc)) <- c("id", "pVal", "FC")

volcano <- ggplot(data = long_p.fc, aes(x= log2FC, 
                                        y=-log10(adj_pVal))) +
#volcano <- ggplot(data = long_p.fc, aes(x= log2(FC), 
#                                        y=-log10(pVal))) +
  #                     # geom_point(aes(color = as.factor(abs(log2FC) > log2(2) & anova.pVal < 0.05 & adj_pVal < 0.05)), size = 1) +
  geom_point(aes(color = as.factor(abs(log2FC) >= log2(re_cutFC) & anova.pVal < cutP & adj_pVal < re_cutP)), size = 3, alpha = 0.5) +  
  scale_color_manual(values = c("grey", "red")) +
  theme(legend.position = "none") +
  xlim(range_x) + 
  ylim(range_y) +
  #scale_x_continuous(breaks = c(-6, -3, -1, 0, 1, 3, 6)) +
  xlab("log2 (fold change)") + ylab("-log10 (p-value)") +
  ggtitle(label = paste0(pair_select, ": Volcano plot at ", re_cutFC, "x fold change and P-value < ", re_cutP)) + 
  theme(plot.title = element_text(lineheight= 0.8, face="bold")) +
  geom_text_repel(data = (subset(long_p.fc, 
                                 abs(log2FC) > cutFC_label & 
                                 -log10(adj_pVal) > cutP_label)),
                  aes(label = gene, size = 0.25),
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines")) 

volcano 

}

#x <- long_p.fc$fold.change
#y <- 2^x
#long_p.fc$fold.change <- y
#write.csv(long_p.fc, file = "long_pfc.csv")
