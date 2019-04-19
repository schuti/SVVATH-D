
swathD_hm_var <- function(x){

#dfmat <- as.matrix(log_ds[ , 2:length(log_ds)])
#med <- apply(t(dfmat), 1, mean)
#dfmat_medScale <- (t(dfmat) - med)

  
#hm_set <<- input$hm_set  
#varTopN <<- input$varTopN
#n_cluster <<- input$n_cluster

#if(hm_set == 'Significant features'){
 # dfvar <- dfmat_medScale_sig
#} else if (hm_set == 'All features'){
 # dfvar <- dfmat_medScale[-anova_pVal, -gene]
#}

rowVar <- apply(dfvar, 1, sd, na.rm=TRUE)

dfvar2 <- data.frame(dfvar, rowVar, gene = rownames(dfvar)) 
dfvar3 <- dfvar2 %>% arrange(desc(rowVar))
rownames(dfvar3) <- dfvar3$gene

#top20p <- ceiling(0.2*nrow(dfvar3))
#dfvar_top20p <- dfvar3 %>% head(n = top20p) %>% dplyr::select(-rowVar, -gene)

#topNp <- ceiling(as.numeric(input$varTopN) * nrow(dfvar3) )

topNp <- ceiling(as.numeric(varTopN/100) * nrow(dfvar3) )
#topNp <<- topNp

dfvar_topNp <- dfvar3 %>% head(n = topNp) %>% dplyr::select(-rowVar, -gene)

#chose_dist <- as.character(input$chose_dist)
#chose_dist <- chose_dist
#chose_link <- as.character(input$chose_link)
#chose_link <- chose_link

hm_var_tmp <- pheatmap(dfvar_topNp, breaks = seq(-(max(round(dfvar_topNp, 0))), max(round(dfvar_topNp, 0)), length.out=101), #breaks = seq(-3, 3, length.out=101), 
                   legend_breaks=seq(-(max(round(dfvar_topNp, 0))), max(round(dfvar_topNp, 0)), length.out=5),
                   #legend_labels = c("-1", "1e+02", "0", "1e+04", "1"),
                   color = colorRampPalette(c("darkblue", "blue", "white", "orangered", "red"))(100),  #color = colorRampPalette(c("darkblue", "cornflowerblue", "white", "orangered", "darkred"))(100),
                   # annotation_col = data.frame(group, row.names = sample_label),
                   annotation_col = data.frame(group = factor(group), row.names = sample_label),
                   #annotation_row = y2,
                   clustering_distance_rows = chose_dist,
                   clustering_distance_cols = chose_dist, #'correlation', 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'
                   clustering_method = chose_link, # 'ward', 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'.
                   #clustering_method_columns = "ward",
                   fontsize_row = 5, fontsize_col = 10, 
                   scale = "none",
                   cutree_rows = n_cluster, 
                   silent = TRUE
                   )
                  # main = paste0(nprot_sig, " proteins with anova_pVal < ", cutP, ' (color in Log2-scale)'))
#hm_var <<- hm_var

hm_var_id <- dfvar_topNp[hm_var_tmp$tree_row$order, ]

ind <- cutree(hm_var_tmp$tree_row, k = n_cluster)

gene_cluster <<- ind

y <- data.frame(cluster = ind)
y2 <- data.frame(cluster = as.factor(paste("cluster", ind)))
rownames(y2) <- rownames(y)


hm_var <- pheatmap(dfvar_topNp, breaks = seq(-(max(round(dfvar_topNp, 0))), max(round(dfvar_topNp, 0)), length.out=101), #breaks = seq(-3, 3, length.out=101), 
                       legend_breaks=seq(-(max(round(dfvar_topNp, 0))), max(round(dfvar_topNp, 0)), length.out=5),
                       #legend_labels = c("-1", "1e+02", "0", "1e+04", "1"),
                       color = colorRampPalette(c("darkblue", "blue", "white", "orangered", "red"))(100),  #color = colorRampPalette(c("darkblue", "cornflowerblue", "white", "orangered", "darkred"))(100),
                       border_color = NA,
                       # annotation_col = data.frame(group, row.names = sample_label),
                       annotation_col = data.frame(group = factor(group), row.names = sample_label),
                       annotation_row = y2,
                       clustering_distance_rows = "correlation",
                       clustering_distance_cols = "correlation", #'correlation', 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'
                       clustering_method = "average", # 'ward', 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'.
                       #clustering_method_columns = "ward",
                       fontsize_row = 5, fontsize_col = 10, 
                       scale = "none",
                       cutree_rows = n_cluster)

hm_var <<- hm_var
}
