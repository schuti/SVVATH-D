library(mixOmics)
library(FactoMineR)
library(NMF)

swathD_ml2enrichr <- function(x){

  
  if(chose_ml == "pca"){

  # PCA from FactoMineR package
    print("run PCA")
    fit_pca <- PCA(log_ds[ , 2:length(log_ds)], graph = FALSE, scale.unit = TRUE, ncp = 3)
    pca <- fit_pca[["var"]][["contrib"]] %>% as.data.frame()
    pca <- data.frame("gene" = rownames(pca), pca, stringsAsFactors = FALSE) 
    com1 <- pca %>% dplyr::arrange(desc(abs(pca$Dim.1))) %>% head(n = 50) %>%
                     dplyr::select("gene") %>% .$gene 
    com2 <- pca %>% dplyr::arrange(desc(abs(pca$Dim.2))) %>% head(n = 50) %>%
                     dplyr::select("gene") %>% .$gene 
    com3 <- pca %>% dplyr::arrange(desc(abs(pca$Dim.3))) %>% head(n = 50) %>%
                     dplyr::select("gene") %>% .$gene
    com_genes <- list(com1 = com1, com2 = com2, com3 = com3)
    com_genes <<- com_genes
  
    
    percentage <- fit_pca$eig[ , 2]
    
    PCs <- data.frame(fit_pca$ind$coord)
    #PCs <- data.frame(fit_pca$ind$cos2)
    PCs$group <- group
    
    plotPCA <- ggplot(data = PCs, aes(x = Dim.1, y = Dim.2)) +
      geom_point(aes(colour = group), size = 3) +
      xlab(paste0('PC1', ' ', '(', round(percentage[1], 2), '%)')) + 
      ylab(paste0('PC2', ' ', '(', round(percentage[2], 2), '%)')) +
      scale_fill_hue(l=40) + 
      coord_fixed(ratio=1, xlim=range(PCs$Dim.1), ylim=range(PCs$Dim.2)) +
      geom_text_repel(label = rownames(PCs)) +
      theme_light()
    
  } else if(chose_ml == 'sipca'){
  #sIPCA from mixOmics package
    print("run sIPCA")
    tmp <- log_ds[ , 2:length(log_ds)] %>% scale()
    tmp[which(is.na(tmp), arr.ind = TRUE)] <- 0
    fit_sipca <- sipca(tmp, mode = "deflation", scale = FALSE, 
                       ncomp = 3, keepX = c(50, 50, 50))
    
    com1 <- selectVar(fit_sipca, comp = 1)$value %>% rownames()
    com2 <- selectVar(fit_sipca, comp = 2)$value %>% rownames()
    com3 <- selectVar(fit_sipca, comp = 3)$value %>% rownames()
    com_genes <- list(com1 = com1, com2 = com2, com3 = com3)
    com_genes <<- com_genes
  
  } else if(chose_ml == 'nmf'){
  # NMF from NMF package
    print("run NMF with 100 iterations")
    tmp <- t(log_ds[ , 2:length(log_ds)]) %>% as.data.frame()
    fit_nmf <- nmf(tmp, rank = 3,  method = 'brunet', seed = 123456, nrun=100) 
    fs <- featureScore(fit_nmf, method = 'kim')
    ef <- extractFeatures(fit_nmf, method = 'kim', format = 'list')
    com1 <- fs[ef[[1]]]
    com1 <- names(com1) #%>% head(n = 100)
    com2 <- fs[ef[[2]]]
    com2 <- names(com2) #%>% head(n = 100)
    com3 <- fs[ef[[3]]]
    com3 <- names(com3) #%>% head(n = 100)
    com_genes <- list(com1 = com1, com2 = com2, com3 = com3)
    com_genes <<- com_genes
  }
  
  mat_com1 <- dfmat_medScale_all[com_genes$com1, ]
  mat_com2 <- dfmat_medScale_all[com_genes$com2, ]
  mat_com3 <-dfmat_medScale_all[com_genes$com3, ]
  mat_coms <- list(mat_com1 = mat_com1, mat_com2 = mat_com2, mat_com3 = mat_com3)
  mat_coms <<- mat_coms
  
  sig_com1 <- dfmat_medScale_sig[com_genes$com1, ]
  ind <- which(!is.na(sig_com1), arr.ind = TRUE)
  sig_com1 <- sig_com1[unique(ind[ ,1]), ]
  sig_com2 <- dfmat_medScale_sig[com_genes$com2, ]
  ind <- which(!is.na(sig_com2), arr.ind = TRUE)
  sig_com2 <- sig_com2[unique(ind[ ,1]), ]
  sig_com3 <-dfmat_medScale_sig[com_genes$com3, ]
  ind <- which(!is.na(sig_com3), arr.ind = TRUE)
  sig_com3 <- sig_com3[unique(ind[ ,1]), ]
  sig_coms <- list(sig_com1 = sig_com1, sig_com2 = sig_com2, sig_com3 = sig_com3)
  sig_coms <<- sig_coms
  
  com_sig_genes <- lapply(sig_coms, rownames)
  com_sig_genes <<- com_sig_genes
  
}



#x_dim1_gn_enrichr <- enrichr(x_dim1_gn, "KEA_2015")
#x_dim1 <- x_dim1_gn_enrichr[["KEA_2015"]]
#x_dim2_gn_enrichr <- enrichr(x_dim2_gn, "KEA_2015")
#x_dim2 <- x_dim2_gn_enrichr[["KEA_2015"]]
#x_dim3_gn_enrichr <- enrichr(x_dim3_gn, "KEA_2015")
#x_dim3 <- x_dim3_gn_enrichr[["KEA_2015"]]

#comVar_enrichr <- enrichr(comVar_gene, "KEA_2015")
#comVar1 <- comVar_enrichr[["KEA_2015"]]
#comVar_enrichr <- enrichr(comVar_gene, "KEA_2015")
#comVar2 <- comVar_enrichr[["KEA_2015"]]
#comVar_enrichr <- enrichr(comVar_gene, "KEA_2015")
#comVar3 <- comVar_enrichr[["KEA_2015"]]

#y3_enr <- enrichr(y3_gn, "KEA_2015")
#z <- y3_enr[["KEA_2015"]]
#y3_2_enr <- enrichr(y3_2_gn, "KEA_2015")
#z2 <- y3_2_enr[["KEA_2015"]]
#y3_3_enr <- enrichr(y3_3_gn, "KEA_2015")
#z3 <- y3_3_enr[["KEA_2015"]]