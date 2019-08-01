library(dplyr)
library(visNetwork)


swathD_ml2net <- function(dat = NULL, dbs = NULL, 
                          cutPval = 0.05, cutScore = 2, cutLinkage = 1){


# node color indicate levels of significant
  # where more red = more significant
# node size indicate the number of genes
# edge width is a proportion to the geneRatio

#chose_ml2net <- com_input[[chose_com]]
  
#chose_ml2net_com <- input$chose_ml2net_com
#chose_ml2net <- com_input[[chose_ml2net_com]]

#ml2net_dbs <- input$ml2net_dbs

genes <- dat  
database <- dbs
adjPval <- cutPval
rScore <- cutScore
nLinkage <- cutLinkage

re_enrichr4 <- enrichr(genes = genes, database = database
                       #database = ml2net_dbs)
                       #databases = list("GO_Biological_Process_2015", 
                        #            "KEGG_2016", "Reactome_2016",
                         #           "KEA_2015", "Transcription_Factor_PPIs")
                         )

re_enrichr4 <<- re_enrichr4

re_enrich4_table <- bind_rows(re_enrichr4) %>% as.data.frame()

# filter by adjusted p value ------------------------------------------------------
if(min(re_enrich4_table$Adjusted.P.value) >= adjPval){
  warning_adjPval <- paste("No matched function/process at the indicated adjusted P-value")
  warning_adjPval <<- warning_adjPval
  print("No matched function/process at the indicated adjusted P-value")
  if(exists("ml2net_vis")){
    rm(ml2net_vis, envir = .GlobalEnv)}
} else if(min(re_enrich4_table$Adjusted.P.value) < adjPval){
  re_enrich4_table <- re_enrich4_table %>% dplyr::filter(Adjusted.P.value < adjPval)
  print("Pass the adjusted P-value filtering")
    if(exists("warning_adjPval")){
      rm(warning_adjPval, envir = .GlobalEnv)}
  
#----------------------------------------------------------------------------------

re_enrich4_table <- tidyr::separate(re_enrich4_table, Overlap, c("geneInput", "geneTotal"), sep = "/")
re_enrich4_table$geneInput <- as.numeric(re_enrich4_table$geneInput)
re_enrich4_table$geneTotal <- as.numeric(re_enrich4_table$geneTotal)
re_enrich4_table <- re_enrich4_table %>% mutate("Pathway_fraction_observed" = (geneInput / geneTotal)) %>% dplyr::arrange(Adjusted.P.value)
#re_enrich4_table <- re_enrich4_table %>% mutate("rankScore" = (geneRatio) * (1/Adjusted.P.value) )
re_enrich4_table <- re_enrich4_table %>% mutate("rankScore" = (100*Pathway_fraction_observed) * (-log10(Adjusted.P.value)) )
re_enrich4_table <<- re_enrich4_table

# filter by rankScore -------------------------------------------------------------
if(sum(re_enrich4_table$rankScore) <= rScore){
  warning_rankScore <- paste("No matched function/process at the indicated rankScore")
  warning_rankScore <<- warning_rankScore
  print("No matched function/process at the indicated rankScore")
  if(exists("ml2net_vis")){
    rm(ml2net_vis, envir = .GlobalEnv)}
} else if(sum(re_enrich4_table$rankScore) > rScore){
  enrich_edge_cutoff <-  re_enrich4_table %>% dplyr::filter(rankScore > rScore)
  print("Pass the rankScore filtering")
  if(exists("warning_rankScore")){
    rm(warning_rankScore, envir = .GlobalEnv)}
    

#----------------------------------------------------------------------------------

#sumScore <- sum(re_enrich4_table$rankScore)
#re_enrich4_table$rankScore <- 100*(re_enrich4_table$rankScore / sumScore)
#re_enrich4_table2 <- re_enrich4_table %>% dplyr::arrange(desc(rankScore))

#edge_cutoff <- median(re_enrich4_table$geneRatio)
#enrich_edge_cutoff <- re_enrich4_table[re_enrich4_table$geneRatio >= median(re_enrich4_table$geneRatio), ]

#enrich_edge_cutoff <- re_enrich4_table %>% dplyr::filter(Adjusted.P.value < 0.05)

#tmp <- split(enrich_edge_cutoff, seq(nrow(enrich_edge_cutoff))) %>% set_names(enrich_edge_cutoff$Term)

tmp_genes <- lapply(enrich_edge_cutoff$Genes, strsplit, ";", fixed = TRUE) #, names = as.character(enrich_edge_cutoff$Term))
tmp_genes <- set_names(tmp_genes, enrich_edge_cutoff$Term)

tmp_term <- enrich_edge_cutoff$Term

#long_tmp <- tmp <- unlist(tmp_genes, recursive = FALSE, use.names = TRUE)

long_tmp <- list()
for(i in 1:length(tmp_genes)){
long_tmp[i] <- list(data.frame("from" = tmp_genes[[i]][[1]] , 
                               "to" = rep(tmp_term[i], length(tmp_genes[[i]][[1]])), 
                               stringsAsFactors = FALSE))
}

long_tmp2 <- bind_rows(long_tmp)

#tmp_pval <- enrich_edge_cutoff$Adjusted.P.value
#tmp_geneRatio <- enrich_edge_cutoff$geneRatio
#tmp_pval <- data_frame("From" = tmp_term, "Adj.P.value" = tmp_pval)

tmp <- enrich_edge_cutoff %>% dplyr::select('Term', 'Adjusted.P.value', 'Pathway_fraction_observed')
colnames(tmp) <- c("to", "score", "width")
#tmp$width <- round(((5*tmp$width) * as.numeric(enrich_edge_cutoff$geneInput))^2, 1)

merge <- left_join(long_tmp2, tmp, by = "to")

tmp_n <- merge %>% dplyr::group_by(to) %>% dplyr::summarize(n_linkage = n()) 
merge <- left_join(merge, tmp_n, by = 'to')

# filter by n_linkage -------------------------------------------------------------
if(max(merge$n_linkage) < nLinkage){
  warning_n_linkage <- paste("No matched function/process at the indicated numbers of linkage")
  warning_n_linkage <<- warning_n_linkage
  print("No matched function/process at the indicated numbers of linkage")
  if(exists("ml2net_vis")){
    rm(ml2net_vis, envir = .GlobalEnv)}
  } else if(max(merge$n_linkage) >= nLinkage){
  ml2net_edge <<- merge %>% dplyr::filter(n_linkage >= nLinkage)
  print("Pass the numbers of linkage filtering")
    if(exists("warning_n_linkage")){
      rm(warning_n_linkage, envir = .GlobalEnv)}

#----------------------------------------------------------------------------------


node <- unique(c(ml2net_edge$from, ml2net_edge$to)) #, chose_ml2net))
x <- node %in% tmp_term
y <- node[x]
z <- node[!x]

group <- toupper(node) %>% replace(. %in% toupper(y), "Function") %>%
  replace(. %in% toupper(z), "Gene involved") #%>%
 # replace(!. %in% c("Function", "Involved gene"), "Not involved")

ml2net_nodes <- data.frame(label = node, id = node, group = group, stringsAsFactors = FALSE) #, title = link_modnodes)

tmp_size <- enrich_edge_cutoff %>% dplyr::select("id" = "Term", "value" = "geneInput", "color" = "Adjusted.P.value")

ml2net_nodes <- left_join(ml2net_nodes, tmp_size, by = 'id')

#ml2net_nodes$value <- (ml2net_nodes$value)

ind <- which(is.na(ml2net_nodes$color), arr.ind = TRUE)

ml2net_nodes_na <- ml2net_nodes[ind, ]
ml2net_nodes_na$color <- "#d9d9d9"
ml2net_nodes_col <- ml2net_nodes[-ind, ]
#ml2net_nodes_rank <- ml2net_nodes_col %>% mutate(rank = rank(ml2net_nodes_col$color))
ml2net_nodes_rank2 <- ml2net_nodes_col %>% mutate(rank = colorRampPalette(c('red', 'yellow'))(length(ml2net_nodes_col$color))[rank(ml2net_nodes_col$color)])
ml2net_nodes_rank2 <- ml2net_nodes_rank2[, -5]
colnames(ml2net_nodes_rank2) <- c("label", "id", "group", "value", "color")
ml2net_nodes <- bind_rows(ml2net_nodes_na, ml2net_nodes_rank2)

ind <- which(is.na(ml2net_nodes$value), arr.ind = TRUE)
ml2net_nodes[ind, 4] <- 1

#ml2net_nodes[ind, 5] <- 1
#ml2net_nodes[ind, 5] <- "#d9d9d9" # input genes from the chosen component 

#ind1 <- suppressWarnings( which(as.numeric(ml2net_nodes$color) >= 0.05) )
#ml2net_nodes[ind1, 5] <- "black"
#ind2 <- suppressWarnings( which(as.numeric(ml2net_nodes$color) < 0.05 & as.numeric(ml2net_nodes$color) >= 0.005) )
#ml2net_nodes[ind2, 5] <- "#ffcccc"
#ind3 <- suppressWarnings( which(as.numeric(ml2net_nodes$color) < 0.005 & as.numeric(ml2net_nodes$color) >= 0.0005) )
#ml2net_nodes[ind3, 5] <- "ffb3b3"
#ind4 <- suppressWarnings( which(as.numeric(ml2net_nodes$color) < 0.0005 & as.numeric(ml2net_nodes$color) >= 0.00005) )
#ml2net_nodes[ind4, 5] <- "ff9999"
#ind5 <- suppressWarnings( which(as.numeric(ml2net_nodes$color) < 0.00005 & as.numeric(ml2net_nodes$color) >= 0.000005) )
#ml2net_nodes[ind5, 5] <- "ff8080"
#ind6 <- suppressWarnings( which(as.numeric(ml2net_nodes$color) < 0.000005) )
#ml2net_nodes[ind6, 5] <- "ff4d4d"

#x <- ml2net_nodes %>% replace(. < 0.05 & . > 0.005,  "white") 

ml2net_nodes <- ml2net_nodes %>% dplyr::group_by(label, id) %>% slice(1)
ml2net_nodes <<- ml2net_nodes

ml2net_vis <- visNetwork(ml2net_nodes, ml2net_edge) %>%
                visNodes(physics = FALSE) %>%
                visEdges(physics= FALSE, smooth = TRUE, color = list(color = "green", highlight = "purple")) %>% #, color = list(color = "#008080", highlight = "purple")) %>%
                visGroups(groupname = "Gene involved", shape = "dot", #size = 1,
                          font =  list(color = "magenta", size = 25),
                          color = list(shadow = TRUE, background = "ddd", border = "black",
                                       highlight = list(background = "feff8f", 
                                                        border = "red")) ) %>%
                visGroups(groupname = "Function", shape = "dot", #size = 25,
                          font =  list(color = "blue", size = 25),
                          color = list(background = "yellow", border = "orange",
                                       highlight = list(background = "feff8f", 
                                                        border = "red")) ) %>%
               # visHierarchicalLayout(enabled = TRUE, levelSeparation = 150, treeSpacing = 100, nodeSpacing = 100, blockShifting = FALSE, edgeMinimization = FALSE) #%>%
                visIgraphLayout(type = "full", randomSeed = 9, layout = "layout_nicely", physics = FALSE, smooth = FALSE) #%>% 
                #visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE, selectedBy = "group") %>%
                #visInteraction(navigationButtons = TRUE)
                
ml2net_vis <<- ml2net_vis

#ml2net_save <- ml2net_vis %>% #visHierarchicalLayout(enabled = TRUE, levelSeparation = 75, treeSpacing = 300, nodeSpacing = 300, blockShifting = FALSE, edgeMinimization = FALSE) %>%
 # visIgraphLayout(type = "square", randomSeed = 9, layout = "layout_nicely", physics = FALSE, smooth = FALSE) %>% 
  #visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE, selectedBy = "group") %>%
  #visInteraction(navigationButtons = TRUE)

#visSave(ml2net_save, file = "~/Desktop/molm13_EV_pca_con2_layoutNicely.html", background = "white")
  
   #}
  } #from n_linkage
  } #from rankScore
  } #from adjPval
}