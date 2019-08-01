
swathD_idConv <- function(dataformat){

#dataformat <- input$dataFormat  

if(species == "Human"){
  chose_species <<- hsapiens_ensembl_012219
} else if(species == "Mouse"){
  chose_species <<- mmusculus_ensembl_012219
} else if(species == "Rat"){
  chose_species <<- rnorvegicus_ensembl_030319
}
    
if(dataformat == 'Peakview SWATH full report'){
  
  areaProt <- read_excel(data_path, sheet = "Area - proteins")
  
  # remove cont and reverse seq, if presented
  ind <- as.character(areaProt$Protein) %>% str_detect("cont")
  areaProt <- areaProt[!ind, ]
  ind <- as.character(areaProt$Protein) %>% str_detect("RRRRR")
  areaProt <- areaProt[!ind, ]
  
  areaProt[ , 2:length(areaProt)] <- lapply(areaProt[ , 2:length(areaProt)], as.numeric)
  
  areaProt <<- areaProt
  
  #df <- areaProt[ , 2:length(areaProt)]
  
  areaPept <- read_excel(data_path, sheet = "Area - peptides")
  ind <- as.character(areaPept$Protein) %>% str_detect("cont")
  areaPept <- areaPept[!ind, ]
  ind <- as.character(areaPept$Protein) %>% str_detect("RRRRR")
  areaPept <- areaPept[!ind, ]
  
  areaPept[ , 3:length(areaPept)] <- lapply(areaPept[ , 3:length(areaPept)], as.numeric)
  
  areaPept <<- areaPept  
  
  df <- areaProt[ ,1] %>% 
    separate(Protein, c("sp", "uniProtID", "entry_name"), sep = "\\|") %>%
    separate(entry_name, c("entry_names", "species"), sep = "_") %>%
    dplyr::select(uniProtID, entry_names, species) 
  
  incProgress(0.1, message = "UniProt ID-Gene mapping")
  
#  ensembl=useMart("ensembl", 
 #                 dataset="hsapiens_gene_ensembl", # will add other species later
  #                host = "useast.ensembl.org",
   #             #  host = "www.ensembl.org", # force to use the main emsembl site to avoid the getBM issues
    #              ensemblRedirect = FALSE)
  
  tmp <- getBM(attributes = c('uniprotswissprot', 'external_gene_name'), #, 'description'), 
               filters = 'uniprotswissprot', 
               values = df$uniProtID, 
            #   mart = ensembl) 
               mart = chose_species)
    
  # detect species before matching IDs by the appropriate dbs
#  species <- df$species[1]
#  if(species == 'MOUSE'){
#    load("~/Desktop/RShiny/SwathShiny/lib/mmusculus_ensembl_012219") # for M.musculus 
#    tmp <- getBM(attributes = c('uniprotswissprot', 'external_gene_name'), 
#                 filters = 'uniprotswissprot', 
#                 values = df$uniProtID, 
#                 mart = mmusculus_ensembl_012219) 
#  } else if(species == 'HUMAN') {
#    load("~/Desktop/RShiny/SwathShiny/lib/hsapiens_ensembl_012219") # for H.sapiens
#    tmp <- getBM(attributes = c('uniprotswissprot', 'external_gene_name'), 
#                 filters = 'uniprotswissprot', 
#                 values = df$uniProtID, 
#                 mart = hsapiens_ensembl_012219) 
#  } else if(species == 'RAT') {
#    load("~/Desktop/RShiny/SwathShiny/lib/rnorvegicus_ensembl_030319") 
#    tmp <- getBM(attributes = c('uniprotswissprot', 'external_gene_name'), # for R.norvegicus
#                 filters = 'uniprotswissprot', 
#                 values = df$uniProtID, 
#                 mart = rnorvegicus_ensembl_012219) 
#  } else if(species == 'CANLF') { 
#    load("~/Desktop/RShiny/SwathShiny/lib/cfamiliaris_ensembl_030319") # for C.familiaris
#    tmp <- getBM(attributes = c('uniprotswissprot', 'external_gene_name'), 
#                 filters = 'uniprotswissprot', 
#                 values = df$uniProtID, 
#                 mart = cfamiliaris_ensembl_030319) 
#  }
  
  colnames(tmp) <- c('uniProtID', "gene.SYMBOL")
  
  df <- left_join(df, tmp[!duplicated(tmp$uniProtID), ], by = "uniProtID") 
  
  # for a given missing gene name, replaced by entry_names 
  ind <- is.na(df$gene.SYMBOL)
  df$gene.SYMBOL[ind] <- df$entry_names[ind]

#-- VERY SLOW -- for a number of UniProt ID that gene names are missing in ensembl, go get them one-by-one from UniProt --------------
#  gn_na <- as.list(unique(df$uniProtID[ind]))
#  url <- "https://www.uniprot.org/uniprot/"
#  api <- lapply(gn_na, function(x){paste0(url, x, ".txt")})
  
#  getGN <- function(x){
#    y <- read_tsv(x, col_names = FALSE)
#    z <- y[str_detect(y$X1, 'GN   Name'), ] %>% as.character() %>% strsplit("=", fixed = TRUE) %>% 
#      sapply("[", 2) %>% strsplit(";", fixed = TRUE) %>% sapply("[", 1) %>% strsplit(" ", fixed = TRUE) %>% sapply("[", 1)
#  }
  
#  x <- suppressMessages(lapply(api, getGN))
#  names(x) <- df$uniProtID[ind]
#  x2 <- unlist(x) %>% data.frame(uniProtID = df$uniProtID[ind], gene.SYMBOL = ., stringsAsFactors = FALSE)
  
  # Finally, if no gene name in both ensembl and UniProt, replace NA with its UniProt ID
#  ind2 <- is.na(x2$gene.SYMBOL)
#  x2$gene.SYMBOL[ind2] <- x2$uniProtID[ind2]
  
#  df$gene.SYMBOL[ind] <- x2$gene.SYMBOL
#-----------------------------------------  
  
  id_all <<- df
  
  df2 <- areaProt[ , 2:length(areaProt)]
  df2[df2 == 0] <- NA
  colnames(df2) <- sample_label
  expr_all <<- df2
  
} else if(dataformat == 'Processed dataset [gene, expr]'){
  
  incProgress(0.1, message = "User's defined dataset")
  
  #tmp <- read_excel(data_path, sheet = 'Area - peptides') 
  tmp <- read.csv(data_path, header = TRUE, stringsAsFactors = FALSE) 
  
  areaProt <- tmp
  colnames(areaProt) <- c('gene.SYMBOL', sample_label)
  
  areaProt[ , 2:length(areaProt)] <- lapply(areaProt[ , 2:length(areaProt)], as.numeric)
  
  areaProt <<- areaProt # use later for the number of peptide/protein
  
  ind2 <- which(areaProt == 0, arr.ind = TRUE)
  areaProt[ind2] <- NA
  
  expr_all <- areaProt[, 2:length(areaProt)]
  colnames(expr_all) <- sample_label
  expr_all <<- expr_all
  
  #id_all <- areaProt[, 1] 
  id_all <- areaProt[,1] %>% data.frame()
  colnames(id_all) <- 'gene.SYMBOL'
  id_all <<- id_all

} else if(dataformat == 'Processed dataset [gene, peptide, expr]'){
  
  incProgress(0.1, message = "User's defined dataset")
  
  #tmp <- read_excel(data_path, sheet = 'Area - peptides') 
  tmp <- read.csv(data_path, header = TRUE, stringsAsFactors = FALSE) 
  
  areaPept <- tmp
  colnames(areaPept) <- c('gene.SYMBOL', 'peptide', sample_label)
  
  areaPept[ , 3:length(areaPept)] <- lapply(areaPept[ , 3:length(areaPept)], as.numeric)
  
  areaPept <<- areaPept # use later for the number of peptide/protein
  
  tmp2 <- areaPept[, c(1, 3:length(areaPept))]
  areaProt <- tmp2 %>% dplyr::group_by(gene.SYMBOL) %>% dplyr::summarise_at(.vars = sample_label, .funs = sum) 
  areaProt <<- areaProt # use later for the number of peptide/protein
  
  ind2 <- which(areaProt == 0, arr.ind = TRUE)
  areaProt[ind2] <- NA
  
  expr_all <- areaProt[, 2:length(areaProt)]
  colnames(expr_all) <- sample_label
  expr_all <<- expr_all
  
  #id_all <- areaProt[, 1] 
  id_all <- areaProt[,1] %>% data.frame()
  colnames(id_all) <- 'gene.SYMBOL'
  id_all <<- id_all  
  
    
} else if(dataformat == "Processed dataset [sp|up|entry_species, peptide, expr]"){

  incProgress(0.1, message = "User's defined dataset")
  
  #tmp <- read_excel(data_path, sheet = 'Area - peptides') 
  tmp <- read.csv(data_path, header = TRUE, stringsAsFactors = FALSE) 
  
  areaPept <- tmp
  colnames(areaPept) <- c('Protein', 'peptide', sample_label)
  
  areaPept[ , 3:length(areaPept)] <- lapply(areaPept[ , 3:length(areaPept)], as.numeric)
  
  areaPept <<- areaPept # use later for the number of peptide/protein
  
  tmp2 <- areaPept[, c(1, 3:length(areaPept))]
  areaProt <- tmp2 %>% dplyr::group_by(Protein) %>% dplyr::summarise_at(.vars = sample_label, .funs = sum) 
  
  
  df <- areaProt[ ,1] %>% 
    separate(Protein, c("sp", "uniProtID", "entry_name"), sep = "\\|") %>%
    separate(entry_name, c("entry_names", "species"), sep = "_") %>%
    dplyr::select(uniProtID, entry_names, species) 
  
  incProgress(0.1, message = "UniProt ID-Gene mapping")
  
  tmp <- getBM(attributes = c('uniprotswissprot', 'external_gene_name'), 
               filters = 'uniprotswissprot', 
               values = df$uniProtID, 
               # mart = ensembl) 
               mart = chose_species)
  
  colnames(tmp) <- c('uniProtID', "gene.SYMBOL")
  df <- left_join(df, tmp[!duplicated(tmp$uniProtID), ], by = "uniProtID") 
  
  # for a given missing gene name, replaced by entry_names 
  ind <- is.na(df$gene.SYMBOL)
  df$gene.SYMBOL[ind] <- df$entry_names[ind]
  
  id_all <<- df
  
  df2 <- areaProt[ , 2:length(areaProt)]
  df2[df2 == 0] <- NA
  colnames(df2) <- sample_label
  expr_all <<- df2
  
} else if(dataformat == "Processed dataset [up, peptide, expr]"){
  
  incProgress(0.1, message = "User's defined dataset")
  
  #tmp <- read_excel(data_path, sheet = 'Area - peptides') 
  tmp <- read.csv(data_path, header = TRUE, stringsAsFactors = FALSE) 
  
  areaPept <- tmp
  colnames(areaPept) <- c('uniProtID', 'peptide', sample_label)
  
  areaPept[ , 3:length(areaPept)] <- lapply(areaPept[ , 3:length(areaPept)], as.numeric)
  
  areaPept <<- areaPept # use later for the number of peptide/protein
  
  tmp2 <- areaPept[, c(1, 3:length(areaPept))]
  areaProt <- tmp2 %>% dplyr::group_by(uniProtID) %>% dplyr::summarise_at(.vars = sample_label, .funs = sum) 
  
  
  df <- areaProt[ ,1]
  
  incProgress(0.1, message = "UniProt ID-Gene mapping")
  
  tmp <- getBM(attributes = c('uniprotswissprot', 'external_gene_name'), 
               filters = 'uniprotswissprot', 
               values = df$uniProtID, 
               # mart = ensembl) 
               mart = chose_species)
 
  colnames(tmp) <- c('uniProtID', "gene.SYMBOL")
  df <- left_join(df, tmp[!duplicated(tmp$uniProtID), ], by = "uniProtID") 
  
  # for a given missing gene name, replaced by uniProtID 
  ind <- is.na(df$gene.SYMBOL)
  df$gene.SYMBOL[ind] <- df$uniProtID[ind]
  
  df$gene.SYMBOL <- make.names(df$gene.SYMBOL, unique = TRUE)
  
  id_all <<- df
  
  df2 <- areaProt[ , 2:length(areaProt)]
  df2[df2 == 0] <- NA
  colnames(df2) <- sample_label
  expr_all <<- df2
  
} else if(dataformat == "Processed dataset [up, expr]"){
  
  incProgress(0.1, message = "User's defined dataset")
  
  #tmp <- read_excel(data_path, sheet = 'Area - peptides') 
  tmp <- read.csv(data_path, header = TRUE, stringsAsFactors = FALSE) 
  
  areaProt <- tmp
  colnames(areaProt) <- c('uniProtID', sample_label)
  
  areaProt[ , 2:length(areaProt)] <- lapply(areaProt[ , 2:length(areaProt)], as.numeric)
  
  areaProt <<- areaProt # use later for the number of peptide/protein
  
  df <- data.frame(uniProtID = areaProt[ ,1])
  
  incProgress(0.1, message = "UniProt ID-Gene mapping")
  
  tmp <- getBM(attributes = c('uniprotswissprot', 'external_gene_name'), 
               filters = 'uniprotswissprot', 
               values = df$uniProtID, 
               # mart = ensembl) 
               mart = chose_species)
  
  colnames(tmp) <- c('uniProtID', "gene.SYMBOL")
  df <- left_join(df, tmp[!duplicated(tmp$uniProtID), ], by = "uniProtID") 
  
  # for a given missing gene name, replaced by uniProtID 
  ind <- is.na(df$gene.SYMBOL)
  df$gene.SYMBOL[ind] <- df$uniProtID[ind]
  
  id_all <<- df
  
  df2 <- areaProt[ , 2:length(areaProt)]
  df2[df2 == 0] <- NA
  colnames(df2) <- sample_label
  expr_all <<- df2

} 

}
  