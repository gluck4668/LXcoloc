#-------R packages ----------------------------------------------
com_packages <- function(){

  packs <- c("BiocManager","devtools","remotes","openxlsx","dplyr","readr","utils","meta",
             "ggplot2","tidyr","stringr","data.table","rlang","coloc","gwasrapidd",
             "MendelianRandomization","VariantAnnotation","foreach","grid","forestploter",
             "VennDiagram","ggvenn","RColorBrewer","purrr",
             "circlize","reshape2","ComplexHeatmap")

  pack_uninst <- packs[!packs %in% installed.packages()[,1]]

  if(length(pack_uninst)>0){
    tryCatch(install.packages(pack_uninst),error=function(e){e})
    tryCatch(BiocManager::install(pack_uninst),error=function(e){e})}

#---------------------------------
  for(i in c(packs)){
    library(i,character.only = T)}
}


