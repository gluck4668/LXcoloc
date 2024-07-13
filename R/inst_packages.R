#-------R packages ----------------------------------------------
inst_packages <- function(){

  packs <- c("BiocManager","devtools","remotes","openxlsx","dplyr","readr",
             "ggplot2","tidyr","stringr","data.table","rlang","coloc",
             "MendelianRandomization","VariantAnnotation","foreach")

  pack_uninst <- packs[!packs %in% installed.packages()[,1]]

  if(length(pack_uninst)>0){
    tryCatch(install.packages(pack_uninst),error=function(e){e})
    tryCatch(BiocManager::install(pack_uninst),error=function(e){e})}

  mr_packs <- c("ieugwasr","MRInstruments", "TwoSampleMR","gwasvcf",
                "gwasglue","plinkbinr","locuscomparer","geni.plots")

  if(!"ieugwasr" %in% installed.packages()[,1])
    install.packages("ieugwasr-master/ieugwasr_1.0.0.tar.gz",
                     repos = NULL, type = "source")

  if(!"MRInstruments" %in% installed.packages()[,1])
    remotes::install_github("MRCIEU/MRInstruments")

  if(!"TwoSampleMR" %in% installed.packages()[,1])
    remotes::install_github("MRCIEU/TwoSampleMR")

  if(!"gwasvcf" %in% installed.packages()[,1])
    remotes::install_github("mrcieu/gwasvcf")

  if(!"gwasglue" %in% installed.packages()[,1])
    remotes::install_github("mrcieu/gwasglue")

  if(!"plinkbinr" %in% installed.packages()[,1])
    devtools::install_github("explodecomputer/plinkbinr")

  if(!"locuscomparer" %in% installed.packages()[,1])
    devtools::install_github("boxiangliu/locuscomparer")

  if(!"geni.plots" %in% installed.packages()[,1])
    remotes::install_github("jrs95/geni.plots", build_vignettes = TRUE)

#---------------------------------
  for(i in c(packs,mr_packs)){
    library(i,character.only = T)}


}

