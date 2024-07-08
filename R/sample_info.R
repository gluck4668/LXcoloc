
# sample_information <- sample_info(eqtlID,samplesize_eqtl_file,outcomeID,samplesize_gwas)

sample_info <- function(eqtlID,samplesize_eqtl_file,outcomeID,samplesize_gwas){

#----R packages--------
# com_packages()

#---eqtl samplesize------------------------------------------------------------
x1 <- str_extract(dir(eqtlID),".*?(?=\\.)")

idx <- data.frame()
for(i in x1){
  if(grepl("_",i))
    x2 <- str_extract(i,".*?(?=_)") %>% data.frame() else
      x2 <- i %>% data.frame()

  idx <- rbind(idx,x2)
}

eqt_list = idx[,1]

#----------------------
if(!is.na(samplesize_eqtl_file)){

  if(file.exists(samplesize_eqtl_file)){
       sam_eqtl <-read.xlsx(samplesize_eqtl_file)

       if(!any(grepl("samples",names(sam_eqtl))))
       warning(paste0(" '",samplesize_eqtl_file,"' did not include the samplesize information. Please check it."))

       not_idx <- paste0(eqt_list[!eqt_list %in% sam_eqtl[,1]],collapse = ", ")
       if(!not_idx=="")
       warning (paste0("'",not_idx,"'"," was not in the samplesize file '",samplesize_eqtl_file,"', Please check it. "))
   } else
  warning(paste0(" The smaplesize file '",samplesize_eqtl_file, "' was not existed. "))

not_smaple_info <- any( !file.exists(samplesize_eqtl_file), !any(grepl("samples",names(sam_eqtl))), !not_idx=="" )

} else
not_smaple_info <-TRUE
#---is.na end--------

#--------------------
if(not_smaple_info){
  warning (paste0("The eQTL samplesize was not provided, and it was obtained from GWAS Catalog online. "))
  data_list=NA
  data_dir=eqtlID
  bat_get_gwasinfo(data_dir,data_list)
  sam_eqtl <-read.xlsx(paste0("gwas catalog information/",data_dir,"_gwas_id_info.xlsx"))
}

#--------------------
sam_eqtl <- sam_eqtl %>% data.frame()
sam_sit <- grep("sample",names(sam_eqtl),ignore.case = T)[1]
names(sam_eqtl)[sam_sit] <- "samplesize"
names(sam_eqtl)[1] <- "ID"
names(sam_eqtl)[2] <- "Trait"

#---------eqtl samplesize end -------------------------------------------------

#---outcome samplesize-----
out_01 <- str_extract(outcomeID,".*?(?=\\.)")
out_02 <- str_extract(out_01,".*?(?=_)")
if(all(is.na(out_02)))
  out_id = out_01 else
    out_id = out_02

cls <- class(samplesize_gwas)

if(cls=="numeric")
  sam_gwas <- samplesize_gwas else {

   if(!is.na(samplesize_gwas)){

      if(file.exists(samplesize_gwas)){
        sam_gwas <-read.xlsx(samplesize_gwas)

        if(!any(grepl("samples",names(sam_gwas))))
          warning(paste0(" '",samplesize_gwas,"' did not include the samplesize information. Please check it."))

        if(!out_id %in% sam_gwas[,1])
          warning (paste0("The outcome '",out_id,"'"," was not in the samplesize file '",samplesize_gwas,"', Please check it. "))

        no_sample <- any(!any(grepl("samples",names(sam_gwas))),!out_id %in% sam_gwas[,1])
        if(no_sample)
        no_out_smaple_info <- TRUE else
          no_out_smaple_info <- FALSE

        } else  #-----if(file.exists(samplesize_gwas))----
        {warning(paste0(" The smaplesize file '",samplesize_gwas, "' was not existed. "))
        no_out_smaple_info <- TRUE}

    }  else #----if(!is.na(samplesize_gwas))-----
    no_out_smaple_info <- TRUE

} #----- samplesize_gwas else -------

if(exists("no_out_smaple_info")){

  if(no_out_smaple_info){
    warning (paste0("The outcome samplesize was not provided, and it was obtained from GWAS Catalog online. "))
    out_id_file <- data.frame(id=out_id)
    write.xlsx(out_id_file,"outcome_sample_info.xlsx")
    data_dir=NA
    data_list="outcome_sample_info.xlsx"
    bat_get_gwasinfo(data_dir=NA,data_list)
    sam_gwas <-read.xlsx("gwas catalog information/outcome_sample_info.xlsx_gwas_info.xlsx")
  }

  sam_gwas <- sam_gwas %>% data.frame()
  sam_sit <- grep("sample",names(sam_gwas),ignore.case = T)[1]
  names(sam_gwas)[sam_sit] <- "samplesize"
  names(sam_gwas)[1] <- "ID"
  names(sam_gwas)[2] <- "Trait"

} #-----if(exists("no_out_smaple_info"))----

  sample_info_list <- list(eqtl_sample=sam_eqtl,
                           outcome_sample=sam_gwas)

  return(sample_info_list)

}
