
#-------------------------------------------
bat_get_gwasinfo <- function(data_dir,data_list){

if(!is.na(data_list) & !is.na(data_dir))
  stop("Choose either 'data_list' or 'data_dir' at a time, not both.")


#---R packages--------
com_packages()

#------------------------
dir_file = "gwas catalog information"
if(!dir.exists(dir_file))
    dir.create(dir_file)
#------------------------
if(!is.na(data_list)) {
  id_list <- read.xlsx(data_list) %>% data.frame()
  }

#------------------------
if(!is.na(data_dir)){
file_list <- dir(data_dir) %>%
            str_extract(".*?(?=\\.)")

if(length(file_list)==0)
  stop(paste0("The '",data_dir,"' is empty."))

file_table <- data.frame()

for(i in file_list){
  if(grepl("_",i))
    x=str_extract(i,".*?(?=_)") else
      x=i
  file_table <- rbind(file_table,cbind(id=x))
}

id_list <- file_table %>% unique()
}

#------------------------
id_list$Trait <- NA
nn <- ncol(id_list)

#------------------------
gwas_data_info <- data.frame()

#-------foreach analysis----------------------
foreach(x=c(1:nrow(id_list)), .errorhandling = "pass") %do% {

print(paste0("It is number ",x, " of ",nrow(id_list)))

gwas_info <- get_studies(study_id =id_list[x,1])

if(grepl("level",gwas_info@studies$reported_trait,ignore.case = T))
 info <- str_extract(gwas_info@studies$reported_trait,".*(?= )") else
   info <- gwas_info@studies$reported_trait

id_list[x,nn] <-info

#samplesize <- str_extract(gwas_info@studies$initial_sample_size,"\\d.*\\d")
samplesize <- gwas_info@ancestries$number_of_individuals
initial_samplesize <- gwas_info@studies$initial_sample_size
study_arry <- gwas_info@studies$study_design_comment

#pop <- str_extract(gwas_info@studies$initial_sample_size,"(?<=\\ ).*?(?=\\ )")
pop <- gwas_info@ancestral_groups$ancestral_group
nsnp=gwas_info@studies$snp_count
pubmed_id <- gwas_info@publications$pubmed_id
date <- gwas_info@publications$publication_date %>% year()

m_info <- data.frame(id=id_list[x,1],
                     trait=info,
                     study_arry,
                     samplesize=samplesize,
                     initial_samplesize,
                     nsnp=nsnp,
                     pop=pop,
                     pubmed_id,
                     date=date)

if(any(duplicated(m_info$id))){

 if(length(table(m_info$samplesize))>1)
    samples_num <- sum(m_info$samplesize) else
       samples_num <- m_info$samplesize[1]

   pop_all <- paste0(m_info$pop[!duplicated(m_info$pop)],collapse = ",")

   m_info <- m_info %>% distinct(id,.keep_all = T)

   m_info$samplesize <- samples_num
   m_info$pop <- pop_all
}

gwas_data_info <- rbind(gwas_data_info,m_info)

p = round(100*x/nrow(id_list),4)
print(paste0("......",p,"% was completed......"))
#-----foreach end-------------
}

if(grepl("\\\\",data_dir))
  data_save <- str_extract(data_dir,"(?<=\\\\)[^\\\\]+$") else
  {if(grepl("/",data_dir))
    data_save <- str_extract(data_dir,"(?<=\\/)[^\\/]+$") else
      data_save <- data_dir
  }

if(!is.na(data_dir)){
  write.xlsx(gwas_data_info,paste0(dir_file,"/",data_save,"_gwas_id_info.xlsx"))
} else{
write.xlsx(id_list, paste0(dir_file,"/",data_list,"_info.xlsx"))
write.xlsx(gwas_data_info,paste0(dir_file,"/",data_list,"_gwas_info.xlsx"))
}

print(paste0("The information of the gwas id can be found in the folder of '",dir_file,"' "))
#----the end----------
}


#bat_get_gwasinfo(data_dir,data_list)


