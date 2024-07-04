
rm(list=ls())
gc()

devtools::load_all()

#---------------------------------------------------
dir_data="data"
eqtl_list <- "meta_list.txt"

#-------gwas data (outcome data) ----------------------------------

#结局数据：可以是id或本地数据，但如果是vcf.gz格式，还需要下载其index文件（tbi格式）
outcomeID="finngen_R10_NAFLD.gz"  #结局数据
type2 = "cc"

#---如果是ieu网站数据（id或本地文件均可）不需要填下面信息；其他来源的数据要填
pval_gwas="pval"
eaf_gwas="af_alt"
beta_gwas="beta"
se_gwas="sebeta"
SNP_gwas="rsids"
chr_gwas="#chrom"
pos_gwas="pos"

samplesize_gwas=412181  #如果表里有，直接用；如果没有，则需要在下载数据时，查看相应信息


#----packages------------------------------------
inst_packages()

#----processing the outcome data----------------
print("processing the gwas data, it will take a long time.......")
gwas_df <- fread(outcomeID)

gwas_df <- gwas_df %>% dplyr::rename(pval.outcome=pval_gwas,
                                     eaf.outcome=eaf_gwas,
                                     beta.outcome=beta_gwas,
                                     se.outcome=se_gwas,
                                     SNP=SNP_gwas,
                                     chr.outcome=chr_gwas,
                                     pos.outcome=pos_gwas)
head(gwas_df)

if(class(samplesize_gwas)=="character")
  gwas_df <- gwas_df %>% dplyr::rename(samplesize.outcome=samplesize_gwas)else
    gwas_df$samplesize.outcome <- samplesize_gwas #添加例数

gwas_df <- gwas_df[!is.na(gwas_df$beta.outcome),]%>%
  .[!duplicated(.$SNP),]


#---metabolites list------------------------------------------
meta_list <- read.table(eqtl_list,header = T)

#------foreach analysis-------------------------------------

result_table <- data.frame()

foreach(x=meta_list$id[1] ,.errorhandling = "pass") %do% {

meta_file=paste0(dir_data,"/",x,".txt")

eqtl_df <- read_table(meta_file) %>% dplyr::arrange(pval.exposure)

#---eqtl共定位区域-------------
geneChr <- eqtl_df[i,1]
geneStart <- top10_snp[i,2]
geneEnd <- top10_snp[i,2]
chrpos <- paste0(eqtl[i,1],":",top10_snp[i,2]-100000,"-",top10_snp[i,2]+100000)
#---eqtl data 取共定位区域----
eqtl_coloc <- eqtl_df %>%
  subset(chr.exposure==geneChr & pos.exposure>geneStart-number & pos.exposure<geneEnd+number)
#rm(eqtl_df)
#write.csv (eqtl_coloc, paste0("analysis results/",str_extract(eqtlID,".*?(?=\\.)"),"_coloc_data.csv"),row.names=F )
#print("Processing the eQTL data is done.")


#---------------------------
#---outcome data 取共定位区域----
gwas_coloc <- gwas_df %>%
  subset(chr.outcome==geneChr & pos.outcome>geneStart-number & pos.outcome<geneEnd+number)
#rm(gwas_df)
#write.csv (gwas_coloc, paste0("analysis results/",str_extract(outcomeID,".*?(?=\\.)"),"_coloc_data.csv"),row.names=F )
#print("Processing the gwas data is done.")


#---整合eqtl和gwas数据-------
merg_data <- merge(eqtl_coloc,gwas_coloc,by="SNP")
if(nrow(merg_data)==0)
  if(nontarget)
    stop("There is no any shared SNP between the eqtl and gwas data.") else
      stop("There is no any shared SNP between the eqtl and gwas data. Thus, the target gene may be unsuitable.")

dataset1 <- list(pvalues=merg_data$pval.exposure,
                 N=merg_data$samplesize.exposure,
                 MAF= ifelse(merg_data$eaf.exposure>0.5,1-merg_data$eaf.exposure,merg_data$eaf.exposure),
                 beta=merg_data$beta.exposure,
                 varbeta= merg_data$se.exposure^2,
                 type=type1,
                 snp=merg_data$SNP,
                 z= merg_data$beta.exposure/merg_data$se.exposure,
                 chr=merg_data$chr.exposure,
                 pos=merg_data$pos.exposure,
                 id=eqtlID)


dataset2 <- list(pvalues=merg_data$pval.outcome,
                 N=merg_data$samplesize.outcome,
                 MAF= ifelse(merg_data$eaf.outcome>0.5,1-merg_data$eaf.outcome,merg_data$eaf.outcome),
                 beta=merg_data$beta.outcome,
                 varbeta= merg_data$se.outcome^2,
                 type=type2,
                 snp=merg_data$SNP,
                 z= merg_data$beta.outcome/merg_data$se.outcome,
                 chr=merg_data$chr.outcome,
                 pos=merg_data$pos.outcome,
                 id=outcomeID)

#------共定位分析--------------------------------------------------------
chr_pos <- paste0(i,": ","chr.exposure:",top10_snp[i,1]," - ","SNP:",top10_snp[i,5]," - ",top10_snp[i,6])

message (chr_pos)

col_result <- coloc::coloc.abf(dataset1, dataset2)

col_data <- col_result$results %>% dplyr::arrange(-SNP.PP.H4)  # -SNP.PP.H4降序

#筛选共定位的位点
#col_snp <- col_data %>% filter(SNP.PP.H4 > SNP_PP_H4)
col_r <- col_data[1,]
col_r$chr_pos <- chr_pos

result=rbind(result,col_r)






#----foreach end-------
}



chr_pos <- function(idchr_pos <- function(dir_file,eqtlID,type1,outcomeID,type2,pval_gwas,eaf_gwas,
         beta_gwas,se_gwas,SNP_gwas,chr_gwasm,pos_gwas,
         samplesize_gwas,geneChr,geneEnd,number,SNP_PP_H4,nontarget,n){
#-----start-------
if(nontarget)
  stop("nontarget is FALSE,and thus the chrpos test cann't be analyzed.")

inst_packages()

#---create dir----------------
if(nontarget)
  dir_file="analysis results (nontarget)" else
    dir_file="analysis results (target)"

  if(!dir.exists(dir_file))
    dir.create(dir_file)

print("processing the eqtl data, it will take a long time.......")
#---ieu id data---
if(!grepl("[.]",eqtlID,ignore.case = T)){
    eqtl_df <- extract_instruments(outcomes=eqtlID,clump=F)
    eqtl_df <- eqtl_df[!is.na(eqtl_df$beta.exposure),] %>%
      .[!duplicated(.$SNP),]} else{ # else 01

#---ieu vcf data---
if(grepl("vcf.gz",eqtlID,ignore.case = T)){
  eqtl_df <- readVcf(eqtlID) %>%
  gwasvcf_to_TwoSampleMR(type = "exposure")
  eqtl_df <- eqtl_df[!is.na(eqtl_df$beta.exposure),] %>%
  .[!duplicated(.$SNP),]} else{  #else 02

#---other format data---
  file_type <- str_extract(eqtlID,"(?<=\\.)[^\\.]+$")
  if(grepl(file_type,"xlsx",ignore.case = T) | grepl(file_type,"csv",ignore.case = T))
  eqtl_df <- eval(str2expression(paste0("read.",file_type,"(eqtlID)"))) else
  eqtl_df <- fread(eqtlID)

  eqtl_df <- eqtl_df %>% dplyr::rename(pval.exposure=pval_eqtl,
                                      eaf.exposure=eaf_eqtl,
                                      beta.exposure=beta_eqtl,
                                      se.exposure=se_eqtl,
                                      SNP=SNP_eqtl,
                                      chr.exposure=chr_eqtl,
                                      pos.exposure=pos_eqtl)

 if(class(samplesize_eqtl)=="character")
  eqtl_df <- eqtl_df %>% dplyr::rename(samplesize.exposure=samplesize_eqtl)else
  eqtl_df$samplesize.exposure <- samplesize_eqtl #添加例数

  eqtl_df <- eqtl_df[!is.na(eqtl_df$beta.exposure),]%>%
            .[!duplicated(.$SNP),]

  } #else 02 end

  }# else 01 end


#------gwas data----------------------------------
  print("processing the gwas data, it will take a long time.......")
  gwas_df <- fread(outcomeID)

  gwas_df <- gwas_df %>% dplyr::rename(pval.outcome=pval_gwas,
                                       eaf.outcome=eaf_gwas,
                                       beta.outcome=beta_gwas,
                                       se.outcome=se_gwas,
                                       SNP=SNP_gwas,
                                       chr.outcome=chr_gwas,
                                       pos.outcome=pos_gwas)
  head(gwas_df)

  if(class(samplesize_gwas)=="character")
    gwas_df <- gwas_df %>% dplyr::rename(samplesize.outcome=samplesize_gwas)else
      gwas_df$samplesize.outcome <- samplesize_gwas #添加例数

  gwas_df <- gwas_df[!is.na(gwas_df$beta.outcome),]%>%
    .[!duplicated(.$SNP),]


#-----------------------------------------------------------------#
top10_snp <- eqtl_df %>%
             dplyr::arrange(pval.exposure) %>% .[1:n,]

write.xlsx(top10_snp,paste0(dir_file,"/eqtl_chr_pos.xlsx"))


#------foreach ------------------
result=data.frame()

foreach (i=1:n,.errorhandling = "pass") %do% {

#---eqtl共定位区域-------------
chrpos <- paste0(top10_snp[i,1],":",top10_snp[i,2]-100000,"-",top10_snp[i,2]+100000)
geneChr <- top10_snp[i,1]
geneStart <- top10_snp[i,2]
geneEnd <- top10_snp[i,2]

#---eqtl data 取共定位区域----
eqtl_coloc <- eqtl_df %>%
  subset(chr.exposure==geneChr & pos.exposure>geneStart-number & pos.exposure<geneEnd+number)
#rm(eqtl_df)
#write.csv (eqtl_coloc, paste0("analysis results/",str_extract(eqtlID,".*?(?=\\.)"),"_coloc_data.csv"),row.names=F )
#print("Processing the eQTL data is done.")


#---------------------------
#---outcome data 取共定位区域----
gwas_coloc <- gwas_df %>%
  subset(chr.outcome==geneChr & pos.outcome>geneStart-number & pos.outcome<geneEnd+number)
#rm(gwas_df)
#write.csv (gwas_coloc, paste0("analysis results/",str_extract(outcomeID,".*?(?=\\.)"),"_coloc_data.csv"),row.names=F )
#print("Processing the gwas data is done.")


#---整合eqtl和gwas数据-------
merg_data <- merge(eqtl_coloc,gwas_coloc,by="SNP")
if(nrow(merg_data)==0)
  if(nontarget)
    stop("There is no any shared SNP between the eqtl and gwas data.") else
      stop("There is no any shared SNP between the eqtl and gwas data. Thus, the target gene may be unsuitable.")

dataset1 <- list(pvalues=merg_data$pval.exposure,
                 N=merg_data$samplesize.exposure,
                 MAF= ifelse(merg_data$eaf.exposure>0.5,1-merg_data$eaf.exposure,merg_data$eaf.exposure),
                 beta=merg_data$beta.exposure,
                 varbeta= merg_data$se.exposure^2,
                 type=type1,
                 snp=merg_data$SNP,
                 z= merg_data$beta.exposure/merg_data$se.exposure,
                 chr=merg_data$chr.exposure,
                 pos=merg_data$pos.exposure,
                 id=eqtlID)


dataset2 <- list(pvalues=merg_data$pval.outcome,
                 N=merg_data$samplesize.outcome,
                 MAF= ifelse(merg_data$eaf.outcome>0.5,1-merg_data$eaf.outcome,merg_data$eaf.outcome),
                 beta=merg_data$beta.outcome,
                 varbeta= merg_data$se.outcome^2,
                 type=type2,
                 snp=merg_data$SNP,
                 z= merg_data$beta.outcome/merg_data$se.outcome,
                 chr=merg_data$chr.outcome,
                 pos=merg_data$pos.outcome,
                 id=outcomeID)

#------共定位分析--------------------------------------------------------
chr_pos <- paste0(i,": ","chr.exposure:",top10_snp[i,1]," - ","SNP:",top10_snp[i,5]," - ",top10_snp[i,6])

message (chr_pos)

col_result <- coloc::coloc.abf(dataset1, dataset2)

col_data <- col_result$results %>% dplyr::arrange(-SNP.PP.H4)  # -SNP.PP.H4降序

#筛选共定位的位点
#col_snp <- col_data %>% filter(SNP.PP.H4 > SNP_PP_H4)
col_r <- col_data[1,]
col_r$chr_pos <- chr_pos

result=rbind(result,col_r)

}

write.xlsx(result,paste0(dir_file,"/chr_pos_coloc.xlsx"))

#-------the end-----------
}


chrpos_test <- function(n=20){

chr_pos (dir_file,eqtlID,type1,outcomeID,type2,pval_gwas,eaf_gwas,
         beta_gwas,se_gwas,SNP_gwas,chr_gwasm,pos_gwas,
         samplesize_gwas,geneChr,geneEnd,number,SNP_PP_H4,nontarget,n)

}




