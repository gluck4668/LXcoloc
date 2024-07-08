
other_format <- function(dir_file,eqtlID,type1,outcomeID,type2,pval_gwas,eaf_gwas,
                        beta_gwas,se_gwas,SNP_gwas,chr_gwasm,pos_gwas,
                        samplesize_gwas,geneChr,geneEnd,number,SNP_PP_H4,nontarget){

#-----eqtl data-----------------------------
print("processing the eQTL data, it will take much time.......")

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


#---共定位区域-------------
if(nontarget){
 top_snp <- eqtl_df %>% dplyr::arrange(pval.exposure) %>% .[1,]
 chrpos <- paste0(top_snp$chr.exposure,":",top_snp$pos.exposure-100000,"-",top_snp$pos.exposure+100000)
 geneChr <- top_snp$chr.exposure
 geneStart <- top_snp$pos.exposure
 geneEnd <- top_snp$pos.exposure
  } else {
 chrpos <- paste0(geneChr, ":", geneStart - number, "-", geneEnd + number)}


#---eqtl data 取共定位区域----
eqtl_coloc <- eqtl_df %>%
subset(chr.exposure==geneChr & pos.exposure>geneStart-number & pos.exposure<geneEnd+number)
#rm(eqtl_df)
#write.csv (eqtl_coloc, paste0("analysis results/",str_extract(eqtlID,".*?(?=\\.)"),"_coloc_data.csv"),row.names=F )
print("Processing the eQTL data is done.")


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
#head(gwas_df)

if(class(samplesize_gwas)=="character")
  gwas_df <- gwas_df %>% dplyr::rename(samplesize.outcome=samplesize_gwas)else
    gwas_df$samplesize.outcome <- samplesize_gwas #添加例数

gwas_df <- gwas_df[!is.na(gwas_df$beta.outcome),]%>%
                  .[!duplicated(.$SNP),]

#---outcome data 取共定位区域----
gwas_coloc <- gwas_df %>%
  subset(chr.outcome==geneChr & pos.outcome>geneStart-number & pos.outcome<geneEnd+number)
#rm(gwas_df)
#write.csv (gwas_coloc, paste0("analysis results/",str_extract(outcomeID,".*?(?=\\.)"),"_coloc_data.csv"),row.names=F )
print("Processing the gwas data is done.")


#---整合eqtl和gwas数据-------
merg_data <- merge(eqtl_coloc,gwas_coloc,by="SNP")
if(nrow(merg_data)==0)
if(nontarget)
  stop("There is no any overlapping SNP between the eqtl and gwas data.") else
  stop("There is no any overlapping SNP between the eqtl and gwas data. Thus, the target gene may be unsuitable.")

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
col_result <- coloc::coloc.abf(dataset1, dataset2)
col_data <- col_result$results %>% dplyr::arrange(-SNP.PP.H4)  # -SNP.PP.H4降序

#筛选共定位的位点
col_snp <- col_data %>% filter(SNP.PP.H4 > SNP_PP_H4)
lead_snp <- col_snp$snp[1]

if(is.na(lead_snp))
  print(paste0("The 'SNP.PP.H4' is ",round(col_data$SNP.PP.H4[1],3),
               ", suggesting that there is no any colocalization snp.")) else
                 print(paste0("The colocalization snp may be '",lead_snp,
                              "' with a 'SNP.PP.H4' value of ",round(col_snp$SNP.PP.H4,4)))

if(nontarget)
  file=paste0(dir_file,"/coloc_SNP (",eqtlID,"-",outcomeID,") ",Sys.Date(),"(nontarget).csv") else
  file=paste0(dir_file,"/coloc_SNP (",eqtlID,"-",outcomeID,") ",Sys.Date(),".csv")
write.csv(col_data, file=file, row.names=F)

#------locuscompare可视化-------------------
fn1_gwas <- data.frame(rsid=merg_data$SNP,pval=merg_data$pval.outcome)
fn2_eqtl <- data.frame(rsid=merg_data$SNP,pval=merg_data$pval.exposure)

coloc_plot <- locuscompare(in_fn1 = fn1_gwas, in_fn2 = fn2_eqtl,
                           title1 = outcomeID, title2 = eqtlID)

coloc_plot

if(nontarget)
  plot01_file=paste0(dir_file,"/coloc_plot_locuscompare (",eqtlID,"-",outcomeID,") ",Sys.Date(),"(nontarget).png") else
    plot01_file=paste0(dir_file,"/coloc_plot_locuscompare (",eqtlID,"-",outcomeID,") ",Sys.Date(),".png")
ggsave(plot01_file,coloc_plot, width = 1000,height = 800,units = "px",dpi = 150)


#---------geni.plots可视化--------------------------------------------------
snp_list=merg_data$SNP
table(duplicated(snp_list))

print("Begin to calculate the correlation coefficient,it may take a long time.......")
corr_list <-ld_snplist(variants=snp_list, with_alleles = TRUE, pop = "EUR",
                       opengwas_jwt = get_opengwas_jwt()) %>%
                       data.frame()

if(ncol(corr_list)<2 | nrow(corr_list)<2)
  stop(paste0("There is no enough corr number between the eqtl and gwas data.
       Maybe the chrompos is inappropriate or the target gene is inappropriate."))

print("Calculating the correlation coefficient is done.")
corr_snp <- rownames(corr_list) %>% str_extract(.,".*?(?=_)")
colnames(corr_list) <- corr_snp
rownames(corr_list) <- corr_snp
corr_list$snp <- corr_snp
corr_list <- corr_list %>% dplyr::arrange(snp)
corr_list <- corr_list[,-ncol(corr_list)]

set_01 <- data.frame(marker=dataset1$snp,pvalue_1=dataset1$pvalues) %>%
  subset(marker %in% corr_snp)

set_02 <- data.frame(marker=dataset2$snp,pvalue_2=dataset2$pvalues) %>%
  subset(marker %in% corr_snp)

mark=data.frame(marker = dataset1$snp,chr = dataset1$chr, pos = dataset1$pos)

assoc <- Reduce(function(x,y)inner_join(x,y,by="marker"),list(mark,set_01,set_02)) %>%
  data.frame() %>% dplyr::arrange(marker)

traits=c(dataset1$id[1],dataset2$id[1])

if(is.na(lead_snp))
  highlights= col_data$snp[1] else
    highlights = lead_snp

library(geni.plots)
assoc$pos <-as.numeric(assoc$pos) %>% as.integer()

plot_geni <- fig_region_stack(data = assoc, corr = corr_list,
                              traits = traits,highlights=highlights,
                              build = 37,title_center = TRUE)
plot_geni

if(nontarget)
  plot02_file=paste0(dir_file,"/coloc_plot_geni.plots (",eqtlID,"-",outcomeID,") ",Sys.Date(),"(nontarget).png") else
  plot02_file=paste0(dir_file,"/coloc_plot_geni.plots (",eqtlID,"-",outcomeID,") ",Sys.Date(),".png")
ggsave(plot02_file,plot_geni,width = 1000,height = 1000,units = "px",dpi = 150)


return(plot_geni)



#----------the end---------------------------------------------
}
