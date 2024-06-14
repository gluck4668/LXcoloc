
ieu_vcf <- function(dir_file,eqtlID,type1,outcomeID,type2,geneChr,
                         geneEnd,number,SNP_PP_H4,nontarget){


#---共定位区域-------------
eqtl_df <-VariantAnnotation::readVcf(eqtlID) %>%
          gwasglue::gwasvcf_to_TwoSampleMR(type = "exposure")

if(nontarget){
  top_snp <- eqtl_df %>% dplyr::arrange(pval.exposure) %>% .[1,]
  chrpos <- paste0(top_snp$chr.exposure,":",top_snp$pos.exposure-100000,
                     "-",top_snp$pos.exposure+100000)} else {
  chrpos <- paste0(geneChr, ":", geneStart - number, "-", geneEnd + number)}

#----gwasvcf_to_coloc---------------
vcf_coloc <- gwasvcf_to_coloc(vcf1 = eqtlID,vcf2 = outcomeID,chrompos = chrpos)

#共定位分析
df1 <- vcf_coloc[[1]] %>% data.frame() %>% na.omit()

dataset1 <- list(pvalues=df1$pvalues,N=df1$N,MAF=df1$MAF,beta=df1$beta,
                   varbeta=df1$varbeta,type=df1$type[1],snp=df1$snp,z=df1$z,
                   chr=df1$chr,pos=df1$pos,id=df1$id)

df2 <- vcf_coloc[[2]] %>% data.frame() %>% na.omit()

dataset2 <- list(pvalues=df2$pvalues,N=df2$N,MAF=df2$MAF,beta=df2$beta,
                   varbeta=df2$varbeta,type=df2$type[1],snp=df2$snp,z=df2$z,
                   chr=df2$chr,pos=df2$pos,id=df2$id)

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

file=paste0(dir_file,"/coloc_SNP (",eqtlID,"-",outcomeID,") ",Sys.Date(),".csv")
write.csv(col_snp, file=file, row.names=F)

#---------geni.plots可视化--------------------------------------------------
dd1 <- data.frame(snp=dataset1$snp,pval=dataset1$pvalues) %>%
  dplyr::arrange(pval)

if(nrow(dd1)>500)
  n_snp=500 else
    n_snp=nrow(dd1)

dd_snp <- dd1[c(1:n_snp),]

snp_list=dd_snp$snp

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
  highlights= dd_snp$snp[1] else
    highlights = lead_snp

library(geni.plots)
plot_geni <- fig_region_stack(data = assoc, corr = corr_list,
                              traits = traits,highlights=highlights,
                              build = 37,title_center = TRUE)
plot_geni

plot01_file=paste0(dir_file,"/coloc_plot_geni.plots (",eqtlID,"-",outcomeID,") ",Sys.Date(),".png")
ggsave(plot01_file,plot_geni,width = 1000,height = 1000,units = "px",dpi = 150)


#------locuscompare可视化-------------------
fn1_gwas <- data.frame(rsid=assoc$marker,pval=assoc$pvalue_2)
fn2_eqtl <- data.frame(rsid=assoc$marker,pval=assoc$pvalue_1)

coloc_plot <- locuscompare(in_fn1 = fn1_gwas, in_fn2 = fn2_eqtl,
                           title1 = 'CAD GWAS', title2 = 'Coronary Artery eQTL')

plot02_file=paste0(dir_file,"/coloc_plot_locuscompare (",eqtlID,"-",outcomeID,") ",Sys.Date(),".png")
ggsave(plot02_file,coloc_plot, width = 1000,height = 800,units = "px",dpi = 150)

print(paste0("The analsyis is done, and the results can be found in the folder '",dir_file,"'."))

coloc_plot

result <- list(plot=coloc_plot)

return(result)
#--------------the end----------------------------

}
