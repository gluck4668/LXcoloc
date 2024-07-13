

ieu_id <- function(dir_file,eqtlID,type1,outcomeID,type2,
                  geneChr,geneEnd,number,SNP_PP_H4,nontarget) {


#-----eqtl data------------
eqtl_df <-extract_instruments(outcomes=eqtlID,clump=F)

#---共定位区域-------------
if(nontarget){
  top_snp <- eqtl_df %>% dplyr::arrange(pval.exposure) %>% .[1,]
  chrpos <- paste0(top_snp$chr.exposure,":",top_snp$pos.exposure-100000,
                   "-",top_snp$pos.exposure+100000)} else {
  chrpos <- paste0(geneChr, ":", geneStart - number, "-", geneEnd + number)}

#提取共定位数据
col_df <- ieugwasr_to_coloc(id1=eqtlID, id2=outcomeID, chrompos=chrpos,
                              type1 =type1,type2 = type2 )

#共定位分析
col_result <- coloc::coloc.abf(col_df[[1]], col_df[[2]])

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

#------locuscompare可视化-------------------
fn1_gwas <- data.frame(rsid=col_df[[2]]$snp,pval=col_df[[2]]$pvalues) %>% dplyr::arrange(pval)
fn2_eqtl <- data.frame(rsid=col_df[[1]]$snp,pval=col_df[[1]]$pvalues) %>% dplyr::arrange(pval)

coloc_plot <- locuscompare(in_fn1 = fn1_gwas, in_fn2 = fn2_eqtl,
                           title1 = 'CAD GWAS', title2 = 'Coronary Artery eQTL')

coloc_plot

plot01_file=paste0(dir_file,"/coloc_plot_locuscompare (",eqtlID,"-",outcomeID,") ",Sys.Date(),".png")
ggsave(plot01_file,coloc_plot, width = 1000,height = 800,units = "px",dpi = 150)

#---------geni.plots可视化--------------------------------------------------
snp_list=col_df[[1]]$snp

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

set_01 <- data.frame(marker=col_df[[1]]$snp,pvalue_1=col_df[[1]]$pvalues) %>%
          subset(marker %in% corr_snp)

set_02 <- data.frame(marker=col_df[[2]]$snp,pvalue_2=col_df[[2]]$pvalues) %>%
       subset(marker %in% corr_snp)

mark=data.frame(marker = col_df[[1]]$snp,chr = col_df[[1]]$chr, pos = col_df[[1]]$pos)

assoc <- Reduce(function(x,y)inner_join(x,y,by="marker"),list(mark,set_01,set_02)) %>%
         data.frame() %>% dplyr::arrange(marker)

traits=c(col_df[[1]]$id,col_df[[2]]$id)

if(is.na(lead_snp))
  highlights= dd_snp$snp[1] else
   highlights = lead_snp

library(geni.plots)
plot_geni <- fig_region_stack(data = assoc, corr = corr_list,
                 traits = traits,highlights=highlights,
                 build = 37,title_center = TRUE)
plot_geni

plot02_file=paste0(dir_file,"/coloc_plot_geni.plots (",eqtlID,"-",outcomeID,") ",Sys.Date(),".png")
ggsave(plot02_file,plot_geni,width = 1000,height = 1000,units = "px",dpi = 150)

print(paste0("The analsyis is done, and the results can be found in the folder '",dir_file,"'."))

return(plot_geni)

#--------------the end----------------------------

}




