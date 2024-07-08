LXcoloc <- function(){
#==========start===================
#----R packages---------------
inst_packages()

#---create dir----------------
if(nontarget)
  dir_file="analysis results (nontarget)" else
    dir_file="analysis results (target)"

if(!dir.exists(dir_file))
  dir.create(dir_file)

#-----ieu gwas id data--------
if(!grepl("[.]",eqtlID,ignore.case = T) & !grepl("[.]",outcomeID,ignore.case = T)){
  ieu_id(dir_file,eqtlID,type1,outcomeID,type2,geneChr,geneEnd,number,SNP_PP_H4,nontarget)} else { # if01 start


#-----gwas vcf data-----------
if(grepl("vcf.gz",eqtlID,ignore.case = T) & grepl("vcf.gz",outcomeID,ignore.case = T)){
  ieu_vcf(dir_file,eqtlID,type1,outcomeID,type2,geneChr,geneEnd,number,SNP_PP_H4,nontarget)} else { # if02 start


  gwas_cata (dir_file,eqtlID,type1,outcomeID,type2,pval_gwas,eaf_gwas,
                 beta_gwas,se_gwas,SNP_gwas,chr_gwasm,pos_gwas,
                 samplesize_gwas,geneChr,geneEnd,number,SNP_PP_H4,nontarget)

  # other_format(dir_file,eqtlID,type1,outcomeID,type2,pval_gwas,eaf_gwas,
  #                   beta_gwas,se_gwas,SNP_gwas,chr_gwasm,pos_gwas,
  #                   samplesize_gwas,geneChr,geneEnd,number,SNP_PP_H4,nontarget)




      } # if02 end
  } # if01 end


#==========the end=================

}
