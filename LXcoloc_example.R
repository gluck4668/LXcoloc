
install.packages("devtools")
library(devtools)

install_github("gluck4668/LXcoloc")

library(LXcoloc)


rm(list=ls())
gc()

#devtools::load_all()

#-----eqtl data (exposure data)-------------------------------------
{

#暴露数据：可以是id或本地数据，但如果是vcf.gz格式，还需要下载其index文件（tbi格式）
#eqtlID="THRB_cis_eQTLs_data.xlsx"
eqtlID="finngen_R10_F5_DEPRESSION_RECURRENT.gz"
type1 = "quant" #数据类型，一般有分类变量"cc"和连续变量"quant"

#---如果是ieu网站数据（id或本地文件均可）不需要填下面信息；其他来源的数据要填
pval_eqtl="pval"
eaf_eqtl="af_alt"
beta_eqtl="beta"
se_eqtl="sebeta"
SNP_eqtl="rsids"
chr_eqtl="#chrom"
pos_eqtl="pos"

samplesize_eqtl=265431 #如果表里有，直接用；如果没有，则需要在下载数据时，查看相应信息


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



#------target gene information------------------------------
target_gene="THRB"  # 用靶基因的染色体范围来筛选暴露因素snp，然后分析eqtl和gwas共定位

geneChr=3                # FASN 基因所在的染色体(https://www.ncbi.nlm.nih.gov/gene)
geneStart=24158644       #染色体起始位置
geneEnd=24537199         #染色体终止位置

number=100000             # 取值前后范围


#------non-target-----------------------------------------
nontarget=T #无靶基因，以暴露data中pvalue最小的snp做为假定“target”来分析两个数据间的共定位


#-----#共定位阈值-----------------------------------------
SNP_PP_H4=0.75  #共定位阈值

}

LXcoloc()

chrpos_test(n=50) # 在进行nontarget=T分析时，批量分析前n个的共定位情况


#============================================================================#

#---查找靶基因的eqtl ID--------------------------------
#devtools::install_github( "rmgpanw/eqtlgentools")
#library(eqtlgentools)
library(dplyr)
library(openxlsx)
library(data.table)

target_gene="THRB"  # 部分基因名，或全基因名

cis_eqtl_data <- fread("2019-12-11-cis-eQTLs (complete data).gz")

gene_eqtl <- cis_eqtl_data %>% filter(grepl(target_gene,GeneSymbol,ignore.case =T))

table(gene_eqtl$GeneSymbol) # 统计显示

write.xlsx(gene_eqtl,paste0(target_gene,"_cis_eQTLs_data.xlsx"))


