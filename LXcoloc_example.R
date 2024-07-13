
install.packages("devtools")
library(devtools)

install_github("gluck4668/LXcoloc")

library(LXcoloc)


rm(list=ls())
gc()

# devtools::load_all()

#-----eqtl data (exposure data)-------------------------------------
{

#暴露数据：可以是id或本地数据，但如果是vcf.gz格式，还需要下载其index文件（tbi格式）
#eqtlID="THRB_cis_eQTLs_data.xlsx"
eqtlID= "exposure"
type1 = "quant" #数据类型，一般有分类变量"cc"和连续变量"quant"

#---查询exposure数据参数
library(data.table)
# head(fread(dir(eqtlID,full.names = T)[1]))

#---如果是ieu网站数据（id或本地文件均可）不需要填下面信息；其他来源的数据要填
pval_eqtl="p_value"
eaf_eqtl="effect_allele_frequency"
beta_eqtl="beta"
se_eqtl="standard_error"
SNP_eqtl="variant_id"
chr_eqtl="chromosome"
pos_eqtl="base_pair_location"

samplesize_eqtl_file= "1400 meta gwas_data_info.xlsx" # "meta_GCST90091033_622_gwas_id_info.xlsx"
# NA # 如果数据来自GWAS Catalog，可以填NA；如果来自其他数据库，则需要提供含有ID和samplesize的excel表格


#-------gwas data (outcome data) ----------------------------------

#结局数据：可以是id或本地数据，但如果是vcf.gz格式，还需要下载其index文件（tbi格式）
outcomeID ="GCST90091033_buildGRCh37.tsv （done).gz"  #结局数据
type2 = "cc"

#---查询outcome数据参数
# head(fread(outcomeID))

#---如果是ieu网站数据（id或本地文件均可）不需要填下面信息；其他来源的数据要填
pval_gwas="p_value"
eaf_gwas="effect_allele_frequency"
beta_gwas="beta"
se_gwas="standard_error"
SNP_gwas="variant_id"
chr_gwas="chromosome"
pos_gwas="base_pair_location"

samplesize_gwas= "NAFLD gwas_info.xlsx"#"outcome_sample_info.xlsx_gwas_info.xlsx"
# 如果数据来自GWAS Catalog，可以填NA；如果来自其他数据库，则需要提供含有ID和samplesize的excel表格

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

LXcoloc()# chrpos_test(n=10) # 在进行nontarget=T分析时，批量分析前n个的共定位情况


#============================================================================#
library(data.table)
head(fread(outcomeID))

#-------------------------

# remotes::install_github('MRCIEU/TwoSampleMR')

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


