library(MRInstruments)
library(MendelianRandomization) 
library(TwoSampleMR) 
library(simex) 
library(data.table)
library(dplyr)
library(VariantAnnotation)
library(VariantAnnotation)
library(gwasglue)
##################################
library(VariantAnnotation)
vcfRT <- readVcf("../ebi-a-GCST006697.vcf.gz")
outcome = gwasvcf_to_TwoSampleMR(vcf = vcfRT, type="outcome")
res_all <- as.data.frame(matrix(ncol = 9))
colnames(res_all) <- c("id.exposure","id.outcome","outcome","exposure","method","nsnp","b","se","pval")
filelist <- filelist("skin/")
for (i in 1:length(filelist)){
  a<- fread(filelist[i])
  a<- dplyr::rename(a, c(
    "chr.exposure"="chromosome",
    "pos.exposure"="base_pair_location",
    "SNP"="rsid",
    "effect_allele.exposure"="effect_allele",
    "other_allele.exposure"="other_allele",
    "beta.exposure"="beta_site1",
    "se.exposure"="standard_error_site1",
    "pval.exposure"="p_value_site1",
    "N"="n_site1"))
  id <- paste0(a$microbial_feature_level[1],"_",a$microbial_feature[1])
  a <- a[,c(1,2,3,4,19,20,21,22,36)]
  a <-a %>%
    dplyr::filter(pval.exposure < 1e-5)
  
  a$eaf.exposure <- id
  a <- a %>%
    dplyr::rename(
      rsid = SNP,        
      pval = pval.exposure,          
      id = eaf.exposure      
    )
  exposure_data_clumped <- ld_clump(
    dat = a,                     
    clump_kb = 1000,               
    clump_r2 = 0.001,              
    clump_p = 1,pop = "EAS",                    
    bfile = "/home/liusai/index/GWAS/EAS",         
    plink_bin = "/home/liusai/miniconda3/envs/GWAS/bin/plink"     
  )
  snp_list <- exposure_data_clumped$rsid
  exp_dat <- format_data(
    as.data.frame(exposure_data_clumped),
    type='exposure',
    snp_col = "rsid",
    beta_col = "beta.exposure",
    se_col = "se.exposure",
    effect_allele_col ="effect_allele.exposure",
    other_allele_col = "other_allele.exposure",
    #eaf_col = "effect_allele_freq",
    pval_col = "pval"
  )
  # 筛选结局数据
  outcome_data_filtered <- outcome[outcome$SNP %in% snp_list, ]
  if (nrow(outcome_data_filtered) != 0){
    exposure_data_clumped <- rename(exposure_data_clumped, c(
      "rsid"="SNP"
    ))
    dat <- harmonise_data(exposure_dat = exp_dat, outcome_dat = outcome)
    res <- mr(dat)
    res$id.exposure <- paste0(filelist[i])
    res_all <- rbind(res_all,res)}else {
      next
    }
}
write.csv(res_all,"res_all_oral.csv")
