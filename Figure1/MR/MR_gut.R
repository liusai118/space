library(MRInstruments)
library(MendelianRandomization) 
library(TwoSampleMR) 
library(simex) 
library(data.table)
library(dplyr)
library(VariantAnnotation)
library(VariantAnnotation)
library(gwasglue)
################################
res_all <- as.data.frame(matrix(ncol = 9))
colnames(res_all) <- c("id.exposure","id.outcome","outcome","exposure","method","nsnp","b","se","pval")
file_list_all <- c("phylum","order","class","family","genus")
library(VariantAnnotation)
vcfRT <- readVcf("../ebi-a-GCST006697.vcf.gz")
outcome = gwasvcf_to_TwoSampleMR(vcf = vcfRT, type="outcome")
for ( j in 1:length(file_list_all)){
filelist <- list.files(paste0("gut",file_list_all[j]))
for (i in 1:length(filelist)){
  a<- fread(paste0(filelist[i],"/"))
  b <- dplyr::rename(a, c(
    "chr.exposure"="chr",
    "pos.exposure"="bp",
    "SNP"="rsID",
    "effect_allele.exposure"="eff.allele",
    "other_allele.exposure"="ref.allele",
    "beta.exposure"="beta",
    "se.exposure"="SE",
    "pval.exposure"="P.weightedSumZ",
    "N"="N"))
  b <- b[,c(2,3,4,5,6,7,8,10,11)]
  b <- b %>%
    dplyr::filter(pval.exposure < 1e-5)
  exp_dat <- format_data(
    as.data.frame(b),
    type='exposure',
    snp_col = "SNP",
    beta_col = "beta.exposure",
    se_col = "se.exposure",
    effect_allele_col ="effect_allele.exposure",
    other_allele_col = "other_allele.exposure",
    #eaf_col = "effect_allele_freq",
    pval_col = "pval.exposure"
  )
  exp_dat$eaf.exposure <- paste0(filelist[i])
  exp_dat <- exp_dat%>%
    dplyr::rename(
      rsid = SNP,        
      pval = pval.exposure,       
      id = eaf.exposure    
    )
  exposure_data_clumped <- ld_clump(
    dat = exp_dat,                 
    clump_kb = 1000,             
    clump_r2 = 0.001,           
    bfile = "/home/liusai/index/GWAS/EUR",       
    plink_bin = "/home/liusai/miniconda3/envs/GWAS/bin/plink"     
  )
  
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
  snp_list <- exposure_data_clumped$rsid
  
  outcome_data_filtered <- outcome[outcome$SNP %in% snp_list, ]
  if (nrow(outcome_data_filtered) != 0){
    
    exposure_data_clumped <- rename(exposure_data_clumped, c(
      "rsid"="SNP"
    ))
    dat <- harmonise_data(exposure_dat = exp_dat, outcome_dat = outcome)
    res <- mr(dat)
    if(nrow(res) != 0){
      res$id.exposure <- paste0(filelist[i])
      res_all <- rbind(res_all,res)
    }
    else{
      next}} 
  else {
    next
  }}
write.csv(res_all,paste0("res_gut.csv"))}
