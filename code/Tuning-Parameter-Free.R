library(optparse)
library(parallel)
library(readr)
library(bigreadr)
library(bigsnpr)
library(data.table)
library(dplyr)
library(scales)
library(stringr) # for str_split


options(stringsAsFactors=F)
option_list = list(
  make_option("--homedir", action = "store", default = NA, type = "character",
              help="Path to save the output folder [Required]"),
  make_option("--PennPRS_path", action = "store", default = NA, type = "character",
              help="Path to the PennPRS folder [Required]"),
  make_option("--input_GWAS_path", action = "store", default = NA, type = "character",
              help="gwas path after qc"),
  # make_option("--userID", action = "store", default = NA, type = "character",
  #             help="User account ID [Required]"),
  make_option("--submissionID", action = "store", default = NA, type = "character",
              help="Job ID [Required]"),
  make_option("--methods", action = "store", default = 'LDpred2-auto', type = "character",
              help="Options: a subset of methods from LDpred2-auto and DBSLMM, divided by comma"),
  make_option("--trait", action = "store", default = NA, type = "character",
              help="trait name [Optional]"),
  make_option("--race", action = "store", default = NA, type = "character",
              help="Race of the training GWAS individuals. Options: EUR (European), AFR (African), 
              AMR (Mixed American, Hispanic/Latio), EAS (East Asian), or SAS (South Asian) [Required]"),
  make_option("--LDrefpanel", action = "store", default = '1kg', type = "character",
              help="LD reference panel. Options: '1kg' (1000 Genomes Project Phase 3) or 'ukbb' (UK Biobank) [Optional]"),

  # -------------------- Parameters in LDpred2-auto -------------------------
  make_option("--coef_shrink", action = "store", default = 1, type = "numeric",
              help="Shrinkage multiplicative coefficient to apply to off-diagonal elements of the correlation matrix.
              Reduce this up to 0.4 if there is (severe) mismatch between GWAS and the LD ref data. 
              Options: any value in [0.4,1] [Optional]"),
  make_option("--allow_jump_sign", action = "store", default = TRUE, type = "logical",
              help="Whether to allow for effects sizes to change sign in consecutive iterations. TRUE: normal sampling. 
              FALSE: to force effects to go through 0 first before changing sign. Setting this parameter to FALSE 
              could be useful to prevent instability (oscillation and ultimately divergence) of the Gibbs sampler. 
              This would also be useful for accelerating convergence of chains with a large initial value for p.
              Options: TRUE or FALSE [Optional]"),
  make_option("--use_MLE", action = "store", default = TRUE, type = "logical",
              help="Whether to use maximum likelihood estimation (MLE) to estimate alpha and the variance component (since v1.11.4), 
              or assume that alpha is -1 and estimate the variance of (scaled) effects as h2/(m*p), 
              as it was done in earlier versions of LDpred2-auto (e.g. in v1.10.8). 
              Default is TRUE, which should provide a better model fit, but might also be less robust.
              Options: TRUE or FALSE [Optional]"),
  make_option("--ensemble", action = "store", default = F, type = "logical",
              help="Whether to train a weighted combination of the single PRS models.
              Options: F [default: %default]"),
  make_option("--verbose", action="store", default=1, type="integer",
              help="Print logfile? 0 = no; 1 = yes [default: %default]"),
  make_option("--NCORES", action = "store", default = '5', type = "numeric",
              help="Number of cores used for parallel computing of PRS-CS across chromosomes.
              Default: 5.
              Options: positive integer [Optional]")
)

opt = parse_args(OptionParser(option_list=option_list))
print(opt)

# Input: ---------
PennPRS_path = opt$PennPRS_path
homedir = opt$homedir
input_GWAS_path = opt$input_GWAS_path
# userID = opt$userID
submissionID = opt$submissionID
methods = str_split(opt$methods,",")[[1]]
trait = opt$trait
race = opt$race
LDrefpanel = opt$LDrefpanel
NCORES = opt$NCORES

ld_path <- paste0(PennPRS_path, '/LD/', race, '/')
PUMAS_path = paste0(PennPRS_path,'/code/')
PRScs_path = paste0(PennPRS_path, 'software/PRScs/')
plink_path = paste0(PennPRS_path, 'software/')

trait_name = paste0(race,'_',trait)
ld_path0 <- paste0(ld_path, 'LD_1kg/') # set to the /LD_1kg folder under /LD/
if (LDrefpanel == '1kg'){
  eval_ld_ref_path <- paste0(ld_path, '/1KGref_plinkfile/') # set to the /1KGref_plinkfile folder under /LD/
  path_precalLD <- paste0(ld_path, '/LDpred2_lassosum2_corr_1kg/') # set to the /LDpred2_lassosum2_corr_1kg folder under /LD/
} 
# Job name/ID: e.g., trait_race_method_userID_submissionID
jobID = paste(c(trait,race, paste0(methods,collapse = '.'), submissionID), collapse = '_')
# Create a job-specific (trait, race, methods, userID, jobID) directory to save all the outputs, set the working directory to this directory
workdir = paste0(homedir,jobID,'/')
suppressWarnings(dir.create(workdir))
setwd(workdir) 


source(paste0(PUMAS_path, 'PennPRS_functions.R')) # please save the PennPRS_functions.R file to the /PUMAS/code/ directory
gwas_path <- paste0(workdir, 'sumdata/')
output_path <- paste0(workdir, 'output/')
input_path <- paste0(workdir, 'input_for_eval/')
PennPRS_finalresults_path <- paste0(workdir, 'PennPRS_results/')
dir.create(gwas_path, showWarnings = F)
dir.create(output_path, showWarnings = F)
dir.create(input_path, showWarnings = F)
dir.create(PennPRS_finalresults_path, showWarnings = F)
# Create a separate directory 'PRS_model_training/' to store input for training PRS models
prsdir0 = paste0(workdir, 'PRS_model_training/')
if (!dir.exists(prsdir0)) dir.create(prsdir0)
output_path_eval = paste0(workdir, 'output_for_eval/')
dir.create(output_path_eval, showWarnings = F)
# Create a separate directory 'PRS_model_training/' to store input for training PRS models
for (method in methods){
  prsdir = paste0(prsdir0, method,'/')
  if (!dir.exists(prsdir)) dir.create(prsdir)
}

# Parameters
if ('LDpred2-auto' %in% methods){
  coef_shrink <- opt$coef_shrink
  allow_jump_sign = opt$allow_jump_sign
  use_MLE = opt$use_MLE
}

# copy the input GWAS summary data, {Ancestry}_{Trait}.txt, to the /sumdata/ folder
system(paste0('cp -r ',input_GWAS_path, trait_name,'.txt ', workdir, 'sumdata/'))



######## QC for GWAS Summary Data:
cat(paste0("\n********************************************"))
cat(paste0("\n**** Step 0: QC for the input GWAS data ****"))
cat(paste0("\n********************************************\n"))

sumraw = bigreadr::fread2(paste0(workdir, 'sumdata/', trait_name, '.txt'))
sumraw$BETA = as.numeric(sumraw$BETA)
sumraw$SE = as.numeric(sumraw$SE)
sumraw$MAF = as.numeric(sumraw$MAF)
sumraw$P = as.numeric(sumraw$P)
# 0. Are there any SNP that have reasonable z-score?
chi2_thr = 30
remaining.SNPs = which(abs(sumraw$BETA/sumraw$SE) < sqrt(chi2_thr))
if (length(remaining.SNPs) < 5){
  stop(paste0("[Terminated] Job is terminated because less than 5 SNPs have z-score < sqrt(30), suggesting issues with the input GWAS data."))
}

n.na = sum(!complete.cases(sumraw))
if (n.na > 0){
  sumraw = sumraw[complete.cases(sumraw), ]
  if (n.na == 1) print(paste0('* 1 SNP has missing GWAS summary-level information and is removed.'))
  if (n.na > 1) print(paste0('* ', n.na, ' SNPs have missing GWAS summary-level information and are removed.'))
}

# 1. Remove SNPs with problematic BETA
beta.thr = 1e3
rm.indx1 = which(abs(sumraw$BETA) > beta.thr)
if (length(rm.indx1) > 0){
  if (length(rm.indx1) == 1) print(paste0('* 1 SNP has problematic GWAS summary statistic with abs(BETA) > ', beta.thr, ' and is removed.'))
  if (length(rm.indx1) > 1) print(paste0('* ', length(rm.indx1), ' SNPs have problematic GWAS summary statistics with abs(BETA) > ', beta.thr, ' and are removed.'))
} 

# 2. Remove SNPs with problematic p-values
rm.indx2 = which( ((sumraw$P) > 1) | (sumraw$P < 0))
if (length(rm.indx2) > 0){
  if (length(rm.indx2) == 1) print(paste0('* 1 SNP has p-value > 1 or < 0 and is removed.'))
  if (length(rm.indx2) > 1) print(paste0('* ', length(rm.indx2), ' SNPs have p-value > 1 or < 0 and are removed.'))
} 

# 3. Remove SNPs with an effective sample size less than 0.67 times the 90th percentile of sample size.
rm.indx3 = numeric()
# N.90percentile = quantile(sumraw$N, 0.9)
# rm.indx3 = which(sumraw$N < N.90percentile)
# if (length(rm.indx3) > 0){
#   if (length(rm.indx3) == 1) print(paste0('* 1 SNP has an effective sample size less than 0.67 times the 90th percentile of the total sample size and is removed.'))
#   if (length(rm.indx3) > 1) print(paste0('* ', length(rm.indx3), ' SNPs have an effective sample size less than 0.67 times the 90th percentile of the total sample size and are removed.'))
# } 

# 4. Remove SNPs with extremely large effect sizes (z^2> 100) 
chi2.thr = 1e3
rm.indx4 = which((sumraw$BETA/sumraw$SE)^2 > chi2.thr)
if (length(rm.indx4) > 0){
  if (length(rm.indx4) == 1) print(paste0('* 1 SNP has an extremely large effect size  (z-score^2 > ', chi2.thr, ') and is removed.'))
  if (length(rm.indx4) > 1) print(paste0('* ', length(rm.indx4), ' SNPs have extremely large effect sizes  (z-score^2 > ', chi2.thr, ') and are removed.'))
} 

# 5. Remove SNPs with zero SE 
rm.indx5 = which(sumraw$SE == 0)
if (length(rm.indx5) > 0){
  if (length(rm.indx5) == 1) print(paste0('* 1 SNP has SE = 0 and is removed.'))
  if (length(rm.indx5) > 1) print(paste0('* ', length(rm.indx5), ' SNPs have SE = 0 and are removed.'))
} 
rm.indx = unique(c(rm.indx1, rm.indx2, rm.indx3, rm.indx4, rm.indx5))

if (length(rm.indx) > 0){
  sumraw = sumraw[-rm.indx, ]
  if (nrow(sumraw) == 0){
    stop(paste0("[Terminated] 0 SNPs remaining after QC. Job terminated.\n * Please check the quality of the input GWAS summary data and make sure the columns are in correct format."))
  }
  if (nrow(sumraw) > 0){
    write_delim(sumraw, file = paste0(workdir, 'sumdata/', trait_name, '.txt'), delim='\t')
    if (length(rm.indx) == 1) print(paste0('* 1 problematic SNP removed. QC step completed.'))
    if (length(rm.indx) > 1) print(paste0('* QC step completed. ', nrow(sumraw), ' SNPs remaining. ', length(rm.indx), ' problematic SNPs removed.'))
  }
}
if (length(rm.indx) == 0) print(paste0('* QC step completed. ', nrow(sumraw), ' SNPs remaining. No SNP was removed.'))



if (('LDpred2-auto' %in% methods) | ('DBSLMM' %in% methods)){
  method = 'LDpred2-auto'
  prsdir = paste0(prsdir0, method,'/')
  if ( opt$verbose >= 1 ){
    print(paste0('************************************************************'))
    print(paste0('****** Start training PRS model based on LDpred2-auto ******'))
    print(paste0('************************************************************'))
  }
  map_ldref <- readRDS(paste0(ld_path, '/map/map_',LDrefpanel,'_ldref.rds'))
  
  sumraw = bigreadr::fread2(paste0(workdir, 'sumdata/', trait_name, '.txt'))
  sumstats = sumraw[,c('CHR','SNP','A1','A2','BETA','SE','P','N', 'MAF')]; rm(sumraw)
  sumstats$P = as.numeric(sumstats$P)
  sumstats$BETA = as.numeric(sumstats$BETA)
  sumstats$SE = as.numeric(sumstats$SE)
  sumstats$N = as.numeric(sumstats$N)
  sumstats$MAF = as.numeric(sumstats$MAF)
  names(sumstats) <- c("chr", "rsid", "a1", "a0", "beta", "beta_se", "p", "n_eff", "a1_sumdata_af")
  
  info_snp <- snp_match(sumstats, map_ldref, strand_flip = T, join_by_pos = F) # important: for real data, strand_flip = T
  info_snp <- tidyr:: drop_na(tibble::as_tibble(info_snp))
  sd_ldref <- with(info_snp, sqrt(2 * a1_af * (1 - a1_af)))
  sd_ss <- with(info_snp, sqrt(2 * a1_sumdata_af * (1 - a1_sumdata_af)))
  is_bad <- sd_ss < (0.5 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0.1 | sd_ldref < 0.05
  df_beta <- info_snp[!is_bad, ]
  
  td = paste0(prsdir0, 'temporary_LDpred2_lassosum2')
  if (!dir.exists(td)) dir.create(td)
  setwd(td)
  tmp <- tempfile(tmpdir = td)
  
  for (chr in 1:22) {
    cat(chr, ".. ", sep = "")
    ## indices in 'df_beta'
    ind.chr <- which(df_beta$chr == chr)
    ## indices in 'map_ldref'
    ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
    ## indices in 'corr0'
    ind.chr3 <- match(ind.chr2, which(map_ldref$chr == chr))
    if (length(ind.chr3) > 0){
      # corr0
      corr0 <- readRDS(paste0(path_precalLD, '/LD_ref_chr', chr, '.rds'))[ind.chr3, ind.chr3]
      if (chr == 1) {
        ld <- Matrix::colSums(corr0^2)
        corr <- as_SFBM(corr0, tmp, compact = TRUE)
      } else {
        if (length(corr0) == 1) corr0 =  as(1, "sparseMatrix") # as(corr0, "sparseMatrix") # as.matrix(corr0, 1, 1)
        ld <- c(ld, Matrix::colSums(corr0^2))
        corr$add_columns(corr0, nrow(corr))
      }
      print(paste0('Complete calculating LD for CHR ', chr))
      rm(corr0)
    }
  }
  
  (ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2,
                                  sample_size = n_eff, blocks = NULL)))
  ldsc_h2_est <- abs(ldsc[["h2"]])
  cat(paste0('Heritability estimate based on LD score regression: ', signif(ldsc[["h2"]], 3)))
  if (ldsc[["h2"]] < 0) cat(paste0('Warning: negative hertability estimate based on LD score regression.'))
  if ('LDpred2-auto' %in% methods){
    set.seed(2024)  # to get the same result every time
    multi_auto <- snp_ldpred2_auto(
      corr, df_beta, h2_init = ldsc_h2_est,
      vec_p_init = seq_log(1e-4, 0.2, length.out = 30), ncores = NCORES,
      use_MLE = use_MLE,  # FALSE if you have convergence issues or when power is low (need v1.11.9)
      allow_jump_sign = allow_jump_sign, shrink_corr = coef_shrink)
    
    # str(multi_auto[[k]])
    # `range` should be between 0 and 2
    range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est)))
    keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE)))
    # To get the final effects / predictions, you should only use chains that pass this filtering:
    if (length(keep) > 0){
      beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
      beta_auto0 = data.frame(df_beta[,c('chr','rsid','a1','a0')], beta_auto)
      colnames(beta_auto0) = c('CHR','SNP','A1','A2', 'BETA')
      # save(multi_auto, file = paste0(temdir,"ldpred2-",trait,"-auto.RData"))
      write_delim(beta_auto0, file = paste0(workdir, trait_name,'.',method,'.PRS.txt'), delim='\t')
      if (opt$verbose >= 1 ){
        print(paste0('****************************************'))
        print(paste0('**** Complete training LDpred2-auto ****'))
        print(paste0('****************************************'))
      }
      rm(corr)
    }
    if (length(keep) == 0){
      print(paste0('[Warning] All 30 chains in ', method, ' were deemed bad chains, and the resulting PRS model may not have sufficient power. \nPotential explanations:\n 1. The trait is not heritable.\n 2. The GWAS have insufficient power (e.g., due to low sample size) to develop a predictive PRS.\n 3. Issues with the input GWAS summary data (e.g., problematic BETA or SE).\n 4. ', method, ' is not powerful for developing PRS for the trait, in which case other methods can be considered.'))
      beta_auto <- rowMeans(sapply(multi_auto, function(auto) auto$beta_est))
      beta_auto0 = data.frame(df_beta[,c('chr','rsid','a1','a0')], beta_auto)
      colnames(beta_auto0) = c('CHR','SNP','A1','A2', 'BETA')
      # save(multi_auto, file = paste0(temdir,"ldpred2-",trait,"-auto.RData"))
      write_delim(beta_auto0, file = paste0(workdir, trait_name,'.',method,'.PRS.txt'), delim='\t')
      if (opt$verbose >= 1 ){
        print(paste0('****************************************'))
        print(paste0('**** Complete training LDpred2-auto ****'))
        print(paste0('****************************************'))
      }
      rm(corr)
    }
  }
  system(paste0('rm -rf ', paste0(prsdir0, 'temporary_LDpred2_lassosum2_ite1')))
}



if ('DBSLMM' %in% methods){
  method = 'DBSLMM'
  prsdir = paste0(prsdir0, method,'/')
  dbslmm_path = paste0(PennPRS_path, 'software/DBSLMM/')
  suppressWarnings(dir.create(paste0(prsdir, 'output/')))
  race.dbslmm = ifelse(race %in% c('EUR', 'AFR', 'EAS'), race, 'EUR')
  if ( opt$verbose >= 1 ){
    print(paste0('************************************************************'))
    print(paste0('********* Start training PRS model based on DBSLMM *********'))
    print(paste0('************************************************************'))
  }
  # ---------------------------------------------------------------------------------
  # --------------------- Step 1: Train PRS models using DBSLMM ---------------------
  # ---------------------------------------------------------------------------------
  # --------------------- Step 1.1: Input preparation for DBSLMM --------------------

  gwasinput = paste0(workdir, 'sumdata/',trait_name,'.txt')
  if (!file.exists(gwasinput)) print(paste0('A valid GWAS summary data file is missing.'))
  if (file.exists(gwasinput)){
    suppressWarnings(dir.create(paste0(prsdir, 'summary_gemma/')))
    sumraw0 = bigreadr::fread2(gwasinput)
    # impute position info:
    ref = bigreadr::fread2(paste0(eval_ld_ref_path, LDrefpanel, '_hm3_', race.dbslmm, '_ref.bim'))[,c(2,4)]
    colnames(ref) = c('SNP', 'ps')
    sumraw0 = merge(sumraw0, ref, by = 'SNP')
    sumraw0$n_mis = max(sumraw0$N) - sumraw0$N
    sumraw0 = sumraw0[,c('CHR', 'SNP', 'ps', 'n_mis', 'N', 'A1', 'A2', 'MAF', 'BETA', 'SE', 'P')] # allele1 <-> A1: REF
    colnames(sumraw0) = c('chr', 'rs',  'ps',  'n_mis',   'n_obs',   'allele1', 'allele0', 'af',  'beta', 'se', 'p_wald')
    for (chr in 1:22){
      sumdat.file = paste0(prsdir, 'summary_gemma/chr', chr, '.assoc.txt')
      sumraw = sumraw0[sumraw0$chr == chr, ]
      write_delim(sumraw, sumdat.file, delim = '\t', col_names = F)
      rm(sumraw)
    }
    print(paste0('Generating input GWAS data for ', method, ': completed.'))
  }
  
  # --------------------- Step 1.2: Run DBSLMM ---------------------
  system(paste0('chmod 777 ', dbslmm_path, 'dbslmm'))
  for (chr in 1:22) {
    # cat(paste0('Rscript ${DBSLMM} --summary ${summf}${chr}.assoc.txt --outPath ${outPath} --plink ${plink} 
    #            --dbslmm ${dbslmm} --ref ${ref}${chr} --model ${model} --n ${n} --nsnp ${m} --block ${blockf}${chr}.bed 
    #            --h2 ${herit} --thread ', NCORES,'\n\n'), file = zz)
    summf = paste0(prsdir, 'summary_gemma/chr')
    outPath = paste0(prsdir, 'output/')
    tem = bigreadr::fread2(paste0(prsdir, 'summary_gemma/chr',chr,'.assoc.txt'))
    n = mean(tem[,4] + tem[,5])
    dbslmmcode = paste(paste0('Rscript ', dbslmm_path, 'software/DBSLMM.R'),
                      paste0('--summary ', summf, chr, '.assoc.txt'), 
                      paste0('--outPath ', outPath),
                      paste0('--plink ', plink_path, 'plink'),
                      paste0('--dbslmm ', dbslmm_path, 'dbslmm'),
                      paste0('--ref ', PennPRS_path, 'LD/', race.dbslmm, '/1KGref_plinkfile/chr', chr),
                      paste0('--model DBSLMM'),
                      paste0('--n ', n),
                      paste0('--nsnp ', nrow(tem)),
                      paste0('--block ', dbslmm_path, 'block_data/', race.dbslmm, '/chr', chr, '.bed'),
                      paste0('--h2 ', ldsc_h2_est),
                      paste0('--thread ', NCORES))
    system(dbslmmcode)
    print(paste0('Complete training ', method, ' for CHR ', chr))
  }
  
  
  # --------------------- Step 1.3: Reformat the trained PRS weight file: ---------------------
  score = NULL
  for(chr in c(1:22)){
    temfile = paste0(prsdir, 'output/chr', chr, '.dbslmm.txt')
    if(file.exists(temfile)){
      scoretemp = bigreadr::fread2(temfile)[, c(1, 2, 4)] # BETA corresponds to A1
      colnames(scoretemp) = c('SNP', 'A1', 'BETA')
      scoretemp$CHR = chr
      ref = bigreadr::fread2(paste0(eval_ld_ref_path, 'chr', chr, '.bim'))[,c(2,5,6)]
      colnames(ref) = c('SNP', 'A1.ref', 'A2.ref')
      scoretemp = merge(scoretemp, ref, by = 'SNP')
      scoretemp$A2 = ifelse(scoretemp$A1 == scoretemp$A1.ref, scoretemp$A2.ref, scoretemp$A1.ref)
      scoretemp = scoretemp[,c(4,1,2,7,3)]
      colnames(scoretemp) = c('CHR','SNP','A1','A2','BETA')
      score = rbind(score, scoretemp)
      rm(scoretemp)
      # print(paste0('Chr ', chr,' Completed'))
    }
    if(!file.exists(temfile)) print(paste0('Chromosome ', chr, ' has zero SNP left in the model.'))
  }
  print(paste0('Combining fitted models across chromosomes: completed.'))
  score = score[,c('CHR','SNP','A1','A2','BETA')]
  
  # Match alleles with GWAS summary data:
  stateval = bigreadr::fread2(paste0(workdir, 'sumdata/',trait_name,'.txt'))
  stateval = stateval[, c('SNP','A1','A2')]
  colnames(stateval) = c('SNP', 'A1.ref','A2.ref')
  stateval = merge(stateval, score, by = 'SNP')
  flipped = which(stateval$A1.ref != stateval$A1)
  print(paste0(length(flipped), ' flipped SNPs.'))
  if (length(flipped) > 0){
    stateval[flipped,'A1'] = stateval[flipped,'A1.ref']
    stateval[flipped,'A2'] = stateval[flipped,'A2.ref']
    stateval[flipped,paste0('BETA')] = - stateval[flipped,paste0('BETA')]
  }
  scores = stateval[,c('CHR','SNP','A1','A2','BETA')] # other files: SNP	CHR	A1	BETA1	BETA2	A2
  write_delim(scores, file = paste0(workdir, trait_name,'.',method,'.PRS.txt'), delim='\t')
  # write.table(scores, paste0(input_path, trait_name,'.',method, '.ite',ite,'.txt'), row.names = F,col.names = T, quote = FALSE, sep = "\t" )
  if (opt$verbose >= 1 ){
    print(paste0('****************************************'))
    print(paste0('******* Complete training DBSLMM *******'))
    print(paste0('****************************************'))
  }
}



# ------------------------------------------------------------------------------------------------------
# ------------------------ Step 5: Write log files and delete intermediate files -----------------------
# ------------------------------------------------------------------------------------------------------
filen<-paste0(workdir, 'PRS_model_training_info_', trait_name, '.txt')
file.create(filen)
zz <- file(filen, "w")

print.title = paste0("Summary of PRS model Training on ",trait, " for ", race)
cat(paste0("\n",paste(rep('*', nchar(print.title)+10),collapse='')), file = zz)
cat(paste0("\n**** ", print.title, " ****"), file = zz)
cat(paste0("\n",paste(rep('*', nchar(print.title)+10),collapse=''),'\n'), file = zz)

if ('PRS-CS-auto' %in% methods){
  cat(paste0("\n******************************************************"), file = zz)
  cat(paste0("\n******************** PRS-CS-auto *********************"), file = zz)
  cat(paste0("\n******************************************************"), file = zz)
  method = 'PRS-CS-auto'
  prsdir = paste0(prsdir0, method,'/')
  tfile = paste0(workdir, trait_name,'.',method,'.PRS.txt')
  if (file.exists(tfile)){
    if (is.na(phi)){
      cat(paste0('\nSince phi is not specified, it was learnt from the data using a fully Bayesian approach. A PRS model was then trained by PRS-CS-auto with the estimated phi value.\n'), file = zz)
    }
    if (!is.na(phi)){
      cat(paste0('\nA PRS model was trained by PRS-CS-auto with a pre-specified phi value ', phi, '.\n'), file = zz)
    }
  }
  if (!file.exists(tfile)){
    cat(paste0('\nNo PRS model was generated. Please check log file for potential issues.\n'), file = zz)
  }
}

if ('LDpred2-auto' %in% methods){
  cat(paste0("\n*****************************************************"), file = zz)
  cat(paste0("\n******** LDpred2-auto (June 8, 2023 Version) ********"), file = zz)
  cat(paste0("\n*****************************************************"), file = zz)
  cat(paste0("\n* Please refer to https://privefl.github.io/bigsnpr/articles/LDpred2.html#ldpred2-auto-automatic-model for detailed implementation of LDpred2-auto."), file = zz)
  method = 'LDpred2-auto'
  prsdir = paste0(prsdir0, method,'/')
  tfile = paste0(workdir, trait_name,'.',method,'.PRS.txt')
  if (file.exists(tfile)){
    cat(paste0("\n* Parameter specifications:"), file = zz)
    cat(paste0('\n  - Shrinkage multiplicative coefficient applied to the off-diagonal elements of the correlation matrix: shrink_corr = ', coef_shrink), file = zz)
    cat(paste0('\n  - Allow for effects sizes to change sign in consecutive iterations? ', ifelse(allow_jump_sign, 'Yes.', 'No.')), file = zz)
    cat(paste0('\n  - ', ifelse(use_MLE, 'Used maximum likelihood estimation (MLE) to estimate alpha and the variance component (since version 1.11.4).', 
                              'Assume alpha = -1 and estimate the variance of (scaled) effects by h2/(m*p), as in earlier versions (e.g. v1.10.8).')), file = zz)
    if (length(keep) == 0){
      cat(paste0('\n##### [Note] All 30 chains in ', method, ' were deemed bad chains, and the resulting PRS model may not have sufficient power. \nPotential explanations:\n 1. The trait is not heritable.\n 2. The GWAS have insufficient power (e.g., due to low sample size) to develop a predictive PRS.\n 3. Issues with the input GWAS summary data (e.g., problematic BETA or SE).\n 4. ', method, ' is not powerful for developing PRS for the trait, in which case other methods can be considered.'), file = zz)
    }
  }
  if (!file.exists(tfile)){
    cat(paste0('\nNo PRS model was generated. Please check log file for potential issues.\n'), file = zz)
  }
}

if ('DBSLMM' %in% methods){
  cat(paste0("\n*****************************************************"), file = zz)
  cat(paste0("\n******************* DBSLMM (V0.3) *******************"), file = zz)
  cat(paste0("\n*****************************************************"), file = zz)
  cat(paste0("\n* Please refer to https://biostat0903.github.io/DBSLMM/Manual.html for detailed implementation of DBSLMM."), file = zz)
  method = 'DBSLMM'
  prsdir = paste0(prsdir0, method,'/')
  tfile = paste0(workdir, trait_name,'.',method,'.PRS.txt')
  if (file.exists(tfile)){
    cat(paste0("\n* DBSLMM default version was implemented with:"), file = zz)
    cat(paste0('\n  - p-value threshold: 1e-06'), file = zz)
    cat(paste0('\n  - LD threshold: 0.2'), file = zz)
  }
  if (!file.exists(tfile)){
    cat(paste0('\nNo PRS model was generated. Please check log file for potential issues.\n'), file = zz)
  }
}

close(zz)



# ---------------------------- README file:
filen<-paste0(workdir,'README.txt')
file.create(filen)
zz <- file(filen, "w")

print.title = paste0("List of Contents:") # ,trait, " for ", race
cat(paste0("\n",paste(rep('*', nchar(print.title)+10),collapse='')), file = zz)
cat(paste0("\n* ", print.title, " ****"), file = zz)
cat(paste0("\n",paste(rep('*', nchar(print.title)+10),collapse=''),'\n'), file = zz)
cat(paste0('\n* PRS models trained by tuning-parameter-free methods:'), file = zz)
for (method in methods){
  prsfile = paste0(workdir, trait_name,'.',method, '.PRS.txt')
  if (file.exists(prsfile)) cat(paste0('\n  ', trait_name,'.',method, '.PRS.txt'), file = zz)
}

cat(paste0('\n\n* Details of the PRS training:'), file = zz)
cat(paste0('\n  PRS_model_training_info_', trait_name, '.txt\n'), file = zz)
close(zz)




# Clean up intermediate files:
unlink(paste0(workdir, 'input_for_eval/'), recursive = TRUE, force = TRUE)
unlink(paste0(workdir, 'sumdata/'), recursive = TRUE, force = TRUE)
unlink(paste0(workdir, 'output/'), recursive = TRUE, force = TRUE)
unlink(paste0(workdir, 'output_for_eval/'), recursive = TRUE, force = TRUE)
unlink(paste0(workdir, 'PRS_model_training/'), recursive = TRUE, force = TRUE)
unlink(paste0(workdir, 'PennPRS_results/'), recursive = TRUE, force = TRUE)

