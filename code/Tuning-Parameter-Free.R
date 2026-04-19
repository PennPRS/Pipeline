library(optparse)
library(parallel)
library(readr)
library(bigreadr)
library(bigsnpr)
library(data.table)
library(dplyr)
library(scales)
library(stringr) # for str_split
# library(plyr)
library(rmarkdown)
library(readxl)

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
              help="Number of cores used for parallel computing of LDpred2-auto.
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


source(paste0(PUMAS_path, 'PennPRS_functions.R'))
source(paste0(PUMAS_path, 'gwas_qc_report_generator.R'))
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

# Read in LD reference file to add SNP position info for output PRS files:
ref.bim = bigreadr::fread2(paste0(eval_ld_ref_path, LDrefpanel,'_hm3_',race,'_ref.bim'))[, c(2,4)]
colnames(ref.bim) = c('SNP', 'chr_position')


######## QC for GWAS Summary Data:


# # copy the input GWAS summary data, {Ancestry}_{Trait}.txt, to the /sumdata/ folder
# system(paste0('cp -r ',input_GWAS_path, trait_name,'.txt ', workdir, 'sumdata/'))


cat(paste0("\n********************************************************"))
cat(paste0("\n********** Step 0: QC for the input GWAS data **********"))
cat(paste0("\n********************************************************\n"))

write_QC_report(race, trait, input_GWAS_path, workdir, PennPRS_path)



if (('LDpred2-auto' %in% methods) | ('DBSLMM' %in% methods)){
  method = 'LDpred2-auto'
  prsdir = paste0(prsdir0, method,'/')
  map_ldref <- readRDS(paste0(ld_path, '/map/map_',LDrefpanel,'_ldref.rds'))
  
  sumraw = bigreadr::fread2(paste0(workdir, 'sumdata/', trait_name, '.txt'))
  sumstats = sumraw[,c('CHR','SNP','A1','A2','BETA','SE','P','N', 'MAF')]; rm(sumraw)
  sumstats$P = as.numeric(sumstats$P)
  sumstats$BETA = as.numeric(sumstats$BETA)
  sumstats$SE = as.numeric(sumstats$SE)
  sumstats$N = as.numeric(sumstats$N)
  sumstats$MAF = as.numeric(sumstats$MAF)
  names(sumstats) <- c("chr", "rsid", "a0", "a1", "beta", "beta_se", "p", "n_eff", "a1_sumdata_af")
  
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
  # dbslmm_h2 <- list()
  
  ld = NULL
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
      if ((chr == 1) | (is.null(ld))) {
        ld.temp <- Matrix::colSums(corr0^2)
        ld <- ld.temp
        corr <- as_SFBM(corr0, tmp, compact = TRUE)
      } else {
        if (length(corr0) == 1) corr0 =  as(1, "sparseMatrix") # as(corr0, "sparseMatrix") # as.matrix(corr0, 1, 1)
        ld.temp <- Matrix::colSums(corr0^2)
        ld <- c(ld, ld.temp)
        corr$add_columns(corr0, nrow(corr))
      }
      # (ldsc.temp <- with(df_beta[df_beta$chr == chr, ], snp_ldsc(ld.temp, length(ld.temp), chi2 = (beta / beta_se)^2,
      #                                 sample_size = n_eff, blocks = NULL)))
      # dbslmm_h2[[chr]] <- abs(ldsc.temp[["h2"]])
      print(paste0('Complete calculating LD for CHR ', chr))
      rm(corr0)
    }
  }
  
  (ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2,
                                  sample_size = n_eff, blocks = NULL)))
  ldsc_h2_est <- abs(ldsc[["h2"]])
  cat(paste0('Heritability estimate based on LD score regression: ', signif(ldsc[["h2"]], 3), '\n'))
  if (ldsc[["h2"]] < 0) cat(paste0('Warning: negative hertability estimate based on LD score regression.\n'))
  if ('LDpred2-auto' %in% methods){
    if ( opt$verbose >= 1 ){
      print(paste0('************************************************************'))
      print(paste0('****** Start training PRS model based on LDpred2-auto ******'))
      print(paste0('************************************************************'))
    }
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
      beta_auto0 = data.frame(df_beta[,c('chr','rsid','a0','a1')], beta_auto)
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
      cat(paste0('[Warning] All 30 chains in ', method, ' were deemed bad chains, and the resulting PRS model may not have sufficient power. \nPotential explanations:\n 1. The trait is not heritable.\n 2. The GWAS have insufficient power (e.g., due to low sample size) to develop a predictive PRS.\n 3. Issues with the input GWAS summary data (e.g., problematic BETA or SE).\n 4. ', method, ' is not powerful for developing PRS for the trait, in which case other methods can be considered.'))
      beta_auto <- rowMeans(sapply(multi_auto, function(auto) auto$beta_est))
      beta_auto0 = data.frame(df_beta[,c('chr','rsid','a0','a1')], beta_auto)
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
  if (!file.exists(gwasinput)){
    print(paste0('A valid GWAS summary data file is missing.'))
    q()
  } 
  if (file.exists(gwasinput)){
    suppressWarnings(dir.create(paste0(prsdir, 'summary_gemma/')))
    sumraw0 = bigreadr::fread2(gwasinput)
    # impute position info:
    ref = bigreadr::fread2(paste0(dbslmm_path, 'LDref/', race.dbslmm, '/merge.bim'))[,c(2,4,5,6)]
    colnames(ref) = c('SNP', 'ps', 'ref', 'alt')
    ref = ref[ref$SNP %in% sumraw0$SNP, ]
    sumraw0 = merge(sumraw0, ref, by = 'SNP')
    sumraw0$n_mis = max(sumraw0$N) - sumraw0$N
    # # Match with reference data:
    # flipped = which(sumraw0$ref != sumraw0$A1)
    # # print(paste0(length(flipped), ' flipped SNPs.'))
    # if (length(flipped) > 0){
    #   sumraw0[flipped,'A1'] = sumraw0[flipped,'ref']
    #   sumraw0[flipped,'A2'] = sumraw0[flipped,'alt']
    #   sumraw0[flipped,paste0('BETA')] = - sumraw0[flipped,paste0('BETA')]
    # }
    sumraw0 = sumraw0[,c('CHR', 'SNP', 'ps', 'n_mis', 'N', 'A1', 'A2', 'MAF', 'BETA', 'SE', 'P')] # allele1 <-> A1: REF: effect allele
    colnames(sumraw0) = c('chr', 'rs',  'ps',  'n_mis',   'n_obs',   'allele1', 'allele0', 'af',  'beta', 'se', 'p_wald')
    # allele1: REF, effect allele
    
    # for (chr in 1:22){
    #   sumdat.file = paste0(prsdir, 'summary_gemma/chr', chr, '.assoc.txt')
    #   sumraw = sumraw0[sumraw0$chr == chr, ]
    #   write_delim(sumraw, sumdat.file, delim = '\t', col_names = F)
    #   rm(sumraw)
    # }
    # # ALL CHROMOSOMES:
    sumdat.file = paste0(prsdir, 'summary_gemma/all.assoc.txt')
    write_delim(sumraw0, sumdat.file, delim = '\t', col_names = T)
    rm(sumraw0)
    print(paste0('Generating input GWAS data for ', method, ': completed.'))
    
    # --------------------- Step 1.2: Run DBSLMM ---------------------
    system(paste0('chmod 777 ', dbslmm_path, 'dbslmm'))
    summf = paste0(prsdir, 'summary_gemma/')
    outPath = paste0(prsdir, 'output/')
    tem = bigreadr::fread2(paste0(summf, 'all.assoc.txt'))
    n = mean(tem[,4] + tem[,5])
    dbslmmcode = paste(paste0('Rscript ', dbslmm_path, 'software/DBSLMM.R'),
                       paste0('--summary ', summf, 'all.assoc.txt'), 
                       paste0('--outPath ', outPath),
                       paste0('--type auto'),
                       paste0('--N ', n),
                       paste0('--dbslmm ', dbslmm_path, 'dbslmm'),
                       # paste0('--ref ', PennPRS_path, 'LD/', race.dbslmm, '/1KGref_plinkfile/1kg_hm3_EUR_ref'),
                       paste0('--reference ', dbslmm_path, 'LDref/', race.dbslmm),
                       paste0('--model DBSLMM'),
                       paste0('--block ', dbslmm_path, 'block_data/', race.dbslmm, '/'))
    system(dbslmmcode)
    print(paste0('Complete training ', method, ' PRS across all chromosomes.'))
    
    
    # --------------------- Step 1.3: Reformat the trained PRS weight file: ---------------------
    ref = ref[, c(1,3,4)]
    colnames(ref) = c('SNP', 'A1.ref', 'A2.ref')
    
    score = NULL
    for(chr in c(1:22)){
      temfile = paste0(prsdir, 'output/all_chr', chr, '.dbslmm.txt')
      if(file.exists(temfile)){
        if (file.info(temfile)$size > 0){
        scoretemp = bigreadr::fread2(temfile)[, c(1, 2, 4)] # BETA corresponds to A1
        colnames(scoretemp) = c('SNP', 'A1', 'BETA')
        scoretemp$CHR = chr
        scoretemp = merge(scoretemp, ref, by = 'SNP')
        scoretemp$A2 = ifelse(scoretemp$A1 == scoretemp$A1.ref, scoretemp$A2.ref, scoretemp$A1.ref)
        scoretemp = scoretemp[,c(4,1,2,7,3)]
        colnames(scoretemp) = c('CHR','SNP','A1','A2','BETA')
        score = rbind(score, scoretemp)
        rm(scoretemp)
        # print(paste0('Chr ', chr,' Completed'))
        }
      }
      if(!file.exists(temfile)) print(paste0('No SNP in Chromosome ', chr, ' was included in the PRS model.'))
    }
    print(paste0('Combining trained PRS models across chromosomes: completed.'))
    
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
}


# Save the final PRS models generated by each single methods to the working directory
for (method in methods){
  tfile = paste0(workdir, trait_name,'.',method,'.PRS.txt')
  if (file.exists(tfile)) {
    SCORE = bigreadr::fread2(tfile)
    # Update output file format according to the pgsc_calc pipeline from the PGS Catalog
    SCORE = merge(SCORE, ref.bim, by = 'SNP')
    SCORE = SCORE[, c('CHR','SNP','chr_position','A1','A2', 'BETA')]
    colnames(SCORE) = c('chr_name','rsid','chr_position','effect_allele','other_allele', 'effect_weight')
    write_delim(SCORE, paste0(workdir, trait_name,'.',method,'.PRS.txt'))
  }
}


# ------------------------------------------------------------------------------------------------------
# ------------------------ Step 5: Write log files and delete intermediate files -----------------------
# ------------------------------------------------------------------------------------------------------
filen<-paste0(workdir, 'PRS_INFO.txt')
file.create(filen)
zz <- file(filen, "wt")

on.exit({
  try(close(zz), silent = TRUE)
}, add = TRUE)

# print.title = paste0("Summary of PRS model Training on ",trait, " for ", race)
# cat(paste0("\n",paste(rep('*', nchar(print.title)+10),collapse='')), file = zz)
# cat(paste0("\n**** ", print.title, " ****"), file = zz)
# cat(paste0("\n",paste(rep('*', nchar(print.title)+10),collapse=''),'\n'), file = zz)

write_line <- function(...) {
  cat(paste0(...), "\n", file = zz, append = TRUE, sep = "")
}

write_section <- function(title) {
  write_line("")
  write_line(paste(rep("=", 70), collapse = ""))
  write_line(title)
  write_line(paste(rep("=", 70), collapse = ""))
}

write_subsection <- function(title) {
  write_line("")
  write_line(title)
  write_line(paste(rep("-", nchar(title)), collapse = ""))
}

# --------------------------------------------------
# Header
# --------------------------------------------------
write_line("PRS INFORMATION REPORT")
write_line("")
write_line(paste0("Trait: ", trait))
write_line(paste0("Ancestry: ", race))
write_line(paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
write_line(paste0("Methods requested: ", paste(methods, collapse = ", ")))

# --------------------------------------------------
# Overview
# --------------------------------------------------
write_section("1. OVERVIEW")
write_line("This file summarizes details of the PRS model training process and provides example")
write_line("commands for calculating PRS using the trained model.")

# --------------------------------------------------
# Training summary
# --------------------------------------------------
write_section("2. SUMMARY OF PRS MODEL TRAINING")
methods.completed = NULL
if ("LDpred2-auto" %in% methods) {
  method <- "LDpred2-auto"
  prsdir <- paste0(prsdir0, method, "/")
  tfile <- paste0(workdir, trait_name, ".", method, ".PRS.txt")
  
  write_subsection("LDpred2-auto (June 8, 2023 Version)")
  write_line("* Reference documentation:")
  write_line("  https://privefl.github.io/bigsnpr/articles/LDpred2.html#ldpred2-auto-automatic-model")
  
  if (file.exists(tfile)) {
    methods.completed = c(methods.completed, method)
    write_line("* Model training status: completed.")
    write_line(paste0("* Output score file: ", trait_name, ".", method, ".PRS.txt"))
    write_line("* Parameter specifications:")
    write_line('  - Shrinkage multiplicative coefficient applied to the off-diagonal elements of the correlation matrix: shrink_corr = ', coef_shrink)
    write_line('  - Allow for effects sizes to change sign in consecutive iterations? ', ifelse(allow_jump_sign, 'Yes.', 'No.'))
    write_line('  - ', ifelse(use_MLE, 'Used maximum likelihood estimation (MLE) to estimate alpha and the variance component (since version 1.11.4).', 
                                'Assume alpha = -1 and estimate the variance of (scaled) effects by h2/(m*p), as in earlier versions (e.g. v1.10.8).'))
    
    # write_line("  - shrink_corr = ", coef_shrink)
    # write_line("  - allow_jump_sign = ", ifelse(isTRUE(allow_jump_sign), "Yes", "No"))
    # write_line(
    #   "  - alpha / variance estimation: ",
    #   ifelse(isTRUE(use_MLE),
    #          "maximum likelihood estimation (MLE) was used",
    #          "alpha was fixed at -1 and variance was estimated as in earlier LDpred2-auto versions")
    # )
    
    if (length(keep) == 0) {
      write_line("")
      write_line("WARNING:")
      write_line("  All 30 LDpred2-auto chains were deemed bad chains.")
      write_line("  The resulting PRS model may have limited predictive power.")
      write_line("  Possible explanations include:")
      write_line("    1. The trait has low heritability.")
      write_line("    2. The GWAS summary statistics have insufficient power (e.g., due to small sample size).")
      write_line("    3. The input GWAS summary data contains problematic values (e.g., BETA or SE).")
      write_line("    4. LDpred2-auto may not be well suited for this trait, and other methods can be considered instead.")
    }
  } else {
    write_line("* Model training status: no PRS model was generated.")
    write_line("* Please check log file for details.")
  }
}

if ("DBSLMM" %in% methods) {
  method <- "DBSLMM"
  prsdir <- paste0(prsdir0, method, "/")
  tfile <- paste0(workdir, trait_name, ".", method, ".PRS.txt")
  
  write_subsection("DBSLMM (V1.0 User Friendly Version)")
  write_line("* Reference documentation:")
  write_line("  https://github.com/biostat0903/DBSLMM")
  
  if (file.exists(tfile)) {
    methods.completed = c(methods.completed, method)
    write_line("* Model training status: completed")
    write_line(paste0("* Output score file: ", trait_name, ".", method, ".PRS.txt"))
    write_line("* DBSLMM default (automatic) version used:")
    write_line("  - p-value threshold = 1e-06")
    write_line("  - LD threshold = 0.2")
  } else {
    write_line("* Model training status: no PRS model was generated.")
    write_line("* Please check log file for details.")
  }
}


# --------------------------------------------------
# Example: compute PRS with PLINK2
# --------------------------------------------------
if (length(methods.completed) > 0){
  write_section("3. EXAMPLE CODE FOR COMPUTING PRS BASED ON THE GENERATED SCORE FILES")
  write_line(paste0("   Example genotype data (prefix): PennPRS/test/evaldir/eval.{bim,bed,fam}"))
  example_score_file = paste0(trait_name, ".", methods.completed[1], ".PRS.txt")
  write_line(paste0("   Example score file: ", example_score_file))
  
  write_section("3.1. Example Command for Comuting PRS using PLINK2")
  # write_line("")
  # write_line("Expected score file columns:")
  # write_line("  1. Variant ID")
  # write_line("  2. Effect allele")
  # write_line("  3. SNP weight")
  write_line("")
  
  # write_line("Example score file:")
  # write_line(example_score_file)
  
  # write_line("")
  write_line("Example PLINK2 command:")
  write_line("")
  write_line("PennPRS/software/plink2 \\")
  write_line("    --bfile PennPRS/test/evaldir/eval \\")
  write_line("    --score ", example_score_file, " 2 3 5 cols=+scoresums,-scoreavgs \\")
    write_line("    --out PRS_", trait_name, ".", methods.completed[1])
    write_line("    --threads 1")
    write_line("")
    write_line("Notes:")
    write_line("  - Replace the score file path and genotype data path with your actual file paths.")
    # write_line("  - Try the example provided in {PennPRS/test/Testing PennPRS with an Example.md}.")
    write_line("  - See https://www.cog-genomics.org/plink/2.0/score for further information.")
    
    # --------------------------------------------------
    # Example: compute PRS with pscs_calc
    # --------------------------------------------------
    write_section("3.2. COMPUTING PRS USING pscs_calc")
    write_line("The pgsc_calc is a tool for calculating PRS using score files published in the PGS Catalog or custom scoring files.")
    write_line("To calculate PRS using pgsc_calc, please install pgsc_calc following the instructions at https://pgsc-calc.readthedocs.io/en/latest/ first.")
    write_line("Prepare required files:")
    write_line("   1. Genotype data (pfile, bfile, or vcf) with prefix: temppath/genotype_data (e.g., PennPRS/test/evaldir/eval.{bim,bed,fam}).")
    write_line("   2. Samplesheet.csv (store genotype data information, save under the same folder: temppath).")
    write_line("   3. Update score file with format required by pgsc_calc: temppath/", trait_name, ".", methods.completed[1], "_for_pgsc_calc.txt (save under the same folder: temppath)")
    
    write_line("")
    write_line("Example pscs_calc command:")
    write_line("")
    write_line("path/to/nextflow run path/to/pgsc_calc \\")
    write_line("    -profile test,docker \\")
    write_line("    --input temppath/samplesheet.csv \\")
    write_line("    --scorefile temppath/", trait_name, ".", methods.completed[1], "_for_pgsc_calc.txt")
    write_line("")
    write_line("Notes:")
    write_line("  - Replace path/to/nextflow, path/to/pgsc_calc, and temppath/ with your actual paths.")
    write_line("  - Generated PRS file can be found in pgsc_calc/results/.")
    # write_line("  - Try the example provided in {PennPRS/test/Testing PennPRS with an Example.md}.")
    write_line("  - See https://pgsc-calc.readthedocs.io/en/latest/ for further information.")
}

write_line("")
# write_line("PRS information report completed.")

invisible(filen)





# ---------------------------- README file:
filen<-paste0(workdir,'README.txt')
file.create(filen)
zz <- file(filen, "w")

print.title = paste0("List of Contents:") # ,trait, " for ", race
cat(paste0("\n",paste(rep('*', nchar(print.title)+10),collapse='')), file = zz)
cat(paste0("\n**** ", print.title, " ****"), file = zz)
cat(paste0("\n",paste(rep('*', nchar(print.title)+10),collapse=''),'\n'), file = zz)

cat(paste0('\n* Report on GWAS QC procedure:'), file = zz)
cat(paste0('\n  QC_report.html\n'), file = zz)

cat(paste0('\n* PRS models trained by tuning-parameter-free methods:'), file = zz)
for (method in methods){
  prsfile = paste0(workdir, trait_name,'.',method, '.PRS.txt')
  if (file.exists(prsfile)) cat(paste0('\n  ', trait_name,'.',method, '.PRS.txt'), file = zz)
}

cat(paste0('\n\n* Details of PRS model training and instructions on PRS computing:'), file = zz)
cat(paste0('\n  PRS_INFO.txt\n'), file = zz)
close(zz)


cat(paste0('Job completed. Results are saved in: ', workdir))


# Clean up intermediate files:
unlink(paste0(workdir, 'input_for_eval/'), recursive = TRUE, force = TRUE)
unlink(paste0(workdir, 'sumdata/'), recursive = TRUE, force = TRUE)
unlink(paste0(workdir, 'output/'), recursive = TRUE, force = TRUE)
unlink(paste0(workdir, 'output_for_eval/'), recursive = TRUE, force = TRUE)
unlink(paste0(workdir, 'PRS_model_training/'), recursive = TRUE, force = TRUE)
unlink(paste0(workdir, 'PennPRS_results/'), recursive = TRUE, force = TRUE)

