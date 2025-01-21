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
  make_option("--methods", action = "store", default = 'C+T,lassosum2,LDpred2', type = "character",
              help="Options: a subset of methods from C+T, lassosum2, and LDpred2, divided by comma"),
  make_option("--trait", action = "store", default = NA, type = "character",
              help="trait name [Optional]"),
  make_option("--race", action = "store", default = NA, type = "character",
              help="Race of the training GWAS individuals. Options: EUR (European), AFR (African), 
              AMR (Mixed American, Hispanic/Latio), EAS (East Asian), or SAS (South Asian) [Required]"),
  make_option("--LDrefpanel", action = "store", default = '1kg', type = "character",
              help="LD reference panel. Options: '1kg' (1000 Genomes Project Phase 3) or 'ukbb' (UK Biobank, will be available soon) [Optional]"),
  make_option("--k", action = "store", default = 2, type = "numeric",
              help = "k-fold Monte Carlo Cross Validation (MCCV) for PUMAS. Options: any integer greater than or equal to 2 [Optional]"),
  
  make_option("--partitions", action = "store", default = '0.8,0.2', type = "character",
              help="Partitions for PUMAS subsampling. 
              Format: '% training, % testing' (equal to 1 - % training), divided by comma [Optional]"),
  
  make_option("--delta", action = "store", default = '0.001,0.01,0.1,1', type = "character",
              help="Candidate values of the shrinkage parameter in L2 regularization. Options: 
              candidate values in (0,Infinity), divided by comma [Optional]"),
  make_option("--nlambda", action = "store", default = 30, type = "numeric",
              help="Number of different candidate values for lambda (shrinkage parameter in the L1 regularization.
              Options: any positive integer [Optional]"),
  make_option("--lambda.min.ratio", action = "store", default = 0.01, type = "numeric",
              help="Ratio between the lowest and highest candidate values of lambda.
              Options: any value in (0,1) [Optional]"),
  
  make_option("--alpha", action = "store", default = '0.7,1.0,1.4', type = "character",
              help="H_2 = alpha * H_20, where H_20 is the heritability estimated by LD score regression.
              Recommended alternative: default values in the LDpred2 algorithm (June 8, 2023 version): 0.3,0.7,1.0,1.4.
              Options: candidate values in (0,Infinity), divided by comma [Optional]"),
  make_option("--p_seq", action = "store", default = '1.0e-05,3.2e-05,1.0e-04,3.2e-04,1.0e-03,3.2e-03,1.0e-02,3.2e-02,1.0e-01,3.2e-01,1.0e+00', 
              type = "character",
              help="A sequence of candidate values for the proportion of causal SNPs.
              Options: candidate values in(0,1], divided by comma [Optional]"),
  make_option("--sparse", action = "store", default = FALSE, type = "character",
              help="Whether we consider a sparse (i.e., the majority of SNPs have effects shrunk to zero) effect size distribution.
              Options: candidate values in {FALSE, TRUE}, divided by comma [Optional]"),
  
  make_option("--NCORES", action = "store", default = 5, type = "numeric",
              help="Number of cores (default: 5) used for parallel computing of LDpred2 and lassosum2.
              Options: positive integer [Optional]"),
  
  make_option("--kb", action = "store", default = 500, type = "numeric",
              help="SNPs within this range of the index SNP are considered for the C (LD Clumping) step.
              Options: any positive integer [Optional]"),
  make_option("--Pvalthr", action = "store", default = '5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01', type = "character",
              help="p-value threshold for the T (p-value Thresholding) step.
              Options: candidate values in(0,1], divided by comma [Optional]"),
  make_option("--R2", action = "store", default = '0.1', type = "character",
              help="SNPs having squared correlation higher than r2 with the index SNPs will be removed. 
              Options: candidate values in (0,1], divided by comma [Optional]"),
  
  make_option("--ensemble", action = "store", default = F, type = "logical",
              help="Whether to train a weighted combination of the single PRS models.
              Options: T, F [default: %default]"),
  make_option("--verbose", action="store", default=1, type="integer",
              help="Print logfile? 0 = no; 1 = yes [default: %default]")
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
# Parameters for subsampling
k = opt$k
NCORES = opt$NCORES
# if multiple methods are selected, also input ensemble into opt
ensemble = FALSE
if (length(methods) > 1) ensemble = opt$ensemble ######## Need Bingxuan's input

# Optional input parameters:
partitions <- opt$partitions

ld_path <- paste0(PennPRS_path, '/LD/', race, '/')
PUMAS_path = paste0(PennPRS_path,'/code/')
plink_path = paste0(PennPRS_path, '/software/')
threads = 1


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


# Parameters
if ('C+T' %in% methods){
  # Parameters for the Clumping step
  kb = opt$kb # SNPs within 500kb of the index SNP are considered for clumping
  p.ldclump = 0.5 # P-value threshold for a SNP to be included as an index SNP
  Pvalthr = as.numeric(str_split(opt$Pvalthr,",")[[1]]) # default: c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01)
  R2 = as.numeric(str_split(opt$R2,",")[[1]])
  params.ct = expand.grid(r2 = R2, pvalthr = Pvalthr)
}
if ('lassosum2' %in% methods){
  delta = as.numeric(str_split(opt$delta,",")[[1]]) # candidate values of the shrinkage parameter in L2 regularization
  nlambda = opt$nlambda # number of different candidate values for lambda (shrinkage parameter in the L1 regularization). Default in lassosum2 pipeline: 30, which may lead to issues when using PUMAS subsampling to tune parameters
  lambda.min.ratio = opt$lambda.min.ratio # Ratio between the lowest and highest candidate values of lambda. Candidate values in (0,Inf), divided by comma
}
if ('LDpred2' %in% methods){
  h2.ratio = as.numeric(str_split(opt$alpha,",")[[1]])
  p_seq <- as.numeric(str_split(opt$p_seq,",")[[1]])  # Default
  sp.temp = toupper(str_split(opt$sparse,",")[[1]])
  sparse.option = ifelse(sp.temp == 'TRUE', TRUE, FALSE) # opt$sparse # If TRUE: generate sparse effect size estimates. Options: subset of {FALSE, TRUE}
}


#source(paste0(PUMAS_path, 'PennPRS_functions.R')) # please save the PennPRS_functions.R file to the /PUMAS/code/ directory
gwas_path <- paste0(workdir, 'sumdata/')
output_path <- paste0(workdir, 'output/')
input_path <- paste0(workdir, 'input_for_eval/')
PennPRS_finalresults_path <- paste0(workdir, 'PennPRS_results/')
eval_ld_ref = paste0(eval_ld_ref_path,LDrefpanel,'_hm3_',race,'_ref') # or hm3+mega
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
if (ensemble) ensemble.methods = methods

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
# N.90percentile = quantile(sumraw$N, 0.1)
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
    write_delim(sumraw, paste0(workdir, 'sumdata/', trait_name, '.txt'), delim = '\t')
    if (length(rm.indx) == 1) print(paste0('* 1 problematic SNP removed. QC step completed.'))
    if (length(rm.indx) > 1) print(paste0('* QC step completed. ', nrow(sumraw), ' SNPs remaining. ', length(rm.indx), ' problematic SNPs removed.'))
  }
}
if (length(rm.indx) == 0) print(paste0('* QC step completed. ', nrow(sumraw), ' SNPs remaining. No SNP was removed.'))


# --------------------------------------------------------------------
# --------------------- Step 1: PUMAS Subsampling ---------------------
# --------------------------------------------------------------------
# Different from PRS-CS (auto) which doesn't have tuning parameters, C+T, lassosum2, and LDpred2 have tuning parameters, 
# and thus we need to conduct k=4 fold MCCV.
# Note!!! For a relatively large tuning sample size (here in our case ~300K), then there is no need to do MCCV and k=2 is sufficient.
# This is usually the case for EUR
# But for non-EUR races, it's very common that the tuning sample size is only several thousands
# If tuning sample size is below 2000 we use k=4 MCCV?, o.w. we just use k=2
if ( opt$verbose >= 1 ) {
  cat(paste0("\n***********************************"))
  cat(paste0("\n**** Step 1: PUMAS Subsampling ****"))
  cat(paste0("\n***********************************\n"))
}

if (T == T){
ld_file <- paste0(ld_path0,race,'_LD_hm3.RData')
rs_file <- paste0(ld_path0,race,'_rs_hm3.RData')
pumascode = paste(paste0('Rscript ', PUMAS_path, 'PUMAS.subsampling.customized.R '),
                  paste0('--k ',k), 
                  paste0('--partitions ',partitions),
                  paste0('--trait_name ',trait_name),
                  paste0('--gwas_path ', gwas_path),
                  paste0('--ld_file ', ld_file),
                  paste0('--rs_file ', rs_file),
                  paste0('--output_path ', output_path),
                  paste0('--threads ', threads))
system(pumascode)
}


# --------------------------------------------------------------------------------------------
# --------------------- Step 2: Train PRS models using different methods ---------------------
# --------------------------------------------------------------------------------------------
if ( opt$verbose >= 1 ) {
  cat(paste0("\n**************************************************************************************"))
  cat(paste0("\n**** Step 2: Train PRS models on pseudo training datasets using different methods ****"))
  cat(paste0("\n**************************************************************************************\n"))
}
if ('C+T' %in% methods){
  cat(paste0("\n*******************************************************"))
  cat(paste0("\n**** Start running C+T on pseudo training datasets ****"))
  cat(paste0("\n*******************************************************\n"))
  
  method = 'C+T'
  prsdir = paste0(prsdir0, method,'/')
  
  # Submit k separate jobs (for ite in 1:k) to the server and run them in parallel. 
  for (ite in 1:k){
    # --------------------- Step 2.1: Preparation for input files for C+T ---------------------
    pumasout = paste0(output_path, trait_name, '.gwas.ite', ite, '.txt')
    if (!file.exists(pumasout)) print(paste0('Subsampling failed to generate ', ite,'-th fold of the MCCV summary data. Rerun pumas.subsampling.R.'))
    if (file.exists(pumasout)){
      sumraw = bigreadr::fread2(pumasout)
      sumstats = sumraw[,c('CHR','SNP','A1','A2','BETA','SE','P','N', 'MAF')]
      colnames(sumstats) <- c("CHR", "SNP", "REF", "ALT", "BETA", "SE", "P", "N", "FRQ") # MAF or REF_FRQ
      rownames(sumstats) = sumstats$SNP
      # REF: effect allele
      print(paste0('Complete Loading GWAS summary data for iteration ', ite))
    }
    
    temdir = paste0(prsdir,'snplist/')
    if (!dir.exists(temdir)){dir.create(temdir)}
    temdir = paste0(prsdir,'clumped/')
    if (!dir.exists(temdir)){dir.create(temdir)}
    temdir = paste0(prsdir,'filtered/')
    if (!dir.exists(temdir)){dir.create(temdir)}
    
    # Create base data (summary statistic) file containing the P-value information:
    sumstats_input <- sumstats[,c('SNP','P')] #,'REF_FRQ'
    fwrite(sumstats_input, file=paste0(prsdir,'snplist/',trait_name,'_ite',ite,'.txt'),row.names = F, quote = F, sep=' ')
    
    # --------------------- Step 2.2: Run C+T ----------------------------
    set.seed(2023)
    for (r2 in R2){
      # --------------------- LD Clumping (C) ---------------------
      ldclumpcode <- paste0(plink_path, 'plink --bfile ', eval_ld_ref,
                            ' --clump ',prsdir,'snplist/',trait_name,'_ite',ite,'.txt',
                            ' --clump-p1 ',p.ldclump,
                            # ' --clump-p2 ',pc,
                            ' --clump-r2 ',r2,
                            ' --clump-kb ',kb,
                            ' --threads ', threads,
                            ' --silent',
                            ' --out ',prsdir,'clumped/',trait_name,'_ite',ite,'_r2=',r2)
      system(ldclumpcode)
      # --------------------- Thresholding (T) ---------------------
      LD <- bigreadr::fread2(paste0(prsdir,'clumped/',trait_name,'_ite',ite,'_r2=',r2,'.clumped'))
      clumped.snp <- LD[,3,drop=F][,1]
      sumstats.clumped <- sumstats[clumped.snp,]
      
      for (pvalthr in Pvalthr){
        keep.SNP = sumstats.clumped[sumstats.clumped$P <= pvalthr,c('SNP')]
        sumstats[,paste0('r2_', r2,'_p_',pvalthr)] = sumstats$BETA
        dump.SNP = which(!sumstats$SNP %in% keep.SNP)
        sumstats[dump.SNP, paste0('r2_', r2,'_p_',pvalthr)] = 0
      }
    }
    SCORE = sumstats[,c('CHR','SNP','REF','ALT',  sapply(1:nrow(params.ct), function(x){paste0('r2_', params.ct[x,'r2'],'_p_',params.ct[x,'pvalthr'])}))]
    colnames(SCORE)[3:4] = c('A1','A2')
    write_delim(SCORE, paste0(prsdir, trait_name,'.',method,'.ite',ite,".txt"), delim = '\t')
  }
  rm(SCORE)
}





if (('lassosum2' %in% methods) | ('LDpred2' %in% methods)){
  cat(paste0("\n****************************************************************************"))
  cat(paste0("\n**** Start running lassosum2 and/or LDpred2 on pseudo training datasets ****"))
  cat(paste0("\n****************************************************************************\n"))
  
  cat(paste0("\n*****************************************************************"))
  cat(paste0("\n**** Preparing LD info for training lassosum2 and/or LDpred2 ****"))
  cat(paste0("\n*****************************************************************\n"))
  map_ldref <- readRDS(paste0(ld_path, '/map/map_',LDrefpanel,'_ldref.rds'))
  
  # Submit k separate jobs (for ite in 1:k) to the server and run them in parallel. 
  # Each job will require < 20G, the memory depends on NCORES (if NCORES = 3 then perhaps 10G is enough).
  for (ite in 1:k){
    # --------------------- Step 2.1: Preparation for input files for each PRS method (here: lassosum2) ---------------------
    # Reformat summary data to use as the input data for lassosum2:
    pumasout = paste0(output_path, trait_name, '.gwas.ite', ite, '.txt')
    if (!file.exists(pumasout)) print(paste0('Subsampling failed to generate ', ite,'-th fold of the MCCV summary data. Rerun pumas.subsampling.R.'))
    if (file.exists(pumasout)){
      sumraw = bigreadr::fread2(pumasout)
      sumstats = sumraw[,c('CHR','SNP','A1','A2','BETA','SE','P','N', 'MAF')]
      names(sumstats) <- c("chr", "rsid", "a1", "a0", "beta", "beta_se", "p", "n_eff", "a1_sumdata_af")
      # a0: effect allele
      
      info_snp <- snp_match(sumstats, map_ldref, strand_flip = T, join_by_pos = F) # important: for real data, strand_flip = T
      info_snp <- tidyr:: drop_na(tibble::as_tibble(info_snp))
      sd_ldref <- with(info_snp, sqrt(2 * a1_af * (1 - a1_af)))
      sd_ss <- with(info_snp, sqrt(2 * a1_sumdata_af * (1 - a1_sumdata_af)))
      is_bad <- sd_ss < (0.5 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0.1 | sd_ldref < 0.05
      df_beta <- info_snp[!is_bad, ]
      print(paste0('Complete pre-processing GWAS summary data for iteration ', ite))
    }
    
    if (ite == 1){
      td = paste0(prsdir0, 'temporary_LDpred2_lassosum2_ite',ite)
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
    }
    
    if ('lassosum2' %in% methods){
      method = 'lassosum2'
      prsdir = paste0(prsdir0, method,'/')
      # Parameters
      delta = as.numeric(str_split(opt$delta,",")[[1]]) # candidate values of the shrinkage parameter in L2 regularization
      
      # --------------------- Step 2.2: Run lassosum2 ----------------------------
      set.seed(2023)
      beta_lassosum2 <- snp_lassosum2(corr, df_beta, # ncores = NCORES, 
                                      delta = delta, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio)
      if (ite == 1){
        params.lassosum2 <- attr(beta_lassosum2, "grid_param")
        write.table(params.lassosum2, file = paste0(prsdir, 'params.lassosum2.txt'), row.names = F)
      }
      
      beta_lassosum2[is.na(beta_lassosum2)] = 0
      # Further fix potential non-convergent issues:
      n.nonconvergent = sapply(1:nrow(params.lassosum2), function(x){sum(abs(beta_lassosum2[,x])>1)})
      indx.nonconvergent = which(n.nonconvergent > 0)
      if (length(indx.nonconvergent)>0) beta_lassosum2[,indx.nonconvergent] = 0
      
      beta_lassosum2 = data.frame(df_beta[,c('chr','rsid','a1','a0')], beta_lassosum2)
      colnames(beta_lassosum2) = c('chr','rsid','a1','a0', paste0('lassosum2_',1:nrow(params.lassosum2)))
      output_lassosum2 = paste0(prsdir, trait_name,'.',method,'.ite',ite,'.txt')
      # write_delim(beta_lassosum2,file = output_lassosum2, delim='\t')
      write_delim(beta_lassosum2, output_lassosum2, delim = '\t')
      
      print(paste0('** Complete training ', method, ' for MCCV ite ', ite, ' **'))
    }
    
    
    if ('LDpred2' %in% methods){
      method = 'LDpred2'
      prsdir = paste0(prsdir0, method,'/')
      
      (ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2,
                                      sample_size = n_eff, blocks = NULL)))
      ldsc_h2_est <- abs(ldsc[["h2"]])
      h2_seq <- round(ldsc_h2_est * h2.ratio, 5); 
      h2_seq[h2_seq == 0] = 1e-5
      h2_seq[duplicated(h2_seq)] = h2_seq[duplicated(h2_seq)] * 1.01
      n.inflated = sum(h2_seq>1)
      if (n.inflated > 0) h2_seq[h2_seq>1] = 0.95 + seq(0,0.01*(n.inflated-1), by = 0.01)
      params.ldpred2 <- expand.grid(p = p_seq, h2 = h2_seq, sparse = sparse.option)
      
      set.seed(2023)
      beta_ldpred2 <- snp_ldpred2_grid(corr, df_beta, params.ldpred2, ncores = NCORES)
      beta_ldpred2[is.na(beta_ldpred2)] = 0
      # Further fix potential non-convergent issues:
      # n.nonconvergent = sapply(1:nrow(params.ldpred2), function(x){sum(abs(beta_ldpred2[,x])>1)})
      # indx.nonconvergent = which(n.nonconvergent > 0)
      # if (length(indx.nonconvergent)>0) beta_ldpred2[,indx.nonconvergent] = 0
      
      beta_ldpred2 = data.frame(df_beta[,c('chr','rsid','a1','a0')], beta_ldpred2)
      colnames(beta_ldpred2) = c('chr','rsid','a1','a0', paste0('LDpred2_',1:nrow(params.ldpred2)))
      output_LDpred2 = paste0(prsdir, trait_name,'.',method,'.ite',ite,'.txt')
      write_delim(beta_ldpred2, output_LDpred2, delim = '\t')
      print(paste0('** Complete training ', method, ' for MCCV ite ', ite, ' **'))
    }
  }
}



# ---------------------------------------------------------------------------------------
# ----------------------------- Step 3: Single Model Tuning -----------------------------
# ---------------------------------------------------------------------------------------
# single_prs() in PUMA-CUBS.evaluation.R
if ( opt$verbose >= 1 ) {
  cat(paste0("\n********************************************************"))
  cat(paste0("\n******* Step 3: Parameter Tuning for Each Method *******"))
  cat(paste0("\n********************************************************\n"))
}


if ('C+T' %in% methods){
  method = 'C+T'
  prsdir = paste0(prsdir0, method,'/')
  
  # This step is also needed for C+T, these input files for evaluation will be used for calculating R2 on testing data and for training ensemble PRS
  for (ite in 1:k){
    output_ct = paste0(prsdir, trait_name,'.',method,'.ite',ite,'.txt')
    if(file.exists(output_ct)){
      score = bigreadr::fread2(output_ct) 
      n.tuning = ncol(score) - 4
      colnames(score) = c('CHR', 'SNP', 'A1', 'A2', paste0('BETA',1:n.tuning))
    }
    
    # Match alleles with GWAS summary data:
    sumstats = bigreadr::fread2(paste0(output_path, trait_name, '.gwas.ite', ite, '.txt'))
    stateval = sumstats[, c('SNP','A1','A2')]
    colnames(stateval) = c('SNP', 'A1.ref','A2.ref')
    stateval = left_join(stateval,score,by="SNP")
    
    na.ind = which(is.na(stateval$A1))
    if (length(na.ind) > 0){
      stateval[na.ind, paste0('BETA',1:n.tuning)] = 0
      stateval[na.ind,'A1'] = stateval[na.ind,'A1.ref']
      stateval[na.ind,'A2'] = stateval[na.ind,'A2.ref']
    }
    flipped = which(stateval$A1.ref != stateval$A1)
    print(paste0(length(flipped), ' flipped SNPs.'))
    if (length(flipped) > 0){
      stateval[flipped,'A1'] = stateval[flipped,'A1.ref']
      stateval[flipped,'A2'] = stateval[flipped,'A2.ref']
      for (t in 1:n.tuning){
        stateval[flipped,paste0('BETA',t)] = - stateval[flipped,paste0('BETA',t)]
      }
    }
    scores = stateval[,c('CHR','SNP','A1','A2',paste0('BETA',1:n.tuning))] # other files: SNP	CHR	A1	BETA1	BETA2	A2
    write_delim(scores, paste0(input_path, trait_name,'.',method, '.ite',ite,'.txt'), delim = '\t')
    # write.table(scores, paste0(input_path, trait_name,'.',method, '.ite',ite,'.txt'), row.names = F,col.names = T, quote = FALSE, sep = "\t" )
    rm(stateval)
  }
  
  xty_path = stats_path = output_path # the "output_path" used for storing output from pumas.subsampling.R
  R2.tuned = cbind(params.ct, matrix(0,length(R2)*length(Pvalthr),k))
  colnames(R2.tuned) = c('r2','pvalthr',paste0('ite',1:k))
  
  tunecode = paste(paste0('Rscript ', PUMAS_path, 'PUMAS.evaluation.customized.R'),
                   paste0('--k ',k), 
                   paste0('--ref_path ', eval_ld_ref),
                   paste0('--trait_name ',trait_name),
                   paste0('--prs_method ', method),
                   paste0('--xty_path ', xty_path),
                   paste0('--stats_path ', stats_path),
                   paste0('--weight_path ', prsdir),
                   paste0('--output_path ', output_path_eval))
  system(tunecode)
  ##### Extract parameter tuning results
  R2.tuned = bigreadr::fread2(paste0(output_path_eval,trait_name, '.', method, '.txt'))
  # Select tuning parameters
  nonzero.ind = which(colMeans(R2.tuned) != 0)
  params.tuned.ct = nonzero.ind[which.max(colMeans(R2.tuned)[nonzero.ind])]
  r2 = params.ct[params.tuned.ct,'r2']; pvalthr = params.ct[params.tuned.ct,'pvalthr']
  r2.ct = r2; pval.ct = pvalthr
  print(paste0('Tuned parameters: r2 = ', r2, ', pval = ', pvalthr))
  tuned.parameters.file = paste0(workdir, 'PRS_model_training/',method,'/tuned_parameters_',trait_name,'.RData')
  save(r2, pvalthr, file = tuned.parameters.file)
}


# Additional filtering for lassosum2 and LDpred2 tuning parameter settings
if ('lassosum2' %in% methods){
  method = 'lassosum2'
  prsdir = paste0(prsdir0, method,'/')
  
  for (ite in 1:k){
    output_lassosum2 = paste0(prsdir, trait_name,'.',method,'.ite',ite,'.txt')
    if(file.exists(output_lassosum2)){
      score = bigreadr::fread2(output_lassosum2) 
      n.tuning = ncol(score) - 4 # as.numeric(strsplit(colnames(score)[ncol(score)],split='_')[[1]][2])
      colnames(score) = c('CHR', 'SNP', 'A1', 'A2', paste0('BETA',1:n.tuning))
    }
    
    # Match alleles with GWAS summary data:
    sumstats = bigreadr::fread2(paste0(output_path, trait_name, '.gwas.ite', ite, '.txt'))
    stateval = sumstats[, c('SNP','A1','A2')]
    colnames(stateval) = c('SNP', 'A1.ref','A2.ref')
    stateval = left_join(stateval,score,by="SNP")
    
    na.ind = which(is.na(stateval$A1))
    if (length(na.ind) > 0){
      stateval[na.ind, paste0('BETA',1:n.tuning)] = 0
      stateval[na.ind,'A1'] = stateval[na.ind,'A1.ref']
      stateval[na.ind,'A2'] = stateval[na.ind,'A2.ref']
    }
    flipped = which(stateval$A1.ref != stateval$A1)
    print(paste0(length(flipped), ' flipped SNPs.'))
    if (length(flipped) > 0){
      stateval[flipped,'A1'] = stateval[flipped,'A1.ref']
      stateval[flipped,'A2'] = stateval[flipped,'A2.ref']
      for (t in 1:n.tuning){
        stateval[flipped,paste0('BETA',t)] = - stateval[flipped,paste0('BETA',t)]
      }
    }
    scores = stateval[,c('CHR','SNP','A1','A2',paste0('BETA',1:n.tuning))] # other files: SNP	CHR	A1	BETA1	BETA2	A2
    write_delim(scores, paste0(input_path, trait_name,'.',method, '.ite',ite,'.txt'), delim = '\t')
    rm(stateval)
  }
  
  xty_path = stats_path = output_path # the "output_path" used forstoring output from pumas.subsampling.R
  
  pumascode = paste(paste0('Rscript ', PUMAS_path, 'PUMAS.evaluation.customized.R'),
                    paste0('--k ',k), 
                    paste0('--ref_path ', eval_ld_ref),
                    paste0('--trait_name ',trait_name),
                    paste0('--prs_method ',method),
                    paste0('--xty_path ', xty_path),
                    paste0('--stats_path ', stats_path),
                    paste0('--weight_path ', input_path),
                    paste0('--output_path ', output_path_eval))
  system(pumascode)
  
  
  ##### Extract parameter tuning results
  r2 = bigreadr::fread2(paste0(output_path_eval,trait_name, '.', method, '.txt'))
  r2.avg = colMeans(r2)
  r2.order = list()
  n.candidates = 20
  for (kk in 1:k) r2.order[[kk]] = order(as.numeric(r2[kk,]),decreasing = T)[1:n.candidates]
  # Select the top parameter settings:
  # params.tuned = as.numeric(substr(names(sort(r2.avg, decreasing = T)[1:5]), 5, 10))
  params.tuned = Reduce(intersect, r2.order)
  nonzero.indx = which(r2.avg != 0)
  params.tuned = unique(params.tuned[(params.tuned <= length(r2.avg)) & (params.tuned %in% nonzero.indx)])
  # print(paste0('Maximum r2 of ',method,': ', max(r2.avg)))
  tuned.parameters.file = paste0(workdir, 'PRS_model_training/',method,'/tuned_parameters_',trait_name,'.RData')
  save(params.tuned, file = tuned.parameters.file)
  if (length(params.tuned) == 0){
    r2 = bigreadr::fread2(paste0(output_path_eval,trait_name, '.', method, '.txt'))
    r2.avg = colMeans(r2)
    nonzero.indx = which(r2.avg != 0)
    r2.order = list()
    n.candidates = length(nonzero.indx)
    for (kk in 1:k) r2.order[[kk]] = order(as.numeric(r2[kk, nonzero.indx]),decreasing = T)[1:n.candidates]
    # Select the top parameter settings:
    # params.tuned = as.numeric(substr(names(sort(r2.avg, decreasing = T)[1:5]), 5, 10))
    params.tuned = Reduce(intersect, r2.order)
    params.tuned = nonzero.indx[params.tuned[!is.na(params.tuned)]]
    params.tuned = unique(params.tuned[(params.tuned <= length(r2.avg)) & (params.tuned %in% nonzero.indx)])
    save(params.tuned, file = tuned.parameters.file)
  }
}


if ('LDpred2' %in% methods){
  method = 'LDpred2'
  prsdir = paste0(prsdir0, method,'/')
  
  for (ite in 1:k){
    output_LDpred2 = paste0(prsdir, trait_name,'.',method,'.ite',ite,'.txt')
    if(file.exists(output_LDpred2)){
      score = bigreadr::fread2(output_LDpred2) 
      n.tuning = ncol(score) - 4 # as.numeric(strsplit(colnames(score)[ncol(score)],split='_')[[1]][2])
      colnames(score) = c('CHR', 'SNP', 'A1', 'A2', paste0('BETA',1:n.tuning))
    }
    
    # Match alleles with GWAS summary data:
    stateval = bigreadr::fread2(paste0(output_path, trait_name, '.gwas.ite', ite, '.txt'))
    stateval = stateval[, c('SNP','A1','A2')]
    colnames(stateval) = c('SNP', 'A1.ref','A2.ref')
    stateval = left_join(stateval,score,by="SNP")
    
    na.ind = which(is.na(stateval$A1))
    if (length(na.ind) > 0){
      stateval[na.ind, paste0('BETA',1:n.tuning)] = 0
      stateval[na.ind,'A1'] = stateval[na.ind,'A1.ref']
      stateval[na.ind,'A2'] = stateval[na.ind,'A2.ref']
    }
    flipped = which(stateval$A1.ref != stateval$A1) # April 14: 0 (already corrected)
    print(paste0(length(flipped), ' flipped SNPs.'))
    if (length(flipped) > 0){
      stateval[flipped,'A1'] = stateval[flipped,'A1.ref']
      stateval[flipped,'A2'] = stateval[flipped,'A2.ref']
      for (t in 1:n.tuning){
        stateval[flipped,paste0('BETA',t)] = - stateval[flipped,paste0('BETA',t)]
      }
    }
    scores = stateval[,c('CHR','SNP','A1','A2',paste0('BETA',1:n.tuning))] # other files: SNP	CHR	A1	BETA1	BETA2	A2
    write_delim(scores, paste0(input_path, trait_name,'.',method, '.ite',ite,'.txt'), delim = '\t')
    # write.table(scores, paste0(input_path, trait_name,'.',method, '.ite', ite, '.txt'), row.names = F,col.names = T, quote = FALSE, sep = "\t" )
    rm(stateval, scores)
  }
  
  xty_path = stats_path = output_path # the "output_path" used forstoring output from pumas.subsampling.R
  
  pumascode = paste(paste0('Rscript ', PUMAS_path, 'PUMAS.evaluation.customized.R'),
                    paste0('--k ',k), 
                    paste0('--ref_path ', eval_ld_ref),
                    paste0('--trait_name ',trait_name),
                    paste0('--prs_method ',method),
                    paste0('--xty_path ', xty_path),
                    paste0('--stats_path ', stats_path),
                    paste0('--weight_path ', input_path),
                    paste0('--output_path ', output_path_eval))
  system(pumascode)
  
  ##### Extract parameter tuning results
  r2 = bigreadr::fread2(paste0(output_path_eval,trait_name, '.', method, '.txt'))
  r2.avg = colMeans(r2)
  r2.order = list()
  n.candidates = 20
  for (kk in 1:k) r2.order[[kk]] = order(as.numeric(r2[kk,]),decreasing = T)[1:n.candidates]
  # Select the top parameter settings:
  # params.tuned = as.numeric(substr(names(sort(r2.avg, decreasing = T)[1:5]), 5, 10))
  params.tuned = Reduce(intersect, r2.order)
  nonzero.indx = which(r2.avg != 0)
  params.tuned = unique(params.tuned[(params.tuned <= length(r2.avg)) & (params.tuned %in% nonzero.indx)])
  # print(paste0('Maximum r2 of ',method,': ', max(r2.avg)))
  tuned.parameters.file = paste0(workdir,'PRS_model_training/',method,'/tuned_parameters_',trait_name,'.RData')
  save(params.tuned, file = tuned.parameters.file)
  if (length(params.tuned) == 0){
    r2 = bigreadr::fread2(paste0(output_path_eval,trait_name, '.', method, '.txt'))
    r2.avg = colMeans(r2)
    nonzero.indx = which(r2.avg != 0)
    r2.order = list()
    n.candidates = length(nonzero.indx)
    for (kk in 1:k) r2.order[[kk]] = order(as.numeric(r2[kk, nonzero.indx]),decreasing = T)[1:n.candidates]
    # Select the top parameter settings:
    # params.tuned = as.numeric(substr(names(sort(r2.avg, decreasing = T)[1:5]), 5, 10))
    params.tuned = Reduce(intersect, r2.order)
    params.tuned = nonzero.indx[params.tuned[!is.na(params.tuned)]]
    params.tuned = unique(params.tuned[(params.tuned <= length(r2.avg)) & (params.tuned %in% nonzero.indx)])
    save(params.tuned, file = tuned.parameters.file)
  }
}



# ---------------------------------------------------------------------------------------
# --------------------- Step 4: Generate best PRS(s) for each method --------------------
# ---------------------------------------------------------------------------------------
# Once we select the "best" PRSs for different models, train on the whole dataset
# For LDpred2 and lassosum, select the best ones that are not all 0 or NA's
# Then come back and train the best PRS model
if ( opt$verbose >= 1 ) {
  cat(paste0("\n************************************************************"))
  cat(paste0("\n******* Step 4: Generate best PRS(s) for each method *******"))
  cat(paste0("\n************************************************************\n"))
}
# ---------------------------------------------------------------------------------------
# ---------- Step 4.1: Train the best PRS(s) from each method on the whole data ---------
# ---------------------------------------------------------------------------------------

sumraw = bigreadr::fread2(paste0(output_path,trait_name,".gwas_matched.txt"))
err.CT = err.lassosum2 = err.LDpred2 = 0

# C+T: find the one best PRS 
if ('C+T' %in% methods){
  method = 'C+T'
  prsdir = paste0(prsdir0, method,'/')
  if ( opt$verbose >= 1 ){
    print(paste0('************************************************************'))
    print(paste0('******* Start training ', method, ' on the original GWAS data *******'))
    print(paste0('************************************************************'))
  }
  # The default number of tuning parameter settings is 12, so NCORES can be 3 (default), 6, or 12. 
  # Parameters for the Clumping step
  z = qnorm(p=p.ldclump/2,lower.tail=FALSE)
  Pvalthr = as.numeric(str_split(opt$Pvalthr,",")[[1]])
  R2 = as.numeric(str_split(opt$R2,",")[[1]])
  params.ct = expand.grid(r2 = R2, pvalthr = Pvalthr)
  
  # --------------------- Step 2.1: Preparation for input files for C+T ---------------------
  sumstats = sumraw[,c('CHR','SNP','A1','A2','BETA','SE','P','N', 'MAF')]
  sumstats$P = as.numeric(sumstats$P)
  sumstats$BETA = as.numeric(sumstats$BETA)
  sumstats$SE = as.numeric(sumstats$SE)
  sumstats$N = as.numeric(sumstats$N)
  sumstats$MAF = as.numeric(sumstats$MAF)
  colnames(sumstats) <- c("CHR", "SNP", "REF", "ALT", "BETA", "SE", "P", "N", "FRQ") # MAF or REF_FRQ
  rownames(sumstats) = sumstats$SNP
  sumstats0 = sumstats
  # REF: effect allele
  print(paste0('** Complete Loading GWAS summary data **'))
  
  tuned.parameters.file = paste0(workdir, 'PRS_model_training/',method,'/tuned_parameters_',trait_name,'.RData')
  load(tuned.parameters.file) # Load tuned r2 and pvalthr
  
  
  # Create base data (summary statistic) file containing the P-value information:
  sumstats_input <- sumstats[,c('SNP','P')] #,'REF_FRQ'
  fwrite(sumstats_input, file=paste0(prsdir,'snplist/',trait_name,'.txt'),row.names = F, quote = F, sep=' ')
  
  # --------------------- Step 2.2: Run C+T ----------------------------
  set.seed(2023)
  # --------------------- LD Clumping (C) ---------------------
  ldclumpcode <- paste0(plink_path, 'plink --bfile ', eval_ld_ref,
                        ' --clump ',prsdir,'snplist/',trait_name,'.txt',
                        ' --clump-p1 ',p.ldclump,
                        # ' --clump-p2 ',pc,
                        ' --clump-r2 ',r2,
                        ' --clump-kb ',kb,
                        ' --threads ', threads,
                        ' --silent',
                        ' --out ',prsdir,'clumped/',trait_name,'_r2=',r2)
  system(ldclumpcode)
  # --------------------- Thresholding (T) ---------------------
  LD <- bigreadr::fread2(paste0(prsdir,'clumped/',trait_name,'_r2=',r2,'.clumped'))
  clumped.snp <- LD[,3,drop=F][,1]
  sumstats.clumped <- sumstats[clumped.snp,]
  
  keep.SNP = sumstats.clumped[sumstats.clumped$P <= pvalthr,c('SNP')]
  dump.SNP = which(!sumstats$SNP %in% keep.SNP)
  sumstats[dump.SNP, 'BETA'] = 0
  
  SCORE = sumstats[,c('CHR','SNP','REF','ALT', 'BETA')]
  colnames(SCORE)[3:4] = c('A1','A2')
  # Remove SNPs with zero effects:
  beta_ct.out = SCORE[SCORE$BETA != 0,]
  
  ### If all effects are 0:
  R2.tuned = bigreadr::fread2(paste0(output_path_eval,trait_name, '.', method, '.txt'))
  nonzero.ind = which(colMeans(R2.tuned) != 0)
  params.tuned.ct = nonzero.ind[which.max(colMeans(R2.tuned)[nonzero.ind])]
  while((nrow(beta_ct.out) == 0) & (length(nonzero.ind) > 0)){
    sumstats = sumstats0
    nonzero.ind = nonzero.ind[-which(nonzero.ind == params.tuned.ct)]
    params.tuned.ct = nonzero.ind[which.max(colMeans(R2.tuned)[nonzero.ind])]
    r2 = params.ct[params.tuned.ct,'r2']; pvalthr = params.ct[params.tuned.ct,'pvalthr']
    keep.SNP = sumstats.clumped[sumstats.clumped$P <= pvalthr,c('SNP')]
    dump.SNP = which(!sumstats$SNP %in% keep.SNP)
    sumstats[dump.SNP, 'BETA'] = 0
    SCORE = sumstats[,c('CHR','SNP','REF','ALT', 'BETA')]
    colnames(SCORE)[3:4] = c('A1','A2')
    # Remove SNPs with zero effects:
    beta_ct.out = SCORE[SCORE$BETA != 0,]
  }
  if (nrow(beta_ct.out) > 0){
    prs_ct_outputfile = paste0(PennPRS_finalresults_path, trait_name,'.',method,'.PRS.txt')
    write_delim(beta_ct.out, prs_ct_outputfile, delim = '\t')
    
    ct.optimal.indx = which((params.ct$r2 == r2) & (params.ct$pvalthr == pvalthr))
    optimal.indx = ct.optimal.indx
    save(optimal.indx, file = paste0(prsdir, trait_name, '.', method, '.optimal.indx.RData'))
    rm(optimal.indx)
    write_delim(params.ct[ct.optimal.indx,], paste0(PennPRS_finalresults_path, trait_name,'.',method,'.optimal_params.txt'), delim = '\t')
    if ( opt$verbose >= 1 ){
      print(paste0('*************************************'))
      print(paste0('******* Complete training ', method, ' *******'))
      print(paste0('*************************************'))
    }
    # If R2<0:
    R2.tuned = bigreadr::fread2(paste0(output_path_eval,trait_name, '.', method, '.txt'))
    if (colMeans(R2.tuned[ct.optimal.indx]) <= 0){
      err.CT = 1
      if (ensemble){
        ensemble.methods = ensemble.methods[which(ensemble.methods != method)]
        print(paste0('[Warning] Trained PRS model based on ', method, ' has an estimated R2 < 0, which indicates that the PRS lacks prediction power and thus we will not incorporate it in the ensemble PRS. \nPotential explanations for R2 < 0 are:\n 1. The trait is not heritable.\n 2. The GWAS have insufficient power (e.g., due to low sample size) to develop a predictive PRS.\n 3. Issues with the input GWAS summary data (e.g., problematic BETA or SE).\n 4. ', method, ' is not powerful for developing PRS for the trait, in which case other methods can be considered.'))
      } 
      if (!ensemble){
        print(paste0('[Warning] Trained PRS model based on ', method, ' has an estimated R2 < 0, which indicates that the PRS lacks prediction power. \nPotential explanations for R2 < 0 are:\n 1. The trait is not heritable.\n 2. The GWAS have insufficient power (e.g., due to low sample size) to develop a predictive PRS.\n 3. Issues with the input GWAS summary data (e.g., problematic BETA or SE).\n 4. ', method, ' is not powerful for developing PRS for the trait, in which case other methods can be considered.'))
      }
    }
  }
  if (nrow(beta_ct.out) == 0){
    err.CT = 2
    if (ensemble){
      ensemble.methods = ensemble.methods[which(ensemble.methods != method)]
      print(paste0('Trained PRS model based on ', method, ' has zero effect estimate for all SNPs and thus we will not incorporate it in the ensemble PRS. Please try other methods.'))
    } 
    if (!ensemble){
      print(paste0('Trained PRS model based on ', method, ' has zero effect estimate for all SNPs. Please try other methods.'))
    }
  }
  rm(SCORE)
}




if (('lassosum2' %in% methods) | ('LDpred2' %in% methods)){
  map_ldref <- readRDS(paste0(ld_path, 'map/map_',LDrefpanel,'_ldref.rds'))
  
  sumstats = sumraw[,c('CHR','SNP','A1','A2','BETA','SE','P','N', 'MAF')]
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
  

  
  (ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2,
                                  sample_size = n_eff, blocks = NULL)))
  ldsc_h2_est <- abs(ldsc[["h2"]])
  cat(paste0('Heritability estimate based on LD score regression: ', signif(ldsc[["h2"]], 3)))
  if (ldsc[["h2"]] < 0) cat(paste0('Warning: negative hertability estimate based on LD score regression.'))
  
  if ('lassosum2' %in% methods){
    method = 'lassosum2'
    prsdir = paste0(prsdir0, method,'/')
    if ( opt$verbose >= 1 ){
      print(paste0('******************************************************************'))
      print(paste0('******* Start training lassosum2 on the original GWAS data *******'))
      print(paste0('******************************************************************'))
    }
    delta = as.numeric(str_split(opt$delta,",")[[1]]) # candidate values of the shrinkage parameter in L2 regularization
    params.lassosum2.grid = expand.grid(lambda = 1:nlambda, delta = delta)
    
    # Load tuned parameters (top 3 choices were saved, in case the top choice gives NA)
    r2 = bigreadr::fread2(paste0(output_path_eval,trait_name, '.', method, '.txt'))
    r2.avg = colMeans(r2)
    tuned.parameters.file = paste0(workdir, 'PRS_model_training/',method,'/tuned_parameters_',trait_name,'.RData')
    load(tuned.parameters.file)
    
    if (length(params.tuned) == 0) cat(paste0('Trained PRS model based on ', method, ' has zero effect estimate for all SNPs. Please try other methods.'))
    if (length(params.tuned) > 0){
      params.lassosum2 = bigreadr::fread2(paste0(prsdir, 'params.lassosum2.txt')) 
      delta = sort(unique(params.lassosum2[params.tuned, 'delta'])) # need to align the order of delta with original order
      params.lassosum2.reduced.grid0 = params.lassosum2.grid[params.tuned,]
      params.lassosum2.reduced.grid = expand.grid(lambda = 1:nlambda, delta = delta)
      # which of these restricted tuning parameter settings are in the top candidate list
      train.indx.lassosum2 = sapply(1:length(params.tuned), function(x){which((params.lassosum2.reduced.grid[, 'lambda'] == params.lassosum2.reduced.grid0[x,'lambda']) & (params.lassosum2.reduced.grid[, 'delta'] == params.lassosum2.reduced.grid0[x,'delta']))})
      
      # --------------------- Step 2.2: Run lassosum2 ----------------------------
      set.seed(2023)
      beta_lassosum2 <- snp_lassosum2(corr, df_beta, # ncores = NCORES,
                                      delta = delta, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio) 
      params.lassosum2.reduced <- attr(beta_lassosum2, "grid_param")
      beta_lassosum2[is.na(beta_lassosum2)] = 0
      # Columns that have nonzero entries: index in params.tuned and train.indx.lassosum2
      # nonzero.indx = which(sapply(train.indx.lassosum2, function(x){sum(beta_lassosum2[,x]!=0)>0}))
      nonzero.indx = which(sapply(train.indx.lassosum2, function(x){sum(abs(beta_lassosum2[,x])>1e-7)>0}))
      if (length(nonzero.indx) > 0){
        if (length(nonzero.indx) == 1) {
          indx.temp = ii = 1; optimal.indx = params.tuned
        }
        if (length(nonzero.indx) > 1) {
          candidates = params.tuned[nonzero.indx]
          stop = 0; ii = 0
          while ((stop == 0) & (ii < length(candidates))){
            ii = ii + 1
            if (((candidates[ii]-1) %in% candidates) | ((candidates[ii]+1) %in% candidates)){
              stop = 1
            }
          }
          if (stop == 0){
            ii = 1
          }
          optimal.indx = candidates[ii]
          indx.temp = ii
        }
        save(optimal.indx, file = paste0(prsdir, trait_name, '.', method, '.optimal.indx.RData'))
        lassosum2.out = signif(params.lassosum2.reduced[train.indx.lassosum2[nonzero.indx][indx.temp],], 4)
        write_delim(lassosum2.out, paste0(PennPRS_finalresults_path, trait_name,'.',method,'.optimal_params.txt'))
        beta_lassosum2 = data.frame(df_beta[,c('chr','rsid','a1','a0')], beta_lassosum2[,train.indx.lassosum2[nonzero.indx][indx.temp]], delim = '\t')
        colnames(beta_lassosum2) = c('CHR','SNP','A1','A2', 'BETA')
        
        nonzero = (sum(beta_lassosum2$BETA!=0)>0)
        if (nonzero){
          # prs_lassosum2_archivefile = paste0(PennPRS_finalresults_path, trait_name,'.',method,'.PRS.fullversion.txt')
          # write_delim(beta_lassosum2, file = prs_lassosum2_archivefile, delim='\t')
          # Remove SNPs with zero effects:
          beta_lassosum2.out = beta_lassosum2[beta_lassosum2$BETA != 0,]
          prs_lassosum2_outputfile = paste0(PennPRS_finalresults_path, trait_name,'.',method,'.PRS.txt')
          write_delim(beta_lassosum2.out, prs_lassosum2_outputfile, delim = '\t')
          print(paste0('Optimal parameter setting: '))
          print(paste0('Delta: ', signif(params.lassosum2[nonzero.indx[indx.temp],'delta'], digits = 3)))
          print(paste0('Lambda: ', signif(params.lassosum2[nonzero.indx[indx.temp],'lambda'], digits = 3)))
          print(paste0('Complete training ', method))
          
          if (r2.avg[optimal.indx] <= 0){
            err.lassosum2 = 1
            if (ensemble){
              ensemble.methods = ensemble.methods[which(ensemble.methods != method)]
              print(paste0('[Warning] Trained PRS model based on ', method, ' has an estimated R2 < 0, which indicates that the PRS lacks prediction power and thus we will not incorporate it in the ensemble PRS. \nPotential explanations for R2 < 0 are:\n 1. The trait is not heritable.\n 2. The GWAS have insufficient power (e.g., due to low sample size) to develop a predictive PRS.\n 3. Issues with the input GWAS summary data (e.g., problematic BETA or SE).\n 4. ', method, ' is not powerful for developing PRS for the trait, in which case other methods can be considered.'))
            } 
            if (!ensemble){
              print(paste0('[Warning] Trained PRS model based on ', method, ' has an estimated R2 < 0, which indicates that the PRS lacks prediction power. \nPotential explanations for R2 < 0 are:\n 1. The trait is not heritable.\n 2. The GWAS have insufficient power (e.g., due to low sample size) to develop a predictive PRS.\n 3. Issues with the input GWAS summary data (e.g., problematic BETA or SE).\n 4. ', method, ' is not powerful for developing PRS for the trait, in which case other methods can be considered.'))
            }
          }
        }
        if (!nonzero){
          err.lassosum2 = 1
          if (ensemble){
            ensemble.methods = ensemble.methods[which(ensemble.methods != method)]
            print(paste0('Trained PRS model based on ', method, ' has zero effect estimate for all SNPs (potentially due to convergence issues) and thus we will not incorporate it in the ensemble PRS. Please consider other methods.'))
          } 
          if (!ensemble){
            print(paste0('Trained PRS model based on ', method, ' has zero effect estimate for all SNPs (potentially due to convergence issues). Please consider other methods.'))
          }
        }
      }
      if (length(nonzero.indx) == 0){
        err.lassosum2 = 0
        if (ensemble){
          ensemble.methods = ensemble.methods[which(ensemble.methods != method)]
          print(paste0('Trained PRS model based on ', method, ' has zero effect estimate for all SNPs (potentially due to convergence issues) and thus we will not incorporate it in the ensemble PRS. Please consider other methods.'))
        } 
        if (!ensemble){
          print(paste0('Trained PRS model based on ', method, ' has zero effect estimate for all SNPs (potentially due to convergence issues). Please consider other methods.'))
        }      
      }
      if ( opt$verbose >= 1 ){
        print(paste0('*************************************'))
        print(paste0('******* Complete training ', method, ' *******'))
        print(paste0('*************************************'))
      }
      suppressWarnings(rm(optimal.indx))
    }
  }
  
  
  if ('LDpred2' %in% methods){
    method = 'LDpred2'
    prsdir = paste0(prsdir0, method,'/')
    
    r2 = bigreadr::fread2(paste0(output_path_eval,trait_name, '.', method, '.txt'))
    r2.avg = colMeans(r2)
    tuned.parameters.file = paste0(workdir,'PRS_model_training/',method,'/tuned_parameters_',trait_name,'.RData')
    load(tuned.parameters.file)
    
    if (length(params.tuned) == 0) cat(paste0('Trained PRS model based on ', method, ' has zero effect estimate for all SNPs. Please try other methods.'))
    if (length(params.tuned) > 0){
      h2_seq <- round(ldsc_h2_est * h2.ratio, 5); 
      h2_seq[h2_seq == 0] = 1e-5
      h2_seq[duplicated(h2_seq)] = h2_seq[duplicated(h2_seq)] * 1.01
      n.inflated = sum(h2_seq>1)
      if (n.inflated > 0) h2_seq[h2_seq>1] = 0.95 + seq(0,0.01*(n.inflated-1), by = 0.01)
      params.ldpred2.train <- expand.grid(p = p_seq, h2 = h2_seq, sparse = sparse.option)
      params.ldpred2.reduced = params.ldpred2.train[params.tuned,]
      
      set.seed(2023)
      beta_ldpred2 <- snp_ldpred2_grid(corr, df_beta, params.ldpred2.reduced, ncores = NCORES)
      beta_ldpred2[is.na(beta_ldpred2)] = 0
      
      # Columns that have nonzero entries:
      # nonzero.indx = which(sapply(1:ncol(beta_ldpred2), function(x){sum(beta_ldpred2[,x]!=0)>0}))
      nonzero.indx = which(sapply(1:ncol(beta_ldpred2), function(x){sum(abs(beta_ldpred2[,x])>1e-10)>0}))
      if (length(nonzero.indx) > 0){
        if (length(nonzero.indx) == 1) {
          ii = 1; optimal.indx = params.tuned
        }
        if (length(nonzero.indx) > 1) {
          candidates = params.tuned[nonzero.indx]
          stop = 0; ii = 0
          while ((stop == 0) & (ii < length(candidates))){
            ii = ii + 1
            if (((candidates[ii]-1) %in% candidates) | ((candidates[ii]+1) %in% candidates)){
              stop = 1
            }
          }
          if (stop == 0){
            ii = 1
          }
          optimal.indx = candidates[ii]
          indx.temp = ii
        }
        # indx.temp = which.max(r2.avg[params.tuned[nonzero.indx]])
        # ldpred2.optimal.indx = params.tuned[nonzero.indx][indx.temp]
        # optimal.indx = ldpred2.optimal.indx
        save(optimal.indx, file = paste0(prsdir, trait_name, '.', method, '.optimal.indx.RData'))
        # params.ldpred2.train.simplified = expand.grid(p = p_seq, h2 = h2.ratio, sparse = sparse.option)
        # write_delim(params.ldpred2.train.simplified[ldpred2.optimal.indx,], 
        #             file = paste0(PennPRS_finalresults_path, trait_name,'.',method,'.optimal_params.txt'))
        ldpred2.out = params.ldpred2.train[optimal.indx,]
        write_delim(ldpred2.out, paste0(PennPRS_finalresults_path, trait_name,'.',method,'.optimal_params.txt'))
        beta_ldpred2 = data.frame(df_beta[,c('chr','rsid','a1','a0')], beta_ldpred2[,nonzero.indx[indx.temp]])
        colnames(beta_ldpred2) = c('CHR','SNP','A1','A2', 'BETA')
        
        nonzero = (sum(beta_ldpred2[,'BETA']!=0)>0)
        if (nonzero){
          # prs_ldpred2_archivefile = paste0(PennPRS_finalresults_path, trait_name,'.',method,'.PRS.fullversion.txt')
          # write_delim(beta_ldpred2, file = prs_ldpred2_archivefile, delim='\t')
          # Remove SNPs with zero effects:
          beta_ldpred2.out = beta_ldpred2[beta_ldpred2$BETA != 0,]
          prs_ldpred2_outputfile = paste0(PennPRS_finalresults_path, trait_name,'.',method,'.PRS.txt')
          write_delim(beta_ldpred2.out, prs_ldpred2_outputfile)
          print(paste0('Optimal parameter setting: '))
          print(paste0('p: ', signif(params.ldpred2.reduced[nonzero.indx[indx.temp],'p'], digits = 3)))
          print(paste0('h2: ', signif(params.ldpred2.reduced[nonzero.indx[indx.temp],'h2'], digits = 3)))
          print(paste0('sparse: ', signif(params.ldpred2.reduced[nonzero.indx[indx.temp],'sparse'], digits = 3)))
          print(paste0('Complete training ', method))
          
          if (r2.avg[optimal.indx] <= 0){
            err.LDpred2 = 1
            if (ensemble){
              ensemble.methods = ensemble.methods[which(ensemble.methods != method)]
              print(paste0('[Warning] Trained PRS model based on ', method, ' has an estimated R2 < 0, which indicates that the PRS lacks prediction power and thus we will not incorporate it in the ensemble PRS. \nPotential explanations for R2 < 0 are:\n 1. The trait is not heritable.\n 2. The GWAS have insufficient power (e.g., due to low sample size) to develop a predictive PRS.\n 3. Issues with the input GWAS summary data (e.g., problematic BETA or SE).\n 4. ', method, ' is not powerful for developing PRS for the trait, in which case other methods can be considered.'))
            } 
            if (!ensemble){
              print(paste0('[Warning] Trained PRS model based on ', method, ' has an estimated R2 < 0, which indicates that the PRS lacks prediction power. \nPotential explanations for R2 < 0 are:\n 1. The trait is not heritable.\n 2. The GWAS have insufficient power (e.g., due to low sample size) to develop a predictive PRS.\n 3. Issues with the input GWAS summary data (e.g., problematic BETA or SE).\n 4. ', method, ' is not powerful for developing PRS for the trait, in which case other methods can be considered.'))
            }
          }
        }
        if (!nonzero){
          err.LDpred2 = 2
          if (ensemble){
            ensemble.methods = ensemble.methods[which(ensemble.methods != method)]
            print(paste0('Trained PRS model based on ', method, ' has zero effect estimate for all SNPs (potentially due to convergence issues) and thus we will not incorporate it in the ensemble PRS. Please consider other methods.'))
          } 
          if (!ensemble){
            print(paste0('Trained PRS model based on ', method, ' has zero effect estimate for all SNPs (potentially due to convergence issues). Please consider other methods.'))
          }
        }
      }
      if (length(nonzero.indx) == 0){
        err.LDpred2 = 2
        if (ensemble){
          ensemble.methods = ensemble.methods[which(ensemble.methods != method)]
          print(paste0('Trained PRS model based on ', method, ' has zero effect estimate for all SNPs (potentially due to convergence issues) and thus we will not incorporate it in the ensemble PRS. Please consider other methods.'))
        } 
        if (!ensemble){
          print(paste0('Trained PRS model based on ', method, ' has zero effect estimate for all SNPs (potentially due to convergence issues). Please consider other methods.'))
        }      
      }
      if (opt$verbose >= 1 ){
        print(paste0('*************************************'))
        print(paste0('******* Complete training ', method, ' *******'))
        print(paste0('*************************************'))
      }
    }
  }
  rm(corr, ld)
  # system(paste0('rm -rf ', paste0(prsdir0, 'temporary_LDpred2_lassosum2_ite1')))
}


if (ensemble) save(ensemble.methods, err.CT, err.lassosum2, err.LDpred2, file = paste0(PennPRS_finalresults_path, 'step1.RData'))
if (!ensemble) save(err.CT, err.lassosum2, err.LDpred2, file = paste0(PennPRS_finalresults_path, 'step1.RData'))


