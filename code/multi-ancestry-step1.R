if(!require(optparse)){
  install.packages("optparse")
  library(optparse)
}
if(!require(parallel)){
  install.packages("parallel")
  library(parallel)
}
if(!require(readr)){
  install.packages("readr")
  library(readr)
}
if(!require(bigreadr)){
  install.packages("bigreadr")
  library(bigreadr)
}
if(!require(bigsnpr)){
  install.packages("bigsnpr")
  library(bigsnpr)
}
if(!require(data.table)){
  install.packages("data.table")
  library(data.table)
}
if(!require(dplyr)){
  install.packages("dplyr")
  library(dplyr)
}
if(!require(scales)){
  install.packages("scales")
  library(scales)
}
if(!require(stringr)){# for str_split
  install.packages("stringr")
  library(stringr)
}
if(!require(caret)){# for findCorrelation to train ensemble PRS
  install.packages("caret")
  library(caret)
}
if(!require(inline)){
  install.packages("inline")
  library(inline)
}
if(!require(devtools)){
  install.packages("devtools")
  library(devtools)
}
if(!require(genio)){
  install.packages("genio")
  library(genio)
}
# if(!require(pryr)){
#   install.packages("pryr")
#   library(pryr)
# }
if(!require(Matrix)){
  install.packages("Matrix")
  library(Matrix)
}
if(!require(lavaan)){
  install.packages("lavaan")
  library(lavaan)
}
if(!require(bigparallelr)){
  install.packages("bigparallelr")
  library(bigparallelr)
}
if(!require(bigmemory)){
  install.packages("bigmemory")
  library(bigmemory)
}
if(!require(Rcpp)){
  install.packages("Rcpp")
  library(Rcpp)
}
if(!require(RcppArmadillo)){
  install.packages("RcppArmadillo")
  library(RcppArmadillo)
}
if(!require(RcppTN)){
  install.packages("RcppTN")
  library(RcppTN)
}
if(!require(doMC)){
  install.packages("doMC")
  library(doMC)
}
if(!require(foreach)){
  install.packages("foreach")
  library(foreach)
}
if(!require(doParallel)){
  install.packages("doParallel")
  library(doParallel)
}
if(!require(MASS)){
  install.packages("MASS")
  library(MASS)
}
if(!require(reshape)){
  install.packages("reshape")
  library(reshape)
}
if(!require(RISCA)){
  install.packages("RISCA")
  library(RISCA)
}


options(stringsAsFactors=F)
option_list = list(
  make_option("--homedir", action = "store", default = NA, type = "character",
              help="Path to save the output folder [Required]"),
  make_option("--PennPRS_path", action = "store", default = NA, type = "character",
              help="Path to the PennPRS folder [Required]"),
  make_option("--input_GWAS_path", action = "store", default = NA, type = "character",
              help="gwas path after qc"),
  make_option("--submissionID", action = "store", default = NA, type = "character",
              help="Job ID [Required]"),
  make_option("--methods", action = "store", default = 'PROSPER', type = "character",
              help="Options: a subset of methods from PROSPER, PRS-CSx, and MUSSEL, divided by comma"),
  make_option("--trait", action = "store", default = NA, type = "character",
              help="trait name [Optional]"),
  make_option("--races", action = "store", default = NA, type = "character",
              help="Races of the training GWAS data. Options: a subset (need to have at least two) of: EUR (European), AFR (African),
              AMR (Mixed American, Hispanic/Latio), EAS (East Asian), or SAS (South Asian), divided by comma [Required]"),
  make_option("--LDrefpanel", action = "store", default = '1kg', type = "character",
              help="LD reference panel. Options: '1kg' (1000 Genomes Project Phase 3) or 'ukbb' (UK Biobank) [Optional]"),
  make_option("--k", action = "store", default = 2, type = "numeric",
              help = "k-fold Monte Carlo Cross Validation (MCCV) for PUMAS. Options: any integer greater than or equal to 2 [Optional]"),
  
  make_option("--partitions", action = "store", default = '0.8,0.2', type = "character",
              help="Partitions for PUMAS subsampling. 
              Format: '% training, % testing' (% testing is equal to 1 - % training), divided by comma [Optional]"),
  
  # Parameters in PROSPER
  make_option("--ndelta", action = "store", default = 5, type = "integer",
              help="Number of candidate values of the shrinkage parameter in L2 regularization. Options: 
              any positive integer number [Optional]"),
  make_option("--nlambda", action = "store", default = 5, type = "integer",
              help="Number of different candidate values for lambda (shrinkage parameter in the L1 regularization.
              Options: any positive integer [Optional]"),
  make_option("--lambda.min.ratio", action = "store", default = 0.01, type = "numeric",
              help="Ratio between the lowest and highest candidate values of lambda.
              Options: any value in (0,1) [Optional]"),
  make_option("--Ll", action="store", default=5, type='integer',
              help="Length of path for the tuning parameter lambda in the PROSPER step [default: %default]"),
  make_option("--Lc", action="store", default=5, type='integer',
              help="Length of path for the tuning parameter c in the PROSPER step [default: %default]"),
  
  # Parameters in PRS-CSx
  make_option("--phi", action = "store", default = '1e-02,1e-04', type = "character",
              help="Global shrinkage parameter phi. For GWAS with limited sample sizes (including most of the current disease GWAS), 
              fixing phi to 1e-2 (for highly polygenic traits) or 1e-4 (for less polygenic traits), or doing a small-scale grid search 
              (e.g., phi=1e-6, 1e-4, 1e-2, 1) to find the optimal phi value in the validation dataset often improves perdictive performance.
              Alternatively, phi can be learnt from the data using a fully Bayesian approach. This works well for polygenic traits 
              with very large GWAS sample sizes (hundreds of thousands of subjects), but is computationally too intensive and thus is not included in our pipeline.
              Default (default values in the PRS-CS algorithm, Nov 21, 2024 version): grid search with candidate values 1e+00,1e-02,1e-04,1e-6.
              Options: candidate values in (0,1], divided by comma [Optional]"),
  
  # Parameters in MUSSEL
  make_option("--p", action="store", default=paste(signif(seq_log(1e-5, 1, length.out = 11), 2), collapse = ','), type='character',
              help="Candidate values for tuning parameter p (causal SNP proportion) [default: %default]"),
  make_option("--H2", action="store", default=paste(c(0.7, 1, 1.4), collapse = ','), type='character',
              help="Candidate values for tuning parameter H2.ratio (heritability = H2.ratio * h2_est from LDSC) [default: %default]"),
  make_option("--sparse", action="store", default='0', type='character',
              help="Whether to consider a sparse model: 0, 1, or 0,1 [default: %default]"),
  make_option("--cors_additional", action="store", default=NA, type='character',
              help="Additional candidate values for tuning parameter: genetic correlation across ancestry groups, example: 3 groups with label 1,2,3, want to add two additional settings: cor_setting1(1,2),cor_setting1(1,3),cor_setting1(2,3);cor_setting2(1,2),cor_setting2(1,3),cor_setting2(2,3) [optional]"),
  make_option("--ps_additional", action="store", default=NA, type='character',
              help="Typically not necessary. Additional candidate values for tuning parameter: ancestry-specific causal SNP proportions, example: 3 groups with label 1,2,3, want to add two additional settings: p1_setting1,p2_setting1,p3_setting1;p1_setting2,p2_setting2,p3_setting2 [optional]"),
  
  make_option("--NCORES", action = "store", default = 5, type = "numeric",
              help="Number of cores (default: 5) used for parallel computing of PROSPER and MUSSEL
              Options: positive integer [Optional]"),
  make_option("--N_THREADS", action = "store", default = 5, type = "numeric",
              help="Number of threads (default: 5) used for parallel computing of PRS-CSx across chromosomes.
              Options: positive integer [Optional]"),
  
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
races = str_split(opt$races,",")[[1]]; K = length(races)
LDrefpanel = opt$LDrefpanel
# Parameters for subsampling
k = opt$k
NCORES = opt$NCORES
N_THREADS = opt$N_THREADS
# ----------------

# Optional input parameters (for PUMAS subsampling):
partitions <- opt$partitions

# ----------------
PUMAS_path = paste0(PennPRS_path,'/code/')
PROSPER_path = paste0(PennPRS_path, '/software/PROSPER')
path_plink = '/dcl01/chatterj/data/jin/software/plink2'
PRScs_path = paste0(PennPRS_path, '/software/PRScs/')
PRScsx_path = paste0(PennPRS_path, '/software/PRScsx/')
MUSSEL_path = paste0(PennPRS_path, '/software/MUSSEL/')
threads = 1

ld_path <- paste0(PennPRS_path, '/LD/', races, '/')
ld_path0 <- paste0(ld_path, 'LD_1kg/') # set to the /LD_1kg folder under /LD/
if (LDrefpanel == '1kg'){
  eval_ld_ref_path <- paste0(ld_path, '/1KGref_plinkfile/') # set to the /1KGref_plinkfile folder under /LD/
  path_precalLD <- paste0(ld_path, '/LDpred2_lassosum2_corr_1kg/') # set to the /LDpred2_lassosum2_corr_1kg folder under /LD/
} 
names(ld_path0) = names(ld_path) = names(eval_ld_ref_path) = names(path_precalLD) = races
# Job name/ID: e.g., trait_race_method_userID_submissionID
jobID = paste(c(trait,paste0(races,collapse = '.'),paste0(methods,collapse = '.'), submissionID), collapse = '_')
# Create a job-specific (trait, race, methods, userID, jobID) directory to save all the outputs, set the working directory to this directory
workdir = paste0(homedir,jobID,'/')
suppressWarnings(dir.create(workdir))
setwd(workdir) 

# ----------------
if ('PROSPER' %in% methods){
  ndelta = opt$ndelta # number of candidate values of the shrinkage parameter in L2 regularization
  nlambda = opt$nlambda # number of different candidate values for lambda (shrinkage parameter in the L1 regularization). Default in lassosum2 pipeline: 30, which may lead to issues when using PUMAS subsampling to tune parameters
  lambda.min.ratio = opt$lambda.min.ratio # Ratio between the lowest and highest candidate values of lambda. Dandidate values in (0,Inf), divided by comma
  Ll = opt$Ll
  Lc = opt$Lc
}
if ('PRS-CSx' %in% methods){
  PHI = opt$phi; phi.vals = as.numeric(str_split(PHI,",")[[1]])
  N_THREADS = opt$N_THREADS
}
if ('MUSSEL' %in% methods){
  # Parameters in the LDpred2 step
  p <- opt$p
  H2 = opt$H2 # as.numeric(strsplit(opt$H2, split = ',')[[1]])
  sparse = opt$sparse
  # Parameters in the MUSS step
  cors_additional = opt$cors_additional
  ps_additional = opt$ps_additional
}

source(paste0(PUMAS_path, 'PennPRS_functions.R')) # please save the PennPRS_functions.R file to the /PUMAS/code/ directory


# please change this directory to the directory specified by jobID
gwas_path <- paste0(workdir, 'sumdata/')
output_path <- paste0(workdir, 'output/')
input_path <- paste0(workdir, 'input_for_eval/')
PennPRS_finalresults_path <- paste0(workdir, 'PennPRS_results/')
dir.create(gwas_path, showWarnings = F)
dir.create(input_path, showWarnings = F)
dir.create(output_path, showWarnings = F)
dir.create(PennPRS_finalresults_path, showWarnings = F)
# Create a separate directory 'PRS_model_training/' to store input for training PRS models
prsdir0 = paste0(workdir, 'PRS_model_training/')
dir.create(prsdir0, showWarnings = F)
for (method in methods) {
  prsdir = paste0(prsdir0, method,'/'); dir.create(prsdir, showWarnings = F)
}
output_path_eval = paste0(workdir, 'output_for_eval/')
dir.create(output_path_eval, showWarnings = F)

if (method %in% 'PROSPER'){
  prsdir = paste0(prsdir0, method,'/')
  summdata = paste0(prsdir, 'summdata/')
  dir.create(summdata, showWarnings = F)
  rscriptsdir = paste0(prsdir, 'rscripts/')
  dir.create(rscriptsdir, showWarnings = F)
  logfiledir = paste0(rscriptsdir, 'logfile/')
  dir.create(logfiledir, showWarnings = F)
  path_out = paste0(prsdir, 'output/')
  dir.create(path_out, showWarnings = F)
  path_out_lassosum2 = paste0(prsdir, 'output/lassosum2')
  dir.create(path_out_lassosum2, showWarnings = F)
  for (race in races) dir.create(paste0(path_out_lassosum2, '/', race), showWarnings = F)
  path_out_PROSPER = paste0(prsdir, 'output/PROSPER/')
  dir.create(path_out_PROSPER, showWarnings = F)
  for (race in races) dir.create(paste0(path_out_PROSPER, '/', race), showWarnings = F)
}
if (method %in% 'MUSSEL'){
  prsdir = paste0(prsdir0, method,'/')
  summdata = paste0(prsdir, 'summdata/')
  dir.create(summdata, showWarnings = F)
  rscriptsdir = paste0(prsdir, 'rscripts/')
  dir.create(rscriptsdir, showWarnings = F)
  logfiledir = paste0(rscriptsdir, 'logfile/')
  dir.create(logfiledir, showWarnings = F)
  path_out = paste0(prsdir, 'output/')
  dir.create(path_out, showWarnings = F)
  path_out_LDpred2 = paste0(prsdir, 'output/LDpred2/')
  dir.create(path_out_LDpred2, showWarnings = F)
  for (race in races) dir.create(paste0(path_out_LDpred2, '/', race), showWarnings = F)
  path_out_MUSSEL = paste0(prsdir, 'output/MUSSEL/')
  dir.create(path_out_MUSSEL, showWarnings = F)
  for (race in races) dir.create(paste0(path_out_MUSSEL, '/', race), showWarnings = F)
}


# copy the input GWAS summary data, {Ancestry}_{Trait}.txt, to the /sumdata/folder
for (race in races){
  trait_name = paste0(race,'_',trait)
  system(paste0('cp -r ',input_GWAS_path, trait_name,'.txt ', workdir, 'sumdata/'))
}


######## QC for GWAS Summary Data:
cat(paste0("\n********************************************"))
cat(paste0("\n**** Step 0: QC for the input GWAS data ****"))
cat(paste0("\n********************************************\n"))

for (race in races){
  trait_name = paste0(race,'_',trait)
  cat(paste0('**** Start QC for GWAS data from ', race, '. ****\n'))
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
      write_delim(sumraw, file = paste0(workdir, 'sumdata/', trait_name, '.txt'), delim='\t')
      if (length(rm.indx) == 1) print(paste0('* 1 problematic SNP removed. QC step completed.'))
      if (length(rm.indx) > 1) print(paste0('* QC step completed. ', nrow(sumraw), ' SNPs remaining. ', length(rm.indx), ' problematic SNPs removed.'))
    }
  }
  if (length(rm.indx) == 0) print(paste0('* QC step completed. ', nrow(sumraw), ' SNPs remaining. No SNP was removed.'))
}

cat(paste0('**** QC for GWAS data completed. ****\n'))


# All inputs/outputs are stored under different folders in the directory, workdir
# --------------------------------------------------------------------
# --------------------- Step 1: PUMAS Subsampling ---------------------
# --------------------------------------------------------------------
for (race in races){
  trait_name = paste0(race,'_',trait)
  ld_file <- paste0(ld_path0[race],race,'_LD_hm3.RData')
  rs_file <- paste0(ld_path0[race],race,'_rs_hm3.RData')
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
  print(paste0('Subsampling completed for ', race))
}



# --------------------------------------------------------------------
# --------------------- Step 2: Train PRS models using PROSPER ---------------------
# --------------------------------------------------------------------
if ('PROSPER' %in% methods){
  method = 'PROSPER'
  path_data = NULL
  for (ite in 1:k){
    Ngwas = numeric()
    path_data[[ite]] = character()
    for (race in races){
      trait_name = paste0(race,'_',trait)
      pumasout = paste0(output_path, trait_name, '.gwas.ite', ite, '.txt')
      if (!file.exists(pumasout)) print(paste0('Subsampling failed to generate ', ite,'-th fold of the MCCV summary data. Rerun pumas.subsampling.customized.R.'))
      if (file.exists(pumasout)){
        sumraw0 = bigreadr::fread2(pumasout)
        sumraw0 = sumraw0[,c('SNP', 'CHR', 'A1', 'A2', 'BETA', 'SE', 'N')]
        colnames(sumraw0) = c('rsid', 'chr', 'a1', 'a0', 'beta', 'beta_se', 'n_eff') # A1/a1: REF
        path_data[[ite]][race] = paste0(summdata, 'gwas.', trait_name, '.ite', ite, '.txt')
        Ngwas[race] = median(sumraw0$n_eff)
        write_delim(sumraw0, path_data[[ite]][race], delim = '\t')
        print(paste0(race, ' ', trait, ': Generating input for iteration ', ite, ' GWAS data for ', method, ' (step 1) completed.'))
      }
    }
    path_data[[ite]] = paste0(path_data[[ite]], collapse = ',')
    
    # --------------------- Step 1: Run lassosum2 ---------------------
    # cat(paste0('\n** Step 2: Run lassosum2 **\n'))
    FILE_sst = path_data[[ite]]
    lassosum2code = paste(paste0('Rscript ', PROSPER_path, '/scripts/lassosum2.R '),#paste0('Rscript ./code/PUMAS.subsampling.R '),
                          paste0('--PATH_package ', PROSPER_path),
                          paste0('--PATH_LD ', PennPRS_path, 'LD/'),
                          paste0('--PATH_out ', path_out_lassosum2),
                          paste0('--Ll ', nlambda),
                          paste0('--Ld ', ndelta),
                          paste0('--lambda.min.ratio ', lambda.min.ratio),
                          paste0('--PATH_plink ',path_plink),
                          paste0('--FILE_sst ', FILE_sst),
                          paste0('--pop ', paste0(races, collapse = ',')),
                          paste0('--chrom 1-22 '),
                          paste0('--NCORES ', NCORES))
    system(lassosum2code)
  }
  
  
  # ---------------------------------------------------------------------------------------
  # ----------- Step 2.3: Train lassosum2 on the whole data ----------
  # ---------------------------------------------------------------------------------------
  path_data_full = character()
  for (race in races){
    trait_name = paste0(race,'_',trait)
    sumraw0 = bigreadr::fread2(paste0(output_path,trait_name,".gwas_matched.txt"))
    sumraw0 = sumraw0[,c('SNP', 'CHR', 'A1', 'A2', 'BETA', 'SE', 'N')]
    colnames(sumraw0) = c('rsid', 'chr', 'a1', 'a0', 'beta', 'beta_se', 'n_eff') # A1/a1: REF
    path_data_full[race] = paste0(summdata, 'gwas_PROSPER_step1.', trait_name, '.full.txt')
    write_delim(sumraw0, path_data_full[race], delim = '\t')
    print(paste0(race, ' ', trait, ': Generating input GWAS data for ', method, ' completed.'))
  }
  path_data_full = paste0(path_data_full, collapse = ',')
  
  # --------------------- Run lassosum2 ---------------------
  if ( opt$verbose >= 1 ){
    print(paste0('******************************************************************'))
    print(paste0('******* Start training lassosum2 on the original GWAS data *******'))
    print(paste0('******************************************************************'))
  }
  lassosum2code = paste(paste0('Rscript ', PROSPER_path, '/scripts/lassosum2.R '),
                        paste0('--PATH_package ', PROSPER_path), 
                        paste0('--PATH_LD ', PennPRS_path, 'LD/'), 
                        paste0('--PATH_out ', path_out_lassosum2),
                        paste0('--Ll ', nlambda),
                        paste0('--Ld ', ndelta),
                        paste0('--lambda.min.ratio ', lambda.min.ratio),
                        paste0('--PATH_plink ',path_plink),
                        paste0('--FILE_sst ', path_data_full),
                        paste0('--pop ', paste0(races, collapse = ',')),
                        paste0('--chrom 1-22 '),
                        paste0('--NCORES ', NCORES))
  system(lassosum2code)
  
  
  # ------ Heritability estimated based the largest ancestry-specific GWAS:
  race.LD = names(Ngwas)[which.max(Ngwas)]
  map_ldref <- readRDS(paste0(ld_path[race.LD], 'map/map_',LDrefpanel,'_ldref.rds'))
  sumstats = bigreadr::fread2(paste0(output_path,race.LD,'_',trait,".gwas_matched.txt"))[,c('CHR','SNP','A1','A2','BETA','SE','P','N', 'MAF')]
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
  
  td = paste0(prsdir0, 'temporary_LD')
  if (!dir.exists(td)) dir.create(td)
  setwd(td)
  tmp <- tempfile(tmpdir = td)

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
      corr0 <- readRDS(paste0(path_precalLD[race.LD], '/LD_ref_chr', chr, '.rds'))[ind.chr3, ind.chr3]
      if ((chr == 1) | (is.null(ld))) {
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
  H2 <- abs(ldsc[["h2"]])
  cat(paste0('Heritability estimate based on LD score regression using the largest ancestry-specific GWAS: ', signif(H2, 3)))
  if (H2 == 0) H2 = 1e-3
  save(H2, file = paste0(input_path, trait,'.',method, '_H2.','.txt'))
  rm(ld, corr, sumstats, info_snp, sumraw0)
  gc()
  
  # --------------------------------------------------------------------
  # -------------- Step 3: PUMAS Evaluation on lassosum2 ---------------
  # --------------------------------------------------------------------
  # --------------------- Step 3.1: Input preparation for pumas.evaluation.R ---------------------
  # Step 2.3: Reformat the trained PRS model (SNP weight file) according to the format of the subsampled summary statistics 
  # from Step 1 (.gwas.omnibus.itei.txt) and use it as the input for pumas.evaluation.R
  for (race in races){
    trait_name = paste0(race,'_',trait)
    for (ite in 1:k){
      output_lassosum2 = paste0(path_out_lassosum2, '/', race, '/', trait_name, '.ite', ite, "_score_file.txt")
      if(file.exists(output_lassosum2)){
        score = bigreadr::fread2(output_lassosum2)  # "rsid"    "a1"      "a0"
        n.tuning = ncol(score) - 3
        colnames(score) = c('SNP', 'A1', 'A2', paste0('BETA',1:n.tuning))
      }
      
      # Match alleles with GWAS summary data:
      stateval = bigreadr::fread2(paste0(output_path, trait_name, '.gwas.ite', ite, '.txt'))
      stateval = stateval[, c('CHR','SNP','A1','A2')]
      colnames(stateval) = c('CHR','SNP', 'A1.ref','A2.ref')
      stateval = left_join(stateval,score,by = 'SNP')
      
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
      write_delim(scores, paste0(input_path, trait_name,'.',method, '_step1.ite',ite,'.txt'), delim='\t')
      rm(stateval)
    }
  }
  
  
  # ---------------------------------------------------------------------------------------
  # ----------- Evaluation + generate the best PRS(s) from lassosum2 on the whole data ----------
  # ---------------------------------------------------------------------------------------
  for (race in races){
    trait_name = paste0(race,'_',trait)
    prsdir = paste0(prsdir0, method,'/')
    xty_path = stats_path = output_path
    eval_ld_ref = paste0(eval_ld_ref_path[race], LDrefpanel,'_hm3_',race,'_ref') # or hm3+mega
    pumascode = paste(paste0('Rscript ', PUMAS_path, 'PUMAS.evaluation.customized.R'),
                      paste0('--k ',k), 
                      paste0('--ref_path ', eval_ld_ref),
                      paste0('--trait_name ',trait_name),
                      paste0('--prs_method ',method,'_step1'),
                      paste0('--xty_path ', xty_path),
                      paste0('--stats_path ', stats_path),
                      paste0('--weight_path ', input_path),
                      paste0('--output_path ', output_path_eval))
    system(pumascode)
    
    
    ##### Model tuning:
    r2 = bigreadr::fread2(paste0(output_path_eval,trait_name, '.', method, '_step1.txt'))
    r2.avg = colMeans(r2)
    # params.tuned = which.max(r2.avg)
    r2.order = list()
    n.candidates = min(15, ncol(r2))
    for (kk in 1:k) r2.order[[kk]] = order(as.numeric(r2[kk,]),decreasing = T)[1:n.candidates]
    # Select the top parameter settings:
    # params.tuned = as.numeric(substr(names(sort(r2.avg, decreasing = T)[1:5]), 5, 10))
    params.tuned = Reduce(intersect, r2.order)
    nonzero.indx = which(r2.avg > 0)
    params.tuned = unique(params.tuned[(params.tuned <= length(r2.avg)) & (params.tuned %in% nonzero.indx)])
    
    if (length(params.tuned) == 0){
      cat(paste0('Warning: all tuning parameter settings for ', race, ' in PROSPER Step 1 gives non-positive R2 values, which may lead to a final PRS model that has insufficient predictive power.'))
      params.tuned = Reduce(intersect, r2.order)
      params.tuned = unique(params.tuned[(params.tuned <= length(r2.avg))])
    }
    # Load in the PRS model trained based on the full GWAS dataset:
    output_lassosum2 = paste0(path_out_lassosum2, '/', race, '/', trait_name, '.full_score_file.txt')
    beta_lassosum2 = bigreadr::fread2(output_lassosum2)  
    beta_lassosum2[is.na(beta_lassosum2)] = 0
    
    # Columns that have nonzero entries: index in params.tuned and train.indx.lassosum2
    # nonzero.indx = which(sapply(params.tuned, function(x){sum(beta_lassosum2[,x]!=0)>0}))
    nonzero.indx = which(sapply(params.tuned, function(x){sum(abs(beta_lassosum2[,paste0('score',x)])>1e-7)>0}))
    if (length(nonzero.indx) > 0){
      candidates = params.tuned[nonzero.indx]
      stop = 0; ii = 0
      while ((stop == 0) & (ii < length(candidates))){
        ii = ii + 1
        if (((candidates[ii]-1) %in% candidates) | ((candidates[ii]+1) %in% candidates)){
          stop = 1
        }
      }
      optimal.indx = candidates[ii]
      indx.temp = ii
    }
    if (length(nonzero.indx) == 0){
      stop(paste0('Trained PRS model for ', race, 'based on ', method, 'has zero effect estimate for all SNPs. Please try other methods.'))
    }
    
    params.tuned = params.tuned[indx.temp]
    # Save optimal parameters for running PROSPER on each CV:
    for (ite in 1:k){
      params = bigreadr::fread2(paste0(path_out_lassosum2, '/', race, '/', trait_name, '.ite', ite, '_score_param.txt'))
      optimal.pars = params[params.tuned,]
      write_delim(optimal.pars, file = paste0(path_out_lassosum2, '/', race, '/', trait_name, '_optimal_param_ite', ite, '.txt'), delim = '\t')
    }
    # Save optimal parameters for running PROSPER on the full GWAS data:
    params = bigreadr::fread2(paste0(path_out_lassosum2, '/', race, '/', trait_name, '.full_score_param.txt'))
    optimal.pars = params[params.tuned,]
    write_delim(optimal.pars, file = paste0(path_out_lassosum2, '/', race, '/', trait_name, '_optimal_param_full.txt'), delim = '\t')
    
    print(paste0('Hyperparameter tuning for ', method,' completed on', race, '.'))
    print(paste0('R2 on tuning dataset for ', race, ': ', signif(r2.avg[params.tuned], 3))) #r2.avg[which.max(r2.avg)]))
    # Create an info/ folder to save files that contain tuned parameter information:
    tuned.parameters.file = paste0(workdir, 'PRS_model_training/',method,'/',method,'_step1_tuned_parameters_',trait_name,'.RData')
    save(params.tuned, file = tuned.parameters.file)
  }
}


if ('PRS-CSx' %in% methods){
  method = 'PRS-CSx'
  prsdir = paste0(prsdir0, method,'/')
  # --------------------- Step 2.1: Input preparation for PRS-CSx ---------------------
  # First, create a "pseudo" .bim file as one of the inputs of PRS-CSx: 
  # since we do not have individual-level validation dataset,  we need to create this pseudo .bim file:
  VALIDATION_BIM_PREFIX = paste0(input_path,'pseudo_validation.', paste(races,collapse='_'), '_',trait) # This is supposed to be a real .bim file for validation individuals, but we create a pseudo validation file instead 
  snp.list = list()
  n.vec = rep(0,K); names(n.vec) = races # the list of GWAS training sample size
  for(i in 1:K){
    race = races[i]
    trait_name = paste0(race,'_',trait)
    valbim = bigreadr::fread2(paste0(output_path,trait_name,'.gwas.ite1.txt')) # the set of SNPs is the same across the k MCCV files.
    n.vec[race] = round(median(valbim$N))
    snp.list[[i]] = valbim %>% mutate(V3 = 0) %>% mutate(BP = 1) 
    snp.list[[i]] = snp.list[[i]][,c('CHR','SNP','V3','BP','A1','A2')] # filter(CHR==chr) %>% 
  }
  snp.file = rbindlist(snp.list) %>% 
    distinct(SNP,.keep_all=TRUE)
  write_delim(snp.file, paste0(VALIDATION_BIM_PREFIX, '.bim'), delim = '\t', col_names = F)
  
  # Then, reformat summary data to use as the input data for PRS-CSx:
  for (race in races){
    trait_name = paste0(race,'_',trait)
    for (ite in 1:k){
      pumasout = paste0(output_path, trait_name, '.gwas.ite', ite, '.txt')
      if (!file.exists(pumasout)) print(paste0('Subsampling failed to generate ', ite,'-th fold of the MCCV summary data. Rerun pumas.subsampling.customized.R.'))
      if (file.exists(pumasout)){
        sumraw0 = bigreadr::fread2(pumasout)
        sumraw0 = sumraw0[,c('CHR', 'SNP', 'A1', 'A2', 'BETA', 'SE')]
        colnames(sumraw0) = c('CHR', 'SNP', 'A1', 'A2', 'BETA', 'SE') # A1: REF
        prscs.sumdat.file = paste0(prsdir,trait_name,'_reformated_gwas.ite', ite, '.txt')
        write_delim(sumraw0[, c('SNP', 'A1', 'A2', 'BETA', 'SE')], prscs.sumdat.file, delim = '\t')
        print(paste0(race, ' ', trait, ': Generating input for iteration ', ite, ' GWAS data for ', method, ' completed.'))
      }
    }
  }
  
  # --------------------- Step 2.2: Run PRS-CSx ---------------------
  PATH_TO_REFERENCE = paste0(PRScs_path,'ref/') # 1000 Genomes reference data
  SEED = 2024
  trait_names = paste0(races,'_',trait)
  chrs = paste0(1:22, collapse = ',')
  for (ite in 1:k){
    training_summary_data_filenames = paste(paste0(prsdir,trait_names,'_reformated_gwas.ite', ite, '.txt'), collapse=',')
    n_gwas = paste(n.vec, collapse=',')
    pop = paste(races, collapse=',')
    out_dir = paste0(prsdir, trait,'.ite',ite) 
    # here different from in PRS-CS, out_dir is a directory, not just a prefix
    # we don't need to specify race in the directory name: PRS-CSx itself will generate files by race
    if (!dir.exists(out_dir)) dir.create(out_dir)
    for (v in 1:length(phi.vals)){
      system(paste0("python ", PRScsx_path, "PRScsx.py",
                    " --ref_dir=", PATH_TO_REFERENCE,
                    " --bim_prefix=", VALIDATION_BIM_PREFIX,
                    " --sst_file=", training_summary_data_filenames,
                    " --n_gwas=", n_gwas,
                    " --pop=", pop,
                    " --chrom=", chrs,
                    " --phi=", phi.vals[v],
                    " --out_dir=", out_dir,
                    " --out_name=", method))
    }
    print(paste0('Complete training ', method, ' for ite ', ite))
  }
  
  
  # ------------------- Parameter Tuning
  for (ite in 1:k){
    for (race in races){
      trait_name = paste0(race,'_',trait)
      out_dir = paste0(prsdir, trait,'.ite',ite, '/', method, '_', race)
      SCORE = NULL
      for (phi in phi.vals){
        score = NULL
        for(chr in c(1:22)){
          temfile = paste0(out_dir, '_pst_eff_a1_b0.5_phi', scientific(phi,digits=2), '_chr',chr,'.txt')
          if(file.exists(temfile)){
            scoretemp = bigreadr::fread2(temfile)[,c(1,2,4,5,6)]; colnames(scoretemp) = c('CHR', 'SNP', 'A1', 'A2', paste0('BETA', which(phi.vals == phi)))
            score = rbind(score, scoretemp)
            rm(scoretemp)
            # print(paste0('Chr ', chr,' Completed'))
          }
          # if(!file.exists(temfile)) print(paste0('Need to rerun ', method, ' on CHR ',chr))
        }
        if (is.null(SCORE)) SCORE = score
        if (!is.null(SCORE)) {
          score = score[,c('SNP','A1',paste0('BETA', which(phi.vals == phi)))]
          flipped = which(score$A1 != SCORE$A1) 
          if (length(flipped)>0) score[,3] = - score[,3]
          SCORE[,paste0('BETA', which(phi.vals == phi))] = score[,paste0('BETA', which(phi.vals == phi))]
        }
      }
      score = SCORE; rm(SCORE)
      
      # Match alleles with GWAS summary data:
      stateval = bigreadr::fread2(paste0(output_path, trait_name, '.gwas.ite', ite, '.txt'))
      stateval = stateval[, c('CHR','SNP','A1','A2')]
      colnames(stateval) = c('CHR.ref','SNP', 'A1.ref','A2.ref')
      stateval = left_join(stateval,score,by="SNP")
      
      na.ind = which(is.na(stateval$A1))
      if (length(na.ind) > 0){
        stateval[na.ind, paste0('BETA',1:length(phi.vals))] = 0
        stateval[na.ind,'A1'] = stateval[na.ind,'A1.ref']
        stateval[na.ind,'A2'] = stateval[na.ind,'A2.ref']
      }
      flipped = which(stateval$A1.ref != stateval$A1) # April 14: 0 (already corrected)
      print(paste0(length(flipped), ' flipped SNPs.'))
      if (length(flipped) > 0){
        stateval[flipped,'A1'] = stateval[flipped,'A1.ref']
        stateval[flipped,'A2'] = stateval[flipped,'A2.ref']
        for (t in 1:length(phi.vals)){
          stateval[flipped,paste0('BETA',t)] = - stateval[flipped,paste0('BETA',t)]
        }
      }
      scores = stateval[,c('CHR.ref','SNP','A1','A2',paste0('BETA',1:length(phi.vals)))] # other files: SNP	CHR	A1	BETA1	BETA2	A2
      colnames(scores) = c('CHR','SNP','A1','A2',paste0('BETA',1:length(phi.vals)))
      # write_delim(scores, paste0(input_path, trait_name,'.',method, '.ite',ite,'.txt'), delim='\t')
      write.table(scores, paste0(input_path, trait_name,'.',method, '.ite', ite, '.txt'), row.names = F,col.names = T, quote = FALSE, sep = "\t" )
      rm(stateval, scores)
    }
  }
  
  for (race in races){
    trait_name = paste0(race,'_',trait)
    prsdir = paste0(prsdir0, method,'/')
    xty_path = stats_path = output_path
    eval_ld_ref = paste0(eval_ld_ref_path[race], LDrefpanel,'_hm3_',race,'_ref') # or hm3+mega
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
    for (kk in 1:k) r2.order[[kk]] = order(as.numeric(r2[kk,]),decreasing = T)[1:length(phi.vals)]
    # Select the top parameter settings:
    params.tuned = Reduce(intersect, r2.order)
    nonzero.indx = which(r2.avg > 0)
    params.tuned = unique(params.tuned[(params.tuned <= length(r2.avg)) & (params.tuned %in% nonzero.indx)])
    tuned.parameters.file = paste0(workdir, 'PRS_model_training/',method,'/tuned_parameters_',trait_name,'.RData')
    save(params.tuned, file = tuned.parameters.file)
    if (length(params.tuned) == 0){
      r2 = bigreadr::fread2(paste0(output_path_eval,trait_name, '.', method, '.txt'))
      r2.avg = colMeans(r2)
      nonzero.indx = which(r2.avg != 0)
      r2.order = list()
      for (kk in 1:k) r2.order[[kk]] = order(as.numeric(r2[kk, nonzero.indx]),decreasing = T)[1:length(phi.vals)]
      # Select the top parameter settings:
      params.tuned = Reduce(intersect, r2.order)
      params.tuned = nonzero.indx[params.tuned[!is.na(params.tuned)]]
      params.tuned = unique(params.tuned[(params.tuned <= length(r2.avg)) & (params.tuned %in% nonzero.indx)])
      save(params.tuned, file = tuned.parameters.file)
    }
  }
}



if ('MUSSEL' %in% methods){
  method = 'MUSSEL'
  prsdir = paste0(prsdir0, method,'/')
  trait_names = paste0(races,'_',trait)
  # --------------------- Step 2.1: Input preparation for MUSSEL ---------------------
  for (race in races){
    trait_name = paste0(race,'_',trait)
    for (ite in 1:k){
      pumasout = paste0(output_path, trait_name, '.gwas.ite', ite, '.txt')
      if (!file.exists(pumasout)) print(paste0('Subsampling failed to generate ', ite,'-th fold of the MCCV summary data. Rerun pumas.subsampling.customized.R.'))
      if (file.exists(pumasout)){
        sumraw0 = bigreadr::fread2(pumasout)
        sumraw0 = sumraw0[,c('SNP', 'CHR', 'A1', 'A2', 'BETA', 'SE', 'N')]
        colnames(sumraw0) = c('rsid',	'chr',	'a1',	'a0',	'beta',	'beta_se', 'n_eff') # a0 (in LDpred2): A1: REF, effective (reference) allele (counted allele in regression), the allele which beta corresponds to.
        suppressWarnings(dir.create(paste0(prsdir, 'ite',ite)))
        suppressWarnings(dir.create(paste0(prsdir, 'ite',ite, '/summdata')))
        mussel.sumdat.file = paste0(prsdir,'ite',ite,'/summdata/', race, '.txt')
        write_delim(sumraw0, mussel.sumdat.file, delim = '\t')
        print(paste0(method, ' - reformatting input summary data for iteration ', ite, ' for ', race, ' completed.'))
      }
    }
  }
  
  # --------------------- Step 2.2: Run MUSSEL on each fold of subsampled training data ---------------------
  SEED = 2024
  for (ite in 1:k){
    path_data = paste0(prsdir,'ite',ite,'/')
    pop = paste(races, collapse=',')
    out_dir = paste0(prsdir, trait,'.ite',ite) 
    if (!dir.exists(out_dir)) dir.create(out_dir)
    if (!dir.exists(out_dir)) for (race in races) dir.create(paste0(out_dir, '/', race))
    system(paste0("Rscript ", MUSSEL_path, "R/LDpred2.R",
                  " --PATH_package=", MUSSEL_path,
                  " --PATH_ref=", paste0(PennPRS_path, '/LD/'),
                  " --PATH_out=", out_dir,
                  " --pop=", pop,
                  paste0(' --p ', opt$p), paste0(' --H2 ', opt$H2), paste0(' --sparse ', opt$sparse),
                  " --FILE_sst=", paste0(paste0(path_data, "summdata/", races, ".txt"), collapse = ','),
                  " --bfile_tuning=", paste0(sapply(1:K, function(x) {paste0(ld_path[races[x]], "1KGref_plinkfile/1kg_hm3_", races[x], "_ref")}), collapse = ','),
                  " --NCORES ", NCORES))
    print(paste0('Complete training ', method, ' for ite ', ite))
  }
  
  # --------------------- Step 2.2: Run MUSSEL on the full GWAS summary data ---------------------
  for (race in races){
    trait_name = paste0(race,'_',trait)
    sumraw0 = bigreadr::fread2(paste0(output_path,trait_name,".gwas_matched.txt"))
    sumraw0 = sumraw0[,c('SNP', 'CHR', 'A1', 'A2', 'BETA', 'SE', 'N')]
    colnames(sumraw0) = c('rsid',	'chr',	'a1',	'a0',	'beta',	'beta_se', 'n_eff') # a1: A1: REF, effective (reference) allele (counted allele in regression), the allele which beta corresponds to.
    suppressWarnings(dir.create(paste0(prsdir, 'full')))
    suppressWarnings(dir.create(paste0(prsdir, 'full/summdata')))
    mussel.sumdat.file = paste0(prsdir,'full/summdata/', race, '.txt')
    write_delim(sumraw0, mussel.sumdat.file, delim = '\t')
    print(paste0(method, ' - reformatting full input summary data for ', race, ' completed.'))
  }
  
  SEED = 2024
  path_data = paste0(prsdir,'full/')
  pop = paste(races, collapse=',')
  out_dir = paste0(prsdir, trait,'.full') 
  if (!dir.exists(out_dir)) dir.create(out_dir)
  if (!dir.exists(out_dir)) for (race in races) dir.create(paste0(out_dir, '/', race))
  system(paste0("Rscript ", MUSSEL_path, "R/LDpred2.R",
                " --PATH_package=", MUSSEL_path,
                " --PATH_ref=", paste0(PennPRS_path, '/LD/'),
                " --PATH_out=", out_dir,
                " --pop=", pop,
                paste0(' --p ', opt$p), paste0(' --H2 ', opt$H2), paste0(' --sparse ', opt$sparse),
                " --FILE_sst=", paste0(paste0(path_data, "summdata/", races, ".txt"), collapse = ','),
                " --bfile_tuning=", paste0(sapply(1:K, function(x) {paste0(ld_path[races[x]], "1KGref_plinkfile/1kg_hm3_", races[x], "_ref")}), collapse = ','),
                " --NCORES ", NCORES))
  print(paste0('Complete training ', method, ' on the original GWAS summary data'))
  
  
  # ------------------------------------------------------------------
  # -------------- Step 3: PUMAS Evaluation on LDpred2 ---------------
  # ------------------------------------------------------------------
  # --------------------- Step 3.1: Input preparation for pumas.evaluation.R ---------------------
  for (race in races){
    trait_name = paste0(race,'_',trait)
    for (ite in 1:k){
      output_LDpred2 = paste0(prsdir, trait,'.ite',ite, '/', race, '/tmp/beta_files/beta_in_all_settings/ldpred2effects.txt')
      if(file.exists(output_LDpred2)){
        score = bigreadr::fread2(output_LDpred2)
        n.tuning = ncol(score) - 4
        score = score[, c('rsid', 'a0', 'a1', paste0('e',1:n.tuning))]
        colnames(score) = c('SNP', 'A1', 'A2', paste0('BETA',1:n.tuning)) # in LDpred2.R, a0 is the effect allele, so we have to flip a1 and a0 here.
      }
      
      # Match alleles with GWAS summary data:
      stateval = bigreadr::fread2(paste0(output_path, trait_name, '.gwas.ite', ite, '.txt'))
      stateval = stateval[, c('CHR','SNP','A1','A2')]
      colnames(stateval) = c('CHR','SNP', 'A1.ref','A2.ref')
      stateval = left_join(stateval,score,by = 'SNP')
      
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
      write_delim(scores, paste0(input_path, trait_name,'.',method, '_step1.ite',ite,'.txt'), delim='\t')
      rm(stateval)
    }
  }
  
  
  # ---------------------------------------------------------------------------------------
  # ----------- Evaluation + generate the best PRS(s) from LDpred2 on the whole data ----------
  # ---------------------------------------------------------------------------------------
  for (race in races){
    trait_name = paste0(race,'_',trait)
    xty_path = stats_path = output_path
    eval_ld_ref = paste0(eval_ld_ref_path[race], LDrefpanel,'_hm3_',race,'_ref') # or hm3+mega
    pumascode = paste(paste0('Rscript ', PUMAS_path, 'PUMAS.evaluation.customized.R'),
                      paste0('--k ',k), 
                      paste0('--ref_path ', eval_ld_ref),
                      paste0('--trait_name ',trait_name),
                      paste0('--prs_method ',method,'_step1'),
                      paste0('--xty_path ', xty_path),
                      paste0('--stats_path ', stats_path),
                      paste0('--weight_path ', input_path),
                      paste0('--output_path ', output_path_eval))
    system(pumascode)
    
    
    ##### Model tuning:
    r2 = bigreadr::fread2(paste0(output_path_eval,trait_name, '.', method, '_step1.txt'))
    r2.avg = colMeans(r2)
    # params.tuned = which.max(r2.avg)
    r2.order = list()
    n.candidates = min(15, ncol(r2))
    for (kk in 1:k) r2.order[[kk]] = order(as.numeric(r2[kk,]),decreasing = T)[1:n.candidates]
    # Select the top parameter settings:
    # params.tuned = as.numeric(substr(names(sort(r2.avg, decreasing = T)[1:5]), 5, 10))
    params.tuned = Reduce(intersect, r2.order)
    nonzero.indx = which(r2.avg > 0)
    params.tuned = unique(params.tuned[(params.tuned <= length(r2.avg)) & (params.tuned %in% nonzero.indx)])
    
    if (length(params.tuned) == 0){
      cat(paste0('Warning: all tuning parameter settings for ', race, ' in MUSSEL Step 1 gives non-positive R2 values, which may lead to a final PRS model that has insufficient predictive power.'))
      params.tuned = Reduce(intersect, r2.order)
      params.tuned = unique(params.tuned[(params.tuned <= length(r2.avg))])
    }
    # Load in the PRS model trained based on the full GWAS dataset:
    output_LDpred2 = paste0(prsdir, trait,'.full/', race, '/tmp/beta_files/beta_in_all_settings/ldpred2effects.txt')
    beta_LDpred2 = bigreadr::fread2(output_LDpred2)  
    beta_LDpred2[is.na(beta_LDpred2)] = 0
    
    # Columns that have nonzero entries: index in params.tuned and train.indx.lassosum2
    # nonzero.indx = which(sapply(params.tuned, function(x){sum(beta_LDpred2[,x]!=0)>0}))
    nonzero.indx = which(sapply(params.tuned, function(x){sum(abs(beta_LDpred2[,paste0('e',x)])>1e-7)>0}))
    if (length(nonzero.indx) > 0){
      candidates = params.tuned[nonzero.indx]
      stop = 0; ii = 0
      while ((stop == 0) & (ii < length(candidates))){
        ii = ii + 1
        if (((candidates[ii]-1) %in% candidates) | ((candidates[ii]+1) %in% candidates)){
          stop = 1
        }
      }
      optimal.indx = candidates[ii]
      indx.temp = ii
    }
    if (length(nonzero.indx) == 0){
      stop(paste0('Trained PRS model for ', race, 'based on ', method, 'has zero effect estimate for all SNPs. Please try other methods.'))
    }
    
    params.tuned = params.tuned[indx.temp]
    # Save optimal parameters for running PROSPER on each CV:
    for (ite in 1:k){
      params = bigreadr::fread2(paste0(prsdir, trait,'.ite',ite, '/', race, '/tmp/beta_files/beta_in_all_settings/params.txt'))
      optimal.pars = params[params.tuned,]; colnames(optimal.pars) = c('p0', 'h20', 'sparse0')
      write_delim(optimal.pars, file = paste0(path_out_LDpred2, '/', race, '/', trait_name, '_optimal_param_ite', ite, '.txt'), delim = '\t')
    }
    # Save optimal parameters for running PROSPER on the full GWAS data:
    params = bigreadr::fread2(paste0(prsdir, trait,'.full/', race, '/tmp/beta_files/beta_in_all_settings/params.txt'))
    optimal.pars = params[params.tuned,]; colnames(optimal.pars) = c('p0', 'h20', 'sparse0')
    write_delim(optimal.pars, file = paste0(path_out_LDpred2, '/', race, '/', trait_name, '_optimal_param_full.txt'), delim = '\t')
    
    print(paste0('Hyperparameter tuning for ', method,' completed on', race, '.'))
    print(paste0('R2 on tuning dataset for ', race, ': ', signif(r2.avg[params.tuned], 3))) #r2.avg[which.max(r2.avg)]))
    # Create an info/ folder to save files that contain tuned parameter information:
    tuned.parameters.file = paste0(workdir, 'PRS_model_training/',method,'/',method,'_step1_tuned_parameters_',trait_name,'.RData')
    save(params.tuned, file = tuned.parameters.file)
  }
}





