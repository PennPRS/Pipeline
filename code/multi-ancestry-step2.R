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
              help="Candidate values for tuning parameter H2 (heritability = H2 * h2_est from LDSC) [default: %default]"),
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
PROSPER_path = paste0(PennPRS_path, 'software/PROSPER')
path_plink = '/dcl01/chatterj/data/jin/software/plink2'
PRScs_path = paste0(PennPRS_path, 'software/PRScs/')
PRScsx_path = paste0(PennPRS_path, 'software/PRScsx/')
MUSSEL_path = paste0(PennPRS_path, 'software/MUSSEL/')
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

source(paste0(PUMAS_path, 'PennPRS_functions.R')) # please save the PennPRS_functions.R file to the /PUMAS/code/ directory

trait_names = paste0(races,'_',trait)
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

if (method %in% 'PROSPER'){
  prsdir = paste0(prsdir0, method,'/')
  summdata = paste0(prsdir, 'summdata/')
  rscriptsdir = paste0(prsdir, 'rscripts/')
  logfiledir = paste0(rscriptsdir, 'logfile/')
  path_out = paste0(prsdir, 'output/')
  path_out_lassosum2 = paste0(prsdir, 'output/lassosum2')
  path_out_PROSPER = paste0(prsdir, 'output/PROSPER/')
}
if (method %in% 'MUSSEL'){
  prsdir = paste0(prsdir0, method,'/')
  summdata = paste0(prsdir, 'summdata/')
  rscriptsdir = paste0(prsdir, 'rscripts/')
  logfiledir = paste0(rscriptsdir, 'logfile/')
  path_out = paste0(prsdir, 'output/')
  path_out_LDpred2 = paste0(prsdir, 'output/LDpred2/')
  path_out_MUSSEL = paste0(prsdir, 'output/MUSSEL/')
}


# --------------------------------------------------------------------
# --------------------- Step 2: Train PRS models ---------------------
# --------------------------------------------------------------------
if ('PROSPER' %in% methods){
  Ngwas = numeric()
  method = 'PROSPER'
  prsdir = paste0(prsdir0, method,'/') 
  load(paste0(input_path, trait,'.',method, '_H2.','.txt'))
  
  path_data = NULL
  for (ite in 1:k){
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
  }
  
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
  
  # --------------------------------------------------------------------
  # -------------- Step 4: PROSPER ---------------
  # --------------------------------------------------------------------
  # ---------------------------------------------------------------------------------------
  # ----------- Step 4.1: Train PROSPER on the subsampled training datasets ----------
  # ---------------------------------------------------------------------------------------
  path_data = lassosum_param = NULL
  for (ite in 1:k){
    path_data[[ite]] = lassosum_param[[ite]] = character()
    for (race in races){
      trait_name = paste0(race,'_',trait)
      pumasout = paste0(output_path, trait_name, '.gwas.ite', ite, '.txt')
      path_data[[ite]][race] = paste0(summdata, 'gwas.', trait_name, '.ite', ite, '.txt')
      lassosum_param[[ite]][race] = paste0(path_out_lassosum2, '/', race, '/', trait_name, '_optimal_param_ite',ite,'.txt')
    }
    path_data[[ite]] = paste0(path_data[[ite]], collapse = ',')
    lassosum_param[[ite]] = paste0(lassosum_param[[ite]], collapse = ',')
    
    # --------------------- Step 2.2: Run PROSPER ---------------------
    # cat(paste0('\n** Step 2: Run lassosum2 **\n'))
    FILE_sst = path_data[[ite]]
    paste0(races, collapse = ',')
    prospercode = paste(paste0('Rscript ', PROSPER_path, '/scripts/PROSPER.R '),
                        paste0('--PATH_package ', PROSPER_path), 
                        paste0('--PATH_LD ', PennPRS_path, 'LD/'), 
                        paste0('--PATH_out ', path_out_PROSPER),
                        paste0('--FILE_sst ', FILE_sst),
                        paste0('--Ll ', Ll),
                        paste0('--Lc ', Lc),
                        paste0('--pop ', paste0(races, collapse = ',')),
                        paste0('--lassosum_param ', lassosum_param),
                        paste0('--ite ', ite),
                        paste0('--chrom 1-22 '),
                        paste0('--NCORES ', NCORES))
    system(prospercode)
    cat(paste0('Completed training PROSPER for MCCV iteration ', ite, '.'))
  }
  
  # ---------------------------------------------------------------------------------------
  # ----------- Step 4.2: Train PROSPER on the whole data ----------
  # ---------------------------------------------------------------------------------------
  lassosum_param_full = character()
  for (race in races){
    trait_name = paste0(race,'_',trait)
    xty_path = stats_path = output_path
    lassosum_param_full[race] = paste0(path_out_lassosum2, '/', race, '/', trait_name, '_optimal_param_full.txt')
  }
  lassosum_param_full = paste0(lassosum_param_full, collapse = ',')
  
  prospercode = paste(paste0('Rscript ', PROSPER_path, '/scripts/PROSPER.R '),
                      paste0('--PATH_package ', PROSPER_path), 
                      paste0('--PATH_LD ', PennPRS_path, 'LD/'), 
                      paste0('--PATH_out ', path_out_PROSPER),
                      paste0('--FILE_sst ', path_data_full),
                      paste0('--Ll ', Ll),
                      paste0('--Lc ', Lc),
                      paste0('--pop ', paste0(races, collapse = ',')),
                      paste0('--lassosum_param ', lassosum_param_full),
                      paste0('--ite full'),
                      paste0('--chrom 1-22 '),
                      paste0('--NCORES ', NCORES))
  system(prospercode)
  cat(paste0('Completed training PROSPER on the original GWAS data.'))
  
  # ------------------------------------------------------------------------------
  # --------------------- Step 5: PUMAS Evaluation on PROSPER + Generate best PRS(s) for PROSPER --------------------
  # ------------------------------------------------------------------------------
  # --------------------- Step 5.1: Input preparation for pumas.evaluation.R ---------------------
  # Step 2.3: Reformat the trained PRS model (SNP weight file) according to the format of the subsampled summary statistics 
  # from Step 1 (.gwas.omnibus.itei.txt) and use it as the input for pumas.evaluation.R
  for (race in races){
    trait_name = paste0(race,'_',trait)
    for (ite in 1:k){
      output_prosper = paste0(path_out_PROSPER, 'before_ensemble/score_file', ite, '.txt')
      if(file.exists(output_prosper)){
        score = bigreadr::fread2(output_prosper)  # "rsid"    "a1"      "a0"
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
      write_delim(scores, paste0(input_path, trait_name,'.',method, '.ite',ite,'.txt'), delim='\t')
      rm(stateval)
    }
  }
  rm(score)
  
  
  # Load in the PRS model trained based on the full GWAS dataset
  # One for all races:
  output_prosper = paste0(path_out_PROSPER, '/before_ensemble/score_filefull.txt')
  beta_prosper = bigreadr::fread2(output_prosper)  
  beta_prosper[is.na(beta_prosper)] = 0
  colnames(beta_prosper)[1:3] = c('SNP', 'A1', ' A2')
  rm.indx = which(sapply(1:(ncol(beta_prosper)-3),function(x){sum(beta_prosper[,3+x]^2)}) > 0.8)
  save(rm.indx, file = paste0(workdir, 'PRS_model_training/',method,'/',method,'_rm.indx_',trait,'.RData'))
  
  stateval0 = NULL
  for (race in races){
    trait_name = paste0(race,'_',trait)
    tem = bigreadr::fread2(paste0(output_path, trait_name, '.gwas.ite1.txt'))
    tem = tem[, c('CHR','SNP')]
    stateval0 = rbind(stateval0, tem)
  }
  stateval0 = stateval0[!duplicated(stateval0),]
  beta_prosper = left_join(beta_prosper, stateval0, by = 'SNP')
  
  for (race in races){
    trait_name = paste0(race,'_',trait)
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
    
    ##### Model tuning:
    r2 = bigreadr::fread2(paste0(output_path_eval,trait_name, '.', method, '.txt'))
    r2.avg = colMeans(r2)
    # params.tuned = which.max(r2.avg)
    r2.order = list()
    n.candidates = ncol(r2) # min(15, ncol(r2))
    for (kk in 1:k) r2.order[[kk]] = order(as.numeric(r2[kk,]),decreasing = T)[1:n.candidates]
    # Select the top parameter settings:
    # params.tuned = as.numeric(substr(names(sort(r2.avg, decreasing = T)[1:5]), 5, 10))
    params.tuned = Reduce(intersect, r2.order)
    nonzero.indx = which(r2.avg > 0)
    params.tuned = unique(params.tuned[(params.tuned <= length(r2.avg)) & (params.tuned %in% nonzero.indx)])
    
    params.tuned = setdiff(params.tuned, rm.indx)
    a = cbind(sapply(1:(ncol(beta_prosper)-4),function(x){sum(beta_prosper[,3+x]^2)}), 
              sapply(1:(ncol(beta_prosper)-4),function(x){sum(abs(beta_prosper[,3+x]) > 0.1)}),
              sapply(1:(ncol(beta_prosper)-4),function(x){sum(abs(beta_prosper[,3+x]) != 0)}),
              t(r2))
    keep.indx = as.numeric(substr(names(a[((a[,1]< 1.5) & (a[,4]< H2)),4]), 5,7))
    if (sum(params.tuned %in% keep.indx) > 0) params.tuned = params.tuned[params.tuned %in% keep.indx]
    # params.tuned = params.tuned[1:(min(5,length(params.tuned)))]
    
    nonzero.indx = which(sapply(params.tuned, function(x){sum(abs(beta_prosper[,paste0('score',x)])>1e-7)>0}))
    if (length(nonzero.indx) > 0){
      params.tuned = params.tuned[nonzero.indx]
      params.tuned = params.tuned[1:(min(5,length(params.tuned)))]
      tuned.parameters.file = paste0(workdir, 'PRS_model_training/',method,'/',method,'_tuned_parameters_',trait_name,'.RData')
      save(params.tuned, file = tuned.parameters.file)
      # Save the R2s on testing data for the selected tuning parameter settings:
      pumas.cor2 = r2[,params.tuned]
      write.table(pumas.cor2,paste0(PennPRS_finalresults_path,trait_name,".",method,".testing.txt"),col.names = T,row.names=F,quote=F,sep="\t")
      
      
      # Save the best PROSPER model:
      params.tuned.full = params.tuned
      params.tuned = params.tuned[1]
      params = bigreadr::fread2(paste0(path_out_PROSPER, 'before_ensemble/score_paramfull.txt'))
      optimal.pars = params[params.tuned,]
      # write_delim(optimal.pars, file = paste0(PennPRS_finalresults_path, trait_name, '_', method, '_single_PRS_optimal_param.txt'), delim = '\t')
      
      beta_prosper_single = data.frame(beta_prosper[,c('CHR', 'SNP', 'A1', ' A2', paste0('score', params.tuned))])
      colnames(beta_prosper_single) = c('CHR', 'SNP','A1', 'A2', 'BETA')
      
      nonzero = (sum(beta_prosper_single$BETA!=0)>0)
      if (nonzero){
        # Remove SNPs with zero effects:
        beta_prosper_single = beta_prosper_single[beta_prosper_single$BETA != 0,]
        prs_prosper_single_outputfile = paste0(PennPRS_finalresults_path, trait_name,'.',method,'.PRS_single_best.txt')
        write_delim(beta_prosper_single, file = prs_prosper_single_outputfile, delim='\t')
        print(paste0('Optimal parameter setting: '))
        optimal.pars[1, c(3,4,5,7)] = signif(optimal.pars[1, c(3,4,5,7)],3)
        print(optimal.pars)
        print(paste0('Complete training ', method, ' for ',race))
      }
      if (!nonzero){
        print(paste0('Trained PRS model based on ', method, 'has zero effect estimate for all SNPs. Please try other methods.'))
      }
    }
    if (length(nonzero.indx) == 0){
      print(paste0('Trained PRS model based on ', method, 'has zero effect estimate for all SNPs. Please try other methods.'))
    }
    
    # ------------------------------------------------------------------------------
    # --------------------- Step 6: Train Ensemble PROSPER PRS ---------------------
    # ------------------------------------------------------------------------------
    pumascode = paste(paste0('Rscript ', PUMAS_path, 'PUMAS.evaluation.customized_MA_single_method_ensemble_2subsamples.R'),
                      paste0('--k ',k), 
                      paste0('--ref_path ', eval_ld_ref),
                      paste0('--trait_name ', trait_name),
                      paste0('--prs_method ', method),
                      paste0('--optimal_params_path ', workdir, 'PRS_model_training/',method,'/'),
                      paste0('--xty_path ', xty_path),
                      paste0('--stats_path ', stats_path),
                      paste0('--weight_path ', input_path),
                      paste0('--output_path_eval ', output_path_eval),
                      paste0('--output_path ', PennPRS_finalresults_path))
    system(pumascode)
    # Generate ensemble PRS file:
    beta_prosper_ensemble = beta_prosper[,c('CHR', 'SNP', 'A1', ' A2', paste0('score',params.tuned.full))]
    colnames(beta_prosper_ensemble) = c('CHR', 'SNP','A1','A2', paste0('BETA',params.tuned.full))
    beta_prosper_ensemble0 = beta_prosper_ensemble
    
    # 1. Ensemble using PUMA-CUBS
    ensemble.weights = bigreadr::fread2(paste0(PennPRS_finalresults_path,trait_name,'.',method,".omnibus.weights.txt"))
    ensemble.weights = colMeans(ensemble.weights)
    ensemble.weights = ensemble.weights/sum(ensemble.weights)
    beta_prosper_ensemble$BETA = as.numeric(as.matrix(beta_prosper_ensemble[,paste0('BETA',params.tuned.full)]) %*% t(t(ensemble.weights)))
    beta_prosper_ensemble = beta_prosper_ensemble[, c('CHR','SNP','A1','A2','BETA')]
    prs_prosper_outputfile = paste0(PennPRS_finalresults_path, trait_name,'.',method,'.PRS.txt')
    write_delim(beta_prosper_ensemble, file = prs_prosper_outputfile, delim='\t')
  }
  
  if ( opt$verbose >= 1 ){
    print(paste0('*************************************'))
    print(paste0('******* Complete training ', method, ' *******'))
    print(paste0('*************************************'))
  }
}


if ('PRS-CSx' %in% methods){
  method = 'PRS-CSx'
  prsdir = paste0(prsdir0, method,'/')
  VALIDATION_BIM_PREFIX = paste0(input_path,'pseudo_validation.', paste(races,collapse='_'), '_',trait)
  n.vec = rep(0,K); names(n.vec) = races # the list of GWAS training sample size
  for(i in 1:K){
    race = races[i]; trait_name = paste0(race,'_',trait)
    valbim = bigreadr::fread2(paste0(output_path,trait_name,'.gwas.ite1.txt')) # the set of SNPs is the same across the k MCCV files.
    n.vec[race] = round(median(valbim$N))
  }
  # Reformat the original summary data to use as the input data for PRS-CSx:
  for (race in races){
    trait_name = paste0(race,'_',trait)
    for (ite in 1:k){
      sumraw0 = bigreadr::fread2(paste0(output_path,trait_name,".gwas_matched.txt"))[,c('SNP', 'A1', 'A2', 'BETA', 'SE')]
      prscs.sumdat.file = paste0(prsdir,trait_name,'_reformated_gwas.txt')
      write_delim(sumraw0, prscs.sumdat.file, delim = '\t')
      print(paste0(race, ' ', trait, ': Generating input GWAS data for ', method, ' completed.'))
    }
  }
  # --------------------- Run PRS-CSx ---------------------
  PATH_TO_REFERENCE = paste0(PRScs_path,'ref/')
  SEED = 2024
  chrs = paste0(1:22, collapse = ',')
  training_summary_data_filenames = paste(paste0(prsdir,trait_names,'_reformated_gwas.txt'), collapse=',')
  n_gwas = paste(n.vec, collapse=',')
  pop = paste(races, collapse=',')
  out_dir = paste0(prsdir, trait,'.full') 
  if (!dir.exists(out_dir)) dir.create(out_dir)
  phis = rep(0, length(races)); names(phis) = races
  for (race in races){
    trait_name = paste0(race,'_',trait)
    prsdir = paste0(prsdir0, method,'/')
    load(paste0(workdir, 'PRS_model_training/',method,'/tuned_parameters_',trait_name,'.RData'))
    phis[race] = phi.vals[params.tuned[1]]
  }
  for (phi in unique(phis)){
    system(paste0("python ", PRScsx_path, "PRScsx.py",
                  " --ref_dir=", PATH_TO_REFERENCE,
                  " --bim_prefix=", VALIDATION_BIM_PREFIX,
                  " --sst_file=", training_summary_data_filenames,
                  " --n_gwas=", n_gwas,
                  " --pop=", pop,
                  " --chrom=", chrs,
                  " --phi=", phi,
                  " --out_dir=", out_dir,
                  " --out_name=", method))
  }
  print(paste0('Complete training ', method, ' on the original GWAS summary data'))
  
  # --------------------- Reformat the trained PRS weight file: ---------------------
  for (race in races){
    trait_name = paste0(race,'_',trait)
    prsdir = paste0(prsdir0, method,'/')
    load(paste0(workdir, 'PRS_model_training/',method,'/tuned_parameters_',trait_name,'.RData'))
    score = NULL
    for(chr in c(1:22)){
      # if (is.na(phi)) temfile = paste0(out_dir, '_pst_eff_a1_b0.5_phiauto_chr',chr,'.txt')
      temfile = paste0(out_dir, '/', method, '_', race, '_pst_eff_a1_b0.5_phi', phi.vals[params.tuned[1]], '_chr',chr,'.txt')
      if(file.exists(temfile)){
        scoretemp = bigreadr::fread2(temfile)[, c(1, 2, 4, 5, 6)]
        colnames(scoretemp) = c('CHR', 'SNP', 'A1', 'A2', 'BETA')
        score = rbind(score, scoretemp)
        rm(scoretemp)
        # print(paste0('Chr ', chr,' Completed'))
      }
      # if(!file.exists(temfile)) print(paste0('Need to rerun ', method, ' on CHR ',chr))
    }
    beta_prscsx = score[,c('CHR','SNP','A1','A2','BETA')]
    beta_prscsx[is.na(beta_prscsx)] = 0
    
    prs_prscsx_outputfile = paste0(PennPRS_finalresults_path, trait_name,'.',method,'.PRS.txt')
    # write_delim(beta_prscsx, prs_prscsx_outputfile, delim='\t')
    write.table(beta_prscsx, prs_prscsx_outputfile, row.names = F, col.names = T, quote = FALSE, sep = "\t" )
  }
  if (opt$verbose >= 1 ){
    print(paste0('*************************************'))
    print(paste0('***** Complete training PRS-CSx *****'))
    print(paste0('*************************************'))
  }
}





if ('MUSSEL' %in% methods){
  method = 'MUSSEL'
  prsdir = paste0(prsdir0, method,'/')
  # -----------------------------------------------------------------------------
  # -------------- Run MUSSEL on the Pseudo Training summary Data ---------------
  # -----------------------------------------------------------------------------
  # ----------------------------------------------
  # ----------- Run MUSS by chromosome -----------
  # ----------------------------------------------
  SEED = 2024
  for (ite in 1:k){
    cat(paste0('Running ', method, ' for MCCV iteration ', ite, '...'))
    path_data = paste0(prsdir,'ite',ite,'/')
    pop = paste(races, collapse=',')
    out_dir = paste0(prsdir, trait,'.ite',ite) 
    if (!dir.exists(out_dir)) dir.create(out_dir)
    if (!dir.exists(out_dir)) for (race in races) dir.create(paste0(out_dir, '/', race))
    for (chr in 1:22){
      system(paste0("Rscript ", MUSSEL_path, "R/MUSS.R",
                    " --PATH_package=", MUSSEL_path,
                    " --PATH_PennPRS=", PennPRS_path, 
                    " --PATH_LDref=", paste0(PennPRS_path, '/LD/'),
                    " --PATH_out=", out_dir,
                    " --pop=", pop,
                    ' --LDpred2_params=', paste(paste0(path_out_LDpred2, '/', races, '/', races,'_',trait, '_optimal_param_ite', ite, '.txt'), collapse = ','),
                    ' --chrom ', chr,
                    ' --cors_additional ', opt$cors_additional,
                    ' --ps_additional ', opt$ps_additional,
                    " --FILE_sst=", paste0(paste0(path_data, "summdata/", races, ".txt"), collapse = ','),
                    " --bfile_tuning=", paste0(sapply(1:K, function(x) {paste0(ld_path[races[x]], "1KGref_plinkfile/1kg_hm3_", races[x], "_ref")}), collapse = ','),
                    " --NCORES ", NCORES))
      cat(paste0('CHR ', chr, ': completed.'))
    }
    cat(paste0(method, ' for MCCV iteration ', ite, ': completed.'))
  }
  
  # ---------------------------------------------------------------
  # ----------- Step 4.2: Train MUSSEL on the whole data ----------
  # ---------------------------------------------------------------
  cat(paste0('Running ', method, ' on the original GWAS summary data...'))
  path_data = paste0(prsdir,'full/')
  pop = paste(races, collapse=',')
  out_dir = paste0(prsdir, trait,'.full') 
  if (!dir.exists(out_dir)) dir.create(out_dir)
  if (!dir.exists(out_dir)) for (race in races) dir.create(paste0(out_dir, '/', race))
  for (chr in c(1:22)){
    system(paste0("Rscript ", MUSSEL_path, "R/MUSS.R",
                  " --PATH_package=", MUSSEL_path,
                  " --PATH_PennPRS=", PennPRS_path, 
                  " --PATH_LDref=", paste0(PennPRS_path, '/LD/'),
                  " --PATH_out=", out_dir,
                  " --pop=", pop,
                  ' --LDpred2_params=', paste(paste0(path_out_LDpred2, '/', races, '/', races,'_',trait, '_optimal_param_full.txt'), collapse = ','),
                  ' --chrom ', chr,
                  ' --cors_additional ', opt$cors_additional,
                  ' --ps_additional ', opt$ps_additional,
                  " --FILE_sst=", paste0(paste0(path_data, "summdata/", races, ".txt"), collapse = ','),
                  " --bfile_tuning=", paste0(sapply(1:K, function(x) {paste0(ld_path[races[x]], "1KGref_plinkfile/1kg_hm3_", races[x], "_ref")}), collapse = ','),
                  " --NCORES ", NCORES))
    cat(paste0('CHR ', chr, ': completed.'))
  }
  cat(paste0('Completed training MUSSEL on the original GWAS data.'))
  
  # ------------------------------------------------------------------------------
  # --------------------- Step 5: PUMAS Evaluation on PROSPER + Generate best PRS(s) for PROSPER --------------------
  # ------------------------------------------------------------------------------
  # --------------------- Step 5.1: Input preparation for pumas.evaluation.R ---------------------
  # Note: we combine all ancestry-specific models
  ref = bigreadr::fread2(paste0(PennPRS_path, 'software/PROSPER/ref_bim.txt'))[,c(1,2,5,6)]
  colnames(ref) = c('CHR.ref', 'SNP', 'A1.ref', 'A2.ref')
  rownames(ref) = ref$SNP; ref0 = ref; rm(ref)

  for (ite in 1:k){
    beta_mussel = list()
    for (race in races){
      trait_name = paste0(race,'_',trait)
      beta = NULL
      for (chr in 1:22){
        output_mussel = paste0(prsdir, trait,'.ite',ite,'/tmp/MUSS_beta_in_all_settings_bychrom/', race, '-chr', chr, '.txt')
        if(file.exists(output_mussel)){
          score = bigreadr::fread2(output_mussel)  # "SNP", "ALT"
          score[is.na(score)] = 0
          n.tuning = ncol(score) - 2
          colnames(score) = c('SNP', 'A2', paste0('BETA',1:n.tuning)) # a1 in MUSS.R (a0 is the effect allele), need to flip here.
          beta = rbind(beta, score); rm(score)
          # print(chr)
        }
      }
      beta_mussel[[race]] = beta; rm(beta)
      stateval0 = NULL
      tem = bigreadr::fread2(paste0(output_path, trait_name, '.gwas.ite1.txt'))
      tem = tem[, c('CHR','SNP')]
      stateval0 = rbind(stateval0, tem)
      stateval0 = stateval0[!duplicated(stateval0),]
      beta_mussel[[race]] = left_join(beta_mussel[[race]], stateval0, by = 'SNP')
    }
    beta = beta_mussel[[1]][, c('CHR', 'SNP', 'A2', paste0('BETA', 1:n.tuning))]; anc = 1
    colnames(beta) = c(paste0('CHR.', anc), 'SNP', paste0('A2.', anc), paste0('BETA', (n.tuning*(anc-1) + 1) : (n.tuning*anc)))
    rownames(beta) = beta$SNP
    for (anc in 2:length(beta_mussel)){
      tem = beta_mussel[[anc]][, c('CHR', 'SNP', 'A2', paste0('BETA', 1:n.tuning))]
      colnames(tem) = c(paste0('CHR.', anc), 'SNP', paste0('A2.', anc), paste0('BETA', (n.tuning*(anc-1) + 1) : (n.tuning*anc)))
      rownames(tem) = tem$SNP
      beta = full_join(beta, tem, by = 'SNP')
    }
    # Match alleles:
    ref = ref0[beta$SNP, c('CHR.ref', 'SNP', 'A1.ref', 'A2.ref')]
    for (anc in 1:length(beta_mussel)){
      na.ind = which(is.na(beta[, paste0('A2.', anc)]))
      if (length(na.ind) > 0){
        beta[na.ind, paste0('BETA', (n.tuning*(anc-1) + 1) : (n.tuning*anc))] = 0
        beta[na.ind, paste0('A1.', anc)] = ref[na.ind, 'A1.ref']
        beta[na.ind, paste0('A2.', anc)] = ref[na.ind, 'A2.ref']
      }
      flipped = which(ref$A2.ref != beta[,paste0('A2.', anc)])
      print(paste0(length(flipped), ' flipped SNPs for ', races[anc]))
      if (length(flipped) > 0){
        beta[flipped, paste0('A1.', anc)] = ref[flipped, 'A1.ref']
        beta[flipped, paste0('A2.', anc)] = ref[flipped, 'A2.ref']
        for (t in (n.tuning*(anc-1) + 1) : (n.tuning*anc)){
          beta[flipped,paste0('BETA',t)] = - beta[flipped,paste0('BETA',t)]
        }
      }
    }
    beta_mussel = cbind(ref[, c('CHR.ref', 'SNP', 'A1.ref', 'A2.ref')], beta[,paste0('BETA',1:(n.tuning*length(races)))]) # other files: SNP	CHR	A1	BETA1	BETA2	A2
    colnames(beta_mussel)[c(1,3,4)] = c('CHR', 'A1', 'A2')
    # ------------------- Match alleles with GWAS summary data for each ancestry separately ------------------- 
    for (race in races){
      trait_name = paste0(race,'_',trait)
      
      stateval = bigreadr::fread2(paste0(output_path, trait_name, '.gwas.ite', ite, '.txt'))
      stateval = stateval[, c('CHR','SNP','A1','A2')]
      colnames(stateval) = c('CHR','SNP', 'A1.ref','A2.ref')
      stateval = left_join(stateval, beta_mussel, by = 'SNP')
      
      na.ind = which(is.na(stateval$A1))
      if (length(na.ind) > 0){
        stateval[na.ind, paste0('BETA',1:(n.tuning*length(races)))] = 0
        stateval[na.ind,'A1'] = stateval[na.ind,'A1.ref']
        stateval[na.ind,'A2'] = stateval[na.ind,'A2.ref']
      }
      flipped = which(stateval$A1.ref != stateval$A1)
      print(paste0(length(flipped), ' flipped SNPs.'))
      if (length(flipped) > 0){
        stateval[flipped,'A1'] = stateval[flipped,'A1.ref']
        stateval[flipped,'A2'] = stateval[flipped,'A2.ref']
        for (t in 1:(n.tuning*length(races))){
          stateval[flipped,paste0('BETA',t)] = - stateval[flipped,paste0('BETA',t)]
        }
      }
      score = stateval[,c('CHR.x','SNP','A1','A2',paste0('BETA',1:(n.tuning*length(races))))] # other files: SNP	CHR	A1	BETA1	BETA2	A2
      colnames(score)[1] = 'CHR'
      write_delim(score, paste0(input_path, trait_name,'.',method, '.ite',ite,'.txt'), delim = '\t')
      rm(stateval,score)
    }
  }

  # Load in the PRS model trained based on the full GWAS dataset
  beta_mussel = list()
  for (race in races){
    trait_name = paste0(race,'_',trait)
    beta = NULL
    for (chr in 1:22){
      output_mussel = paste0(prsdir, trait,'.full/tmp/MUSS_beta_in_all_settings_bychrom/', race, '-chr', chr, '.txt')
      if(file.exists(output_mussel)){
        score = bigreadr::fread2(output_mussel) # "SNP", "ALT"
        score[is.na(score)] = 0
        n.tuning = ncol(score) - 2
        colnames(score) = c('SNP', 'A2', paste0('BETA',1:n.tuning)) # a1 in MUSS.R (a0 is the effect allele), need to flip here.
        beta = rbind(beta, score); rm(score)
        # print(chr)
      }
    }
    beta_mussel[[race]] = beta; rm(beta)
    stateval0 = NULL
    tem = bigreadr::fread2(paste0(output_path, trait_name, '.gwas.ite1.txt'))
    tem = tem[, c('CHR','SNP')]
    stateval0 = rbind(stateval0, tem)
    stateval0 = stateval0[!duplicated(stateval0),]
    beta_mussel[[race]] = left_join(beta_mussel[[race]], stateval0, by = 'SNP')
    rm.indx = which(sapply(1:(ncol(beta_mussel[[race]])-2),function(x){sum(beta_mussel[[race]][,2+x]^2)}) > 0.8)
    save(rm.indx, file = paste0(workdir, 'PRS_model_training/', method, '/', trait_name, '_', method, '_rm.indx.RData'))
  }

  beta = beta_mussel[[1]][, c('CHR', 'SNP', 'A2', paste0('BETA', 1:n.tuning))]; anc = 1
  colnames(beta) = c(paste0('CHR.', anc), 'SNP', paste0('A2.', anc), paste0('BETA', (n.tuning*(anc-1) + 1) : (n.tuning*anc)))
  rownames(beta) = beta$SNP
  for (anc in 2:length(beta_mussel)){
    tem = beta_mussel[[anc]][, c('CHR', 'SNP', 'A2', paste0('BETA', 1:n.tuning))]
    colnames(tem) = c(paste0('CHR.', anc), 'SNP', paste0('A2.', anc), paste0('BETA', (n.tuning*(anc-1) + 1) : (n.tuning*anc)))
    rownames(tem) = tem$SNP
    beta = full_join(beta, tem, by = 'SNP')
  }
  # Match alleles:
  ref = ref0[beta$SNP, c('CHR.ref', 'SNP', 'A1.ref', 'A2.ref')]
  for (anc in 1:length(beta_mussel)){
    na.ind = which(is.na(beta[, paste0('A2.', anc)]))
    if (length(na.ind) > 0){
      beta[na.ind, paste0('BETA', (n.tuning*(anc-1) + 1) : (n.tuning*anc))] = 0
      beta[na.ind, paste0('A1.', anc)] = ref[na.ind, 'A1.ref']
      beta[na.ind, paste0('A2.', anc)] = ref[na.ind, 'A2.ref']
    }
    flipped = which(ref$A2.ref != beta[,paste0('A2.', anc)])
    print(paste0(length(flipped), ' flipped SNPs for ', races[anc]))
    if (length(flipped) > 0){
      beta[flipped, paste0('A1.', anc)] = ref[flipped, 'A1.ref']
      beta[flipped, paste0('A2.', anc)] = ref[flipped, 'A2.ref']
      for (t in (n.tuning*(anc-1) + 1) : (n.tuning*anc)){
        beta[flipped,paste0('BETA',t)] = - beta[flipped,paste0('BETA',t)]
      }
    }
  }
  beta_mussel = cbind(ref[, c('CHR.ref', 'SNP', 'A1.ref', 'A2.ref')], beta[,paste0('BETA',1:(n.tuning*length(races)))]) # other files: SNP	CHR	A1	BETA1	BETA2	A2
  colnames(beta_mussel)[c(1,3,4)] = c('CHR', 'A1', 'A2')
  
  # ------------------- Match alleles with GWAS summary data for each ancestry separately ------------------- 
  for (race in races){
    trait_name = paste0(race,'_',trait)
    
    stateval = bigreadr::fread2(paste0(output_path, trait_name, '.gwas.ite', ite, '.txt'))
    stateval = stateval[, c('CHR','SNP','A1','A2')]
    colnames(stateval) = c('CHR','SNP', 'A1.ref','A2.ref')
    stateval = left_join(stateval, beta_mussel, by = 'SNP')
    
    na.ind = which(is.na(stateval$A1))
    if (length(na.ind) > 0){
      stateval[na.ind, paste0('BETA',1:(n.tuning*length(races)))] = 0
      stateval[na.ind,'A1'] = stateval[na.ind,'A1.ref']
      stateval[na.ind,'A2'] = stateval[na.ind,'A2.ref']
    }
    flipped = which(stateval$A1.ref != stateval$A1)
    print(paste0(length(flipped), ' flipped SNPs.'))
    if (length(flipped) > 0){
      stateval[flipped,'A1'] = stateval[flipped,'A1.ref']
      stateval[flipped,'A2'] = stateval[flipped,'A2.ref']
      for (t in 1:(n.tuning*length(races))){
        stateval[flipped,paste0('BETA',t)] = - stateval[flipped,paste0('BETA',t)]
      }
    }
    score = stateval[,c('CHR.x','SNP','A1','A2',paste0('BETA',1:(n.tuning*length(races))))] # other files: SNP	CHR	A1	BETA1	BETA2	A2
    colnames(score)[1] = 'CHR'
    write_delim(score, paste0(input_path, trait_name,'.',method, '.full.txt'), delim = '\t')
    rm(stateval,score)
  }

  for (race in races){
    trait_name = paste0(race,'_',trait)
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
    
    ##### Model tuning:
    r2 = bigreadr::fread2(paste0(output_path_eval,trait_name, '.', method, '.txt'))
    r2.avg = colMeans(r2)
    # params.tuned = which.max(r2.avg)
    r2.order = list()
    n.candidates = ncol(r2) # min(15, ncol(r2))
    for (kk in 1:k) r2.order[[kk]] = order(as.numeric(r2[kk,]),decreasing = T)[1:n.candidates]
    # Select the top parameter settings:
    # params.tuned = as.numeric(substr(names(sort(r2.avg, decreasing = T)[1:5]), 5, 10))
    params.tuned = Reduce(intersect, r2.order)
    nonzero.indx = which(r2.avg != 0)
    params.tuned = unique(params.tuned[(params.tuned <= length(r2.avg)) & (params.tuned %in% nonzero.indx)])
    
    # load rm.indx
    load(paste0(workdir, 'PRS_model_training/', method, '/', trait_name, '_', method, '_rm.indx.RData'))
    params.tuned = setdiff(params.tuned, rm.indx)
    a = cbind(sapply(1:(ncol(beta_mussel)-4),function(x){sum(beta_mussel[,paste0('BETA', x)]^2)}), 
              sapply(1:(ncol(beta_mussel)-4),function(x){sum(abs(beta_mussel[,paste0('BETA', x)]) > 0.1)}),
              sapply(1:(ncol(beta_mussel)-4),function(x){sum(abs(beta_mussel[,paste0('BETA', x)]) != 0)}),
              t(r2))
    keep.indx = as.numeric(substr(names(a[((a[,1]< 1.5) & (a[,4]< 1)),4]), 5,7))
    if (sum(params.tuned %in% keep.indx) > 0) params.tuned = params.tuned[params.tuned %in% keep.indx]
    
    nonzero.indx = which(sapply(params.tuned, function(x){sum(abs(beta_mussel[,paste0('BETA',x)])>1e-7)>0}))
    if (length(nonzero.indx) > 0){
      params.tuned = params.tuned[nonzero.indx]
      params.tuned = params.tuned[1:(min(5,length(params.tuned)))]
      tuned.parameters.file = paste0(workdir, 'PRS_model_training/',method,'/',method,'_tuned_parameters_',trait_name,'.RData')
      save(params.tuned, file = tuned.parameters.file)
      # Save the R2s on testing data for the selected tuning parameter settings:
      pumas.cor2 = r2[,params.tuned]
      write.table(pumas.cor2,paste0(PennPRS_finalresults_path,trait_name,".",method,".testing.txt"),col.names = T,row.names=F,quote=F,sep="\t")
      
      
      # Save the best PROSPER model:
      params.tuned.full = params.tuned
      params.tuned = params.tuned[1]
      params = NULL
      for (ra in races){
        tem = cbind(bigreadr::fread2(paste0(prsdir, trait,'.full/tmp/MUSS_beta_in_all_settings_bychrom/settings_1.txt')), ra)
        colnames(tem)[ncol(tem)] = 'ancestry'
        params = rbind(params, tem)
      }
      optimal.pars = params[params.tuned,]
      # write_delim(optimal.pars, file = paste0(PennPRS_finalresults_path, trait_name, '_', method, '_single_PRS_optimal_param.txt'), delim = '\t')
      
      beta_mussel_single = data.frame(beta_mussel[,c('CHR', 'SNP', 'A1', 'A2', paste0('BETA', params.tuned))])
      colnames(beta_mussel_single) = c('CHR', 'SNP','A1', 'A2', 'BETA')
      
      nonzero = (sum(beta_mussel_single$BETA!=0)>0)
      if (nonzero){
        # Remove SNPs with zero effects:
        beta_mussel_single = beta_mussel_single[beta_mussel_single$BETA != 0,]
        prs_mussel_single_outputfile = paste0(PennPRS_finalresults_path, trait_name,'.',method,'.PRS_single_best.txt')
        write_delim(beta_mussel_single, file = prs_mussel_single_outputfile, delim='\t')
        print(paste0('Optimal parameter setting: '))
        # optimal.pars[1,] = signif(optimal.pars[1,],3)
        print(optimal.pars)
        print(paste0('Complete training ', method, ' for ',race))
      }
      if (!nonzero){
        print(paste0('Trained PRS model based on ', method, 'has zero effect estimate for all SNPs. Please try other methods.'))
      }
    }
    if (length(nonzero.indx) == 0){
      print(paste0('Trained PRS model based on ', method, 'has zero effect estimate for all SNPs. Please try other methods.'))
    }
    
    # ------------------------------------------------------------------------------
    # --------------------- Step 6: Train Ensemble PROSPER PRS ---------------------
    # ------------------------------------------------------------------------------
    pumascode = paste(paste0('Rscript ', PUMAS_path, 'PUMAS.evaluation.customized_MA_single_method_ensemble_2subsamples.R'),
                      paste0('--k ',k), 
                      paste0('--ref_path ', eval_ld_ref),
                      paste0('--trait_name ', trait_name),
                      paste0('--prs_method ', method),
                      paste0('--optimal_params_path ', workdir, 'PRS_model_training/',method,'/'),
                      paste0('--xty_path ', xty_path),
                      paste0('--stats_path ', stats_path),
                      paste0('--weight_path ', input_path),
                      paste0('--output_path_eval ', output_path_eval),
                      paste0('--output_path ', PennPRS_finalresults_path))
    system(pumascode)
    # Generate ensemble PRS file:
    beta_prosper_ensemble = beta_mussel[,c('CHR', 'SNP', 'A1', ' A2', paste0('score',params.tuned.full))]
    colnames(beta_prosper_ensemble) = c('CHR', 'SNP','A1','A2', paste0('BETA',params.tuned.full))
    beta_prosper_ensemble0 = beta_prosper_ensemble
    
    # 1. Ensemble using PUMA-CUBS
    ensemble.weights = bigreadr::fread2(paste0(PennPRS_finalresults_path,trait_name,'.',method,".omnibus.weights.txt"))
    ensemble.weights = colMeans(ensemble.weights)
    ensemble.weights = ensemble.weights/sum(ensemble.weights)
    beta_prosper_ensemble$BETA = as.numeric(as.matrix(beta_prosper_ensemble[,paste0('BETA',params.tuned.full)]) %*% t(t(ensemble.weights)))
    beta_prosper_ensemble = beta_prosper_ensemble[, c('CHR','SNP','A1','A2','BETA')]
    prs_prosper_outputfile = paste0(PennPRS_finalresults_path, trait_name,'.',method,'.PRS.txt')
    write_delim(beta_prosper_ensemble, file = prs_prosper_outputfile, delim='\t')
  }
  
  if ( opt$verbose >= 1 ){
    print(paste0('*************************************'))
    print(paste0('******* Complete training ', method, ' *******'))
    print(paste0('*************************************'))
  }
}







# Save the final PRS models generated by different methods to the working directory
for (method in methods){
  for (race in races){
    trait_name = paste0(race,'_',trait)
    system(paste('cp -r', paste0(PennPRS_finalresults_path, trait_name,'.',method,'.PRS.txt'), workdir))
    if (method == 'PROSPER') system(paste('cp -r', paste0(PennPRS_finalresults_path, trait_name,'.',method,'.PRS_single_best.txt'), workdir))
  }
}

filen<-paste0(workdir, 'PRS_model_training_info.txt')
file.create(filen)
zz <- file(filen, "w")

print.title = paste0("Summary of PRS model Training on ",trait, " for ", race)
cat(paste0("\n",paste(rep('*', nchar(print.title)+10),collapse='')), file = zz)
cat(paste0("\n**** ", print.title, " ****"), file = zz)
cat(paste0("\n",paste(rep('*', nchar(print.title)+10),collapse=''),'\n'), file = zz)


if ('PROSPER' %in% methods){
  method = 'PROSPER'
  prsdir = paste0(prsdir0, method,'/') 
  cat(paste0("\n\n****************************************************"), file = zz)
  cat(paste0("\n******** PROSPER (January 14, 2024 Version) ********"), file = zz)
  cat(paste0("\n****************************************************"), file = zz)
  cat(paste0("\n* Please refer to https://github.com/Jingning-Zhang/PROSPER for details of PROSPER."), file = zz)
  cat(paste0("\n\n************ Step 1: lassosum2 by race *************"), file = zz)
  cat(paste0('\n************ Tuning parameter settings: ************'), file = zz)
  cat(paste0('\nnLambda = ', opt$nlambda, ' (number of candidate values for the shrinkage parameter in the L1 regularization)'), file = zz)
  cat(paste0('\nnDelta = ', opt$ndelta, ' (number of candidate values for the shrinkage parameter in the L2 regularization)'), file = zz)
  cat(paste0('\nlambda.min.ratio = ', opt$lambda.min.ratio, ' (ratio between the lowest and highest candidate values of lambda)'), file = zz)
  cat(paste0('\n************** Tuned parameter values: *************'), file = zz)
  for (race in races){
    cat(paste0('\n********************** ',race,' *********************'), file = zz)
    trait_name = paste0(race,'_',trait)
    pars.lassosum2 = bigreadr::fread2(paste0(path_out_lassosum2, '/', race, '/', trait_name, '_optimal_param_full.txt'))
    cat(paste0('\nLambda = ', signif(as.numeric(pars.lassosum2[1,'lambda']),3)), file = zz)
    cat(paste0('\nDelta = ', signif(as.numeric(pars.lassosum2[1,'delta']),3)), file = zz)
  }
  
  cat(paste0("\n\n****************** Step 2: PROSPER *****************"), file = zz)
  cat(paste0('\n************ Tuning parameter settings: ************'), file = zz)
  cat(paste0('\nnL1 = ', 5, ' (number of candidate values for the shrinkage parameter in the L1 regularization for SNP effect sizes)'), file = zz)
  cat(paste0('\nnL2 = ', 5, ' (number of candidate values for the shrinkage parameter in the L2 regularization for inducing similarity across ancestries)'), file = zz)
  cat(paste0('\nNumber of top-performing PRS models combined in the final PROSPER PRS models: ', 15), file = zz)
}


if ('PRS-CSx' %in% methods){
  method = 'PRS-CSx'
  prsdir = paste0(prsdir0, method,'/') 
  cat(paste0("\n\n**********************************************"), file = zz)
  cat(paste0("\n******** PRS-CSx (May 14, 2024 Version) ********"), file = zz)
  cat(paste0("\n************************************************"), file = zz)
  cat(paste0("\n* Please refer to https://github.com/getian107/PRScsx for details of PRS-CSx."), file = zz)
  cat(paste0('\n************ Tuning parameter settings: ************'), file = zz)
  cat(paste0('\nphi = ', opt$phi, ' (global shrinkage parameter)'), file = zz)
  cat(paste0('\n************** Tuned parameter values: *************'), file = zz)
  for (race in races){
    cat(paste0('\n********************** ',race,' *********************'), file = zz)
    trait_name = paste0(race,'_',trait)
    load(paste0(workdir, 'PRS_model_training/',method,'/tuned_parameters_',trait_name,'.RData'))
    cat(paste0('\nphi = ', phi.vals[params.tuned[1]]), file = zz)
  }
}

if ('MUSSEL' %in% methods){
  method = 'MUSSEL'
  prsdir = paste0(prsdir0, method,'/') 
  cat(paste0("\n\n***********************************************"), file = zz)
  cat(paste0("\n******** MUSSEL (April 10, 2024 Version) ********"), file = zz)
  cat(paste0("\n*************************************************"), file = zz)
  cat(paste0("\n* Please refer to https://github.com/Jin93/MUSSEL for details of MUSSEL"), file = zz)
  cat(paste0('\n************ Tuning parameter settings: ************'), file = zz)
  cat(paste0('\nnLambda = ', opt$nlambda, ' (number of candidate values for the shrinkage parameter in the L1 regularization)'), file = zz)
  cat(paste0('\nnDelta = ', opt$ndelta, ' (number of candidate values for the shrinkage parameter in the L2 regularization)'), file = zz)
  cat(paste0('\nlambda.min.ratio = ', opt$lambda.min.ratio, ' (ratio between the lowest and highest candidate values of lambda)'), file = zz)
  cat(paste0('\n************** Tuned parameter values: *************'), file = zz)
  for (race in races){
    cat(paste0('\n********************** ',race,' *********************'), file = zz)
    trait_name = paste0(race,'_',trait)
    pars.lassosum2 = bigreadr::fread2(paste0(path_out_lassosum2, '/', race, '/', trait_name, '_optimal_param_full.txt'))
    cat(paste0('\nLambda = ', signif(as.numeric(pars.lassosum2[1,'lambda']),3)), file = zz)
    cat(paste0('\nDelta = ', signif(as.numeric(pars.lassosum2[1,'delta']),3)), file = zz)
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
cat(paste0('\n* PRS models trained by multi-ancestry methods:'), file = zz)
for (method in methods){
  cat(paste0('\n* ', method, ':'), file = zz)
  for (race in races){
    cat(paste0('\n********************** ',race,' *********************'), file = zz)
    trait_name = paste0(race,'_',trait)
    prsfile = paste0(workdir, trait_name,'.',method, '.PRS.txt')
    if (file.exists(prsfile)) cat(paste0('\n  ', trait_name,'.',method, '.PRS.txt'), file = zz)
  }
}
cat(paste0('\n\n* Details of the PRS training:'), file = zz)
cat(paste0('\n  PRS_model_training_info.txt\n'), file = zz)
close(zz)


# Clean up intermediate files:
system(paste0('rm -rf ', paste0(workdir, 'input_for_eval/')))
system(paste0('rm -rf ', paste0(workdir, 'sumdata/')))
system(paste0('rm -rf ', paste0(workdir, 'output/')))
system(paste0('rm -rf ', paste0(workdir, 'output_for_eval/')))
system(paste0('rm -rf ', paste0(workdir, 'PRS_model_training/')))
system(paste0('rm -rf ', paste0(workdir, 'PennPRS_results/')))







