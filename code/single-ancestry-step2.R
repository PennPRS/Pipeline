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
# userID = opt$userID
submissionID = opt$submissionID
methods = str_split(opt$methods,",")[[1]]
trait = opt$trait
race = opt$race
LDrefpanel = opt$LDrefpanel
# Parameters for subsampling
k = opt$k
# if multiple methods are selected, also input ensemble into opt

# ----------------
ensemble = FALSE
if (length(methods) > 1) ensemble = opt$ensemble

# Optional input parameters:
partitions <- opt$partitions

ld_path <- paste0(PennPRS_path, '/LD/', race, '/')
PUMAS_path = paste0(PennPRS_path,'/code/')
trait_name = paste0(race,'_',trait)
if (LDrefpanel == '1kg'){
  eval_ld_ref_path <- paste0(ld_path, '/1KGref_plinkfile/') # set to the /1KGref_plinkfile folder under /LD/
} 
jobID = paste(c(trait,race, paste0(methods,collapse = '.'), submissionID), collapse = '_')
workdir = paste0(homedir,jobID,'/')
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
eval_ld_ref = paste0(eval_ld_ref_path,LDrefpanel,'_hm3_',race,'_ref')
# Create a separate directory 'PRS_model_training/' to store input for training PRS models
prsdir0 = paste0(workdir, 'PRS_model_training/')
output_path_eval = paste0(workdir, 'output_for_eval/')

ensemble.methods = NULL
load(paste0(PennPRS_finalresults_path, 'step1.RData'))

# Read in LD reference file to add SNP position info for output PRS files:
ref.bim = bigreadr::fread2(paste0(eval_ld_ref, '.bim'))[, c(2,4)]
colnames(ref.bim) = c('SNP', 'chr_position')

# ---------------------------------------------------------------------------------------
# --------------------- Step 5: Train Ensemble PRS and Benchmarking ---------------------
# ---------------------------------------------------------------------------------------
# If more than one method is selected, can train an ensemble PRS
if ((ensemble) & (length(ensemble.methods) > 1)){
  xty_path = stats_path = output_path
  pumascode = paste(paste0('Rscript ', PUMAS_path, 'PUMAS.evaluation.customized_single_ans_ensemble_2subsamples.R'),
                    paste0('--k ',k), 
                    paste0('--ref_path ', eval_ld_ref),
                    paste0('--trait_name ', trait_name),
                    paste0('--prs_method ', paste0(ensemble.methods, collapse = ',')),
                    paste0('--optimal_params_path ', prsdir0),
                    paste0('--xty_path ', xty_path),
                    paste0('--stats_path ', stats_path),
                    paste0('--weight_path ', input_path),
                    paste0('--output_path_eval ', output_path_eval),
                    paste0('--output_path ', PennPRS_finalresults_path))
  system(pumascode)
  
  # Generate ensemble PRS file:
  n.snps = rep(0, length(methods))
  score = list()
  for (m in 1:length(methods)){
    method = methods[m]
    score.file = paste0(PennPRS_finalresults_path, trait_name,'.',method,'.PRS.txt')
    if (file.exists(score.file)){
      score[[m]] = bigreadr::fread2(score.file) # only contains non-zero-effect SNPs
      n.snps[m] = nrow(score[[m]])
    }
    if (!file.exists(score.file)) score[[m]] = matrix()
  }
  m.indx = which.max(n.snps)
  if (max(n.snps) == 0) cat(paste0('None of the selected methods generated a PRS model with nonzero SNP effects. Please consider the possibility that the outcome is not heritable.'))
  if (max(n.snps) > 0){
    SCORE = score[[m.indx]]; colnames(SCORE)[which(colnames(SCORE) == 'BETA')] = paste0('BETA.',methods[m.indx])
    for (m in c(1:length(methods))[-m.indx]){
      if (sum(is.na(score[[m]])) == 0) {
        tem = score[[m]][,c('SNP','BETA')]; colnames(tem)[2] = paste0('BETA.',methods[m])
        SCORE = left_join(SCORE, tem, by="SNP")
      }
      if (sum(is.na(score[[m]])) != 0){
        SCORE = cbind(SCORE, 0)
        colnames(SCORE)[ncol(SCORE)] = paste0('BETA.',methods[m])
      }    
    }
    SCORE[is.na(SCORE)] = 0
    SCORE = SCORE[,c('CHR','SNP','A1','A2', paste0('BETA.',ensemble.methods))]
    SCORE0 = SCORE
    
    # 1. Ensemble using PUMA-CUBS
    ensemble.weights = bigreadr::fread2(paste0(PennPRS_finalresults_path, trait_name, '.', paste0(ensemble.methods, collapse = '.'), '.omnibus.weights.txt'))
    ensemble.weights = colMeans(ensemble.weights)
    ensemble.weights = ensemble.weights/sum(ensemble.weights)
    ensemble.weights[is.na(ensemble.weights)] = 0
    ensemble.weights1 = ensemble.weights
    SCORE$BETA = as.numeric(as.matrix(SCORE[,paste0('BETA.',ensemble.methods)]) %*% t(t(ensemble.weights)))
    SCORE = SCORE[, c('CHR','SNP','A1','A2','BETA')]
    
    # Update output file format according to the pgsc_calc pipeline from the PGS Catalog
    SCORE = merge(SCORE, ref.bim, by = 'SNP')
    SCORE = SCORE[, c('CHR','chr_position','A1','A2', 'BETA')]
    colnames(SCORE) = c('chr_name','chr_position','effect_allele','other_allele', 'effect_weight')
    
    prs_ensemble_outputfile = paste0(workdir, trait_name,'.ensemble_PRS.txt')
    write_delim(SCORE, prs_ensemble_outputfile, delim='\t')
    
    # 2. Ensemble using an alternative approach
    ensemble.weights = bigreadr::fread2(paste0(PennPRS_finalresults_path, trait_name, '.', paste0(ensemble.methods, collapse = '.'), '.omnibus.weights.alternative.txt'))
    ensemble.weights = colMeans(ensemble.weights)
    ensemble.weights = ensemble.weights/sum(ensemble.weights)
    ensemble.weights[is.na(ensemble.weights)] = 0
    ensemble.weights2 = ensemble.weights
    SCORE = SCORE0
    SCORE$BETA = as.numeric(as.matrix(SCORE[,paste0('BETA.',ensemble.methods)]) %*% t(t(ensemble.weights)))
    SCORE = SCORE[, c('CHR','SNP','A1','A2','BETA')]
    
    # Update output file format according to the pgsc_calc pipeline from the PGS Catalog
    SCORE = merge(SCORE, ref.bim, by = 'SNP')
    SCORE = SCORE[, c('CHR','chr_position','A1','A2', 'BETA')]
    colnames(SCORE) = c('chr_name','chr_position','effect_allele','other_allele', 'effect_weight')
    
    prs_ensemble_outputfile = paste0(workdir, trait_name,'.ensemble_PRS_by_model_averaging.txt')
    write_delim(SCORE, prs_ensemble_outputfile, delim='\t')
    
    if ( opt$verbose >= 1 ){
      print.title = paste0('Complete training Ensemble PRS combining single PRS trained by ', paste0(ensemble.methods,collapse = ', '))
      cat(paste0("\n",paste(rep('*', nchar(print.title)+10),collapse='')))
      cat(paste0("\n**** ", print.title, " ****"))
      cat(paste0("\n",paste(rep('*', nchar(print.title)+10),collapse=''),'\n'))
    }
  }
}

# Save the final PRS models generated by each single methods to the working directory
for (method in methods){
  tfile = paste0(PennPRS_finalresults_path, trait_name,'.',method,'.PRS.txt')
  if (file.exists(tfile)) {
    SCORE = bigreadr::fread2(tfile)
    # Update output file format according to the pgsc_calc pipeline from the PGS Catalog
    SCORE = merge(SCORE, ref.bim, by = 'SNP')
    SCORE = SCORE[, c('CHR','SNP','chr_position','A1','A2', 'BETA')]
    colnames(SCORE) = c('chr_name','rsid','chr_position','effect_allele','other_allele', 'effect_weight')
    write_delim(SCORE, paste0(workdir, trait_name,'.',method,'.PRS.txt'))
  }
}





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
write_line("commands for calculating PRS using the trained models.")

# --------------------------------------------------
# Training summary
# --------------------------------------------------
write_section("2. SUMMARY OF PRS MODEL TRAINING")
methods.completed = NULL

if ('C+T' %in% methods){
  method = 'C+T'
  prsdir = paste0(prsdir0, method,'/')
  tfile = paste0(workdir, trait_name,'.',method,'.PRS.txt')
  
  write_subsection("C+T")
  write_line("* Reference documentation:")
  write_line("  https://pmc.ncbi.nlm.nih.gov/articles/PMC6904799/")
  
  if (file.exists(tfile)){
    methods.completed = c(methods.completed, method)
    write_line("* Model training status: completed.")
    write_line(paste0("* Output score file: ", trait_name, ".", method, ".PRS.txt"))
    write_line("* Tuning parameter settings:")
    
    tuning.file = paste0(workdir, 'PRS_model_training/',method,'/tuned_parameters_',trait_name,'.RData')
    if (file.exists(tuning.file)){
      load(tuning.file)
      r2.ct = r2; pval.ct = pvalthr
      write_line('  - r2 = ', opt$R2, ' (SNPs having squared correlation higher than r2 with the index SNPs will be removed)')
      write_line('  - P-value Threshold = ', opt$Pvalthr, ' (p-value threshold for the thresholding step)')
      write_line('* Tuned parameter values:')
      write_line('  - r2 = ', r2.ct)
      write_line('  - P-value Threshold = ', pval.ct)
    }
    if (!file.exists(tuning.file)){
      write_line(' ', method, ' failed to run, no PRS model was generated.\n')
    }
  }
  if (!file.exists(tfile)){
    write_line("* Model training status: no PRS model was generated.")
    write_line("* Please check log file for details.")
  }
  if (err.CT == 1){
    if (ensemble){
      write_line('  [Warning] Trained PRS model has an estimated R2 < 0, indicating that the PRS lacks prediction power and is thus not incorporated in the ensemble PRS.')
    } 
    if (!ensemble){
      write_line('  [Warning] Trained PRS model has an estimated R2 < 0, indicating that the PRS lacks prediction power.')
    }
    write_line("  Possible explanations include:")
    write_line("    1. The trait has low heritability.")
    write_line("    2. The GWAS summary statistics have insufficient power (e.g., due to small sample size).")
    write_line("    3. The input GWAS summary data contains problematic values (e.g., BETA or SE).")
    write_line("    4. C+T is not powerful for developing PRS for the trait, and other methods can be considered instead.")
  }
  if (err.CT == 2){
    write_line('  [Warning] C+T failed to select any set of independent SNPs with p-value < ',max(Pvalthr), '.')
  }
}


if ('lassosum2' %in% methods){
  method = 'lassosum2'
  prsdir <- paste0(prsdir0, method, "/")
  tfile <- paste0(workdir, trait_name, ".", method, ".PRS.txt")
  
  write_subsection("lassosum2 (June 8, 2023 Version)")
  write_line("* Reference documentation:")
  write_line("  https://privefl.github.io/bigsnpr/articles/LDpred2.html")
  
  tuning.file = paste0(PennPRS_finalresults_path, trait_name,'.',method,'.optimal_params.txt')
  if ((file.exists(tfile)) & (file.exists(tuning.file))){
    methods.completed = c(methods.completed, method)
    write_line("* Model training status: completed")
    write_line(paste0("* Output score file: ", trait_name, ".", method, ".PRS.txt"))
    write_line("* Tuning parameter settings:")
    
    lassosum2.out = bigreadr::fread2(tuning.file)
    write_line('  - nlambda = ', opt$nlambda, ' (number of candidate values for the shrinkage parameter in the L1 regularization)')
    write_line('  - Delta = ', opt$delta, ' (shrinkage parameter in the L2 regularization)')
    write_line('  - lambda.min.ratio = ', opt$lambda.min.ratio, ' (ratio between the lowest and highest candidate values of lambda)')
    write_line('* Tuned parameter values:')
    write_line('  - Lambda = ', as.numeric(lassosum2.out[1,1]))
    write_line('  - Delta = ', as.numeric(lassosum2.out[1,2]),'\n')
  }
  if (!file.exists(tuning.file)){
    write_line("* Model training status: no PRS model was generated.")
    write_line("* Please check log file for details.")
  }
  if (err.lassosum2 == 1){
    if (ensemble) write_line('  [Warning] Trained PRS model has an estimated R2 < 0, indicating that the PRS lacks prediction power and is thus not incorporated in the ensemble PRS.')
    if (!ensemble) write_line('  [Warning] Trained PRS model has an estimated R2 < 0, indicating that the PRS lacks prediction power.')
    write_line("  Possible explanations include:")
    write_line("    1. The trait has low heritability.")
    write_line("    2. The GWAS summary statistics have insufficient power (e.g., due to small sample size).")
    write_line("    3. The input GWAS summary data contains problematic values (e.g., BETA or SE).")
    write_line("    4. C+T is not powerful for developing PRS for the trait, and other methods can be considered instead.")
  }
  if (err.lassosum2 == 2){
    write_line('  No tuning parameter setting led to a PRS model with nonzero SNP effect estimates.')
    write_line('  This is possibly due to convergence issues of the ', method, ' algorithm on the input GWAS data.')
  }
}


if ('LDpred2' %in% methods){
  method = 'LDpred2'
  prsdir <- paste0(prsdir0, method, "/")
  tfile <- paste0(workdir, trait_name, ".", method, ".PRS.txt")
  
  write_subsection("LDpred2 (June 8, 2023 Version)")
  write_line("* Reference documentation:")
  write_line("  https://privefl.github.io/bigsnpr/articles/LDpred2.html")
  
  tuning.file = paste0(PennPRS_finalresults_path, trait_name,'.',method,'.optimal_params.txt')
  if ((file.exists(tfile)) & (file.exists(tuning.file))){
    methods.completed = c(methods.completed, method)
    write_line("* Model training status: completed")
    write_line(paste0("* Output score file: ", trait_name, ".", method, ".PRS.txt"))
    write_line("* Tuning parameter settings:")
    
    ldpred2.out = bigreadr::fread2(tuning.file)
    write_line('  - h2.ratio = ', opt$alpha, ' (heritability = h2.ratio * H_20, where H_20 is the heritability estimated by LD score regression)')
    write_line('  - p = ', opt$p_seq, ' (causal SNP proportion)')
    write_line('  - sparse = ', opt$sparse, ' (whether to generate a sparse PRS model)')
    write_line('* Tuned parameter values:')
    write_line('  - p = ', ldpred2.out[1,1])
    write_line('  - H2 = ', ldpred2.out[1,2], ' (heritability)')
    write_line('  - sparse = ', ldpred2.out[1,3])
  }
  if (!file.exists(tuning.file)){
    write_line("* Model training status: no PRS model was generated.")
    write_line("* Please check log file for details.")
  }
  if (err.LDpred2 == 1){
    if (ensemble) write_line('  [Warning] Trained PRS model has an estimated R2 < 0, indicating that the PRS lacks prediction power and is thus not incorporated in the ensemble PRS.')
    if (!ensemble) write_line('  [Warning] Trained PRS model has an estimated R2 < 0, indicating that the PRS lacks prediction power.')
    write_line("  Possible explanations include:")
    write_line("    1. The trait has low heritability.")
    write_line("    2. The GWAS summary statistics have insufficient power (e.g., due to small sample size).")
    write_line("    3. The input GWAS summary data contains problematic values (e.g., BETA or SE).")
    write_line("    4. C+T is not powerful for developing PRS for the trait, and other methods can be considered instead.")
  }
  if (err.LDpred2 == 2){
    write_line('  No tuning parameter setting led to a PRS model with nonzero SNP effect estimates.')
    write_line('  This is possibly due to convergence issues of the ', method, ' algorithm on the input GWAS data.')
  }
}



if ((ensemble) & (length(ensemble.methods) > 1)){
  write_subsection("Ensemble PRS")
  write_line("* Weights of PRS by each method:")
  
  for (m in 1:length(ensemble.methods)){
    method = ensemble.methods[m]
    write_line(method, ': ', signif(ensemble.weights1[m],3))
  }
  if (sum(ensemble.weights1 != 0) == 0) write_line('##### Note: weights of all PRS are zero due to insufficient prediction power of all PRS.')
  
  write_subsection("Alternative: Optimal PRS by Model Averaging")
  write_line("* Weights of PRS by each method:")
  for (m in 1:length(ensemble.methods)){
    method = ensemble.methods[m]
    write_line(method, ': ', signif(ensemble.weights2[m],3))
  }
  if (sum(n.snps == 0)>0){
    write_line('##### Note: ', paste(methods[which(n.snps == 0)], collapse = ', '), ' did not generate a PRS model that has non-zero SNP weights or sufficient prediction power and thus did not contribute to the ensemble PRS.')
    write_line('#####       This can be due to low heritability of the trait, insufficient GWAS training sample size, or convergence issue of the algorithms.')
  }
  if (sum(ensemble.weights2 != 0) == 0) write_line('##### Note: weights of all PRS are zero due to minimal prediction power of all PRS.')
}

if ((ensemble) & (length(ensemble.methods) == 1)){
  write_subsection("Ensemble PRS")
  write_line('##### Note: no ensemble PRS model was generated because only one of the selected methods generated a PRS model that has non-zero SNP weights or sufficient prediction power.')
}

if ((ensemble) & (length(ensemble.methods) == 0)){
  write_subsection("Ensemble PRS")
  write_line('##### Note: no ensemble PRS model was generated because none of the selected methods generated a PRS model that has non-zero SNP weights or sufficient prediction power.')
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
  # write_line("Example PLINK2 command:")
  # write_line("")
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
cat(paste0('\n  QC_report.txt\n'), file = zz)

cat(paste0('\n* PRS models trained by single-ancestry methods:'), file = zz)
for (method in methods){
  prsfile = paste0(workdir, trait_name,'.',method, '.PRS.txt')
  if (file.exists(prsfile)) cat(paste0('\n  ', trait_name,'.',method, '.PRS.txt'), file = zz)
}
if ((ensemble) & (length(ensemble.methods) > 1)){
  cat(paste0('\n\n* Ensemble PRS models combining PRS models trained by different methods:'), file = zz)
  ensembleprsfile1 = paste0(workdir, trait_name, '.ensemble_PRS.txt')
  if (file.exists(ensembleprsfile1)) cat(paste0('\n  ', trait_name, '.ensemble_PRS.txt'), file = zz)
  ensembleprsfile2 = paste0(workdir, trait_name, '.ensemble_PRS_by_model_averaging.txt')
  if (file.exists(ensembleprsfile2)) cat(paste0('\n  ', trait_name, '.ensemble_PRS_by_model_averaging.txt\n'), file = zz)
} 

cat(paste0('\n\n* Details of the PRS training:'), file = zz)
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
