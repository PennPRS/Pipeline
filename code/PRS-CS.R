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
  make_option("--methods", action = "store", default = 'PRS-CS', type = "character",
              help="Options: PRS-CS"),
  make_option("--type", action = "store", default = 'auto', type = "character",
              help="Method type: auto (PRS-CS auto, default) or grid (PRS-CS) [Optional]"),
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

  make_option("--phi", action = "store", default = NA, type = "character",
              help="Global shrinkage parameter phi. In PRS-CS-auto, phi is not specified and will be learnt from the data using a fully Bayesian approach. 
              This usually works well for polygenic traits with large GWAS sample sizes (hundreds of thousands of subjects). 
              For GWAS with limited sample sizes (including most of the current disease GWAS), using PRS-CS-auto with phi set to 1e-2 
              (for highly polygenic traits) or 1e-4 (for less polygenic traits), or using PRS-CS with a small-scale grid search 
              to find the optimal phi value in the validation dataset often improves perdictive performance.
              Recommended alternative (default values in the PRS-CS algorithm, May 14, 2024 version): 1e+00,1e-02,1e-04,1e-6.
              Options: candidate values in (0,1], divided by comma [Optional]"),
  make_option("--N_THREADS", action = "store", default = 5, type = "numeric",
              help="Number of threads used for parallel computing of PRS-CS across chromosomes.
              Default: 5.
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
type = str_split(opt$type,",")[[1]]
trait = opt$trait
race = opt$race
LDrefpanel = opt$LDrefpanel
# Parameters for subsampling
k = opt$k

# Optional input parameters:
partitions <- opt$partitions

ld_path <- paste0(PennPRS_path, '/LD/', race, '/')
PUMAS_path = paste0(PennPRS_path,'/code/')
PRScs_path = paste0(PennPRS_path, 'software/PRScs/')
plink_path = paste0(PennPRS_path, 'software/')
PATH_TO_REFERENCE = paste0(PennPRS_path,'/LD/', race, '/ldblk_1kg_',tolower(race)) # We can use 1000 Genomes reference data for now
threads = 1


trait_name = paste0(race,'_',trait)
ld_path0 <- paste0(ld_path, 'LD_1kg/') # set to the /LD_1kg folder under /LD/
if (LDrefpanel == '1kg'){
  eval_ld_ref_path <- paste0(ld_path, '/1KGref_plinkfile/') # set to the /1KGref_plinkfile folder under /LD/
} 
# Job name/ID: e.g., trait_race_method_userID_submissionID
jobID = paste(c(trait,race, paste0(methods,collapse = '.'), submissionID), collapse = '_')
# Create a job-specific (trait, race, methods, userID, jobID) directory to save all the outputs, set the working directory to this directory
workdir = paste0(homedir,jobID,'/')
suppressWarnings(dir.create(workdir))
setwd(workdir) 


# Parameters
if ('PRS-CS' %in% methods){
  PHI = opt$phi; phi.vals = as.numeric(str_split(PHI,",")[[1]])
  N_THREADS = as.numeric(opt$N_THREADS)
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

if ('grid' %in% type){
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
    cat(paste0("\n**************************************************************"))
    cat(paste0("\n**** Step 2: Train PRS models on pseudo training datasets ****"))
    cat(paste0("\n**************************************************************\n"))
  }
  
  if ('PRS-CS' %in% methods){
    cat(paste0("\n**********************************************************"))
    cat(paste0("\n**** Start running PRS-CS on pseudo training datasets ****"))
    cat(paste0("\n**********************************************************\n"))
    method = 'PRS-CS'
    prsdir = paste0(prsdir0, method,'/')
    # --------------------- Input preparation for PRS-CS ---------------------
    # Create a "pseudo" .bim file: since we do not have individual-level validation dataset, 
    # we need to create this .bim file as one of the input files for PRS-CS:
    VALIDATION_BIM_PREFIX = paste0(input_path,'pseudo_validation_',trait_name) # This is supposed to be a real .bim file for validation individuals, but we create a pseudo validation file instead 
    valbim = bigreadr::fread2(paste0(output_path,trait_name,'.gwas.ite1.txt')) # the set of SNPs is the same across the k MCCV files.
    GWAS_SAMPLE_SIZE = as.character(round(mean(valbim$N)))
    valbim = cbind(valbim[,c('CHR', 'SNP')], 0, 1, valbim[,c('A1', 'A2')])
    write_delim(valbim, paste0(VALIDATION_BIM_PREFIX, '.bim'), delim = '\t', col_names = F)
    
    # Reformat summary data to use as the input data for PRS-CS:
    chrs = paste0(1:22, collapse = ',')
    for (ite in 1:k){
      pumasout = paste0(output_path, trait_name, '.gwas.ite', ite, '.txt')
      if (!file.exists(pumasout)) print(paste0('Subsampling failed to generate ', ite,'-th fold of the MCCV summary data. Rerun pumas.subsampling.R.'))
      if (file.exists(pumasout)){
        sumraw0 = bigreadr::fread2(pumasout)
        sumraw0 = sumraw0[,c('CHR', 'SNP', 'A1', 'A2', 'BETA', 'SE')]
        colnames(sumraw0) = c('CHR', 'SNP', 'A1', 'A2', 'BETA', 'SE') # A1: REF
        prscs.sumdat.file = paste0(prsdir, trait_name, '_reformated_gwas.ite', ite, '.txt')
        write_delim(sumraw0[, c('SNP', 'A1', 'A2', 'BETA', 'SE')], prscs.sumdat.file, delim = '\t')
        # for (chr in 1:22){
        #   prscs.sumdat.file = paste0(prsdir, trait_name, '_reformated_gwas_chr', chr, '.ite', ite, '.txt')
        #   sumraw = sumraw0[sumraw0$CHR == chr, c('SNP', 'A1', 'A2', 'BETA', 'SE')]
        #   write_delim(sumraw, prscs.sumdat.file, delim = '\t')
        # }
        print(paste0('Reformatting the input GWAS data for iteration ', ite, ' for ', method, ': completed.'))
        
        out_dir = paste0(prsdir, trait_name,'.',method,'.ite',ite)
        SUM_STATS_FILE = paste0(prsdir, trait_name, '_reformated_gwas.ite', ite, '.txt')
        SEED = 2024
        for (phi in phi.vals){
          prscscode = paste(paste0('python ', PRScs_path, 'PRScs.py'),
                            paste0('--ref_dir=', PATH_TO_REFERENCE),
                            paste0('--bim_prefix=', VALIDATION_BIM_PREFIX),
                            paste0('--sst_file=', SUM_STATS_FILE),
                            paste0('--n_gwas=', GWAS_SAMPLE_SIZE),
                            paste0('--out_dir=', out_dir),
                            paste0('--chrom=', chrs),
                            paste0('--phi=', phi),
                            paste0('--seed=', SEED))
          system(prscscode)
          print(paste0('Complete training ', method, ' for ite ', ite, ' with phi = ', phi))
        }
      }
    }
  }
  
  
  # ---------------------------------------------------------------------------------------
  # ----------------------------- Step 3: Single Model Tuning -----------------------------
  # ---------------------------------------------------------------------------------------
  # single_prs() in PUMA-CUBS.evaluation.R
  if ( opt$verbose >= 1 ) {
    cat(paste0("\n****************************************"))
    cat(paste0("\n******* Step 3: Parameter Tuning *******"))
    cat(paste0("\n****************************************\n"))
  }
  
  if ('PRS-CS' %in% methods){
    method = 'PRS-CS'
    prsdir = paste0(prsdir0, method,'/')
    n.tuning = length(phi.vals)
    
    for (ite in 1:k){
      out_dir = paste0(prsdir, trait_name,'.',method,'.ite',ite)
      SCORE = NULL
      for (phi in phi.vals){
        score = NULL
        for(chr in c(1:22)){
          temfile = paste0(out_dir, '_pst_eff_a1_b0.5_phi', strsplit(PHI,split=',')[[1]][which(phi.vals == phi)], '_chr',chr,'.txt')
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
      scores = stateval[,c('CHR.ref','SNP','A1','A2',paste0('BETA',1:n.tuning))] # other files: SNP	CHR	A1	BETA1	BETA2	A2
      write_delim(scores, paste0(input_path, trait_name,'.',method, '.ite',ite,'.txt'), delim = '\t')
      # write.table(scores, paste0(input_path, trait_name,'.',method, '.ite', ite, '.txt'), row.names = F,col.names = T, quote = FALSE, sep = "\t" )
      rm(stateval, scores)
    }
    
    xty_path = stats_path = output_path # the "output_path" used for storing output from pumas.subsampling.R
    
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
    for (kk in 1:k) r2.order[[kk]] = order(as.numeric(r2[kk,]),decreasing = T)[1:n.tuning]
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
      for (kk in 1:k) r2.order[[kk]] = order(as.numeric(r2[kk, nonzero.indx]),decreasing = T)[1:n.tuning]
      # Select the top parameter settings:
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
    cat(paste0("\n****************************************************************"))
    cat(paste0("\n******* Step 4: Generate best PRS(s) on the orignal GWAS *******"))
    cat(paste0("\n****************************************************************\n"))
  }
  # ---------------------------------------------------------------------------------------
  # ---------- Step 4.1: Train the best PRS(s) from each method on the whole data ---------
  # ---------------------------------------------------------------------------------------
  
  sumraw = bigreadr::fread2(paste0(output_path,trait_name,".gwas_matched.txt"))
  err.PRSCS = 0
  
  if ('PRS-CS' %in% methods){
    method = 'PRS-CS'
    prsdir = paste0(prsdir0, method,'/')
    if ( opt$verbose >= 1 ){
      print(paste0('***************************************************************'))
      print(paste0('******* Start training PRS-CS on the original GWAS data *******'))
      print(paste0('***************************************************************'))
    }
    
    # Reformat summary data to use as the input data for PRS-CS-auto:
    sumraw0 = sumraw[,c('CHR', 'SNP', 'A1', 'A2', 'BETA', 'SE')] # A1: REF
    prscs.sumdat.file = paste0(prsdir,trait_name,'_reformated_gwas.txt')
    write_delim(sumraw[, c('SNP', 'A1', 'A2', 'BETA', 'SE')], prscs.sumdat.file)
    print(paste0('Reformatting the full GWAS data for ', method, ': completed.'))
    
    out_dir = paste0(prsdir, trait_name,'.',method)
    SUM_STATS_FILE = prscs.sumdat.file
    load(paste0(workdir, 'PRS_model_training/',method,'/tuned_parameters_',trait_name,'.RData'))
    SEED = 2024
    phi = phi.vals[params.tuned[1]]
    chrs = paste0(1:22, collapse = ',')
    prscscode = paste(paste0('python ', PRScs_path, 'PRScs.py'),
                      paste0('--ref_dir=', PATH_TO_REFERENCE),
                      paste0('--bim_prefix=', VALIDATION_BIM_PREFIX),
                      paste0('--sst_file=', SUM_STATS_FILE),
                      paste0('--n_gwas=', GWAS_SAMPLE_SIZE),
                      paste0('--out_dir=', out_dir),
                      paste0('--chrom=', chrs),
                      paste0('--phi=', phi),
                      paste0('--seed=', SEED))
    system(prscscode)
    print(paste0('Complete training ', method, ' on the full GWAS summary data with phi = ', phi))
    
    
    # Reformat output PRS file:
    n.tuning = length(phi.vals)
    out_dir = paste0(prsdir, trait_name,'.',method)
    score = NULL
    for(chr in c(1:22)){
      temfile = paste0(out_dir, '_pst_eff_a1_b0.5_phi', strsplit(PHI,split=',')[[1]][which(phi.vals == phi)], '_chr',chr,'.txt')
      if(file.exists(temfile)){
        scoretemp = bigreadr::fread2(temfile)[,c(1,2,4,5,6)]; colnames(scoretemp) = c('CHR', 'SNP', 'A1', 'A2', 'BETA')
        score = rbind(score, scoretemp)
        rm(scoretemp)
        # print(paste0('Chr ', chr,' Completed'))
      }
      if(!file.exists(temfile)) print(paste0('Need to rerun ', method, ' on CHR ',chr))
    }
    
    beta_prscs.out = score[score$BETA != 0,]
    
    ### If all effects are 0:
    r2 = bigreadr::fread2(paste0(output_path_eval,trait_name, '.', method, '.txt'))
    nonzero.ind = which(colMeans(r2) != 0)
    params.tuned.prscs = nonzero.ind[which.max(colMeans(r2)[nonzero.ind])]
    while((nrow(beta_prscs.out) == 0) & (length(nonzero.ind) > 0)){
      sumstats = sumstats0
      nonzero.ind = nonzero.ind[-which(nonzero.ind == params.tuned.prscs)]
      params.tuned.prscs = nonzero.ind[which.max(colMeans(r2)[nonzero.ind])]
      keep.SNP = sumstats.clumped[sumstats.clumped$P <= pvalthr,c('SNP')]
      dump.SNP = which(!sumstats$SNP %in% keep.SNP)
      sumstats[dump.SNP, 'BETA'] = 0
      SCORE = sumstats[,c('CHR','SNP','REF','ALT', 'BETA')]
      colnames(SCORE)[3:4] = c('A1','A2')
      # Remove SNPs with zero effects:
      beta_prscs.out = SCORE[SCORE$BETA != 0,]
    }
    if (nrow(beta_prscs.out) > 0){
      prs_prscs_outputfile = paste0(PennPRS_finalresults_path, trait_name,'.',method,'.PRS.txt')
      write_delim(beta_prscs.out, prs_prscs_outputfile)
      
      optimal.indx = params.tuned.prscs
      save(optimal.indx, file = paste0(prsdir, trait_name, '.', method, '.optimal.indx.RData'))
      rm(optimal.indx)
      write_delim(data.frame(phi = phi.vals[params.tuned.prscs]), paste0(PennPRS_finalresults_path, trait_name,'.',method,'.optimal_params.txt'))
      if ( opt$verbose >= 1 ){
        print(paste0('*************************************'))
        print(paste0('******* Complete training ', method, ' *******'))
        print(paste0('*************************************'))
      }
      # If R2<0:
      r2 = bigreadr::fread2(paste0(output_path_eval,trait_name, '.', method, '.txt'))
      if (colMeans(r2[params.tuned.prscs]) <= 0){
        err.PRSCS = 1
        print(paste0('[Warning] Trained PRS model based on ', method, ' has an estimated R2 < 0, which indicates that the PRS lacks prediction power. \nPotential explanations for R2 < 0 are:\n 1. The trait is not heritable.\n 2. The GWAS have insufficient power (e.g., due to low sample size) to develop a predictive PRS.\n 3. Issues with the input GWAS summary data (e.g., problematic BETA or SE).\n 4. ', method, ' is not powerful for developing PRS for the trait, in which case other methods can be considered.'))
      }
    }
    if (nrow(beta_prscs.out) == 0){
      err.PRSCS = 2
      print(paste0('Trained PRS model based on ', method, ' has zero effect estimate for all SNPs. Please try other methods.'))
    }
    rm(SCORE)
  }
  
  
  # Save the final PRS models generated by each single methods to the working directory
  for (method in methods){
    tfile = paste0(PennPRS_finalresults_path, trait_name,'.',method,'.PRS.txt')
    if (file.exists(tfile)) system(paste('cp -r', tfile, workdir))
  }
  
  
  filen<-paste0(workdir, 'PRS_model_training_info_', trait_name, '.txt')
  file.create(filen)
  zz <- file(filen, "w")
  
  print.title = paste0("Summary of PRS model Training on ",trait, " for ", race)
  cat(paste0("\n",paste(rep('*', nchar(print.title)+10),collapse='')), file = zz)
  cat(paste0("\n**** ", print.title, " ****"), file = zz)
  cat(paste0("\n",paste(rep('*', nchar(print.title)+10),collapse=''),'\n'), file = zz)
  
  if ('PRS-CS' %in% methods){
    cat(paste0("\n****************************************************"), file = zz)
    cat(paste0("\n********** PRS-CS (May 14, 2024 Version) **********"), file = zz)
    cat(paste0("\n****************************************************"), file = zz)
    cat(paste0("\n* Please refer to https://github.com/getian107/PRScs for detailed implementation of PRS-CS."), file = zz)
    method = 'PRS-CS'
    prsdir = paste0(prsdir0, method,'/') 
    tuning.file = paste0(PennPRS_finalresults_path, trait_name,'.',method,'.optimal_params.txt')
    if (file.exists(tuning.file)){
      prscs.out = bigreadr::fread2(tuning.file)
      cat(paste0('\n************ Tuning parameter settings: ************'), file = zz)
      cat(paste0('\nphi = ', opt$phi, '  (global shrinkage parameter)'), file = zz)
      cat(paste0('\n************** Tuned parameter values: *************'), file = zz)
      cat(paste0('\nphi = ', prscs.out[1,1], '\n'), file = zz)
    }
    if (!file.exists(tuning.file)){
      cat(paste0('\nNo PRS model was generated.\n'), file = zz)
    }
    if (err.PRSCS == 1){
      cat(paste0('[Warning] Trained PRS model has an estimated R2 < 0, indicating that the PRS lacks prediction power. \nPotential explanations for R2 < 0:\n 1. The trait is not heritable.\n 2. The GWAS have insufficient power (e.g., due to low sample size) to develop a predictive PRS.\n 3. Issues with the input GWAS summary data (e.g., problematic BETA or SE).\n 4. ', method, ' is not powerful for developing PRS for the trait.\n'), file = zz)
    }
    if (err.PRSCS == 2){
      cat(paste0('No tuning parameter setting led to a PRS model with nonzero SNP effect estimates. This is possibly due to convergence issues of the ', method, ' algorithm on the input GWAS data.\n'), file = zz)
    }
  }
}



if ('auto' %in% type){
  if ('PRS-CS' %in% methods){
    method = 'PRS-CS-auto'
    prsdir = paste0(prsdir0, method,'/')
    PRScs_path = '/dcs04/nilanjan/data/jjin/software/PRScs/'
    if ( opt$verbose >= 1 ){
      print(paste0('************************************************************'))
      print(paste0('****** Start training PRS model based on PRS-CS-auto *******'))
      print(paste0('************************************************************'))
    }
    # --------------------------------------------------------------------------------------
    # --------------------- Step 1: Train PRS models using PRS-CS-auto ---------------------
    # --------------------------------------------------------------------------------------
    # Create a separate directory 'PRS_model_training/' to store input for training PRS models
    # Within the directory 'PRS_model_training/': create a directory for each PRS method, in this case, lassosum2:
    
    # --------------------- Step 1.1: Input preparation for PRS-CS-auto ---------------------
    # First, create a "pseudo" .bim file: since we do not have individual-level validation dataset, 
    # we need to create this .bim file as one of the input files for PRS-CS-auto:
    VALIDATION_BIM_PREFIX = paste0(input_path,'pseudo_validation_',trait_name) # This is supposed to be a real .bim file for validation individuals, but we create a pseudo validation file instead 
    
    valbim = bigreadr::fread2(paste0(workdir, 'sumdata/',trait_name,'.txt')) 
    GWAS_SAMPLE_SIZE = as.character(round(mean(valbim$N)))
    valbim = cbind(valbim[,c('CHR', 'SNP')], 0, 1, valbim[,c('A1', 'A2')])
    write_delim(valbim, paste0(VALIDATION_BIM_PREFIX, '.bim'), delim = '\t', col_names = F)
    
    
    # Reformat summary data to use as the input data for PRS-CS-auto:
    gwasinput = paste0(workdir, 'sumdata/',trait_name,'.txt')
    if (!file.exists(gwasinput)) print(paste0('A valid GWAS summary data file is missing.'))
    if (file.exists(gwasinput)){
      sumraw0 = bigreadr::fread2(gwasinput)
      sumraw0 = sumraw0[,c('CHR', 'SNP', 'A1', 'A2', 'BETA', 'SE')]
      colnames(sumraw0) = c('CHR', 'SNP', 'A1', 'A2', 'BETA', 'SE') # A1: REF
      for (chr in 1:22){
        prscs.sumdat.file = paste0(prsdir,trait_name,'_reformated_gwas_chr', chr, '.txt')
        sumraw = sumraw0[sumraw0$CHR == chr, c('SNP', 'A1', 'A2', 'BETA', 'SE')]
        write_delim(sumraw, prscs.sumdat.file, delim = '\t')
      }
      print(paste0('Generating input GWAS data for ', method, ': completed.'))
    }
    
    # --------------------- Step 1.2: Run PRS-CS-auto ---------------------
    # PATH_TO_REFERENCE = paste0(PRScs_path, 'LDref/ldblk_1kg_',tolower(race)) # We can use 1000 Genomes reference data for now
    SEED = 2023
    # Can we submit 22 separate jobs to server and run them in parallel? 
    # Each of the 22 prscscode job should only require <1.5G memory.
    for (chr in 1:22) {
      out_dir = paste0(prsdir, trait_name,'.',method)
      SUM_STATS_FILE = paste0(prsdir, trait_name, '_reformated_gwas_chr', chr, '.txt')
      if (is.na(phi)){
        prscscode = paste(paste0('python ', PRScs_path, 'PRScs.py'),
                          paste0('--ref_dir=', PATH_TO_REFERENCE),
                          paste0('--bim_prefix=', VALIDATION_BIM_PREFIX),
                          paste0('--sst_file=', SUM_STATS_FILE),
                          paste0('--n_gwas=', GWAS_SAMPLE_SIZE),
                          paste0('--out_dir=', out_dir),
                          paste0('--chrom=', chr),
                          #paste0('--phi=1e+00,1e-02,1e-04,1e-6'), # this is for PRS-CS-auto grid search
                          paste0('--seed=', SEED))
      }
      if (!is.na(phi)){
        prscscode = paste(paste0('python ', PRScs_path, 'PRScs.py'),
                          paste0('--ref_dir=', PATH_TO_REFERENCE),
                          paste0('--bim_prefix=', VALIDATION_BIM_PREFIX),
                          paste0('--sst_file=', SUM_STATS_FILE),
                          paste0('--n_gwas=', GWAS_SAMPLE_SIZE),
                          paste0('--out_dir=', out_dir),
                          paste0('--chrom=', chr),
                          paste0('--phi=', phi.vals),
                          paste0('--seed=', SEED))
      }
      system(prscscode)
      print(paste0('Complete training ', method, ' for CHR ', chr))
    }
    
    
    # --------------------- Step 2.3: Reformat the trained PRS weight file: ---------------------
    score = NULL
    for(chr in c(1:22)){
      if (is.na(phi)) temfile = paste0(out_dir, '_pst_eff_a1_b0.5_phiauto_chr',chr,'.txt')
      if (!is.na(phi)) temfile = paste0(out_dir, '_pst_eff_a1_b0.5_phi', phi, '_chr',chr,'.txt')
      if(file.exists(temfile)){
        scoretemp = bigreadr::fread2(temfile)[, c(1, 2, 4, 5, 6)]
        colnames(scoretemp) = c('CHR', 'SNP', 'A1', 'A2', 'BETA')
        score = rbind(score, scoretemp)
        rm(scoretemp)
        print(paste0('Chr ', chr,' Completed'))
      }
      # if(!file.exists(temfile)) print(paste0('Need to rerun ', method, ' on CHR ',chr))
    }
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
      print(paste0('*****************************************'))
      print(paste0('***** Complete training PRS-CS-auto *****'))
      print(paste0('*****************************************'))
    }
  }
  
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
}
close(zz)
rm(zz)

# ---------------------------- README file:
filen<-paste0(workdir,'README.txt')
file.create(filen)
zz <- file(filen, "w")

print.title = paste0("List of Contents:") # ,trait, " for ", race
cat(paste0("\n",paste(rep('*', nchar(print.title)+10),collapse='')), file = zz)
cat(paste0("\n* ", print.title, " ****"), file = zz)
cat(paste0("\n",paste(rep('*', nchar(print.title)+10),collapse=''),'\n'), file = zz)
cat(paste0('\n* PRS models trained by single-ancestry methods:'), file = zz)
for (method in methods){
  prsfile = paste0(workdir, trait_name,'.',method, '.PRS.txt')
  if (file.exists(prsfile)) cat(paste0('\n  ', trait_name,'.',method, '.PRS.txt'), file = zz)
}
cat(paste0('\n\n* Details of the PRS training:'), file = zz)
cat(paste0('\n  PRS_model_training_info_', trait_name, '.txt\n'), file = zz)
close(zz)

cat(paste0('Results are saved in: ', workdir))

# Clean up intermediate files:
system(paste0('rm -rf ', paste0(workdir, 'input_for_eval/')))
system(paste0('rm -rf ', paste0(workdir, 'sumdata/')))
system(paste0('rm -rf ', paste0(workdir, 'output/')))
system(paste0('rm -rf ', paste0(workdir, 'output_for_eval/')))
system(paste0('rm -rf ', paste0(workdir, 'PRS_model_training/')))
system(paste0('rm -rf ', paste0(workdir, 'PennPRS_results/')))



