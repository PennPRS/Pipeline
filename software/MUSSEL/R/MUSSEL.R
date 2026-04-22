rm(list=ls())
suppressMessages(library(RISCA))
suppressMessages(library(optparse))
suppressMessages(library(bigreadr))
suppressMessages(library(bigsnpr))
suppressMessages(library(bigparallelr))
suppressMessages(library(bigmemory))
suppressMessages(library(stringr))
suppressMessages(library(caret))
suppressMessages(library(Rcpp))
suppressMessages(library(RcppArmadillo))
suppressMessages(library(RcppTN))
suppressMessages(library(inline))
suppressMessages(library(doMC))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
suppressMessages(library(data.table))
suppressMessages(library(readr))
suppressMessages(library(MASS)) # for mvrnorm and ginv
suppressMessages(library(reshape)) # for melt
suppressMessages(library(parallel))
suppressMessages(library(devtools))
suppressMessages(library(genio)) # for read_plink
suppressMessages(library(dplyr))
# suppressMessages(library(pryr))
suppressMessages(library(Matrix))
suppressMessages(library(lavaan))
suppressMessages(library(xtable))
suppressMessages(library(SuperLearner))

option_list = list(
  make_option("--PATH_package", action="store", default=NA, type='character',
              help="Path to the directory where the downloaded files (decompressed) are saved [required]"),
  make_option("--PATH_out", action="store", default=NA, type='character',
              help="Path to the output directory where the results are saved [required]"),
  make_option("--PATH_plink", action="store", default=NA, type='character',
              help="Path to plink2 [required]"),
  make_option("--pop", action="store", default=NA, type='character',
              help="Populations of the GWAS samples, separated by comma [required]"),
  make_option("--chrom", action="store", default="1-22", type='character',
              help="The chromosome on which the model is fitted, input in the format of 1-22 or 1,2,3 [required]"),
  make_option("--SL_library", action="store", default="SL.glmnet,SL.ridge,SL.lm", type='character',
              help="The base learners implemented in SuperLearner, separated by comma [default: %default]"),
  make_option("--linear_score", action="store", default=T, type='logical',
              help="Whether the trained linear models will be saved. If not, only the Super Learner model will be saved. Note: some models in SL_library are non-linear. In this case, linear score file cannot be generated [default: %default]"),
  make_option("--target_pop", action="store", default=NA, type='character',
              help="Target population (used to save output), separated by comma [required]"),
  
  make_option("--trait_type", action="store", default='continuous', type='character',
              help="Type of phenotype, continuous or binary [required]"),
  make_option("--bfile_tuning", action="store", default=NA, type='character',
              help="Path to PLINK binary input file prefix (minus .bed/.bim/.fam) for tuning, save by chromosome [required]"),
  make_option("--pheno_tuning", action="store", default=NA, type='character',
              help="Path to phenotype file (PLINK format) for tuning, separated by comma [optional, taken from bfile otherwise]"),
  make_option("--covar_tuning", action="store", default=NA, type='character',
              help="Path to quantitative covariates (PLINK format) for tuning, separated by comma [optional]"),
  
  make_option("--testing", action="store", default=T, type='logical',
              help="Whether to perform testing in seperate dataset [required]"),
  make_option("--bfile_testing", action="store", default=NA, type='character',
              help="Path to PLINK binary input file prefix (.bed/.bim/.fam) for testing, save by chromosome [required]"),
  make_option("--pheno_testing", action="store", default=NA, type='character',
              help="Path to phenotype file (PLINK format) for testing, separated by comma [optional, taken from bfile otherwise]"),
  make_option("--covar_testing", action="store", default=NA, type='character',
              help="Path to quantitative covariates (PLINK format) for testing, separated by comma [optional]"),
  
  make_option("--verbose", action="store", default=1, type="integer",
              help="How much chatter to print: 0=nothing; 1=minimal; 2=all [default: %default]"),
  make_option("--cleanup", action="store", default=T, type="logical",
              help="Cleanup temporary files or not [default: %default]"),
  make_option("--NCORES", action="store", default=1, type="integer",
              help="How many cores to use [default: %default]")
)
opt = parse_args(OptionParser(option_list=option_list))

if ( opt$verbose == 2 ) { SYS_PRINT = F } else { SYS_PRINT = T }

NCORES <- opt$NCORES
races = str_split(opt$pop,",")[[1]]; K <- length(races)
target = str_split(opt$target_pop,",")[[1]]; K.target <- length(target)

out_paths <- paste0(opt$PATH_out,"/",races)
opt$chrom <- gsub("-",":",opt$chrom)
eval(parse(text=paste0("chrs = c(",opt$chrom,")")))
path_plink=opt$PATH_plink

trait_type <- opt$trait_type
bfile_tuning_vec <- str_split(opt$bfile_tuning,",")[[1]]
pheno_tuning_vec <- str_split(opt$pheno_tuning,",")[[1]]
covar_tuning_vec <- str_split(opt$covar_tuning,",")[[1]]
bfile_testing_vec <- str_split(opt$bfile_testing,",")[[1]]
pheno_testing_vec <- str_split(opt$pheno_testing,",")[[1]]
covar_testing_vec <- str_split(opt$covar_testing,",")[[1]]


# Perform i/o checks here:

files <- NULL
for(mmm in 1:K){ files <- c(files, paste(bfile_tuning_vec[mmm],c(".bed",".bim",".fam"),sep='')) }
suppressWarnings(if ( sum(is.na(pheno_tuning_vec)) == 0 ) { files <- c(files, pheno_tuning_vec) })
suppressWarnings(if ( sum(is.na(covar_tuning_vec)) == 0 ) { files <- c(files, covar_tuning_vec) })
if(opt$testing){
  if(is.na(bfile_testing_vec)[1]){
    cat( "ERROR: Please provide testing bfile\n" , sep='', file=stderr() )
    q()
  }
  for(mmm in 1:K){ files <- c(files, paste(bfile_testing_vec[mmm],c(".bed",".bim",".fam"),sep='')) }
  suppressWarnings(if ( sum(is.na(pheno_testing_vec)) == 0 ) { files <- c(files, pheno_testing_vec) })
  suppressWarnings(if ( sum(is.na(covar_testing_vec)) == 0 ) { files <- c(files, covar_testing_vec) })
}

for ( f in files ) {
  if ( !file.exists(f) ){
    cat( "ERROR: ", f , " input file does not exist\n" , sep='', file=stderr() )
    q()
  }
}
rm(list="files")


source(paste0(opt$PATH_package,"/R/MUSS-functions.R"))


# First, check if LDpred2 outputs were saved for all chromosomes.
outputstatus = sapply(1:K, function(x){sapply(chrs, function(y){file.exists(paste0(opt$PATH_out,'/tmp/MUSS_beta_in_all_settings_bychrom/', races[x], '-chr', y,'.txt'))})})

if (sum(outputstatus) < K * length(chrs)){
  rerun = which(rowSums(outputstatus) < 2)
  cat(paste0('\n** Terminated: need to rerun MUSS for the following chromosomes: ', paste(rerun, collapse = ','), ' first. **\n'))
  cat(paste0('\n** Rerun MUSS_jobs.R with --chrom ', paste(rerun, collapse = ','), ' **\n'))
}

settings <- bigreadr::fread2(paste0(opt$PATH_out,'/tmp/MUSS_beta_in_all_settings_bychrom/settings_',chrs[length(chrs)],".txt"));

## Combine estimated SNP effect sizes across all chromosomes
if (sum(outputstatus) == K * length(chrs)){
  
  unregister_dopar()
  
  allsnps = NULL
  betafiles = list()
  ############
  for(mmm in 1:K){
    race <- races[mmm]
    out_path <- out_paths[mmm]
    
    bfile_tuning <- bfile_tuning_vec[mmm]
    pheno_tuning <- pheno_tuning_vec[mmm]
    covar_tuning <- covar_tuning_vec[mmm]
    bfile_testing <- bfile_testing_vec[mmm]
    pheno_testing <- pheno_testing_vec[mmm]
    covar_testing <- covar_testing_vec[mmm]
    
    ncpu <- NCORES #detectCores()
    cl <- makeCluster(ncpu)
    registerDoMC(ncpu)
    ## Combine all chromosomes
    score <- foreach(j = 1:length(chrs), .combine='rbind') %dopar% {
      chr <- chrs[j]
      betas <- bigreadr::fread2(paste0(opt$PATH_out, '/tmp/MUSS_beta_in_all_settings_bychrom/', race, '-chr',chr,'.txt'))
      return(betas)
    }
    # registerDoMC(1)
    tmp <- apply(score[,-c(1:2)], MARGIN=1, function(x){sum(x!=0)}); m <- !(tmp==0)
    score <- score[m,,drop=F]
    colnames(score) <- c("rsid","a1",paste0("score",1:(ncol(score)-2)))
    allsnps = unique(c(allsnps, score$rsid))
    betafiles[[mmm]] = score
    bigreadr::fwrite2(score, paste0(opt$PATH_out,"/MUSS/beta_file_", race, ".txt"), col.names = T, sep="\t")#, nThread=NCORES)
    bigreadr::fwrite2(settings, paste0(opt$PATH_out,"/MUSS/beta_settings.txt"), col.names = T, sep="\t")#, nThread=NCORES)
  }
}

snpunion = NULL
for (k in 1:K){
  if (k == 1){
    snptemp = betafiles[[k]]$rsid
    scoretemp = data.frame(rsid = allsnps); rownames(scoretemp) = scoretemp$rsid
    scoretemp[snptemp,'a1'] = betafiles[[k]]$a1
    scoretemp[snptemp,paste0('score',1:nrow(settings), '_', races[k])] = betafiles[[k]][,paste0('score',1:nrow(settings))]
    scoretemp[is.na(scoretemp)] = 0
    scoreall = scoretemp
  }
  if (k > 1){
    snptemp = betafiles[[k]]$rsid
    scoretemp = data.frame(rsid = allsnps); rownames(scoretemp) = scoretemp$rsid
    scoretemp[snptemp,'a1'] = betafiles[[k]]$a1
    scoretemp[snptemp,paste0('score',1:nrow(settings), '_', races[k])] = betafiles[[k]][,paste0('score',1:nrow(settings))]
    scoretemp[is.na(scoretemp)] = 0
    snpint = intersect(snptemp, snpunion)
    
    matched = which((scoreall[snpint, 'a1'] == scoretemp[snpint, 'a1']))
    flipped = which((scoreall[snpint, 'a1'] != scoretemp[snpint, 'a1']))
    if (length(flipped) > 0){
      scoretemp[snpint[flipped],'a1'] = scoreall[snpint[flipped],'a1']
      scoretemp[snpint[flipped],paste0('score',1:nrow(settings), '_', races[k])] = - scoretemp[snpint[flipped],paste0('score',1:nrow(settings), '_', races[k])]
    }
    scoreall[snptemp, 'a1'] = scoretemp[snptemp, 'a1']
    scoreall[snptemp, paste0('score',1:nrow(settings), '_', races[k])] = scoretemp[snptemp, paste0('score',1:nrow(settings), '_', races[k])]
    scoreall[is.na(scoreall)] = 0
  }
  snpunion = unique(c(snpunion, snptemp))
}
bigreadr::fwrite2(scoreall, paste0(opt$PATH_out,"/MUSS/beta_file_all.txt"), col.names = T, sep="\t")#, nThread=NCORES)

if ( opt$verbose >= 1 ) cat(paste0("\n** PRSs in all tuning parameter settings are saved in ", opt$PATH_out,"/MUSS/beta_file.txt **\n"))
if ( opt$verbose >= 1 ) cat(paste0("\n** Their corresponding tuning parameter settings are saved in ", opt$PATH_out,"/MUSS/beta_settings.txt **\n"))


############

SL_library <- str_split(opt$SL_library,",")[[1]]

if(("SL.nnet" %in% SL_library) & opt$linear_score){
  question1 <- readline("nnet is non-linear and a linear score text file is unavailable. Continue? (Y/N)")
  if(regexpr(question1, 'Y', ignore.case = TRUE) == 1){
    continue <- TRUE
    opt$linear_score <- FALSE
  }else{
    q()
  }
}



########################################################################
########################################################################
settings <- fread2(paste0(opt$PATH_out,"/MUSS/beta_settings.txt"), sep="\t")#, nThread=NCORES); 
nprs <- nrow(settings) * K

for (mmm in 1:K.target){
  race <- races[mmm]
  out_path <- out_paths[mmm]
  
  bfile_tuning <- bfile_tuning_vec[mmm]
  pheno_tuning <- pheno_tuning_vec[mmm]
  covar_tuning <- covar_tuning_vec[mmm]
  bfile_testing <- bfile_testing_vec[mmm]
  pheno_testing <- pheno_testing_vec[mmm]
  covar_testing <- covar_testing_vec[mmm]
  
  # suppressWarnings(dir.create(paste0(opt$PATH_out, "/MUSS/",target[mmm])))
  
  if ( opt$verbose >= 1 ) cat("\n** Step 3. Calculate cross-ancestry PRS under all parameter settings for tuning samples **\n")
  
  ############
  ## Step 3.1.  Load data
  
  # Make/fetch the phenotype file
  fam <- read.table(paste(bfile_tuning,".fam",sep=''),as.is=T)
  if ( !is.na(pheno_tuning) ) {
    pheno <- read.table(pheno_tuning, as.is = T)
    # Match up data
    m <- match( paste(fam[,1],fam[,2]) , paste(pheno[,1],pheno[,2]) )
    m.keep <- !is.na(m)
    fam <- fam[m.keep,,drop=F]
    pheno <- pheno[m[m.keep],,drop=F]
  }else {
    pheno <- fam[,c(1,2,6)]
  }
  m <- is.na(pheno[,3]) # Remove samples with missing phenotype
  fam <- fam[!m,,drop=F]
  pheno <- pheno[!m,,drop=F]
  
  # Load in the covariates if needed
  if ( !is.na(covar_tuning) ) {
    covar <- ( read.table(covar_tuning,as.is=T,head=T) )
    if ( opt$verbose >= 1 ) cat( "Loaded",ncol(covar)-2,"covariates\n")
    # Match up data
    m <- match( paste(fam[,1],fam[,2]) , paste(covar[,1],covar[,2]) )
    m.keep <- !is.na(m)
    fam <- fam[m.keep,]
    pheno <- pheno[m.keep,]
    covar <- covar[m[m.keep],]
    if (trait_type == 'continuous'){
      reg <- summary(lm( pheno[,3] ~ as.matrix(covar[,3:ncol(covar)]) ))
      if ( opt$verbose >= 1 ) cat( reg$r.sq , "variance in phenotype explained by covariates in tuning samples \n" )
      pheno[,3] <- scale(reg$resid)
    }
  }
  
  ############
  ## Step 3.2. Calculate scores for all tuning parameter settings on tuning samples
  
  if ( opt$verbose == 2 ) cat("Calculating MUSS PRS for tuning samples\n")
  
  suppressWarnings(dir.create(paste0(opt$PATH_out,"/tmp/PRS")))
  
  
  SCORE = NULL
  
  arg <- paste0(opt$PATH_plink ," --threads 1",#NCORES,
                " --bfile ", bfile_tuning,
                " --score ", opt$PATH_out,"/MUSS/beta_file_all.txt header-read",
                " cols=+scoresums,-scoreavgs --score-col-nums 3-",nprs+2,
                " --out ",opt$PATH_out,"/tmp/PRS/MUSS_tuning_",target[mmm])
  system( arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)
  
  SCORE <- fread2(paste0(opt$PATH_out,"/tmp/PRS/MUSS_tuning_",target[mmm],".sscore"))
  
  m <- match( paste(fam[,1],fam[,2]) , paste(SCORE[,1],SCORE[,2]) )
  m.keep <- !is.na(m)
  fam <- fam[m.keep,]
  pheno <- pheno[m.keep,]
  SCORE <- SCORE[m[m.keep],]
  SCORE_id <- SCORE[,1:4]
  SCORE <- SCORE[,-1:-4]
  colnames(SCORE) <- paste0("score",1:ncol(SCORE))
  
  ############
  ## Step 3.3. Run Super Learner to obtain the ensemble PRS
  
  if ( opt$verbose >= 1 ) cat(paste0("Running Super Learner with base learners ",opt$SL_library," to obtain the final ensemble PRS \n"))
  
  # Remove constant scores (marked by score_drop)
  score_sd <- apply(SCORE,2,sd)
  score_drop <- which(is.na(score_sd) | score_sd==0)
  if(length(score_drop)>0){ SCORE <- SCORE[,-score_drop,drop=F] }
  
  set.seed(2020)
  if ( opt$verbose == 2 ){
    sl <- SuperLearner(Y = pheno[,3],
                       X = SCORE,
                       family = gaussian(),
                       SL.library = SL_library)
  }else{
    suppressWarnings(sl <- SuperLearner(Y = pheno[,3],
                                        X = SCORE,
                                        family = gaussian(),
                                        SL.library = SL_library))
  }
  suppressWarnings(dir.create(paste0(opt$PATH_out, "/MUSSEL/")))
  
  save(sl, score_drop, file = paste0(opt$PATH_out,"/MUSSEL/superlearner_function_", target[mmm], ".RData"))
  
  if ( opt$verbose >= 1 ) cat(paste0("Superlearner model of ensemble PRS saved in ", opt$PATH_out,"/superlearner_function_", target[mmm], ".RData \n"))
  
  # Predictions of ensembled scores from PROSPER on tuning samples
  after_ensemble_tuning <- cbind(pheno[,1:2], ensemble_score = predict(sl, SCORE, onlySL = TRUE)[[1]])
  fwrite2(after_ensemble_tuning, paste0(opt$PATH_out,"/tmp/PRS/after_ensemble_tuning_", target[mmm],".txt"), col.names = T, sep="\t", nThread=NCORES)
  if ( opt$verbose == 2 ) cat(paste0("MUSSEL PRS for tuning samples saved in ", opt$PATH_out,"/tmp/PRS/after_ensemble_tuning_", target[mmm],".txt \n"))
  
  # Get tuning R2
  if (trait_type == 'continuous'){
    fit <- lm(pheno[,3]~after_ensemble_tuning$ensemble_score)
    R2 <- summary(fit)$r.square
    R2_res <- data.frame(tuning_R2=R2)
    fwrite2(R2_res, paste0(opt$PATH_out,"/MUSSEL/R2_", target[mmm],".txt"), col.names = T, sep="\t", nThread=NCORES)
  }
  if (trait_type == 'binary'){
    if (!is.na(covar_tuning)){
      n.covar = ncol(covar)-2
      dat.temp = as.data.frame(cbind(pheno[,3], after_ensemble_tuning$ensemble_score, covar[,-(1:2)]))
      colnames(dat.temp) = c('y', 'prs', paste0('covar',1:n.covar))
      roc_obj = roc.binary(status = 'y',
                           variable = 'prs',
                           confounders = formula(paste0(paste('', paste(c(paste0('covar',1:n.covar)),collapse="+"), sep='~'))),
                           data = dat.temp,
                           precision=seq(0.05,0.95, by=0.05))
    }
    if (is.na(covar_tuning)){
      dat.temp = as.data.frame(cbind(pheno[,3], after_ensemble_tuning$ensemble_score))
      colnames(dat.temp) = c('y', 'prs')
      roc_obj = roc.binary(status = 'y',
                           variable = 'prs', confounders = '~1',
                           data = dat.temp,
                           precision=seq(0.05,0.95, by=0.05))
    }
    AUC <- roc_obj$auc
    rm(roc_obj)
  }
  rm(list=c("after_ensemble_tuning"))
  
  
  if(opt$linear_score){
    
    score <- fread2(paste0(opt$PATH_out,"/MUSS/beta_file_all.txt"), sep="\t", nThread=NCORES)
    if(length(score_drop)>0){ score <- score[,-score_drop,drop=F] }
    
    # Get weights of all scores (wi) by predicting on an identity matrix
    tmp <- rbind(rep(0,ncol(SCORE)), diag(ncol(SCORE)))
    colnames(tmp) <- colnames(data.frame(SCORE))
    tmp <- predict(sl, tmp, onlySL = TRUE)[[1]]
    coef <- (tmp-tmp[1])[-1]
    
    # Get weights of variants (wi * snpj)
    ensemble_score <- data.frame(score[,1:2], weight = as.matrix(score[,-(1:2)]) %*% matrix(coef, ncol=1))
    ensemble_score <- ensemble_score[ensemble_score$weight!=0,]
    
    fwrite2(ensemble_score, paste0(opt$PATH_out,"/MUSSEL/MUSSEL_beta_file_", target[mmm], ".txt"), col.names = T, sep="\t", nThread=NCORES)
    
    if ( opt$verbose >= 1 ) cat(paste0("Ensembled PROSPER model is saved in ", opt$PATH_out,"/MUSSEL/MUSSEL_prs_file_", target[mmm], ".txt \n"))
    
    rm(list=c("score","ensemble_score","coef","tmp"))
    
  }
  
  if ( (opt$verbose >= 1) & !(opt$testing)) cat(paste0("** COMPLETED! R2 on tuning samples saved in ", opt$PATH_out,"/MUSSEL/R2_", target[mmm],".txt \n"))
  
  rm(list=c("SCORE","pheno","fam"))
  
  
  ########################################################################
  ########################################################################
  
  if(opt$testing){
    
    if ( opt$verbose >= 1 ) cat("\n** Step 4. Testing **\n")
    
    ############
    ## Step 4.1. Load data
    
    # Make/fetch the phenotype file
    fam <- read.table(paste(bfile_testing,".fam",sep=''),as.is=T)
    if ( !is.na(pheno_testing) ) {
      pheno <- read.table(pheno_testing, as.is=T)
      # Match up data
      m <- match( paste(fam[,1],fam[,2]) , paste(pheno[,1],pheno[,2]) )
      m.keep <- !is.na(m)
      fam <- fam[m.keep,,drop=F]
      pheno <- pheno[m[m.keep],,drop=F]
    } else {
      pheno <- fam[,c(1,2,6)]
    }
    m <- is.na(pheno[,3]) # Remove samples with missing phenotype
    fam <- fam[!m,,drop=F]
    pheno <- pheno[!m,,drop=F]
    
    # Load in the covariates if needed
    if ( !is.na(covar_testing) ) {
      covar <- ( read.table(covar_testing,as.is=T,head=T) )
      if ( opt$verbose >= 1 ) cat( "Loaded",ncol(covar)-2,"covariates\n")
      # Match up data
      m <- match( paste(fam[,1],fam[,2]) , paste(covar[,1],covar[,2]) )
      m.keep <- !is.na(m)
      fam <- fam[m.keep,]
      pheno <- pheno[m.keep,]
      covar <- covar[m[m.keep],]
      if (trait_type == 'continuous'){
        reg <- summary(lm( pheno[,3] ~ as.matrix(covar[,3:ncol(covar)]) ))
        if ( opt$verbose >= 1 ) cat( reg$r.sq , "variance in phenotype explained by covariates in testing samples \n" )
        pheno[,3] <- scale(reg$resid)
      }
    }
    
    ############
    ## Step 4.2. Calculate scores for all tuning parameter settings on tuning samples
    
    arg <- paste0(opt$PATH_plink ," --threads 1",#NCORES,
                  " --bfile ", bfile_testing,
                  " --score ", opt$PATH_out,"/MUSS/beta_file_all.txt header-read",
                  " cols=+scoresums,-scoreavgs --score-col-nums 3-",nprs+2,
                  " --out ",opt$PATH_out,"/tmp/PRS/MUSS_testing_",target[mmm])
    system( arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)
    
    SCORE <- fread2(paste0(opt$PATH_out,"/tmp/PRS/MUSS_testing_",target[mmm],".sscore"))
    
    m <- match( paste(fam[,1],fam[,2]) , paste(SCORE[,1],SCORE[,2]) )
    m.keep <- !is.na(m)
    fam <- fam[m.keep,]
    pheno <- pheno[m.keep,]
    SCORE <- SCORE[m[m.keep],]
    SCORE_id <- SCORE[,1:4]
    SCORE <- SCORE[,-1:-4]
    colnames(SCORE) <- paste0("score",1:ncol(SCORE))
    
    if(length(score_drop)>0){ SCORE <- SCORE[,-score_drop,drop=F] }
    
    # Predictions of ensembled scores from PROSPER on testing samples
    after_ensemble_testing <- cbind(pheno[,1:2], ensemble_score = predict(sl, newdata = SCORE, onlySL = TRUE)[[1]])
    fwrite2(after_ensemble_testing, paste0(opt$PATH_out,"/tmp/PRS/after_ensemble_testing_", target[mmm],".txt"), col.names = T, sep="\t", nThread=NCORES)
    if ( opt$verbose == 2 ) cat(paste0("Predicted PROSPER scores for testing samples is saved in ", opt$PATH_out,"/tmp/sample_scores_",opt$prefix,"/after_ensemble_testing_", target[mmm],".txt \n"))
    
    # Get testing R2
    if (trait_type == 'continuous'){
      fit <- lm(pheno[,3]~after_ensemble_testing$ensemble_score)
      R2 <- summary(fit)$r.square
      R2_res <- cbind(R2_res,data.frame(testing_R2=R2))
      fwrite2(R2_res, paste0(opt$PATH_out,"/MUSSEL/R2_", target[mmm], ".txt"), col.names = T, sep="\t", nThread=NCORES)
      if ( opt$verbose >= 1 ) cat(paste0("** COMPLETED! MUSSEL R2 saved in ", opt$PATH_out,"/MUSSEL/R2_", target[mmm], ".txt \n"))
    }
    if (trait_type == 'binary'){
      if (!is.na(covar_tuning)){
        n.covar = ncol(covar)-2
        dat.temp = as.data.frame(cbind(pheno[,3], after_ensemble_testing$ensemble_score, covar[,-(1:2)]))
        colnames(dat.temp) = c('y', 'prs', paste0('covar',1:n.covar))
        roc_obj = roc.binary(status = 'y',
                             variable = 'prs',
                             confounders = formula(paste0(paste('', paste(c(paste0('covar',1:n.covar)),collapse="+"), sep='~'))),
                             data = dat.temp,
                             precision=seq(0.05,0.95, by=0.05))
      }
      if (is.na(covar_tuning)){
        dat.temp = as.data.frame(cbind(pheno[,3], after_ensemble_testing$ensemble_score))
        colnames(dat.temp) = c('y', 'prs')
        roc_obj = roc.binary(status = 'y',
                             variable = 'prs', confounders = '~1',
                             data = dat.temp,
                             precision=seq(0.05,0.95, by=0.05))
      }
      AUC <- roc_obj$auc
      AUC_res <- cbind(AUC_res,data.frame(testing_AUC=AUC))
      fwrite2(AUC_res, paste0(opt$PATH_out,"/MUSSEL/AUC_", target[mmm], ".txt"), col.names = T, sep="\t", nThread=NCORES)
      if ( opt$verbose >= 1 ) cat(paste0("** COMPLETED! MUSSEL AUC saved in ", opt$PATH_out,"/MUSSEL/AUC_", target[mmm], ".txt \n"))
    }
  }
  
  # if(opt$cleanup){
  #   arg = paste0("rm -rf " , opt$PATH_out, "/tmp")
  #   system(arg)
  # }
}



