rm(list=ls())
suppressMessages(library(RISCA))
suppressMessages(library(optparse))
suppressMessages(library(bigsnpr))
suppressMessages(library(bigreadr))
suppressMessages(library(readr))
suppressMessages(library(stringr))
suppressMessages(library(caret))
suppressMessages(library(Rcpp))
suppressMessages(library(RcppArmadillo))
suppressMessages(library(inline))
suppressMessages(library(doMC))
suppressMessages(library(foreach))
suppressMessages(library(parallel))
suppressMessages(library(doParallel))


option_list = list(
  make_option("--PATH_package", action="store", default=NA, type='character',
              help="Path to the directory where the downloaded files (decompressed) are saved [required]"),
  make_option("--PATH_out", action="store", default=NA, type='character',
              help="Path to the output directory where the results are saved [required]"),
  make_option("--PATH_plink", action="store", default=NA, type='character',
              help="Path to plink2 [required]"),
  make_option("--FILE_sst", action="store", default=NA, type='character',
              help="Paths followed by file names of the population-specific GWAS summary statistics, separated by comma [required] [required columns: chr, rsid, pos, a0, a1, beta, beta_se, n_eff]"),
  make_option("--pop", action="store", default=NA, type='character',
              help="Populations of the GWAS samples, separated by comma [required]"),
  make_option("--chrom", action="store", default="1-22", type='character',
              help="The chromosome on which the model is fitted, input in the format of 1-22 or 1,2,3 [required]"),
  
  make_option("--bfile_tuning", action="store", default=NA, type='character',
              help="Path to PLINK binary input file prefix (minus .bed/.bim/.fam) for tuning, save by chromosome [required]"),
  make_option("--pheno_tuning", action="store", default=NA, type='character',
              help="Path to phenotype file (PLINK format) for tuning, separated by comma [optional, taken from bfile otherwise]"),
  make_option("--covar_tuning", action="store", default=NA, type='character',
              help="Path to quantitative covariates (PLINK format) for tuning, separated by comma [optional]"),
  
  make_option("--trait_type", action="store", default='continuous', type='character',
              help="Type of phenotype, continuous or binary [required]"),
  make_option("--testing", action="store", default=F, type='logical',
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

suppressWarnings(dir.create(opt$PATH_out))

NCORES <- opt$NCORES
races = str_split(opt$pop,",")[[1]]; K <- length(races)
sumdata_paths = str_split(opt$FILE_sst,",")[[1]]
out_paths <- paste0(opt$PATH_out,"/",races)
path_plink <- opt$PATH_plink

opt$chrom <- gsub("-",":",opt$chrom)
eval(parse(text=paste0("chrs = c(",opt$chrom,")")))

trait_type <- opt$trait_type
bfile_tuning_vec <- str_split(opt$bfile_tuning,",")[[1]]
pheno_tuning_vec <- str_split(opt$pheno_tuning,",")[[1]]
covar_tuning_vec <- str_split(opt$covar_tuning,",")[[1]]
bfile_testing_vec <- str_split(opt$bfile_testing,",")[[1]]
pheno_testing_vec <- str_split(opt$pheno_testing,",")[[1]]
covar_testing_vec <- str_split(opt$covar_testing,",")[[1]]


# Perform i/o check:
files <- NULL
files <- c(files, sumdata_paths)
for(mmm in 1:K){ files <- c(files, paste(bfile_tuning_vec[mmm],c(".bed",".bim",".fam"),sep='')) }
n.pheno_tuning_vec = length(pheno_tuning_vec)
n.covar_tuning_vec = length(covar_tuning_vec)
suppressWarnings(if ( n.pheno_tuning_vec == K ) { files <- c(files, pheno_tuning_vec) })
suppressWarnings(if ( n.covar_tuning_vec == K ) { files <- c(files, covar_tuning_vec) })
if(opt$testing){
  if(is.na(bfile_testing_vec)[1]){
    cat( "ERROR: Please provide testing bfile\n" , sep='', file=stderr() )
    q()
  }
  for(mmm in 1:K){ files <- c(files, paste(bfile_testing_vec[mmm],c(".bed",".bim",".fam"),sep='')) }
  n.pheno_testing_vec = length(pheno_testing_vec)
  n.covar_testing_vec = length(covar_testing_vec)
  suppressWarnings(if ( n.pheno_testing_vec == K ) { files <- c(files, pheno_testing_vec) })
  suppressWarnings(if ( n.covar_testing_vec == K ) { files <- c(files, covar_testing_vec) })
  # suppressWarnings(if ( !is.na(pheno_testing_vec) ) { files <- c(files, pheno_testing_vec) })
  # suppressWarnings(if ( !is.na(covar_testing_vec) ) { files <- c(files, covar_testing_vec) })
}

for ( f in files ) {
  if ( !file.exists(f) ){
    cat( "ERROR: ", f , " input file does not exist\n" , sep='', file=stderr() )
    q()
  }
}
rm(list="files")


source(paste0(opt$PATH_package,"/R/source-functions.R"))


# First, check if LDpred2 outputs were saved for all chromosomes.

outputstatus = sapply(1:K, function(x){sapply(1:22, function(y){file.exists(paste0(out_paths[x],'/tmp/beta_files/beta_in_all_settings/ldpred2effects.txt'))})})

if (sum(outputstatus) < (K * 22)){
  rownames(outputstatus) = 1:22
  colnames(outputstatus) = races
  rerun = list()
  RERUN = character()
  for (race in races){
    rerun[[race]] = NULL
    for (chromo in 1:22) {
      if (!outputstatus[chromo, race]) rerun[[race]] = c(rerun[[race]], chromo)
    }
    if (!is.null(rerun[[race]])) RERUN[race] = paste0(race,': CHR=', paste0(rerun[[race]], collapse = ','))
  }
  rerun = which(rowSums(outputstatus) < 2)
  
  cat(paste0('\n** Terminated: need to rerun LDpred2 for: ', paste0(RERUN, collapse='; '), '. **\n'))
  cat(paste0('\n** Rerun LDpred2_jobs.R with --chrom ', paste(rerun, collapse = ','), ' **\n'))
}

if (sum(outputstatus) == (K * 22)){
  
  unregister_dopar()
  
  ############
  for(mmm in 1:K){
    race <- races[mmm]
    sumdata_path <- sumdata_paths[mmm]
    out_path <- out_paths[mmm]
    
    bfile_tuning <- bfile_tuning_vec[mmm]
    pheno_tuning <- pheno_tuning_vec[mmm]
    covar_tuning <- covar_tuning_vec[mmm]
    bfile_testing <- bfile_testing_vec[mmm]
    pheno_testing <- pheno_testing_vec[mmm]
    covar_testing <- covar_testing_vec[mmm]
    
    suppressWarnings(dir.create(paste0(out_path)))
    suppressWarnings(dir.create(paste0(out_path, "/tmp")))
    suppressWarnings(dir.create(paste0(out_path, "/tmp/beta_files")))
    suppressWarnings(dir.create(paste0(out_path, "/tmp/beta_files/beta_in_all_settings")))

    ncpu <- NCORES #detectCores()
    # cl <- makeCluster(ncpu)
    # registerDoMC(ncpu)
    ## Combine all chromosomes
    # score <- foreach(j = 1:length(chrs), .combine='rbind') %dopar% {
    #   chr <- chrs[j]
    #   betas <- bigreadr::fread2(paste0(out_path,"/tmp/beta_files/beta_in_all_settings/ldpred2effects.txt"))
    #   return(betas)
    # }
    score <- bigreadr::fread2(paste0(out_path,"/tmp/beta_files/beta_in_all_settings/ldpred2effects.txt"))
    # registerDoMC(1)
    params <- bigreadr::fread2(paste0(out_path,"/tmp/beta_files/beta_in_all_settings/params.txt")); nprs <- nrow(params)
    # params$sparsity <- apply(score[,-c(1:2)], MARGIN = 2, FUN = function (x){mean(x!=0)})
    tmp <- apply(score[,-c(1:4)], MARGIN=1, function(x){sum(x!=0)}); m <- !(tmp==0)
    score <- score[m,,drop=F]
    score <- score[,c(c('rsid', 'a0'),paste0('e',1:nrow(params)))]
    colnames(score) <- c("rsid", "a0", paste0("score",1:(ncol(score)-2)))
    bigreadr::fwrite2(score, paste0(out_path,"/tmp/beta_files/beta_file.txt"), col.names = T, sep="\t")#, nThread=NCORES)
    bigreadr::fwrite2(params, paste0(out_path,"/tmp/beta_files/beta_params.txt"), col.names = T, sep="\t")#, nThread=NCORES)
    
    
    ########################################################################
    ########################################################################
    
    if ( opt$verbose >= 1 ) cat("\n** Step 3: LDpred2 tuning **\n")
    
    ############
    ## Step 3.1. Load data
    
    # Make/fetch the phenotype file
    fam <- read.table(paste(bfile_tuning,".fam",sep=''),as.is=T)
    if ( !is.na(pheno_tuning) ) {
      pheno <- read.table(pheno_tuning, as.is=T)
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
    if ( !is.na(covar_tuning) ) {
      covar <- read.table(covar_tuning,as.is=T,head=T)
      if ( opt$verbose >= 1 ) cat(ncol(covar)-2,"covariates loaded. \n")
      # Match up data
      m <- match( paste(fam[,1],fam[,2]) , paste(covar[,1],covar[,2]) )
      m.keep <- !is.na(m)
      fam <- fam[m.keep,]
      pheno <- pheno[m.keep,]
      covar <- covar[m[m.keep],]
      if (trait_type == 'continuous'){
        reg <- summary(lm( pheno[,3] ~ as.matrix(covar[,3:ncol(covar)]) ))
        if ( opt$verbose == 2 ) cat( reg$r.sq , "variance in phenotype explained by covariates in tuning samples \n" )
        pheno[,3] <- scale(reg$resid)
      }
    }
    
    ############
    ## Step 3.2. Calculate PRS under all tuning parameter settings on the tuning samples
    
    if ( opt$verbose == 1 ) cat("Calculating LDpred2 PRS for tuning samples\n")
    
    suppressWarnings(dir.create(paste0(out_path, "/tmp/PRS")))
    
    arg <- paste0(opt$PATH_plink," --threads 1",#NCORES,
                  " --bfile ", bfile_tuning,
                  " --score ", out_path,"/tmp/beta_files/beta_file.txt header-read",
                  " cols=+scoresums,-scoreavgs --score-col-nums 3-",nprs+2,
                  " --out ", out_path,"/tmp/PRS/LDpred2_PRS_tuning")
    system(arg, ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)
    
    SCORE <- bigreadr::fread2(paste0(out_path,"/tmp/PRS/LDpred2_PRS_tuning.sscore"))
    
    m <- match( paste(fam[,1],fam[,2]) , paste(SCORE[,1],SCORE[,2]) )
    m.keep <- !is.na(m)
    fam <- fam[m.keep,]
    pheno <- pheno[m.keep,]
    SCORE <- SCORE[m[m.keep],]
    SCORE_id <- SCORE[,1:4]
    SCORE <- SCORE[,-1:-4]
    colnames(SCORE) <- paste0("score",1:ncol(SCORE))
    
    ############
    ## Step 3.3. Find optimal tuning parameters
    
    suppressWarnings(dir.create(paste0(opt$PATH_out, "/LDpred2/")))
    set.seed(1)
    
    if (trait_type == 'continuous'){
      R2 <- numeric(length = ncol(SCORE))
      for (i in 1:ncol(SCORE)){
        fit <- lm( pheno[,3] ~ SCORE[,i] )
        R2[i] <- summary(fit)$r.square
      }
      indx <- which.max(R2)
    }
    if (trait_type == 'binary'){
      AUC <- numeric(length = ncol(SCORE))
      for (i in 1:ncol(SCORE)){
        if (!is.na(covar_tuning)){
          n.covar = ncol(covar)-2
          dat.temp = as.data.frame(cbind(pheno[,3], SCORE[,i], covar[,-(1:2)]))
          colnames(dat.temp) = c('y', 'prs', paste0('covar',1:n.covar))
          roc_obj = roc.binary(status = 'y',
                               variable = 'prs',
                               confounders = formula(paste0(paste('', paste(c(paste0('covar',1:n.covar)),collapse="+"), sep='~'))),
                               data = dat.temp,
                               precision=seq(0.05,0.95, by=0.05))
        }
        if (is.na(covar_tuning)){
          dat.temp = as.data.frame(cbind(pheno[,3], SCORE[,i]))
          colnames(dat.temp) = c('y', 'prs')
          roc_obj = roc.binary(status = 'y',
                               variable = 'prs', confounders = '~1',
                               data = dat.temp,
                               precision=seq(0.05,0.95, by=0.05))
        }
        AUC[i] <- roc_obj$auc
        rm(roc_obj)
      }
      indx <- which.max(AUC)
    }
    
    # Save optimal tuning parameters
    p0 <- params$p[indx]
    h20 <- params$h2[indx]
    sparse0 <- params$sparse[indx]
    optim_params <- data.frame(p0 = p0, h20 = h20, sparse0 = sparse0)
    bigreadr::fwrite2(optim_params, paste0(opt$PATH_out, '/LDpred2/', race, "_optim_params.txt"), col.names = T, sep="\t")#, nThread=NCORES)
    
    # Get tuning R2/AUC
    if (trait_type == 'continuous'){
      R2_res <- data.frame(tuning_R2=R2[indx])
      bigreadr::fwrite2(R2_res, paste0(opt$PATH_out, '/LDpred2/', race, "_LDpred2_best_R2_tuning.txt"), col.names = T, sep="\t")#, nThread=NCORES)
      if ( opt$verbose >= 1 ) cat(paste0(race," LDpred2 optimal tuning parameters saved in ", opt$PATH_out, '/LDpred2/', race, "_optimal_param.txt \n"))
    }
    if (trait_type == 'binary'){
      AUC_res <- data.frame(tuning_AUC=AUC[indx])
      bigreadr::fwrite2(AUC_res, paste0(opt$PATH_out, '/LDpred2/', race, "_LDpred2_best_AUC_tuning.txt"), col.names = T, sep="\t")#, nThread=NCORES)
      if ( opt$verbose >= 1 ) cat(paste0(race," LDpred2 optimal tuning parameters saved in ", opt$PATH_out, '/LDpred2/', race, "_optimal_param.txt \n"))
    }
    
    # Save estimated SNP effect sizes from LDpred2
    best_score <- score[,c(1,2,indx+2)]
    colnames(best_score)[3] <- "weight"
    best_score <- best_score[best_score$weight!=0,]
    bigreadr::fwrite2(best_score, paste0(opt$PATH_out, '/LDpred2/', race, "_LDpred2_beta.txt"), col.names = T, sep="\t")#, nThread=NCORES)
    if ( opt$verbose >= 1 ) cat(paste0("LDpred2 model is saved in ", opt$PATH_out, '/LDpred2/', race, "_LDpred2_beta.txt \n"))
    if ( (opt$verbose >= 1) & !(opt$testing)) cat(paste0("***** Completed! R2 on tuning sample is saved in ", opt$PATH_out, '/LDpred2/', race, "_LDpred2_best_R2_tuning.txt *****\n"))
    
    
    ########################################################################
    ########################################################################
    
    if(opt$testing){
      
      if ( opt$verbose >= 1 ) cat("\n** Step 4: LDpred2 testing **\n")
      
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
          if ( opt$verbose >= 1 ) cat( round(reg$r.sq, 4) , "variance in phenotype explained by covariates in testing samples \n" )
          pheno[,3] <- scale(reg$resid)
        }
      }
      
      ############
      ## Step 4.2. Calculate scores for all tuning parameter settings on tuning samples
      
      arg <- paste0(opt$PATH_plink," --threads 1 --silent",#NCORES,
                    " --bfile ", bfile_testing,
                    " --score ", opt$PATH_out, '/LDpred2/', race, "_LDpred2_beta.txt header-read",
                    " cols=+scoresums,-scoreavgs --score-col-nums 3",
                    " --out ", out_path,"/tmp/PRS/LDpred2_PRS_testing")
      
      system( arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)
      
      SCORE <- bigreadr::fread2(paste0(out_path,"/tmp/PRS/LDpred2_PRS_testing.sscore"))
      
      m <- match( paste(fam[,1],fam[,2]) , paste(SCORE[,1],SCORE[,2]) )
      m.keep <- !is.na(m)
      fam <- fam[m.keep,]
      pheno <- pheno[m.keep,]
      SCORE <- SCORE[m[m.keep],]
      
      # Get testing R2
      if (trait_type == 'continuous'){
        fit <- lm(pheno[,3]~SCORE[,5])
        R2 <- summary(fit)$r.square
        R2_res <- cbind(R2_res,data.frame(testing_R2=R2))
        bigreadr::fwrite2(R2_res, paste0(opt$PATH_out, '/LDpred2/', race, "_LDpred_R2_testing.txt"), col.names = T, sep="\t")#, nThread=NCORES)
        if ( opt$verbose >= 1 ) cat(paste0("** !COMPLETED! R2 is saved in ", opt$PATH_out, '/LDpred2/', race, "_LDpred_R2_testing.txt \n"))
      }
      if (trait_type == 'binary'){
        if (!is.na(covar_tuning)){
          n.covar = ncol(covar)-2
          dat.temp = as.data.frame(cbind(pheno[,3], SCORE[,5], covar[,-(1:2)]))
          colnames(dat.temp) = c('y', 'prs', paste0('covar',1:n.covar))
          roc_obj = roc.binary(status = 'y',
                               variable = 'prs',
                               confounders = formula(paste0(paste('', paste(c(paste0('covar',1:n.covar)),collapse="+"), sep='~'))),
                               data = dat.temp,
                               precision=seq(0.05,0.95, by=0.05))
        }
        if (is.na(covar_tuning)){
          dat.temp = as.data.frame(cbind(pheno[,3], SCORE[,5]))
          colnames(dat.temp) = c('y', 'prs')
          roc_obj = roc.binary(status = 'y',
                               variable = 'prs', confounders = '~1',
                               data = dat.temp,
                               precision=seq(0.05,0.95, by=0.05))
        }
        AUC <- roc_obj$auc
        AUC_res <- cbind(AUC_res,data.frame(testing_AUC=AUC))
        bigreadr::fwrite2(AUC_res, paste0(opt$PATH_out, '/LDpred2/', race, "_LDpred_AUC_testing.txt"), col.names = T, sep="\t")#, nThread=NCORES)
        if ( opt$verbose >= 1 ) cat(paste0("** !COMPLETED! AUC is saved in ", opt$PATH_out, '/LDpred2/', race, "_LDpred_AUC_testing.txt \n"))
      }
    }
    
    # if(opt$cleanup){
    #   arg = paste0("rm -rf " , out_path, "/tmp")
    #   system(arg)
    # }
  }
}



