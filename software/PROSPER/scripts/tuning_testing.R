
rm(list=ls())
suppressMessages(library("optparse"))
suppressMessages(library("bigreadr"))
suppressMessages(library("stringr"))
suppressMessages(library("SuperLearner"))
suppressMessages(library("glmnet"))
suppressMessages(library("MASS"))


option_list = list(
  make_option("--PATH_out", action="store", default=NA, type='character',
              help="Output directory of the prs matrix and corresponding parameter settings [required]"),
  make_option("--PATH_plink", action="store", default=NA, type='character',
              help="Path to plink2 executable [required]"),
  make_option("--prefix", action="store", default=NA, type='character',
              help="Prefix for target population (used to save output and results) [required]"),

  make_option("--SL_library", action="store", default="SL.glmnet,SL.ridge,SL.lm", type='character',
              help="Models input to SL.library in SuperLearner function, separated by comma [default: %default]"),
  make_option("--linear_score", action="store", default=TRUE, type='logical',
              help="Whether to output linear score text file. If not, only superlearner model is saved. Note some models available for SL.library are non-linear models. In this case, linear score file cannot be output [default: %default]"),

  make_option("--bfile_tuning", action="store", default=NA, type='character',
              help="Path to PLINK binary input file prefix (minus bed/bim/fam) for tuning [required]"),
  make_option("--pheno_tuning", action="store", default=NA, type='character',
              help="Path to phenotype file for tuning (PLINK format) [optional, taken from bfile otherwise]"),
  make_option("--covar_tuning", action="store", default=NA, type='character',
              help="Path to quantitative covariates for tuning (PLINK format) [optional]"),

  make_option("--testing", action="store", default=F, type='logical',
              help="Whether to perform testing in seperate dataset [required]"),
  make_option("--bfile_testing", action="store", default=NA, type='character',
              help="Path to PLINK binary input file prefix (minus bed/bim/fam) for testing [required]"),
  make_option("--pheno_testing", action="store", default=NA, type='character',
              help="Path to phenotype file for testing (PLINK format) [optional, taken from bfile otherwise]"),
  make_option("--covar_testing", action="store", default=NA, type='character',
              help="Path to quantitative covariates for testing (PLINK format) [optional]"),

  make_option("--verbose", action="store", default=1, type="integer",
              help="How much chatter to print: 0=nothing; 1=minimal; 2=all [default: %default]"),
  make_option("--cleanup", action="store", default=T, type="logical",
              help="Cleanup temporary files or not [default: %default]"),
  make_option("--NCORES", action="store", default=1, type="integer",
              help="How many cores to use [default: %default]")

)
opt = parse_args(OptionParser(option_list=option_list))


# Perform i/o checks here:
files <- paste(opt$bfile_tuning,c(".bed",".bim",".fam"),sep='')
if ( !is.na(opt$pheno_tuning)[1] ) files <- c(files,opt$pheno_tuning)
if ( !is.na(opt$covar_tuning)[1] ) files <- c(files,opt$covar_tuning)
if(opt$testing){
  if(is.na(opt$bfile_testing)){
    cat( "ERROR: Please provide testing bfile\n" , sep='', file=stderr() )
    q()
  }
  files <- c(files, paste(opt$bfile_testing,c(".bed",".bim",".fam"),sep=''))
  if ( !is.na(opt$pheno_testing)[1] ) files <- c(files,opt$pheno_testing)
  if ( !is.na(opt$covar_testing)[1] ) files <- c(files,opt$covar_testing)
}
for ( f in files ) {
  if ( !file.exists(f) ){
    cat( "ERROR: ", f , " input file does not exist\n" , sep='', file=stderr() )
    q()
  }
}
rm(list="files")

if ( opt$verbose == 2 ) { SYS_PRINT = F } else { SYS_PRINT = T }

NCORES <- opt$NCORES

SL_library <- str_split(opt$SL_library,",")[[1]]

if("SL.nnet" %in% SL_library & opt$linear_score){
  question1 <- readline("nnet is non-linear so linear score text file cannot be output. Continue? (Y/N)")
  if(regexpr(question1, 'Y', ignore.case = TRUE) == 1){
    continue <- TRUE
    opt$linear_score <- FALSE
  }else{
    q()
  }
}

suppressWarnings(dir.create(paste0(opt$PATH_out)))

if(! dir.exists(opt$PATH_out)){
  cat( "ERROR: output path does not exist\n" )
  q()
}

suppressWarnings(dir.create(paste0(opt$PATH_out, "/after_ensemble_",opt$prefix)))


########################################################################
########################################################################

if ( opt$verbose >= 1 ) cat("\n** Step 3. Tuning **\n")


############
## Step 3.1.  Load data

# Make/fetch the phenotype file
fam <- read.table(paste(opt$bfile_tuning,".fam",sep=''),as.is=T)
if ( !is.na(opt$pheno_tuning) ) {
  pheno <- read.table(opt$pheno_tuning, as.is = T)
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
if ( !is.na(opt$covar_tuning) ) {
  covar <- ( read.table(opt$covar_tuning,as.is=T,head=T) )
  if ( opt$verbose >= 1 ) cat( "Loaded",ncol(covar)-2,"covariates\n")
  # Match up data
  m <- match( paste(fam[,1],fam[,2]) , paste(covar[,1],covar[,2]) )
  m.keep <- !is.na(m)
  fam <- fam[m.keep,]
  pheno <- pheno[m.keep,]
  covar <- covar[m[m.keep],]
  reg <- summary(lm( pheno[,3] ~ as.matrix(covar[,3:ncol(covar)]) ))
  if ( opt$verbose >= 1 ) cat( reg$r.sq , "variance in phenotype explained by covariates in tuning samples \n" )
  pheno[,3] <- scale(reg$resid)
}

############
## Step 3.2. Calculate scores for all tuning parameter settings on tuning samples

if ( opt$verbose == 2 ) cat("Calculating lassosum2 scores for tuning samples\n")

suppressWarnings(dir.create(paste0(opt$PATH_out,"/tmp/sample_scores_",opt$prefix)))

param <- fread2(paste0(opt$PATH_out,"/before_ensemble/score_param.txt"), sep="\t", nThread=NCORES); nprs <- nrow(param)

arg <- paste0(opt$PATH_plink ," --threads ",NCORES,
              " --bfile ",opt$bfile_tuning,
              " --score ", opt$PATH_out,"/before_ensemble/score_file.txt header-read",
              " cols=+scoresums,-scoreavgs --score-col-nums 4-",nprs+3,
              " --out ",opt$PATH_out,"/tmp/sample_scores_",opt$prefix,"/before_ensemble_tuning")
system( arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)

SCORE <- fread2(paste0(opt$PATH_out,"/tmp/sample_scores_",opt$prefix,"/before_ensemble_tuning.sscore"))

m <- match( paste(fam[,1],fam[,2]) , paste(SCORE[,1],SCORE[,2]) )
m.keep <- !is.na(m)
fam <- fam[m.keep,]
pheno <- pheno[m.keep,]
SCORE <- SCORE[m[m.keep],]
SCORE_id <- SCORE[,1:4]
SCORE <- SCORE[,-1:-4]
colnames(SCORE) <- paste0("score",1:ncol(SCORE))

############
## Step 3.3. Perform super learning for ensemble PRS

if ( opt$verbose >= 1 ) cat(paste0("Performing super learning for ensemble PRS using ",opt$SL_library," \n"))

# Remove constant scores (marked by score_drop)
score_sd <- apply(SCORE,2,sd)
score_drop <- which(is.na(score_sd) | score_sd==0)
if(length(score_drop)>0){ SCORE <- SCORE[,-score_drop,drop=F] }

set.seed(1)
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
save(sl, score_drop, file = paste0(opt$PATH_out,"/after_ensemble_",opt$prefix,"/superlearner_function.RData"))

if ( opt$verbose >= 1 ) cat(paste0("Superlearner model of ensemble PRS is saved in ", opt$PATH_out,"/after_ensemble_",opt$prefix,"/superlearner_function.RData \n"))

# Predictions of ensembled scores from PROSPER on tuning samples
after_ensemble_tuning <- cbind(pheno[,1:2], ensemble_score = predict(sl, SCORE, onlySL = TRUE)[[1]])
fwrite2(after_ensemble_tuning, paste0(opt$PATH_out,"/tmp/sample_scores_",opt$prefix,"/after_ensemble_tuning.txt"), col.names = T, sep="\t", nThread=NCORES)
if ( opt$verbose == 2 ) cat(paste0("Predicted PROSPER scores for tuning samples is saved in ", opt$PATH_out,"/tmp/sample_scores_",opt$prefix,"/after_ensemble_tuning.txt \n"))

# Get tuning R2
fit <- lm(pheno[,3]~after_ensemble_tuning$ensemble_score)
R2 <- summary(fit)$r.square
R2_res <- data.frame(tuning_R2=R2)
fwrite2(R2_res, paste0(opt$PATH_out,"/after_ensemble_",opt$prefix,"/R2.txt"), col.names = T, sep="\t", nThread=NCORES)

rm(list=c("after_ensemble_tuning"))



if(opt$linear_score){

  # read coefficients of all candidate PRS models (row: variants; column: PRS models)
  score <- fread2(paste0(opt$PATH_out,"/before_ensemble/score_file.txt"), sep="\t", nThread=NCORES)
  if(length(score_drop)>0){ score <- score[,-score_drop,drop=F] }

  # Get weights of all scores (wi) by predicting on an identity matrix
  tmp <- rbind(rep(0,ncol(SCORE)), diag(ncol(SCORE)))
  colnames(tmp) <- colnames(data.frame(SCORE))
  tmp <- predict(sl, tmp, onlySL = TRUE)[[1]]
  coef <- (tmp-tmp[1])[-1]

  # Get weights of variants (wi * betaj)
  ensemble_score <- data.frame(score[,1:3], weight = as.matrix(score[,-(1:3)]) %*% matrix(coef, ncol=1))
  ensemble_score <- ensemble_score[ensemble_score$weight!=0,]

  fwrite2(ensemble_score, paste0(opt$PATH_out,"/after_ensemble_",opt$prefix,"/PROSPER_prs_file.txt"), col.names = T, sep="\t", nThread=NCORES)

  if ( opt$verbose >= 1 ) cat(paste0("Ensembled PROSPER model is saved in ", opt$PATH_out,"/after_ensemble_",opt$prefix,"/PROSPER_prs_file.txt \n"))

  rm(list=c("score","ensemble_score","coef","tmp"))

}

if ( (opt$verbose >= 1) & !(opt$testing)) cat(paste0("** !COMPLETED! R2 is saved in ", opt$PATH_out,"/after_ensemble_",opt$prefix,"/R2.txt \n"))

rm(list=c("SCORE","pheno","fam"))


########################################################################
########################################################################

if(opt$testing){

  if ( opt$verbose >= 1 ) cat("\n** Step 4. Testing **\n")

  ############
  ## Step 4.1. Load data

  # Make/fetch the phenotype file
  fam <- read.table(paste(opt$bfile_testing,".fam",sep=''),as.is=T)
  if ( !is.na(opt$pheno_testing) ) {
      pheno <- read.table(opt$pheno_testing, as.is=T)
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
  if ( !is.na(opt$covar_testing) ) {
      covar <- ( read.table(opt$covar_testing,as.is=T,head=T) )
      if ( opt$verbose >= 1 ) cat( "Loaded",ncol(covar)-2,"covariates\n")
      # Match up data
      m <- match( paste(fam[,1],fam[,2]) , paste(covar[,1],covar[,2]) )
      m.keep <- !is.na(m)
      fam <- fam[m.keep,]
      pheno <- pheno[m.keep,]
      covar <- covar[m[m.keep],]
      reg <- summary(lm( pheno[,3] ~ as.matrix(covar[,3:ncol(covar)]) ))
      if ( opt$verbose >= 1 ) cat( reg$r.sq , "variance in phenotype explained by covariates in testing samples \n" )
      pheno[,3] <- scale(reg$resid)
  }

  ############
  ## Step 4.2. Calculate scores for all tuning parameter settings on tuning samples

  arg <- paste0(opt$PATH_plink ," --threads ",NCORES,
                " --bfile ",opt$bfile_testing,
                " --score ", opt$PATH_out,"/before_ensemble/score_file.txt header-read",
                " cols=+scoresums,-scoreavgs --score-col-nums 4-",nprs+3,
                " --out ",opt$PATH_out,"/tmp/sample_scores_",opt$prefix,"/before_ensemble_testing")
  system( arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)

  SCORE <- fread2(paste0(opt$PATH_out,"/tmp/sample_scores_",opt$prefix,"/before_ensemble_testing.sscore"))

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
  after_ensemble_testing <- cbind(pheno[,1:2], ensemble_score = predict(sl, SCORE, onlySL = TRUE)[[1]])
  fwrite2(after_ensemble_testing, paste0(opt$PATH_out,"/tmp/sample_scores_",opt$prefix,"/after_ensemble_testing.txt"), col.names = T, sep="\t", nThread=NCORES)
  if ( opt$verbose == 2 ) cat(paste0("Predicted PROSPER scores for testing samples is saved in ", opt$PATH_out,"/tmp/sample_scores_",opt$prefix,"/after_ensemble_testing.txt \n"))

  # Get testing R2
  fit <- lm(pheno[,3]~after_ensemble_testing$ensemble_score)
  R2 <- summary(fit)$r.square
  R2_res <- cbind(R2_res,data.frame(testing_R2=R2))

  fwrite2(R2_res, paste0(opt$PATH_out,"/after_ensemble_",opt$prefix,"/R2.txt"), col.names = T, sep="\t", nThread=NCORES)

  if ( opt$verbose >= 1 ) cat(paste0("** !COMPLETED! R2 is saved in ", opt$PATH_out,"/after_ensemble_",opt$prefix,"/R2.txt \n"))

}

if(opt$cleanup){
  arg = paste0("rm -rf " , opt$PATH_out, "/tmp")
  system(arg)
}

