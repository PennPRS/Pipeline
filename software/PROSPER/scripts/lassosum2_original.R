
rm(list=ls())
suppressMessages(library("optparse"))
suppressMessages(library("bigreadr"))
suppressMessages(library("readr"))
suppressMessages(library("stringr"))
suppressMessages(library("caret"))

suppressMessages(library("Rcpp"))
suppressMessages(library("RcppArmadillo"))
suppressMessages(library("inline"))

suppressMessages(library("doMC"))
suppressMessages(library("foreach"))
progBar <- function(ii, N, per = 10) {
  #ii is current iteration.
  #N is total number of iterations to perform.
  #per is step (percent) when to update the progress bar. We need this as when multiple iterations are being performed progress bar gets updated too often and output is messed up.
  if (ii %in% seq(1, N, per)) {
    x <- round(ii * 100 / N)
    message("[ ",
            paste(rep("=", x), collapse = ""),
            paste(rep("-", 100 - x), collapse = ""),
            " ] ", x, "%", "\r",
            appendLF = FALSE)
    if (ii == N) cat("\r")
  }
}

option_list = list(
  make_option("--PATH_package", action="store", default=NA, type='character',
              help="Full path to the directory for the downloaded file after decompression [required]"),
  make_option("--PATH_out", action="store", default=NA, type='character',
              help="Output directory of the prs matrix and corresponding parameter settings [required]"),
  make_option("--PATH_plink", action="store", default=NA, type='character',
              help="Path to plink2 executable [required]"),

  make_option("--FILE_sst", action="store", default=NA, type='character',
              help="Full path and the file name of the GWAS summary statistics, separated by comma [required] [must have columns: rsid, chr, beta, beta_se, a1, a0; alternative columns: n_eff, and a1_af]"),
  make_option("--pop", action="store", default=NA, type='character',
              help="Population of the GWAS sample, separated by comma [required]"),
  make_option("--chrom", action="store", default="1-22", type='character',
              help="The chromosome on which the model is fitted, separated by comma or dash for consecutive chromosomes [required]"),

  make_option("--Ll", action="store", default=5, type='integer',
              help="Length of path for the tuning parameter lambda [default: %default]"),
  make_option("--Ld", action="store", default=5, type='integer',
              help="Length of path for the tuning parameter delta [default: %default]"),

  make_option("--bfile_tuning", action="store", default=NA, type='character',
              help="Path to PLINK binary input file prefix (minus bed/bim/fam) for tuning by chromosome [required]"),
  make_option("--pheno_tuning", action="store", default=NA, type='character',
              help="Path to phenotype file for tuning (PLINK format) [optional, taken from bfile otherwise]"),
  make_option("--covar_tuning", action="store", default=NA, type='character',
              help="Path to quantitative covariates for tuning (PLINK format) [optional]"),

  make_option("--testing", action="store", default=F, type='logical',
              help="Whether to perform testing in seperate dataset [required]"),
  make_option("--bfile_testing", action="store", default=NA, type='character',
              help="Path to PLINK binary input file prefix (minus bed/bim/fam) for testing by chromosome [required]"),
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

if ( opt$verbose == 2 ) { SYS_PRINT = F } else { SYS_PRINT = T }

suppressWarnings(dir.create(opt$PATH_out))

if(! dir.exists(opt$PATH_out)){
  cat( "ERROR: output path does not exist\n" )
  q()
}

ethnic_vec = str_split(opt$pop,",")[[1]]; M <- length(ethnic_vec)
sumdata_path_vec = str_split(opt$FILE_sst,",")[[1]]
out_path_vec <- paste0(opt$PATH_out,"/",ethnic_vec)

opt$chrom <- gsub("-",":",opt$chrom)
eval(parse(text=paste0("allchrom = c(",opt$chrom,")")))
Ll <- opt$Ll
Ld <- opt$Ld

bfile_tuning_vec <- str_split(opt$bfile_tuning,",")[[1]]
pheno_tuning_vec <- str_split(opt$pheno_tuning,",")[[1]]
covar_tuning_vec <- str_split(opt$covar_tuning,",")[[1]]
bfile_testing_vec <- str_split(opt$bfile_testing,",")[[1]]
pheno_testing_vec <- str_split(opt$pheno_testing,",")[[1]]
covar_testing_vec <- str_split(opt$covar_testing,",")[[1]]


# Perform i/o checks here:
files <- NULL
files <- c(files, sumdata_path_vec)
for(mmm in 1:M){ files <- c(files, paste(bfile_tuning_vec[mmm],c(".bed",".bim",".fam"),sep='')) }
suppressWarnings(if ( !is.na(pheno_tuning_vec)[1] ) { files <- c(files, pheno_tuning_vec) })
suppressWarnings(if ( !is.na(covar_tuning_vec)[1] ) { files <- c(files, covar_tuning_vec) })
if(opt$testing){
  if(is.na(bfile_testing_vec)[1]){
    cat( "ERROR: Please provide testing bfile\n" , sep='', file=stderr() )
    q()
  }
  for(mmm in 1:M){ files <- c(files, paste(bfile_testing_vec[mmm],c(".bed",".bim",".fam"),sep='')) }
  suppressWarnings(if ( !is.na(pheno_testing_vec)[1] ) { files <- c(files, pheno_testing_vec) })
  suppressWarnings(if ( !is.na(covar_testing_vec)[1] ) { files <- c(files, covar_testing_vec) })
}

for ( f in files ) {
  if ( !file.exists(f) ){
    cat( "ERROR: ", f , " input file does not exist\n" , sep='', file=stderr() )
    q()
  }
}
rm(list="files")

NCORES <- opt$NCORES

sourceCpp(paste0(opt$PATH_package,"/scripts/lassosum2.cpp"))

ref <- fread2(paste0(opt$PATH_package,"/ref_bim.txt"))

for(mmm in 1:M){

  ethnic <- ethnic_vec[mmm]
  sumdata_path <- sumdata_path_vec[mmm]
  out_path <- out_path_vec[mmm]

  suppressWarnings(dir.create(paste0(out_path)))
  suppressWarnings(dir.create(paste0(out_path, "/tmp")))
  suppressWarnings(dir.create(paste0(out_path, "/tmp/PRS_files")))
  suppressWarnings(dir.create(paste0(out_path, "/tmp/PRS_files/PRS_in_all_settings_bychrom")))

  if ( opt$verbose >= 1 ) {
    cat(paste0("\n*************************"))
    cat(paste0("\n*** LASSOSUM2 FOR ",ethnic," ***"))
    cat(paste0("\n*************************\n"))
  }

  ########################################################################
  ########################################################################
  if ( opt$verbose >= 1 ) cat("\n** Step 1. Preprocessing data **\n")

  ############
  ## Step 1.1. Loading summary data and matching to reference data

  # Match summary statistics to reference data
  df_beta <- fread2(sumdata_path)
  df_beta <- df_beta[!is.na(df_beta$beta) & !is.na(df_beta$beta_se), ]
  df_beta$n_eff[is.na(df_beta$n_eff)] <- mean(df_beta$n_eff, na.rm=T)
  df_beta <- df_beta[df_beta$rsid %in% ref$V2,]
  ref_tmp <- ref[match(df_beta$rsid, ref$V2),]

  # Match effect alleles to that in reference data
  tmp1 <- paste0(df_beta$a1, ":", df_beta$a0); tmp2 <- paste0(df_beta$a0, ":", df_beta$a1)
  tmp0 <- paste0(ref_tmp$V5, ":", ref_tmp$V6)
  flip <-  tmp2 == tmp0
  keep <- (tmp1 == tmp0 | tmp2 == tmp0)
  if ( opt$verbose == 2 ) cat( paste0(sum(flip)," SNPs are flipped and corrected \n"))

  # Adjust direction of effect for flipped variants
  if(sum(flip)>0){
    df_beta$beta[flip] <- -1 * df_beta$beta[flip]
  }
  df_beta <- df_beta[keep,,drop=F]


  # Get scale factor for standardized effect size

  N0 <- max(df_beta$n_eff, na.rm=T)
  df_beta$n_eff[is.na(df_beta$n_eff)] <- N0
  df_beta$snps_scale <- sqrt(df_beta$n_eff * df_beta$beta_se^2 + df_beta$beta^2)
  df_beta$beta_hat <- df_beta$beta / df_beta$snps_scale

  # Get maximum beta_hat to determine parameter path for lambda in the algorithm (marked by summ_max0)
  summ_max0 <- max(abs(df_beta$beta_hat))

  rm(list = c("ref_tmp","tmp0","tmp1","tmp2","flip","keep"))

  ############
  ## Step 1.2. Set parameter path

  # Parameter path for lambda
  lambdapath <- l_path(summmax = summ_max0, nlambda = Ll)
  # Parameter path for delta
  deltapath <- d_path(max=100, min=0.5, ndelta = Ld)
  # Total number of grid search
  Ngridsearch <- Ll*Ld

  ################################################################
  ################################################################

  if ( opt$verbose >= 1 ) cat("\n** Step 2. Fitting models by chromosome **\n")

  registerDoMC(NCORES)

  # Run algorithm parallelled by chromosomes
  ff <- foreach(j = 1:length(allchrom), ii = icount(), .final = function(x) NULL) %dopar% {

    chr <- allchrom[j]

    ############
    ## Step 2.1. Extract variants in both provided summary statistics and reference data

    load(paste0(opt$PATH_package,"/",ethnic,"/chr",chr,"_LD.RData"))

    # Mark overlapped variants
    m <- lapply(snps_list, FUN=function (x){x %in% df_beta$rsid})
    tmpLD <- LD_list
    tmpSNP <- snps_list
    for(i in 1:length(m)){
      if(length(tmpLD[[i]])>0){
        # Subset reference data by overlapped variants
        tmpLD[[i]] <- tmpLD[[i]][m[[i]],m[[i]],drop=F]; tmpLD[[i]][is.nan(tmpLD[[i]])] <- 1
        tmpSNP[[i]] <- tmpSNP[[i]][m[[i]]]

        # Remove variants due to homozygosity or perfect correlations
        if(nrow(tmpLD[[i]])>1){ drop = findCorrelation(tmpLD[[i]],cutoff = 0.99999) }else{ drop <- integer(0) }
        if(length(drop)>0){
          tmpLD[[i]] <- tmpLD[[i]][-drop, -drop, drop=F]
          tmpSNP[[i]] <- tmpSNP[[i]][-drop]
        }
      }
    }
    LD_list0 <- tmpLD
    snps_list0 <- tmpSNP
    Nsnps0 <- unlist(lapply(tmpSNP, length))

    # Match standardized effect size and scale factors
    tmp <- lapply(snps_list0, FUN=function (x){ df_beta[match(x, df_beta$rsid),] } )
    summ_list0 <- lapply(tmp, FUN=function (x){ x$beta_hat } )
    snps_scale0 <- lapply(tmp, FUN=function (x){ x$snps_scale } )

    if ( opt$verbose == 2 ) cat(paste0("* CHR",chr,": ", sum(unlist(m))," SNPs are included in the analysis \n"))

    rm(list = c("i","LD_list","Nsnps","snps_list","tmp","tmpLD","tmpSNP","m"))

    nblock <- length(LD_list0)

    summ_list <- summ_list0
    snps_scale <- snps_scale0
    LD_list <- LD_list0

    snp_list <- snps_list0
    Nsnps <- Nsnps0

    rm(list=c("Nsnps0","snps_list0","snps_scale0","summ_list0","LD_list0"))

    ############
    ## Step 2.2. Run algorithm

    if ( opt$verbose == 2 ) cat(paste0("Starting model fitting for CHR",chr,".. "))

    res <- enet_singlethnic(summ=summ_list, R=LD_list,
                            deltapath=deltapath, lambdapath=lambdapath,
                            verbose=opt$verbose)

    ############
    ## Step 2.3. Clean PRSs into a matrix (#variant X #grid_search)

    for (i in 1:Ngridsearch){
      for (bl in 1:nblock){
        tmp1 <- res$b[[i]][[bl]]
        tmp2 <- snps_scale[[bl]]; tmp2[is.na(tmp2)] <- 0
        if(bl==1){ b_tmp <- tmp1 * tmp2 }else{ b_tmp <- rbind(b_tmp, tmp1 * tmp2) }
      }
      if(i==1){ prs <- b_tmp }else{ prs <- cbind(prs, b_tmp) }
    }
    prs[is.na(prs)] <- 0; prs[prs > 10] <- 0; prs[prs < -10] <- 0

    ############
    ## Step 2.4. Summarize tuning parameter setting for each grid search

    param <- matrix(nrow=Ngridsearch, ncol = 3)
    param[,1] <- res$delta
    param[,2] <- res$lambda
    param[,3] <- apply(prs, MARGIN = 2, FUN = function (x){mean(x!=0)})
    colnames(param) <- c("delta", "lambda", "sparsity_nonzero_percentage")
    param <- data.frame(param)

    ############
    ## Step 2.5. Save files
    #
    # Note: 1. In the final prs file, the columns are: (rsid, a1: effect allele, a0: reference allele, PRSs...)
    #       2. For the param file, the order of its rows is same as the order of columns for PRSs. The param file indicate the tuning parameters and score source of the PRSs.

    snps <- unlist(snp_list)
    ref_tmp <- ref[match(snps, ref$V2),]
    df <- data.frame(rsid = snps, a1= ref_tmp$V5, a0= ref_tmp$V6, prs, stringsAsFactors=F)

    fwrite2(df, paste0(out_path,"/tmp/PRS_files/PRS_in_all_settings_bychrom/prs_chr",chr,".txt"), col.names = F, sep="\t", nThread=1)
    fwrite2(param, paste0(out_path,"/tmp/PRS_files/PRS_in_all_settings_bychrom/param_chr",chr,".txt"), col.names = T, sep="\t", nThread=1)

    rm(list=c("ref_tmp","b_tmp","res","tmp1","tmp2","bl","i","prs","param","df","snps",
              "nblock","summ_list","snps_scale","LD_list","snp_list","Nsnps"))

    progBar(ii, length(allchrom), per=5)

  }

  ############
  ## Step 2.6. Combine all chromosomes

  score <- foreach(j = 1:length(allchrom), .combine='rbind') %dopar% {
    chr <- allchrom[j]
    prs <- fread2(paste0(out_path,"/tmp/PRS_files/PRS_in_all_settings_bychrom/prs_chr",chr,".txt"))
    return(prs)
  }
  registerDoMC(1)
  param <- fread2(paste0(out_path,"/tmp/PRS_files/PRS_in_all_settings_bychrom/param_chr",allchrom[1],".txt")); nprs <- nrow(param)
  param[,ncol(param)] <- apply(score[,-1:-3], MARGIN = 2, FUN = function (x){mean(x!=0)})
  tmp <- apply(score[,-1:-3], MARGIN=1, function(x){sum(x!=0)}); m <- !(tmp==0)
  score <- score[m,,drop=F]
  colnames(score) <- c("rsid","a1","a0",paste0("score",1:(ncol(score)-3)))
  fwrite2(score, paste0(out_path,"/tmp/PRS_files/score_file.txt"), col.names = T, sep="\t", nThread=NCORES)
  fwrite2(param, paste0(out_path,"/tmp/PRS_files/score_param.txt"), col.names = T, sep="\t", nThread=NCORES)


  ########################################################################
  ########################################################################

  if ( opt$verbose >= 1 ) cat("\n** Step 3. Tuning **\n")

  bfile_tuning <- bfile_tuning_vec[mmm]
  pheno_tuning <- pheno_tuning_vec[mmm]
  covar_tuning <- covar_tuning_vec[mmm]
  bfile_testing <- bfile_testing_vec[mmm]
  pheno_testing <- pheno_testing_vec[mmm]
  covar_testing <- covar_testing_vec[mmm]

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
      if ( opt$verbose >= 1 ) cat( "Loaded",ncol(covar)-2,"covariates\n")
      # Match up data
      m <- match( paste(fam[,1],fam[,2]) , paste(covar[,1],covar[,2]) )
      m.keep <- !is.na(m)
      fam <- fam[m.keep,]
      pheno <- pheno[m.keep,]
      covar <- covar[m[m.keep],]
      reg <- summary(lm( pheno[,3] ~ as.matrix(covar[,3:ncol(covar)]) ))
      if ( opt$verbose == 2 ) cat( reg$r.sq , "variance in phenotype explained by covariates in tuning samples \n" )
      pheno[,3] <- scale(reg$resid)
  }

  ############
  ## Step 3.2. Calculate scores for all tuning parameter settings on tuning samples

  if ( opt$verbose == 1 ) cat("Calculating lassosum2 scores for tuning samples\n")

  suppressWarnings(dir.create(paste0(out_path, "/tmp/sample_scores")))

  arg <- paste0(opt$PATH_plink ," --threads ",NCORES,
                " --bfile ", bfile_tuning,
                " --score ", out_path,"/tmp/PRS_files/score_file.txt header-read",
                " cols=+scoresums,-scoreavgs --score-col-nums 4-",nprs+3,
                " --out ", out_path,"/tmp/sample_scores/lassosum2_tuning")
  system( arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)

  SCORE <- fread2(paste0(out_path,"/tmp/sample_scores/lassosum2_tuning.sscore"))

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

  set.seed(1)

  R2 <- numeric(length = ncol(SCORE))
  for (i in 1:ncol(SCORE)){
    fit <- lm( pheno[,3] ~ SCORE[,i] )
    R2[i] <- summary(fit)$r.square
  }
  best <- which.max(R2)

  # Save optimal tuning parameters
  delta0 <- param$delta[best]
  lambda0 <- param$lambda[best]
  tmp <- data.frame(delta0=delta0, lambda0=lambda0)
  fwrite2(tmp, paste0(out_path,"/optimal_param.txt"), col.names = T, sep="\t", nThread=NCORES)

  # Get tuning R2
  R2_res <- data.frame(tuning_R2=R2[best])
  fwrite2(R2_res, paste0(out_path,"/R2.txt"), col.names = T, sep="\t", nThread=NCORES)
  if ( opt$verbose >= 1 ) cat(paste0("lassosum2 optimal tuning parameters are saved in ", out_path,"/optimal_param.txt \n"))

  # Save lassosum2 model
  best_score <- score[,c(1,2,3,best+3)]
  colnames(best_score)[4] <- "weight"
  best_score <- best_score[best_score$weight!=0,]
  fwrite2(best_score, paste0(out_path,"/lassosum2_prs_file.txt"), col.names = T, sep="\t", nThread=NCORES)
  if ( opt$verbose >= 1 ) cat(paste0("lassosum2 model is saved in ", out_path,"/lassosum2_prs_file.txt \n"))
  if ( (opt$verbose >= 1) & !(opt$testing)) cat(paste0("** !COMPLETED! R2 is saved in ", out_path,"/R2.txt \n"))


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
        reg <- summary(lm( pheno[,3] ~ as.matrix(covar[,3:ncol(covar)]) ))
        if ( opt$verbose >= 1 ) cat( round(reg$r.sq, 4) , "variance in phenotype explained by covariates in testing samples \n" )
        pheno[,3] <- scale(reg$resid)
    }

    ############
    ## Step 4.2. Calculate scores for all tuning parameter settings on tuning samples

    arg <- paste0(opt$PATH_plink ," --threads ",NCORES,
                  " --bfile ", bfile_testing,
                  " --score ", out_path,"/lassosum2_prs_file.txt header-read",
                  " cols=+scoresums,-scoreavgs --score-col-nums 4",
                  " --out ", out_path,"/tmp/sample_scores/lassosum2_testing")
    system( arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)

    SCORE <- fread2(paste0(out_path,"/tmp/sample_scores/lassosum2_testing.sscore"))

    m <- match( paste(fam[,1],fam[,2]) , paste(SCORE[,1],SCORE[,2]) )
    m.keep <- !is.na(m)
    fam <- fam[m.keep,]
    pheno <- pheno[m.keep,]
    SCORE <- SCORE[m[m.keep],]

    # Get testing R2
    fit <- lm(pheno[,3]~SCORE[,5])
    R2 <- summary(fit)$r.square
    R2_res <- cbind(R2_res,data.frame(testing_R2=R2))
    fwrite2(R2_res, paste0(out_path,"/R2.txt"), col.names = T, sep="\t", nThread=NCORES)
    if ( opt$verbose >= 1 ) cat(paste0("** !COMPLETED! R2 is saved in ", out_path,"/R2.txt \n"))

  }

  if(opt$cleanup){
    arg = paste0("rm -rf " , out_path, "/tmp")
    system(arg)
  }

}



