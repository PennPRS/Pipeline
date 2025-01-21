
rm(list=ls())
suppressMessages(library("optparse"))
suppressMessages(library("bigreadr"))
suppressMessages(library("stringr"))
suppressMessages(library("caret"))
suppressMessages(library("Rcpp"))
suppressMessages(library("RcppArmadillo"))
suppressMessages(library("inline"))

suppressMessages(library("doMC"))
suppressMessages(library("foreach"))

## progress bar function from https://stackoverflow.com/questions/51213293/is-it-possible-to-get-a-progress-bar-with-foreach-and-a-multicore-kind-of-back
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
              help="Output directory of results [required]"),

  make_option("--FILE_sst", action="store", default=NA, type='character',
              help="Full path and the file name of the GWAS summary statistics, separated by comma [required] [must have columns: rsid, chr, beta, beta_se, a1, a0, n_eff]"),
  make_option("--pop", action="store", default=NA, type='character',
              help="Population of the GWAS sample, separated by comma [required]"),
  make_option("--lassosum_param", action="store", default=NA, type='character',
              help="Directory of the optimal tuning parameters in lassosum (delta and lambda), separated by comma [required]"),
  make_option("--chrom", action="store", default="1-22", type='character',
              help="The chromosome on which the model is fitted, separated by comma or dash for consecutive chromosomes [required]"),

  make_option("--Ll", action="store", default=5, type='integer',
              help="Length of path for the tuning parameter lambda [default: %default]"),
  make_option("--Lc", action="store", default=5, type='integer',
              help="Length of path for the tuning parameter c [default: %default]"),

  make_option("--verbose", action="store", default=1, type="integer",
              help="How much chatter to print: 0=nothing; 1=minimal; 2=all [default: %default]"),
  make_option("--NCORES", action="store", default=1, type="integer",
              help="How many cores to use [default: %default]")

)
opt = parse_args(OptionParser(option_list=option_list))


NCORES <- opt$NCORES

ethnic = str_split(opt$pop,",")[[1]]; M <- length(ethnic)
sumdata_path = str_split(opt$FILE_sst,",")[[1]]
lassosum_param_path = str_split(opt$lassosum_param,",")[[1]]

# Perform i/o checks here:
for ( f in sumdata_path ) {
  if ( !file.exists(f) ){
    cat( "ERROR: ", f , " input file does not exist\n" , sep='', file=stderr() )
    q()
  }
}

opt$chrom <- gsub("-",":",opt$chrom)
eval(parse(text=paste0("allchrom = c(",opt$chrom,")")))
Ll <- opt$Ll
Lc <- opt$Lc
Ngridsearch <- Ll*Lc


suppressWarnings(dir.create(opt$PATH_out))

if(! dir.exists(opt$PATH_out)){
  cat( "ERROR: output path does not exist\n" )
  q()
}

suppressWarnings(dir.create(paste0(opt$PATH_out, "/tmp")))
suppressWarnings(dir.create(paste0(opt$PATH_out, "/tmp/PRS_in_all_settings_bychrom")))
suppressWarnings(dir.create(paste0(opt$PATH_out, "/before_ensemble")))

sourceCpp(paste0(opt$PATH_package,"/scripts/PROSPER.cpp"))


########################################################################
########################################################################

if ( opt$verbose >= 1 ) cat("\n** Step 1. Preprocessing data **\n")

############
## Step 1.1. Loading summary data and matching to reference data

ref <- fread2(paste0(opt$PATH_package,"/ref_bim.txt"))

df_beta_list <- vector("list", length = M)
summ_max <- numeric(length = M)
N0 <- integer(length = M)
for (l in 1:M){

  # Match summary statistics to reference data
  df_beta <- fread2(sumdata_path[l])
  df_beta <- df_beta[!is.na(df_beta$beta) & !is.na(df_beta$beta_se), ]
  df_beta$n_eff[is.na(df_beta$n_eff)] <- mean(df_beta$n_eff, na.rm=T)
  df_beta <- df_beta[df_beta$rsid %in% ref$V2,]
  ref_tmp <- ref[match(df_beta$rsid, ref$V2),]

  # Match effect alleles to that in reference data
  tmp1 <- paste0(df_beta$a1, ":", df_beta$a0); tmp2 <- paste0(df_beta$a0, ":", df_beta$a1)
  tmp0 <- paste0(ref_tmp$V5, ":", ref_tmp$V6)
  flip <-  tmp2 == tmp0
  keep <- (tmp1 == tmp0 | tmp2 == tmp0)
  if ( opt$verbose == 2 ){
    messageflip <- paste0(". ",sum(flip)," SNPs are flipped and corrected \n")
  }else{
    messageflip <- paste0(" \n")
  }

  # Adjust direction of effect for flipped variants
  if(sum(flip)>0){
    df_beta$beta[flip] <- -1 * df_beta$beta[flip]
  }
  df_beta <- df_beta[keep,,drop=F]


  # Get scale factor for standardized effect size
  N0[l] <- max(df_beta$n_eff, na.rm=T)
  df_beta$n_eff[is.na(df_beta$n_eff)] <- N0[l]
  df_beta$snps_scale <- sqrt(df_beta$n_eff * df_beta$beta_se^2 + df_beta$beta^2)
  df_beta$beta_hat <- df_beta$beta / df_beta$snps_scale

  # Get maximum beta_hat to determine parameter path for lambda in the algorithm (marked by summ_max)
  summ_max[l] <- max(abs(df_beta$beta_hat))
  df_beta_list[[l]] <- df_beta

  rm(list = c("df_beta","ref_tmp","tmp0","tmp1","tmp2","flip","keep"))
}

############
## Step 1.2. Load optimal tuning parameter from the single-ancestry analaysis

delta <- numeric(length = M)
lambda <- numeric(length = M)
for (l in 1:M){
  tmp <- fread2(lassosum_param_path[l])
  delta[l] <- tmp$delta0
  lambda[l] <- tmp$lambda0
  rm(list=c("tmp"))
}

############
## Step 1.3. Set parameter path

# Parameter path for lambda
lambdapath <- matrix(nrow = M, ncol = Ll)
for (l in 1:M){
  lambdapath[l,] <- r_path(rmax = min(summ_max/lambda), rmin = min(0.001/lambda), nr=Ll ) * lambda[l]
}

# Parameter path for c
cpath_tmp <- c_path(maxc=100, minc=2, nc=Lc)
cpath <- list()
for (lc in 1:Lc){
  cpath[[lc]] <- matrix(cpath_tmp[lc],ncol = M, nrow = M)
}


########################################################################
########################################################################

if ( opt$verbose >= 1 ) cat("\n** Step 2. Fitting models by chromosome **\n")

registerDoMC(NCORES)

# Run algorithm parallelled by chromosomes

ff <- foreach(j = 1:length(allchrom), ii = icount(), .final = function(x) NULL) %dopar% {

  chr <- allchrom[j]

  ############
  ## Step 2.1. Extract variants in both provided summary statistics and reference data

  Nsnps0 <- vector("list", length = M)
  snps_list0 <- vector("list", length = M)
  LD_list0 <- vector("list", length = M)
  summ_list0 <- vector("list", length = M)
  snps_scale0 <- vector("list", length = M)
  for (l in 1:M){

    load(paste0(opt$PATH_package,"/",ethnic[l],"/chr",chr,"_LD.RData"))
    df_beta <- df_beta_list[[l]]

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
    LD_list0[[l]] <- tmpLD
    snps_list0[[l]] <- tmpSNP
    Nsnps0[[l]] <- unlist(lapply(tmpSNP, length))

    # Match standardized effect size and scale factors
    tmp <- lapply(snps_list0[[l]], FUN=function (x){ df_beta[match(x, df_beta$rsid),] } )
    summ_list0[[l]] <- lapply(tmp, FUN=function (x){ x$beta_hat } )
    snps_scale0[[l]] <- lapply(tmp, FUN=function (x){ x$snps_scale } )

    if ( opt$verbose == 2 ) cat(paste0(sum(Nsnps0[[l]])," SNPs are included in the analysis of ",ethnic[l]," on CHR",chr,messageflip))

    rm(list = c("i","LD_list","Nsnps","snps_list","tmp","tmpLD","tmpSNP","m","df_beta"))
  }

  nblock <- length(LD_list0[[l]])

  ############
  ## Step 2.2. Transform to standard data format

  # Organize data in a structure fit in the algorithm
  indx_block1 <- integer(length = nblock)
  snp_list1 <- vector("list", length = nblock)
  Nsnps1 <- integer(length = nblock)
  indx1 <- vector("list", length = nblock)
  summ_list1 <- vector("list", length = nblock)
  snps_scale1 <- vector("list", length = nblock)
  LD_list1 <- vector("list", length = nblock)

  for (bl in 1:nblock){

    ## snp_list1, Nsnps1
    snp_list_tmp <- vector("list", length = M)
    tmp <- character()
    for (l in 1:M){
      if(is.null(snps_list0[[l]][[bl]])) { next }
      snp_list_tmp[[l]] <- snps_list0[[l]][[bl]]
      tmp <- c(tmp, snp_list_tmp[[l]])
    }
    tmp <- unique(tmp)
    Nsnps1[bl] <- length(tmp)
    if(Nsnps1[bl]==0){ indx_block1[bl] <- 0; next }
    snp_list1[[bl]] <- tmp
    indx_block1[bl] <- 1

    ## indx1: the position of ref SNP in original summ_list
    ## summ_list1: summ stat matched to reference snp list (set to 0 for snps not in a certain ethnic group)
    ## LD_list1: LD correlations matched to reference snp list (set to 0 for snps not in a certain ethnic group)
    indx_tmp <- matrix(nrow = Nsnps1[bl], ncol = M)
    summ_list_tmp <- matrix(nrow = Nsnps1[bl], ncol = M)
    snps_scale_tmp <- matrix(nrow = Nsnps1[bl], ncol = M)
    LD_list_tmp <- vector("list", length = M)
    for (l in 1:M){
      m <- match(snp_list1[[bl]], snp_list_tmp[[l]]); m1 <- m; m1[is.na(m1)] <- 0; indx_tmp[,l] <- m1
      m1 <- summ_list0[[l]][[bl]][m]; summ_list_tmp[,l] <- m1
      m1 <- snps_scale0[[l]][[bl]][m]; snps_scale_tmp[,l] <- m1
      m1 <- as.matrix(LD_list0[[l]][[bl]][m,m]); m1[is.na(m1)] <- 0; diag(m1)[is.na(m)] <- 1; LD_list_tmp[[l]] <- m1
    }
    indx1[[bl]] <- indx_tmp
    summ_list1[[bl]] <- summ_list_tmp
    snps_scale1[[bl]] <- snps_scale_tmp
    LD_list1[[bl]] <- LD_list_tmp

    rm(list=c("indx_tmp","summ_list_tmp","snps_scale_tmp","LD_list_tmp","snp_list_tmp","tmp","m1"))
  }

  summ_list <- summ_list1
  snps_scale <- snps_scale1
  LD_list <- LD_list1
  indx <- indx1
  indx_block <- indx_block1
  snp_list <- snp_list1
  Nsnps <- Nsnps1
  N <- N0

  rm(list=c("summ_list1","snps_scale1","LD_list1","indx1","indx_block1","snp_list1","Nsnps1",
            "Nsnps0","snps_list0","snps_scale0","summ_list0","N0","l","LD_list0","tmp"))

  ############
  ## Step 2.3. Run algorithm

  res <- enet_multiethnic(summ=summ_list, R=LD_list,
                          M=M, indx=indx,
                          indx_block=indx_block,
                          delta=delta, lambdapath=lambdapath, cpath=cpath,
                          verbose=opt$verbose)

  rm(list=c("summ_list","LD_list","Nsnps","indx","indx_block"))

  ############
  ## Step 2.4. Clean PRSs into a matrix (#variant X #grid_search)

  for (i in 1:Ngridsearch){
    for (bl in 1:nblock){
      tmp1 <- res$b[[i]][[bl]]
      tmp2 <- snps_scale[[bl]]; tmp2[is.na(tmp2)] <- 0
      if(bl==1){ b_tmp <- tmp1 * tmp2 }else{ b_tmp <- rbind(b_tmp, tmp1 * tmp2) }
    }
    if(i==1){ prs <- b_tmp }else{ prs <- cbind(prs, b_tmp) }
  }
  prs[is.na(prs)] <- 0; prs[prs > 10] <- 0; prs[prs < -10] <- 0
  rm(list=c("tmp1","tmp2","b_tmp"))

  ############
  ## Step 2.5. Summarize tuning parameter setting for each grid search

  param <- matrix(nrow=Ngridsearch*M, ncol = 2*M+2)
  for (m in 1:M){ param[,m] <- delta[m] }
  for (m in 1:M){ param[,(M+m)] <- res$lambda[m,] }
  param[,(2*M+1)] <- rep(unlist(lapply(res$c, FUN = function (x){x[1,2]})), each=M)
  param[,(2*M+2)] <- apply(prs, MARGIN = 2, FUN = function (x){mean(x!=0)})
  colnames(param) <- c(paste0("delta_",ethnic), paste0("lambda_",ethnic), "c","sparsity_nonzero_percentage")
  param <- data.frame(score_origin = rep(ethnic, Ngridsearch), param)

  ############
  ## Step 2.6. Save files
  #
  # Note: 1. In the final prs file, the columns are: (rsid, a1: effect allele, a0: reference allele, PRSs...)
  #       2. For the param file, the order of its rows is same as the order of columns for PRSs. The param file indicate the tuning parameters and score source of the PRSs.

  snps <- unlist(snp_list)
  ref_tmp <- ref[match(snps, ref$V2),]
  df <- data.frame(rsid = snps, a1= ref_tmp$V5, a0= ref_tmp$V6, prs, stringsAsFactors=F)

  fwrite2(df, paste0(opt$PATH_out,"/tmp/PRS_in_all_settings_bychrom/prs_chr",chr,".txt"), col.names = F, sep="\t", nThread=1)
  fwrite2(param, paste0(opt$PATH_out,"/tmp/PRS_in_all_settings_bychrom/param_chr",chr,".txt"), col.names = T, sep="\t", nThread=1)

  rm(list=c("snps","snp_list","res", "df","param","prs"))

  progBar(ii, length(allchrom), per=5)
}

rm(list=c("df_beta_list"))

############
## Step 2.7. Combine all chromosomes

score <- foreach(j = 1:length(allchrom), .combine='rbind') %dopar% {
  chr <- allchrom[j]
  prs <- fread2(paste0(opt$PATH_out,"/tmp/PRS_in_all_settings_bychrom/prs_chr",chr,".txt"))
  return(prs)
}
registerDoMC(1)
param <- fread2(paste0(opt$PATH_out,"/tmp/PRS_in_all_settings_bychrom/param_chr",allchrom[1],".txt"))
param[,ncol(param)] <- apply(score[,-1:-3], MARGIN = 2, FUN = function (x){mean(x!=0)})
tmp <- apply(score[,-1:-3], MARGIN=1, function(x){sum(x!=0)}); m <- !(tmp==0)
score <- score[m,,drop=F]
colnames(score) <- c("rsid","a1","a0",paste0("score",1:(ncol(score)-3)))
fwrite2(score, paste0(opt$PATH_out,"/before_ensemble/score_file.txt"), col.names = T, sep="\t", nThread=NCORES)
fwrite2(param, paste0(opt$PATH_out,"/before_ensemble/score_param.txt"), col.names = T, sep="\t", nThread=NCORES)

if ( opt$verbose >= 1 ) cat(paste0("PRSs in all tuning parameter settings are saved in ", opt$PATH_out,"/before_ensemble/score_file.txt \n"))
if ( opt$verbose >= 1 ) cat(paste0("Their corresponding tuning parameter settings are saved in ", opt$PATH_out,"/before_ensemble/score_param.txt \n"))

