#!/s/bin/R35

if(!require(data.table)){
  install.packages("data.table")
  library(data.table)
}
if(!require(BEDMatrix)){
  install.packages("BEDMatrix")
  library(BEDMatrix)
}
if(!require(optparse)){
  install.packages("optparse")
  library(optparse)
}
if(!require(caret)){
  install.packages("caret")
  library(caret)
}

# Read the argument into R
options(stringsAsFactors=F)
option_list = list(
  make_option("--k", action = "store", default = NA, type = "numeric"),
  make_option("--ref_path", action = "store", default = NA, type = "character"),
  make_option("--trait_name", action = "store", default = NA, type = "character"),
  make_option("--prs_method", action = "store", default = NA, type = "character"),
  make_option("--optimal_params_path", action = "store", default = NA, type = "character"),
  make_option("--xty_path", action = "store", default = NA, type = "character"),
  make_option("--stats_path", action = "store", default = NA, type = "character"),
  make_option("--weight_path", action = "store", default = NA, type = "character"),
  make_option("--output_path_eval", action = "store", default = NA, type = "character"),
  make_option("--output_path", action = "store", default = NA, type = "character")
)

opt = parse_args(OptionParser(option_list=option_list))
k <- opt$k
ref_path <- opt$ref_path
trait_name <- opt$trait_name
prs_method <- as.character(unlist(lapply((strsplit(opt$prs_method, ',')),trimws)))
xty_path <- opt$xty_path
stats_path <- opt$stats_path
weight_path <- opt$weight_path
output_path <- opt$output_path
output_path_eval <- opt$output_path_eval
optimal_params_path <- opt$optimal_params_path


### functions

## calculate sumstats-based R2 for a single PRS method
get_sumR2 <- function(XtY.t, weight, var.Y, N.t, X.ref){
  weight <- t((t(weight) - colMeans(weight)) / apply(weight, 2, sd)) # idk why we cannot use colSds (from matrixStats)
  cov.Y_Y.Hat <- colSums(weight*(XtY.t/N.t))
  Y.Hat <- X.ref %*% weight
  var.Y.Hat <- apply(Y.Hat, 2, var)
  sum.R2 <- (unlist(cov.Y_Y.Hat))^2/(var.Y*unlist(var.Y.Hat))
  sum.R2[is.na(sum.R2)] <-0
  return(sum.R2)
}

## calculate sumstats-based omnibus R2
get_omnibus_sumR2 <- function(XtY.vt, weights, prs.weights, var.Y, N.vt, X.ref){
  
  weights.std <- apply(weights,2,function(s){return((s-mean(s))/sd(s))})
  Y.Hats <- X.ref %*% weights.std
  
  # calculate omnibus R2
  cov.Y_Y.Hat.vt <- t(prs.weights) %*% t(weights.std) %*% XtY.vt / N.vt
  var.Y.Hat <- var(as.numeric(Y.Hats %*% prs.weights))
  sum.R2 <- cov.Y_Y.Hat.vt^2/(var.Y*unlist(var.Y.Hat))
  
  # return values
  return(list(sum.R2=sum.R2))
}

get_omnibus_weights <- function(XtY.vtr, weights, N.vtr, X.ref){
  weights.std <- apply(weights,2,function(s){return((s-mean(s))/sd(s))})
  na.indx = which(is.na(X.ref))
  if (length(na.indx)>0) X.ref[na.indx] = 0
  Y.Hats <- X.ref %*% weights.std
  
  # calculate PRS weights
  cov.Y_Y.Hat.vtr <- t(weights.std) %*% XtY.vtr # W^T x^T  y
  # tmp.cor = cor(Y.Hats)
  # drop = findCorrelation(tmp.cor,cutoff = 0.995)
  # if (length(drop) > 0){
  #   weights.std = weights.std[,-drop]
  #   Y.Hats <- X.ref %*% weights.std
  #   cov.Y_Y.Hat.vtr <- t(weights.std) %*% XtY.vtr # W^T x^T  y
  # }
  Sigma.Y.hats <- N.vtr*cov(Y.Hats) # z^T z
  nan.indx = which(is.na(cov.Y_Y.Hat.vtr))
  if (length(nan.indx) == 0){
    prs.weights <- as.numeric(solve(Sigma.Y.hats, tol = 1e-80) %*% cov.Y_Y.Hat.vtr)
    prs.weights[prs.weights<0] <- 0
  }
  if (length(nan.indx) > 0){
    prs.weights = numeric(nrow(Sigma.Y.hats))
    prs.weights[-nan.indx] <- as.numeric(solve(Sigma.Y.hats[-nan.indx, -nan.indx], tol = 1e-80) %*% cov.Y_Y.Hat.vtr[-nan.indx,])
    prs.weights[prs.weights<0] <- 0
  }
  # prs.weights0 = rep(0,ncol(tmp.cor))
  # prs.weights0[-drop] = prs.weights
  # return values
  # return(list(prs.weights=prs.weights0))
  return(list(prs.weights=prs.weights))
}


get_omnibus_weights_new <- function(method){
  r2.ent = bigreadr::fread2(paste0(output_path,trait_name,".",method,".testing.txt"))
  w.ent = log( r2.ent )
  w.ent[w.ent == '-Inf'] = 0
  w.ent[is.na(w.ent)] = 0
  for (kk in 1:k){
    remove.indx = which(w.ent[kk,] == 0)
    if (length(remove.indx) > 0){
      w.ent[kk, remove.indx] = 0
      w.ent[kk, -remove.indx] = exp(w.ent[kk,-remove.indx]-suppressWarnings(max(w.ent[kk,-remove.indx])))
      w.ent[kk, -remove.indx] = w.ent[kk, -remove.indx]/sum(w.ent[kk, -remove.indx])
    }else{
      w.ent[kk,] = exp(w.ent[kk,]-max(w.ent[kk,]))
      w.ent[kk,] = w.ent[kk,]/sum(w.ent[kk,])
    }
  }
  return(list(prs.weights=w.ent))
}


## main function for omnibus R2
omnibus_prs <- function(method){
  # first find the best subset of tuning parameters for the PRS method
  # method_r2 <- as.matrix(fread(paste0(output_path,trait_name,".",method,".tuning.txt"),header=T))
  load(paste0(optimal_params_path, method,'_tuned_parameters_',trait_name,'.RData'))
  n.ensemble = length(params.tuned) # min(length(params.tuned), 15)
  params.tuned = params.tuned[1:n.ensemble]
  best_param <- paste0('BETA',params.tuned)
  
  # coordinate genotype and stats
  xty <- bigreadr::fread2(paste0(xty_path,trait_name,".xty.ite1.txt"))
  stats.pumas <- bigreadr::fread2(paste0(stats_path,trait_name,".forEVAL.txt"))
  geno.ref <- geno_ref[,match(xty$SNP,ref.geno0$SNP)]
  
  # calculate PRS weights
  pumas.weights <- c()
  pumas.weights.new <- c()
  for (j in 1:k) {
    xty <- bigreadr::fread2(paste0(xty_path,trait_name,".xty.ite",j,".txt"))
    
    # make a combined m*l matrix of m SNPs' weights for l methods
    snp.w <- bigreadr::fread2(paste0(weight_path,trait_name,".",method,".ite",j,".txt"),select=best_param)
    snp.w[snp.w==Inf|snp.w==-Inf] <- 0
    snp.w[is.na(snp.w)] <- 0

    # Match alleles between the reference data and the snp.w data (snp.w already matched with the summary data):
    ref.geno = ref.geno0[xty$SNP,]
    flipped = which(ref.geno$A1 != xty$A1)
    if (length(flipped) > 0){
      snp.w[flipped,] = - snp.w[flipped,]
      xty[flipped, 'test'] = - xty[flipped, 'test']
    }
    pumas.tmp <- get_omnibus_weights(XtY.vtr=xty$test, weights=snp.w, N.vtr=stats.pumas$N.t, X.ref=geno.ref)
    pumas.weights <- rbind(pumas.weights,pumas.tmp$prs.weights) # alpha's
    rm(ref.geno)
  }
  pumas.weights.avg <- colMeans(pumas.weights)
  pumas.weights.new <- get_omnibus_weights_new(method)$prs.weights
  pumas.weights.avg.new <- colMeans(pumas.weights.new)
  
  colnames(pumas.weights) <- paste0('BETA',params.tuned)
 
  # store results
  write.table(pumas.weights, paste0(output_path,trait_name,'.',method,".omnibus.weights.txt"),col.names = T,row.names=F,quote=F,sep="\t")
  write.table(pumas.weights.new, paste0(output_path,trait_name,'.',method,".omnibus.weights.alternative.txt"),col.names = T,row.names=F,quote=F,sep="\t")
}


## executable
# read genotype for later usage
geno_ref <- BEDMatrix(ref_path)
# geno.SNP <- fread(paste0(ref_path,".bim"),header=F)$V2

# tem = sapply(1:ncol(geno_ref), function(x){strsplit(colnames(geno_ref)[x], split='_')[[1]][1]})
# colnames(geno_ref) = tem

ref.geno0 <- bigreadr::fread2(paste0(ref_path,".bim"))
colnames(ref.geno0) = c('CHR', 'SNP', 'NA', 'POS', 'A1', 'A2')
rownames(ref.geno0) = ref.geno0$SNP

# get single PRS method r2 separately, testing
# single_prs_test(prs_method) # not necessary, already done

# get omnibus r2
omnibus_prs(prs_method)


