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
  # cov.Y_Y.Hat[cov.Y_Y.Hat<0] = 0
  Y.Hat <- X.ref %*% weight
  var.Y.Hat <- apply(Y.Hat, 2, var)
  sum.R2 <- (unlist(cov.Y_Y.Hat))^2/(var.Y*unlist(var.Y.Hat))
  sum.R2[is.na(sum.R2)] <-0
  return(sum.R2)
}

get_omnibus_weights <- function(XtY.vtr, weights, N.vtr, X.ref){
  weights = as.data.frame(weights)
  weights.std <- apply(weights,2,function(s){return((s-mean(s))/sd(s))})
  n.nonzero = sapply(1:ncol(weights), function(x){sum(weights[,x]!=0)})
  if (sum(n.nonzero == 1)>1){
    idx.adjust = which(n.nonzero == 1)
    weights.std[,idx.adjust[-1]] = 0
  }
  # for (ii in 1:length(n.nonzero)){
  #   if (n.nonzero[ii] == 1){
  #     indx.nonzero = which(weights[,ii]!=0)
  #     weights.std[,ii] = #(weights[indx.nonzero,ii]-mean(weights[indx.nonzero,ii]))/sd(weights[indx.nonzero,ii])
  #   }
  # }
  # rm(ii)
  weights.std[is.na(weights.std)] = 0
  na.indx = which(is.na(X.ref))
  if (length(na.indx)>0) X.ref[na.indx] = 0
  Y.Hats <- X.ref %*% weights.std
  
  # calculate PRS weights
  cov.Y_Y.Hat.vtr <- t(weights.std) %*% XtY.vtr # W^T x^T  y
  # cov.Y_Y.Hat.vtr[cov.Y_Y.Hat.vtr<0] = 0
  Sigma.Y.hats <- N.vtr*cov(Y.Hats) # z^T z
  keep.indx = which(sapply(1:ncol(Sigma.Y.hats), function(x){sum(Sigma.Y.hats[,x] != 0) > 0}))
  if (length(keep.indx) == ncol(Sigma.Y.hats)){
    prs.weights <- as.numeric(solve(Sigma.Y.hats, tol = 1e-100) %*% cov.Y_Y.Hat.vtr)
  } 
  if (length(keep.indx) < ncol(Sigma.Y.hats)){
    weights.matrix = solve(Sigma.Y.hats[keep.indx, keep.indx], tol = 1e-200) %*% as.matrix(cov.Y_Y.Hat.vtr[keep.indx,],ncol = 1)
    prs.weights = rep(0, ncol(weights))
    prs.weights[keep.indx] <- as.numeric(weights.matrix)
  }
  prs.weights[prs.weights<0] <- 0
  
  # return values
  return(list(prs.weights=prs.weights))
}


get_omnibus_weights_new <- function(methods){
  r2.ent = matrix(0, k, length(methods))
  for (m in 1:length(methods)){
    method = methods[m]
    # r2.ent[,m] = as.numeric(bigreadr::fread2(paste0(output_path,trait_name,".",method,".ensembletraining.txt"))[,1])
    testing.file = paste0(output_path,trait_name,".",method,".testing.txt")
    if (file.exists(testing.file)) r2.ent[,m] = as.numeric(bigreadr::fread2(testing.file)[,1])
  }
  w.ent = suppressWarnings(log(r2.ent))
  w.ent[w.ent == '-Inf'] = 0
  w.ent[is.nan(w.ent)] = 0
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
  # prs.weights=colMeans(w.ent)
  
  # return values
  # return(list(prs.weights=prs.weights))
  return(list(prs.weights=w.ent))
}


single_prs_test <- function(single_method){
  optim.file = paste0(optimal_params_path, single_method, '/', trait_name, '.', single_method, '.optimal.indx.RData')
  if (file.exists(optim.file)){
    load(optim.file)
    r2 = bigreadr::fread2(paste0(output_path_eval,trait_name,".",single_method,".txt"))
    pumas.cor2 = data.frame(BETA = r2[,optimal.indx])
    colnames(pumas.cor2) = paste0('BETA',optimal.indx)
    # store R2
    write.table(pumas.cor2,paste0(output_path,trait_name,".",single_method,".testing.txt"),col.names = T,row.names=F,quote=F,sep="\t")
  }
  # if (!file.exists(optim.file)){
  #   cat(paste0('\nTrained PRS model based on ', single_method, ' has zero effect estimate for all SNPs and will not be incorporated in the ensemble PRS.'))
  # }
}

single_prs_ensembletraining <- function(single_method){
  stats.pumas <- as.data.frame(fread(paste0(stats_path,trait_name,".forEVAL.txt"),header=T))
  xty <- as.data.frame(fread(paste0(xty_path,trait_name,".xty.ite1.txt"),header=T))
  load(paste0(optimal_params_path, single_method, '/', trait_name, '.', single_method, '.optimal.indx.RData'))
  ## add a function here to coordinate A1/A2 between XtY and SNP weights data ##
  
  ## end ##
  
  geno.ref <- geno_ref[,match(xty$SNP,ref.geno0$SNP)]
  pumas.cor2 <- c()
  for (j in 1:k) {
    xty <- as.data.frame(fread(paste0(xty_path,trait_name,".xty.ite",j,".txt"),header=T))
    snp.w <- bigreadr::fread2(paste0(weight_path,trait_name,".",single_method,".ite",j,".txt"))
    suffix <- colnames(snp.w)[-c(1:4)][optimal.indx]
    
    snp.weight <- snp.w[, suffix]
    snp.weight[snp.weight==Inf|snp.weight==-Inf] <- 0
    snp.weight[is.na(snp.weight)] <- 0
    if (length(suffix) == 1) snp.weight = data.frame(BETA = snp.weight)
    
    # Match alleles between the reference data and the snp.w data (snp.w already matched with the summary data):
    ref.geno = ref.geno0[xty$SNP,]
    flipped = which(ref.geno$A1 != xty$A1)
    if (length(flipped) > 0){
      snp.weight[flipped,] = - snp.weight[flipped,]
      xty[flipped, 'validation_train'] = - xty[flipped, 'validation_train']
    }
    
    pumas.cor2.tmp <- get_sumR2(XtY.t=xty$validation_train, weight=snp.weight, var.Y=stats.pumas$var.Y, N.t=stats.pumas$N.vtr, X.ref=geno.ref)
    pumas.cor2 <- rbind(pumas.cor2,pumas.cor2.tmp)
    rm(ref.geno)
  }
  colnames(pumas.cor2) <- suffix
  # store R2
  write.table(pumas.cor2,paste0(output_path,trait_name,".",single_method,".ensembletraining.txt"),col.names = T,row.names=F,quote=F,sep="\t")
}




## main function for omnibus R2
omnibus_prs <- function(prs_methods){
  # first find the best tuning parameter for each PRS method
  best_param <- c()
  for (method in prs_methods){
    optim.file = paste0(optimal_params_path, method, '/', trait_name, '.', method, '.optimal.indx.RData')
    if (file.exists(optim.file)){
      load(optim.file)
      best_param <- c(best_param,  paste0('BETA',optimal.indx))
    }
    if (!(file.exists(optim.file))){
      best_param <- c(best_param,  NA)
    }
  }
  
  ## add a function here to coordinate A1/A2 between XtY and SNP weights data ##
  
  ## end ##
  
  # coordinate genotype and stats
  xty <- bigreadr::fread2(paste0(xty_path,trait_name,".xty.ite1.txt"))
  stats.pumas <- bigreadr::fread2(paste0(stats_path,trait_name,".forEVAL.txt"))
  geno.ref <- geno_ref[,match(xty$SNP,ref.geno0$SNP)]
  
  # calculate PRS weights
  pumas.weights <- pumas.weights.new <- c()
  for (j in 1:k) {
    xty <- bigreadr::fread2(paste0(xty_path,trait_name,".xty.ite",j,".txt"))
    
    # make a combined m*l matrix of m SNPs' weights for l methods
    snp.w = NULL
    if (!is.na(best_param[1])) snp.w <- bigreadr::fread2(paste0(weight_path,trait_name,".",prs_methods[1],".ite",j,".txt"),select=best_param[1])
    for (i in 2:length(prs_methods)){
      if (!is.na(best_param[i])) snp.w <- cbind(snp.w, bigreadr::fread2(paste0(weight_path,trait_name,".",prs_methods[i],".ite",j,".txt"),select=best_param[i])[,1])
      if (is.na(best_param[i])) snp.w = cbind(snp.w, 0)
    }
    if (is.na(best_param[1])) snp.w = cbind(0, snp.w)
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
  pumas.weights.new <- get_omnibus_weights_new(prs_methods)$prs.weights
  pumas.weights.avg.new <- colMeans(pumas.weights.new)
  
  colnames(pumas.weights) <- prs_methods
  
  # store results
  write.table(pumas.weights, paste0(output_path,trait_name,'.',paste0(prs_methods,collapse = '.'),".omnibus.weights.txt"),col.names = T,row.names=F,quote=F,sep="\t")
  write.table(pumas.weights.new, paste0(output_path,trait_name,'.',paste0(prs_methods,collapse = '.'),".omnibus.weights.alternative.txt"),col.names = T,row.names=F,quote=F,sep="\t")
}


## executable
cat(paste0('Start training an ensemble PRS combining the PRSs trained by ', paste0(prs_method, collapse = ', '), '.'))
# read genotype for later usage
geno_ref <- BEDMatrix(ref_path)
# geno.SNP <- fread(paste0(ref_path,".bim"),header=F)$V2

# tem = sapply(1:ncol(geno_ref), function(x){strsplit(colnames(geno_ref)[x], split='_')[[1]][1]})
# colnames(geno_ref) = tem

ref.geno0 <- bigreadr::fread2(paste0(ref_path,".bim"))
colnames(ref.geno0) = c('CHR', 'SNP', 'NA', 'POS', 'A1', 'A2')
rownames(ref.geno0) = ref.geno0$SNP


# get single PRS method r2 separately on the testing data
for (single_method in prs_method) {
  single_prs_test(single_method)
}

# get omnibus r2
omnibus_prs(prs_method)
