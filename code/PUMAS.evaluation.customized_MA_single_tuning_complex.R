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
  make_option("--xty_path", action = "store", default = NA, type = "character"),
  make_option("--stats_path", action = "store", default = NA, type = "character"),
  make_option("--weight_path", action = "store", default = NA, type = "character"),
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
  Y.Hats <- X.ref %*% weights.std
  
  # calculate PRS weights
  cov.Y_Y.Hat.vtr <- t(weights.std) %*% XtY.vtr
  Sigma.Y.hats <- N.vtr*cov(Y.Hats)
  prs.weights <- as.numeric(solve(Sigma.Y.hats) %*% cov.Y_Y.Hat.vtr)
  prs.weights[prs.weights<0] <- 0
  
  # return values
  return(list(prs.weights=prs.weights))
}

## main function for single method R2
single_prs <- function(prs_method){
  stats.pumas <- as.data.frame(fread(paste0(stats_path,trait_name,".omnibus.forEVAL.txt"),header=T))
  
  # snp.w <- as.data.frame(fread(paste0(weight_path,trait_name,".",prs_method,".ite1.txt"),header=T))
  # suffix <- colnames(snp.w)[-c(1:4)]
  # xty <- as.data.frame(fread(paste0(xty_path,trait_name,".xty.omnibus.ite1.txt"),header=T))
  
  ## add a function here to coordinate A1/A2 between XtY and SNP weights data ##
  
  ## end ##
  
  # geno.ref <- geno_ref[,match(xty$SNP,rs.geno)]
  pumas.cor2 <- c()
  for (j in 1:k) {
    xty <- as.data.frame(fread(paste0(xty_path,trait_name,".xty.omnibus.ite",j,".txt"),header=T))
    snp.w <- as.data.frame(fread(paste0(weight_path,trait_name,".",prs_method,".ite",j,".txt"),header=T))
    suffix <- colnames(snp.w)[-c(1:4)]
    
    # Added: find intersection of SNPs:
    snps = intersect(snp.w$SNP, xty$SNP)
    snps = intersect(snps, colnames(geno_ref0))
    geno_ref = geno_ref0[,snps]
    rownames(snp.w) = snp.w$SNP
    snp.w = snp.w[snps,]
    rownames(xty) = xty$SNP
    xty = xty[snps,]
    snp.weight <- snp.w[, suffix]
    snp.weight[snp.weight==Inf|snp.weight==-Inf] <- 0
    snp.weight[is.na(snp.weight)] <- 0
    if (length(suffix) == 1) snp.weight = data.frame(BETA = snp.weight)
    
    # Match alleles between the reference data and the snp.w data (snp.w already matched with the summary data):
    ref.geno = ref.geno0[snps,]
    flipped = which(ref.geno$A1 == xty$A1)
    if (length(flipped) > 0){
      snp.weight[flipped,] = - snp.weight[flipped,]
      xty[flipped, 'test'] = - xty[flipped, 'test']
    }
    
    pumas.cor2.tmp <- get_sumR2(XtY.t=xty$test, weight=snp.weight, var.Y=stats.pumas$var.Y, N.t=stats.pumas$N.t, X.ref=geno_ref)
    pumas.cor2 <- rbind(pumas.cor2,pumas.cor2.tmp)
    rm(ref.geno)
  }
  colnames(pumas.cor2) <- suffix
  # store R2
  write.table(pumas.cor2,paste0(output_path,trait_name,".",prs_method,".tuning.txt"),col.names = T,row.names=F,quote=F,sep="\t")
}


## executable
# read genotype for later usage
geno_ref0 <- BEDMatrix(ref_path)
# rs.geno <- fread(paste0(ref_path,".bim"),header=F)$V2
tem = sapply(1:ncol(geno_ref0), function(x){strsplit(colnames(geno_ref0)[x], split='_')[[1]][1]})
colnames(geno_ref0) = tem

ref.geno0 <- bigreadr::fread2(paste0(ref_path,".bim"))
colnames(ref.geno0) = c('CHR', 'SNP', 'NA', 'POS', 'A1', 'A2')
rownames(ref.geno0) = ref.geno0$SNP


# get single PRS method r2 separately, model tuning
for (single_method in prs_method) {
  single_prs(single_method)
}

# get single PRS method r2 separately, testing
# for (single_method in prs_method) {
#   single_prs_test(single_method)
# }

# get omnibus r2
# omnibus_prs(prs_method)



