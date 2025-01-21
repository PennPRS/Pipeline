#! /usr/bin/env Rscript

########################################################################
# Deterministic Bayesian Sparse Linear Mixed Model (DBSLMM)            #
# Copyright (C) 2019  Sheng Yang and Xiang Zhou                        #
#                                                                      #
# This program is free software: you can redistribute it and/or modify #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or    #
# (at your option) any later version.                                  #
#                                                                      #
# This program is distributed in the hope that it will be useful,      #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
# GNU General Public License for more details.                         #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with this program. If not, see <http://www.gnu.org/licenses/>. #
########################################################################

library(bigreadr)
library(optparse)

## Parameter setting
args_list <- list(
  make_option("--extsumm", type = "character", default = NULL,
              help = "INPUT: perfix of summary statistics (LDSC format: SNP, N, Z, A1, A2)",
              metavar = "character"),
  make_option("--esteff", type = "character", default = NULL,
              help = "INPUT: perfix estimated effect (PLINK format: snp, effect allele, beta)", 
              metavar = "character"),
  make_option("--LDpath", type = "character", default = NULL,
              help = "INPUT: path of LD (the result from MAKELD.R)", 
              metavar = "character"), 
  make_option("--outpath", type = "character", default = NULL,
              help = "OUTPUT: log file", metavar = "character")
)

opt_parser <- OptionParser(option_list=args_list)
opt <- parse_args(opt_parser)

# output the options
cat("Acccept Options: \n")
cat("--extsumm:    ", opt$extsumm, "\n")
cat("--esteff:     ", opt$esteff, "\n")
cat("--LDpath:     ", opt$LDpath, "\n")
cat("--outpath:    ", opt$outpath, "\n")

# check the options
if (!file.exists(opt$extsumm)){
  cat(paste0("ERROR: ", opt$extsumm, " does not exist! Please check!\n"))
  q()
}
esteff_str <- unlist(strsplit(opt$esteff, ","))
if (length(esteff_str) != 4 & length(esteff_str) != 5){
  cat(paste0("ERROR: ", opt$esteff, " is not correct! Please check!\n"))
  q() 
}
if (!file.exists(esteff_str[1])){
  cat(paste0("ERROR: ", opt$esteff, " does not exist! Please check!\n"))
  q()
}

LD_str <- paste0(opt$LDpath, "chr", c(1: 22), ".RData")
if (any(file.exists(LD_str)) == F){
  cat(paste0("ERROR: ", opt$LDpath, " dose not include LD matrix for 22 chromosomes! Please check!\n"))
  q() 
}

# estimate for each chromosome
est_chr_r2 <- function(LD_path, summ_dat, est_dat, chr){
  
  require(plyr)
  comp_str <- "/net/mulan/disk2/yasheng/comparisonProject/"
  load(paste0(LD_path, "/chr", chr, ".RData"))
  r2_chr <- matrix(NA, length(LD_list), 2)
  count <- 0
  for (b in 1: length(LD_list)){
    summ_dat <- summ_dat[!is.na(summ_dat[, 3]), ]
    est_dat <- est_dat[!is.na(est_dat[, 3]), ]
    snp_inter <- Reduce(intersect, list(summ_dat[, 1], est_dat[, 1],
                                        LD_list[[b]][2][[1]]))
    if (length(snp_inter) == 0){
      
      r2_chr[b, 1] <- r2_chr[b, 2] <- NA
    } else {
      
      summ_inter <- summ_dat[match(snp_inter, summ_dat[, 1]), ]
      est_inter <- est_dat[match(snp_inter, est_dat[, 1]), ]
      summ_inter[, 3] <- ifelse(summ_inter[, 2] == est_inter[, 2], 
                                summ_inter[, 3], -summ_inter[, 3])
      
      LD_inter <- LD_list[[b]][1][[1]][match(snp_inter, LD_list[[b]][2][[1]]), 
                                       match(snp_inter, LD_list[[b]][2][[1]])]
      r2_chr[b, 1] <- as.numeric(t(summ_inter[, 3]) %*% est_inter$beta)
      r2_chr[b, 2] <- as.numeric(t(est_inter$beta) %*% LD_inter %*% est_inter$beta)
    }
    if(r2_chr[b, 1] < 0 & is.na(r2_chr[b, 1]) == FALSE){
      count <- count + 1
    }
    
  }
  return(r2_chr)
}

# process external and esteff data 
extsumm_dat <- fread2(opt$extsumm)
extsumm_std <- data.frame(SNP = extsumm_dat$SNP, 
                          A1 = extsumm_dat$A1, 
                          EFFSTD = extsumm_dat$Z / sqrt(extsumm_dat$N))
extsumm_std <- extsumm_std[!is.na(extsumm_std[, 3]), ]
esteff_dat_str <- esteff_str[1]
esteff_col_str <- as.numeric(esteff_str[-1])
esteff_dat <- fread2(esteff_dat_str)[, esteff_col_str]
if (length(esteff_col_str == 3)){
  esteff_std <- esteff_dat
  colnames(esteff_std) <- c("SNP", "A1", "beta_s")
} else {
  maf_std <- sqrt(2*(1-esteff_dat[, esteff_col_str[4]]) * esteff_dat[, esteff_col_str[4]])
  esteff_std <- data.frame(SNP = esteff_dat[, esteff_col_str[1]], 
                           A1 = esteff_dat[, esteff_col_str[2]], 
                           beta_s = esteff_dat[, esteff_col_str[3]]*maf_std)
}

## external test
ext_path <- "/home/pgsProject/Docker-pgsAlgorithm/Algorithms/package/external/"
r2_list <- list()
for (chr in 1: 22){
  r2_list[[chr]] <- est_chr_r2(LD_path = opt$LDpath, 
                               summ_dat = extsumm_std, 
                               est_dat = esteff_std, 
                               chr = chr)
  cat("chr", chr, "is finished.\n")
}
r2_dat <- ldply(r2_list, function (a) a)
r2 <- sum(r2_dat[, 1], na.rm = T)^2 / sum(r2_dat[, 2], na.rm = T)
cat("r2: ", r2, "\n")
write.table(r2, file = paste0(opt$outpath, "/r2.log"), 
            row.names = F, col.names = F, quote = F)
