#! /usr/bin/env Rscript

########################################################################
# Deterministic Bayesian Sparse Linear Mixed Model (DBSLMM v1.0)       #
# Copyright (C) 2024  Sheng Yang and Xiang Zhou                        #
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

# Load parameters
library(bigsnpr)
library(bigreadr)
library(tidyverse)
library(Metrics)
require(scoring)
require(pROC)
library(optparse)

## Set parameter
args_list <- list(
  make_option("--dbslmm_eff", type = "character", default = NULL,
              help = "INPUT: ", metavar = "character"),
  make_option("--validation_g", type = "character", default = NULL,
              help = "INPUT: the predicted phenotype (PLINK output)", metavar = "character"),
  make_option("--validation_p", type = "character", default = NULL,
              help = "INPUT: the validation phenotype file and column number", metavar = "character"),
  make_option("--h2f", type = "character", default = NULL,
              help = "INPUT: the fold of heritability of trait", metavar = "character"),
  make_option("--cov", type = "character", default = NULL,
              help = "INPUT: covariate", metavar = "character")
)
opt_parser <- OptionParser(option_list=args_list)
opt <- parse_args(opt_parser)

# opt <- list(dbslmm_eff = "/public/home/biostat03/biosoft/DBSLMM/tmp_output/PGC_h2f",
#             validation_g = "/public/home/biostat03/project/pgsfusionProject/ref_geno/ref", 
#             validation_p = "/public/home/biostat03/project/pgsfusionProject/ref_geno/ref_p.txt", 
#             h2f = "1,0.7", 
#             cov = NULL)

# Check parameters
val_p <- fread2(opt$validation_p)[, 1]
val_g_fam <- fread2(paste0(opt$validation_g, ".fam"))
if(length(val_p) != nrow(val_g_fam)){
  
  cat(paste0("ERROR: Sample size of ${valdaition_g} and ${valdaition_p} does not the same! Please check!\n"))
  q()
}

na_idx <- which(!is.na(val_p))
cat(paste0(length(na_idx), " samples involves in validation datasets.\n"))

data_type <- "binary"
if(length(unique(val_p)) != 2){
  
  dat_type <- "qunatitative"
}


# Load data
h_vec <- unlist(strsplit(opt$h2f, ","))
val_dir <- dirname(opt$validation_g)

# Select parameters
index <- matrix(NA, length(h_vec), 2)
count <- 1
for (hh in h_vec){

  dbslmm_eff <- fread2(paste0(opt$dbslmm_eff, "_h2f", hh, ".dbslmm.txt"))
  colnames(dbslmm_eff) <- c("chr", "rs", "pos", "a0", "a1", "beta", 
                            "beta_ss", "large")
  if (!file.exists(paste0(opt$validation_g, ".rds")))
    snp_readBed(paste0(opt$validation_g, ".bed"))
  val_g <- snp_attach(paste0(opt$validation_g, ".rds"))
  map <- transmute(val_g$map,
                   chr = chromosome, pos = physical.pos, 
                   a0 = allele2, a1 = allele1)
  df_beta <- snp_match(dbslmm_eff, map)
  val_g_imp <- snp_fastImputeSimple(val_g$genotypes)
  pred_p <- big_prodVec(val_g_imp, df_beta$beta, ind.col = df_beta[["_NUM_ID_"]])
  
  ## Include covariables
  if(!is.null(opt$cov)){

    cov_df <- data.frame(y = val_p, cov)
    cov_model <- lm(y ~ cov, data = cov_df)
    pred_c <- predict(cov_model, cov_df)
    pred_p <- pred_p + pred_c
  } 
  
  if (dat_type == "qunatitative"){
    
    index[count, 1] <- cor(val_p[na_idx], pred_p[na_idx])^2
    index[count, 2] <- mse(val_p[na_idx], pred_p[na_idx])
  } else {
    
    index[count, 1] <- as.numeric(auc(roc(val_p[na_idx], pred_p[na_idx], 
                                          levels = c(0, 1))))
    index[count, 2] <- mean(brierscore(val_p[na_idx] ~ pred_p[na_idx]))
  }
  count <- count + 1
}
best_h2f <- h_vec[which.max(index[, 1])]

# Process file
paste0("mv ", opt$dbslmm_eff, "_h2f", best_h2f, ".dbslmm.txt ", 
       opt$dbslmm_eff, "_best.dbslmm.txt ") %>% system()
# paste0("rm -rf ", opt$dbslmm_eff, "_h2f*") %>% system()
