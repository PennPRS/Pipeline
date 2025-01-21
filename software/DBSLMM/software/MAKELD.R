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
library(plyr)
library(optparse)

## Parameter setting
args_list <- list(
  make_option("--bfile", type = "character", default = NULL,
              help = "INPUT: the prefix of reference panel (plink file)", metavar = "character"),
  make_option("--blockfile", type = "character", default = NULL,
              help = "INPUT: the prefix of block information", 
              metavar = "character"),
  make_option("--chr", type = "character", default = NULL,
              help = "INPUT: chromosome", metavar = "character"), 
  make_option("--plink", type = "character", default = NULL,
              help = "INPUT: prefix of PLINK", metavar = "character"), 
  make_option("--outpath", type = "character", default = NULL,
              help = "OUTPUT: outpath of block LD matrix", 
              metavar = "character"))

opt_parser <- OptionParser(option_list=args_list)
opt <- parse_args(opt_parser)

## output the options
cat("Acccept Options: \n")
cat("--bfile:      ", opt$bfile, "\n")
cat("--blockfile:  ", opt$blockfile, "\n")
cat("--chr:        ", opt$chr, "\n")
cat("--plink:      ", opt$plink, "\n")
cat("--outpath:    ", opt$outpath, "\n")

## check
bfile_str <- paste0(opt$bfile, c(".bed", ".bim", ".fam"))
if(any(file.exists(bfile_str)) == F){
  cat(paste0("ERROR: ", opt$bfile, " does not exist! Please check!\n"))
  q()
}
if(!file.exists(opt$blockfile)){
  cat(paste0("ERROR: ", opt$blockfile, " does not exist! Please check!\n"))
  q()
}

## load data
chr_bim <- fread2(paste0(opt$bfile, ".bim"), select = c(2, 4, 5, 6))
blk_info <- fread2(paste0(opt$blockfile), select = c(2, 3))
blk_snp_dat <- aaply(c(1: nrow(blk_info)), 1, function(a) 
  ifelse(chr_bim[, 2]>blk_info[a, 1] & chr_bim[, 2]<=blk_info[a, 2], a, 0))

## calculate LD matrix
chr_snp <- vector()
LD_list <- list()
for (b in 1: nrow(blk_info)){
  LD_mat <- list()
  ## get LD matrix
  blk_b_str <- paste0(opt$outpath, "/chr", opt$chr, "_blk", b)
  snp_blk <- chr_bim[blk_snp_dat[b, ] == b, 1]
  if (length(snp_blk) != 0){
    write.table(snp_blk, file = paste0(blk_b_str, ".snplist"),
                row.names = F, col.names = F, quote = F)
    LD_cmd <- paste0(opt$plink, " --silent --bfile ", opt$bfile, " --r square ",
                     " --extract ", blk_b_str, ".snplist --out ", blk_b_str)
    system(LD_cmd)
    chr_snp <- c(chr_snp, snp_blk)
    cat ("blk", b, "snp number: ", length(snp_blk), "\n")
    
    LD_mat[[1]] <- as.matrix(fread2(paste0(blk_b_str, ".ld")))
    LD_mat[[2]] <- as.matrix(fread2(paste0(blk_b_str, ".snplist"), header = F))
    LD_list[[b]] <- LD_mat
    
    ## remove files
    system(paste0("rm ", blk_b_str, ".log"))
    system(paste0("rm ", blk_b_str, ".snplist"))
    system(paste0("rm ", blk_b_str, ".ld"))
    system(paste0("rm ", blk_b_str, ".nosex"))
  }
}
cat(opt$outpath, "/chr", opt$chr, ".RData\n")
save(LD_list, file = paste0(opt$outpath, "/chr", opt$chr, ".RData"))
