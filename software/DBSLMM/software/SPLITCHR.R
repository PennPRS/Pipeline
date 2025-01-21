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
  make_option("--summary", type = "character", default = NULL,
              help = "INPUT: the summary statistics (gemma output format)", 
              metavar = "character"))

opt_parser <- OptionParser(option_list=args_list)
opt <- parse_args(opt_parser)

s <- strsplit(opt$summary, "\\/")[[1]]
dir <- s[-length(s)]
dir <- paste0(dir, collapse="/")
file_name <- s[length(s)]
prefix_file <- strsplit(file_name, "\\.")[[1]]
len_prefix_file <- length(prefix_file)
prefix_file <- paste(prefix_file[-c((len_prefix_file-1):len_prefix_file)], collapse = ".")

summstats <- fread2(opt$summary)

for (chr in 1: 22){
  sum_chr <- summstats[summstats[,1]==chr, ]
  write.table(sum_chr, 
              file = paste0(dir, "/", prefix_file, "_chr", chr, ".assoc.txt"), 
              row.names = F, quote = F, col.names = F, 
              sep = "\t")
}

  