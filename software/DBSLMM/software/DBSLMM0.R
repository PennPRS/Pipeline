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

# Load packages
library(bigsnpr)
library(bigreadr)
library(dplyr)
library(plyr)
library(optparse)

# Set parameters
args_list <- list(
  ## Required parameters
  make_option("--summary", type = "character", default = NULL,
              help = "INPUT: the summary statistics (gemma output format)", metavar = "character"),
  make_option("--dbslmm", type = "character", default = NULL,
              help = "INPUT: the perfix of dbslmm software", metavar = "character"),
  make_option("--type", type = "character", default = "auto",
              help="INPUT: type of DBSLMM (default: default version)", metavar="character"),
  make_option("--model", type = "character", default = NULL,
              help = "INPUT: the chioce of LMM and DBSLMM", metavar = "character"),
  make_option("--reference", type = "character", default = NULL,
              help = "INPUT: the perfix of reference panel", metavar = "character"),
  make_option("--block", type = "character", default = NULL,
              help = "INPUT: the block information (Berisa and Pickrell 2015)", metavar = "character"),
  make_option("--N", type = "character", default = NULL,
              help="INPUT: sample size for summary statistics", metavar="character"),
  make_option("--outPath", type = "character", default = NULL,
              help="INPUT: the output path", metavar = "character"),
  ## Non-required parameters
  make_option("--h2f", type = "character", default = NULL,
              help = "INPUT: the fold of heritability of trait", metavar = "character"),
  make_option("--mafMax", type = "character", default = "0.2",
              help = "INPUT: the maximium of the difference between reference panel and summary data (default:0.2)", 
              metavar = "character"),
  make_option("--thread", type = "character", default = "1",
              help = "INPUT: the number of threads (default: 1)", 
              metavar = "character")
)
opt_parser <- OptionParser(option_list=args_list)
opt <- parse_args(opt_parser)

# opt <- list(summary = "", 
#             dbslmm = "/public/home/biostat03/biosoft/DBSLMM/scr/dbslmm",
#             type = "auto", 
#             model = "LMM", 
#             block = "/public/home/biostat03/biosoft/DBSLMM/block_data/EAS/",
#             reference = "/public/home/biostat03/project/pgsfusionProject/ref_geno/", 
#             N = "4479,75725",
#             outPath = "/public/home/biostat03/biosoft/DBSLMM/tmp_output/", 
#             mafMax = "0.2", 
#             thread = "1")

# Check the software
if (!file.exists(opt$dbslmm)){
  cat(paste0("ERROR: ", opt$dbslmm, " does not exist! Please check!\n"))
  q()
}

# Check files
if (!file.exists(opt$block)){
  cat(paste0("ERROR: ", opt$block, " does not exist! Please check!\n"))
  q()
}
if (!file.exists(opt$summary)){
  cat(paste0("ERROR: ", opt$summary, " does not exist! Please check!\n"))
  q()
}
if (!file.exists(opt$outPath)){
  cat(paste0("ERROR: ", opt$outPath, " does not exist! Please check!\n"))
  q()
}
if (!file.exists(opt$reference)){
  cat(paste0("ERROR: ", opt$outPath, " does not exist! Please check!\n"))
  q()
}

# Check settings
if (opt$type == "tuning"){
  
  if(is.null(opt$h2f)){
    cat(paste0("ERROR: Fold does not set! Please check!\n"))
    q()
  }
}

# Set parameters
R2 <- 0.2
PVAL <- 1e-6
WIN <- 1000

# Estimate heritability for whole genome
map <- transmute(readRDS(paste0(opt$reference,"/map_hm3_plus.rds")),
                 chr = as.integer(chr), pos = pos, rsid, af_val = af_ref,
                 a0 = a0, a1 = a1, ldsc = ldsc)
sumstats <- fread2(opt$summary, header = T)
sumstats <- sumstats[order(sumstats[, 1], sumstats[, 3]), ]
sumstats_gemma <- sumstats
colnames(sumstats_gemma) <- c("chr", "rs", "ps", "n_mis", "n_obs", "allele1", 
                              "allele0", "af", "beta", "se", "p_wald")
## Calculate effective sample size
N_str <- strsplit(opt$N, ",")[[1]] %>% as.numeric
if (length(N_str) == 1){
  
  n_eff <- N_str
} else {
  
  n_eff <- 4 / (1 / N_str[1] + 1 / N_str[2])
}
sumstats <- sumstats[, c(1, 3, 2, 6, 7, 8, 9, 10, 11)]
sumstats$n_eff <- n_eff
colnames(sumstats) <- c("chr", "pos", "rsid", "a1", "a0", "freq", "beta",
                        "beta_se", "pval", "n_eff")
## Coordinate summary statistics and map
df_beta <- as_tibble(snp_match(sumstats, map, strand_flip = FALSE, 
                               return_flip_and_rev = TRUE))
ldsc <- with(df_beta, snp_ldsc(ldsc, length(ldsc), chi2 = (beta / beta_se)^2,
                               sample_size = n_eff, blocks = NULL))
cat("heritability:", ldsc[[2]], "\n")
cat("number of SNPs: ", nrow(df_beta), "\n")

prefix_file <- strsplit(basename(opt$summary), "\\.")[[1]][1]

# DBSLMM-lmm
if (opt$model == "LMM"){

  lmm <- alply(c(1: 22), 1, function(CHR){

    sumstats_chr <- sumstats_gemma[sumstats_gemma$chr == CHR, ]
    chr_str <- paste0(opt$outPath, prefix_file, "_chr", CHR, ".assoc.txt")
    fwrite2(sumstats_chr, file = chr_str, col.names = F, sep = "\t")
    lmm_cmd <- paste0(opt$dbslmm,
                      " -s ",      chr_str,
                      " -r ",      paste0(opt$reference, "/merge"),
                      " -nsnp ",   nrow(df_beta),
                      " -n ",      as.integer(n_eff),
                      " -b ",      paste0(opt$block, "chr", CHR, ".bed"),
                      " -mafMax ", opt$mafMax,
                      " -t ",      opt$thread)
    if (opt$type == "auto"){
      ## auto model
      r_lmm_cmd <- paste0(lmm_cmd,
                          " -h ",      ldsc[[2]],
                          " -eff ",    opt$outPath, prefix_file, "_chr", CHR, ".dbslmm") %>% system()
      paste0("rm -rf ", opt$outPath, prefix_file, "_chr", CHR, ".dbslmm.badsnps") %>% system()
    } else {
      ## tuning model
      h2_vec <- as.numeric(unlist(strsplit(opt$h2f, ",")))
      for (hh in h2_vec) {

        # fwrite2(sumstats_chr, file = chr_str, col.names = F, sep = "\t")
        r_lmm_cmd <- paste0(lmm_cmd,
                            " -h ",      ldsc[[2]] * hh,
                            " -eff ",    opt$outPath, prefix_file, "_chr", CHR,  "_h2f", hh, ".dbslmm") %>% system()
        paste0("rm -rf ", opt$outPath, prefix_file, "_chr", CHR,  "_h2f", hh, ".dbslmm.badsnps") %>% system()
      }
    }

    ## remove temporary files
    system(paste0("rm -rf ", chr_str))
    return(CHR)
  })
}

# DBSLMM model
if (opt$model == "DBSLMM"){

  ## C+T
  df_beta_sig <- df_beta[df_beta$pval < PVAL, ]
  if(!file.exists(paste0(opt$reference, "/merge.bk")))
    snp_readBed(paste0(opt$reference, "/merge.bed"))
  ref_bed <- snp_attach(paste0(opt$reference, "/merge.rds"))
  ref_sub_str <- paste0(opt$outPath, "/merge_sub-", as.numeric(as.POSIXlt(Sys.time())))
  ref_sub_bed <- snp_attach(snp_subset(ref_bed,
                                       ind.col = df_beta_sig$`_NUM_ID_`,
                                       backingfile = ref_sub_str))
  # ref_G <- snp_fastImputeSimple(ref_sub_bed$genotypes)
  ref_G <- ref_sub_bed$genotypes
  lp_val <- -log10(df_beta_sig$pval)
  all_keep <- snp_grid_clumping(ref_G,
                                grid.thr.r2 = R2,
                                grid.base.size = WIN,
                                infos.chr = df_beta_sig$chr,
                                infos.pos = df_beta_sig$pos,
                                lpS = lp_val,
                                ncores = 1)
  system(paste0("rm -rf ", ref_sub_str, "*"))
  l_chr <- data.frame(chr = c(1: 22),
                      l_b = c(1: 22) %in% names(all_keep),
                      count = NA)
  l_chr$count[l_chr$l_b] <- c(1: sum(l_chr$l_b))
  ### Fit model
  dbslmm <- alply(c(1: 22), 1, function(CHR){

    dbslmm_cmd <- paste0(opt$dbslmm,
                         " -r ",      paste0(opt$reference, "/merge"),
                         " -nsnp ",   nrow(df_beta),
                         " -n ",      as.integer(n_eff),
                         " -b ",      paste0(opt$block, "chr", CHR, ".bed"),
                         " -mafMax ", opt$mafMax,
                         " -t ",      opt$thread)

    if(CHR %in% names(all_keep)){

      l_rsid_chr <- df_beta_sig$rsid.ss[all_keep[[l_chr$count[l_chr$chr == CHR]]][[1]]]
      l_summ_chr <- sumstats_gemma[sumstats_gemma$rs %in% l_rsid_chr, ]
      l_chr_str <- paste0(opt$outPath, "l_", prefix_file, "_chr", CHR, ".assoc.txt")
      s_summ_chr <- sumstats_gemma %>%
        filter(chr == CHR & !sumstats_gemma$rs %in% l_rsid_chr)
      s_chr_str <- paste0(opt$outPath, "s_", prefix_file, "_chr", CHR, ".assoc.txt")
      fwrite2(l_summ_chr, file = l_chr_str, col.names = F, sep = "\t")
      fwrite2(s_summ_chr, file = s_chr_str, col.names = F, sep = "\t")
      if (opt$type == "auto"){
        ## auto model
        r_dbslmm_cmd <- paste0(dbslmm_cmd,
                               " -l ",     l_chr_str,
                               " -s ",     s_chr_str,
                               " -h ",     ldsc[[2]],
                               " -eff ",   opt$outPath, prefix_file, "_chr", CHR, ".dbslmm") %>% system()
        system(paste0("rm -rf ", opt$outPath, prefix_file, "_chr", CHR, ".dbslmm.badsnps"))
      } else {
        ## tuning model
        h2_vec <- as.numeric(unlist(strsplit(opt$h2f, ",")))
        for (hh in h2_vec) {

          r_dbslmm_cmd <- paste0(dbslmm_cmd,
                                 " -l ",      l_chr_str,
                                 " -s ",      s_chr_str,
                                 " -h ",      ldsc[[2]] * hh,
                                 " -eff ",    opt$outPath, prefix_file, "_chr", CHR,  "_h2f", hh, ".dbslmm") %>% system()
          system(paste0("rm -rf ", opt$outPath, prefix_file, "_chr", CHR, "_h2f", hh, ".dbslmm.badsnps"))
        }
      }

      system(paste0("rm -rf ", l_chr_str))
    }else{

      sumstats_chr <- sumstats_gemma[sumstats_gemma$chr == CHR, ]
      s_chr_str <- paste0(opt$outPath, "s_", prefix_file, "_chr", CHR, ".assoc.txt")
      fwrite2(sumstats_chr, file = s_chr_str, col.names = F, sep = "\t")
      if (opt$type == "auto"){
        ## auto model
        r_dbslmm_cmd <- paste0(dbslmm_cmd,
                               " -s ",     s_chr_str,
                               " -h ",     ldsc[[2]],
                               " -eff ",   opt$outPath, prefix_file, "_chr", CHR, ".dbslmm") %>% system()
        system(paste0("rm -rf ", opt$outPath, prefix_file, "_chr", CHR, ".dbslmm.badsnps"))
      } else {
        ## tuning model
        h2_vec <- as.numeric(unlist(strsplit(opt$h2f, ",")))
        for (hh in h2_vec) {

          r_dbslmm_cmd <- paste0(dbslmm_cmd,
                                 " -s ",      s_chr_str,
                                 " -h ",      ldsc[[2]] * hh,
                                 " -eff ",    opt$outPath, prefix_file, "_chr", CHR,  "_h2f", hh, ".dbslmm") %>% system()
          system(paste0("rm -rf ", opt$outPath, prefix_file, "_chr", CHR, "_h2f", hh, ".dbslmm.badsnps"))
        }
      }
    }

    ## remove temporary files
    system(paste0("rm -rf ", s_chr_str))
    return(CHR)
  })
}

# # Format files
# paste0(opt$outPath, prefix_file, "_chr", 1:22, ".dbslmm.txt") %>% paste(., collapse = " ") %>%
#     paste("cat",. , ">", paste0(opt$outPath, prefix_file, ".dbslmm.txt")) %>% system()
# for (chr in 1: 22)
#   paste0("rm -rf ", opt$outPath, prefix_file, "_chr", chr, ".dbslmm") %>% system()
