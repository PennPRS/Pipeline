rm(list=ls())
suppressMessages(library(optparse))
suppressMessages(library(bigreadr))
suppressMessages(library(bigsnpr))
suppressMessages(library(bigparallelr))
suppressMessages(library(bigmemory))
suppressMessages(library(stringr))
suppressMessages(library(caret))
suppressMessages(library(Rcpp))
suppressMessages(library(RcppArmadillo))
suppressMessages(library(RcppTN))
suppressMessages(library(inline))
suppressMessages(library(doMC))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))

suppressMessages(library(data.table))
suppressMessages(library(readr))
suppressMessages(library(MASS)) # for mvrnorm and ginv
suppressMessages(library(reshape)) # for melt
suppressMessages(library(parallel))
suppressMessages(library(devtools))
# suppressMessages(library(mvnfast))
suppressMessages(library(genio)) # for read_plink
suppressMessages(library(dplyr))
# suppressMessages(library(gdata)) 
# suppressMessages(library(R.utils)) # for gzip
# suppressMessages(library(pROC))

suppressMessages(library(pryr))
suppressMessages(library(Matrix))
suppressMessages(library(lavaan))
suppressMessages(library(xtable))
# library(corpcor) #for pseudoinverse
# suppressMessages(library(DescTools))
suppressMessages(library(SuperLearner))

option_list = list(
  make_option("--PATH_package", action="store", default=NA, type='character',
              help="Path to the directory where the downloaded files (decompressed) are saved [required]"),
  make_option("--PATH_data", action="store", default=NA, type='character',
              help="Path to the directory where the training data are saved [required]"),
  make_option("--PATH_LDref", action="store", default=NA, type='character',
              help="Path to the directory where the LD reference data by ancestry group and chromosome are saved [required]"),
  make_option("--PATH_out", action="store", default=NA, type='character',
              help="Path to the output directory where the results are saved [required]"),
  make_option("--FILE_sst", action="store", default=NA, type='character',
              help="Paths followed by file names of the population-specific GWAS summary statistics, separated by comma [required] [required columns: chr, rsid, pos, a0, a1, beta, beta_se, n_eff]"),
  make_option("--pop", action="store", default=NA, type='character',
              help="Populations of the GWAS samples, separated by comma [required]"),
  make_option("--LDpred2_params", action="store", default=NA, type='character',
              help="Path to the directory where the tuned LDpred2 parameters (population-specific causal SNP proportions, heritability and whether or not a sparse model is used) are saved, separated by comma [required]"),
  make_option("--chrom", action="store", default="1-22", type='character',
              help="The chromosome on which the model is fitted, input in the format of 
              1-22 or 1,2,3 [required]"),
  
  make_option("--cors_additional", action="store", default=NA, type='character',
              help="Additional candidate values for tuning parameter: genetic correlation across ancestry groups, example: 3 groups with label 1,2,3, want to add two additional settings: cor_setting1(1,2),cor_setting1(1,3),cor_setting1(2,3);cor_setting2(1,2),cor_setting2(1,3),cor_setting2(2,3) [optional]"),
  make_option("--ps_additional", action="store", default=NA, type='character',
              help="Typically not necessary. Additional candidate values for tuning parameter: ancestry-specific causal SNP proportions, example: 3 groups with label 1,2,3, want to add two additional settings: p1_setting1,p2_setting1,p3_setting1,p1_setting2,p2_setting2,p3_setting2 [optional]"),
  
  make_option("--bfile_tuning", action="store", default=NA, type='character',
              help="Path to PLINK binary input file prefix (minus .bed/.bim/.fam) for tuning, save by chromosome [required]"),
  make_option("--verbose", action="store", default=1, type="integer",
              help="How much chatter to print: 0=nothing; 1=minimal; 2=all [default: %default]"),
  make_option("--cleanup", action="store", default=T, type="logical",
              help="Cleanup temporary files or not [default: %default]"),
  make_option("--NCORES", action="store", default=5, type="integer",
              help="How many cores to use [default: %default]")
)
opt = parse_args(OptionParser(option_list=option_list))

if ( opt$verbose == 2 ) { SYS_PRINT = F } else { SYS_PRINT = T }

NCORES <- opt$NCORES
races = str_split(opt$pop,",")[[1]]; K <- length(races)
sumdata_paths = str_split(opt$FILE_sst,",")[[1]]
out_paths <- paste0(opt$PATH_out,"/",races)
opt$chrom <- gsub("-",":",opt$chrom)
eval(parse(text=paste0("chrs = c(",opt$chrom,")")))
ldpred2_params_path = str_split(opt$LDpred2_params,",")[[1]]
bfile_tuning_vec <- str_split(opt$bfile_tuning,",")[[1]]


# Perform i/o checks here:
for ( f in sumdata_paths ) {
  if ( !file.exists(f) ){
    cat( "ERROR: input GWAS summary data files", f , " do not exist \n" , sep='', file=stderr() )
    q()
  }
}

files <- NULL
files <- c(files, sumdata_paths)
for(mmm in 1:K){ files <- c(files, paste(bfile_tuning_vec[mmm],c(".bed",".bim",".fam"),sep='')) }

for ( f in files ) {
  if ( !file.exists(f) ){
    cat( "ERROR: ", f , " files for tuning data do not exist\n" , sep='', file=stderr() )
    q()
  }
}
rm(list="files")

NCORES <- opt$NCORES

source(paste0(opt$PATH_package,"/R/source-functions.R"))

# First, write rscript files for running LDpred2 by chromosome
rscripts_path = paste0(opt$PATH_data, '/rscripts/')
suppressWarnings(dir.create(rscripts_path))
suppressWarnings(dir.create(paste0(rscripts_path,'logfile')))


for (chr in chrs){
  filen<-paste0(rscripts_path, 'MUSS_rscript_chr', chr, ".sh")
  system(paste0('rm -rf ',filen))
  file.create(filen)
  zz <- file(filen, "w")
  cat("#!/bin/bash", "", file = zz, sep = "\n")
  cat("\n", file=zz)
  cat("#$ -cwd", "", file = zz, sep = "\n")
  cat(paste0('module load R'), file = zz, sep = "\n")
  cat("\n", file=zz)
  cat(paste0('#$ -o ',rscripts_path,'logfile'), file = zz, sep = "\n")
  cat(paste0('#$ -e ',rscripts_path,'logfile'), file = zz, sep = "\n")
  cat("\n", file=zz)
  cat(paste0('Rscript ',opt$PATH_package,'/R/MUSS.R '), file = zz, sep = " ")
  cat(paste0(' --PATH_package ', opt$PATH_package), file = zz, sep = " ")
  cat(paste0(' --PATH_LDref ', opt$PATH_LDref), file = zz, sep = " ")
  cat(paste0(' --PATH_out ', opt$PATH_out), file = zz, sep = " ")
  cat(paste0(' --FILE_sst ', opt$FILE_sst), file = zz, sep = " ")
  cat(paste0(' --pop ', opt$pop), file = zz, sep = " ")
  cat(paste0(' --LDpred2_params ', opt$LDpred2_params), file = zz, sep = " ")
  cat(paste0(' --chrom ', chr), file = zz, sep = " ")
  if (!is.na(opt$cors_additional)) cat(paste0(' --cors_additional ', opt$cors_additional), file = zz, sep = " ")
  if (!is.na(opt$ps_additional)) cat(paste0(' --ps_additional ', opt$ps_additional), file = zz, sep = " ")
  cat(paste0(' --bfile_tuning ', opt$bfile_tuning), file = zz, sep = " ")
  cat(paste0(' --NCORES ', opt$NCORES), file = zz, sep = "\n")
  cat("\n", file=zz)
  cat(paste0('wait'), file = zz, sep = " ")
}

for (chr in chrs){
  if ((chr >= 1)&(chr < 13)) system(paste0('sbatch --mem=30G ' ,rscripts_path, 'MUSS_rscript_chr', chr, ".sh"))
  if ((chr >= 13)&(chr < 23)) system(paste0('sbatch --mem=15G ' ,rscripts_path, 'MUSS_rscript_chr', chr, ".sh"))
}

print(paste0('R scripts submitted for running MUSS by chromosome in parallel.'))

cat(paste0('\n** Wait until all following files have been saved: **\n'))
cat(paste0('\n** ', opt$PATH_out, '/tmp/MUSS_beta_in_all_settings_bychrom/{', paste(races,collapse = ','), '}-chr{1..22}.txt **\n'))
cat(paste0('\n** Done. **\n'))








