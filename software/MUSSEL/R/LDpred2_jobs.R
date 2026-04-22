rm(list=ls())
suppressMessages(library("optparse"))
suppressMessages(library("bigsnpr"))
suppressMessages(library("bigreadr"))
suppressMessages(library("readr"))
suppressMessages(library("stringr"))
suppressMessages(library("caret"))
suppressMessages(library("Rcpp"))
suppressMessages(library("RcppArmadillo"))
suppressMessages(library("inline"))
suppressMessages(library("doMC"))
suppressMessages(library("foreach"))

option_list = list(
  make_option("--PATH_package", action="store", default=NA, type='character',
              help="Path to the directory where the downloaded files (decompressed) are saved [required]"),
  make_option("--PATH_data", action="store", default=NA, type='character',
              help="Path to the directory where the training data by ancestry group are saved [required]"),
  make_option("--PATH_LDref", action="store", default=NA, type='character',
              help="Path to the directory where the LD reference data by ancestry group and chromosome are saved [required]"),
  make_option("--PATH_out", action="store", default=NA, type='character',
              help="Path to the output directory where the results are saved [required]"),
  make_option("--FILE_sst", action="store", default=NA, type='character',
              help="Paths followed by file names of the population-specific GWAS summary statistics, separated by comma [required] [required columns: chr, rsid, pos, a0, a1, beta, beta_se, n_eff]"),
  make_option("--pop", action="store", default=NA, type='character',
              help="Populations of the GWAS samples, separated by comma [required]"),
  make_option("--p", action="store", default=paste(signif(seq_log(1e-5, 1, length.out = 21), 2), collapse = ','), type='character',
              help="Candidate values for tuning parameter p (causal SNP proportion) [default: %default]"),
  make_option("--H2", action="store", default=paste(c(0.3, 0.7, 1, 1.4), collapse = ','), type='character',
              help="Candidate values for tuning parameter H2 (heritability = H2 * h2_est from LDSC) [default: %default]"),
  make_option("--sparse", action="store", default='0', type='character',
              help="Whether to consider a sparse model: 0, 1, or 0,1 [default: %default]"),
  make_option("--bfile_tuning", action="store", default=NA, type='character',
              help="Path to PLINK binary input file prefix (minus .bed/.bim/.fam) for tuning, save by chromosome [required]"),
  
  make_option("--verbose", action="store", default=1, type="integer",
              help="How much chatter to print: 0=nothing; 1=minimal; 2=all [default: %default]"),
  make_option("--cleanup", action="store", default=T, type="logical",
              help="Cleanup temporary files or not [default: %default]"),
  make_option("--NCORES", action="store", default=12, type="integer",
              help="How many cores to use [default: %default]")
)
opt = parse_args(OptionParser(option_list=option_list))

if ( opt$verbose == 2 ) { SYS_PRINT = F } else { SYS_PRINT = T }

suppressWarnings(dir.create(opt$PATH_out))

races = str_split(opt$pop,",")[[1]]; K <- length(races)
sumdata_paths = str_split(opt$FILE_sst,",")[[1]]

bfile_tuning_vec <- str_split(opt$bfile_tuning,",")[[1]]


# Perform i/o check:
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

filen<-paste0(rscripts_path, 'LDpred2_rscript.sh')
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
cat(paste0('Rscript ',opt$PATH_package,'/R/LDpred2.R '), file = zz, sep = " ")
cat(paste0(' --PATH_package ', opt$PATH_package), file = zz, sep = " ")
cat(paste0(' --PATH_ref ', opt$PATH_LDref), file = zz, sep = " ")
cat(paste0(' --PATH_out ', opt$PATH_out), file = zz, sep = " ")
cat(paste0(' --FILE_sst ', opt$FILE_sst), file = zz, sep = " ")
cat(paste0(' --pop ', opt$pop), file = zz, sep = " ")
cat(paste0(' --p ', opt$p), file = zz, sep = " ")
cat(paste0(' --H2 ', opt$H2), file = zz, sep = " ")
cat(paste0(' --sparse ', opt$sparse), file = zz, sep = " ")
cat(paste0(' --bfile_tuning ', opt$bfile_tuning), file = zz, sep = " ")
cat(paste0(' --NCORES ', opt$NCORES), file = zz, sep = "\n")
cat("\n", file=zz)
cat(paste0('wait'), file = zz, sep = " ")

system(paste0('sbatch --mem=25G ' ,rscripts_path, 'LDpred2_rscript.sh'))

print(paste0('R scripts submitted for running LDpred2.'))

cat(paste0('\n** Wait until the following files have been saved: **\n'))
cat(paste0('\n** ', opt$PATH_out, '/{', paste(races,collapse = ','), '}', '/tmp/beta_files/beta_in_all_settings/params.txt **\n'))
cat(paste0('\n** Done. **\n'))
