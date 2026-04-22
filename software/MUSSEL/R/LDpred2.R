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
  make_option("--PATH_ref", action="store", default=NA, type='character',
              help="Path to the directory where the LD reference data by ancestry group and chromosome are saved [required]"),
  make_option("--PATH_out", action="store", default=NA, type='character',
              help="Path to the output directory where the results are saved [required]"),
  make_option("--FILE_sst", action="store", default=NA, type='character',
              help="Paths followed by file names of the population-specific GWAS summary statistics, separated by comma [required] [required columns: chr, rsid, pos, a0, a1, beta, beta_se, n_eff]"),
  make_option("--pop", action="store", default=NA, type='character',
              help="Populations of the GWAS samples, separated by comma [required]"),
  
  make_option("--p", action="store", default=paste(signif(seq_log(1e-5, 1, length.out = 11), 2), collapse = ','), type='character',
              help="Candidate values for tuning parameter p (causal SNP proportion) [default: %default]"),
  make_option("--H2", action="store", default=paste(c(0.7, 1, 1.4), collapse = ','), type='character',
              help="Candidate values for tuning parameter H2 (heritability = H2 * h2_est from LDSC) [default: %default]"),
  make_option("--sparse", action="store", default='0', type='character',
              help="Whether to consider a sparse model: 0, 1, or 0,1 [default: %default]"),
  make_option("--bfile_tuning", action="store", default=NA, type='character',
              help="Path to PLINK binary input file prefix (minus .bim) for tuning, save by chromosome [required]"),

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
ref_paths <- paste0(opt$PATH_ref, '/', races, "/1KGref_plinkfile/")
precalLD_paths <- paste0(opt$PATH_ref, '/', races, '/LDpred2_lassosum2_corr_1kg/')
map_paths <- paste0(opt$PATH_ref, '/', races, '/map')
out_paths <- paste0(opt$PATH_out,"/",races)

bfile_tuning_vec <- str_split(opt$bfile_tuning,",")[[1]]

# Perform i/o check:
files <- NULL
files <- c(files, sumdata_paths)
for(mmm in 1:K){ files <- c(files, paste(bfile_tuning_vec[mmm],c(".bim"),sep='')) }

for ( f in files ) {
  if ( !file.exists(f) ){
    cat( "ERROR: ", f , " files for tuning data do not exist\n" , sep='', file=stderr() )
    q()
  }
}
rm(list="files")

NCORES <- opt$NCORES

temp <- commandArgs(TRUE)
race = races[as.numeric(temp[1])]
chr =  as.numeric(temp[2])

source(paste0(opt$PATH_package,"/R/source-functions.R"))
ldr = 3/1000 # default ld radius in LDpred2
ref_bim_path = strsplit(opt$PATH_package, '/')[[1]]
ref_bim_path = paste0(ref_bim_path[-length(ref_bim_path)], collapse = '/')
ref <- fread2(paste0(ref_bim_path,"/PROSPER/ref_bim.txt"))
p_seq <- as.numeric(strsplit(opt$p, split = ',')[[1]])
H2s = as.numeric(strsplit(opt$H2, split = ',')[[1]])
sparse = as.numeric(strsplit(opt$sparse, split = ',')[[1]])



for(mmm in 1:K){
  race <- races[mmm]
  sumdata_path <- sumdata_paths[mmm]
  ref_path <- ref_paths[mmm]
  precalLD_path <- precalLD_paths[mmm]
  map_path <- map_paths[mmm]
  out_path <- out_paths[mmm]
  
  bfile_tuning <- bfile_tuning_vec[mmm]

  suppressWarnings(dir.create(paste0(out_path)))
  suppressWarnings(dir.create(paste0(out_path, "/tmp")))
  suppressWarnings(dir.create(paste0(out_path, "/tmp/beta_files")))
  suppressWarnings(dir.create(paste0(out_path, "/tmp/beta_files/beta_in_all_settings")))

  if ( opt$verbose >= 1 ) {
    cat(paste0("\n*************************"))
    cat(paste0("\n**** LDpred2 on ",race," ****"))
    cat(paste0("\n*************************\n"))
  }
    
  ########################################################################
  ########################################################################
  if ( opt$verbose >= 1 ) cat(paste0('\n** Step 1: data preparation **\n'))
  
  ############
  map_ldref <- readRDS(paste0(map_path, '/map_1kg_ldref.rds'))
  
  ## Load GWAS summary data, only keep SNPs that are present in both LD reference samples and tuning samples (according to rsid)
  sum.raw <-bigreadr::fread2(sumdata_path)
  sum.raw <- sum.raw[sum.raw$rsid %in% ref$V2,]
  tun <- bigreadr::fread2(paste0(bfile_tuning,'.bim'))
  sum.raw <- sum.raw[sum.raw$rsid %in% tun$V2,]
  ref_tmp <- ref[match(sum.raw$rsid, ref$V2),]
  
  sumstats = sum.raw[,c("chr", "rsid", "a1", "a0", "beta", "beta_se", "n_eff")]
  names(sumstats) <- c("chr", "rsid", "a0", "a1", "beta", "beta_se", "n_eff") # a0: ref allele
  info_snp <- snp_match(sumstats, map_ldref, strand_flip = T, join_by_pos = F) # important: for real data, strand_flip = T
  info_snp <- tidyr:: drop_na(tibble::as_tibble(info_snp))
  df_beta <- info_snp
  print(paste0('Complete pre-processing GWAS summary data.'))
  
  td = paste0(out_path, '/intermediate/')
  if (!dir.exists(td)) dir.create(td)
  setwd(td)
  tmp <- tempfile(tmpdir = td)
  ld = NULL
  for (chr in 1:22) {
    cat(chr, ".. ", sep = "")
    ## indices in 'df_beta'
    ind.chr <- which(df_beta$chr == chr)
    ## indices in 'map_ldref'
    ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
    ## indices in 'corr0'
    ind.chr3 <- match(ind.chr2, which(map_ldref$chr == chr))
    if (length(ind.chr3) > 0){
      # corr0
      corr0 <- readRDS(paste0(precalLD_path, '/LD_ref_chr', chr, '.rds'))[ind.chr3, ind.chr3]
      if ((chr == 1) | (is.null(ld))) {
        ld <- Matrix::colSums(corr0^2)
        corr <- as_SFBM(corr0, tmp, compact = TRUE)
      } else {
        if (length(corr0) == 1) corr0 =  as(1, "sparseMatrix") # as(corr0, "sparseMatrix") # as.matrix(corr0, 1, 1)
        ld <- c(ld, Matrix::colSums(corr0^2))
        corr$add_columns(corr0, nrow(corr))
      }
      print(paste0('Complete summarizing LD matrix info for CHR ', chr))
      rm(corr0)
    }
  }
  rm(sum.raw, tun, ref_tmp, sumstats)
  
  # --------------------- Step 2.2: Run LDpred2 ---------------------
  (ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2,
                                  sample_size = n_eff, blocks = NULL)))
  ldsc_h2_est <- ldsc[["h2"]]
  h2_seq <- round(ldsc_h2_est * H2s, 4)
  params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = sparse)
  
  if ( opt$verbose >= 1 ) cat("\n** Step 2: run LDpred2 under all tuning parameter settings **\n")
  set.seed(2023)
  beta_ldpred2 <- snp_ldpred2_grid(corr, df_beta, params, ncores = NCORES)
  beta_ldpred2[is.na(beta_ldpred2)] = 0; beta_ldpred2[abs(beta_ldpred2) > 5] = 0
  beta_ldpred2 = data.frame(df_beta[,c('chr','rsid','a0','a1')], beta_ldpred2)
  rownames(beta_ldpred2) = df_beta$rsid
  colnames(beta_ldpred2) = c('chr','rsid','a0','a1', paste0('e',1:nrow(params)))
  output_LDpred2 = paste0(out_path, '/tmp/beta_files/beta_in_all_settings/ldpred2effects.txt')
  # write_delim(beta_ldpred2,file = output_LDpred2, delim='\t')
  write.table(beta_ldpred2,file = output_LDpred2, col.names = T,row.names = F,quote=F)
  
  # The rows are in the same order as the columns in beta_grid.
  params <- expand.grid(p = p_seq, h2 = H2s, sparse = sparse)
  fwrite2(params, paste0(out_path, '/tmp/beta_files/beta_in_all_settings/params.txt'), col.names = T, sep="\t", nThread=1)
  
  if ( opt$verbose >= 1 ) cat(paste0("\n** Completed: results saved in ", paste0(out_path, "/tmp/beta_files/beta_in_all_settings/params.txt","**\n")))
  ############
  rm(list=c("corr"))
}

  

