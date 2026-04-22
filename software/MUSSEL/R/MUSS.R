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
suppressMessages(library(genio)) # for read_plink
suppressMessages(library(dplyr))

# suppressMessages(library(pryr))
suppressMessages(library(Matrix))
suppressMessages(library(lavaan))
suppressMessages(library(xtable))
suppressMessages(library(SuperLearner))

option_list = list(
  make_option("--PATH_package", action="store", default=NA, type='character',
              help="Path to the directory where the downloaded files (decompressed) are saved [required]"),
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
  make_option("--chrom", action="store", default=NA, type='integer',
              help="The chromosome on which the model is fitted [required]"),
  
  make_option("--cors_additional", action="store", default=NA, type='character',
              help="Additional candidate values for tuning parameter: genetic correlation across ancestry groups, example: 3 groups with label 1,2,3, want to add two additional settings: cor_setting1(1,2),cor_setting1(1,3),cor_setting1(2,3);cor_setting2(1,2),cor_setting2(1,3),cor_setting2(2,3) [optional]"),
  make_option("--ps_additional", action="store", default=NA, type='character',
              help="Typically not necessary. Additional candidate values for tuning parameter: ancestry-specific causal SNP proportions, example: 3 groups with label 1,2,3, want to add two additional settings: p1_setting1,p2_setting1,p3_setting1;p1_setting2,p2_setting2,p3_setting2 [optional]"),
  
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
ref_paths <- paste0(opt$PATH_LDref, "/", races, "/1KGref_plinkfile/") # raw LD reference genotype data by ancestry group and chromosome
LD_paths <- paste0(opt$PATH_LDref, "/", races, "/LD") # precalculated block-wise reference LD matrices by ancestry group and chromosome (generated based on the raw reference genotype data)
map_paths <- paste0(opt$PATH_LDref,"/",races, "/map")
out_paths <- paste0(opt$PATH_out,"/",races)
chr = opt$chrom
ldpred2_params_path = str_split(opt$LDpred2_params,",")[[1]]
bfile_tuning_vec <- str_split(opt$bfile_tuning,",")[[1]]

# Perform i/o checks here:
for ( f in sumdata_paths ) {
  if ( !file.exists(f) ){
    cat( "ERROR: input GWAS summary data files", f , " do not exist \n" , sep='', file=stderr() )
    q()
  }
}

# cors_additional = ps_additional = NA
if (opt$cors_additional == 'NA') opt$cors_additional = NA
if (!is.na(opt$cors_additional)){
  tem = strsplit(opt$cors_additional, split = ';')[[1]]
  cors_additional = lapply(1:2, function(x){as.numeric(strsplit(tem[[x]], split = ',')[[1]])})
  validinput = sum(sapply(1:length(cors_additional), function(x){length(cors_additional[[x]]) == (K * (K-1))/2}))
  if (validinput != length(cors_additional)) cat(paste0('Error: wrong input format for --cors_additional \n'))
}
if (is.na(opt$cors_additional)) cors_additional = opt$cors_additional
if (opt$ps_additional == 'NA') opt$ps_additional = NA
if (!is.na(opt$ps_additional)){
  tem = strsplit(opt$ps_additional, split = ';')[[1]]
  ps_additional = lapply(1:2, function(x){as.numeric(strsplit(tem[[x]], split = ',')[[1]])})
  validinput = sum(sapply(1:length(ps_additional), function(x){length(ps_additional[[x]]) == K}))
  if (validinput != length(ps_additional)) cat(paste0('Error: wrong input format for --ps_additional \n'))
}
if (is.na(opt$ps_additional)) ps_additional = opt$ps_additional

suppressWarnings(dir.create(opt$PATH_out))
suppressWarnings(dir.create(paste0(opt$PATH_out, "/tmp")))
suppressWarnings(dir.create(paste0(opt$PATH_out, "/tmp/MUSS_beta_in_all_settings_bychrom")))
suppressWarnings(dir.create(paste0(opt$PATH_out, "/MUSS/"))) # before ensemble by SL

sourceCpp(paste0(opt$PATH_package,"/src/MUSS.cpp"))
source(paste0(opt$PATH_package,"/R/source-functions.R"))

########################################################################
########################################################################

if ( opt$verbose >= 1 ) cat("\n** Step 1: data preparation. **\n")

############
## Step 1.1. Load GWAS summary data, only keep SNPs that are present in both LD reference samples and tuning samples (according to rsid)

# ref <- fread2(paste0(opt$PATH_package,"/ref_bim.txt")) # reference SNP information combined across all chromosomes
ref_bim_path = strsplit(opt$PATH_package, '/')[[1]]
ref_bim_path = paste0(ref_bim_path[-length(ref_bim_path)], collapse = '/')
ref <- fread2(paste0(ref_bim_path,"/PROSPER/ref_bim.txt"))

sum.raw = tun = list()
for (k in 1:K){
  sum.raw[[k]] = bigreadr::fread2(sumdata_paths[k])
  sum.raw[[k]] = sum.raw[[k]][,c('rsid','chr','a1','a0','beta','beta_se','n_eff')]
  colnames(sum.raw[[k]]) = c('SNP_ID','CHR','REF','ALT','BETA','SE','N') # REF: allele corresponding to BETA
  sum.raw[[k]] = sum.raw[[k]] %>% filter((SNP_ID %in% ref$V2) & (CHR == chr))
  tun[[k]] <- bigreadr::fread2(paste0(bfile_tuning_vec[k],'.bim'))
  tun[[k]] = tun[[k]] %>% filter((!duplicated(V2)) & (V1 == chr))
  colnames(tun[[k]]) = c('chrom','sid','na','pos','nt1','nt2')
  snps = intersect(tun[[k]]$sid, sum.raw[[k]]$SNP_ID)
  sum.raw[[k]] = sum.raw[[k]] %>% filter(SNP_ID %in% snps)
  tun[[k]] <- tun[[k]] %>% filter(sid %in% snps)
}


## Step 1.2. Delete SNPs with ambiguous alleles:
snp.remove = list()
l=1
for (k1 in 1:(K-1)){
  for (k2 in (k1+1):K){
    bim1 = tun[[k1]]; bim2 = tun[[k2]]
    sharedsnp = intersect(sum.raw[[k1]]$SNP_ID, sum.raw[[k2]]$SNP_ID)
    bim1 = bim1 %>% filter(sid %in% sharedsnp)
    bim1$sid = as.character(bim1$sid)
    bim1 = bim1[,c('sid','nt1','nt2')]
    colnames(bim1)[1] = 'SNP_ID'
    
    bim2 = bim2 %>% filter(sid %in% sharedsnp)
    bim2$sid = as.character(bim2$sid)
    bim2 = bim2[,c('sid','nt1','nt2')]
    colnames(bim2)[1] = 'SNP_ID'
    
    bim = merge(bim1, bim2, by='SNP_ID')
    rownames(bim) = bim$SNP_ID
    bim = bim[sharedsnp,]
    
    ## --------- scenario 2: A1 A2 flipped between the biomarker and the outcome also modifiable strand-ambiguous:
    alleles = cbind(bim$nt1.x, bim$nt2.x, bim$nt1.y, bim$nt2.y)
    combine.alleles = function(x) paste(x,collapse='')
    alleles = apply(alleles,1,combine.alleles)
    inds.ambiguous= which(alleles %in% c('ACGT','AGCT','TCGA','TGCA','CATG','CTAG','GATC','GTAC'))
    inds.ambiguous.keep = which(alleles %in% c('ACTG','AGTC','TCAG','TGAC','CAGT','CTGA','GACT','GTCA'))
    matched = which((bim$nt1.x == bim$nt1.y)&(bim$nt2.x == bim$nt2.y))
    flipped = which((bim$nt1.x == bim$nt2.y)&(bim$nt2.x == bim$nt1.y))
    snp.remove[[l]] = sharedsnp[-c(inds.ambiguous, inds.ambiguous.keep, matched, flipped)]
    rm(sharedsnp)
    l = l + 1
  }
}
snp.remove = unlist(snp.remove)
for (k in 1:K) sum.raw[[k]] = sum.raw[[k]] %>% filter(!(SNP_ID %in% snp.remove))

## Step 1.3. Match alleles between GWAS summary data and reference panel, summarize LD information by ancestry group
snps = summ = snps.qc = C.sfbm = beta_inf = NSNPS = SNPS = LD = GENO = INDX = list()
hsq = numeric()
for (k in 1:K){
  ref_path <- ref_paths[k]
  LD_path <- LD_paths[k]
  map_path <- map_paths[k]
  out_path <- out_paths[k]
  
  map <- readRDS(paste0(map_path,'/map_MUSS_chr',chr,'.rds'))
  ## Load LD information (saved in a list) by LD block:
  # Load snps_list (rsid), Nsnps: 
  load(paste0(LD_path, '/standard_data/chr',chr,'_snps.RData'))
  # Load LD_list
  load(paste0(LD_path, '/standard_data/chr',chr,'_LD.RData'))
  remove.indx =  which(sapply(1:length(snps_list), function(x){length(snps_list[[x]])}) == 0)
  if (length(remove.indx) > 0){
    snps_list = snps_list[-remove.indx]
    Nsnps = Nsnps[-remove.indx]
  } 
  snpslist = unlist(snps_list)
  
  set.seed(2020)
  sumstats = sum.raw[[k]][sum.raw[[k]]$CHR == chr,c('CHR', 'SNP_ID', 'REF', 'ALT', 'BETA', 'SE', 'N')]
  names(sumstats) <- c("chr", "rsid", "a0", "a1", "beta", "beta_se", "n_eff")
  rsid_int = intersect(map$rsid, snpslist)
  sumstats = sumstats[sumstats$rsid %in% rsid_int,]
  info_snp <- snp_match(sumstats, map, strand_flip = T, join_by_pos = F)
  rownames(info_snp) = info_snp$rsid
  
  ind.chr <- which(info_snp$chr == chr)
  df_beta <- info_snp[ind.chr, c("beta", "beta_se", "n_eff")]
  ## indices in G
  ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]

  # Remove SNPs that do not exist in snps_list:
  remove.indx =  which(sapply(1:length(LD_list), function(x){length(LD_list[[x]])}) == 0)# c(1,length(LD_list))
  if (length(remove.indx) > 0) LD_list = LD_list[-remove.indx]
  allele_folder = paste0(LD_path, '/tmp/byblock/chr',chr,'/')
  # if (length(remove.indx) == 2) num_files <- (length(list.files(allele_folder))-2)/4
  # if (length(remove.indx) == 1) num_files <- (length(list.files(allele_folder))-1)/4
  # if (length(remove.indx) == 0) num_files <- (length(list.files(allele_folder)))/4
  allele_files = list.files(allele_folder)
  allele_files = allele_files[which(str_detect(allele_files, c('.bim')))]
  allele_list = list()
  for (l in 1:length(allele_files)){
    allele_list[[l]] = bigreadr::fread2(paste0(allele_folder,allele_files[l]))
    rownames(allele_list[[l]]) = allele_list[[l]][,2]
  }
  
  ## Summarize LD and corresponding SNP information by LD block
  snps_list0 = snps_list; snps_list = list()
  geno_list = list()
  for (l in 1:length(snps_list0)){
    geno_list[[l]] = character()
    if (!is.null(snps_list0[[l]])){
      indx = which(snps_list0[[l]] %in% rownames(df_beta))
      Nsnps[l] = length(indx)
      snps_list[[l]] = snps_list0[[l]][indx]
      
      LD_list[[l]] = LD_list[[l]][indx, indx]
      temrsid = snps_list[[l]]
      for (l2 in 1:length(allele_files)){
        te = intersect(temrsid, allele_list[[l2]][,2])
        if (length(te)>0){
          geno_list[[l]] = allele_list[[l2]][snps_list[[l]],]
        }
      }
      if (is.null(geno_list[[l]])) print('Error: no SNP remaining after filtering.')
    }
  }
  
  ## Construct sparse LD correlation matrix, C.sfbm[[k]]:
  df_beta0 = df_beta; df_beta0$indx = 1:nrow(df_beta); indx = list(); xindx = yindx = entries = numeric()
  for (l in 1:length(LD_list)){
    # if (!is.null(snps_list[[l]])){
    if (length(snps_list[[l]]) > 0){
      indx[[l]] = df_beta0[snps_list[[l]],'indx']
      te = expand.grid(x = indx[[l]], y = indx[[l]])
      xindx = c(xindx, te$x); yindx = c(yindx, te$y)
      entries = c(entries, as.numeric(LD_list[[l]]))
      # flipped snps:
      matched = which((info_snp[indx[[l]],'a0'] == geno_list[[l]][,5])&(info_snp[indx[[l]],'a1'] == geno_list[[l]][,6]))
      flipped = which((info_snp[indx[[l]],'a0'] == geno_list[[l]][,6])&(info_snp[indx[[l]],'a1'] == geno_list[[l]][,5]))
      info_snp[indx[[l]][flipped],'a0'] = geno_list[[l]][flipped,5]; info_snp[indx[[l]][flipped],'a1'] = geno_list[[l]][flipped,6];
      info_snp[indx[[l]][flipped],'beta'] = -info_snp[indx[[l]][flipped],'beta']; 
      df_beta[indx[[l]][flipped],'beta'] = -df_beta[indx[[l]][flipped],'beta']
    }
  }
  corr0 <- sparseMatrix(i=xindx,j=yindx,x=entries)
  C.sfbm[[k]] <- bigsparser::as_SFBM(as(corr0, "dgCMatrix"))
  
  ldsc <- snp_ldsc2(corr0, df_beta)
  hsq[k] <- abs(ldsc[["h2"]])
  
  # Obtain initial values of beta on the standardized scale for MCMC:
  Ns = info_snp$n_eff
  diag.add = ncol(C.sfbm[[k]]) / (hsq[k] * Ns)
  beta_inf[[k]] = bigsparser::sp_solve_sym(C.sfbm[[k]], info_snp$beta, add_to_diag = diag.add) # on the standardized scale
  
  summ[[k]] = info_snp; rownames(summ[[k]]) = summ[[k]]$rsid
  snps.qc[[k]] = summ[[k]]$rsid # # identical to info_snp$rsid
  
  colnames(summ[[k]]) = c('CHR','SNP_ID','REF','ALT','BETA','SE','N','NUM_ID_SS','pos','NUM_ID')
  summ[[k]] = summ[[k]][,c('CHR','SNP_ID','REF','ALT','BETA','SE','N')]
  rm(corr0)
  NSNPS[[k]] = Nsnps; SNPS[[k]] = snps_list; LD[[k]] = LD_list; GENO[[k]] = geno_list; INDX[[k]] = indx
  if ( opt$verbose >= 1 ) cat(paste0('\n** ',races[k],' LD information preparation completed. **\n'))
}
names(hsq) = races; names(C.sfbm) = paste0('C',1:K,'.sfbm')
rm(sum.raw,sumstats)

# ---------------- Standardize BETA and SE ----------------
for (k in 1:K){
  summ[[k]]$c = (summ[[k]]$SE*sqrt(summ[[k]]$N)) #c: scale
  summ[[k]]$Bh1 = summ[[k]]$BETA/(summ[[k]]$SE*sqrt(summ[[k]]$N))
  summ[[k]]$SE = summ[[k]]$SE/(summ[[k]]$SE*sqrt(summ[[k]]$N))
  rownames(summ[[k]]) = summ[[k]]$SNP_ID
  colnames(summ[[k]])[which(colnames(summ[[k]]) == 'Bh1')] = paste0('Bh',k)
}

# ---------------- redefine SNP indices ----------------
snps = unique(unlist(snps.qc))
snpinfo = data.frame(rsid = snps, ind0 = 1:length(snps))
snpinfo$rsid = as.character(snpinfo$rsid)
rownames(snpinfo) = as.character(snpinfo$ind0)
temmat = as.data.frame(matrix(0,nrow(snpinfo),K))
colnames(temmat) = paste0('ind',1:K)
snpinfo = cbind(snpinfo, temmat)

# ---------------- qc ----------------
for (k in 1:K){
  tem = which(!(snpinfo$rsid %in% snps.qc[[k]]))
  snpinfo[,paste0('ind',k)][tem] = 0
  indtemp = lapply(c(1:nrow(snpinfo))[-tem],function(x) {which(snps.qc[[k]] == snpinfo$rsid[x])})
  snpinfo[,paste0('ind',k)][c(1:nrow(snpinfo))[-tem]] = unlist(indtemp)
  # print(paste0('Complete generating indtemp ',k))
}
rm(indtemp,tem,snps)

# ---------------- Summarize remaining information required:
Ns = sigmasq = Mt = tem = Bh = beta_init = list()
for (k in 1:K){
  Ns[[k]] = summ[[k]]$N
  sigmasq[[k]] = 1/Ns[[k]]
  Mt[[k]] = nrow(summ[[k]])
  tem[[k]] = sapply(1:Mt[[k]], function(x){which(snpinfo[,paste0('ind',k)]==x)})
  beta_init[[k]] = beta_inf[[k]][snpinfo[tem[[k]],paste0('ind',k)]]/summ[[k]][snpinfo$rsid[tem[[k]]],'c']
  Bh[[k]] = summ[[k]][snpinfo$rsid[tem[[k]]],paste0('Bh',k)]
  names(Bh[[k]]) = snpinfo$rsid[tem[[k]]]
}

# Define SNP indices: snp.indexabcde: SNPs that exist in all ancestry groups in (a,b,c,d,e):
indxinfo = snp.index.def(tem)
M = indxinfo$M; Mt = indxinfo$Mt; snp.index = indxinfo$snp.index; rm(indxinfo)

# Multiply genetic correlation by (-1) for SNP that have flipped alleles between ancestry groups 
rsigns = genetic.cor.signs(tem,snpinfo,ref_paths)
r.sign = rsigns$r.sign; rm(rsigns)

############ Step 1.4. Set up tuning parameter settings
## Load optimal tuning parameter from the single-ancestry analaysis
optim_ldpred2_params = lapply(1:K, function(x){fread2(ldpred2_params_path[x])})
p_est <- sapply(1:K, function(x){optim_ldpred2_params[[x]][,'p0']})
H2 <- sapply(1:K, function(x){abs(hsq[races[x]]) * optim_ldpred2_params[[x]][,'h20']})
sparse <- sapply(1:K, function(x){optim_ldpred2_params[[x]][,'sparse0']})

rs = tuning.params.setup(races, cors_additional)$rs
rs.dup = duplicated(rs)
if (sum(rs.dup) > 1) rs = rs[-duplicated(rs)]

settings0 = expand.grid(r.indx = 1:length(rs))
settings.p = as.data.frame(rbind(p_est, t(sapply(1:K, function(x){rep(p_est[x],K)}))))
if (sum(is.na(ps_additional)) == 0){
  temp.p = matrix(unlist(ps_additional),byrow=T,ncol = K)
  colnames(temp.p) = colnames(settings.p)
  settings.p = rbind(settings.p, temp.p)
}
colnames(settings.p) = paste0('p.causal',1:K); rownames(settings.p) = NULL
settings = do.call(rbind, lapply(1:nrow(settings.p), function(x){(cbind(settings0, settings.p[rep(x,nrow(settings0)),]))}))
if (sum(duplicated(settings)) > 0) settings = settings[-which(duplicated(settings)),]


indmat = as.matrix(snpinfo[,3:(K+2)])
snpinfo_rsid = snpinfo$rsid
n.burnin = 4e1 * K
niter = 1e2 * K

########################################################################
########################################################################

if ( opt$verbose >= 1 ) cat(paste0("\n** Step 2. Running MUSS on chromosome ", chr, ". ** \n"))
ncores = ifelse(NCORES >= nrow(settings), nrow(settings), NCORES)

registerDoMC(ncores)

# Run algorithm in parallel by chromosome

if (opt$verbose >= 1) cat(paste0('\n** MCMC started under all parameter settings. ** \n'))
ff <- foreach(ss = 1:nrow(settings), .combine='c', .multicombine=TRUE) %dopar% {
  # ii = icount(), .final = function(x) NULL
  r = rs[[settings[ss,'r.indx']]]
  p.causal = as.numeric(settings[ss,paste0('p.causal',1:K)]);
  
  MUSSout = MUSS(ss, chain=1, n.burnin, niter, settings, r, r.sign,
                    M, Mt, indmat, tem, beta_init, sigmasq, Bh, C.sfbm,
                    H2, snp.index, sparse, snpinfo)
  # progBar(ii, nrow(settings), per=1)
  list(MUSSout)
}

############
## Step 2.4. Combine estimated SNP effect sizes across all parameter settings and all ancestry groups
for (k in 1:K){
  prs.file = data.frame(SNP = summ[[k]]$SNP_ID, ALT = summ[[k]]$ALT)
  # Obtain estimated SNP effect sizes on the original scale
  for (ss in 1:nrow(settings)){
    tempbeta = ff[[ss]][[k]][rownames(summ[[k]])]*summ[[k]]$c
    tempbeta[is.na(tempbeta)] <- 0; tempbeta[tempbeta > 5] <- 0; tempbeta[tempbeta < -5] <- 0
    prs.file[,paste0('BETA',ss)] = tempbeta
  }
  write.table(prs.file,file = paste0(opt$PATH_out, '/tmp/MUSS_beta_in_all_settings_bychrom/', races[k], '-chr',chr,'.txt'), col.names = T,row.names = F,quote=F)
}
fwrite2(settings, paste0(opt$PATH_out,'/tmp/MUSS_beta_in_all_settings_bychrom/settings_',chr,".txt"), col.names = T, sep="\t", nThread=1)

if ( opt$verbose >= 1 ) cat(paste0('\n** Completed! Estimated SNP effect sizes and tuning parameter settings saved for chromosome ', chr, ' **\n'))


