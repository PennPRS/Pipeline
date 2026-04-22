library("optparse")
library(caret)
gbuild = 'hg19'
eth = c('EUR', 'AFR', 'EAS', 'AMR', 'SAS'); 
et = c('EUR', 'AFR', 'ASN', 'EUR', 'EUR')

for (race in eth[3:5]){
  a = NULL
  for (chr in 1:22){
    tem = bigreadr::fread2(paste0('/dcl01/chatterj/data/jzhang2/1000G/GRCh37/',race,'/chr',chr,'.bim'))
    tm = duplicated(tem[,2])
    if (sum(tm) > 0){
      a = c(a, unique(tem[tm,2]))
    }
  }
}


# Step 1: extract block information

# /tmp/byblock/
for (k in 1:length(eth)){
  race = eth[k];
  blockdir = paste0('/dcs04/nilanjan/data/jjin/prs/LDblocks/',race)
  if (!dir.exists(blockdir)){dir.create(blockdir)}
  temdir = paste0('/dcs04/nilanjan/data/jjin/prs/LDblocks/',race,'/tmp/')
  if (!dir.exists(temdir)){dir.create(temdir)}
  temdir = paste0('/dcs04/nilanjan/data/jjin/prs/LDblocks/',race,'/tmp/byblock/')
  if (!dir.exists(temdir)){dir.create(temdir)}
  ldblocks = bigreadr::fread2(paste0('/users/jjin/R/4.1.x/lassosum/data/Berisa.',et[k],'.',gbuild,'.bed')) # this LD block information is from the R package "lassosum".
  for (chr in 1:22){
    temdir = paste0('/dcs04/nilanjan/data/jjin/prs/LDblocks/',race,'/tmp/byblock/chr',chr,'/')
    if (!dir.exists(temdir)){dir.create(temdir)}
    ldblock = ldblocks[ldblocks$chr == paste0('chr',chr),]
    for (bl in 1:nrow(ldblock)){
      plinkcode = paste(
        "/dcl01/chatterj/data/jin/software/plink2",
        paste0('--bfile /dcs04/nilanjan/data/jjin/UKBval/',race,'/genotype/chr',chr),
        paste0('--from-bp ', ldblock[bl, 'start']),
        paste0('--to-bp ', ldblock[bl, 'stop']-1), 
        paste0('--chr ',chr),
        paste0("--make-bed"),
        paste0('--threads 1'),
        paste0('--out /dcs04/nilanjan/data/jjin/prs/LDblocks/',race,'/tmp/byblock/chr',chr,'/chr',chr,'_',ldblock[bl, 'start'],'_',ldblock[bl, 'stop'])
      )
      system(plinkcode)
    }
    plinkcode = paste(
      "/dcl01/chatterj/data/jin/software/plink2",
      paste0('--bfile /dcs04/nilanjan/data/jjin/UKBval/',race,'/genotype/chr',chr),
      paste0('--to-bp ', ldblock[1, 'start']),
      paste0('--chr ',chr),
      paste0("--make-bed"),
      paste0('--threads 1'),
      paste0('--out /dcs04/nilanjan/data/jjin/prs/LDblocks/',race,'/tmp/byblock/chr',chr,'/chr',chr,'_start_',ldblock[1, 'start'])
    )
    system(plinkcode)
    plinkcode = paste(
      "/dcl01/chatterj/data/jin/software/plink2",
      paste0('--bfile /dcs04/nilanjan/data/jjin/UKBval/',race,'/genotype/chr',chr),
      paste0('--from-bp ', ldblock[nrow(ldblock), 'stop']),
      paste0('--chr ',chr),
      paste0("--make-bed"),
      paste0('--threads 1'),
      paste0('--out /dcs04/nilanjan/data/jjin/prs/LDblocks/',race,'/tmp/byblock/chr',chr,'/chr',chr,'_',ldblock[nrow(ldblock), 'stop'],'_end')
    )
    system(plinkcode)
    print(paste0('Complete ',race,' chr ',chr))
  }
}


# /tmp/LD/:
for (k in 2){
  race = eth[k];
  blockdir = paste0('/dcs04/nilanjan/data/jjin/prs/LDblocks/',race)
  if (!dir.exists(blockdir)){dir.create(blockdir)}
  temdir = paste0('/dcs04/nilanjan/data/jjin/prs/LDblocks/',race,'/tmp/')
  if (!dir.exists(temdir)){dir.create(temdir)}
  temdir = paste0('/dcs04/nilanjan/data/jjin/prs/LDblocks/',race,'/tmp/LD/')
  if (!dir.exists(temdir)){dir.create(temdir)}
  ldblocks = bigreadr::fread2(paste0('/users/jjin/R/4.1.x/lassosum/data/Berisa.',et[k],'.',gbuild,'.bed'))
  for (chr in 1:22){
    temdir = paste0('/dcs04/nilanjan/data/jjin/prs/LDblocks/',race,'/tmp/LD/chr',chr,'/')
    if (!dir.exists(temdir)){dir.create(temdir)}
    ldblock = ldblocks[ldblocks$chr == paste0('chr',chr),]
    for (bl in 1:nrow(ldblock)){
      plinkcode = paste(
        # "/dcl01/chatterj/data/jin/software/plink",
        '/dcs04/nilanjan/data/jzhang2/MEPRS/PRSdir/plink1',
        paste0('--bfile /dcs04/nilanjan/data/jjin/prs/LDblocks/',race,'/tmp/byblock/chr',chr,'/chr',chr,'_',ldblock[bl, 'start'],'_',ldblock[bl, 'stop']),
        paste0('--keep-allele-order'),
        paste0('--r bin4'), 
        paste0('--threads 1'),
        paste0('--out /dcs04/nilanjan/data/jjin/prs/LDblocks/',race,'/tmp/LD/chr',chr,'/chr',chr,'_',ldblock[bl, 'start'],'_',ldblock[bl, 'stop'])
      )
      system(plinkcode)
    }
    plinkcode = paste(
      # "/dcl01/chatterj/data/jin/software/plink",
      '/dcs04/nilanjan/data/jzhang2/MEPRS/PRSdir/plink1',
      paste0('--bfile /dcs04/nilanjan/data/jjin/prs/LDblocks/',race,'/tmp/byblock/chr',chr,'/chr',chr,'_start_',ldblock[1, 'start']),
      paste0('--keep-allele-order'),
      paste0('--r bin4'), 
      paste0('--threads 1'),
      paste0('--out /dcs04/nilanjan/data/jjin/prs/LDblocks/',race,'/tmp/LD/chr',chr,'/chr',chr,'_start_',ldblock[1, 'start'])
    )
    system(plinkcode)
    plinkcode = paste(
      # "/dcl01/chatterj/data/jin/software/plink",
      '/dcs04/nilanjan/data/jzhang2/MEPRS/PRSdir/plink1',
      paste0('--bfile /dcs04/nilanjan/data/jjin/prs/LDblocks/',race,'/tmp/byblock/chr',chr,'/chr',chr,'_',ldblock[nrow(ldblock), 'stop'],'_end'),
      paste0('--keep-allele-order'),
      paste0('--r bin4'), 
      paste0('--threads 1'),
      paste0('--out /dcs04/nilanjan/data/jjin/prs/LDblocks/',race,'/tmp/LD/chr',chr,'/chr',chr,'_',ldblock[nrow(ldblock), 'stop'],'_end')
    )
    system(plinkcode)
    print(paste0('Complete ',race,' chr ',chr))
  }
}

# num_files <- (length(list.files(workdir))-2)/4
# log_files = list.files(workdir)
# log_files = log_files[which(str_detect(log_files, c('.log')))]
# log_list = list()
# for (l in 1:length(log_files)){
#   log_list[[l]] = bigreadr::fread2(paste0(workdir,log_files[l]))
#   rownames(log_list[[l]]) = log_list[[l]][,2]
# }
# pos_files = list.files(workdir)
# pos_files = pos_files[which(str_detect(pos_files, c('.log')))]
# pos_list = list()
# for (l in 1:length(pos_files)){
#   pos_list[[l]] = bigreadr::fread2(paste0(workdir,pos_files[l]))
#   rownames(pos_list[[l]]) = pos_list[[l]][,2]
# }


# ---------------------------- step 2 ----------------------------
# 2_reformat_to_RData_by_chr.R
suppressMessages(library("optparse"))
option_list = list(
  make_option("--packagedir", action="store", default=NA, type='character',
              help=" [required]"),
  make_option("--refgeno", action="store", default=NA, type='character',
              help=" [required]"),
  make_option("--hg", action="store", default=NA, type='character',
              help=" [required]"),
  make_option("--chr", action="store", default=NA, type='character',
              help=" [required]"),
  make_option("--workdir0", action="store", default=NA, type='character',
              help=" [required]"),
  make_option("--workdir", action="store", default=NA, type='character',
              help=" [required]")
)
opt = parse_args(OptionParser(option_list=option_list))


####################################################################################################
####################################################################################################

## reformat to standard data
library(caret)

for (k in 1:length(eth)){ # length(eth)
  race = eth[k];
  workdir = paste0('/dcs04/nilanjan/data/jjin/prs/LDblocks/',race,'/standard_data')
  if (!dir.exists(workdir)) dir.create(paste0(workdir))
  block_info <- read.table(paste0('/users/jjin/R/4.1.x/lassosum/data/Berisa.',et[k],'.',gbuild,'.bed'),header = T)
  workdir0 = paste0("/dcs04/nilanjan/data/jjin/prs/LDblocks/", race)
  for (chr in 1:22){
    block_info_tmp <- block_info[block_info$chr== paste0('chr',chr),]
    
    Nsnps <- integer(length = nrow(block_info_tmp)+2)
    snps_list <- vector("list", length = nrow(block_info_tmp)+2)
    LD_list <- vector("list", length = nrow(block_info_tmp)+2)
    
    #### start
    snps <- character()
    tmpfile = paste0(workdir0, "/tmp/byblock/chr", chr,"/chr", chr,"_start_",block_info_tmp$start[1],".bim")
    if (file.exists(tmpfile)){
      tmp.snps <- try(read.table(tmpfile, stringsAsFactors = F), silent=TRUE)
      if ('try-error' %in% class(tmp.snps)) {
        Nsnps[1] <- 0
      }else{
        n.snp.tmp <- nrow(tmp.snps)
        tmp.LD <- readBin(paste0(workdir0, "/tmp/LD/chr", chr,"/chr",
                                 chr,"_start_", block_info_tmp$start[1],".ld.bin"),
                          what="numeric", size=4, n=(n.snp.tmp)^2)
        cat(paste0("Total #of SNP is ",n.snp.tmp,", and ",sum(is.nan(tmp.LD))," is nan.\n"))
        tmp.LD[is.nan(tmp.LD)] <- 1; tmp.LD <- matrix(tmp.LD, ncol=n.snp.tmp)
        drop = findCorrelation(tmp.LD,cutoff = 0.999999)
        
        Nsnps[1] <- n.snp.tmp - length(drop)
        snps_list[[1]] <- tmp.snps$V2[-drop]
        LD_list[[1]] <- tmp.LD[-drop, -drop]
      }
    }
    #### Median
    for (i in 1:nrow(block_info_tmp)){
      snps <- character()
      tmp.snps <- try(read.table(paste0(workdir0, "/tmp/byblock/chr", chr,"/chr",
                                        chr,"_",block_info_tmp$start[i],"_",block_info_tmp$stop[i],".bim"), stringsAsFactors = F), silent=TRUE)
      if ('try-error' %in% class(tmp.snps)) {
        Nsnps[i+1] <- 0
      }else{
        n.snp.tmp <- nrow(tmp.snps)
        tmp.LD <- readBin(paste0(workdir0, "/tmp/LD/chr", chr,"/chr",
                                 chr,"_", block_info_tmp$start[i],"_",block_info_tmp$stop[i],".ld.bin"),
                          what="numeric", size=4, n=(n.snp.tmp)^2)
        print(paste0("Total #of SNP is ",n.snp.tmp,", and ",sum(is.nan(tmp.LD))," is nan.\n"))
        tmp.LD[is.nan(tmp.LD)] <- 1; tmp.LD <- matrix(tmp.LD, ncol=n.snp.tmp)
        drop = findCorrelation(tmp.LD,cutoff = 0.999999) # 0.999999
        
        Nsnps[i+1] <- n.snp.tmp - length(drop)
        snps_list[[i+1]] <- tmp.snps$V2[-drop]
        LD_list[[i+1]] <- tmp.LD[-drop, -drop]
      }
      cat(paste0(chr,": ",i,"/",nrow(block_info_tmp),"\n"))
    }
    #warnings()
    
    #### end
    snps <- character()
    tmpfile = paste0(workdir0, "/tmp/byblock/chr", chr,"/chr", chr,"_",block_info_tmp$stop[i],"_end.bim")
    if (file.exists(tmpfile)){
      tmp.snps <- try(read.table(tmpfile, stringsAsFactors = F), silent=TRUE)
      if ('try-error' %in% class(tmp.snps)) {
        Nsnps[i+2] <- 0
      }else{
        n.snp.tmp <- nrow(tmp.snps)
        tmp.LD <- readBin(paste0(workdir0, "/tmp/LD/chr", chr,"/chr",
                                 chr,"_",block_info_tmp$stop[i],"_end.ld.bin"),
                          what="numeric", size=4, n=(n.snp.tmp)^2)
        cat(paste0("Total #of SNP is ",n.snp.tmp,", and ",sum(is.nan(tmp.LD))," is nan.\n"))
        tmp.LD[is.nan(tmp.LD)] <- 1; tmp.LD <- matrix(tmp.LD, ncol=n.snp.tmp)
        drop = findCorrelation(tmp.LD,cutoff = 0.999999)
        
        Nsnps[i+2] <- n.snp.tmp - length(drop)
        snps_list[[i+2]] <- tmp.snps$V2[-drop]
        LD_list[[i+2]] <- tmp.LD[-drop, -drop]
      }
    }
    
    cat(paste0("Saving standard data for ",chr ,"...\n"))
    save(Nsnps, snps_list,
         file = paste0(workdir,"/chr",chr,"_snps.RData"))
    
    save(Nsnps, snps_list, LD_list,
         file = paste0(workdir,"/chr",chr,"_LD.RData"))
    
    cat(paste0(chr, " completed.\n"))
  }
}





