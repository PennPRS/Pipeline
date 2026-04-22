library("optparse")
library(caret)
# First, download the R package lassosum, directory to the package:
dir.lassosum = '/users/jjin/R/4.1.x/lassosum/'
eth = c('EUR','AFR','AMR','EAS','SAS') # race
et = c('EUR', 'AFR', 'EUR', 'ASN', 'EUR') # corresponding file names in the lassosum package.
gbuild = 'hg19'
dir.plink2 = '/dcl01/chatterj/data/jin/software/plink2' # path to the PLINK 2.0 software
dir.plink = '/dcl01/chatterj/data/jin/software/plink' # path to the PLINK 1.9 software

# directory where you want to save the generated LD matrix info by LD block
LDdir = '/dcs04/nilanjan/data/jjin/prs/testLDblocks/' 

# Split genotype data files by LD block
# /tmp/byblock/:
for (k in 1:length(eth)){ # Generate LD info for each ancestry k
  race = eth[k];
  blockdir = paste0(LDdir,race) 
  if (!dir.exists(blockdir)){dir.create(blockdir)}
  temdir = paste0(blockdir, '/tmp/')
  if (!dir.exists(temdir)){dir.create(temdir)}
  temdir = paste0(blockdir, '/tmp/byblock/')
  if (!dir.exists(temdir)){dir.create(temdir)}
  
  ldblocks = bigreadr::fread2(paste0(dir.lassosum, 'data/Berisa.',et[k],'.',gbuild,'.bed')) # load LD block position information
  refdir = paste0('/dcs04/nilanjan/data/jjin/prs/1KGref_MEGA/GRCh37/',race) # directory to the reference genotype data 
  for (chr in 1:22){
    temdir = paste0(blockdir,'/tmp/byblock/chr',chr,'/')
    if (!dir.exists(temdir)){dir.create(temdir)}
    ldblock = ldblocks[ldblocks$chr == paste0('chr',chr),]
    for (bl in 1:nrow(ldblock)){
      plinkcode = paste(
        "/dcl01/chatterj/data/jin/software/plink2",
        paste0('--bfile ', refdir, '/chr',chr),
        paste0('--from-bp ', ldblock[bl, 'start']),
        paste0('--to-bp ', ldblock[bl, 'stop']-1),
        paste0('--chr ',chr),
        paste0("--make-bed"),
        paste0('--threads 1'),
        paste0('--out ', blockdir, '/tmp/byblock/chr',chr,'/chr',chr,'_',ldblock[bl, 'start'],'_',ldblock[bl, 'stop'])
      )
      system(plinkcode)
    }
    plinkcode = paste(
      "/dcl01/chatterj/data/jin/software/plink2",
      paste0('--bfile ', refdir, '/chr',chr),
      paste0('--to-bp ', ldblock[1, 'start']),
      paste0('--chr ',chr),
      paste0("--make-bed"),
      paste0('--threads 1'),
      paste0('--out ', blockdir, '/tmp/byblock/chr',chr,'/chr',chr,'_start_',ldblock[1, 'start'])
    )
    system(plinkcode)
    plinkcode = paste(
      "/dcl01/chatterj/data/jin/software/plink2",
      paste0('--bfile ', refdir, '/chr',chr),
      paste0('--from-bp ', ldblock[nrow(ldblock), 'stop']),
      paste0('--chr ',chr),
      paste0("--make-bed"),
      paste0('--threads 1'),
      paste0('--out ', blockdir, '/tmp/byblock/chr',chr,'/chr',chr,'_',ldblock[nrow(ldblock), 'stop'],'_end')
    )
    system(plinkcode)
    print(paste0('Genotype files split by LD block for ',race,' chr ',chr))
  }
}


# Generate LD information by LD block
# /tmp/LD/:
for (k in 1:length(eth)){
  race = eth[k];
  blockdir = paste0(LDdir,race) 
  if (!dir.exists(blockdir)){dir.create(blockdir)}
  temdir = paste0(blockdir, '/tmp/LD/')
  if (!dir.exists(temdir)){dir.create(temdir)}
  ldblocks = bigreadr::fread2(paste0(dir.lassosum, 'data/Berisa.',et[k],'.',gbuild,'.bed')) # load LD block position information
  for (chr in 1:22){
    temdir = paste0(blockdir,'/tmp/LD/chr',chr,'/')
    if (!dir.exists(temdir)){dir.create(temdir)}
    ldblock = ldblocks[ldblocks$chr == paste0('chr',chr),]
    for (bl in 1:nrow(ldblock)){
      plinkcode = paste(
        dir.plink,
        paste0('--bfile ', blockdir, '/tmp/byblock/chr',chr,'/chr',chr,'_',ldblock[bl, 'start'],'_',ldblock[bl, 'stop']),
        paste0('--keep-allele-order'),
        paste0('--r bin4'),
        paste0('--threads 1'),
        paste0('--out ', blockdir, '/tmp/LD/chr',chr,'/chr',chr,'_',ldblock[bl, 'start'],'_',ldblock[bl, 'stop'])
      )
      system(plinkcode)
    }
    plinkcode = paste(
      dir.plink,
      paste0('--bfile ', blockdir, '/tmp/byblock/chr',chr,'/chr',chr,'_start_',ldblock[1, 'start']),
      paste0('--keep-allele-order'),
      paste0('--r bin4'),
      paste0('--threads 1'),
      paste0('--out ', blockdir, '/tmp/LD/chr',chr,'/chr',chr,'_start_',ldblock[1, 'start'])
    ) # this step may print error message like "Error: No variants remaining after main filters.", which is fine.
    system(plinkcode)
    plinkcode = paste(
      dir.plink,
      paste0('--bfile ', blockdir, '/tmp/byblock/chr',chr,'/chr',chr,'_',ldblock[nrow(ldblock), 'stop'],'_end'),
      paste0('--keep-allele-order'),
      paste0('--r bin4'),
      paste0('--threads 1'),
      paste0('--out ', blockdir, '/tmp/LD/chr',chr,'/chr',chr,'_',ldblock[nrow(ldblock), 'stop'],'_end')
    ) # this step may print error message like "Error: No variants remaining after main filters.", which is fine.
    system(plinkcode)
    print(paste0('Complete ',race,' chr ',chr))
  }
}





####################################################################################################
####################################################################################################

## Reformat to standard data
library(caret)

# /standard_data:
for (k in 1:length(eth)){
  race = eth[k];
  workdir = paste0(LDdir,race,'/standard_data')
  if (!dir.exists(workdir)) dir.create(paste0(workdir))
  block_info <- read.table(paste0(dir.lassosum, 'data/Berisa.',et[k],'.',gbuild,'.bed'),header = T)
  blockdir = paste0(LDdir, race)
  for (chr in 1:22){
    block_info_tmp <- block_info[block_info$chr== paste0('chr',chr),]
    
    Nsnps <- integer(length = nrow(block_info_tmp)+2)
    snps_list <- vector("list", length = nrow(block_info_tmp)+2)
    LD_list <- vector("list", length = nrow(block_info_tmp)+2)
    
    #### start
    snps <- character()
    tmpfile = paste0(blockdir, "/tmp/byblock/chr", chr,"/chr", chr,"_start_", block_info_tmp$start[1],".bim")
    if (file.exists(tmpfile)){
      tmp.snps <- try(read.table(tmpfile, stringsAsFactors = F), silent=TRUE)
      if ('try-error' %in% class(tmp.snps)) {
        Nsnps[1] <- 0
      }else{
        n.snp.tmp <- nrow(tmp.snps)
        tmp.LD <- readBin(paste0(blockdir, "/tmp/LD/chr", chr,"/chr",
                                 chr,"_start_", block_info_tmp$start[1],".ld.bin"),
                          what="numeric", size=4, n=(n.snp.tmp)^2)
        cat(paste0("Total #of SNP is ",n.snp.tmp,", and ",sum(is.nan(tmp.LD))," is nan.\n"))
        tmp.LD[is.nan(tmp.LD)] <- 1; tmp.LD <- matrix(tmp.LD, ncol=n.snp.tmp)
        if (n.snp.tmp > 1){
          drop = findCorrelation(tmp.LD,cutoff = 0.999999)
          Nsnps[1] <- n.snp.tmp - length(drop)
          snps_list[[1]] <- tmp.snps$V2[-drop]
          LD_list[[1]] <- tmp.LD[-drop, -drop]
        } 
      }
    }
    #### Median
    for (i in 1:nrow(block_info_tmp)){
      snps <- character()
      tmp.snps <- try(read.table(paste0(blockdir, "/tmp/byblock/chr", chr,"/chr",
                                        chr, "_", block_info_tmp$start[i], "_", block_info_tmp$stop[i],".bim"), stringsAsFactors = F), silent=TRUE)
      if ('try-error' %in% class(tmp.snps)) {
        Nsnps[i+1] <- 0
      }else{
        n.snp.tmp <- nrow(tmp.snps)
        tmp.LD <- readBin(paste0(blockdir, "/tmp/LD/chr", chr,"/chr",
                                 chr,"_", block_info_tmp$start[i], "_", block_info_tmp$stop[i],".ld.bin"),
                          what="numeric", size=4, n=(n.snp.tmp)^2)
        print(paste0("Total #of SNP is ",n.snp.tmp,", and ",sum(is.nan(tmp.LD))," is nan.\n"))
        tmp.LD[is.nan(tmp.LD)] <- 1; tmp.LD <- matrix(tmp.LD, ncol=n.snp.tmp)
        if (n.snp.tmp > 1){
          drop = findCorrelation(tmp.LD,cutoff = 0.999999)
          Nsnps[i+1] <- n.snp.tmp - length(drop)
          snps_list[[i+1]] <- tmp.snps$V2[-drop]
          LD_list[[i+1]] <- tmp.LD[-drop, -drop]
        } 
      }
      cat(paste0(chr,": ",i,"/",nrow(block_info_tmp),"\n"))
    }
    #warnings()
    
    #### end
    snps <- character()
    tmpfile = paste0(blockdir, "/tmp/byblock/chr", chr, "/chr", chr, "_", block_info_tmp$stop[i],"_end.bim")
    if (file.exists(tmpfile)){
      tmp.snps <- try(read.table(tmpfile, stringsAsFactors = F), silent=TRUE)
      if ('try-error' %in% class(tmp.snps)) {
        Nsnps[i+2] <- 0
      }else{
        n.snp.tmp <- nrow(tmp.snps)
        tmp.LD <- readBin(paste0(blockdir, "/tmp/LD/chr", chr, "/chr",
                                 chr, "_", block_info_tmp$stop[i], "_end.ld.bin"),
                          what="numeric", size=4, n=(n.snp.tmp)^2)
        cat(paste0("Total #of SNP is ",n.snp.tmp,", and ",sum(is.nan(tmp.LD))," is nan.\n"))
        tmp.LD[is.nan(tmp.LD)] <- 1; tmp.LD <- matrix(tmp.LD, ncol=n.snp.tmp)
        if (n.snp.tmp > 1){
          drop = findCorrelation(tmp.LD,cutoff = 0.999999)
          Nsnps[i+2] <- n.snp.tmp - length(drop)
          snps_list[[i+2]] <- tmp.snps$V2[-drop]
          LD_list[[i+2]] <- tmp.LD[-drop, -drop]
        } 
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

