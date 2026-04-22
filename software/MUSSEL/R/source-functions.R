unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

## progress bar function: https://stackoverflow.com/questions/51213293/is-it-possible-to-get-a-progress-bar-with-foreach-and-a-multicore-kind-of-back
progBar <- function(ii, N, per = 10) {
  # ii: index of the current iteration.
  # N: total number of iterations.
  # per: per how may iterations do we want to update the progress bar. We need this as when multiple iterations are being performed progress bar gets updated too often and output is messed up.
  if (ii %in% seq(1, N, per)) {
    x <- round(ii * 100 / N)
    message("[ ",
            paste(rep("-", x), collapse = ""),
            paste(rep(" ", 100 - x), collapse = ""),
            " ] ", x, "%", "\r", 
            appendLF = FALSE)
    if (ii == N) cat("\r")
  }
}

combine.alleles = function(x) paste(x,collapse='')

snp.index.def = function(tem){
  K = length(tem)
  if (K == 2){
    snp.index1 = setdiff(tem[[1]],Reduce(union, tem[-1]))
    snp.index2 = setdiff(tem[[2]],Reduce(union, tem[-2]))
    
    snp.index12 = setdiff(Reduce(intersect, tem[c(1,2)]),Reduce(union, tem[-c(1,2)]))
    snp.index1pool = tem[[1]]; snp.index2pool = tem[[2]]
    snp.indexpool = unique(unlist(tem))
    
    M1 = length(snp.index1); M2 = length(snp.index2); M12 = length(snp.index12)
    M1t = length(snp.index1pool); M2t = length(snp.index2pool)
    M = nrow(snpinfo); Mt = c(M1t, M2t)
    snp.index = list(snp.index1 = snp.index1, snp.index2 = snp.index2, snp.index12 = snp.index12)
  }
  if (K == 3){
    snp.index1 = setdiff(tem[[1]],Reduce(union, tem[-1]))
    snp.index2 = setdiff(tem[[2]],Reduce(union, tem[-2]))
    snp.index3 = setdiff(tem[[3]],Reduce(union, tem[-3]))
    
    snp.index12 = setdiff(Reduce(intersect, tem[c(1,2)]),Reduce(union, tem[-c(1,2)]))
    snp.index13 = setdiff(Reduce(intersect, tem[c(1,3)]),Reduce(union, tem[-c(1,3)]))
    snp.index23 = setdiff(Reduce(intersect, tem[c(2,3)]),Reduce(union, tem[-c(2,3)]))
    snp.index123 = setdiff(Reduce(intersect, tem[c(1,2,3)]),Reduce(union, tem[-c(1,2,3)]))
    
    snp.index1pool = tem[[1]]; snp.index2pool = tem[[2]]; snp.index3pool = tem[[3]]
    snp.indexpool = unique(unlist(tem))
    
    M1 = length(snp.index1); M2 = length(snp.index2); M3 = length(snp.index3)
    M12 = length(snp.index12); M13 = length(snp.index13); M23 = length(snp.index23);
    M123 = length(snp.index123); 
    M1t = length(snp.index1pool); M2t = length(snp.index2pool); M3t = length(snp.index3pool)
    M = nrow(snpinfo); Mt = c(M1t, M2t, M3t)
    snp.index = list(snp.index1 = snp.index1, snp.index2 = snp.index2, snp.index3 = snp.index3, 
                     snp.index12 = snp.index12, snp.index13 = snp.index13,
                     snp.index23 = snp.index23, 
                     snp.index123 = snp.index123)
  }
  if (K == 4){
    snp.index1 = setdiff(tem[[1]],Reduce(union, tem[-1]))
    snp.index2 = setdiff(tem[[2]],Reduce(union, tem[-2]))
    snp.index3 = setdiff(tem[[3]],Reduce(union, tem[-3]))
    snp.index4 = setdiff(tem[[4]],Reduce(union, tem[-4]))
    
    snp.index12 = setdiff(Reduce(intersect, tem[c(1,2)]),Reduce(union, tem[-c(1,2)]))
    snp.index13 = setdiff(Reduce(intersect, tem[c(1,3)]),Reduce(union, tem[-c(1,3)]))
    snp.index14 = setdiff(Reduce(intersect, tem[c(1,4)]),Reduce(union, tem[-c(1,4)]))
    snp.index23 = setdiff(Reduce(intersect, tem[c(2,3)]),Reduce(union, tem[-c(2,3)]))
    snp.index24 = setdiff(Reduce(intersect, tem[c(2,4)]),Reduce(union, tem[-c(2,4)]))
    snp.index34 = setdiff(Reduce(intersect, tem[c(3,4)]),Reduce(union, tem[-c(3,4)]))
    
    snp.index123 = setdiff(Reduce(intersect, tem[c(1,2,3)]),Reduce(union, tem[-c(1,2,3)]))
    snp.index124 = setdiff(Reduce(intersect, tem[c(1,2,4)]),Reduce(union, tem[-c(1,2,4)]))
    snp.index134 = setdiff(Reduce(intersect, tem[c(1,3,4)]),Reduce(union, tem[-c(1,3,4)]))
    snp.index234 = setdiff(Reduce(intersect, tem[c(2,3,4)]),Reduce(union, tem[-c(2,3,4)]))
    
    snp.index1234 = setdiff(Reduce(intersect, tem[c(1,2,3,4)]),Reduce(union, tem[-c(1,2,3,4)]))
    
    snp.index1pool = tem[[1]]; snp.index2pool = tem[[2]]; snp.index3pool = tem[[3]]; snp.index4pool = tem[[4]]
    snp.indexpool = unique(unlist(tem))
    
    M1 = length(snp.index1); M2 = length(snp.index2); M3 = length(snp.index3); M4 = length(snp.index4)
    M12 = length(snp.index12); M13 = length(snp.index13); M14 = length(snp.index14)
    M23 = length(snp.index23); M24 = length(snp.index24); M34 = length(snp.index34);
    
    M123 = length(snp.index123); M124 = length(snp.index124); M134 = length(snp.index134); M234 = length(snp.index234); M1234 = length(snp.index1234)
    M1t = length(snp.index1pool); M2t = length(snp.index2pool); M3t = length(snp.index3pool); M4t = length(snp.index4pool)
    M = nrow(snpinfo); Mt = c(M1t, M2t, M3t, M4t)
    snp.index = list(snp.index1 = snp.index1, snp.index2 = snp.index2, snp.index3 = snp.index3, snp.index4 = snp.index4,
                     snp.index12 = snp.index12, snp.index13 = snp.index13, snp.index14 = snp.index14,
                     snp.index23 = snp.index23, snp.index24 = snp.index24, snp.index34 = snp.index34,
                     snp.index123 = snp.index123, snp.index124 = snp.index124, snp.index134 = snp.index134, 
                     snp.index234 = snp.index234, snp.index1234 = snp.index1234)
    
  }
  if (K == 5){
    snp.index1 = setdiff(tem[[1]],Reduce(union, tem[-1]))
    snp.index2 = setdiff(tem[[2]],Reduce(union, tem[-2]))
    snp.index3 = setdiff(tem[[3]],Reduce(union, tem[-3]))
    snp.index4 = setdiff(tem[[4]],Reduce(union, tem[-4]))
    snp.index5 = setdiff(tem[[5]],Reduce(union, tem[-5]))
    
    snp.index12 = setdiff(Reduce(intersect, tem[c(1,2)]),Reduce(union, tem[-c(1,2)]))
    snp.index13 = setdiff(Reduce(intersect, tem[c(1,3)]),Reduce(union, tem[-c(1,3)]))
    snp.index14 = setdiff(Reduce(intersect, tem[c(1,4)]),Reduce(union, tem[-c(1,4)]))
    snp.index15 = setdiff(Reduce(intersect, tem[c(1,5)]),Reduce(union, tem[-c(1,5)]))
    snp.index23 = setdiff(Reduce(intersect, tem[c(2,3)]),Reduce(union, tem[-c(2,3)]))
    snp.index24 = setdiff(Reduce(intersect, tem[c(2,4)]),Reduce(union, tem[-c(2,4)]))
    snp.index25 = setdiff(Reduce(intersect, tem[c(2,5)]),Reduce(union, tem[-c(2,5)]))
    snp.index34 = setdiff(Reduce(intersect, tem[c(3,4)]),Reduce(union, tem[-c(3,4)]))
    snp.index35 = setdiff(Reduce(intersect, tem[c(3,5)]),Reduce(union, tem[-c(3,5)]))
    snp.index45 = setdiff(Reduce(intersect, tem[c(4,5)]),Reduce(union, tem[-c(4,5)]))
    snp.index123 = setdiff(Reduce(intersect, tem[c(1,2,3)]),Reduce(union, tem[-c(1,2,3)]))
    snp.index124 = setdiff(Reduce(intersect, tem[c(1,2,4)]),Reduce(union, tem[-c(1,2,4)]))
    snp.index125 = setdiff(Reduce(intersect, tem[c(1,2,5)]),Reduce(union, tem[-c(1,2,5)]))
    snp.index134 = setdiff(Reduce(intersect, tem[c(1,3,4)]),Reduce(union, tem[-c(1,3,4)]))
    snp.index135 = setdiff(Reduce(intersect, tem[c(1,3,5)]),Reduce(union, tem[-c(1,3,5)]))
    snp.index145 = setdiff(Reduce(intersect, tem[c(1,4,5)]),Reduce(union, tem[-c(1,4,5)]))
    snp.index234 = setdiff(Reduce(intersect, tem[c(2,3,4)]),Reduce(union, tem[-c(2,3,4)]))
    snp.index235 = setdiff(Reduce(intersect, tem[c(2,3,5)]),Reduce(union, tem[-c(2,3,5)]))
    snp.index245 = setdiff(Reduce(intersect, tem[c(2,4,5)]),Reduce(union, tem[-c(2,4,5)]))
    snp.index345 = setdiff(Reduce(intersect, tem[c(3,4,5)]),Reduce(union, tem[-c(3,4,5)]))
    
    snp.index1234 = setdiff(Reduce(intersect, tem[c(1,2,3,4)]),Reduce(union, tem[-c(1,2,3,4)]))
    snp.index1235 = setdiff(Reduce(intersect, tem[c(1,2,3,5)]),Reduce(union, tem[-c(1,2,3,5)]))
    snp.index1245 = setdiff(Reduce(intersect, tem[c(1,2,4,5)]),Reduce(union, tem[-c(1,2,4,5)]))
    snp.index1345 = setdiff(Reduce(intersect, tem[c(1,3,4,5)]),Reduce(union, tem[-c(1,3,4,5)]))
    snp.index2345 = setdiff(Reduce(intersect, tem[c(2,3,4,5)]),Reduce(union, tem[-c(2,3,4,5)]))
    
    snp.index12345 = Reduce(intersect, tem[c(1,2,3,4,5)])
    
    snp.index1pool = tem[[1]]; snp.index2pool = tem[[2]]; snp.index3pool = tem[[3]]; snp.index4pool = tem[[4]]; snp.index5pool = tem[[5]]
    snp.indexpool = unique(unlist(tem))
    
    M1 = length(snp.index1); M2 = length(snp.index2); M3 = length(snp.index3); M4 = length(snp.index4); M5 = length(snp.index5)
    M12 = length(snp.index12); M13 = length(snp.index13); M14 = length(snp.index14); M15 = length(snp.index15)
    M23 = length(snp.index23); M24 = length(snp.index24); M25 = length(snp.index25); 
    M34 = length(snp.index34); M35 = length(snp.index35); M45 = length(snp.index45);
    
    M123 = length(snp.index123); M124 = length(snp.index124); M125 = length(snp.index125); 
    M134 = length(snp.index134); M135 = length(snp.index135); M145 = length(snp.index145); 
    M234 = length(snp.index234); M235 = length(snp.index235); M245 = length(snp.index245); M345 = length(snp.index345);
    M1234 = length(snp.index1234); M1235 = length(snp.index1235); M1245 = length(snp.index1245);
    M1345 = length(snp.index1345); M2345 = length(snp.index2345);
    M1t = length(snp.index1pool); M2t = length(snp.index2pool); M3t = length(snp.index3pool); M4t = length(snp.index4pool); M5t = length(snp.index5pool)
    M = nrow(snpinfo); Mt = c(M1t, M2t, M3t, M4t, M5t)
    snp.index = list(snp.index1 = snp.index1, snp.index2 = snp.index2, snp.index3 = snp.index3, snp.index4 = snp.index4, snp.index5 = snp.index5,
                     snp.index12 = snp.index12, snp.index13 = snp.index13, snp.index14 = snp.index14, snp.index15 = snp.index15,
                     snp.index23 = snp.index23, snp.index24 = snp.index24, snp.index25 = snp.index25,
                     snp.index34 = snp.index34, snp.index35 = snp.index35, snp.index45 = snp.index45,
                     snp.index123 = snp.index123, snp.index124 = snp.index124, snp.index125 = snp.index125,
                     snp.index134 = snp.index134, snp.index135 = snp.index135, snp.index145 = snp.index145,
                     snp.index234 = snp.index234, snp.index235 = snp.index235, snp.index245 = snp.index245, snp.index345 = snp.index345,
                     snp.index1234 = snp.index1234, snp.index1235 = snp.index1235, snp.index1245 = snp.index1245,
                     snp.index1345 = snp.index1345, snp.index2345 = snp.index2345,
                     snp.index12345 = snp.index12345)
  }
  return(list(M = M, Mt = Mt, snp.index = snp.index))
}

genetic.cor.signs = function(tem,snpinfo,ref_paths){
  K = length(tem)
  for (k1 in 1:(K-1)){
    for (k2 in (k1+1):K){
      snpindx = intersect(tem[[k1]],tem[[k2]])
      sharedsnp = snpinfo$rsid[snpindx]
      bim1 = bigreadr::fread2(paste0(ref_paths[k1],'/chr',chr,'.bim'),header=F)
      bim1 = bim1 %>% filter(!duplicated(V2))
      colnames(bim1) = c('chrom','sid','na','pos','nt1','nt2')
      bim1 = bim1 %>% filter(sid %in% sharedsnp)
      bim1$sid = as.character(bim1$sid)
      bim1 = bim1[,c('sid','nt1','nt2')]
      colnames(bim1)[1] = 'SNP_ID'
      
      bim2 = bigreadr::fread2(paste0(ref_paths[k2],'/chr',chr,'.bim'),header=F)
      bim2 = bim2 %>% filter(!duplicated(V2))
      colnames(bim2) = c('chrom','sid','na','pos','nt1','nt2')
      bim2 = bim2 %>% filter(sid %in% sharedsnp)
      bim2$sid = as.character(bim2$sid)
      bim2 = bim2[,c('sid','nt1','nt2')]
      colnames(bim2)[1] = 'SNP_ID'
      
      ###
      bim = merge(bim1, bim2, by='SNP_ID')
      rownames(bim) = bim$SNP_ID
      bim = bim[snpinfo$rsid[snpindx],]
      
      alleles = cbind(bim$nt1.x, bim$nt2.x, bim$nt1.y, bim$nt2.y)
      combine.alleles = function(x) paste(x,collapse='')
      alleles = apply(alleles,1,combine.alleles)
      inds.ambiguous= which(alleles %in% c('ACGT','AGCT','TCGA','TGCA','CATG','CTAG','GATC','GTAC'))
      inds.ambiguous.keep = which(alleles %in% c('ACTG','AGTC','TCAG','TGAC','CAGT','CTGA','GACT','GTCA'))
      matched = which((bim$nt1.x == bim$nt1.y)&(bim$nt2.x == bim$nt2.y))
      flipped = which((bim$nt1.x == bim$nt2.y)&(bim$nt2.x == bim$nt1.y))
      
      if ((k1 == 1)&(k2 == 2)){
        r12.sign = rep(0,M)
        r12.sign[snpindx[matched]] = 1
        r12.sign[snpindx[flipped]] = -1
      }
      if ((k1 == 1)&(k2 == 3)){
        r13.sign = rep(0,M)
        r13.sign[snpindx[matched]] = 1
        r13.sign[snpindx[flipped]] = -1
      }
      if ((k1 == 1)&(k2 == 4)){
        r14.sign = rep(0,M)
        r14.sign[snpindx[matched]] = 1
        r14.sign[snpindx[flipped]] = -1
      }
      if ((k1 == 1)&(k2 == 5)){
        r15.sign = rep(0,M)
        r15.sign[snpindx[matched]] = 1
        r15.sign[snpindx[flipped]] = -1
      }
      if ((k1 == 2)&(k2 == 3)){
        r23.sign = rep(0,M)
        r23.sign[snpindx[matched]] = 1
        r23.sign[snpindx[flipped]] = -1
      }
      if ((k1 == 2)&(k2 == 4)){
        r24.sign = rep(0,M)
        r24.sign[snpindx[matched]] = 1
        r24.sign[snpindx[flipped]] = -1
      }
      if ((k1 == 2)&(k2 == 5)){
        r25.sign = rep(0,M)
        r25.sign[snpindx[matched]] = 1
        r25.sign[snpindx[flipped]] = -1
      }
      if ((k1 == 3)&(k2 == 4)){
        r34.sign = rep(0,M)
        r34.sign[snpindx[matched]] = 1
        r34.sign[snpindx[flipped]] = -1
      }
      if ((k1 == 3)&(k2 == 5)){
        r35.sign = rep(0,M)
        r35.sign[snpindx[matched]] = 1
        r35.sign[snpindx[flipped]] = -1
      }
      if ((k1 == 4)&(k2 == 5)){
        r45.sign = rep(0,M)
        r45.sign[snpindx[matched]] = 1
        r45.sign[snpindx[flipped]] = -1
      }
    }
  }
  if (K == 2){
    r.sign = list(r12.sign = r12.sign)
  }
  if (K == 3){
    r.sign = list(r12.sign = r12.sign, r13.sign = r13.sign, r23.sign = r23.sign)
  }
  if (K == 4){
    r.sign = list(r12.sign = r12.sign, r13.sign = r13.sign, r14.sign = r14.sign, 
                  r23.sign = r23.sign, r24.sign = r24.sign, r34.sign = r34.sign)
  }
  if (K == 5){
    r.sign = list(r12.sign = r12.sign, r13.sign = r13.sign, r14.sign = r14.sign, r15.sign = r15.sign,
                  r23.sign = r23.sign, r24.sign = r24.sign, r25.sign = r25.sign, 
                  r34.sign = r34.sign, r35.sign = r35.sign, r45.sign = r45.sign)
  }
  return(list(r.sign = r.sign))
}

tuning.params.setup = function(races, cors_additional){
  cor_with_afr = 0.75
  cor_wo_afr = 0.9
  K = length(races)
  
  r12 = c(0.7,0.8,0.9,0.95)
  rs = lapply(1:length(r12), function(x){rep(r12[x],K*(K-1)/2)})
  rs.unequal = rep(cor_wo_afr,K*(K-1)/2)
  pair.indx = matrix(unlist(sapply(1:(K-1), function(x){as.numeric(sapply((x+1):K, function(y){c(x,y)}))})),byrow = T,ncol = 2)
  afr.indx = sapply(1:nrow(pair.indx), function(x){2 %in% pair.indx[x,]})
  rs.unequal[afr.indx] = cor_with_afr
  rs[[length(rs) + 1]] = rs.unequal
  if (sum(is.na(cors_additional)) == 0) rs[(length(rs)+1):(length(rs)+length(cors_additional))] = cors_additional
  return(list(rs = rs))
}


MUSS = function(ss, chain, n.burnin, niter, settings, r, r.sign, 
                    M, Mt, indmat, tem, beta_init, sigmasq, Bh, C.sfbm,
                    H2, snp.index, sparse, snpinfo){
  K = ncol(indmat)
  set.seed(2020)
  if (K == 2){
    sigmasq1 = sigmasq[[1]]; sigmasq2 = sigmasq[[2]]
    p.causal = settings[ss,paste0('p.causal',1:K)]; p.causal = as.numeric(p.causal)
    p.causal1 = p.causal[1]; p.causal2 = p.causal[2]
    h2 = numeric()
    for (k in 1:K) h2[k] = H2[k]/(Mt[k]*p.causal[k])
    p.causal.common = 1 # settings[ss,'p.causal.common']
    rs12 = r[1] * r.sign$r12.sign
    
    C1.sfbm = C.sfbm$C1.sfbm; C2.sfbm = C.sfbm$C2.sfbm
    # initial values
    hsqh1 = hsqh2 = numeric()
    pm1 = rep(0,Mt[1]); pm2 = rep(0,Mt[2])
    hsqh1[1:niter] = h2[1]; hsqh2[1:niter] = h2[2]
    bjh1 = beta_init[[1]]; bjh2 = beta_init[[2]]
    
    p_12_11 = min(p.causal1, p.causal2)*p.causal.common; p_12_10 = p.causal1-p_12_11; p_12_01 = p.causal2-p_12_11
    p_12_00 = 1 - p_12_11 - p_12_10 - p_12_01
    
    for (s in 2:niter){
      # update by block
      res12 = update_beta2(snp.index$snp.index12, indmat[,c(1,2)], C1.sfbm, C2.sfbm, Bh[[1]], Bh[[2]],
                           bjh1, bjh2, sigmasq1, sigmasq2, hsqh1[s], hsqh2[s], rs12, p_12_11, p_12_10, p_12_01, sparse[c(1,2)])
      bjh1[indmat[snp.index$snp.index12,1]] = res12[[1]][1,]; bjh2[indmat[snp.index$snp.index12,2]] = res12[[1]][2,];
      
      res1 = update_beta1(indmat[snp.index$snp.index1,1], C1.sfbm, Bh[[1]], bjh1, sigmasq1, hsqh1[s], p.causal1, sparse[1])
      bjh1[indmat[snp.index$snp.index1,1]] = res1[1,]; 
      res2 = update_beta1(indmat[snp.index$snp.index2,2], C2.sfbm, Bh[[2]], bjh2, sigmasq2, hsqh2[s], p.causal2, sparse[2])
      bjh2[indmat[snp.index$snp.index2,2]] = res2[1,];
      
      if (s>n.burnin){
        pm1[indmat[snp.index$snp.index12,1]] = pm1[indmat[snp.index$snp.index12,1]] + res12[[2]][1,]; pm2[indmat[snp.index$snp.index12,2]] = pm2[indmat[snp.index$snp.index12,2]] + res12[[2]][2,]
        pm1[indmat[snp.index$snp.index1,1]] = pm1[indmat[snp.index$snp.index1,1]] + res1[2,]; pm2[indmat[snp.index$snp.index2,2]] = pm2[indmat[snp.index$snp.index2,2]] + res2[2,]
      }
      progBar(s, niter, per=20)
    }
    pm1 = pm1/(niter-n.burnin); pm2 = pm2/(niter-n.burnin)
    names(pm1) = snpinfo$rsid[tem[[1]]]; names(pm2) = snpinfo$rsid[tem[[2]]]

    return(list(pm1 = pm1, pm2 = pm2))
  }
  if (K == 3){
    sigmasq1 = sigmasq[[1]]; sigmasq2 = sigmasq[[2]]; sigmasq3 = sigmasq[[3]]
    p.causal = settings[ss,paste0('p.causal',1:K)]; p.causal = as.numeric(p.causal)
    p.causal1 = p.causal[1]; p.causal2 = p.causal[2]; p.causal3 = p.causal[3]
    h2 = numeric()
    for (k in 1:K) h2[k] = H2[k]/(Mt[k]*p.causal[k])
    p.causal.common = 1 # settings[ss,'p.causal.common']
    rs12 = r[1] * r.sign$r12.sign; rs13 = r[2] * r.sign$r13.sign; rs23 = r[3] * r.sign$r23.sign
    
    C1.sfbm = C.sfbm$C1.sfbm; C2.sfbm = C.sfbm$C2.sfbm; C3.sfbm = C.sfbm$C3.sfbm
    # initial values
    hsqh1 = hsqh2 = hsqh3 = numeric()
    pm1 = rep(0,Mt[1]); pm2 = rep(0,Mt[2]); pm3 = rep(0,Mt[3])
    #set.seed(ss+chain+2020)
    hsqh1[1:niter] = h2[1]; hsqh2[1:niter] = h2[2]; hsqh3[1:niter] = h2[3]
    bjh1 = beta_init[[1]]; bjh2 = beta_init[[2]]; bjh3 = beta_init[[3]]
    
    p_12_11 = min(p.causal1, p.causal2)*p.causal.common; p_12_10 = p.causal1-p_12_11; p_12_01 = p.causal2-p_12_11
    p_12_00 = 1 - p_12_11 - p_12_10 - p_12_01
    p_13_11 = min(p.causal1, p.causal3)*p.causal.common; p_13_10 = p.causal1-p_13_11; p_13_01 = p.causal3-p_13_11
    p_13_00 = 1 - p_13_11 - p_13_10 - p_13_01
    p_23_11 = min(p.causal2, p.causal3)*p.causal.common; p_23_10 = p.causal2-p_23_11; p_23_01 = p.causal3-p_23_11
    p_23_00 = 1 - p_23_11 - p_23_10 - p_23_01
    
    p123_111 = min(p.causal1, p.causal2, p.causal3)*(p.causal.common^2)
    p123_110 = p_12_11 - p123_111; p123_101 = p_13_11 - p123_111; p123_011 = p_23_11 - p123_111
    p123_100 = p_13_10 - p123_110; p123_010 = p_12_01 - p123_011; p123_001 = p_13_01 - p123_011
    p123_000 = p_12_00 - p123_001
    
    for (s in 2:niter){
      # update by block
      res123 = update_beta3(snp.index$snp.index123, indmat[,1:3], C1.sfbm, C2.sfbm, C3.sfbm,
                            Bh[[1]], Bh[[2]], Bh[[3]], bjh1, bjh2, bjh3,
                            sigmasq1, sigmasq2, sigmasq3, hsqh1[s], hsqh2[s], hsqh3[s], rs12, rs13, rs23,
                            p123_111, p123_110, p123_101, p123_011, p123_100, p123_010, p123_001, p123_000, sparse[1:3])
      bjh1[indmat[snp.index$snp.index123,1]] = res123[[1]][1,]; bjh2[indmat[snp.index$snp.index123,2]] = res123[[1]][2,]; bjh3[indmat[snp.index$snp.index123,3]] = res123[[1]][3,]
      
      res12 = update_beta2(snp.index$snp.index12, indmat[,c(1,2)], C1.sfbm, C2.sfbm, Bh[[1]], Bh[[2]],
                           bjh1, bjh2, sigmasq1, sigmasq2, hsqh1[s], hsqh2[s], rs12, p_12_11, p_12_10, p_12_01, sparse[c(1,2)])
      bjh1[indmat[snp.index$snp.index12,1]] = res12[[1]][1,]; bjh2[indmat[snp.index$snp.index12,2]] = res12[[1]][2,];
      res13 = update_beta2(snp.index$snp.index13, indmat[,c(1,3)], C1.sfbm, C3.sfbm, Bh[[1]], Bh[[3]],
                           bjh1, bjh3, sigmasq1, sigmasq3, hsqh1[s], hsqh3[s], rs13, p_13_11, p_13_10, p_13_01, sparse[c(1,3)])
      bjh1[indmat[snp.index$snp.index13,1]] = res13[[1]][1,]; bjh3[indmat[snp.index$snp.index13,3]] = res13[[1]][2,];
      res23 = update_beta2(snp.index$snp.index23, indmat[,c(2,3)], C2.sfbm, C3.sfbm, Bh[[2]], Bh[[3]],
                           bjh2, bjh3, sigmasq2, sigmasq3, hsqh2[s], hsqh3[s], rs23, p_23_11, p_23_10, p_23_01, sparse[c(2,3)])
      bjh2[indmat[snp.index$snp.index23,2]] = res23[[1]][1,]; bjh3[indmat[snp.index$snp.index23,3]] = res23[[1]][2,];
      
      res1 = update_beta1(indmat[snp.index$snp.index1,1], C1.sfbm, Bh[[1]], bjh1, sigmasq1, hsqh1[s], p.causal1, sparse[1])
      bjh1[indmat[snp.index$snp.index1,1]] = res1[1,]; 
      res2 = update_beta1(indmat[snp.index$snp.index2,2], C2.sfbm, Bh[[2]], bjh2, sigmasq2, hsqh2[s], p.causal2, sparse[2])
      bjh2[indmat[snp.index$snp.index2,2]] = res2[1,];
      res3 = update_beta1(indmat[snp.index$snp.index3,3], C3.sfbm, Bh[[3]], bjh3, sigmasq3, hsqh3[s], p.causal3, sparse[3])
      bjh3[indmat[snp.index$snp.index3,3]] = res3[1,]; 
      
      if (s>n.burnin){
        pm1[indmat[snp.index$snp.index123,1]] = pm1[indmat[snp.index$snp.index123,1]] + res123[[2]][1,];
        pm2[indmat[snp.index$snp.index123,2]] = pm2[indmat[snp.index$snp.index123,2]] + res123[[2]][2,]
        pm3[indmat[snp.index$snp.index123,3]] = pm3[indmat[snp.index$snp.index123,3]] + res123[[2]][3,]
        
        pm1[indmat[snp.index$snp.index12,1]] = pm1[indmat[snp.index$snp.index12,1]] + res12[[2]][1,]; pm2[indmat[snp.index$snp.index12,2]] = pm2[indmat[snp.index$snp.index12,2]] + res12[[2]][2,]
        pm1[indmat[snp.index$snp.index13,1]] = pm1[indmat[snp.index$snp.index13,1]] + res13[[2]][1,]; pm3[indmat[snp.index$snp.index13,3]] = pm3[indmat[snp.index$snp.index13,3]] + res13[[2]][2,]
        pm2[indmat[snp.index$snp.index23,2]] = pm2[indmat[snp.index$snp.index23,2]] + res23[[2]][1,]; pm3[indmat[snp.index$snp.index23,3]] = pm3[indmat[snp.index$snp.index23,3]] + res23[[2]][2,]
        
        pm1[indmat[snp.index$snp.index1,1]] = pm1[indmat[snp.index$snp.index1,1]] + res1[2,]; pm2[indmat[snp.index$snp.index2,2]] = pm2[indmat[snp.index$snp.index2,2]] + res2[2,]
        pm3[indmat[snp.index$snp.index3,3]] = pm3[indmat[snp.index$snp.index3,3]] + res3[2,];
      }
      progBar(s, niter, per=20)
    }
    pm1 = pm1/(niter-n.burnin); pm2 = pm2/(niter-n.burnin); pm3 = pm3/(niter-n.burnin)
    names(pm1) = snpinfo$rsid[tem[[1]]]; names(pm2) = snpinfo$rsid[tem[[2]]]; names(pm3) = snpinfo$rsid[tem[[3]]]
    
    return(list(pm1 = pm1, pm2 = pm2, pm3 = pm3))
  }
  if (K == 4){
    sigmasq1 = sigmasq[[1]]; sigmasq2 = sigmasq[[2]]; sigmasq3 = sigmasq[[3]]; sigmasq4 = sigmasq[[4]]
    p.causal = settings[ss,paste0('p.causal',1:K)]; p.causal = as.numeric(p.causal)
    p.causal1 = p.causal[1]; p.causal2 = p.causal[2]; p.causal3 = p.causal[3]; p.causal4 = p.causal[4]
    h2 = numeric()
    for (k in 1:K) h2[k] = H2[k]/(Mt[k]*p.causal[k])
    p.causal.common = 1 # settings[ss,'p.causal.common']
    rs12 = r[1] * r.sign$r12.sign; rs13 = r[2] * r.sign$r13.sign; rs14 = r[3] * r.sign$r14.sign
    rs23 = r[4] * r.sign$r23.sign; rs24 = r[5] * r.sign$r24.sign; rs34 = r[6] * r.sign$r34.sign
    
    C1.sfbm = C.sfbm$C1.sfbm; C2.sfbm = C.sfbm$C2.sfbm; C3.sfbm = C.sfbm$C3.sfbm; C4.sfbm = C.sfbm$C4.sfbm;
    # initial values
    hsqh1 = hsqh2 = hsqh3 = hsqh4 = numeric()
    pm1 = rep(0,Mt[1]); pm2 = rep(0,Mt[2]); pm3 = rep(0,Mt[3]); pm4 = rep(0,Mt[4]);
    #set.seed(ss+chain+2020)
    hsqh1[1:niter] = h2[1]; hsqh2[1:niter] = h2[2]; hsqh3[1:niter] = h2[3]; hsqh4[1:niter] = h2[4];
    bjh1 = beta_init[[1]]; bjh2 = beta_init[[2]]; bjh3 = beta_init[[3]]; bjh4 = beta_init[[4]];
    
    p_12_11 = min(p.causal1, p.causal2)*p.causal.common; p_12_10 = p.causal1-p_12_11; p_12_01 = p.causal2-p_12_11
    p_12_00 = 1 - p_12_11 - p_12_10 - p_12_01
    p_13_11 = min(p.causal1, p.causal3)*p.causal.common; p_13_10 = p.causal1-p_13_11; p_13_01 = p.causal3-p_13_11
    p_13_00 = 1 - p_13_11 - p_13_10 - p_13_01
    p_14_11 = min(p.causal1, p.causal4)*p.causal.common; p_14_10 = p.causal1-p_14_11; p_14_01 = p.causal4-p_14_11
    p_14_00 = 1 - p_14_11 - p_14_10 - p_14_01
    p_23_11 = min(p.causal2, p.causal3)*p.causal.common; p_23_10 = p.causal2-p_23_11; p_23_01 = p.causal3-p_23_11
    p_23_00 = 1 - p_23_11 - p_23_10 - p_23_01
    p_24_11 = min(p.causal2, p.causal4)*p.causal.common; p_24_10 = p.causal2-p_24_11; p_24_01 = p.causal4-p_24_11
    p_24_00 = 1 - p_24_11 - p_24_10 - p_24_01
    p_34_11 = min(p.causal3, p.causal4)*p.causal.common; p_34_10 = p.causal3-p_34_11; p_34_01 = p.causal4-p_34_11
    p_34_00 = 1 - p_34_11 - p_34_10 - p_34_01
    
    p123_111 = min(p.causal1, p.causal2, p.causal3)*(p.causal.common^2)
    p123_110 = p_12_11 - p123_111; p123_101 = p_13_11 - p123_111; p123_011 = p_23_11 - p123_111
    p123_100 = p_13_10 - p123_110; p123_010 = p_12_01 - p123_011; p123_001 = p_13_01 - p123_011
    p123_000 = p_12_00 - p123_001
    p124_111 = min(p.causal1, p.causal2, p.causal4)*(p.causal.common^2)
    p124_110 = p_12_11 - p124_111; p124_101 = p_14_11 - p124_111; p124_011 = p_24_11 - p124_111
    p124_100 = p_14_10 - p124_110; p124_010 = p_12_01 - p124_011; p124_001 = p_14_01 - p124_011
    p124_000 = p_12_00 - p124_001
    p134_111 = min(p.causal1, p.causal3, p.causal4)*(p.causal.common^2)
    p134_110 = p_13_11 - p134_111; p134_101 = p_14_11 - p134_111; p134_011 = p_34_11 - p134_111
    p134_100 = p_14_10 - p134_110; p134_010 = p_13_01 - p134_011; p134_001 = p_14_01 - p134_011
    p134_000 = p_13_00 - p134_001
    p234_111 = min(p.causal2, p.causal3, p.causal4)*(p.causal.common^2)
    p234_110 = p_23_11 - p234_111; p234_101 = p_24_11 - p234_111; p234_011 = p_34_11 - p234_111
    p234_100 = p_24_10 - p234_110; p234_010 = p_23_01 - p234_011; p234_001 = p_24_01 - p234_011
    p234_000 = p_23_00 - p234_001
    
    p1234_1111 = min(p.causal1, p.causal2, p.causal3,p.causal4)*(p.causal.common^3)
    p1234_1110 = p123_111 - p1234_1111; p1234_1101 = p124_111 - p1234_1111; p1234_1011 = p134_111 - p1234_1111; p1234_0111 = p234_111 - p1234_1111
    p1234_1100 = p123_110 - p1234_1101; p1234_1010 = p134_110 - p1234_1110; p1234_0110 = p234_110 - p1234_1110; 
    p1234_1001 = p124_101 - p1234_1011; p1234_0101 = p124_011 - p1234_0111; p1234_0011 = p134_011 - p1234_0111; 
    p1234_1000 = p134_100 - p1234_1100; p1234_0100 = p234_100 - p1234_1100; p1234_0010 = p123_001 - p1234_0011; p1234_0001 = p124_001 - p1234_0011; 
    p1234_0000 = p234_000 - p1234_1000 
    
    for (s in 2:niter){
      # update by block
      res1234 = update_beta4(snp.index$snp.index1234, indmat[,c(1,2,3,4)], C1.sfbm, C2.sfbm, C3.sfbm, C4.sfbm, Bh[[1]], Bh[[2]], Bh[[3]], Bh[[4]], bjh1, bjh2, bjh3, bjh4,
                             sigmasq1, sigmasq2, sigmasq3, sigmasq4, hsqh1[s], hsqh2[s], hsqh3[s], hsqh4[s],
                             rs12, rs13, rs14, rs23, rs24, rs34,
                             p1234_1111, p1234_1110, p1234_1101, p1234_1011, p1234_0111,
                             p1234_1100, p1234_1010, p1234_1001, p1234_0110, p1234_0101, p1234_0011,
                             p1234_1000, p1234_0100, p1234_0010, p1234_0001, p1234_0000, sparse)
      bjh1[indmat[snp.index$snp.index1234,1]] = res1234[[1]][1,]; bjh2[indmat[snp.index$snp.index1234,2]] = res1234[[1]][2,]
      bjh3[indmat[snp.index$snp.index1234,3]] = res1234[[1]][3,]; bjh4[indmat[snp.index$snp.index1234,4]] = res1234[[1]][4,]
      
      res123 = update_beta3(snp.index$snp.index123, indmat[,1:3], C1.sfbm, C2.sfbm, C3.sfbm,
                            Bh[[1]], Bh[[2]], Bh[[3]],  bjh1, bjh2, bjh3,
                            sigmasq1, sigmasq2, sigmasq3, hsqh1[s], hsqh2[s], hsqh3[s], rs12, rs13, rs23,
                            p123_111, p123_110, p123_101, p123_011, p123_100, p123_010, p123_001, p123_000, sparse[1:3])
      bjh1[indmat[snp.index$snp.index123,1]] = res123[[1]][1,]; bjh2[indmat[snp.index$snp.index123,2]] = res123[[1]][2,]; bjh3[indmat[snp.index$snp.index123,3]] = res123[[1]][3,]
      res124 = update_beta3(snp.index$snp.index124, indmat[,c(1,2,4)], C1.sfbm, C2.sfbm, C4.sfbm,
                            Bh[[1]], Bh[[2]], Bh[[4]],  bjh1, bjh2, bjh4,
                            sigmasq1, sigmasq2, sigmasq4, hsqh1[s], hsqh2[s], hsqh4[s], rs12, rs14, rs24,
                            p124_111, p124_110, p124_101, p124_011, p124_100, p124_010, p124_001, p124_000, sparse[c(1,2,4)])
      bjh1[indmat[snp.index$snp.index124,1]] = res124[[1]][1,]; bjh2[indmat[snp.index$snp.index124,2]] = res124[[1]][2,]; bjh4[indmat[snp.index$snp.index124,4]] = res124[[1]][3,]
      res134 = update_beta3(snp.index$snp.index134, indmat[,c(1,3,4)], C1.sfbm, C3.sfbm, C4.sfbm,
                            Bh[[1]], Bh[[3]], Bh[[4]],  bjh1, bjh3, bjh4,
                            sigmasq1, sigmasq3, sigmasq4, hsqh1[s], hsqh3[s], hsqh4[s], rs13, rs14, rs34,
                            p134_111, p134_110, p134_101, p134_011, p134_100, p134_010, p134_001, p134_000, sparse[c(1,3,4)])
      bjh1[indmat[snp.index$snp.index134,1]] = res134[[1]][1,]; bjh3[indmat[snp.index$snp.index134,3]] = res134[[1]][2,]; bjh4[indmat[snp.index$snp.index134,4]] = res134[[1]][3,]
      res234 = update_beta3(snp.index$snp.index234, indmat[,c(2,3,4)], C2.sfbm, C3.sfbm, C4.sfbm,
                            Bh[[2]], Bh[[3]], Bh[[4]],  bjh2, bjh3, bjh4,
                            sigmasq2, sigmasq3, sigmasq4, hsqh2[s], hsqh3[s], hsqh4[s], rs23, rs24, rs34,
                            p234_111, p234_110, p234_101, p234_011, p234_100, p234_010, p234_001, p234_000, sparse[c(2,3,4)])
      bjh2[indmat[snp.index$snp.index234,2]] = res234[[1]][1,]; bjh3[indmat[snp.index$snp.index234,3]] = res234[[1]][2,]; bjh4[indmat[snp.index$snp.index234,4]] = res234[[1]][3,]
      
      res12 = update_beta2(snp.index$snp.index12, indmat[,c(1,2)], C1.sfbm, C2.sfbm, Bh[[1]], Bh[[2]],
                           bjh1, bjh2, sigmasq1, sigmasq2, hsqh1[s], hsqh2[s], rs12, p_12_11, p_12_10, p_12_01, sparse[c(1,2)])
      bjh1[indmat[snp.index$snp.index12,1]] = res12[[1]][1,]; bjh2[indmat[snp.index$snp.index12,2]] = res12[[1]][2,];
      res13 = update_beta2(snp.index$snp.index13, indmat[,c(1,3)], C1.sfbm, C3.sfbm, Bh[[1]], Bh[[3]],
                           bjh1, bjh3, sigmasq1, sigmasq3, hsqh1[s], hsqh3[s], rs13, p_13_11, p_13_10, p_13_01, sparse[c(1,3)])
      bjh1[indmat[snp.index$snp.index13,1]] = res13[[1]][1,]; bjh3[indmat[snp.index$snp.index13,3]] = res13[[1]][2,];
      res14 = update_beta2(snp.index$snp.index14, indmat[,c(1,4)], C1.sfbm, C4.sfbm, Bh[[1]], Bh[[4]],
                           bjh1, bjh4, sigmasq1, sigmasq4, hsqh1[s], hsqh4[s], rs14, p_14_11, p_14_10, p_14_01, sparse[c(1,4)])
      bjh1[indmat[snp.index$snp.index14,1]] = res14[[1]][1,]; bjh4[indmat[snp.index$snp.index14,4]] = res14[[1]][2,];
      res23 = update_beta2(snp.index$snp.index23, indmat[,c(2,3)], C2.sfbm, C3.sfbm, Bh[[2]], Bh[[3]],
                           bjh2, bjh3, sigmasq2, sigmasq3, hsqh2[s], hsqh3[s], rs23, p_23_11, p_23_10, p_23_01, sparse[c(2,3)])
      bjh2[indmat[snp.index$snp.index23,2]] = res23[[1]][1,]; bjh3[indmat[snp.index$snp.index23,3]] = res23[[1]][2,];
      res24 = update_beta2(snp.index$snp.index24, indmat[,c(2,4)], C2.sfbm, C4.sfbm, Bh[[2]], Bh[[4]],
                           bjh2, bjh4, sigmasq2, sigmasq4, hsqh2[s], hsqh4[s], rs24, p_24_11, p_24_10, p_24_01, sparse[c(2,4)])
      bjh2[indmat[snp.index$snp.index24,2]] = res24[[1]][1,]; bjh4[indmat[snp.index$snp.index24,4]] = res24[[1]][2,];
      res34 = update_beta2(snp.index$snp.index34, indmat[,c(3,4)], C3.sfbm, C4.sfbm, Bh[[3]], Bh[[4]],
                           bjh3, bjh4, sigmasq3, sigmasq4, hsqh3[s], hsqh4[s], rs34, p_34_11, p_34_10, p_34_01, sparse[c(3,4)])
      bjh3[indmat[snp.index$snp.index34,3]] = res34[[1]][1,]; bjh4[indmat[snp.index$snp.index34,4]] = res34[[1]][2,];
      
      res1 = update_beta1(indmat[snp.index$snp.index1,1], C1.sfbm, Bh[[1]], bjh1, sigmasq1, hsqh1[s], p.causal1, sparse[1])
      bjh1[indmat[snp.index$snp.index1,1]] = res1[1,]; 
      res2 = update_beta1(indmat[snp.index$snp.index2,2], C2.sfbm, Bh[[2]], bjh2, sigmasq2, hsqh2[s], p.causal2, sparse[2])
      bjh2[indmat[snp.index$snp.index2,2]] = res2[1,];
      res3 = update_beta1(indmat[snp.index$snp.index3,3], C3.sfbm, Bh[[3]], bjh3, sigmasq3, hsqh3[s], p.causal3, sparse[3])
      bjh3[indmat[snp.index$snp.index3,3]] = res3[1,]; 
      res4 = update_beta1(indmat[snp.index$snp.index4,4], C4.sfbm, Bh[[4]], bjh4, sigmasq4, hsqh4[s], p.causal4, sparse[4])
      bjh4[indmat[snp.index$snp.index4,4]] = res4[1,];
      
      if (s>n.burnin){
        pm1[indmat[snp.index$snp.index1234,1]] = pm1[indmat[snp.index$snp.index1234,1]] + res1234[[2]][1,];
        pm2[indmat[snp.index$snp.index1234,2]] = pm2[indmat[snp.index$snp.index1234,2]] + res1234[[2]][2,];
        pm3[indmat[snp.index$snp.index1234,3]] = pm3[indmat[snp.index$snp.index1234,3]] + res1234[[2]][3,];
        pm4[indmat[snp.index$snp.index1234,4]] = pm4[indmat[snp.index$snp.index1234,4]] + res1234[[2]][4,];
        
        pm1[indmat[snp.index$snp.index123,1]] = pm1[indmat[snp.index$snp.index123,1]] + res123[[2]][1,];
        pm2[indmat[snp.index$snp.index123,2]] = pm2[indmat[snp.index$snp.index123,2]] + res123[[2]][2,]
        pm3[indmat[snp.index$snp.index123,3]] = pm3[indmat[snp.index$snp.index123,3]] + res123[[2]][3,]
        pm1[indmat[snp.index$snp.index124,1]] = pm1[indmat[snp.index$snp.index124,1]] + res124[[2]][1,];
        pm2[indmat[snp.index$snp.index124,2]] = pm2[indmat[snp.index$snp.index124,2]] + res124[[2]][2,]
        pm4[indmat[snp.index$snp.index124,4]] = pm4[indmat[snp.index$snp.index124,4]] + res124[[2]][3,]
        pm1[indmat[snp.index$snp.index134,1]] = pm1[indmat[snp.index$snp.index134,1]] + res134[[2]][1,];
        pm3[indmat[snp.index$snp.index134,3]] = pm3[indmat[snp.index$snp.index134,3]] + res134[[2]][2,]
        pm4[indmat[snp.index$snp.index134,4]] = pm4[indmat[snp.index$snp.index134,4]] + res134[[2]][3,]
        pm2[indmat[snp.index$snp.index234,2]] = pm2[indmat[snp.index$snp.index234,2]] + res234[[2]][1,];
        pm3[indmat[snp.index$snp.index234,3]] = pm3[indmat[snp.index$snp.index234,3]] + res234[[2]][2,]
        pm4[indmat[snp.index$snp.index234,4]] = pm4[indmat[snp.index$snp.index234,4]] + res234[[2]][3,]
        
        pm1[indmat[snp.index$snp.index12,1]] = pm1[indmat[snp.index$snp.index12,1]] + res12[[2]][1,]; pm2[indmat[snp.index$snp.index12,2]] = pm2[indmat[snp.index$snp.index12,2]] + res12[[2]][2,]
        pm1[indmat[snp.index$snp.index13,1]] = pm1[indmat[snp.index$snp.index13,1]] + res13[[2]][1,]; pm3[indmat[snp.index$snp.index13,3]] = pm3[indmat[snp.index$snp.index13,3]] + res13[[2]][2,]
        pm1[indmat[snp.index$snp.index14,1]] = pm1[indmat[snp.index$snp.index14,1]] + res14[[2]][1,]; pm4[indmat[snp.index$snp.index14,4]] = pm4[indmat[snp.index$snp.index14,4]] + res14[[2]][2,]
        pm2[indmat[snp.index$snp.index23,2]] = pm2[indmat[snp.index$snp.index23,2]] + res23[[2]][1,]; pm3[indmat[snp.index$snp.index23,3]] = pm3[indmat[snp.index$snp.index23,3]] + res23[[2]][2,]
        pm2[indmat[snp.index$snp.index24,2]] = pm2[indmat[snp.index$snp.index24,2]] + res24[[2]][1,]; pm4[indmat[snp.index$snp.index24,4]] = pm4[indmat[snp.index$snp.index24,4]] + res24[[2]][2,]
        pm3[indmat[snp.index$snp.index34,3]] = pm3[indmat[snp.index$snp.index34,3]] + res34[[2]][1,]; pm4[indmat[snp.index$snp.index34,4]] = pm4[indmat[snp.index$snp.index34,4]] + res34[[2]][2,]
        
        pm1[indmat[snp.index$snp.index1,1]] = pm1[indmat[snp.index$snp.index1,1]] + res1[2,]; pm2[indmat[snp.index$snp.index2,2]] = pm2[indmat[snp.index$snp.index2,2]] + res2[2,]
        pm3[indmat[snp.index$snp.index3,3]] = pm3[indmat[snp.index$snp.index3,3]] + res3[2,]; pm4[indmat[snp.index$snp.index4,4]] = pm4[indmat[snp.index$snp.index4,4]] + res4[2,]
      }
      progBar(s, niter, per=20)
    }
    pm1 = pm1/(niter-n.burnin); pm2 = pm2/(niter-n.burnin); pm3 = pm3/(niter-n.burnin); pm4 = pm4/(niter-n.burnin)
    names(pm1) = snpinfo$rsid[tem[[1]]]; names(pm2) = snpinfo$rsid[tem[[2]]]
    names(pm3) = snpinfo$rsid[tem[[3]]]; names(pm4) = snpinfo$rsid[tem[[4]]]
    
    return(list(pm1 = pm1, pm2 = pm2, pm3 = pm3, pm4 = pm4))
  }
  if (K == 5){
    sigmasq1 = sigmasq[[1]]; sigmasq2 = sigmasq[[2]]; sigmasq3 = sigmasq[[3]]; sigmasq4 = sigmasq[[4]]; sigmasq5 = sigmasq[[5]]
    p.causal = settings[ss,paste0('p.causal',1:K)]; p.causal = as.numeric(p.causal)
    p.causal1 = p.causal[1]; p.causal2 = p.causal[2]; p.causal3 = p.causal[3]; p.causal4 = p.causal[4]; p.causal5 = p.causal[5]
    h2 = numeric()
    for (k in 1:K){
      h2[k] = H2[k]/(Mt[k]*p.causal[k])
    }
    p.causal.common = 1
    rs12 = r[1] * r.sign$r12.sign; 
    rs13 = r[2] * r.sign$r13.sign; rs14 = r[3] * r.sign$r14.sign; rs15 = r[4] * r.sign$r15.sign;
    rs23 = r[5] * r.sign$r23.sign; rs24 = r[6] * r.sign$r24.sign; rs25 = r[7] * r.sign$r25.sign; 
    rs34 = r[8] * r.sign$r34.sign; rs35 = r[9] * r.sign$r35.sign; rs45 = r[10] * r.sign$r45.sign
    
    C1.sfbm = C.sfbm$C1.sfbm; C2.sfbm = C.sfbm$C2.sfbm; C3.sfbm = C.sfbm$C3.sfbm; C4.sfbm = C.sfbm$C4.sfbm; C5.sfbm = C.sfbm$C5.sfbm
    # initial values
    hsqh1 = hsqh2 = hsqh3 = hsqh4 = hsqh5 = numeric()
    pm1 = rep(0,Mt[1]); pm2 = rep(0,Mt[2]); pm3 = rep(0,Mt[3]); pm4 = rep(0,Mt[4]); pm5 = rep(0,Mt[5]);
    
    hsqh1[1:niter] = h2[1]; hsqh2[1:niter] = h2[2]; hsqh3[1:niter] = h2[3]; hsqh4[1:niter] = h2[4]; hsqh5[1:niter] = h2[5]
    bjh1 = beta_init[[1]]; bjh2 = beta_init[[2]]; bjh3 = beta_init[[3]]; bjh4 = beta_init[[4]]; bjh5 = beta_init[[5]]; 
    
    p_12_11 = min(p.causal1, p.causal2)*p.causal.common; p_12_10 = p.causal1-p_12_11; p_12_01 = p.causal2-p_12_11
    p_12_00 = 1 - p_12_11 - p_12_10 - p_12_01
    p_13_11 = min(p.causal1, p.causal3)*p.causal.common; p_13_10 = p.causal1-p_13_11; p_13_01 = p.causal3-p_13_11
    p_13_00 = 1 - p_13_11 - p_13_10 - p_13_01
    p_14_11 = min(p.causal1, p.causal4)*p.causal.common; p_14_10 = p.causal1-p_14_11; p_14_01 = p.causal4-p_14_11
    p_14_00 = 1 - p_14_11 - p_14_10 - p_14_01
    p_15_11 = min(p.causal1, p.causal5)*p.causal.common; p_15_10 = p.causal1-p_15_11; p_15_01 = p.causal5-p_15_11
    p_15_00 = 1 - p_15_11 - p_15_10 - p_15_01
    p_23_11 = min(p.causal2, p.causal3)*p.causal.common; p_23_10 = p.causal2-p_23_11; p_23_01 = p.causal3-p_23_11
    p_23_00 = 1 - p_23_11 - p_23_10 - p_23_01
    p_24_11 = min(p.causal2, p.causal4)*p.causal.common; p_24_10 = p.causal2-p_24_11; p_24_01 = p.causal4-p_24_11
    p_24_00 = 1 - p_24_11 - p_24_10 - p_24_01
    p_25_11 = min(p.causal2, p.causal5)*p.causal.common; p_25_10 = p.causal2-p_25_11; p_25_01 = p.causal5-p_25_11
    p_25_00 = 1 - p_25_11 - p_25_10 - p_25_01
    p_34_11 = min(p.causal3, p.causal4)*p.causal.common; p_34_10 = p.causal3-p_34_11; p_34_01 = p.causal4-p_34_11
    p_34_00 = 1 - p_34_11 - p_34_10 - p_34_01
    p_35_11 = min(p.causal3, p.causal5)*p.causal.common; p_35_10 = p.causal3-p_35_11; p_35_01 = p.causal5-p_35_11
    p_35_00 = 1 - p_35_11 - p_35_10 - p_35_01
    p_45_11 = min(p.causal4, p.causal5)*p.causal.common; p_45_10 = p.causal4-p_45_11; p_45_01 = p.causal5-p_45_11
    p_45_00 = 1 - p_45_11 - p_45_10 - p_45_01
    
    p123_111 = min(p.causal1, p.causal2, p.causal3)*(p.causal.common^2)
    p123_110 = p_12_11 - p123_111; p123_101 = p_13_11 - p123_111; p123_011 = p_23_11 - p123_111
    p123_100 = p_13_10 - p123_110; p123_010 = p_12_01 - p123_011; p123_001 = p_13_01 - p123_011
    p123_000 = p_12_00 - p123_001
    p124_111 = min(p.causal1, p.causal2, p.causal4)*(p.causal.common^2)
    p124_110 = p_12_11 - p124_111; p124_101 = p_14_11 - p124_111; p124_011 = p_24_11 - p124_111
    p124_100 = p_14_10 - p124_110; p124_010 = p_12_01 - p124_011; p124_001 = p_14_01 - p124_011
    p124_000 = p_12_00 - p124_001
    p125_111 = min(p.causal1, p.causal2, p.causal5)*(p.causal.common^2)
    p125_110 = p_12_11 - p125_111; p125_101 = p_15_11 - p125_111; p125_011 = p_25_11 - p125_111
    p125_100 = p_15_10 - p125_110; p125_010 = p_12_01 - p125_011; p125_001 = p_15_01 - p125_011
    p125_000 = p_12_00 - p125_001
    p134_111 = min(p.causal1, p.causal3, p.causal4)*(p.causal.common^2)
    p134_110 = p_13_11 - p134_111; p134_101 = p_14_11 - p134_111; p134_011 = p_34_11 - p134_111
    p134_100 = p_14_10 - p134_110; p134_010 = p_13_01 - p134_011; p134_001 = p_14_01 - p134_011
    p134_000 = p_13_00 - p134_001
    p135_111 = min(p.causal1, p.causal3, p.causal5)*(p.causal.common^2)
    p135_110 = p_13_11 - p135_111; p135_101 = p_15_11 - p135_111; p135_011 = p_35_11 - p135_111
    p135_100 = p_15_10 - p135_110; p135_010 = p_13_01 - p135_011; p135_001 = p_15_01 - p135_011
    p135_000 = p_13_00 - p135_001
    p145_111 = min(p.causal1, p.causal4, p.causal5)*(p.causal.common^2)
    p145_110 = p_14_11 - p145_111; p145_101 = p_15_11 - p145_111; p145_011 = p_45_11 - p145_111
    p145_100 = p_15_10 - p145_110; p145_010 = p_14_01 - p145_011; p145_001 = p_15_01 - p145_011
    p145_000 = p_14_00 - p145_001
    p234_111 = min(p.causal2, p.causal3, p.causal4)*(p.causal.common^2)
    p234_110 = p_23_11 - p234_111; p234_101 = p_24_11 - p234_111; p234_011 = p_34_11 - p234_111
    p234_100 = p_24_10 - p234_110; p234_010 = p_23_01 - p234_011; p234_001 = p_24_01 - p234_011
    p234_000 = p_23_00 - p234_001
    p235_111 = min(p.causal2, p.causal3, p.causal5)*(p.causal.common^2)
    p235_110 = p_23_11 - p235_111; p235_101 = p_25_11 - p235_111; p235_011 = p_35_11 - p235_111
    p235_100 = p_25_10 - p235_110; p235_010 = p_25_01 - p235_011; p235_001 = p_25_01 - p235_011
    p235_000 = p_23_00 - p235_001
    p245_111 = min(p.causal2, p.causal4, p.causal5)*(p.causal.common^2)
    p245_110 = p_24_11 - p245_111; p245_101 = p_25_11 - p245_111; p245_011 = p_45_11 - p245_111
    p245_100 = p_25_10 - p245_110; p245_010 = p_25_01 - p245_011; p245_001 = p_25_01 - p245_011
    p245_000 = p_24_00 - p245_001
    p345_111 = min(p.causal3, p.causal4, p.causal5)*(p.causal.common^2)
    p345_110 = p_34_11 - p345_111; p345_101 = p_35_11 - p345_111; p345_011 = p_45_11 - p345_111
    p345_100 = p_35_10 - p345_110; p345_010 = p_35_01 - p345_011; p345_001 = p_35_01 - p345_011
    p345_000 = p_34_00 - p345_001
    
    p1234_1111 = min(p.causal1, p.causal2, p.causal3,p.causal4)*(p.causal.common^3)
    p1234_1110 = p123_111 - p1234_1111; p1234_1101 = p124_111 - p1234_1111; p1234_1011 = p134_111 - p1234_1111; p1234_0111 = p234_111 - p1234_1111
    p1234_1100 = p123_110 - p1234_1101; p1234_1010 = p134_110 - p1234_1110; p1234_0110 = p234_110 - p1234_1110; 
    p1234_1001 = p124_101 - p1234_1011; p1234_0101 = p124_011 - p1234_0111; p1234_0011 = p134_011 - p1234_0111; 
    p1234_1000 = p134_100 - p1234_1100; p1234_0100 = p234_100 - p1234_1100; p1234_0010 = p123_001 - p1234_0011; p1234_0001 = p124_001 - p1234_0011; 
    p1234_0000 = p234_000 - p1234_1000 
    
    p1235_1111 = min(p.causal1, p.causal2, p.causal3,p.causal5)*(p.causal.common^3)
    p1235_1110 = p123_111 - p1235_1111; p1235_1101 = p125_111 - p1235_1111; p1235_1011 = p135_111 - p1235_1111; p1235_0111 = p235_111 - p1235_1111
    p1235_1100 = p123_110 - p1235_1101; p1235_1010 = p135_110 - p1235_1110; p1235_0110 = p235_110 - p1235_1110; 
    p1235_1001 = p125_101 - p1235_1011; p1235_0101 = p125_011 - p1235_0111; p1235_0011 = p135_011 - p1235_0111; 
    p1235_1000 = p135_100 - p1235_1100; p1235_0100 = p235_100 - p1235_1100; p1235_0010 = p123_001 - p1235_0011; p1235_0001 = p125_001 - p1235_0011; 
    p1235_0000 = p235_000 - p1235_1000 
    
    p1245_1111 = min(p.causal1, p.causal2, p.causal4,p.causal5)*(p.causal.common^3)
    p1245_1110 = p124_111 - p1245_1111; p1245_1101 = p125_111 - p1245_1111; p1245_1011 = p145_111 - p1245_1111; p1245_0111 = p245_111 - p1245_1111
    p1245_1100 = p124_110 - p1245_1101; p1245_1010 = p145_110 - p1245_1110; p1245_0110 = p245_110 - p1245_1110; 
    p1245_1001 = p125_101 - p1245_1011; p1245_0101 = p125_011 - p1245_0111; p1245_0011 = p145_011 - p1245_0111; 
    p1245_1000 = p145_100 - p1245_1100; p1245_0100 = p245_100 - p1245_1100; p1245_0010 = p124_001 - p1245_0011; p1245_0001 = p125_001 - p1245_0011; 
    p1245_0000 = p245_000 - p1245_1000 
    
    p1345_1111 = min(p.causal1, p.causal3, p.causal4,p.causal5)*(p.causal.common^3)
    p1345_1110 = p134_111 - p1345_1111; p1345_1101 = p135_111 - p1345_1111; p1345_1011 = p145_111 - p1345_1111; p1345_0111 = p345_111 - p1345_1111
    p1345_1100 = p134_110 - p1345_1101; p1345_1010 = p145_110 - p1345_1110; p1345_0110 = p345_110 - p1345_1110; 
    p1345_1001 = p135_101 - p1345_1011; p1345_0101 = p135_011 - p1345_0111; p1345_0011 = p145_011 - p1345_0111; 
    p1345_1000 = p145_100 - p1345_1100; p1345_0100 = p345_100 - p1345_1100; p1345_0010 = p134_001 - p1345_0011; p1345_0001 = p135_001 - p1345_0011; 
    p1345_0000 = p345_000 - p1345_1000 
    
    p2345_1111 = min(p.causal2, p.causal3, p.causal4,p.causal5)*(p.causal.common^3)
    p2345_1110 = p234_111 - p2345_1111; p2345_1101 = p235_111 - p2345_1111; p2345_1011 = p245_111 - p2345_1111; p2345_0111 = p345_111 - p2345_1111
    p2345_1100 = p234_110 - p2345_1101; p2345_1010 = p245_110 - p2345_1110; p2345_0110 = p345_110 - p2345_1110; 
    p2345_1001 = p235_101 - p2345_1011; p2345_0101 = p235_011 - p2345_0111; p2345_0011 = p245_011 - p2345_0111; 
    p2345_1000 = p245_100 - p2345_1100; p2345_0100 = p345_100 - p2345_1100; p2345_0010 = p234_001 - p2345_0011; p2345_0001 = p235_001 - p2345_0011; 
    p2345_0000 = p345_000 - p2345_1000 
    
    p12345_11111 = min(p.causal1, p.causal2, p.causal3,p.causal4,p.causal5)*(p.causal.common^4)
    p12345_11110 = p1234_1111 - p12345_11111; p12345_11101 = p1235_1111 - p12345_11111; p12345_11011 = p1245_1111 - p12345_11111; p12345_10111 = p1345_1111 - p12345_11111; p12345_01111 = p2345_1111 - p12345_11111;
    p12345_11100 = p1234_1110 - p12345_11101; p12345_11010 = p1234_1101 - p12345_11011; p12345_10110 = p1234_1011 - p12345_10111; p12345_01110 = p1234_0111 - p12345_01111; 
    p12345_11001 = p1235_1101 - p12345_11011; p12345_10101 = p1235_1011 - p12345_10111; p12345_01101 = p1235_0111 - p12345_01111; 
    p12345_10011 = p1245_1011 - p12345_10111; p12345_01011 = p1245_0111 - p12345_01111; 
    p12345_00111 = p1345_0111 - p12345_01111; 
    p12345_00011 = p1245_0011 - p12345_00111; p12345_00101 = p1235_0011 - p12345_00111; p12345_00110 = p1234_0011 - p12345_00111; 
    p12345_01001 = p1245_0101 - p12345_01101; p12345_01010 = p1245_0110 - p12345_01110; p12345_01100 = p2345_1100 - p12345_11100; 
    p12345_10001 = p1345_1001 - p12345_11001; p12345_10010 = p1345_1010 - p12345_11010; p12345_10100 = p1345_1100 - p12345_11100; p12345_11000 = p1245_1100 - p12345_11100;
    
    p12345_00001 = p1345_0001 - p12345_01001; p12345_00010 = p1345_0010 - p12345_01010; p12345_00100 = p1234_0010 - p12345_00101;
    p12345_01000 = p1234_0100 - p12345_01001; p12345_10000 = p1234_1000 - p12345_10001; 
    p12345_00000 = p1234_0000 - p12345_00001; 
    
    pall = c(p12345_11111, p12345_11110, p12345_11101, p12345_11011, p12345_10111,
             p12345_11100, p12345_11010, p12345_11001, p12345_10110, p12345_10101, p12345_10011,
             p12345_11000, p12345_10100, p12345_10010, p12345_10001, p12345_10000, 
             p12345_01111, p12345_01110, p12345_01101, p12345_01011, p12345_00111, 
             p12345_01100, p12345_01010, p12345_01001, p12345_00110, p12345_00101, p12345_00011,
             p12345_01000, p12345_00100, p12345_00010, p12345_00001, p12345_00000)
    for (s in 2:niter){
      res12345 = update_beta5(snp.index$snp.index12345, indmat, C1.sfbm, C2.sfbm, C3.sfbm, C4.sfbm, C5.sfbm, 
                              Bh[[1]], Bh[[2]], Bh[[3]], Bh[[4]], Bh[[5]], bjh1, bjh2, bjh3, bjh4, bjh5,
                              sigmasq1, sigmasq2, sigmasq3, sigmasq4, sigmasq5, hsqh1[s], hsqh2[s], hsqh3[s], hsqh4[s], hsqh5[s],
                              rs12, rs13, rs14, rs15, rs23, rs24, rs25, rs34, rs35, rs45, 
                              pall, sparse)
      bjh1[indmat[snp.index$snp.index12345,1]] = res12345[[1]][1,]; bjh2[indmat[snp.index$snp.index12345,2]] = res12345[[1]][2,]
      bjh3[indmat[snp.index$snp.index12345,3]] = res12345[[1]][3,]; bjh4[indmat[snp.index$snp.index12345,4]] = res12345[[1]][4,]
      bjh5[indmat[snp.index$snp.index12345,5]] = res12345[[1]][5,]
      # 
      res1234 = update_beta4(snp.index$snp.index1234, indmat[,c(1,2,3,4)], C1.sfbm, C2.sfbm, C3.sfbm, C4.sfbm, Bh[[1]], Bh[[2]], Bh[[3]], Bh[[4]], bjh1, bjh2, bjh3, bjh4,
                             sigmasq1, sigmasq2, sigmasq3, sigmasq4, hsqh1[s], hsqh2[s], hsqh3[s], hsqh4[s],
                             rs12, rs13, rs14, rs23, rs24, rs34,
                             p1234_1111, p1234_1110, p1234_1101, p1234_1011, p1234_0111,
                             p1234_1100, p1234_1010, p1234_1001, p1234_0110, p1234_0101, p1234_0011,
                             p1234_1000, p1234_0100, p1234_0010, p1234_0001, p1234_0000, sparse)
      bjh1[indmat[snp.index$snp.index1234,1]] = res1234[[1]][1,]; bjh2[indmat[snp.index$snp.index1234,2]] = res1234[[1]][2,]
      bjh3[indmat[snp.index$snp.index1234,3]] = res1234[[1]][3,]; bjh4[indmat[snp.index$snp.index1234,4]] = res1234[[1]][4,]
      
      res1235 = update_beta4(snp.index$snp.index1235, indmat[,c(1,2,3,5)], C1.sfbm, C2.sfbm, C3.sfbm, C5.sfbm, Bh[[1]], Bh[[2]], Bh[[3]], Bh[[5]], bjh1, bjh2, bjh3, bjh5,
                             sigmasq1, sigmasq2, sigmasq3, sigmasq5, hsqh1[s], hsqh2[s], hsqh3[s], hsqh5[s],
                             rs12, rs13, rs15, rs23, rs25, rs35,
                             p1235_1111, p1235_1110, p1235_1101, p1235_1011, p1235_0111,
                             p1235_1100, p1235_1010, p1235_1001, p1235_0110, p1235_0101, p1235_0011,
                             p1235_1000, p1235_0100, p1235_0010, p1235_0001, p1235_0000, sparse)
      bjh1[indmat[snp.index$snp.index1235,1]] = res1235[[1]][1,]; bjh2[indmat[snp.index$snp.index1235,2]] = res1235[[1]][2,]
      bjh3[indmat[snp.index$snp.index1235,3]] = res1235[[1]][3,]; bjh5[indmat[snp.index$snp.index1235,5]] = res1235[[1]][4,]
      
      res1245 = update_beta4(snp.index$snp.index1245, indmat[,c(1,2,4,5)], C1.sfbm, C2.sfbm, C4.sfbm, C5.sfbm, Bh[[1]], Bh[[2]], Bh[[4]], Bh[[5]], bjh1, bjh2, bjh4, bjh5,
                             sigmasq1, sigmasq2, sigmasq4, sigmasq5, hsqh1[s], hsqh2[s], hsqh4[s], hsqh5[s],
                             rs12, rs14, rs15, rs24, rs25, rs45,
                             p1245_1111, p1245_1110, p1245_1101, p1245_1011, p1245_0111,
                             p1245_1100, p1245_1010, p1245_1001, p1245_0110, p1245_0101, p1245_0011,
                             p1245_1000, p1245_0100, p1245_0010, p1245_0001, p1245_0000, sparse)
      bjh1[indmat[snp.index$snp.index1245,1]] = res1245[[1]][1,]; bjh2[indmat[snp.index$snp.index1245,2]] = res1245[[1]][2,]
      bjh4[indmat[snp.index$snp.index1245,4]] = res1245[[1]][3,]; bjh5[indmat[snp.index$snp.index1245,5]] = res1245[[1]][4,]
      
      res1345 = update_beta4(snp.index$snp.index1345, indmat[,c(1,3,4,5)], C1.sfbm, C3.sfbm, C4.sfbm, C5.sfbm, Bh[[1]], Bh[[3]], Bh[[4]], Bh[[5]], bjh1, bjh3, bjh4, bjh5,
                             sigmasq1, sigmasq3, sigmasq4, sigmasq5, hsqh1[s], hsqh3[s], hsqh4[s], hsqh5[s],
                             rs13, rs14, rs15, rs34, rs35, rs45,
                             p1345_1111, p1345_1110, p1345_1101, p1345_1011, p1345_0111,
                             p1345_1100, p1345_1010, p1345_1001, p1345_0110, p1345_0101, p1345_0011,
                             p1345_1000, p1345_0100, p1345_0010, p1345_0001, p1345_0000, sparse[c(1,3,4,5)])
      bjh1[indmat[snp.index$snp.index1345,1]] = res1345[[1]][1,]; bjh3[indmat[snp.index$snp.index1345,3]] = res1345[[1]][2,]
      bjh4[indmat[snp.index$snp.index1345,4]] = res1345[[1]][3,]; bjh5[indmat[snp.index$snp.index1345,5]] = res1345[[1]][4,]
      
      res2345 = update_beta4(snp.index$snp.index2345, indmat[,c(2,3,4,5)], C2.sfbm, C3.sfbm, C4.sfbm, C5.sfbm, Bh[[2]], Bh[[3]], Bh[[4]], Bh[[5]], bjh2, bjh3, bjh4, bjh5,
                             sigmasq2, sigmasq3, sigmasq4, sigmasq5, hsqh2[s], hsqh3[s], hsqh4[s], hsqh5[s],
                             rs23, rs24, rs25, rs34, rs35, rs45,
                             p2345_1111, p2345_1110, p2345_1101, p2345_1011, p2345_0111,
                             p2345_1100, p2345_1010, p2345_1001, p2345_0110, p2345_0101, p2345_0011,
                             p2345_1000, p2345_0100, p2345_0010, p2345_0001, p2345_0000, sparse)
      bjh2[indmat[snp.index$snp.index2345,2]] = res2345[[1]][1,]; bjh3[indmat[snp.index$snp.index2345,3]] = res2345[[1]][2,]
      bjh4[indmat[snp.index$snp.index2345,4]] = res2345[[1]][3,]; bjh5[indmat[snp.index$snp.index2345,5]] = res2345[[1]][4,]
      
      res123 = update_beta3(snp.index$snp.index123, indmat[,1:3], C1.sfbm, C2.sfbm, C3.sfbm,
                            Bh[[1]], Bh[[2]], Bh[[3]],  bjh1, bjh2, bjh3,
                            sigmasq1, sigmasq2, sigmasq3, hsqh1[s], hsqh2[s], hsqh3[s], rs12, rs13, rs23,
                            p123_111, p123_110, p123_101, p123_011, p123_100, p123_010, p123_001, p123_000, sparse[1:3])
      bjh1[indmat[snp.index$snp.index123,1]] = res123[[1]][1,]; bjh2[indmat[snp.index$snp.index123,2]] = res123[[1]][2,]; bjh3[indmat[snp.index$snp.index123,3]] = res123[[1]][3,]
      res124 = update_beta3(snp.index$snp.index124, indmat[,c(1,2,4)], C1.sfbm, C2.sfbm, C4.sfbm,
                            Bh[[1]], Bh[[2]], Bh[[4]],  bjh1, bjh2, bjh4,
                            sigmasq1, sigmasq2, sigmasq4, hsqh1[s], hsqh2[s], hsqh4[s], rs12, rs14, rs24,
                            p124_111, p124_110, p124_101, p124_011, p124_100, p124_010, p124_001, p124_000, sparse[c(1,2,4)])
      bjh1[indmat[snp.index$snp.index124,1]] = res124[[1]][1,]; bjh2[indmat[snp.index$snp.index124,2]] = res124[[1]][2,]; bjh4[indmat[snp.index$snp.index124,4]] = res124[[1]][3,]
      res125 = update_beta3(snp.index$snp.index125, indmat[,c(1,2,5)], C1.sfbm, C2.sfbm, C5.sfbm,
                            Bh[[1]], Bh[[2]], Bh[[5]],  bjh1, bjh2, bjh5,
                            sigmasq1, sigmasq2, sigmasq5, hsqh1[s], hsqh2[s], hsqh5[s], rs12, rs15, rs25,
                            p125_111, p125_110, p125_101, p125_011, p125_100, p125_010, p125_001, p125_000, sparse[c(1,2,5)])
      bjh1[indmat[snp.index$snp.index125,1]] = res125[[1]][1,]; bjh2[indmat[snp.index$snp.index125,2]] = res125[[1]][2,]; bjh5[indmat[snp.index$snp.index125,5]] = res125[[1]][3,]
      res134 = update_beta3(snp.index$snp.index134, indmat[,c(1,3,4)], C1.sfbm, C3.sfbm, C4.sfbm,
                            Bh[[1]], Bh[[3]], Bh[[4]],  bjh1, bjh3, bjh4,
                            sigmasq1, sigmasq3, sigmasq4, hsqh1[s], hsqh3[s], hsqh4[s], rs13, rs14, rs34,
                            p134_111, p134_110, p134_101, p134_011, p134_100, p134_010, p134_001, p134_000, sparse[c(1,3,4)])
      bjh1[indmat[snp.index$snp.index134,1]] = res134[[1]][1,]; bjh3[indmat[snp.index$snp.index134,3]] = res134[[1]][2,]; bjh4[indmat[snp.index$snp.index134,4]] = res134[[1]][3,]
      res135 = update_beta3(snp.index$snp.index135, indmat[,c(1,3,5)], C1.sfbm, C3.sfbm, C5.sfbm,
                            Bh[[1]], Bh[[3]], Bh[[5]],  bjh1, bjh3, bjh5,
                            sigmasq1, sigmasq3, sigmasq5, hsqh1[s], hsqh3[s], hsqh5[s], rs13, rs15, rs35,
                            p135_111, p135_110, p135_101, p135_011, p135_100, p135_010, p135_001, p135_000, sparse[c(1,3,5)])
      bjh1[indmat[snp.index$snp.index135,1]] = res135[[1]][1,]; bjh3[indmat[snp.index$snp.index135,3]] = res135[[1]][2,]; bjh5[indmat[snp.index$snp.index135,5]] = res135[[1]][3,]
      res145 = update_beta3(snp.index$snp.index145, indmat[,c(1,4,5)], C1.sfbm, C4.sfbm, C5.sfbm,
                            Bh[[1]], Bh[[4]], Bh[[5]],  bjh1, bjh4, bjh5,
                            sigmasq1, sigmasq4, sigmasq5, hsqh1[s], hsqh4[s], hsqh5[s], rs14, rs15, rs45,
                            p145_111, p145_110, p145_101, p145_011, p145_100, p145_010, p145_001, p145_000, sparse[c(1,4,5)])
      bjh1[indmat[snp.index$snp.index145,1]] = res145[[1]][1,]; bjh4[indmat[snp.index$snp.index145,4]] = res145[[1]][2,]; bjh5[indmat[snp.index$snp.index145,5]] = res145[[1]][3,]
      res234 = update_beta3(snp.index$snp.index234, indmat[,c(2,3,4)], C2.sfbm, C3.sfbm, C4.sfbm,
                            Bh[[2]], Bh[[3]], Bh[[4]],  bjh2, bjh3, bjh4,
                            sigmasq2, sigmasq3, sigmasq4, hsqh2[s], hsqh3[s], hsqh4[s], rs23, rs24, rs34,
                            p234_111, p234_110, p234_101, p234_011, p234_100, p234_010, p234_001, p234_000, sparse[c(2,3,4)])
      bjh2[indmat[snp.index$snp.index234,2]] = res234[[1]][1,]; bjh3[indmat[snp.index$snp.index234,3]] = res234[[1]][2,]; bjh4[indmat[snp.index$snp.index234,4]] = res234[[1]][3,]
      res235 = update_beta3(snp.index$snp.index235, indmat[,c(2,3,5)], C2.sfbm, C3.sfbm, C5.sfbm,
                            Bh[[2]], Bh[[3]], Bh[[5]],  bjh2, bjh3, bjh5,
                            sigmasq2, sigmasq3, sigmasq5, hsqh2[s], hsqh3[s], hsqh5[s], rs23, rs25, rs35,
                            p235_111, p235_110, p235_101, p235_011, p235_100, p235_010, p235_001, p235_000, sparse[c(2,3,5)])
      bjh2[indmat[snp.index$snp.index235,2]] = res235[[1]][1,]; bjh3[indmat[snp.index$snp.index235,3]] = res235[[1]][2,]; bjh5[indmat[snp.index$snp.index235,5]] = res235[[1]][3,]
      res245 = update_beta3(snp.index$snp.index245, indmat[,c(2,4,5)], C2.sfbm, C4.sfbm, C5.sfbm,
                            Bh[[2]], Bh[[4]], Bh[[5]],  bjh2, bjh4, bjh5,
                            sigmasq2, sigmasq4, sigmasq5, hsqh2[s], hsqh4[s], hsqh5[s], rs24, rs25, rs45,
                            p245_111, p245_110, p245_101, p245_011, p245_100, p245_010, p245_001, p245_000, sparse[c(2,4,5)])
      bjh2[indmat[snp.index$snp.index245,2]] = res245[[1]][1,]; bjh4[indmat[snp.index$snp.index245,4]] = res245[[1]][2,]; bjh5[indmat[snp.index$snp.index245,5]] = res245[[1]][3,]
      res345 = update_beta3(snp.index$snp.index345, indmat[,c(3,4,5)], C3.sfbm, C4.sfbm, C5.sfbm,
                            Bh[[3]], Bh[[4]], Bh[[5]],  bjh3, bjh4, bjh5,
                            sigmasq3, sigmasq4, sigmasq5, hsqh3[s], hsqh4[s], hsqh5[s], rs34, rs35, rs45,
                            p345_111, p345_110, p345_101, p345_011, p345_100, p345_010, p345_001, p345_000, sparse[c(3,4,5)])
      bjh3[indmat[snp.index$snp.index345,3]] = res345[[1]][1,]; bjh4[indmat[snp.index$snp.index345,4]] = res345[[1]][2,]; bjh5[indmat[snp.index$snp.index345,5]] = res345[[1]][3,]
      
      res12 = update_beta2(snp.index$snp.index12, indmat[,c(1,2)], C1.sfbm, C2.sfbm, Bh[[1]], Bh[[2]],
                           bjh1, bjh2, sigmasq1, sigmasq2, hsqh1[s], hsqh2[s], rs12, p_12_11, p_12_10, p_12_01, sparse[c(1,2)])
      bjh1[indmat[snp.index$snp.index12,1]] = res12[[1]][1,]; bjh2[indmat[snp.index$snp.index12,2]] = res12[[1]][2,];
      res13 = update_beta2(snp.index$snp.index13, indmat[,c(1,3)], C1.sfbm, C3.sfbm, Bh[[1]], Bh[[3]],
                           bjh1, bjh3, sigmasq1, sigmasq3, hsqh1[s], hsqh3[s], rs13, p_13_11, p_13_10, p_13_01, sparse[c(1,3)])
      bjh1[indmat[snp.index$snp.index13,1]] = res13[[1]][1,]; bjh3[indmat[snp.index$snp.index13,3]] = res13[[1]][2,];
      res14 = update_beta2(snp.index$snp.index14, indmat[,c(1,4)], C1.sfbm, C4.sfbm, Bh[[1]], Bh[[4]],
                           bjh1, bjh4, sigmasq1, sigmasq4, hsqh1[s], hsqh4[s], rs14, p_14_11, p_14_10, p_14_01, sparse[c(1,4)])
      bjh1[indmat[snp.index$snp.index14,1]] = res14[[1]][1,]; bjh4[indmat[snp.index$snp.index14,4]] = res14[[1]][2,];
      res15 = update_beta2(snp.index$snp.index15, indmat[,c(1,5)], C1.sfbm, C5.sfbm, Bh[[1]], Bh[[5]],
                           bjh1, bjh5, sigmasq1, sigmasq5, hsqh1[s], hsqh5[s], rs15, p_15_11, p_15_10, p_15_01, sparse[c(1,5)])
      bjh1[indmat[snp.index$snp.index15,1]] = res15[[1]][1,]; bjh5[indmat[snp.index$snp.index15,5]] = res15[[1]][2,];
      res23 = update_beta2(snp.index$snp.index23, indmat[,c(2,3)], C2.sfbm, C3.sfbm, Bh[[2]], Bh[[3]],
                           bjh2, bjh3, sigmasq2, sigmasq3, hsqh2[s], hsqh3[s], rs23, p_23_11, p_23_10, p_23_01, sparse[c(2,3)])
      bjh2[indmat[snp.index$snp.index23,2]] = res23[[1]][1,]; bjh3[indmat[snp.index$snp.index23,3]] = res23[[1]][2,];
      res24 = update_beta2(snp.index$snp.index24, indmat[,c(2,4)], C2.sfbm, C4.sfbm, Bh[[2]], Bh[[4]],
                           bjh2, bjh4, sigmasq2, sigmasq4, hsqh2[s], hsqh4[s], rs24, p_24_11, p_24_10, p_24_01, sparse[c(2,4)])
      bjh2[indmat[snp.index$snp.index24,2]] = res24[[1]][1,]; bjh4[indmat[snp.index$snp.index24,4]] = res24[[1]][2,];
      res25 = update_beta2(snp.index$snp.index25, indmat[,c(2,5)], C2.sfbm, C5.sfbm, Bh[[2]], Bh[[5]],
                           bjh2, bjh5, sigmasq2, sigmasq5, hsqh2[s], hsqh5[s], rs25, p_25_11, p_25_10, p_25_01, sparse[c(2,5)])
      bjh2[indmat[snp.index$snp.index25,2]] = res25[[1]][1,]; bjh5[indmat[snp.index$snp.index25,5]] = res25[[1]][2,];
      res34 = update_beta2(snp.index$snp.index34, indmat[,c(3,4)], C3.sfbm, C4.sfbm, Bh[[3]], Bh[[4]],
                           bjh3, bjh4, sigmasq3, sigmasq4, hsqh3[s], hsqh4[s], rs34, p_34_11, p_34_10, p_34_01, sparse[c(3,4)])
      bjh3[indmat[snp.index$snp.index34,3]] = res34[[1]][1,]; bjh4[indmat[snp.index$snp.index34,4]] = res34[[1]][2,];
      res35 = update_beta2(snp.index$snp.index35, indmat[,c(3,5)], C3.sfbm, C5.sfbm, Bh[[3]], Bh[[5]],
                           bjh3, bjh5, sigmasq3, sigmasq5, hsqh3[s], hsqh5[s], rs35, p_35_11, p_35_10, p_35_01, sparse[c(3,5)])
      bjh3[indmat[snp.index$snp.index35,3]] = res35[[1]][1,]; bjh5[indmat[snp.index$snp.index35,5]] = res35[[1]][2,];
      res45 = update_beta2(snp.index$snp.index45, indmat[,c(4,5)], C4.sfbm, C5.sfbm, Bh[[4]], Bh[[5]],
                           bjh4, bjh5, sigmasq4, sigmasq5, hsqh4[s], hsqh5[s], rs45, p_45_11, p_45_10, p_45_01, sparse[c(4,5)])
      bjh4[indmat[snp.index$snp.index45,4]] = res45[[1]][1,]; bjh5[indmat[snp.index$snp.index45,5]] = res45[[1]][2,];
      
      res1 = update_beta1(indmat[snp.index$snp.index1,1], C1.sfbm, Bh[[1]], bjh1, sigmasq1, hsqh1[s], p.causal1, sparse[1])
      bjh1[indmat[snp.index$snp.index1,1]] = res1[1,]; 
      res2 = update_beta1(indmat[snp.index$snp.index2,2], C2.sfbm, Bh[[2]], bjh2, sigmasq2, hsqh2[s], p.causal2, sparse[2])
      bjh2[indmat[snp.index$snp.index2,2]] = res2[1,];
      res3 = update_beta1(indmat[snp.index$snp.index3,3], C3.sfbm, Bh[[3]], bjh3, sigmasq3, hsqh3[s], p.causal3, sparse[3])
      bjh3[indmat[snp.index$snp.index3,3]] = res3[1,]; 
      res4 = update_beta1(indmat[snp.index$snp.index4,4], C4.sfbm, Bh[[4]], bjh4, sigmasq4, hsqh4[s], p.causal4, sparse[4])
      bjh4[indmat[snp.index$snp.index4,4]] = res4[1,];
      res5 = update_beta1(indmat[snp.index$snp.index5,5], C5.sfbm, Bh[[5]], bjh5, sigmasq5, hsqh5[s], p.causal5, sparse[5])
      bjh5[indmat[snp.index$snp.index5,5]] = res5[1,];
      
      if (s>n.burnin){
        pm1[indmat[snp.index$snp.index12345,1]] = pm1[indmat[snp.index$snp.index12345,1]] + res12345[[2]][1,];
        pm2[indmat[snp.index$snp.index12345,2]] = pm2[indmat[snp.index$snp.index12345,2]] + res12345[[2]][2,];
        pm3[indmat[snp.index$snp.index12345,3]] = pm3[indmat[snp.index$snp.index12345,3]] + res12345[[2]][3,];
        pm4[indmat[snp.index$snp.index12345,4]] = pm4[indmat[snp.index$snp.index12345,4]] + res12345[[2]][4,];
        pm5[indmat[snp.index$snp.index12345,5]] = pm5[indmat[snp.index$snp.index12345,5]] + res12345[[2]][5,];
        
        pm1[indmat[snp.index$snp.index1234,1]] = pm1[indmat[snp.index$snp.index1234,1]] + res1234[[2]][1,];
        pm2[indmat[snp.index$snp.index1234,2]] = pm2[indmat[snp.index$snp.index1234,2]] + res1234[[2]][2,];
        pm3[indmat[snp.index$snp.index1234,3]] = pm3[indmat[snp.index$snp.index1234,3]] + res1234[[2]][3,];
        pm4[indmat[snp.index$snp.index1234,4]] = pm4[indmat[snp.index$snp.index1234,4]] + res1234[[2]][4,];
        
        pm1[indmat[snp.index$snp.index1235,1]] = pm1[indmat[snp.index$snp.index1235,1]] + res1235[[2]][1,];
        pm2[indmat[snp.index$snp.index1235,2]] = pm2[indmat[snp.index$snp.index1235,2]] + res1235[[2]][2,];
        pm3[indmat[snp.index$snp.index1235,3]] = pm3[indmat[snp.index$snp.index1235,3]] + res1235[[2]][3,];
        pm5[indmat[snp.index$snp.index1235,5]] = pm5[indmat[snp.index$snp.index1235,5]] + res1235[[2]][4,];
        
        pm1[indmat[snp.index$snp.index1245,1]] = pm1[indmat[snp.index$snp.index1245,1]] + res1245[[2]][1,];
        pm2[indmat[snp.index$snp.index1245,2]] = pm2[indmat[snp.index$snp.index1245,2]] + res1245[[2]][2,];
        pm4[indmat[snp.index$snp.index1245,4]] = pm4[indmat[snp.index$snp.index1245,4]] + res1245[[2]][3,];
        pm5[indmat[snp.index$snp.index1245,5]] = pm5[indmat[snp.index$snp.index1245,5]] + res1245[[2]][4,];
        
        pm1[indmat[snp.index$snp.index1345,1]] = pm1[indmat[snp.index$snp.index1345,1]] + res1345[[2]][1,];
        pm3[indmat[snp.index$snp.index1345,3]] = pm3[indmat[snp.index$snp.index1345,3]] + res1345[[2]][2,];
        pm4[indmat[snp.index$snp.index1345,4]] = pm4[indmat[snp.index$snp.index1345,4]] + res1345[[2]][3,];
        pm5[indmat[snp.index$snp.index1345,5]] = pm5[indmat[snp.index$snp.index1345,5]] + res1345[[2]][4,];
        
        pm2[indmat[snp.index$snp.index2345,2]] = pm2[indmat[snp.index$snp.index2345,2]] + res2345[[2]][1,];
        pm3[indmat[snp.index$snp.index2345,3]] = pm3[indmat[snp.index$snp.index2345,3]] + res2345[[2]][2,];
        pm4[indmat[snp.index$snp.index2345,4]] = pm4[indmat[snp.index$snp.index2345,4]] + res2345[[2]][3,];
        pm5[indmat[snp.index$snp.index2345,5]] = pm5[indmat[snp.index$snp.index2345,5]] + res2345[[2]][4,];
        
        pm1[indmat[snp.index$snp.index123,1]] = pm1[indmat[snp.index$snp.index123,1]] + res123[[2]][1,];
        pm2[indmat[snp.index$snp.index123,2]] = pm2[indmat[snp.index$snp.index123,2]] + res123[[2]][2,]
        pm3[indmat[snp.index$snp.index123,3]] = pm3[indmat[snp.index$snp.index123,3]] + res123[[2]][3,]
        pm1[indmat[snp.index$snp.index124,1]] = pm1[indmat[snp.index$snp.index124,1]] + res124[[2]][1,];
        pm2[indmat[snp.index$snp.index124,2]] = pm2[indmat[snp.index$snp.index124,2]] + res124[[2]][2,]
        pm4[indmat[snp.index$snp.index124,4]] = pm4[indmat[snp.index$snp.index124,4]] + res124[[2]][3,]
        pm1[indmat[snp.index$snp.index125,1]] = pm1[indmat[snp.index$snp.index125,1]] + res125[[2]][1,];
        pm2[indmat[snp.index$snp.index125,2]] = pm2[indmat[snp.index$snp.index125,2]] + res125[[2]][2,]
        pm5[indmat[snp.index$snp.index125,5]] = pm5[indmat[snp.index$snp.index125,5]] + res125[[2]][3,]
        
        pm1[indmat[snp.index$snp.index134,1]] = pm1[indmat[snp.index$snp.index134,1]] + res134[[2]][1,];
        pm3[indmat[snp.index$snp.index134,3]] = pm3[indmat[snp.index$snp.index134,3]] + res134[[2]][2,]
        pm4[indmat[snp.index$snp.index134,4]] = pm4[indmat[snp.index$snp.index134,4]] + res134[[2]][3,]
        pm1[indmat[snp.index$snp.index135,1]] = pm1[indmat[snp.index$snp.index135,1]] + res135[[2]][1,];
        pm3[indmat[snp.index$snp.index135,3]] = pm3[indmat[snp.index$snp.index135,3]] + res135[[2]][2,]
        pm5[indmat[snp.index$snp.index135,5]] = pm5[indmat[snp.index$snp.index135,5]] + res135[[2]][3,]
        pm1[indmat[snp.index$snp.index145,1]] = pm1[indmat[snp.index$snp.index145,1]] + res145[[2]][1,];
        pm4[indmat[snp.index$snp.index145,4]] = pm4[indmat[snp.index$snp.index145,4]] + res145[[2]][2,]
        pm5[indmat[snp.index$snp.index145,5]] = pm5[indmat[snp.index$snp.index145,5]] + res145[[2]][3,]
        
        pm2[indmat[snp.index$snp.index234,2]] = pm2[indmat[snp.index$snp.index234,2]] + res234[[2]][1,];
        pm3[indmat[snp.index$snp.index234,3]] = pm3[indmat[snp.index$snp.index234,3]] + res234[[2]][2,]
        pm4[indmat[snp.index$snp.index234,4]] = pm4[indmat[snp.index$snp.index234,4]] + res234[[2]][3,]
        pm2[indmat[snp.index$snp.index235,2]] = pm2[indmat[snp.index$snp.index235,2]] + res235[[2]][1,];
        pm3[indmat[snp.index$snp.index235,3]] = pm3[indmat[snp.index$snp.index235,3]] + res235[[2]][2,]
        pm5[indmat[snp.index$snp.index235,5]] = pm5[indmat[snp.index$snp.index235,5]] + res235[[2]][3,]
        pm2[indmat[snp.index$snp.index245,2]] = pm2[indmat[snp.index$snp.index245,2]] + res245[[2]][1,];
        pm4[indmat[snp.index$snp.index245,4]] = pm4[indmat[snp.index$snp.index245,4]] + res245[[2]][2,]
        pm5[indmat[snp.index$snp.index245,5]] = pm5[indmat[snp.index$snp.index245,5]] + res245[[2]][3,]
        pm3[indmat[snp.index$snp.index345,3]] = pm3[indmat[snp.index$snp.index345,3]] + res345[[2]][1,];
        pm4[indmat[snp.index$snp.index345,4]] = pm4[indmat[snp.index$snp.index345,4]] + res345[[2]][2,]
        pm5[indmat[snp.index$snp.index345,5]] = pm5[indmat[snp.index$snp.index345,5]] + res345[[2]][3,]
        
        pm1[indmat[snp.index$snp.index12,1]] = pm1[indmat[snp.index$snp.index12,1]] + res12[[2]][1,]; pm2[indmat[snp.index$snp.index12,2]] = pm2[indmat[snp.index$snp.index12,2]] + res12[[2]][2,]
        pm1[indmat[snp.index$snp.index13,1]] = pm1[indmat[snp.index$snp.index13,1]] + res13[[2]][1,]; pm3[indmat[snp.index$snp.index13,3]] = pm3[indmat[snp.index$snp.index13,3]] + res13[[2]][2,]
        pm1[indmat[snp.index$snp.index14,1]] = pm1[indmat[snp.index$snp.index14,1]] + res14[[2]][1,]; pm4[indmat[snp.index$snp.index14,4]] = pm4[indmat[snp.index$snp.index14,4]] + res14[[2]][2,]
        pm1[indmat[snp.index$snp.index15,1]] = pm1[indmat[snp.index$snp.index15,1]] + res15[[2]][1,]; pm5[indmat[snp.index$snp.index15,5]] = pm5[indmat[snp.index$snp.index15,5]] + res15[[2]][2,]
        pm2[indmat[snp.index$snp.index23,2]] = pm2[indmat[snp.index$snp.index23,2]] + res23[[2]][1,]; pm3[indmat[snp.index$snp.index23,3]] = pm3[indmat[snp.index$snp.index23,3]] + res23[[2]][2,]
        pm2[indmat[snp.index$snp.index24,2]] = pm2[indmat[snp.index$snp.index24,2]] + res24[[2]][1,]; pm4[indmat[snp.index$snp.index24,4]] = pm4[indmat[snp.index$snp.index24,4]] + res24[[2]][2,]
        pm2[indmat[snp.index$snp.index25,2]] = pm2[indmat[snp.index$snp.index25,2]] + res25[[2]][1,]; pm5[indmat[snp.index$snp.index25,5]] = pm5[indmat[snp.index$snp.index25,5]] + res25[[2]][2,]
        pm3[indmat[snp.index$snp.index34,3]] = pm3[indmat[snp.index$snp.index34,3]] + res34[[2]][1,]; pm4[indmat[snp.index$snp.index34,4]] = pm4[indmat[snp.index$snp.index34,4]] + res34[[2]][2,]
        pm3[indmat[snp.index$snp.index35,3]] = pm3[indmat[snp.index$snp.index35,3]] + res35[[2]][1,]; pm5[indmat[snp.index$snp.index35,5]] = pm5[indmat[snp.index$snp.index35,5]] + res35[[2]][2,]
        pm4[indmat[snp.index$snp.index45,4]] = pm4[indmat[snp.index$snp.index45,4]] + res45[[2]][1,]; pm5[indmat[snp.index$snp.index45,5]] = pm5[indmat[snp.index$snp.index45,5]] + res45[[2]][2,]
        
        pm1[indmat[snp.index$snp.index1,1]] = pm1[indmat[snp.index$snp.index1,1]] + res1[2,]; pm2[indmat[snp.index$snp.index2,2]] = pm2[indmat[snp.index$snp.index2,2]] + res2[2,]
        pm3[indmat[snp.index$snp.index3,3]] = pm3[indmat[snp.index$snp.index3,3]] + res3[2,]; pm4[indmat[snp.index$snp.index4,4]] = pm4[indmat[snp.index$snp.index4,4]] + res4[2,]
        pm5[indmat[snp.index$snp.index5,5]] = pm5[indmat[snp.index$snp.index5,5]] + res5[2,]
      }
      progBar(s, niter, per=20)
    }
    pm1 = pm1/(niter-n.burnin); pm2 = pm2/(niter-n.burnin); pm3 = pm3/(niter-n.burnin); pm4 = pm4/(niter-n.burnin); pm5 = pm5/(niter-n.burnin)
    names(pm1) = snpinfo$rsid[tem[[1]]]; names(pm2) = snpinfo$rsid[tem[[2]]]
    names(pm3) = snpinfo$rsid[tem[[3]]]; names(pm4) = snpinfo$rsid[tem[[4]]]
    names(pm5) = snpinfo$rsid[tem[[5]]]
    
    return(list(pm1 = pm1, pm2 = pm2, pm3 = pm3, pm4 = pm4, pm5 = pm5))
  }
}














