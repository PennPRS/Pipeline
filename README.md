# PennPRS Offline Pipeline 

[PennPRS](https://pennprs.org/) is a cloud-based platform dedicated to online PRS model training without requiring individual-level data for parameter optimization. On this Github page, we provide an offline version of PennPRS. 

To use the tool, please follow the instructions in **[the Wiki page](https://github.com/PennPRS/Pipeline/wiki)**.
</br>



## Version History
- [ ] __January, 2025:__  The PennPRS offline pipeline was made available on Github.
</br>



## Getting Started

To install the PennPRS Offline Pipeline, please clone the Github repository by `git clone https://github.com/PennPRS/Pipeline.git` and rename the unzipped folder as `/PennPRS/`.

Download LD reference data files for different populations and save the uncompressed folder(s) in `/PennPRS/LD/`.

[EUR LD information](https://www.dropbox.com/scl/fi/r3cwscxycfbgaxb4slh1d/EUR.tar.gz?rlkey=2um75zag5sgzpb82xbr504qhc&st=34n1hx12&dl=0) (~30.46G), decompress by `tar -zxvf EUR.tar.gz`

[AFR LD information](https://www.dropbox.com/scl/fi/zhlpeuaiqjbt1azx0r67h/AFR.tar.gz?rlkey=zflny7tra9bku3e24xehlcak3&st=f7it0p9o&dl=0) (~39.68G), decompress by `tar -zxvf AFR.tar.gz`

[AMR LD information](https://www.dropbox.com/scl/fi/54uxowqs5qhkbe9t776rc/AMR.tar.gz?rlkey=vqw0j78tyrqo6jgymiwevh9q9&st=vtqur2ws&dl=0) (~36.36G), decompress by `tar -zxvf AMR.tar.gz`

[EAS LD information](https://www.dropbox.com/scl/fi/s0a6mqpi14qdqvop871mg/EAS.tar.gz?rlkey=jeodg6upmbi2kifuk9iijrjvg&st=3djfh3fx&dl=0) (~23.20G), decompress by `tar -zxvf EAS.tar.gz`

[SAS LD information](https://www.dropbox.com/scl/fi/5b8937g2wb25q2gomvplr/SAS.tar.gz?rlkey=c9c6v7kadansbdee2xnyr297l&st=h3i4di6d&dl=0) (~24.94G), decompress by `tar -zxvf SAS.tar.gz`


Before running the pipeline, please consider the following quality control (QC) steps for the GWAS summary data:

- Only keep the biallelic [HapMap3 SNPs](https://www.dropbox.com/scl/fi/sktcg9u52jw1clvlj9qwx/hapmap3rsid.txt?rlkey=bwfqpqf9br4ptniee4wjd92c4&st=kefhjw6g&dl=0) to avoid troubles caused by reading huge files (e.g., > 8 million SNPs) in R.
- Remove SNPs with minor allele frequencies (MAF) lower than 1% in all populations.
- The genetic ancestry for each input GWAS summary data needs to be identified. If the GWAS training samples consist of multiple ancestry groups, please choose the ancestry group with the largest sample size.

## Notes

1. If you encounter errors regarding installing/loading R packages when running the pipeline, please manually install the following R packages first.

```
install.packages(c('RISCA','optparse','bigreadr','bigsnpr','bigparallelr', 'bigmemory','stringr','caret','scales','Rcpp', 'RcppArmadillo','RcppTN','inline','doMC','foreach','doParallel','data.table','readr','MASS','reshape','parallel',
'devtools','genio','dplyr','pryr','Matrix','lavaan'))
```

2. If PLINK or PLINK2 in `/PennPRS/software/` is not working, please follow the tutorials for [PLINK1.9](https://www.cog-genomics.org/plink/) and [PLINK2](https://www.cog-genomics.org/plink/2.0/) to re-install them under the same directory.
<be>

## PRS method options
PennPRS supports the following PRS pseudo-training and tuning-parameter-free methods. Please navigate to the Wiki page for the implementation of each of the methods.

- [Single-Ancestry PRS Modeling](https://github.com/PennPRS/Pipeline/wiki/2.-Single-Ancestry-PRS-Modeling)<summary>
    1. C+T-pseudo
    2. Lassosum2-pseudo
    3. LDpred2-pseudo
    4. Ensemble PRS combining PRS trained by different methods
    5. PRS-CS
    6. PRS-CS-auto
    7. LDpred2-auto
    8. DBSLMM
    
- [Multi-Ancestry PRS Modeling](https://github.com/PennPRS/Pipeline/wiki/3.-Multi-Ancestry-PRS-Modeling)<summary>

    9. PROSPER 
    10. MUSSEL 
    11. PRS-CSx 

## Contact
Please report questions and bugs on the Issues page or contact us at pennprs@googlegroups.com.


## Citation


