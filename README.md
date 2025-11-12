# PennPRS Offline Pipeline 

[PennPRS](https://pennprs.org/) is a cloud-based platform dedicated to online PRS model training without requiring individual-level data for parameter optimization. On this Github page, we provide the offline version of PennPRS that can be downloaded to local servers for large-scale PRS training with various method options. 

To use the tool, please follow the instructions in **[the Wiki page](https://github.com/PennPRS/Pipeline/wiki)**.
</br>



## Version History
- [ ] __November 2025:__  Updated code/Tuning-Parameter-Free.R.
- [ ] __January 2025:__  The PennPRS offline pipeline was made available on Github.
</br>



## Getting Started

To install the PennPRS Offline Pipeline, please clone the Github repository by `git clone https://github.com/PennPRS/Pipeline.git` and rename the unzipped folder as `/PennPRS/`.

Download LD reference data files for different populations and save the uncompressed folder(s) in `/PennPRS/LD/`.

[EUR LD information](https://www.dropbox.com/scl/fi/h3sv4l0wh36ki2lrmektl/EUR.tar.gz?rlkey=4ndd32swtbx1uo2awjv79a9mm&st=t3169p1q&dl=0) (~35.02G), decompress by `tar -zxvf EUR.tar.gz`

[AFR LD information](https://www.dropbox.com/scl/fi/ljmyncadxpehnx7j1scli/AFR.tar.gz?rlkey=13bb3qer2zt7s95cb377yexd7&st=tbjcnf4a&dl=0) (~44.12G), decompress by `tar -zxvf AFR.tar.gz`

[AMR LD information](https://www.dropbox.com/scl/fi/8f2i8l7f49tuarpfmsmzq/AMR.tar.gz?rlkey=lgxm7gr5sekedqx7ku1sw3yg0&st=7rvphfcj&dl=0) (~40.17G), decompress by `tar -zxvf AMR.tar.gz`

[EAS LD information](https://www.dropbox.com/scl/fi/3zg2zdv9o8txmhbj2vzdb/EAS.tar.gz?rlkey=56r4dieiqzu52knnjlphxkjll&st=wchdvytz&dl=0) (~25.68G), decompress by `tar -zxvf EAS.tar.gz`

[SAS LD information](https://www.dropbox.com/scl/fi/ki5ar39uzfgbqjor5hy1b/SAS.tar.gz?rlkey=3fcqio7n4w1lmr7c52wjue4ua&st=4e6uzvua&dl=0) (~30.56G), decompress by `tar -zxvf SAS.tar.gz`


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

[Single-Ancestry PRS Modeling](https://github.com/PennPRS/Pipeline/wiki/2.-Single%E2%80%90Ancestry-PRS-Modeling)
  1. C+T-pseudo
  2. Lassosum2-pseudo
  3. LDpred2-pseudo
  4. PRS-CS-pseudo
  5. PRS-CS-auto
  6. LDpred2-auto
  7. DBSLMM
  
  (We also provide the option to train an ensemble PRS combining PRS trained by a subset of C+T-pseudo, Lassosum2-pseudo, and LDpred2-pseudo)
  
[Multi-Ancestry PRS Modeling](https://github.com/PennPRS/Pipeline/wiki/3.-Multi%E2%80%90Ancestry-PRS-Modeling-with-Pseudo%E2%80%90Training-Methods)

  8. PROSPER-pseudo 
  9. MUSSEL-pseudo 
  10. PRS-CSx-pseudo 

## Demo and Runtime Information
We have provided example GWAS summary datasets and the corresponding outputs can be found in Sections 2.1 - 2.4 in **[the Wiki page](https://github.com/PennPRS/Pipeline/wiki)**.
The average run time for completing a job that runs C+T-pseudo, Lassosum2-pseudo, LDpred2-pseudo, and ensemble PRS for ~1.2 million HapMap3 SNPs using 2 CPUs (with 30 GB RAM) is approximately 2.5 hours, while increasing to 4 CPUs reduced the run time to approximately two hours.

## Contact
Please report questions and bugs on the Issues page or contact us at pennprs@googlegroups.com.


## Citation
Jin, J., Li, B., Wang, X., Yang, X., Li, Y., Wang, R., Ye, C., Shu, J., Fan, Z., Xue, F. and Ge, T., 2025. PennPRS: a centralized cloud computing platform for efficient polygenic risk score training in precision medicine. medRxiv, 2025-02. [Link](https://www.medrxiv.org/content/10.1101/2025.02.07.25321875v1)


