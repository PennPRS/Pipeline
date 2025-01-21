# Simplified version of lassosum2

This is a command line tool based on R programming language for lassosum2. This is a simplified version of lassosum2 implemented in PROSPER. Please be careful if using it separately. Please contact Jingning Zhang (jingningzhang238@gmail.com) if you have any questions.


## Getting Started

The required scripts, LD reference panels, R packages, and inputs are same as [PROSPER]( https://github.com/Jingning-Zhang/PROSPER ). 

## Using lassosum2

A toy example is provided also in the last section of [README.md](https://github.com/Jingning-Zhang/PROSPER/blob/main/README.md) on the website of [PROSPER]( https://github.com/Jingning-Zhang/PROSPER).

Example codes for running lassosum2 in both EUR and AFR populations:

```
package='/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/try_from_github/PROSPER'
path_example='/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/try_from_github/PROSPER/example/'
path_result='/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/try_from_github/PROSPER/PROSPER_example_results/'
path_plink='/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2'

mkdir ${path_result}

Rscript ${package}/scripts/lassosum2.R \
--PATH_package ${package} \
--PATH_out ${path_result}/lassosum2 \
--PATH_plink ${path_plink} \
--FILE_sst ${path_example}/summdata/EUR.txt,${path_example}/summdata/AFR.txt \
--pop EUR,AFR \
--chrom 1-22 \
--bfile_tuning ${path_example}/sample_data/EUR/tuning_geno,${path_example}/sample_data/AFR/tuning_geno \
--pheno_tuning ${path_example}/sample_data/EUR/pheno.fam,${path_example}/sample_data/AFR/pheno.fam \
--bfile_testing ${path_example}/sample_data/EUR/testing_geno,${path_example}/sample_data/AFR/testing_geno \
--pheno_testing ${path_example}/sample_data/EUR/pheno.fam,${path_example}/sample_data/AFR/pheno.fam \
--testing TRUE \
--NCORES 5

```


<!-- 

```
directory = '/dcs04/nilanjan/data/jzhang2/PRS-epr'
target_pop = 'AFR'
path_out = '/dcs04/nilanjan/data/jzhang2/example/PRS-epr'
path_sumdata = '/dcs04/nilanjan/data/jzhang2/example/summdata'
path_lassosum2 = '/dcs04/nilanjan/data/jzhang2/example/lassosum2'
path_plink = '/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2'
path_geno = '/dcs04/nilanjan/data/jzhang2/UKBB/genotype'
path_pheno = '/dcs04/nilanjan/data/jzhang2/UKBB/phenotype'
path_covar = '/dcs04/nilanjan/data/jzhang2/UKBB/covariate'


Rscript ${directory}/scripts/PRS-epr.R \
--PATH_package ${directory} \
--PATH_out ${path_out} \
--FILE_sst ${path_sumdata}/EUR.txt,${path_sumdata}/AFR.txt,${path_sumdata}/AMR.txt \
--pop EUR,AFR,AMR \
--lassosum_param ${path_lassosum2}/EUR/optimal_param.txt,${path_lassosum2}/AFR/optimal_param.txt,${path_lassosum2}/AMR/optimal_param.txt \
--chrom 1-22 \
--Ll 10 \
--Lc 10 \
--verbose 1 \
--NCORES 22 

Rscript ${directory}/scripts/tuning_testing.R \
--PATH_plink ${path_plink} \
--PATH_out ${path_out} \
--prefix ${target_pop} \
--testing T \
--bfile_tuning ${path_geno}/${target_pop}/allchrom_tuning \
--pheno_tuning ${path_pheno}/${target_pop}/tuning.txt \
--covar_tuning ${path_covar}/${target_pop}/tuning.txt \
--bfile_testing ${path_geno}/${target_pop}/allchrom_testing \
--pheno_testing ${path_pheno}/${target_pop}/testing.txt \
--covar_testing ${path_covar}/${target_pop}/testing.txt \
--NCORES 50
```

 - **PATH_package** (required): Full path to the directory mentioned above that contains: 1. a folder scripts downloaded from github 2. LD reference panels downloaded from Google Drive.

 - **PATH_plink** (required): Full path to plink2 executable

 - **PATH_out** (required): Full path to the directory for output results

 - **FILE_sst** (required): Full path and the file name of the GWAS summary statistics from multiple populations, separated by comma. An example of the format
```
    rsid           chr     a1     a0    beta       beta_se    n_eff
    rs3131972      1       G      A     0.00235    0.01269    15000
    rs3131969      1       G      A     -0.01228   0.01206    15000
    rs1048488      1       T      C     0.00989    0.01162    15000

```
Required columns in the files:
  1. rsid: SNP ID in the format of rsXXXX.
  2. chr: chromosome number in the format of 1,2,...,22.
  3. beta: SNP effect. Note that for binary traits, beta is the coefficient in logistic regression, i.e. log(OR). 
  4. beta_se: Standard error of beta.
  5. a1: effective allele (counted allele in regression).
  6. a0: alternative allele (non-A1 allele).
  7. n_eff: Sample size per variant. Note that for binary traits, it is `effective sample sizes = 4 / (1 / N_control + 1 / N_case)`; and for continuous traits, it is simply the sample size.

Note that the summary statistics files are suggested to be cleaned as follows before using:
  1. Keep variants in reference panels to avoid troubles caused by reading huge files. The rsid of variants in reference panels can be found in `ref_bim.txt`.
  2. Alleles are suggested to match to reference panels to avoid flipping strands. The alleles of variants in reference panels can be found in `ref_bim.txt`.
  3. For each population, remove variants whose allele frequencies are less than 1%. If your summary statistics does not contain information of population-specific allele frequencies, we recommanded to use that in 1000G, which can be found in `ref_af.txt`.
  4. Remove variants with small sample size (90% of median sample size per variant).

 - **pop** (required): Population of the GWAS sample (AFR, AMR, EAS, EUR and SAS), in the same order of the GWAS summary statistics files, separated by comma.

 - **lassosum_param** (required): Full path and the file name of the optimal tuning parameters in lassosum2, separated by comma. The optimal parameters can be obtained by running the official [lassosum2]( https://privefl.github.io/bigsnpr-extdoc/polygenic-scores-pgs.html#example-ldpred2-and-lassosum2 "lassosum2"), or used a Simplified version [here]( https://drive.google.com/file/d/1aGwAGdIeVmTaTskSuUdEPE6Z9-evSlUG/view?usp=sharing "lassosum2") which is faster and requires less computation resource. The file must have the following format
```
    delta0    lambda0
    49.8      0.002
```
 - **chrom** (optional): The chromosome on which the model is fitted, separated by comma or dash for consecutive chromosomes. Default is 1-22.

 - **Ll** (optional): Length of tuning parameter path for the tuning parameter `lambda`. Default is 5.

 - **Lc** (optional): Length of tuning parameter path for the tuning parameter `c`. Default is 5.

 - **prefix** (required): Output filename prefix. Suggested value is name of target population. 

 - **bfile_tuning** (required): Full path to PLINK binary input file prefix (minus bed/bim/fam) for tuning.

 - **pheno_tuning** (optional): Full path to phenotype file for tuning (PLINK format). This is optional, taken from `bfile_tuning` otherwise.

 - **covar_tuning**  (optional): Full path to quantitative covariates for tuning (PLINK format).

 - **testing** (optional): Whether to perform testing in a seperate dataset. If TRUE, tuning datasets should be provided (bfile_testing, pheno_testing, covar_testing). Default is FALSE.

 - **bfile_tuning** (optional): Should be provided if `testing=TRUE`. Full path to PLINK binary input file prefix (minus bed/bim/fam) for testing.

 - **pheno_tuning** (optional): Should be provided if `testing=TRUE`. Full path to phenotype file for testing (PLINK format). This is optional, taken from `bfile_tuning` otherwise.

 - **covar_tuning**  (optional): Should be provided if `testing=TRUE`. Full path to quantitative covariates for testing (PLINK format).

 - **SL_library** (optional): Models input to SL.library in SuperLearner function, separated by comma. Default is "SL.glmnet,SL.ridge,SL.lm".

 - **linear_score** (optional): Whether to output PRS text file in the format of `rsid, a1, a0, weight`. Default is TRUE; otherwise, only superlearner model will be saved. Note some models available for SL.library are non-linear models, such as SL.nnet. In this case, PRS text file cannot be output.

 - **verbose** (optional): How much chatter to print: 0=nothing; 1=minimal; 2=all. Default is 1.

 - **NCORES** (optional): Number of cores used for parallel computation. Default is 1. Parallel is used on chromosomes in the script PRS-epr.R, and used on plink for PRS calculation in the script tuning_testing.R. So the suggested value is 22 for script PRS-epr.R, and as many cores as possible for tuning_testing.R. 


## Output -->


## Support

Please direct any problems or questions to Jingning Zhang (jingningzhang238@gmail.com).




