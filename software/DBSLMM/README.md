# DBSLMM
## Overview
There are two versions of DBSLMM: the tuning version and the deterministic version. The tuning version examines three different heritability choices and requires a validation data to tune the heritability hyper-parameter. The deterministic version uses one heritability estimate and directly fit the model in the training data without a separate validation data. Both versions requires a reference data to compute the SNP correlation matrix. <br>
For binary traits, especially for psychiatry diseases (i.e. major depression), DBSLMM usually outperforms some existing PGS construction methods.

## Relative papers in my group
### PGS tools
[Method paper](https://linkinghub.elsevier.com/retrieve/pii/S0002-9297(20)30109-9):  <em><strong>Yang S</strong></em>\#, Zhou X. Accurate and Scalable Construction of Polygenic Scores in Large Biobank Data Sets.  <em>Am J Hum Genet</em>. 2020 May 7;106(5):679-693. <br>
[Bechmarking paper](https://academic.oup.com/bib/article/23/2/bbac039/6534383?login=false): <em><strong>Yang S</strong></em>\#\*, Zhou X. PGS-server: accuracy, robustness and transferability of polygenic score methods for biobank scale studies. <em>Brief Bioinform</em>. 2022 Mar 10;23(2):bbac039. <br>
[Database paper](https://academic.oup.com/nar/article/52/D1/D963/7416385): Cao C, Zhang S, Wang J, Tian M, Ji X, Huang D\*, <em><strong>Yang S</strong></em>\*, Gu N. PGS-Depot: a comprehensive resource for polygenic scores constructed by summary statistics based methods. <em>Nucleic Acids Res</em>. 2024 Jan 5;52(D1):D963-D971. <br>
[Webserver paper](https://www.biorxiv.org/content/10.1101/2024.08.05.606619v1): <em><strong>Yang S</strong></em>\#\*, Ye X\#, Ji X\#, Li Z, Tian M, Huang P, Cao C. PGSFusion streamlines polygenic score construction and epidemiological applications in biobank-scale cohorts. <em>Genome Med</em>. 2025 Jul 14;17(1):77. <br>
### Applications
[Psycho-metabolic nexus](https://www.sciencedirect.com/science/article/pii/S2352396424005668): Guo X, Feng Y, Ji X, Jia N, Mainaiti A, Lai J, Wang Z\*, <em><strong>Yang S</strong></em>\*, Hu S\*. Shared genetic architecture and bidirectional clincial risks with psycho-metabolic nexus. <em>EBioMedicine</em>. 2025 Jan;111:105530. doi: 10.1016/j.ebiom.2024.105530.   <br>
[PTRS for stroke](https://link.springer.com/article/10.1007/s00439-024-02717-7): Cai X, Li H, Cao X, Ma X, Zhu W, Xu L, <em><strong>Yang S</em></strong>, Yu R, Huang P. Integrating transcriptomic and polygenic risk scores to enhance predictive accuracy for ischemic stroke subtypes. <em>Hum Genet</em>. 2024 Nov 18. doi: 10.1007/s00439-024-02717-7. 

## Update log
### v1.0 User-friendly DBSLMM
We update the `software/DBSLMM.R` and `software/TUNE.R`.
* integrates [<em>bigsnpr</em>](https://privefl.github.io/bigsnpr/) to estiamte the heritability, select large effect SNPs by clumping, and calculate the $\hat{y}$ for the validation set
* provides <em>'map_hm3_plus.rds'</em> and corresponding PLINK files for three ancestries from 1000 Genomes Project
* fits the whole genomoe for three DBSLMM models with one R script
* provides the online server [PGSFsuion](http://www.pgsfusion.net/) constructing PGS and performing epidemiological application in UKBB
### v0.3 Online Server
* provides the online version https://www.pgs-server.com/.
### v0.3
* fixes the bug of model fiting
### v0.21
* fits DBSLMM in different heritabiliies with one clumping
* uses a script to make the usage of DBSLMM-tuning version more easier.
### v0.2 
* validates the flods of heritability: 0.8, 1 and 1.2
* fits the LMM, when the chromosome without any large effect SNPs
* fits the external validation all by R code. 

## Tutorial for DBSLMM (v1.0)
In this version, we support DBSLMM and LMM models in automatic and tuning versions. 
### Prepartion
* install R packages: <em>bigsnpr, bigreadr, dplyr, plyr, optparse</em>
* download execution file [dbslmm](https://drive.google.com/file/d/1eAbEyhF8rO_faOFL3jqRo9LmfgJNRH6K/view?usp=sharing)
* download the reference panel of European, East Asian and African populations from [1000 Genome Project](https://drive.google.com/drive/folders/1tC5dT6f2otpY0iXMPRzIxfihERHyURr0?usp=sharing)
* ensure the summary statistics in [GEMMA format](https://github.com/genetics-statistics/GEMMA/blob/master/doc/manual.pdf)
### DBSLMM-auto
````bash
Summary_stat=/summary/statistics/in/GEMMA/fromat.txt
anc=EUR
BLOCK=${PACK_DIR}/block_data/${anc}/
outpath=${DATA_DIR}output/
model=DBSLMM
Rscript ${PACK_DIR}/software/DBSLMM.R --summary ${Summary_stat} --dbslmm ${PACK_DIR}/dbslmm --type auto --model ${model} \
                                      --reference ${ref_panel} --block ${BLOCK} --N ${N} --outPath ${outpath}/
````
<em><strong>Tips</strong></em>
* If you set `model=LMM`, we can fit DBSLMM-LMM model.
* If the summary statistics is not EUR, you can set `anc=AFR` or `anc=EAS`.
* If the trait is quantitative, you can set `N=n_sample`. 
* If the trait is binary, you can set `N=${n_case},${n_control}`. DBSLMM will use effective sample size to estimate heritability. You CAN also regard the binary trait as the quantitative trait. 


### DBSLMM-tuning
````bash
Summary_stat=/summary/statistics/in/GEMMA/fromat.txt
anc=EUR
BLOCK=${PACK_DIR}/block_data/${anc}
outpath=${DATA_DIR}output/
model=DBSLMM
Rscript ${PACK_DIR}/software/DBSLMM.R --summary ${Summary_stat} --dbslmm ${PACK_DIR}/dbslmm --type tuning --model DBSLMM
                                      --reference ${ref_panel} --block ${BLOCK} --N ${N} --outPath ${outpath}/ --h2f 0.8,1,1.2 		   
Summary_prefix=`echo "$Summary_stat" | awk -F'.txt' '{print $1}'`
Rscript ${PACK_DIR}/software/TUNE.R --dbslmm_eff ${Summary_prefix} --h2f 0.8,1,1.2 \
                                    --validation_g ${val_genotype} --validation_p ${val_phenotype} 
````
<em><strong>Tips</strong></em>
* `h2f` in `DBSLMM.R` and `TUNE.R` MUST be the same.
* If you have covarites, you can specify `cov` flag.

## Tutorial for external test only using summary statistics
### Make SNP correlation matrix
We should use reference panel to construct SNP correlation matrix. You can treat [1000 Genome Project](https://www.internationalgenome.org/data#download) as reference panel. We also need the [block information](http://bitbucket.org/nygcresearch/ldetect-data). The code is `MAKELD.R`. 
````bash
mkld=/your/path/MAKELD.R
# construct SNP corrlation matrix for each chromosome
for chr in `seq 1 22`
do
bfile=/your/path/chr${chr}
blockfile=/your/path/chr${chr}.bed
plink=/your/path/plink-1.9
outpath=/your/path/outpath/
Rscript ${mkld} --bfile ${bfile} --blockfile ${blockfile} --plink ${plink} --outpath ${outpath} --chr ${chr}
done
````
### Estimate R<sup>2
The format of `extsumm` is the same as that of `LDSC`. 
````bash
SNP N Z A1 A2
rs1983865 253135 3.79310344827586 T C
rs1983864 251364 -4.51612903225806 T G
rs12411954 253213 0.413793103448276 T C
rs7077266 250092 -0.92 T G
````
We should use the standardized SNP effect size. You have two options. First, the variant ID is read from column 1, an allele code is read from the following column, the standardized effect size associated with the named allele is read from the column after the allele column. E.g.
````bash
esteff=/your/path/esteff.txt,1,2,3
````
Second, the variant ID is read from column 1, an allele code is read from the following column, the standardized effect size associated with the named allele is read from the column after the allele column, the MAF with the named allele is read from the last column. E.g.
````bash
esteff=/your/path/esteff.txt,1,2,3,4
````
The example of `EXTERNAL.R` is as following:
````bash
external=/your/path/EXTERNAL.R
extsumm=/your/path/extsumm.ldsc.gz
esteff=/your/path/esteff.txt,1,2,3
LDpath=/your/path/LDpath/
outpath=/your/path/outpath/
Rscript ${external} --extsumm ${extsumm} --esteff ${esteff} --LDpath ${LDpath} --outpath ${outpath}
````
