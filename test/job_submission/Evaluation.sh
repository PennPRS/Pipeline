#!/bin/bash

# Load modules
module load r

PennPRS_path=$1
homedir=$2
phenofile=$3
bfile=$4
PRSdir=$5
ID_col_num=$6
pheno_col_num=$7
covar_col_nums=$8

Rscript ${PennPRS_path}/code/Evaluation_with_individual_level_dataset.R \
--PennPRS_path ${PennPRS_path} \
--homedir ${homedir} \
--phenofile ${phenofile} \
--bfile ${bfile} \
--PRSdir ${PRSdir} \
--ID_col_num ${ID_col_num} \
--pheno_col_num ${pheno_col_num} \
--covar_col_nums ${covar_col_nums} \