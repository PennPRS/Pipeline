#!/bin/bash

# Load modules
module load r

PennPRS_path=$1
homedir=$2
input_GWAS_path=$3
submissionID=$4
methods=$5
trait=$6
race=$7
NCORES=$8

chmod 777 ${PennPRS_path}/software/DBSLMM/dbslmm 
Rscript ${PennPRS_path}/code/Tuning-Parameter-Free.R \
--PennPRS_path ${PennPRS_path} \
--homedir ${homedir} \
--input_GWAS_path ${input_GWAS_path} \
--submissionID ${submissionID} \
--methods ${methods} \
--trait ${trait} \
--race ${race} \
--NCORES ${NCORES} \
