#!/bin/bash

# Load modules
module load r
module load anaconda

PennPRS_path=$1
homedir=$2
input_GWAS_path=$3
submissionID=$4
trait=$5
race=$6
N_THREADS=$7
type=$8
phi=$9

export MKL_NUM_THREADS=$N_THREADS
export NUMEXPR_NUM_THREADS=$N_THREADS
export OMP_NUM_THREADS=$N_THREADS

Rscript ${PennPRS_path}/code/PRS-CS.R \
--PennPRS_path ${PennPRS_path} \
--homedir ${homedir} \
--input_GWAS_path ${input_GWAS_path} \
--submissionID ${submissionID} \
--trait ${trait} \
--race ${race} \
--N_THREADS ${N_THREADS} \
--type ${type} \
--phi ${phi} \