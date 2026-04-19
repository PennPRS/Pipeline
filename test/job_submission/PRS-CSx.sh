#!/bin/bash

# Load modules
module load r
module load anaconda

PennPRS_path=$1
homedir=$2
input_GWAS_path=$3
submissionID=$4
methods=$5
trait=$6
races=$7
N_THREADS=$8
phi=${9}

export MKL_NUM_THREADS=$N_THREADS
export NUMEXPR_NUM_THREADS=$N_THREADS
export OMP_NUM_THREADS=$N_THREADS

Rscript ${PennPRS_path}code/multi-ancestry-step1.R \
--PennPRS_path "${PennPRS_path}" \
--homedir "${homedir}" \
--input_GWAS_path "${input_GWAS_path}" \
--submissionID "${submissionID}" \
--methods "${methods}" \
--trait "${trait}" \
--races "${races}" \
--N_THREADS ${N_THREADS} \
--phi "${phi}" 

if [ $? -eq 0 ]; then
Rscript ${PennPRS_path}code/multi-ancestry-step2.R \
--PennPRS_path "${PennPRS_path}" \
--homedir "${homedir}" \
--submissionID "${submissionID}" \
--methods "${methods}" \
--trait "${trait}" \
--races "${races}" \
--N_THREADS ${N_THREADS} \
--phi "${phi}" 
else
  # If the first script failed, print an error message and do not proceed
  echo "Error: First script failed. Second script will not run."
fi