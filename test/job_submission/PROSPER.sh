#!/bin/bash

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

# Load modules
module load r

PennPRS_path=$1
homedir=$2
input_GWAS_path=$3
submissionID=$4
methods=$5
trait=$6
races=$7
NCORES=$8

Rscript ${PennPRS_path}code/multi-ancestry-step1.R \
--PennPRS_path "${PennPRS_path}" \
--homedir "${homedir}" \
--input_GWAS_path "${input_GWAS_path}" \
--submissionID "${submissionID}" \
--methods "${methods}" \
--trait "${trait}" \
--races "${races}" \
--NCORES "${NCORES}" \

if [ $? -eq 0 ]; then
Rscript ${PennPRS_path}code/multi-ancestry-step2.R \
--PennPRS_path "${PennPRS_path}" \
--homedir "${homedir}" \
--submissionID "${submissionID}" \
--methods "${methods}" \
--trait "${trait}" \
--races "${races}" \
--NCORES "${NCORES}"
else
  # If the first script failed, print an error message and do not proceed
  echo "Error: First script failed. Second script will not run."
fi