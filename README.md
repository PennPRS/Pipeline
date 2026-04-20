# PennPRS Offline Pipeline 

[PennPRS](https://pennprs.org/) is a cloud-based platform dedicated to online PRS model training without requiring individual-level data for parameter optimization. On this Github page, we provide the offline version of PennPRS that can be downloaded to local servers for large-scale PRS training with various method options. 

To use the tool, please follow the instructions in **[the Wiki page](https://github.com/PennPRS/Pipeline/wiki)**.
</br>

A test runner configuration is provided for users to easily test the pipeline in their local environment. To conduct an automated testing, please follow the instructions  
[here](./test-runner/).
</br>

## Version History
- [ ] __April 2026:__  Fixed bugs and made the following major updates: 

    * Added a new [mode](https://github.com/PennPRS/Pipeline/wiki/4.-Model-Evaluation-with-Individual%E2%80%90Level-Data) for evaluating trained PRS models based on user-provided individual-level data.<br>
    * Added a new [function](https://github.com/PennPRS/Pipeline?tab=readme-ov-file#query-data-with-pennprs), which allows users to directly query harmonized GWAS summary datasets from the GWAS Catalog or FinnGen.<br>
    * Added [Environment setup](https://github.com/PennPRS/Pipeline?tab=readme-ov-file#enviroment-set-up) to ensure all dependencies are available.<br>
    * Added a [test runner configuration](https://github.com/PennPRS/Pipeline/tree/main/test-runner) which allows users enable users to easily run and validate the pipeline in a local environment.<br>
    * Reorganized structure of the tutorial on the wiki page with [examples](https://github.com/PennPRS/Pipeline/wiki/4.-Test-Examples) for all supported methods.<br>
    * Updated the DBSLMM pipeline to the latest V1.0 Version, which requires downloading an additional LD folder (see [instructions](https://github.com/PennPRS/Pipeline?tab=readme-ov-file#getting-started)).<br>
    
- [ ] __November 2025:__  Updated code/Tuning-Parameter-Free.R.
- [ ] __January 2025:__  The PennPRS offline pipeline was made available on Github.
</br>



## Getting Started

To install the PennPRS Offline Pipeline, please clone the Github repository by `git clone https://github.com/PennPRS/Pipeline.git` and rename the unzipped folder as `/PennPRS/`.

Download LD reference data files for different populations and save the uncompressed folder(s) in `PennPRS/LD/`.

[EUR LD information](https://www.dropbox.com/scl/fi/h3sv4l0wh36ki2lrmektl/EUR.tar.gz?rlkey=4ndd32swtbx1uo2awjv79a9mm&st=t3169p1q&dl=0) (~32.62G), download by `curl -L -o EUR.tar.gz "https://www.dropbox.com/scl/fi/h3sv4l0wh36ki2lrmektl/EUR.tar.gz?rlkey=4ndd32swtbx1uo2awjv79a9mm&st=ui9t5mxo&dl=1"`, then decompress by `tar -zxvf EUR.tar.gz`

[AFR LD information](https://www.dropbox.com/scl/fi/ljmyncadxpehnx7j1scli/AFR.tar.gz?rlkey=13bb3qer2zt7s95cb377yexd7&st=tbjcnf4a&dl=0) (~41.09G), download by `curl -L -o AFR.tar.gz "https://www.dropbox.com/scl/fi/ljmyncadxpehnx7j1scli/AFR.tar.gz?rlkey=13bb3qer2zt7s95cb377yexd7&st=mmjuvsws&dl=1"`, then decompress by `tar -zxvf AFR.tar.gz`

[AMR LD information](https://www.dropbox.com/scl/fi/8f2i8l7f49tuarpfmsmzq/AMR.tar.gz?rlkey=lgxm7gr5sekedqx7ku1sw3yg0&st=7rvphfcj&dl=0) (~37.41G), download by `curl -L -o AMR.tar.gz "https://www.dropbox.com/scl/fi/8f2i8l7f49tuarpfmsmzq/AMR.tar.gz?rlkey=lgxm7gr5sekedqx7ku1sw3yg0&e=1&st=7rvphfcj&dl=1"`, then decompress by `tar -zxvf AMR.tar.gz`

[EAS LD information](https://www.dropbox.com/scl/fi/4er74cwemmaqh7cj796wr/EAS.tar.gz?rlkey=j3r5zsc421kizuri8cati9iuz&st=y2ou3lds&dl=0) (~25.66G), download by `curl -L -o EAS.tar.gz "https://www.dropbox.com/scl/fi/4er74cwemmaqh7cj796wr/EAS.tar.gz?rlkey=j3r5zsc421kizuri8cati9iuz&st=y2ou3lds&dl=1"`, then decompress by `tar -zxvf EAS.tar.gz`

[SAS LD information](https://www.dropbox.com/scl/fi/ki5ar39uzfgbqjor5hy1b/SAS.tar.gz?rlkey=3fcqio7n4w1lmr7c52wjue4ua&st=4e6uzvua&dl=0) (~28.46G), download by `curl -L -o SAS.tar.gz "https://www.dropbox.com/scl/fi/ki5ar39uzfgbqjor5hy1b/SAS.tar.gz?rlkey=3fcqio7n4w1lmr7c52wjue4ua&st=szbmp67j&dl=1"`, then decompress by `tar -zxvf SAS.tar.gz`

If you intend to run PRS-CSx, please download and save an additional [SNP information file](https://www.dropbox.com/scl/fi/0i74j1kpz24unfy82itsj/snpinfo_mult_1kg_hm3?rlkey=mhzcfm83v0jxdoemlsw8we716&e=1&dl=0) to `/PennPRS/LD/` by `curl -L -o snpinfo_mult_1kg_hm3 "https://www.dropbox.com/scl/fi/0i74j1kpz24unfy82itsj/snpinfo_mult_1kg_hm3?rlkey=mhzcfm83v0jxdoemlsw8we716&e=1&dl=1"` 

If you intend to run DBSLMM, please download the LD reference data for the ancestries of interest provided by the [DBSLMM team](https://drive.google.com/drive/folders/1tC5dT6f2otpY0iXMPRzIxfihERHyURr0) and save the folder(s) in `PennPRS/software/DBSLMM/LDref/`.

Before running the pipeline, please consider the following quality control (QC) steps for the input GWAS summary data:

- Only keep the biallelic [HapMap3 SNPs](https://www.dropbox.com/scl/fi/sktcg9u52jw1clvlj9qwx/hapmap3rsid.txt?rlkey=bwfqpqf9br4ptniee4wjd92c4&st=kefhjw6g&dl=0) to avoid troubles caused by reading huge files (e.g., > 8 million SNPs) in R.
- The genetic ancestry for each input GWAS summary data needs to be identified. If the GWAS training samples consist of multiple ancestry groups, please choose the ancestry group with the largest sample size.


## Enviroment Set Up

In `/PennPRS/`, create and activate the conda environment so all dependencies are available:

```bash
conda env create -f environment.yml
conda activate pennprs
```

## Notes

1. If you encounter errors regarding installing/loading R packages when running the pipeline, please manually install the following R packages first.

```
install.packages(c('pROC', 'readxl','optparse','bigreadr','bigsnpr','bigparallelr', 'bigmemory','stringr','caret','scales','Rcpp', 'RcppArmadillo','RcppTN','inline','doMC','foreach','doParallel','data.table','readr','MASS','reshape','parallel',
'devtools','genio','dplyr','pryr','Matrix','lavaan','BEDMatrix','ROCnReg'))
```
We did not set automatic installation of R packages because on servers or HPC systems, install.packages() may fail if you do not have write permission to the default R library. In that case, you may need to set a personal library path first and then install the R packages.

2. If PLINK or PLINK2 in `/PennPRS/software/` is not working, please follow the tutorials for [PLINK1.9](https://www.cog-genomics.org/plink/) and [PLINK2](https://www.cog-genomics.org/plink/2.0/) to re-install them under the same directory.
<be>


## PRS method options
PennPRS supports the following PRS pseudo-training and tuning-parameter-free methods. Please navigate to **[the Wiki page](https://github.com/PennPRS/Pipeline/wiki)** for the implementation of each method.

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


Each output folder contains the following contents:

1. README.txt
   List of contents in the output folder.
2. QC_report.txt - QC process for the input GWAS summary statistics file, including input data format check (detecting required columns) and standard QC steps applied to GWAS summary data for PRS training
3. PRS_INFO.txt - a summary report for PRS training, including:
   (1) method versions, tuning parameter settings, optimized tuning parameter values
   (2) example code for calculating PRS based on the trained PRS models using [PLINK2](https://www.cog-genomics.org/plink/2.0/score) or [pgsc_calc](https://pgsc-calc.readthedocs.io/en/latest/).
4. Trained PRS models:
   SNP weight files for the trained PRS models (`{ancestry}_{trait}_{method}.txt`)
   
   

## Query Data with PennPRS

We provide the option to directy query public, harmonized GWAS summary data files from the following two GWAS databases:
  1. [The GWAS Catalog](https://www.ebi.ac.uk/gwas/)
  2. [FinnGen](https://www.finngen.fi/en/access_results) (Note: FinnGen requires filling out an online form before downloading)

To query data for offline usage, run `query_data.py` to download the data file to your local server, then run PRS training pipelines with the queried data. 

```bash
module load anaconda
cd PennPRS
conda activate pennprs
python code/query_data.py <source> <trait_id> [options]
```

**Sources**

- `gwas` — GWAS Catalog: EBI summary statistics; trait ID (study accession, e.g., `GCST...`) must appear in the harmonised list (see [PennPRS data](https://pennprs.org/data)).
- `finngen` — FinnGen: R12 EUR files; path pattern `EUR_finngen_R12_<phenocode>.txt`. Phenocode must be in the queryable list (see [PennPRS data](https://pennprs.org/data)).

**Examples**

```
bash
module load anaconda
cd PennPRS
conda activate pennprs

# Resolve GWAS Catalog URL only (no download)
python code/query_data.py gwas GCST009979 --resolve-only

# Resolve and download GWAS Catalog file into a data dir (${gwas-data-dir})
python code/query_data.py gwas GCST009979 --download --gwas-data-dir=test/inputGWAS/

# Resolve FinnGen local path (prints path; exit 1 if file missing)
python code/query_data.py finngen F5_ALZHDEMENT
```

**CLI options**

| Option | Description |
|-----------|-------------|
| `--resolve-only` | Only print URL; do not download (GWAS only). |
| `--download` | For GWAS: download the file into the data dir. |
| `--output-dir DIR` | Directory for logs (e.g. `query_data.log`). Default: current directory. |
| `--harmonised-file PATH` | Path to `harmonised.txt` (GWAS only). Overrides env. |
| `--gwas-data-dir DIR` | Directory for saving GWAS summary data file downloaded from the GWAS Catalog. Overrides env. |
| `--finngen-data-dir DIR` | Base directory for saving GWAS summary data file downloaded from FinnGen. Overrides env. |

**Exit codes**

- `0` — Success (URL/path printed, or file downloaded).
- `1` — Trait not found or download failed.



## Evaluation of the Trained PRS models with Individual-Level Data 

We provide a pipeline for evaluating the performance of the trained PRS models on an individual-level dataset provided by the user. 

### Example

Prepare input files:

1. PRS model files. Save all the trained PRS models (not necessarily trained from the same job) in `${PRSdir}`. Do not change file names, keep the original file names as they are.

2. Individual-level data for evaluation purpose:

  (1) genotype data in PLINK format: `$PennPRS/test/Evaluation/geno/validation.{bim,bed,fam}` <br>
  (2) phenotype data: $PennPRS/test/Evaluation/PRSdir/pheno.txt.
    
    # Input Arguments
    PennPRS_path='PennPRS/'
    homedir="$PennPRS_path/test/Evaluation/output/"
    phenofile="$PennPRS_path/test/Evaluation/pheno/pheno.txt"
    bfile="$PennPRS_path/test/Evaluation/geno/validation"
    PRSdir="$PennPRS_path/test/Evaluation/PRSdir/"
    ID_col_num='1'
    pheno_col_num='2'
    covar_col_nums='3-44'
    
    # Job Submission
    sbatch test/job_submission/Evaluation.sh ${PennPRS_path} ${homedir} ${phenofile} ${bfile} ${PRSdir} ${ID_col_num} ${pheno_col_num} ${covar_col_nums}

**CLI options**

| Option | Description |
|--------|-------------|
| `‑‑PennPRS_path` | Path to PennPRS (required). |
| `‑‑homedir` | Folder where the output results are saved (required). |
| `‑‑phenofile` | Path to the individual-level phenotype data file for PRS model evaluation. The file can be in either .txt, .tsv, .csv, or .xlsx format, with required columns for individual ID and phenotype value, and optional columns for covariate information (required). |
| `‑‑bfile` | Path to the individual-level genotype data in PLINK format ([prefix].{bim,bed,fam}). Only prefix is needed (required). |
| `‑‑PRSdir` | Path to the folder the PRS model files are saved in. Do not change file names, keep the original file names generated by PennPRS as they are (required). |
| `‑‑ID_col_num` | Column number for individual ID (required). |
| `‑‑pheno_col_num` | Column number for phenotype value (required). |
| `‑‑covar_col_nums` | Column number(s) for covariate information (Optional). |

Please refer to the [tutorial](https://github.com/PennPRS/Pipeline/wiki/4.-Model-Evaluation-with-Individual%E2%80%90Level-Data#appendix) for details.


## Automated Testing

A SLURM-driven test harness for the worked examples in [Wiki § 5. Test Examples](https://github.com/PennPRS/Pipeline/wiki/5.-Test-Examples)
is provided under [`test-runner/`](./test-runner). The harness submits each example (5.1 – 5.7) via `sbatch`, monitors job completion, and verifies that
the expected output artifacts (per-method weight files, evaluation results, PLINK2 `.sscore` files) are produced and non-empty. This enables end-to-end
regression testing of the pipeline on any SLURM-managed cluster.

Typical usage:

```bash
cd test-runner
./run_all.sh --list               # enumerate configured tests
./run_all.sh --only 5.1a          # run a single example
./run_all.sh --parallel 4         # run all examples, up to 4 concurrently
```

Configuration (cluster paths, partitions, wall-time, memory) is centralized
in [`test-runner/tests.yaml`](./test-runner/tests.yaml). 

To conduct an automated pipeline testing, please go to 
[`test-runner/README.md`](./test-runner/README.md) for full setup
instructions, the complete test matrix, details of the verification
model, and test examples.






## Demo and Memory & Runtime Information
We have provided example GWAS summary datasets and the corresponding outputs can be found in Sections 2.1 - 2.4 in **[the Wiki page](https://github.com/PennPRS/Pipeline/wiki)**.
The average run time for completing a job that runs C+T-pseudo, Lassosum2-pseudo, LDpred2-pseudo, and ensemble PRS for ~1.2 million HapMap3 SNPs using 2 CPUs (with 30 GB RAM) is approximately 2.5 hours, while increasing to 4 CPUs reduced the run time to approximately two hours.
With ~1.2 million SNPs, single-ancestry analysis pipelines typically require 30GB memory, while for multi-ancestry analysis pipelines, it is recommended that a 25GB * #ancestries is requested to ensure job completion.

Note: fitting the following models with > 1 million SNPs may generate large temporary files (> 20GB per job), and please make sure you have enough storage space to run multiple jobs in parallel before submitting jobs. 
  
  LDpred2-pseudo
  LDpred2-auto
  lassosum2-pseudo
  DBSLMM
  
The temporary/intermediate files in the output folder will be cleaned up if a job is completed successfully. However, when a job unexpectedly fails, the large temporary files should be manually deleted to free up space, especially the subfolder `/PRS_model_training/`.


## Contact
Please report questions and bugs on the Issues page or contact us at pennprs@googlegroups.com.


## Citation
Jin, J., Li, B., Wang, X., Yang, X., Li, Y., Wang, R., Ye, C., Shu, J., Fan, Z., Xue, F. and Ge, T., 2025. PennPRS: a centralized cloud computing platform for efficient polygenic risk score training in precision medicine. medRxiv, 2025-02. [Link](https://www.medrxiv.org/content/10.1101/2025.02.07.25321875v1)


