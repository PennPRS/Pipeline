# PennPRS Test Runner

A SLURM-driven harness that submits test jobs for the worked examples we provided in
([Test Examples 5.1 – 5.7](https://github.com/PennPRS/Pipeline/wiki/5.-Test-Examples))
against a local clone of [PennPRS/Pipeline](https://github.com/PennPRS/Pipeline).

## Layout

```
test-runner/
├── run_all.sh           # orchestrator
├── tests.yaml           # test matrix — edit here to tweak tests
├── lib/
│   ├── common.sh        # logging + yq detection + ${VAR} expansion
│   ├── parse.sh         # YAML → bash record files
│   ├── submit.sh        # sbatch wrapper (supports --wait and parallel modes)
│   └── verify.sh        # existence + schema checks
├── logs/                # Slurm stdout/stderr per job (auto-created)
└── results/             # per-test records, jobids, verify status
```
## Preparation

Before running the test, please follow the instructions on [pipeline installation](https://github.com/PennPRS/Pipeline/wiki/1.-Installation) and ensure the following steps are completed:

1. LD and supported file download: make sure the LD information for EUR and EAS were successfully downloaded and saved in `PennPRS/LD/`.
2. Environment setup

If needed, please modify the command `module load r` and `module load anaconda` in all the job submission scripts (.sh files) in `test/job_submission/` to match the module configuration on your server.
To see what R and conda versions are available to load in your system, please run 
```
module avail
```


## Requirements

- A SLURM cluster (the runner submits every example with `sbatch`). If you use other types of clusters (e.g., bsub, qsub), please adjust job submission command accordingly (see [Setup](https://github.com/PennPRS/Pipeline/tree/main/test-runner#setup))
- `yq` — either [mikefarah/yq](https://github.com/mikefarah/yq) (Go) or
  [kislyuk/yq](https://github.com/kislyuk/yq) (Python). Install with
  any of:
  ```bash
     cd test-runner/
     pip install yq --user                    # python-based
     conda install -c conda-forge yq          # if pip install does not work
     # or download a binary from the mikefarah/yq releases page
  ```
- The pipeline itself, installed following
  [wiki § 1.-Installation](https://github.com/PennPRS/Pipeline/wiki/1.-Installation).


## Setup

1. Edit `tests.yaml`:
   - Change `PennPRS_path` in line 19 to the path to your local cloned PennPRS/ folder.
   - If your cluster uses a different partition / qos than `cpu` /
     `normal`, adjust `globals.slurm` (lines 32 - 36).
   - (Optional) Tweak per-test `slurm.time`, `cpus`, or `mem_per_cpu` if your
     cluster is faster or tighter on memory than the wiki's reference
     (mid-range Intel Xeon).

2. Please modify the command `module load r` and `module load anaconda`) inside each
  `.sh` script under `$test/job_submission/` to match the module configuration on your server.
  
3. The job submission command may vary depending on your computing environment or scheduler. Please modify the  
command in `sbatch_args` (line 39 in `test-runner/lib/submit.sh`) as needed for your server. For example, replace `myaccount` in `test-runner/lib/submit.sh`, with your own account name.

    

**Usage**

```bash
cd test-runner/

# Run all examples sequentially.
./run_all.sh

# Run all examples, submit up to 4 jobs at once (uses sbatch without --wait; polls squeue).
./run_all.sh --parallel 4

# Just list what is configured.
./run_all.sh --list

# Run a single test.
./run_all.sh --only 5.1a

# Run every test under section 5.5 and section 5.7.
./run_all.sh --only 5.5,5.7

# Print the sbatch command the runner would issue, without submitting.
./run_all.sh --dry-run --only 5.2a

# Skip submission; re-run existence/schema checks on jobs that already ran.
./run_all.sh --verify-only
```


**Exit codes**

- `0` — every selected test was submitted **and** passed verification
- `1` — any failure.
- `2` — invalid invocation.


**Test matrix**

| ID    | Section | Description                                         | CPUs | Mem / cpu | Wall-time |
| :---- | :------ | :-------------------------------------------------- | :--: | :-------: | :-------: |
| 5.1a  | 5.1     | SA Pseudo-Training — EUR continuous                 | 11   | 1 G       | 2 h       |
| 5.1b  | 5.1     | SA Pseudo-Training — EUR binary                     | 11   | 1 G       | 2 h       |
| 5.2a  | 5.2     | SA PRS-CS (grid) — EUR continuous                   | 11   | 1.5 G     | 4 h       |
| 5.2b  | 5.2     | SA PRS-CS (grid) — EUR binary                       | 11   | 1.5 G     | 4 h       |
| 5.2c  | 5.2     | SA PRS-CS (auto) — EUR continuous                   | 11   | 1.5 G     | 4 h       |
| 5.2d  | 5.2     | SA PRS-CS (auto) — EUR binary                       | 11   | 1.5 G     | 4 h       |
| 5.3a  | 5.3     | SA Tuning-Parameter-Free — EUR continuous           | 11   | 1 G       | 2 h       |
| 5.3b  | 5.3     | SA Tuning-Parameter-Free — EUR binary               | 11   | 1 G       | 2 h       |
| 5.4a  | 5.4     | PRS-CSx-pseudo — EUR+EAS continuous                 | 11   | 2 G       | 12 h      |
| 5.4b  | 5.4     | PRS-CSx-pseudo — EUR+EAS binary                     | 11   | 2 G       | 12 h      |
| 5.5a  | 5.5     | MUSSEL-pseudo — EUR+EAS continuous                  | 11   | 2 G       | 10 h      |
| 5.5b  | 5.5     | MUSSEL-pseudo — EUR+EAS binary                      | 11   | 2 G       | 10 h      |
| 5.6a  | 5.6     | PROSPER-pseudo — EUR+EAS continuous                 | 11   | 2 G       | 5 h       |
| 5.6b  | 5.6     | PROSPER-pseudo — EUR+EAS binary                     | 11   | 2 G       | 5 h       |
| 5.7   | 5.7     | Model Evaluation with individual-level data         | 1    | 5 G       | 30 m      |



**Outputs**

Each test has:

```
results/<id>.rec          # parsed record (bash-sourceable)
results/<id>.submit.log   # sbatch stdout
results/<id>.jobid        # Slurm job id
results/<id>.state        # final Slurm state (parallel mode only)
results/<id>.outdir       # resolved output folder (set after verify)
results/<id>.verify       # PASS or FAIL
logs/<id>_<jobid>.out     # Slurm job stdout
logs/<id>_<jobid>.err     # Slurm job stderr
```

A terminal summary is printed at the end of `run_all.sh`.



## Test examples 

To test examples in [wiki](https://github.com/PennPRS/Pipeline/wiki/5.-Test-Examples) across all supported methods and modes, recommend running single tests on each category:

```
bash
cd test-runner/

# Run all examples in parallel (recommended), which submits up to 4 jobs at once, takes approximately 2 - 3 hours to complete.
./run_all.sh --parallel 4

# Or run all examples sequentially, takes approximately 8 hours to complete.
./run_all.sh
```

### Test Example Outputs

Reference outputs for all test cases are provided on [Dropbox](https://www.dropbox.com/scl/fo/84rvif06ffausgkczorpw/AELVtRppWb5ewlVp1Hf7U0o?rlkey=beavrxi18raym5w2j920j74ve&st=a3dosg5n&dl=0).

The output folder for every test case in 5.1 - 5.6 will be saved in `$PennPRS_path/test/PennPRSoutput/`; the output for test 5.7 will be saved in `$PennPRS_path/test/Evaluation/output/`. After running the tests, you can compare your local results against these to confirm
correctness.

<table>
  <thead>
    <tr>
      <th>Test</th>
      <th>Output folder (under <code>$PennPRS_path/test/PennPRSoutput/</code>)</th>
      <th>Example output</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>5.1a</td>
      <td nowrap><code>testcontinuous_EUR_C+T.lassosum2.LDpred2_SApseudo_continuous/</code></td>
      <td><a href="https://www.dropbox.com/scl/fo/5rjdc1jy5bycj248r0u8c/ADevbdWjBLn2_h7njvD3HVo?rlkey=ny13v8exjw9i4hnr75iekkeq9&st=qf54lbeb&dl=0">Link</a></td>
    </tr>
    <tr>
      <td>5.1b</td>
      <td nowrap><code>testbinary_EUR_C+T.lassosum2.LDpred2_SApseudo_binary/</code></td>
      <td><a href="https://www.dropbox.com/scl/fo/jhnhe0p4ymz8v81zrpb5c/AIsBVeDgVTb6cG-l5OLBtuc?rlkey=sv0fi31zn6z6rmziqa08n57um&st=3qg1m9rd&dl=0">Link</a></td>
    </tr>
    <tr>
      <td>5.2a</td>
      <td nowrap><code>testcontinuous_EUR_PRS-CS_PRSCSpseudo_continuous/</code></td>
      <td><a href="https://www.dropbox.com/scl/fo/nu9s2vnvnix9ci7cpty6o/AC-LGV2SntbtW3KQ5Ztm-9U?rlkey=ookz2ga2uw6tfz23badp0bqtu&st=n2ukajdg&dl=0">Link</a></td>
    </tr>
    <tr>
      <td>5.2b</td>
      <td nowrap><code>testbinary_EUR_PRS-CS_PRSCSpseudo_binary/</code></td>
      <td><a href="https://www.dropbox.com/scl/fo/wym7ksuqivjqrfugll97n/AJfYTccn61A_0MgBEcgXo_c?rlkey=9h92sgt2gfswgozd6k4auwegd&st=8c0qq4tx&dl=0">Link</a></td>
    </tr>
    <tr>
      <td>5.2c</td>
      <td nowrap><code>testcontinuous_EUR_PRS-CS_PRSCSauto_fullyBayes_continuous/</code></td>
      <td><a href="https://www.dropbox.com/scl/fo/gdlvtcpgeqmtjzra7cnwe/AKYou2FZIsYKiG-zLUlY9p4?rlkey=e2pjr8q6w0ercnoyg2p8ub9w2&st=fieu7r6b&dl=0">Link</a></td>
    </tr>
    <tr>
      <td>5.2d</td>
      <td nowrap><code>testbinary_EUR_PRS-CS_PRSCSauto_fullyBayes_binary/</code></td>
      <td><a href="https://www.dropbox.com/scl/fo/qf3pcwlndcb8s7nyyzxqn/AFUopXhW6xCUSNk7XCTfoDw?rlkey=nrd1wohos79mww4kkj8fioycl&st=6p45a34v&dl=0">Link</a></td>
    </tr>
    <tr>
      <td>5.3a</td>
      <td nowrap><code>testcontinuous_EUR_LDpred2-auto.DBSLMM_tuningfree_continuous/</code></td>
      <td><a href="https://www.dropbox.com/scl/fo/o2xxkfqbziufdj3lxbztl/AB9Y7OZvQAz9ijhRPge0Y6k?rlkey=8tp04mqs63dwvzoz23khns4nn&st=7688vq13&dl=0">Link</a></td>
    </tr>
    <tr>
      <td>5.3b</td>
      <td nowrap><code>testbinary_EUR_LDpred2-auto.DBSLMM_tuningfree_binary/</code></td>
      <td><a href="https://www.dropbox.com/scl/fo/ci7oqkjgro4rzv2g62rk9/AOK3k3Ctk3uu0O_nVP5Mn4I?rlkey=gygxxwczwpkwtw17iljvgnx0v&st=fylft36g&dl=0">Link</a></td>
    </tr>
    <tr>
      <td>5.4a</td>
      <td nowrap><code>testcontinuous_EUR.EAS_PRS-CSx_PRSCSXpseudo_continuous/</code></td>
      <td><a href="https://www.dropbox.com/scl/fo/y6vkg7zu11s81lx1mbtr0/AH84TzjdQbdln5z3ziWGRYE?rlkey=9jlyz58r8uy1lt0x2ciru4oep&st=de61z3ph&dl=0">Link</a></td>
    </tr>
    <tr>
      <td>5.4b</td>
      <td nowrap><code>testbinary_EUR.EAS_PRS-CSx_PRSCSxpseudo_binary/</code></td>
      <td><a href="https://www.dropbox.com/scl/fo/nkb85wo5rnfq06906xpjl/AJTEPF72001hzlVIh7hBPys?rlkey=8dgodr0hscc9qym4zsr51s6be&st=ct3vmdqm&dl=0">Link</a></td>
    </tr>
    <tr>
      <td>5.5a</td>
      <td nowrap><code>testcontinuous_EUR.EAS_MUSSEL_MUSSELpseudo_continuous/</code></td>
      <td><a href="https://www.dropbox.com/scl/fo/7aqm8dig0zpptggm37itq/ABZGD3Qm4-9cLgdigGcIjr8?rlkey=x8tm2023tdy5s6coatrltcvnl&st=snrvhgce&dl=0">Link</a></td>
    </tr>
    <tr>
      <td>5.5b</td>
      <td nowrap><code>testbinary_EUR.EAS_MUSSEL_MUSSELpseudo_binary/</code></td>
      <td><a href="https://www.dropbox.com/scl/fo/300ytqgaiv7jen59nius9/AL7ZpsnWWerdzgtiWPCCnx8?rlkey=jj31okcf3najx0gb2j041pdye&st=0h8pqy3x&dl=0">Link</a></td>
    </tr>
    <tr>
      <td>5.6a</td>
      <td nowrap><code>testcontinuous_EUR.EAS_PROSPER_PROSPERpseudo_continuous/</code></td>
      <td><a href="https://www.dropbox.com/scl/fo/xwergagb2y77vk1hz6h7z/AKEmhVCyZiUYlYy_e4HJsyA?rlkey=4c6j6vw8qd89gc8dgywmbt4uq&st=w7heuycm&dl=0">Link</a></td>
    </tr>
    <tr>
      <td>5.6b</td>
      <td nowrap><code>testbinary_EUR.EAS_PROSPER_PROSPERpseudo_binary/</code></td>
      <td><a href="https://www.dropbox.com/scl/fo/d63h5pekb9oa2twmorzlb/ANV-XuGc57dZH0BHWM9b6hU?rlkey=64wm4mgelvtm2qowohujogbzs&st=g69shyu8&dl=0">Link</a></td>
    </tr>
    <tr>
      <td>5.7</td>
      <td nowrap><code>$PennPRS/test/Evaluation/output/</code> (contains <code>evaluation_results.txt</code>)</td>
      <td><a href="https://www.dropbox.com/scl/fo/5q29yxbt8jt4dj5u08d0d/AH3SPeQcUOggTshFAeIihzY?rlkey=n4tu29lh6w1ki026qeqk0athm&st=1378juy6&dl=0">Link</a></td>
    </tr>
  </tbody>
</table>

## Runtime and memory

- Running `./run_all.sh --parallel 4` (run 4 jobs parallel) should take less than 3 hours to complete.
- With less than 2000 SNPs per example GWAS data file, the average run time for completing a job using 11 cores is less than 30 minutes for each single-ancestry job and less than one hour for a multi-ancestry job.

