import pandas as pd
for pheno_idx in range(1,676):
  print("Start to reformat for pheno_"+str(pheno_idx))
  path = "/your/fastgwa/path/pheno_"+str(pheno_idx)+".fastGWA"
  data = pd.read_csv(path, sep="\t", header=None, names=["chr", "rs", "ps", "allele1", "allele0", "n_obs", "af", "beta", "se", "p_wald"])
  data = data.drop(data.index[0])
  data["chr"] = data["chr"].astype(str)
  data["n_obs"] = data["n_obs"].astype(int)
  total = max(data["n_obs"])
  data["n_mis"] = total - data["n_obs"]
  data = data[["chr", "rs", "ps", "n_mis", "n_obs","allele1", "allele0", "af", "beta", "se", "p_wald"]]
  print("Save pheno_"+str(pheno_idx))
  for i in range(1,23):
    temp = data[data["chr"] == str(i)]
    temp.to_csv("/your/dbslmm/input/path/pheno_"+str(pheno_idx)+"_chr_"+str(i)+".assoc.txt", sep="\t", index=False)

