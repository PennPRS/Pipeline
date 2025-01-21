#!/bin/bash

while getopts "D:p:B:s:m:H:n:G:P:l:T:c:i:t:o:" opt; do
  case $opt in
    D) software_path="$OPTARG"
    ;;
    p) plink="$OPTARG"
    ;;
    B) block_prefix="$OPTARG"
    ;;
    s) summary_file_prefix="$OPTARG"
    ;;
    m) model="$OPTARG"
    ;;
    T) type="$OPTARG"
    ;;
    H) herit="$OPTARG"
    ;;
    G) val_geno_prefix="$OPTARG"
    ;;
    o) outpath="$OPTARG"
    ;;
    P) val_pheno="$OPTARG"
    ;;
    l) col="$OPTARG"
    ;;
    c) cov="$OPTARG"
    ;;
    i) index="$OPTARG"
    ;;
    t) thread="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

printf "\033[33mArgument software_path is %s  \033[0m\n" "$software_path"
printf "\033[33mArgument plink is %s  \033[0m\n" "$plink"
printf "\033[33mArgument block_prefix is %s  \033[0m\n" "$block_prefix"
printf "\033[33mArgument summary_file_prefix is %s  \033[0m\n" "$summary_file_prefix"
printf "\033[33mArgument model is %s  \033[0m\n" "$model"
printf "\033[33mArgument herit is %s  \033[0m\n" "$herit"
printf "\033[33mArgument val_geno_prefix is %s  \033[0m\n" "$val_geno_prefix"
printf "\033[33mArgument type is %s  \033[0m\n" "$type"
if [[ "${type}" == "t" ]]
then
	printf "\033[33mArgument valid_pheno is %s  \033[0m\n" "$val_pheno"
	printf "\033[33mArgument col is %s  \033[0m\n" "$col"
	printf "\033[33mArgument index is %s  \033[0m\n" "$index"
fi
if [ -n "$cov" ]; then 
	printf "\033[33mArgument cov is %s  \033[0m\n" "$cov"
fi
printf "\033[33mArgument thread is %s  \033[0m\n" "$thread"
printf "\033[33mArgument outpath is %s  \033[0m\n" "$outpath"

DBSLMM=${software_path}DBSLMM/software/DBSLMM.R
TUNE=${software_path}DBSLMM/software/TUNE.R
dbslmm=${software_path}/DBSLMM/scr/dbslmm

# LDSC: heritability and number of SNP
nsnp=`sed -n '24p' ${herit} | cut -d ',' -f 2 | cut -d ' ' -f 2`
h2=`sed -n '26p' ${herit} | cut -d ":" -f 2 | cut -d '(' -f 1 | cut -d " " -f 2`

# DBSLMM: tuning version
if [[ "$type" == "t" ]]
then
	for chr in `seq 1 22`
	do

		BLOCK=${block_prefix}${chr}
		summchr=${summary_file_prefix}${chr}
		nobs=`sed -n "2p" ${summchr}.assoc.txt | awk '{print $5}'`
		nmis=`sed -n "2p" ${summchr}.assoc.txt | awk '{print $4}'`
		n=$(echo "${nobs}+${nmis}" | bc -l)
		val_geno=${val_geno_prefix}${chr}
		Rscript ${DBSLMM} --summary ${summchr}.assoc.txt --outPath ${outpath} --plink ${plink} --model ${model}\
						  --dbslmm ${dbslmm} --ref ${val_geno} --n ${n} --type ${type} --nsnp ${nsnp} --block ${BLOCK}.bed\
						  --h2 ${h2} --h2f 0.8,1,1.2 --thread ${thread}
		for h2f in 0.8 1 1.2
		do
			if [ -f "${outpath}${summchr_prefix}_h2f${h2f}.dbslmm.badsnps" ];then
				rm ${outpath}${summchr_prefix}_h2f${h2f}.dbslmm.badsnps
			fi
			summchr_prefix=`echo ${summchr##*/}`
			${plink}  --silent --bfile ${val_geno} --score ${outpath}${summchr_prefix}_h2f${h2f}.dbslmm.txt 1 2 4 sum\
					  --out ${outpath}${summchr_prefix}_h2f${h2f}
			rm ${outpath}${summchr_prefix}_h2f${h2f}.log
			if [ -f "${outpath}${summchr_prefix}_h2f${h2f}.nopred" ];then
				rm ${outpath}${summchr_prefix}_h2f${h2f}.nopred
			fi
		done
	done

	summchr_prefix2=`echo ${summchr_prefix%_*}`
	if [[ ! -n "$cov" ]]
	then 
	Rscript ${TUNE} --phenoPred ${outpath}${summchr_prefix2} --phenoVal ${val_pheno},${col} \
		   --h2Range 0.8,1,1.2 --index ${index}
	else 
	Rscript ${TUNE} --phenoPred ${outpath}${summchr_prefix2} --phenoVal ${val_pheno},${col} \
		   --h2Range 0.8,1,1.2 --index ${index} --cov ${cov}
	fi

	hbest=`cat ${outpath}${summchr_prefix2}_hbest.${index}`
	for chr in `seq 1 22`
	do
		mv ${outpath}${summchr_prefix2}_chr${chr}_h2f${hbest}.dbslmm.txt ${outpath}${summchr_prefix2}_chr${chr}_best.dbslmm.txt
		rm ${outpath}${summchr_prefix2}_chr${chr}_h2f*
	done

fi

## DBSLMM default version
if [[ "$type" == "d" ]]
then
for chr in `seq 1 22` 
do
	BLOCK=${block_prefix}${chr}
	summchr=${summary_file_prefix}${chr}
	val_geno=${val_geno_prefix}${chr}
	summchr=${summary_file_prefix}${chr}
	nobs=`sed -n "2p" ${summchr}.assoc.txt | awk '{print $5}'`
	nmis=`sed -n "2p" ${summchr}.assoc.txt | awk '{print $4}'`
	n=$(echo "${nobs}+${nmis}" | bc -l)
	echo ${model}
	Rscript ${DBSLMM} --summary ${summchr}.assoc.txt --outPath ${outpath} --plink ${plink} --model ${model}\
					  --dbslmm ${dbslmm} --ref ${val_geno} --n ${n} --nsnp ${nsnp} --block ${BLOCK}.bed\
					  --h2 ${h2} --thread ${thread}
	summchr_prefix=`echo ${summchr##*/}`
	rm ${outpath}${summchr_prefix}.dbslmm.badsnps

done
fi
