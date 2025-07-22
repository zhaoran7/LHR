#!/bin/bash
source /work/sph-zhaor/miniconda3/etc/profile.d/conda.sh
conda activate ldsc # python2 environment

gwasdir=/work/sph-zhaor/analysis/ldsc/data/raw
refdir=/work/sph-zhaor/analysis/ldsc/data/ref/hm3

#---------------------------------------
# covert gwas data to ldsc format (standardization)
#---------------------------------------
dats=`cd $gwasdir; ls *.gz | sed 's/\.gz//g'`
for dat in $dats; do
	if [[ -f $dat.sumstats.gz ]]; then echo $dat already run; continue; fi
	echo $dat
	SNP="SNP";EA="EA";NEA="NEA";BETA="BETA";P="P";N="N"           
	datf=$gwasdir/$dat.gz

	if [[ -z "$N" ]]; then
		n_str="--N 100000"
	else
		n_str="--N-col $N"
	fi

	bad_p=`zcat $datf | awk -v p=$P 'NR>1 && $p<=0' | wc -l`
	if [[ $bad_p -gt 0 ]]; then
		zcat $datf | awk -v p=$P '{if($p<=1e-300) $p=1e-300; print $0}' | sed 's/ /\t/g' | gzip -f > $dat.gz
		datf=$dat.gz
	fi

#---------------------------------------
# run munge_sumstats.py 🏃
#---------------------------------------
	python $dir0/apps/ldsc/munge_sumstats.py \
		--chunksize 10000 \
		--sumstats $datf \
		--merge-alleles $refdir/w_hm3.snplist \
		--out $dir0/analysis/ldsc/data/clean/$dat \
		--snp $SNP \
		--a1 $EA \
		--a2 $NEA \
		$n_str \
		--signed-sumstats $BETA,0 \
		--p $P \
		--ignore SNPID,OR
done