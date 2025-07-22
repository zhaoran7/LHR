#!/bin/bash
source /work/sph-zhaor/miniconda3/etc/profile.d/conda.sh
conda activate ldsc  # Python2 environment

dir0=/work/sph-zhaor
refdir=$dir0/analysis/ldsc/data/ref/hm3
cleandir=$dir0/analysis/ldsc/data/clean
logdir=$dir0/analysis/ldsc/log
resdir=$dir0/analysis/ldsc/result/height2cmd

#---------------------------------------
# run genetic correlation analysis using LDSC 🏃
#---------------------------------------
dats=($(cd $cleandir; ls *.sumstats.gz | sed 's/\.sumstats\.gz//')) # standardized data

for ((i=0; i<${#dats[@]}-1; i++)); do
    first=${dats[i]}
    other=("${dats[@]:i+1}")

    input_files="$cleandir/$first.sumstats.gz"
    for trait in "${other[@]}"; do
        input_files+=",$cleandir/$trait.sumstats.gz"
    done

    # run LDSC
    python2 $dir0/apps/ldsc/ldsc.py --rg $input_files \
        --out $logdir/$first.rg \
        --ref-ld-chr $refdir/eur_w_ld_chr/ --w-ld-chr $refdir/eur_w_ld_chr/

    # extract all results and save
    beginl=$(awk '$1=="Summary" {printf NR}' $logdir/$first.rg.log)
    awk -v s=$beginl 'FNR > s' $logdir/$first.rg.log | head -n -3 | sed 's/.sumstats.gz//g' \
        > $resdir/$first.rg.txt
done

# all in one
awk 'NR==1 || FNR>1' $resdir/*.rg.txt | sed 's/  */ /g; s/^ //' \
    > $resdir/all.rg.res
