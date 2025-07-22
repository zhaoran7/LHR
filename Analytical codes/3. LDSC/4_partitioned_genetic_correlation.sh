#!/bin/bash

# here, we set tissue-specific partitioned genetic association analysis as example
# for analysis of gene functional categories, you can just change the ldscore files
# make ldscore files of 53 tissues on GTEx is not needed, you can just download here (require payment):
# https://console.cloud.google.com/storage/browser/_details/broad-alkesgroup-public-requester-pays/LDSCORE/LDSC_SEG_ldscores/biorxiv/GTEx_1000Gv3.tgz

dir0=/work/sph-zhaor
ldsc_path=$dir0/apps/ldsc/ldsc.py
cleandir=$dir0/analysis/ldsc/data/clean
refdir=$dir0/analysis/ldsc/data/ref/GTEx/GTEx_ldscores
frqdir=$dir0/analysis/ldsc/data/ref/1000G/1000G.EUR.QC
wlddir=$dir0/analysis/ldsc/data/ref/1000G/1000G_Phase3_weights_hm3_no_MHC
resdir=$dir0/analysis/ldsc/result/height2cmd/GTEx
mkdir -p $resdir

# define exposure and outcome
exposures=(height LHR)
outcomes=(AF VTE CM AA HF CKD HS IS PAD AS T2D TIA CAD HT MI)

cd $refdir
func_list=($(ls GTEx.*.*.l2.ldscore.gz | cut -d. -f1,2 | sort -u)) #list all ldscore files of 53 tissues

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# step1: run partitioned genetic correlation analysis 🏃
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for func in "${func_list[@]}"; do
  jobfile=$resdir/run_${func}.cmd
  echo "#!/bin/bash" > $jobfile
  echo "source /work/sph-zhaor/miniconda3/etc/profile.d/conda.sh" >> $jobfile
  echo "conda activate ldsc" >> $jobfile

  refpath=$refdir/${func}.
  outdir=$resdir/$func
  mkdir -p $outdir

  for outcome in "${outcomes[@]}"; do
    for exposure in "${exposures[@]}"; do
      log=$outdir/${outcome}_vs_${exposure}.log
      echo "python2 $ldsc_path \
        --ref-ld-chr $refpath \
        --out $outdir/tmp_${outcome}_${exposure} \
        --overlap-annot \
        --rg $cleandir/${outcome}.sumstats.gz,$cleandir/${exposure}.sumstats.gz \
        --frqfile-chr $frqdir/1000G.EUR.QC. \
        --w-ld-chr $wlddir/weights.hm3_noMHC. \
        > $log 2>&1" >> $jobfile

      echo "beginl=\$(awk '\$1==\"Summary\" {print NR}' $log)" >> $jobfile
      echo "awk -v s=\$beginl 'FNR > s' $log | head -n -3 | awk -v f=$func -v o=$outcome -v e=$exposure '{print f\" \"o\" \"e\" \"\$3\" \"\$4\" \"\$5\" \"\$6\" \"\$7\" \"\$8}' >> $outdir/${func}.rg.txt" >> $jobfile
      echo "rm -f $outdir/tmp_${outcome}_${exposure}.*" >> $jobfile
    done
  done

  bsub -q short -n 40 -J rg.$func -o $resdir/log_$func.out -e $resdir/log_$func.err < $jobfile # submit jobs to SUSTECH HPC
done

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# step2: all in one 1️⃣
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#!/bin/bash

resdir=/work/sph-zhaor/analysis/ldsc/result/height2cmd/GTEx
outfile=$resdir/all_partitioned_rg.res
echo -e "Annot Trait1 Trait2 rg se z p h2_1 h2_2" > $outfile

find $resdir -type f -name "*.rg.txt" | while read file; do
  awk 'NR==1 && FNR==1 { next } 1' "$file" >> $outfile
done

