#!/bin/bash

# reference files of s-LDSC can be downloaded at https://zenodo.org/records/10515792
baseline_dir=/work/sph-zhaor/analysis/ldsc/data/ref/1000G/1000G_Phase3_baselineLD_v2.2_ldscores
bfile_dir=/work/sph-zhaor/analysis/ldsc/data/ref/1000G/1000G_Phase3_plink
output_dir=/work/sph-zhaor/analysis/ldsc/data/ref/1000G/1000G_Phase3_functional_ldscores
ldsc_path=/work/sph-zhaor/apps/ldsc/ldsc.py
snplist=/work/sph-zhaor/analysis/ldsc/data/ref/1000G/hm3_no_MHC.list.txt

# obtain functional categroies (start at col5)
func_list=$(zcat $baseline_dir/baselineLD.1.annot.gz | head -n 1 | awk '{for (i=5;i<=NF;i++) print $i}')

#---------------------------------------
# calculate LDscores for every functional category 🧮
#---------------------------------------
for func in $func_list; do
  script=${func}.cmd
  log_out=${func}.out
  log_err=${func}.err
#creat command file for every category
  cat > $script <<EOF
#!/bin/bash
#BSUB -J $func
#BSUB -q short
#BSUB -n 40
#BSUB -o $log_out
#BSUB -e $log_err

source /work/sph-zhaor/miniconda3/etc/profile.d/conda.sh
conda activate ldsc

mkdir -p $output_dir/$func

for chr in {1..22}; do
  zcat $baseline_dir/baselineLD.\${chr}.annot.gz | \
    awk -v col="$func" '
      BEGIN { FS=OFS="\t"; }
      NR==1 {
        for (i=1; i<=NF; i++) {
          if (\$i==col) target_col=i;
        }
      }
      NR>=1 {
        print \$1, \$2, \$3, \$4, \$(target_col);
      }
    ' > $output_dir/$func/${func}.\${chr}.annot

  python2 $ldsc_path \\
    --l2 \\
    --bfile $bfile_dir/1000G.EUR.QC.\${chr} \\
    --ld-wind-cm 1 \\
    --annot $output_dir/$func/${func}.\${chr}.annot \\
    --out $output_dir/$func/${func}.\${chr} \\
    --print-snps $snplist
done
EOF

  bsub < $script # submit jobs to SUSTECH HPC 
done
