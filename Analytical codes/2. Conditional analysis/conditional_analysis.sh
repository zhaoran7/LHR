#!/bin/bash

raw_dir="/work/sph-zhaor/data/gwas/main"
ma_dir="/work/sph-zhaor/analysis/gcta.cojo/data"
bfile_dir="/work/sph-zhaor/data/ukb/gen/typ"
snplist="/work/sph-zhaor/analysis/gcta.cojo/output/lhr/lhr.snplist"

gwas_list=("height.gz" "LHR.gz")

#---------------------------------------
# covert gwas data to GCTA-COJO format (standardization)
#---------------------------------------
for gwas_file in "${gwas_list[@]}"; do
  base="${gwas_file%.gz}"
  in_file="$raw_dir/$gwas_file"
  ma_file="$ma_dir/${base}.ma"

  if [ ! -f "$ma_file" ]; then
    header=$(zcat "$in_file" | head -n1)
    for col in SNP EA NEA EAF BETA SE P N; do
      eval col_${col}=$(echo "$header" | tr '\t' '\n' | grep -n -w "$col" | cut -d: -f1)
    done

    zcat "$in_file" | awk -v OFS="\t" \
      -v snp=$col_SNP -v ea=$col_EA -v nea=$col_NEA -v eaf=$col_EAF \
      -v beta=$col_BETA -v se=$col_SE -v p=$col_P -v n=$col_N \
      'NR==1 {print "SNP","EA","NEA","EAF","BETA","SE","P","N"}
       NR>1 && $snp && $ea && $nea && $eaf+0==$eaf && $beta+0==$beta && $se+0==$se && $p+0==$p && $n+0==$n {
         print $snp, $ea, $nea, $eaf, $beta, $se, $p, $n
       }' > "$ma_file"
  fi

#---------------------------------------
# run conditional analysis 🏃
#---------------------------------------
  out_dir="/work/sph-zhaor/analysis/gcta.cojo/output/${base}"
  cmd_dir="/work/sph-zhaor/analysis/gcta.cojo/cmd/${base}"

  mkdir -p "$out_dir" "$cmd_dir"

  for chr in {1..22}; do
    job="${base}.joint2.chr${chr}"
    bfile="${bfile_dir}/chr${chr}"
    out="${out_dir}/${job}"
    cmd="${cmd_dir}/${job}.cmd"
    log="${cmd_dir}/${job}.log"
    err="${cmd_dir}/${job}.err"

    cat <<EOF > "$cmd"
#!/bin/bash
gcta64 --bfile $bfile --chr $chr --maf 0.01 --cojo-file $ma_file --cojo-slct --out $out
EOF

    chmod +x "$cmd"
    bsub -q short -n 40 -J "$job" -o "$log" -e "$err" < "$cmd"
  done
done



