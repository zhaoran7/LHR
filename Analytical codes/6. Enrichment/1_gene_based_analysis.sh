#!/bin/bash

gwas_height="/data/GWAS/height"
gwas_cmd="/data/GWAS/CMD"
magmadir="/data/magma"
datadir="$magmadir/data"
annodir="$magmadir/anno"
outdir="$magmadir/output"
ref_bfile="$magmadir/ref/g1000_eur/g1000_eur"
gene_loc="$magmadir/ref/NCBI38.gene.loc"

mkdir -p "$datadir" "$annodir" "$outdir"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# step 1: process gwas files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
process_gwas() {
    f=$1
    name=$2
    chr_pos_col=$3
    p_col=$4

    if [ ! -f "$datadir/${name}.snp.p.txt" ]; then
        echo "Processing $name..."
        gunzip -c "$f" | awk -v c="$chr_pos_col" 'NR>1 {print $c, $1, $2}' > "$datadir/${name}.snp.chr.pos.txt"
        gunzip -c "$f" | awk -v p="$p_col" 'NR>1 {val=($p+0 < 1e-300 ? 1e-300 : $p+0); print $1, val}' > "$datadir/${name}.snp.p.txt"
    else
        echo "$name files exist, skip."
    fi
}

# process height and LHR
process_gwas "$gwas_height/height.gz" "height" 2 10
process_gwas "$gwas_height/LHR.gz" "LHR" 2 10

# process CMD
for f in "$gwas_cmd"/*.gz; do
    base=$(basename "$f" .gz)
    process_gwas "$f" "$base" 5 7
done

# annotation
annot="$annodir/genes_annot.genes.annot"
if [ ! -f "$annot" ]; then
    echo "Running MAGMA annotation..."
    magma --annotate window=35,10 --gene-loc "$gene_loc" --snp-loc "$datadir/height.snp.chr.pos.txt" --out "${annot%.genes.annot}"
else
    echo "MAGMA annotation exists, skip."
fi

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# step 2: run gene-base analysis 🏃
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
run_magma() {
    file=$1
    name=$(basename "$file" .snp.p.txt)
    nval=$2
    echo "Running MAGMA gene analysis: $name"
    magma --bfile "$ref_bfile" --pval "$file" N=$nval --gene-annot "$annot" --out "$outdir/$name"
}

[ -f "$datadir/height.snp.p.txt" ] && run_magma "$datadir/height.snp.p.txt" 457892
[ -f "$datadir/LHR.snp.p.txt" ] && run_magma "$datadir/LHR.snp.p.txt" 457381

for f in "$datadir"/*.snp.p.txt; do
    base=$(basename "$f" .snp.p.txt)
    [[ $base == "height" || $base == "LHR" ]] && continue
    run_magma "$f" 520210
done

echo "All done."
