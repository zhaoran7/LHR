#!/bin/bash

threads=40
impdir=/work/sph-zhaor/data/ukb/gen/imp
label=height; binary=N
pheno_f=/work/sph-zhaor/data/ukb/phe/ukb.pheno; phenos="height leg LHR"
covar_f=$pheno_f; covars="age,sex,center,PC1-PC40"

if [[ $binary == Y ]]; then col_str="+a1freqcc,+a1countcc,+totallelecc"; phe_str="--1"; postfix="logistic.hybrid"
else col_str="+a1freq,+a1count,+totallelecc"; phe_str="--1 --pheno-quantile-normalize"; postfix="linear"; fi
		
#---------------------------------------
# step1: run PLINK 🏃‍
#---------------------------------------
outdir=/work/sph-zhaor/data/gwas/self/$label; mkdir -p $outdir
for chr in {1..22} X; do
	x_str=""; if [[ $chr == X ]]; then x_str="--xchr-model 2"; fi
	
	echo "#!/bin/bash
	source phe.f.sh
	if [[ ! -f $impdir/chr$chr.extract ]]; then
		plink2 --memory 12000 --threads $threads --pfile $impdir/chr$chr --freq counts --out $impdir/chr$chr
		awk '\$5>100 {print \$2}' $impdir/chr$chr.acount | sed '1 s/ID/SNP/' > $impdir/chr$chr.extract
	fi
	plink2 --memory 12000 --threads $threads --pfile $impdir/chr$chr --extract $impdir/chr$chr.extract --glm cols=+omitted,+nobs,+beta,$col_str hide-covar allow-no-covars no-x-sex --pheno $pheno_f --no-psam-pheno --pheno-name $phenos $phe_str $x_str --covar $covar_f --covar-name $covars --out chr$chr --no-input-missing-phenotype
	for t in $phenos; do 
		sed -i '1 s/^#//; 1 s/?//g' chr$chr.\$t.glm.$postfix
		header_names chr$chr.\$t.glm.$postfix # 🏮
		cat chr$chr.\$t.glm.$postfix | awk -v snp_col=\$SNP_col -v chr_col=\$CHR_col -v pos_col=\$POS_col -v ea_col=\$EA_col -v nea_col=\$NEA_col -v eaf_col=\$EAF_col -v n_col=\$N_col -v beta_col=\$BETA_col -v se_col=\$SE_col -v p_col=\$P_col '{if (NR==1) print \"SNP CHR POS EA NEA EAF N BETA SE P\"; else print \$snp_col,\$chr_col,\$pos_col, \$ea_col,\$nea_col,\$eaf_col,\$n_col, \$beta_col,\$se_col,\$p_col}' | grep -vwE \"nan|inf|NA\" | sed 's/ /\t/g' | sort -k 2,2n -k 3,3n | gzip -f > chr$chr.\$t.glm.$postfix.gz
	#	rm chr$chr.\$t.glm.$postfix
	done
	" > $outdir/chr$chr.cmd 
	
	cd $outdir # -q smp -n 10 -R "span[ptile=192]"
	bsub -q short -n 40 -J $label.chr$chr -o chr$chr.LOG -e chr$chr.ERR < chr$chr.cmd #submit jobs to SUSTECH HPC
done

#---------------------------------------
# step2: all in one 1️⃣
#---------------------------------------
cat <<EOF > $outdir/step2.cmd
#!/bin/bash
for t in $phenos; do
    ls -v chr*.\$t.glm.$postfix.gz | xargs zcat | awk 'NR==1 || \$8 - /[0-9]/' | gzip -f > ../\$t.gz
done
EOF

cd $outdir
bsub -q short -n 40 -J step2 -o step2.LOG -e step2.ERR < step2.cmd # submit jobs to SUSTECH HPC