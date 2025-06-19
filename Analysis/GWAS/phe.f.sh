header_names() {

	Arr1=("SNP" "CHR" "POS" "EA" "NEA" "EAF" "N" "BETA" "SE" "P")
	Arr2=("^snp$|^id$|^rsid|^variant_idX$" "^chr$|^chrom$|^chromosome$" "^pos$|^bp$|^base_pair$" "^a1$|^ea$|^eff.allele$|^effect_allele$|^allele1$" "^OMITTED$|^nea$|^ref.allele$|^other_allele$|^allele0$" "^eaf$|^a1freq$|^a1_freq$|^effect_allele_frequency" "^n$|^obs_ct$|^Neff$" "^beta$" "^se$|^standard_error" "^p$|^pval$|^p_value$|^p_bolt_lmm$")

	datf=$1
	if [[ $datf == *.gz ]]; then
		head_row=$(zcat $datf | head -1 | sed 's/\t/ /g')
	else
		head_row=$(cat $datf | head -1 | sed 's/\t/ /g')
	fi
	head_row=$(echo $head_row | dos2unix) # 🏮

	for i in ${!Arr1[@]}; do
		eval ${Arr1[$i]}=`echo $head_row | tr ' ' '\n' | grep -Einw ${Arr2[$i]} | sed 's/\w*://'`
		eval ${Arr1[$i]}_col=`echo $head_row | tr ' ' '\n' | grep -Einw ${Arr2[$i]} | sed 's/:\w*//'`
	done
	echo dat $datf, snp $SNP $SNP_col, chr $CHR $CHR_col, pos $POS $POS_col, ea $EA $EA_col, nea $NEA $NEA_col, eaf $EAF $EAF_col, n $N $N_col, beta $BETA $BETA_col, se $SE $SE_col, p $P $P_col

}
