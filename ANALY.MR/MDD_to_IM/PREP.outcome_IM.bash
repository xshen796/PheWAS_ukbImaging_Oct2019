cd /exports/igmm/eddie/GenScotDepression/shen/MendelianRandomisation/MR_IM_20k/run_MR/MDD_to_IM

while read var
list_outcome_gwas='../../GWAS/'$var'.bgenie.QC_forLDSC'

do

if test -z "$var" 
then
      break
else
         echo -n "" > ./outcome_stats/$var.outcome_match
	 cat ./exposure_stats/MDD.exposure_dat | while read line
	 do
		  a=($line) 
		  grep -w ${a[0]} $list_outcome_gwas >> ./outcome_stats/$var.outcome_match
	 done

fi

echo $var outcome matched to exposure dat   done >> XS.log

done < 'ls.IMtraits'

