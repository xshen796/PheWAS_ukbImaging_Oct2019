

while read var
list_exposure='./exposure_stats/'$var'.exposure_dat'

do

if test -z "$var" 
then
      break
else
         echo -n "" > ./outcome_stats/$var.MDDoutcome_match
	 cat $list_exposure | while read line
	 do
		  a=($line) 
		  grep -w ${a[0]} /exports/igmm/eddie/GenScotDepression/shen/MDD_meta/mdd_meta_Shen/23andmePGCUKBnoMRI_20181028_shen_3cohorts.meta >> ./outcome_stats/$var.MDDoutcome_match
	 done

fi

echo $var outcome matched to exposure dat   done >> XS.log

done < 'ls.IMtraits'

