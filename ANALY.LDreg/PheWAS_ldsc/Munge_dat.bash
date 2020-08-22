##### set up the environment for LDSC : https://github.com/bulik/ldsc

# conda env create --file environment.yml
# conda activate ldsc
# sort -rn -k13,13 ampN10.bgenie.QC | head 

if [ ! -d './LDready.stats' ]; then
	 mkdir LDready.stats
fi

for fname in ./bgenie.QC.summstats/*.bgenie.QC
do

### prep the stats
## >0.9maximum N, >0.005MAF, info>0.1, phwe>1e-6  --->  for all imaging traits
awk '{ if($6 >=  0.005) { print }}' $fname | awk '{ if($6 <=  0.995) { print }}' | awk '{ if($7 >= 0.1) { print }}' | awk '{ if($8 >= 1e-6) { print }}'| awk '{ if($13>13126){ print }}' >  tmpdat.txt
awk '{print $1,$2,$3,$5,$4,$6,$7,$8,$9,$10,$11,$12,$13}' tmpdat.txt > tmpdat.middle
echo -e "chr SNP pos A1 A2 af info hwe BETA se T p N" | cat - tmpdat.middle > $fname'_forLDSC'

### munge data
munge_sumstats.py \
--sumstats $fname'_forLDSC' \
--N 14585 \
--out $fname'_forLDSCmunged' \
--merge-alleles ../ldsc/w_hm3.snplist

echo $fname done
 
done

mv ./bgenie.QC.summstats/*_forLDSC* ./LDready.stats

if [ ! -d './LDready.stats/logs' ]; then
	 mkdir LDready.stats/logs
fi

mv *.log ./LDready.stats/logs
rm tmpdat.*
