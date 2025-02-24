#qlogin -l h_vmem=64G

 #module add igmm/apps/R/3.2.4
 #module add R


### prepare log file
if [ -f 'XS.log' ]; then
    rm XS.log
fi
    echo $(date)  > XS.log
    echo XShen script >> XS.log
    echo  >> XS.log
 


### prepare exposure
	if [ ! -d './exposure_stats' ]; then
	 mkdir exposure_stats
	fi
R CMD BATCH PREP.exposure_IM.R
echo '=> All exposure prep   done'
echo  >> XS.log

### prepare outcome
	if [ ! -d './outcome_stats' ]; then
	 mkdir outcome_stats
	fi
sh PREP.outcome_MDD.bash
R CMD BATCH PREP.outcome_MDD.R
echo '=> All exposure prep   done'
echo  >> XS.log

### analysis
R CMD BATCH ANALY.MR_IMtoMDD.R
echo '=> Analyses done: check ../results/IM_to_MDD'
echo  >> XS.log

### clean up redundant files
echo Finished at: $(date)  >> XS.log
rm *.Rout
