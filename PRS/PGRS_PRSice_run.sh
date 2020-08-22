#!/bin/sh
#$ -N MEGApgrs
#$ -cwd
#$ -m beas
#$ -l h_vmem=128G
#$ -l h_rt=24:00:00
. /etc/profile.d/modules.sh

module add igmm/apps/R/3.2.4


R --file=PRSice_v1.25.R -q --args \
plink plink_linux_x86_64/plink \
base 3cohortsmeta.txt \
target autosome \
clump.p1 1 \
clump.p2 1 \
clump.r2 0.25 \
clump.kb 250 \
fastscore T \
barchart.levels 0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1 \
report.individual.scores T \
report.best.score.only F \
covary F \
no.regression F \
debug.mode T \
plink.silent F \
binary.target T \
figname PGRS \
pheno.file dummy.pheno.MDDpgrs \
wd prsice.output/

