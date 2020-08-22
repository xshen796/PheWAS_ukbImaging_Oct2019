library('dplyr')
library("data.table")

gwas_summarystat = fread("3cohorts.meta.txt",
                         header=T,colClasses=NULL)
#colnames(gwas_summarystat) = c('MarkerName','Allele1','Allele2','Freq1','FreqSE','MinFreq','MaxFreq','Effect','Std','P-value','Direction')
colnames(gwas_summarystat)[grep('P-value',colnames(gwas_summarystat))]='P'
SigniSNP=filter(gwas_summarystat,Freq1>0.005,Freq1<0.995,P< 5*10^-8,(Direction=='---')|(Direction=='+++'))
#SigniSNP=filter(gwas_summarystat,Freq1>0.005,Freq1<0.995,P< 5*10^-8)

#harmonization
SigniSNP$NewAllele1=SigniSNP$Allele1
SigniSNP$NewAllele2=SigniSNP$Allele2
SigniSNP$NewFreq1=SigniSNP$Freq1
SigniSNP$NewEffect=SigniSNP$Effect

lst=c(which(SigniSNP$Effect<0))
    SigniSNP$NewEffect[lst]=-1*SigniSNP$Effect[lst]
    SigniSNP$NewFreq1[lst]=1-1*SigniSNP$Freq1[lst]
    SigniSNP$NewAllele1[lst]=SigniSNP$Allele2[lst]
    SigniSNP$NewAllele2[lst]=SigniSNP$Allele1[lst]
SigniSNP=data.frame(SigniSNP)
SigniSNP=SigniSNP[,c('MarkerName','NewAllele1','NewAllele2','NewFreq1','NewEffect','StdErr','P')]

colnames(SigniSNP)=c("SNP","effect_allele","other_allele","eaf", "beta","se","pval")
SigniSNP=data.frame(Phenotype='MDD_23_UKB_PGC',units='unit',N=785581,SigniSNP)

write.table(SigniSNP,file="./exposure_stats/MDD.tops",col.names=T,row.names = F,sep="\t",quot=F)

#####prepare Broad depression as exposure 

if (!require(TwoSampleMR)) {
      install.packages("devtools")
      library(devtools)
      install_github("MRCIEU/TwoSampleMR",lib.loc="/exports/igmm/eddie/GenScotDepression/shen/Rlibrary")
}

library(TwoSampleMR)
#read broad as exposure
exposure_dat <- read_exposure_data("./exposure_stats/MDD.tops",sep = "\t")
#clump
exposure_dat <- clump_data(exposure_dat,clump_kb=3000,clump_r2=0.001,clump_p1=1, clump_p2=1)
exposure_dat <- exposure_dat[,c('SNP','effect_allele.exposure','other_allele.exposure',
                                'eaf.exposure','beta.exposure','se.exposure','pval.exposure',
                                'units.exposure','exposure','mr_keep.exposure','pval_origin.exposure',
                                'units.exposure_dat','id.exposure','data_source.exposure')]
write.table(exposure_dat,file="./exposure_stats/MDD.exposure_dat",col.names=T,row.names = F,sep="\t",quot=F)


system('echo MDD as an exposure prep done >> XS.log')

rm(list=ls())