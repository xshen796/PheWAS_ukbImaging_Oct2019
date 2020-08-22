library('dplyr')
library("data.table")


ls.IM=as.vector(read.table('ls.IMtraits',header = F,stringsAsFactors = F)[,1])

for (f in ls.IM){
      
      summ.stats.fname=paste0('../../GWAS/',f,'.bgenie.QC_forLDSC')
      gwas_summarystat = fread(summ.stats.fname,header=T,colClasses=NULL)
      traitA=f
      
      ### find top SNPs
      SigniSNP=filter(gwas_summarystat,af>0.005,af<0.995,info>0.1,hwe>10^-6,N>=5000,p<8*10^-6)
      colnames(SigniSNP)=c('CHR','SNP','pos','Allele1','Allele2','Freq1','info','hwe','Effect','se','t','p','N')
      
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
      SigniSNP=SigniSNP[,c('SNP','NewAllele1','NewAllele2','NewFreq1','NewEffect','se','p','N')]
      
      colnames(SigniSNP)=c("SNP","effect_allele","other_allele","eaf", "beta","se","pval",'samplesize')
      SigniSNP=data.frame(Phenotype=traitA,units='unit',SigniSNP)
      
      write.table(SigniSNP,file=paste0('./exposure_stats/',traitA,".tops"),col.names=T,row.names = F,sep="\t",quot=F)
      
      
      ### clump exposure dat
      
      if (!require(TwoSampleMR)) {
            install.packages("devtools")
            library(devtools)
            install_github("MRCIEU/TwoSampleMR",lib.loc="/exports/igmm/eddie/GenScotDepression/shen/Rlibrary")
      }
      
      
      library(TwoSampleMR)
      #read broad as exposure
      exposure_dat <- read_exposure_data(paste0('./exposure_stats/',traitA,".tops"),sep='\t')
      #clump
      exposure_dat <- clump_data(exposure_dat,clump_kb=3000,clump_r2=0.001,clump_p1=1, clump_p2=1)
      #exposure_dat <- filter(exposure_dat,eaf.exposure>0.1,eaf.exposure<0.90,abs(beta.exposure)<0.038)
      exposure_dat <- exposure_dat[,c('SNP','effect_allele.exposure','other_allele.exposure',
                                      'eaf.exposure','beta.exposure','se.exposure','pval.exposure','samplesize.exposure',
                                      'units.exposure','exposure','mr_keep.exposure','pval_origin.exposure',
                                      'units.exposure_dat','id.exposure','data_source.exposure')]
      write.table(exposure_dat,file=paste0('./exposure_stats/',traitA,".exposure_dat"),col.names=T,row.names = F,sep="\t",quot=F)
      
      system(paste0('echo ',grep(f,ls.IM),' ',traitA,' as an exposure prep done >> XS.log'))
      system(paste0('wc -l ./exposure_stats/',traitA,'.exposure_dat >> XS.log'))
      
}

rm(list=ls())