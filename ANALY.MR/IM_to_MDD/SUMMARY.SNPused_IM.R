library("data.table")
library('dplyr')
library('knitr')
library('TwoSampleMR')
library('MRPRESSO')
library("ggplot2")
options(bitmapType='cairo')

ls.IM='tMD_ATR'

for (f in ls.IM){
      # file names
      exposure.fname=paste0('exposure_stats/',f,'.exposure_dat')
      traitA=f
      traitB='MDD'
      
      ## exposure dat
      exposure_dat=read.table(exposure.fname,header=T,sep='\t')
      ## outcome dat
      outcome_dat <- read_outcome_data(file=paste0('outcome_stats/',traitA,'.',traitB,".outcome_dat"),sep='\t')
      
      # Harmonise the exposure and outcome data
      dat <- harmonise_data(exposure_dat, outcome_dat)
      
      dat = dat[,c('SNP','effect_allele.exposure','other_allele.exposure','eaf.exposure','beta.exposure',
                   'se.exposure','pval.exposure','beta.outcome','se.outcome','pval.outcome')]
      dat$beta.exposure = paste0(format(dat$beta.exposure, scientific=T, digits = 2),' (',
                                 format(dat$se.exposure, scientific=T, digits = 2),')')
      dat$beta.outcome = paste0(format(dat$beta.outcome, scientific=T, digits = 2),' (',
                                 format(dat$se.outcome, scientific=T, digits = 2),')')
      dat$pval.exposure = format(dat$pval.exposure, scientific=T, digits = 2)
      dat$pval.outcome = format(dat$pval.outcome, scientific=T, digits = 2)
      
      dat=dat[,c('SNP','effect_allele.exposure','other_allele.exposure','eaf.exposure','beta.exposure',
                 'pval.exposure','beta.outcome','pval.outcome')]
      colnames(dat)[grep('exposure',colnames(dat))]=gsub('exposure','depression',colnames(dat)[grep('exposure',colnames(dat))])
      colnames(dat)[grep('outcome',colnames(dat))]=gsub('outcome',traitB,colnames(dat)[grep('outcome',colnames(dat))])
      
      
     
      # summarise stats
      if (f==ls.IM[1]){
            SNPused.dat=dat
      }else{
            SNPused.dat=cbind(SNPused.dat,dat[,grep(traitB,colnames(dat))])
      }
     
}
write.table(SNPused.dat,file="../results/IM_to_MDD/SNPused_IM.summary",col.names = T,row.names=F,sep="\t",quot=F)
rm(list=ls())