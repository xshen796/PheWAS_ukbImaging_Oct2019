library("data.table")
library('dplyr')

ls.IM=as.vector(read.table('ls.IMtraits',header = F,stringsAsFactors = F)[,1])

for (f in ls.IM){
      # file names
      outcome.fname=paste0('./outcome_stats/',f,'.outcome_match')
      traitA='MDD'
      traitB=f

      ## outcome dat
      outcome_match=fread(outcome.fname,header=T,colClasses=NULL)
      outcome_match=outcome_match[,c(2,4,5,6,9,10,12,13)]
      colnames(outcome_match)=c("SNP","effect_allele","other_allele","eaf","beta","se","pval",'samplesize')
      outcome_match$Phenotype=traitB
      
      write.table(outcome_match,file=paste0('outcome_stats/',traitB,".outcome_dat"),col.names=T,row.names = F,sep="\t",quot=F)
      system(paste0('echo ',grep(f,ls.IM),' ',traitA,'as an outcome prep   done >> XS.log'))
      
}

rm(list=ls())