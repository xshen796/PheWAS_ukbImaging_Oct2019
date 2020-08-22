library("data.table")
library('dplyr')

ls.IM=as.vector(read.table('ls.IMtraits',header = F,stringsAsFactors = F)[,1])

for (f in ls.IM){
      # file names
      outcome.fname=paste0('./outcome_stats/',f,'.MDDoutcome_match')
      traitA=f
      traitB='MDD'

      ## outcome dat
      outcome_match=fread(outcome.fname,header=F,colClasses=NULL,sep='\t')
      outcome_match=outcome_match[,c(1:4,8:11)]
      colnames(outcome_match)=c("SNP","effect_allele","other_allele","eaf","beta","se","pval",'direction')
      outcome_match$Phenotype=traitB
      outcome_match$samplesize=785581
      #outcome_match=filter(outcome_match,(direction=='+++')|(direction=='---'))
      
      write.table(outcome_match,file=paste0('outcome_stats/',traitA,'.',traitB,".outcome_dat"),col.names=T,row.names = F,sep="\t",quot=F)
      system(paste0('echo ',grep(f,ls.IM),' ','MDD as an outcome for ',traitA,' prep   done >> XS.log'))
      
}

rm(list=ls())