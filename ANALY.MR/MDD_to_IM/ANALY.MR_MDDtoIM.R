library("data.table")
library('dplyr')
library('knitr')
library('TwoSampleMR')
library('MRPRESSO')
library("ggplot2")
options(bitmapType='cairo')

ls.IM=as.vector(read.table('ls.IMtraits',header = F,stringsAsFactors = F)[,1])

for (f in ls.IM){
      # file names
      exposure.fname=paste0('exposure_stats/MDD.exposure_dat')
      traitA='MDD'
      traitB=f
      
      ## exposure dat
      exposure_dat=read.table(exposure.fname,header=T)
      ## outcome dat
      outcome_dat <- read_outcome_data(file=paste0('outcome_stats/',traitB,".outcome_dat"),sep='\t')
      
      # Harmonise the exposure and outcome data
      dat <- harmonise_data(exposure_dat, outcome_dat)
      
      # analyse
      MR=mr(dat = dat,method_list =c("mr_ivw","mr_egger_regression","mr_weighted_median"))
      res_single = mr_singlesnp(dat)
      res_loo = mr_leaveoneout(dat)
      MR.presso=mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                          SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                          OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
                          data = dat, NbDistribution = 3000,  SignifThreshold = 0.05)
      if (MR.presso$`MR-PRESSO results`$`Global Test`$Pvalue>=0.05){
            tmp.mr_res=data.frame(MR,mr_pleiotropy_test(dat)[5:7],mr_heterogeneity(dat)[2,6:8],
                                  MRpresso_RSSobs=MR.presso$`MR-PRESSO results`$`Global Test`$RSSobs,
                                  MRpresso_Pval=as.character(MR.presso$`MR-PRESSO results`$`Global Test`$Pvalue),
                                  MRpresso_outliers=NA,MRpresso_distortion_pval=NA)
      }else{
            n.outlier=length(MR.presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)
            MR=rbind(MR,MR[1,])
            MR$method=as.character(MR$method)
            MR$method[4]='MR-Presso_outlier-corrected'
            MR$nsnp[4]=MR$nsnp[4]-n.outlier
            MR[4,7:9]=MR.presso$`Main MR results`[2,c(3,4,6)]
            tmp.mr_res=data.frame(MR,mr_pleiotropy_test(dat)[5:7],mr_heterogeneity(dat)[2,6:8],
                                  MRpresso_RSSobs=MR.presso$`MR-PRESSO results`$`Global Test`$RSSobs,
                                  MRpresso_Pval=as.character(MR.presso$`MR-PRESSO results`$`Global Test`$Pvalue),
                                  MRpresso_outliers=paste0(MR.presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`,collapse=','),
                                  MRpresso_distortion_pval=MR.presso$`MR-PRESSO results`$`Distortion Test`$Pvalue)
      }
      # summarise stats
      if (f==ls.IM[1]){
            MRsummary=tmp.mr_res
      }else{
            MRsummary=rbind(MRsummary,tmp.mr_res)
      }
      
      # visualise
      mr_report(dat,output_path='../results/MDD_to_IM/')
      
      fig.scatter = mr_scatter_plot(MR, dat)
      fig.forest = mr_forest_plot(res_single)
      fig.loo = mr_leaveoneout_plot(res_loo)
      fig.funnel = mr_funnel_plot(res_single)
      
      ggsave(fig.scatter[[1]], file=paste0("../results/MDD_to_IM/figure/",traitA,'_',traitB,"_scatter_plot.png"), width=7, height=7)
      ggsave(fig.forest[[1]], file=paste0("../results/MDD_to_IM/figure/",traitA,'_',traitB,"_forestr_plot.png"), width=7, height=7)
      ggsave(fig.loo[[1]], file=paste0("../results/MDD_to_IM/figure/",traitA,'_',traitB,"_loo_plot.png"), width=7, height=7)
      ggsave(fig.funnel[[1]], file=paste0("../results/MDD_to_IM/figure/",traitA,'_',traitB,"_funnel_plot.png"), width=7, height=7)     
      system(paste0('echo ',grep(f,ls.IM),' ',traitA,' against ',traitB,'analysis  done >> XS.log'))
}
write.table(MRsummary,file="../results/MDD_to_IM/outcomes.MRsummary",col.names = T,row.names=F,sep="\t",quot=F)
rm(list=ls())