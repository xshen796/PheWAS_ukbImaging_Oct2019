# 20190803 XShen

# Basic setups ------------------------------------------------------------
library(dplyr)


# Load data ---------------------------------------------------------------

IDP.phewas=readRDS('data/IDP_phewas_meta.rds')
IDP.long_phewas=readRDS('data/dat_long_phewas_meta.rds')
IDP.long_to_short=(IDP.long_phewas[1:(nrow(IDP.long_phewas)/2),2:75]+
                         IDP.long_phewas[(nrow(IDP.long_phewas)/2+1):nrow(IDP.long_phewas),2:75])/2


# Define functions --------------------------------------------------------

source('FUNs/reg_phewasStyle_interaction.R')

# Define global vars ------------------------------------------------------

targetdata=cbind(IDP.phewas,IDP.long_to_short)
dat_long=IDP.long_phewas

# E factors
ls.Env=c('sex')

# dependent variables
ls.dep.behav=colnames(targetdata)[grep('^early|^sociodemo|^lifestyle|^physical|^cog\\.|^mh\\.',colnames(targetdata))]
ls.dep.wm=c(colnames(targetdata)[grep('^g.FA|^g.MD|^g.ICVF|^g.ISOVF|^g.OD',colnames(targetdata))],
            colnames(dat_long)[grep('^FA.wm',colnames(dat_long))],            
            colnames(targetdata)[grep('^bl.FA.wm',colnames(targetdata))],
            colnames(dat_long)[grep('^MD.wm',colnames(dat_long))],            
            colnames(targetdata)[grep('^bl.MD.wm',colnames(targetdata))],
            colnames(dat_long)[grep('^ICVF.wm',colnames(dat_long))],            
            colnames(targetdata)[grep('^bl.ICVF.wm',colnames(targetdata))],
            colnames(dat_long)[grep('^ISOVF.wm',colnames(dat_long))],            
            colnames(targetdata)[grep('^bl.ISOVF.wm',colnames(targetdata))],
            colnames(dat_long)[grep('^OD.wm',colnames(dat_long))],            
            colnames(targetdata)[grep('^bl.OD.wm',colnames(targetdata))])
ls.dep.subcor=c(colnames(targetdata)[grep('^ICV$',colnames(targetdata))],
                colnames(dat_long)[grep('^vol.sub.',colnames(dat_long))],
                colnames(targetdata)[grep('^bl.vol.sub.brain_stem',colnames(targetdata))])
ls.T2star=c(colnames(targetdata)[grep('^Total_volume_of_white_matter_hyperintensities',colnames(targetdata))],
            colnames(dat_long)[grep('T2.star.',colnames(dat_long))])
ls.dep.rsconn=colnames(targetdata)[grep('^rsconn',colnames(targetdata))]
ls.dep.rsamp=colnames(targetdata)[grep('^rsamp',colnames(targetdata))]

ls.dep.all=c(ls.dep.behav,ls.T2star,ls.dep.subcor,ls.dep.wm,ls.dep.rsconn,ls.dep.rsamp)

# factors
ls.factor=paste0('meta3cohorts_pT_',c('0.0005','0.001','0.005','0.01','0.05','0.1','0.5','1'))

# combine the two
ls.dep.factor.combo=expand.grid(ls.dep.all,ls.factor,stringsAsFactors = F)


create_covs <- function(ls.keyvars,targetdata,env_var) {
      
      dep = ls.keyvars[1]
      gene = ls.keyvars[2]
      cov.type = as.numeric(ls.keyvars[3])
      
      ls.covs = switch(
            cov.type,
            colnames(targetdata)[c(4:9,11:26)],
            c(colnames(targetdata)[c(4:9,11:26)]),
            colnames(targetdata)[c(4:9,11:26)],
            c(colnames(targetdata)[c(4:9,11:26)],'ICV'),
            c(colnames(targetdata)[c(4:9,11:26)],'ICV'),
            c(colnames(targetdata)[c(4:9,11:26)],'ICV'),
            c(colnames(targetdata)[c(4:9,11:26)],'r.motion'),
            colnames(targetdata)[c(4,8:9,11:26)]
            )
      
      covs.w_g = paste0(ls.covs,'*',gene,collapse = '+')
      covs.w_env = paste0(ls.covs,'*',env_var,collapse = '+')
      covs.sum = paste0(covs.w_g,'+',covs.w_env)
      return(covs.sum)
}

load('script/ANALY.PRS_by_E/ls.sigpheno.replication.RData')


# Analysis ----------------------------------------------------------------
library('pbapply')
library('nlme')
for (i in ls.Env){
      ls.models=data.frame(dependent=ls.dep.factor.combo$Var1,
                           factor=ls.dep.factor.combo$Var2,
                           cov_type=999,covs='',stringsAsFactors = F)
      #wm
      ls.models$cov_type[grep('^bl.FA|^bl.MD|^bl.ICVF|^bl.ISOVF|^bl.OD|^g\\.',ls.models$dependent)]=1
      ls.models$cov_type[grep('^FA\\.wm|^MD\\.wm|^ICVF\\.wm|^ISOVF\\.wm|^OD\\.wm',ls.models$dependent)]=2
      #subcor
      ls.models$cov_type[grep('^ICV$',ls.models$dependent)]=3
      ls.models$cov_type[grep('^vol\\.sub\\.',ls.models$dependent)]=4
      ls.models$cov_type[grep('^bl.vol\\.sub\\.',ls.models$dependent)]=5
      #T2star
      ls.models$cov_type[grep('^T2\\.star\\.',ls.models$dependent)]=6
      #rsfMRI
      ls.models$cov_type[grep('^rsconn|^rsamp',ls.models$dependent)]=7
      # the rest - behavioural
      ls.models$cov_type[ls.models$cov_type==999]=8
      
      ls.models=ls.models[!(ls.models$dependent %in% ls.Env),]
      ls.models=ls.models[ls.models$dependent %in% ls.sigpheno.replication,]
      
      covs_formodel=apply(X = ls.models[,1:3],MARGIN = 1,FUN = create_covs, 
                          targetdata=targetdata, env_var=i)
      covs_formodel=as.data.frame(covs_formodel)
      ls.models$covs=covs_formodel$covs_formodel
      
      # model.est
      ls.models$model.est=''
      ls.models$model.est[grep('hemi',ls.models$covs)]='lme'
      ls.models$model.est[ls.models$model.est=='']='glm'
      
      ls.models=ls.models[,!grepl('cov_type',colnames(ls.models))]
      
      result.tmp=reg_phewasStyle_interaction(ls.models,dat_short=targetdata,dat_long=NA,interact_var = i)
      eval(parse(text=paste0('result.',i,'=result.tmp')))
      cat(paste0(i,'\n'))
      obj.name=paste0('result.',i)
      save(list=obj.name,file=paste0('result/iv.PRS_by_E/PRS_X_E_',i,'.RData'))
}