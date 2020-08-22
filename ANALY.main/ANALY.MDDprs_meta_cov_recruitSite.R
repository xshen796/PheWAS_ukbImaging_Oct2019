# 20190716 XShen

# Basic setups ------------------------------------------------------------

library(dplyr)


# Load data ---------------------------------------------------------------

IDP.phewas=readRDS('data/IDP_phewas_meta.rds')
count=data.frame(pheno=colnames(IDP.phewas),N=colSums(!is.na(IDP.phewas)))
IDP.phewas=IDP.phewas[,count$N>=2000]

IDP.long_phewas=readRDS('data/dat_long_phewas_meta.rds')

recruitment = readRDS('/sdata/images/projects/UKBIOBANK/data/phenotypes/fields/2018-10-phenotypes-ukb24262/Recruitment.rds')
assess_centre_instance0 = data.frame(f.eid=recruitment$f.eid, assessment_centre_instance0 = recruitment$f.54.0.0)
assess_centre_instance0$assessment_centre_instance0 = as.factor(assess_centre_instance0$assessment_centre_instance0)

IDP.phewas=merge(IDP.phewas,assess_centre_instance0,by='f.eid',all.x=T)
IDP.long_phewas=merge(IDP.long_phewas,assess_centre_instance0,by='f.eid',all.x=T)
rm(assess_centre_instance0)

# Define functions --------------------------------------------------------

source('FUNs/reg_phewasStyle.R')

# Define global vars ------------------------------------------------------

targetdata=IDP.phewas
dat_long=IDP.long_phewas

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
ls.dep.subcor=c(colnames(targetdata)[grep('^ICV',colnames(targetdata))],
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

# covs
ls.models=data.frame(dependent=ls.dep.factor.combo$Var1,
                     factor=ls.dep.factor.combo$Var2,
                     covs='',stringsAsFactors = F)
#wm
ls.models$covs[grep('^bl.FA|^bl.MD|^bl.ICVF|^bl.ISOVF|^bl.OD|^g\\.',ls.models$dependent)]=paste0(c(colnames(targetdata)[4:26],'assessment_centre_instance0'),collapse = '+')
ls.models$covs[grep('^FA\\.wm|^MD\\.wm|^ICVF\\.wm|^ISOVF\\.wm|^OD\\.wm',ls.models$dependent)]=paste0(c(colnames(targetdata)[4:26],'hemi',
                                                                                                       'assessment_centre_instance0'),collapse = '+')
#subcor
ls.models$covs[grep('^ICV$',ls.models$dependent)]=paste0(c(colnames(targetdata)[4:26],'assessment_centre_instance0'),collapse = '+')
ls.models$covs[grep('^vol\\.sub\\.',ls.models$dependent)]=paste0(c(colnames(targetdata)[c(4:26)],'ICV','hemi',
                                                                   'assessment_centre_instance0'),collapse = '+')
ls.models$covs[grep('^bl.vol\\.sub\\.',ls.models$dependent)]=paste0(c(colnames(targetdata)[c(4:26)],'ICV',
                                                                      'assessment_centre_instance0'),collapse = '+')
#T2star
ls.models$covs[grep('^T2\\.star\\.',ls.models$dependent)]=paste0(c(colnames(targetdata)[c(4:26)],'ICV','hemi',
                                                                   'assessment_centre_instance0'),collapse = '+')
#rsfMRI
ls.models$covs[grep('^rsconn|^rsamp',ls.models$dependent)]=paste0(c(colnames(targetdata)[4:26],'r.motion',
                                                                    'assessment_centre_instance0'),collapse = '+')
# the rest - behavioural
ls.models$covs[ls.models$covs=='']=paste0(c(colnames(targetdata)[c(4,8:26)],'assessment_centre_instance0'),collapse = '+')

# model.est
ls.models$model.est=''
ls.models$model.est[grep('hemi',ls.models$covs)]='lme'
ls.models$model.est[ls.models$model.est=='']='glm'


# Analysis ----------------------------------------------------------------
library('pbapply')
library('nlme')
result.cov_recruitSite=reg_phewasStyle(ls.models,dat_short=targetdata,dat_long=dat_long)

save(result.cov_recruitSite,file='result/i.Main_result/result_meta_covRecruitmentCentre.RData')
