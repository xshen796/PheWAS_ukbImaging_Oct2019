# 20190803 XShen

# Basic setups ------------------------------------------------------------



library(dplyr)


# Load data ---------------------------------------------------------------

IDP.phewas=readRDS('data/IDP_phewas_meta.rds')
IDP.long_phewas=readRDS('data/dat_long_phewas_meta.rds')

# Define functions --------------------------------------------------------

source('FUNs/reg_phewasStyle.R')

# Define global vars ------------------------------------------------------

targetdata=IDP.phewas
dat_long=IDP.long_phewas
load('script/ANALY.main/ls.models.RData')

#### list of dependent variables
# list of significant phenotypes replicated in both samples
load('result/i.Main_result/result_replication.RData')
targetResult=result.replication
mat.p=matrix(targetResult$p.corrected,ncol = 8)
colnames(mat.p)=paste0('mddPRS_pT_',c('0.0005','0.001','0.005','0.01','0.05','0.1','0.5','1'))
rownames(mat.p)=unique(targetResult$dependent)
ls.sigpheno.replication=rowSums(mat.p<0.05)>=4
ls.sigpheno.replication=rownames(mat.p)[ls.sigpheno.replication]
# remove MDD related phenotypes
ls.dep.var=ls.sigpheno.replication[!grepl('^mh.MDD_|^mh.Depressi|PHQ9|MDD.Severity',ls.sigpheno.replication)]

#### new ls.models
ls.models=ls.models[ls.models$dependent %in% ls.dep.var,]
ls.models$covs=paste0(ls.models$covs,'+mh.PHQ9.Severity+mh.CIDI.MDD.Severity+mh.Depressive_symptoms_current_PHQ4')


# Analysis ----------------------------------------------------------------
library('pbapply')
library('nlme')
result.cov_depre_symptoms=reg_phewasStyle(ls.models,dat_short=targetdata,dat_long=dat_long)



# Compare beta between results --------------------------------------------
load('result//i.Main_result//result_meta.RData')
result.no_cov=result.meta[result.meta$dependent %in% ls.dep.var,]
mat.beta.no_cov=matrix(result.no_cov$beta,ncol = 8)
mat.beta.cov_depre_symptoms=matrix(result.cov_depre_symptoms$beta,ncol = 8)
mat.p.cov=matrix(result.cov_depre_symptoms$p.corrected,ncol = 8)


deltaBeta=data.frame(dep=unique(result.no_cov$dependent),meanBeta_noCov=rowMeans(mat.beta.no_cov),
                     meanBeta_Cov_DepreSymptoms=rowMeans(mat.beta.cov_depre_symptoms),
                     p.new=rowSums(mat.p.cov)>=4)
deltaBeta$deltaBeta=(deltaBeta$meanBeta_noCov-deltaBeta$meanBeta_Cov_DepreSymptoms)/deltaBeta$meanBeta_noCov

save(deltaBeta,file='result/x.supplementary_materials/deltaBeta.RData')
save(result.cov_depre_symptoms,file='result//x.supplementary_materials/result.cov_depre_symp.RData')
