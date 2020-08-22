# 20190803 XShen

# Basic setups ------------------------------------------------------------



library(dplyr)


# Load data ---------------------------------------------------------------

IDP.phewas=readRDS('data/IDP_phewas_meta.rds')
IDP.long_phewas=readRDS('data/dat_long_phewas_meta.rds')

# Define functions --------------------------------------------------------

source('FUNs/reg_phewasStyle_depreCorr.R')

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
ls.models$covs=gsub(pattern = paste0('\\+genotyping.array\\+',paste0('pc',1:15,collapse = '\\+')),'',
                    ls.models$covs)
ls.models=ls.models[!duplicated(ls.models$dependent),]
ls.models$factor='mh.CIDI.MDD.Severity'
ls.models.tmp=ls.models
ls.models.tmp$factor='mh.PHQ9.Severity'
ls.models=rbind(ls.models,ls.models.tmp)
ls.models.tmp$factor='mh.Depressive_symptoms_current_PHQ4'
ls.models=rbind(ls.models,ls.models.tmp)

# Analysis ----------------------------------------------------------------
library('pbapply')
library('nlme')
result.corr_depre_symptoms=reg_phewasStyle(ls.models,dat_short=targetdata,dat_long=dat_long)

save(result.corr_depre_symptoms,file='result/x.supplementary_materials/Depre_symptoms_corr.RData')
