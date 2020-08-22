# 20190627 XShen

# Basic setups ------------------------------------------------------------



library(dplyr)


# Load data ---------------------------------------------------------------

IDP.phewas=readRDS('data/IDP_phewas_replication.rds')
IDP.long_phewas=readRDS('data/dat_long_phewas_replication.rds')

# Define functions --------------------------------------------------------

source('FUNs/reg_phewasStyle.R')

# Define global vars ------------------------------------------------------

targetdata=IDP.phewas
dat_long=IDP.long_phewas
load('script/ANALY.main/ls.models.RData')

# list of significant phenotypes from the discovery sample
load('result/i.Main_result/result_discovery.RData')
targetResult=result.discovery
      mat.p=matrix(targetResult$p.corrected,ncol = 8)
      colnames(mat.p)=paste0('mddPRS_pT_',c('0.0005','0.001','0.005','0.01','0.05','0.1','0.5','1'))
      rownames(mat.p)=unique(targetResult$dependent)
p.discovery=mat.p
ls.sigpheno.discovery=rowSums(mat.p<0.05)>=4
ls.sigpheno.discovery=rownames(mat.p)[ls.sigpheno.discovery]

ls.models=ls.models[ls.models$dependent %in% ls.sigpheno.discovery,]


# Analysis ----------------------------------------------------------------
library('pbapply')
library('nlme')
result.replication=reg_phewasStyle(ls.models,dat_short=targetdata,dat_long=dat_long)

save(result.replication,file='result/i.Main_result/result_replication.RData')
