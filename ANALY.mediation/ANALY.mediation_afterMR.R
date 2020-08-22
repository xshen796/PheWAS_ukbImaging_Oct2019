library(lavaan)
library(dplyr)

mediation.dat=readRDS('data/mediation.dat.rds')
targetdata=filter(mediation.dat,!is.na(relaQC))

load('script/ANALY.mediation/ls.sig.MR.RData')

# scale pos and r.motion
targetdata$sex=as.numeric(targetdata$sex)
targetdata$genotyping.array=as.numeric(targetdata$genotyping.array)
targetdata$MRI_site=as.numeric(targetdata$MRI_site)
targetdata[,c(4:27)]=scale(targetdata[,c(4:27)])

############################################################################
##########################  MDDpgrs - MDD - IM  ############################
############################################################################

###### basic lists
ls.covs=3:25

x.ls=colnames(targetdata)[grep('meta3cohorts',colnames(targetdata))]
m.ls=colnames(targetdata)[c(36:39)]  #147
y.ls=colnames(targetdata)[c(30,31,28,33)]#46


###### make mediation model list
ls.total=data.frame(predictor=rep(x.ls,each=length(m.ls)),
                    medi=m.ls)
ls.total=data.frame(ls.total,outc=rep(y.ls,each=(nrow(ls.total))))

targetdata$mh.MDD_CIDI=as.numeric(targetdata$mh.MDD_CIDI)

for (c.ls in 1:nrow(ls.total)){
      
      x.n = ls.total[c.ls,'predictor']
      m.n = ls.total[c.ls,'medi']
      y.n = ls.total[c.ls,'outc']
      
      expre1=paste0('X=targetdata$',x.n)
      expre2=paste0('M=targetdata$',m.n)
      expre3=paste0('Y=targetdata$',y.n)
      eval(parse(text=expre1))
      eval(parse(text=expre2))
      eval(parse(text=expre3))
      
      set.seed(1234)
      Data <- data.frame(X = X, Y = Y, M = M,targetdata[,ls.covs])
      Data$genotyping.array=as.numeric(Data$genotyping.array)
      Data$sex=as.numeric(Data$sex)
      Data$age=Data$MRI_age.calculated
      Data=scale(Data)
      Data=Data[complete.cases(Data),]
      
      
      if (length(grep('FA\\.|MD\\.',m.n))==1){
            m.model='M~age+I(age^2)+sex+pos.x+pos.y+pos.z'
      }else if (length(grep('^N',m.n))==1){
            m.model='M~age+I(age^2)+sex+pos.x+pos.y+pos.z+r.motion'
      }else {
            m.model='M~age+I(age^2)+sex'
      }
      
      if (length(grep('FA\\.|MD\\.',y.n))==1){
            y.model='Y~age+I(age^2)+sex+pos.x+pos.y+pos.z'
      }else if (length(grep('^N',y.n))==1){
            y.model='Y~age+I(age^2)+sex+pos.x+pos.y+pos.z+r.motion'
      }else {
            y.model='Y~age+I(age^2)+sex'
      }
      
      model <- paste0(' # direct effect
                      Y ~ c*X
                      # mediator
                      M ~ a*X 
                      Y ~ b*M
                      # total effect
                      indirect := a*b
                      direct := c
                      total := c + (a*b)
                      
                      # control for age and sex
                      X ~ age+I(age^2)+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+genotyping.array
                      \n',
                      m.model,'\n',
                      y.model,'\n'
      )
      fit <- sem(model, data = Data)
      fits=fitMeasures(fit,c('cfi','tli','rmsea','rmsea.ci.lower','rmsea.ci.upper','rmsea.pvalue'))
      para.fits=parameterEstimates(fit)
      fits=c((para.fits[nrow(para.fits)-2,5]/para.fits[nrow(para.fits),5]),fits,para.fits[3,8])
      fits=cbind(rbind(ls.total[c.ls,],rep(NA,3),rep(NA,3)),
                 para.fits[(nrow(para.fits)-2):(nrow(para.fits)),3:ncol(para.fits)],
                 rbind(fits,rep(NA,8),rep(NA,8)))
      
      if (c.ls==1){
            fits.total=fits
      }else{
            fits.total=rbind(fits.total,fits)
      }
      cat(paste0(c.ls,'\n'))
}
rownames(fits.total)=NULL
tmp.res=fits.total
tmp.res=tmp.res[!is.na(tmp.res$predictor),]
tmp.res=filter(tmp.res,cfi>=0.9,tli>=0.9,rmsea.pvalue>0.05)#,V8<0.05


##### correct p values
ls.correct=data.frame(predictor=rep(x.ls[order(x.ls)],each=length(m.ls)),
                      medi=m.ls[order(m.ls)],stringsAsFactors = F)
tmp.res.tocorrect=tmp.res[order(tmp.res$predictor,tmp.res$medi),]

for (i in 1:nrow(ls.correct)){
  p=ls.correct$predictor[i]
  m=ls.correct$medi[i]
    tmp.block=tmp.res.tocorrect[grep(p,tmp.res.tocorrect$predictor),]
    tmp.block=tmp.block[grep(m,tmp.block$medi),]
    p.g=p.adjust(tmp.block$pvalue[grep('^g.MD',tmp.block$outc)],method='fdr')
    p.t=p.adjust(tmp.block$pvalue[grep('^ICVF',tmp.block$outc)],method='fdr')
    p.rs=p.adjust(tmp.block$pvalue[grep('^rsamp',tmp.block$outc)],method='fdr')
  tmp.p.all=c(p.g,p.t,p.rs)
  
  if (i==1){p.corrected=tmp.p.all}else{p.corrected=c(p.corrected,tmp.p.all)}
}

tmp.res.tocorrect$p.corrected=p.corrected
MDDpgrs_MDD_IM=tmp.res.tocorrect[,c(1:3,6:9,20,12:15,18)]
colnames(MDDpgrs_MDD_IM)[grep('V1',colnames(MDDpgrs_MDD_IM))]='C_change'
MDDpgrs_MDD_IM$predictor=gsub('meta3cohorts','Depression-PRS',MDDpgrs_MDD_IM$predictor)
MDDpgrs_MDD_IM$outc=gsub('^MD.wm.','MD in ',MDDpgrs_MDD_IM$outc)
MDDpgrs_MDD_IM$outc=gsub('^g.MD.','g.MD-',MDDpgrs_MDD_IM$outc)
MDDpgrs_MDD_IM$outc=gsub('^N14','Amplitude of Salience Network (N14)',MDDpgrs_MDD_IM$outc)
MDDpgrs_MDD_IM$medi=gsub('^mh.','',MDDpgrs_MDD_IM$medi)
MDDpgrs_MDD_IM$medi=gsub('MDD_CIDI','CIDI depression',MDDpgrs_MDD_IM$medi)
MDDpgrs_MDD_IM$medi=gsub('CIDI.MDD.Severity','Severity of depression (CIDI)',MDDpgrs_MDD_IM$medi)
MDDpgrs_MDD_IM$medi=gsub('Depressive_symptoms_current_PHQ4','Current depressive symptoms',MDDpgrs_MDD_IM$medi)
MDDpgrs_MDD_IM[,6:12]=round(MDDpgrs_MDD_IM[,6:12],3)





save(MDDpgrs_MDD_IM,file='result/iii.mediation/MEDIATION.afterMR.MDDpgrs_MDD_IM.RData')
#write.table(MDDpgrs_MDD_IM,file='table/MEDIATION.afterMR.MDDpgrs_MDD_IM.txt',quote = F,sep = '\t',row.names = F)





############################################################################
##########################  MDDpgrs - IM - MDD  ############################
############################################################################

###### basic lists
ls.covs=3:25

x.ls=colnames(targetdata)[grep('meta3cohorts',colnames(targetdata))]
y.ls=colnames(targetdata)[c(36:39)]  #147
m.ls=colnames(targetdata)[c(32,28,29)]#46


###### make mediation model list
ls.total=data.frame(predictor=rep(x.ls,each=length(m.ls)),
                    medi=m.ls)
ls.total=data.frame(ls.total,outc=rep(y.ls,each=(nrow(ls.total))))

targetdata$mh.MDD_CIDI=as.numeric(targetdata$mh.MDD_CIDI)

for (c.ls in 1:nrow(ls.total)){
  
  x.n = ls.total[c.ls,'predictor']
  m.n = ls.total[c.ls,'medi']
  y.n = ls.total[c.ls,'outc']
  
  expre1=paste0('X=targetdata$',x.n)
  expre2=paste0('M=targetdata$',m.n)
  expre3=paste0('Y=targetdata$',y.n)
  eval(parse(text=expre1))
  eval(parse(text=expre2))
  eval(parse(text=expre3))
  
  set.seed(1234)
  Data <- data.frame(X = X, Y = Y, M = M,targetdata[,ls.covs])
  Data$genotyping.array=as.numeric(Data$genotyping.array)
  Data$sex=as.numeric(Data$sex)
  Data$age=Data$MRI_age.calculated
  Data=scale(Data)
  Data=Data[complete.cases(Data),]
  
  
  if (length(grep('FA\\.|MD\\.',m.n))==1){
    m.model='M~age+I(age^2)+sex+pos.x+pos.y+pos.z'
  }else if (length(grep('^N',m.n))==1){
    m.model='M~age+I(age^2)+sex+pos.x+pos.y+pos.z+r.motion'
  }else {
    m.model='M~age+I(age^2)+sex'
  }
  
  if (length(grep('FA\\.|MD\\.',y.n))==1){
    y.model='Y~age+I(age^2)+sex+pos.x+pos.y+pos.z'
  }else if (length(grep('^N',y.n))==1){
    y.model='Y~age+I(age^2)+sex+pos.x+pos.y+pos.z+r.motion'
  }else {
    y.model='Y~age+I(age^2)+sex'
  }
  
  model <- paste0(' # direct effect
                  Y ~ c*X
                  # mediator
                  M ~ a*X 
                  Y ~ b*M
                  # total effect
                  indirect := a*b
                  direct := c
                  total := c + (a*b)
                  
                  # control for age and sex
                  X ~ age+I(age^2)+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+genotyping.array
                  \n',
                  m.model,'\n',
                  y.model,'\n'
  )
  fit <- sem(model, data = Data)
  fits=fitMeasures(fit,c('cfi','tli','rmsea','rmsea.ci.lower','rmsea.ci.upper','rmsea.pvalue'))
  para.fits=parameterEstimates(fit)
  fits=c((para.fits[nrow(para.fits)-2,5]/para.fits[nrow(para.fits),5]),fits,para.fits[3,8])
  fits=cbind(rbind(ls.total[c.ls,],rep(NA,3),rep(NA,3)),
             para.fits[(nrow(para.fits)-2):(nrow(para.fits)),3:ncol(para.fits)],
             rbind(fits,rep(NA,8),rep(NA,8)))
  
  if (c.ls==1){
    fits.total=fits
  }else{
    fits.total=rbind(fits.total,fits)
  }
  cat(paste0(c.ls,'\n'))
}
rownames(fits.total)=NULL
tmp.res=fits.total
tmp.res=tmp.res[!is.na(tmp.res$predictor),]
tmp.res=filter(tmp.res,cfi>=0.9,tli>=0.9,rmsea.pvalue>0.05)#,V8<0.05


##### correct p values
ls.correct=data.frame(predictor=rep(x.ls[order(x.ls)],each=length(m.ls)),
                      medi=m.ls[order(m.ls)],stringsAsFactors = F)
tmp.res.tocorrect=tmp.res[order(tmp.res$predictor,tmp.res$outc),]

for (i in 1:nrow(ls.correct)){
      p=ls.correct$predictor[i]
      m=ls.correct$medi[i]
      tmp.block=tmp.res.tocorrect[grep(p,tmp.res.tocorrect$predictor),]
      tmp.block=tmp.block[grep(m,tmp.block$medi),]
      p.g=p.adjust(tmp.block$pvalue[grep('^MD',tmp.block$medi)],method='fdr')
      p.t=p.adjust(tmp.block$pvalue[grep('^ICVF',tmp.block$medi)],method='fdr')
      p.rs=p.adjust(tmp.block$pvalue[grep('^g.ICVF',tmp.block$medi)],method='fdr')
      tmp.p.all=c(p.g,p.t,p.rs)
      
      if (i==1){p.corrected=tmp.p.all}else{p.corrected=c(p.corrected,tmp.p.all)}
}
# 
tmp.res.tocorrect$p.corrected=p.corrected
MDDpgrs_IM_MDD=tmp.res.tocorrect[,c(1:3,6:9,20,12:15,18)]
colnames(MDDpgrs_IM_MDD)[grep('V1',colnames(MDDpgrs_IM_MDD))]='C_change'

MDDpgrs_IM_MDD$predictor=gsub('meta3cohorts','Depression-PRS',MDDpgrs_IM_MDD$predictor)
MDDpgrs_IM_MDD$medi=gsub('^MD.wm.','MD in ',MDDpgrs_IM_MDD$medi)
MDDpgrs_IM_MDD$medi=gsub('^g.MD.','g.MD-',MDDpgrs_IM_MDD$medi)
MDDpgrs_IM_MDD$medi=gsub('^N14','Amplitude of Salience Network (N14)',MDDpgrs_IM_MDD$medi)
MDDpgrs_IM_MDD$outc=gsub('^mh.','',MDDpgrs_IM_MDD$outc)
MDDpgrs_IM_MDD$outc=gsub('MDD_CIDI','CIDI depression',MDDpgrs_IM_MDD$outc)
MDDpgrs_IM_MDD$outc=gsub('CIDI.MDD.Severity','Severity of depression (CIDI)',MDDpgrs_IM_MDD$outc)
MDDpgrs_IM_MDD$outc=gsub('Depressive_symptoms_current_PHQ4','Current depressive symptoms',MDDpgrs_IM_MDD$outc)
MDDpgrs_IM_MDD[,6:12]=round(MDDpgrs_IM_MDD[,6:12],3)


save(MDDpgrs_IM_MDD,file='result/iii.mediation/MEDIATION.afterMR.MDDpgrs_IM_MDD.RData')
#write.table(MDDpgrs_IM_MDD,file='table/MEDIATION.afterMR.MDDpgrs_IM_MDD.txt',quote = F,sep = '\t',row.names = F)
