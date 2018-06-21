############ Prereqs ############
oldPar<-par()
options(stringsAsFactors=FALSE)
library(tidyverse)
library(rstan)
library(brms)

############ Import annotation and metabolite data ############
setwd("~/gdrive/AthroMetab/Data")
spec<-read.csv("AthroACSRawSpectra.csv")
key<-read.csv("metabolite_key2.csv")
key[grepl("X - ",key$biochemical),]$biochemical<-key[grepl("X - ",key$biochemical),]$Unknown.Name

key<-spec %>% left_join(key,c("comp_id"="MEBID"))

getMZ<-function(x){
  x<-strsplit(x," ")[[1]]
  return(unlist(strsplit(x[grepl(":100",x)],":"))[1])
}
key$peakMZ<-unlist(lapply(lapply(key$spectra,getMZ),function(x) ifelse(is.null(x),NA,x)))

metab<-read.csv("scaled.csv")
metab$ptid<-as.character(metab$ptid)
pheno<-metab %>% filter(timepoint=="T0") %>% select(ptid,group)
metab<-metab %>% filter(timepoint=="T0") %>% select(-ptid,-group,-timepoint)
metab<-metab %>% log2()
rownames(metab)<-pheno$ptid

############ Load feature selection data ############
load("~/gdrive/AthroMetab/WoAC/outs.RData")
weights<-c()
for(i in 1:length(outs)){
  nams<-outs[[i]]$vars
  experts<-outs[[i]]$pop@varInclude[order(outs[[i]]$pop@cost)[1]]
  experts<-sapply(experts,as.numeric)
  rownames(experts)<-nams
  blanks<-matrix(0,nrow=ncol(metab),ncol=1)
  rownames(blanks)<-colnames(metab)
  blanks[match(rownames(experts),rownames(blanks)),]<-experts
  weights<-cbind(weights,apply(blanks,1,mean))
}
key<-key %>% arrange(as.numeric(gsub("M","",id)))
weights<-apply(weights,1,function(x) mean(x))
key$weights<-0
key$weights[match(names(weights),key$id)]<-weights
include<-key %>% filter(weights>.025)
metab2<-metab[,names(metab) %in% include$id]

############ Model ############
metab2$ptid<-rownames(metab2)
metab2<-metab2 %>% left_join(pheno)

priors<-get_prior(group~.,data=metab2[,names(metab2)!="ptid"],family="categorical")

ptm<-proc.time()
brm1<-brm(group~.,data=metab2[,names(metab2)!="ptid"],
        family="categorical",chains=2,iter=5000,algorithm="sampling")
proc.time()-ptm
launch_shinystan(brm1)

ptm<-proc.time()
prior2<-c(prior(normal(0,10),class=b),
          prior(normal(0,1),class=b,coef=muType2MI_M133),
          prior(normal(1,2),class=Intercept))
priors$prior[5]<-"exponential(1)"
brm2<-brm(group~.,data=metab2[,names(metab2)!="ptid"],
          family="categorical",chains=2,iter=500,algorithm="sampling",
          prior=priors)
proc.time()-ptm
launch_shinystan(brm1)
