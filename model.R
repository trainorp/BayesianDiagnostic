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
include<-key %>% filter(weights>.035)
metab2<-metab[,names(metab) %in% include$id]
rm(outs)

############ Model ############
setwd("~/gdrive/Dissertation/Aim3")
metab2$ptid<-rownames(metab2)
metab2<-metab2 %>% left_join(pheno)

ptm<-proc.time()
priors<-c(prior(normal(0,1),class=b),
          prior(normal(0,4),class=Intercept))
brm1<-brm(group~.,data=metab2[,names(metab2)!="ptid"],
          family="categorical",chains=4,iter=10000,algorithm="sampling",
          prior=priors,seed=33)
proc.time()-ptm

# Summary and predicted probabilities
summary(brm1)
predBrm1<-predict(brm1,newdata=metab2[,!names(metab2)%in%c("group","ptid")])
predBrm1<-as.data.frame(predBrm1)
pheno2<-cbind(pheno,predBrm1)

# Shiny stan
launch_shinystan(brm1)

# Coefficients:
coefDf50<-as.data.frame(posterior_interval(brm1,prob=.50))
coefDf95<-as.data.frame(posterior_interval(brm1,prob=.95))
coefDfMean<-as.data.frame(posterior_summary(brm1))
coefDfMean<-coefDfMean %>% select(Mean=Estimate)
coefDf<-cbind(coefDfMean,coefDf50,coefDf95)
coefDf$Parameter<-rownames(coefDf)
coefDf$Metabolite<-str_split(coefDf$Parameter,"_",simplify=TRUE)[,3]
coefDf$Group<-gsub("mu","",str_split(coefDf$Parameter,"_",simplify=TRUE)[,2])
coefDf<-coefDf %>% filter(Metabolite!="") 
coefDf$Group<-factor(coefDf$Group)
levels(coefDf$Group)<-c("Thrombotic MI","Non-Thrombotic MI")
coefDf<-coefDf %>% filter(Metabolite!="Intercept")
coefDf$Metabolite<-key$biochemical[match(coefDf$Metabolite,key$id)]

# Colors
png(file="brm1Coef.png",height=4.5,width=8.5,units="in",res=400)
ggplot(data=coefDf,aes(x=Metabolite,y=Mean,ymin=`25%`,ymax=`75%`,color=Group))+
  geom_pointrange()+geom_hline(yintercept=0,lty=2)+
  geom_errorbar(aes(ymin=`25%`,ymax=`75%`),width=0.5)+
  coord_flip()+theme_bw()+ylab("Estimate")+
  scale_color_manual(values=c(rgb(255,51,51,max=255,alpha=125),
        rgb(0,0,153,max=255,alpha=125)))
dev.off()

############ Cross-validation error ############
set.seed(3)
cvF<-cvTools::cvFolds(n=nrow(metab2),K=nrow(metab2))
cvF<-data.frame(fold=cvF$which,id=cvF$subsets)

phenoFolds<-data.frame()
for(k in 1:nrow(metab2)){
  brmFold<-brm(group~.,data=metab2[cvF$id[cvF$fold!=k],names(metab2)!="ptid"],
               family="categorical",chains=2,iter=5000,algorithm="sampling",
               prior=priors,seed=k+3)
  predBrmFold<-predict(brmFold,
            newdata=metab2[cvF$id[cvF$fold==k],!names(metab2)%in%c("group","ptid")])
  phenoFold<-cbind(pheno[cvF$id[cvF$fold==k],],predBrmFold)
  phenoFolds<-rbind(phenoFolds,phenoFold)
}

########### Horseshoe prior LDA ###########
metab3<-model.matrix(~group,data=metab2)[,2:3]
colnames(metab3)<-c("Type1","Type2")
metab3<-cbind(metab2[,!names(metab2)%in%c("ptid","group")],metab3)
prior3<-set_prior(horseshoe(1))
ptm<-proc.time()
brm3<-brm(cbind(Type1,Type2)~.,data=metab3,chains=4,iter=10000,algorithm="sampling",
          prior=prior3,control=list(adapt_delta=.9))
proc.time()-ptm
launch_shinystan(brm3)
summary(brm3)
predict(brm3)
