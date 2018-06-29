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
set.seed(333)
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
# write.csv(coefDf,file="coefDf.csv")

# Coefficient forest plot:
# png(file="brm1Coef2.png",height=4.5,width=8.5,units="in",res=400)
ggplot(data=coefDf,aes(x=Metabolite,y=Mean,ymin=`25%`,ymax=`75%`,color=Group))+
  geom_pointrange()+geom_hline(yintercept=0,lty=2)+
  geom_errorbar(aes(ymin=`25%`,ymax=`75%`),width=0.5)+
  coord_flip()+theme_bw()+ylab("Estimate")+
  scale_color_manual(values=c(rgb(255,51,51,max=255,alpha=125),
        rgb(0,0,153,max=255,alpha=125)))
# dev.off()

# MCMC chain:
exMCMCChain<-rbind(data.frame(Value=brm1$fit@sim$samples[[1]]$b_muType1MI_M110,
                 Iteration=1:10000,Group="Thrombotic MI"),
      data.frame(Value=brm1$fit@sim$samples[[1]]$b_muType2MI_M110,
                 Iteration=1:10000,Group="Non-Thrombotic MI"))
exMCMCChain$Group<-factor(exMCMCChain$Group,levels=c("Thrombotic MI","Non-Thrombotic MI"))

# png(file="exMCMCChain.png",height=4,width=6,units="in",res=300)
ggplot(exMCMCChain,aes(x=Iteration,y=Value,color=Group))+
  geom_path()+geom_hline(yintercept=0,lty=2)+
  scale_color_manual(values=c(rgb(255,51,51,max=255,alpha=150),
                                          rgb(0,0,153,max=255,alpha=150)))+
  theme_bw()+xlim(5000,7000)+ylim(-2.5,3.5)
# dev.off()

# png(file="exHist.png",height=4,width=4,units="in",res=300)
ggplot(exMCMCChain,aes(x=Value,fill=Group))+geom_density(alpha=.3)+
  theme_bw()+ylab("Density")+
  scale_fill_manual(values=c(rgb(255,51,51,max=255,alpha=150),
                              rgb(0,0,153,max=255,alpha=150)))+
  coord_flip()
# dev.off()

############ Frequentist model ############
freqModel<-nnet::multinom(group~.,
                          data=metab2[,names(metab2)!="ptid"])
summary(freqModel)

############ Posterior predictive distribution ############
brm1Coefs<-brm1$fit@sim$samples[[1]]
brm1Coefs<-brm1Coefs[names(brm1Coefs)!="lp__"]

# Get coefficients from simulated posterior
brm1Coefs2<-data.frame()
for(i in 1:length(brm1Coefs)){
  x<-brm1Coefs[i]
  brm1Coefs2<-rbind(brm1Coefs2,data.frame(Parameter=names(x),
         iteration=1:length(x[[1]]),Estimate=x[[1]]))
}

# Posterior predictive distribution for 2010
brm1Preds<-list()
for(i in 1:10000){
  brm1Samp1<-metab2[6,]
  brm1Coefs3<-brm1Coefs2[brm1Coefs2$iteration==i,]
  brm1Int1<-matrix(1,ncol=2)
  colnames(brm1Int1)<-c("b_muType1MI_Intercept","b_muType2MI_Intercept")
  brm1Samp1a<-brm1Samp1b<-brm1Samp1[!names(brm1Samp1)%in%c("ptid","group")]
  names(brm1Samp1a)<-paste0("b_muType1MI_",names(brm1Samp1a))
  names(brm1Samp1b)<-paste0("b_muType2MI_",names(brm1Samp1b))
  brm1Samp1<-cbind(brm1Int1,brm1Samp1a,brm1Samp1b)
  brm1Coefs3$x<-t(brm1Samp1[,match(brm1Coefs3$Parameter,names(brm1Samp1))])
  brm1Coefs3$prod<-c(brm1Coefs3$x)*c(brm1Coefs3$Estimate)
  brm1Coefs3$group<-c("Type1","Type2")[(!grepl("Type1",brm1Coefs3$Parameter))+1]
  brm1Coefs3<-brm1Coefs3 %>% group_by(group) %>% 
    summarise(sum=sum(prod)) %>% as.data.frame()
  brm1Coefs3<-rbind(brm1Coefs3,data.frame(group="sCAD",sum=1))
  brm1Pred<-data.frame(t(exp(brm1Coefs3$sum)/sum(exp(brm1Coefs3$sum))))
  names(brm1Pred)<-c("Type1","Type2","sCAD")
  brm1Preds[[i]]<-brm1Pred
  print(i)
}
brm1Preds<-do.call("rbind",brm1Preds)

# Reshape:
brm1Preds$Iteration<-1:nrow(brm1Preds)
brm1Preds<-brm1Preds %>% gather(key="Group",value="Probability",-Iteration)
brm1Preds$Group<-factor(brm1Preds$Group)
levels(brm1Preds$Group)<-c("sCAD","Thrombotic MI","Non-Thrombotic MI")
brm1Preds$Group<-factor(brm1Preds$Group,
        levels=levels(brm1Preds$Group)[c(2,3,1)])

# PTID 2010 MCMC
# png(file="ptid2010MCMC.png",height=3.75,width=7,units="in",res=300)
ggplot(brm1Preds,aes(x=Iteration,y=Probability,color=Group))+
  geom_path()+xlim(5000,5250)+theme_bw()+
  scale_color_manual(values=c(rgb(255,51,51,max=255,alpha=170),
                              rgb(0,0,153,max=255,alpha=170),
                              rgb(50,190,0,max=255,alpha=200)))
# dev.off()

# PTID 2010 Histogram:
# png(file="ptid2010Hist.png",height=3,width=7.5,units="in",res=300)
ggplot(brm1Preds,aes(x=Probability,fill=Group,y=..density..))+
  geom_histogram(bins=45)+facet_wrap(~Group,scales="free_y")+theme_bw()+
  scale_fill_manual(values=c(rgb(255,51,51,max=255,alpha=255),
                              rgb(0,0,153,max=255,alpha=255),
                              rgb(50,190,0,max=255,alpha=255)))+
  ylab("Density")
# dev.off()

############ Cross-validation error ############
set.seed(3)
cvF<-cvTools::cvFolds(n=nrow(metab2),K=nrow(metab2))
cvF<-data.frame(fold=cvF$which,id=cvF$subsets)

phenoFoldsFreq<-phenoFolds<-data.frame()
for(k in 1:nrow(metab2)){
  brmFold<-brm(group~.,data=metab2[cvF$id[cvF$fold!=k],names(metab2)!="ptid"],
               family="categorical",chains=4,iter=5000,algorithm="sampling",
               prior=priors,seed=k+3)
  brmFoldFreq<-nnet::multinom(group~.,
                data=metab2[cvF$id[cvF$fold!=k],names(metab2)!="ptid"])
  predBrmFold<-predict(brmFold,
            newdata=metab2[cvF$id[cvF$fold==k],!names(metab2)%in%c("group","ptid")])
  predBrmFoldFreq<-predict(brmFoldFreq,
            newdata=metab2[cvF$id[cvF$fold==k],!names(metab2)%in%c("group","ptid")],
            type="probs")
  phenoFold<-cbind(pheno[cvF$id[cvF$fold==k],],predBrmFold)
  phenoFoldFreq<-cbind(pheno[cvF$id[cvF$fold==k],],as.data.frame(t(predBrmFoldFreq)))
  phenoFolds<-rbind(phenoFolds,phenoFold)
  phenoFoldsFreq<-rbind(phenoFoldsFreq,phenoFoldFreq)
  print(k)
}

names(phenoFolds)[3:5]<-c("sCAD","Type 1 MI","Type 2 MI")
names(phenoFoldsFreq)[3:5]<-c("sCAD","Type 1 MI","Type 2 MI")
phenoFolds$Predicted<-names(phenoFolds)[3:5][apply(phenoFolds[,3:5],1,which.max)]
phenoFoldsFreq$Predicted<-names(phenoFoldsFreq)[3:5][apply(phenoFoldsFreq[,3:5],1,which.max)]
xtabs(~group+Predicted,data=phenoFolds)
xtabs(~group+Predicted,data=phenoFoldsFreq)

########### Add troponin in: ###########
oxPL<-read.csv("~/gdrive/Athro/oxPL6/wide_data_20150529.csv")
trop<-oxPL %>% select(ptid,tropT0)
trop$ptid<-as.character(trop$ptid)
metab3<-metab2 %>% left_join(trop)

ptm<-proc.time()
set.seed(333)
brm2<-brm(group~.,data=metab3[,names(metab3)!="ptid"],
          family="categorical",chains=4,iter=10000,algorithm="sampling",
          prior=priors,seed=33)
proc.time()-ptm

# Summary and predicted probabilities
summary(brm2)
predBrm2<-predict(brm2,newdata=metab3[,!names(metab3)%in%c("group","ptid")])
predBrm2<-as.data.frame(predBrm2)
pheno3<-cbind(pheno,predBrm2)

# Shiny stan
launch_shinystan(brm2)

# Coefficients:
coefDf50_2<-as.data.frame(posterior_interval(brm2,prob=.50))
coefDf95_2<-as.data.frame(posterior_interval(brm2,prob=.95))
coefDfMean_2<-as.data.frame(posterior_summary(brm2))
coefDfMean_2<-coefDfMean_2 %>% select(Mean=Estimate)
coefDf_2<-cbind(coefDfMean_2,coefDf50_2,coefDf95_2)
coefDf_2$Parameter<-rownames(coefDf_2)
coefDf_2$Metabolite<-str_split(coefDf_2$Parameter,"_",simplify=TRUE)[,3]
coefDf_2$Group<-gsub("mu","",str_split(coefDf_2$Parameter,"_",simplify=TRUE)[,2])
coefDf_2<-coefDf_2 %>% filter(Metabolite!="") 
coefDf_2$Group<-factor(coefDf_2$Group)
levels(coefDf_2$Group)<-c("Thrombotic MI","Non-Thrombotic MI")
coefDf_2<-coefDf_2 %>% filter(Metabolite!="Intercept")
coefDf_2$Metabolite<-key$biochemical[match(coefDf_2$Metabolite,key$id)]
# write.csv(coefDf,file="coefDf_2.csv")

# Coefficient forest plot:
coefDf_2$Metabolite[is.na(coefDf_2$Metabolite)]<-"Troponin"
# png(file="brm2Coef.png",height=4.5,width=8.5,units="in",res=400)
ggplot(data=coefDf_2,aes(x=Metabolite,y=Mean,ymin=`25%`,ymax=`75%`,color=Group))+
  geom_pointrange()+geom_hline(yintercept=0,lty=2)+
  geom_errorbar(aes(ymin=`25%`,ymax=`75%`),width=0.5)+
  coord_flip()+theme_bw()+ylab("Estimate")+
  scale_color_manual(values=c(rgb(255,51,51,max=255,alpha=125),
                              rgb(0,0,153,max=255,alpha=125)))
# dev.off()

# Model comparison
compare_ic(WAIC(brm1),WAIC(brm2))

############ Cross-validation error w/ Troponin ############
set.seed(3)
cvF<-cvTools::cvFolds(n=nrow(metab3),K=nrow(metab3))
cvF<-data.frame(fold=cvF$which,id=cvF$subsets)

phenoFolds2<-data.frame()
for(k in 1:nrow(metab3)){
  brmFold<-brm(group~.,data=metab3[cvF$id[cvF$fold!=k],names(metab3)!="ptid"],
               family="categorical",chains=4,iter=5000,algorithm="sampling",
               prior=priors,seed=k+3)
  print(k)
  predBrmFold<-predict(brmFold,
                       newdata=metab3[cvF$id[cvF$fold==k],
                                      !names(metab3)%in%c("group","ptid")])
  phenoFold2<-cbind(pheno[cvF$id[cvF$fold==k],],predBrmFold)
  phenoFolds2<-rbind(phenoFolds2,phenoFold2)
}

names(phenoFolds2)[3:5]<-c("sCAD","Type 1 MI","Type 2 MI")
phenoFolds2$Predicted<-names(phenoFolds2)[3:5][apply(phenoFolds2[,3:5],1,which.max)]
xtabs(~group+Predicted,data=phenoFolds2)

############ Combined probability estimates ############
names(phenoFolds)[names(phenoFolds) %in% c("sCAD","Type 1 MI","Type 2 MI")]<-
  paste0(names(phenoFolds)[names(phenoFolds) %in% c("sCAD","Type 1 MI","Type 2 MI")],"_1")

names(phenoFolds2)[names(phenoFolds2) %in% c("sCAD","Type 1 MI","Type 2 MI")]<-
  paste0(names(phenoFolds2)[names(phenoFolds2) %in% c("sCAD","Type 1 MI","Type 2 MI")],"_2")

phenoFolds3<-phenoFolds %>% left_join(phenoFolds2,by=c("ptid"))
# write.csv(phenoFolds3,file="phenoFolds.csv")

############ Coefficient correlations w/ Troponin ############
brm2Coefs<-brm2$fit@sim$samples[[1]]
#brm2Coefs<-brm2Coefs[names(brm2Coefs)!="lp__"]

# Get coefficients from simulated posterior
brm2Coefs2<-data.frame()
for(i in 1:length(brm2Coefs)){
  x<-brm2Coefs[i]
  brm2Coefs2<-rbind(brm2Coefs2,data.frame(Parameter=names(x),
                                          iteration=1:length(x[[1]]),Estimate=x[[1]]))
}
brm2Coefs3<-brm2Coefs2 %>% spread(key="iteration",value="Estimate")
tempParams<-brm2Coefs3[,1]
brm2Coefs3<-t(brm2Coefs3[,-1])
colnames(brm2Coefs3)<-tempParams
# png(file="corplot.png",height=6,width=6,units="in",res=300)
corrplot::corrplot(cor(brm2Coefs3))
# dev.off()

plot(brm2Coefs3[,c("b_muType1MI_M25","b_muType1MI_M72")])
plot(brm2Coefs3[,c("b_muType2MI_M25","b_muType2MI_M72")])

brm2Coefs3<-brm2Coefs3[brm2Coefs3[,'lp__']>-100,]
# png(file="coefPost.png",height=6,width=6,units="in",res=300)
plot3D::scatter3D(brm2Coefs3[,'b_muType1MI_M25'],brm2Coefs3[,'b_muType1MI_M72'],
                  brm2Coefs3[,'lp__'],zlim=c(-75,-40),phi=20,theta=30,
                  xlab="1-linoleoylglycerol (18:2)",ylab="2-linoleoylglycerol (18:2)",
                  zlab="Log posterior")
# dev.off()

brm2Coefs4<-brm2Coefs2 %>% 
  filter(Parameter %in% c('iteration','b_muType1MI_M25','b_muType1MI_M72','lp__'))
brm2Coefs4<-brm2Coefs4 %>% spread(key='Parameter',value='Estimate')
names(brm2Coefs4)<-c("Iteration","1-linoleoylglycerol (18:2)",
                     "2-linoleoylglycerol (18:2)","Log posterior")

# png(file="coefPost2.png",height=4,width=6,units="in",res=300)
ggplot(brm2Coefs4 %>% filter(`Log posterior`>-100) %>% sample_frac(.75),
       aes(x=`1-linoleoylglycerol (18:2)`,y=`2-linoleoylglycerol (18:2)`,
                      color=`Log posterior`))+
  geom_point(size=.5)+theme_bw()+scale_color_gradient2(high="red",mid="purple",
                                      low="blue",midpoint=-60)
# dev.off()

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
