rm(list=ls())
set.seed(1)

#import data
setwd('U:\\GIT_models\\FC_Model')
source('gibbs_FC_functions.R')
source('gibbs_FC.R')
dat=read.csv('fake data.csv',as.is=T)
n=nrow(dat)

#create covariate matrix for model
trat2=matrix(0,n,2)
trat2[,1]=ifelse(dat$trat=='Pred',1,0)
trat2[,2]=ifelse(dat$trat=='Descrip',1,0)
colnames(trat2)=c('trat1','trat2')
prescore=dat$prescore
interact=cbind(trat2[,1]*prescore,trat2[,2]*prescore)
colnames(interact)=c('int.trat1','int.trat2')
xmat=cbind(1,trat2,prescore,interact)
colnames(xmat)[1]='Intercept'

#other settings
nparam=ncol(xmat)
lo=0
hi=7
ngibbs=10000
var.betas=c(10,rep(1,nparam-1))
a.sig2=b.sig2=1

mod=gibbs_FC(xmat=xmat,y=dat$postscore,lo=lo,hi=hi,ngibbs=ngibbs,var.betas=var.betas,
             a.sig2=a.sig2,b.sig2=b.sig2)
nburn=ngibbs/2
betas.resumo=summarize.param(betas=mod$betas,nburn=nburn,ngibbs=ngibbs,nomes=colnames(xmat))