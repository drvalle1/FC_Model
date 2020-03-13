rm(list=ls())
set.seed(3)

#basic settings
n=90
lo=0
hi=7

#create covariates
trat=sample(c('Pred','Descrip','Control'),size=n,replace=T)
trat2=matrix(0,n,2)
trat2[,1]=ifelse(trat=='Pred',1,0)
trat2[,2]=ifelse(trat=='Descrip',1,0)
colnames(trat2)=c('trat1','trat2')
prescore=runif(n,min=lo,max=hi)
interact=cbind(trat2[,1]*prescore,trat2[,2]*prescore)
colnames(interact)=c('inter1','inter2')
xmat=cbind(1,trat2,prescore,interact)

#parameters
betas.true=betas=c(1,-1,-0.5,0.8,0.7,-0.8)
sig2.true=sig2=0.4^2
z=rnorm(n,mean=xmat%*%betas,sd=sqrt(sig2))
y=z; y[z>hi]=hi; y[z<lo]=lo
plot(prescore,y,pch=19,col='grey',ylim=range(z))
points(prescore,z)

#create final dataset
postscore=y
fim=cbind(postscore,trat,prescore)
setwd('U:\\GIT_models\\FC_Model')
write.csv(fim,'fake data.csv',row.names=F)