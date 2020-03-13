nburn=ngibbs/2
betas.estim=apply(store.betas[nburn:ngibbs,],2,quantile,c(0.025,0.5,0.975))

rango=range(rbind(betas.estim,betas.true))
plot(betas.true,betas.estim[2,],ylim=rango,xlim=rango)
lines(rango,rango)
for (i in 1:nparam){
  x=rep(betas.true[i],2)
  y=betas.estim[-2,i]
  lines(x,y,col='grey')
}

plot(store.sig2,type='l')
abline(h=sig2.true,col='red')

