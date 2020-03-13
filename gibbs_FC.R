gibbs_FC=function(xmat,y,lo,hi,ngibbs,var.betas,a.sig2,b.sig2){
  xmat=data.matrix(xmat)
  xtx=t(xmat)%*%xmat
  
  #initial values
  nparam=ncol(xmat)
  betas=rep(0,nparam)
  sig2=1
  z=y
  n=length(y)
  
  #priors
  invT=diag(1/var.betas)
  
  #MCMC stuff
  store.betas=matrix(NA,ngibbs,nparam)
  store.sig2=matrix(NA,ngibbs,1)
  
  for (i in 1:ngibbs){
    print(i)
    betas=sample.betas(xtx=xtx,xmat=xmat,sig2=sig2,invT=invT,z=z)
    sig2=sample.sig2(n=n,a.sig2=a.sig2,b.sig2=b.sig2,z=z,xmat=xmat,betas=betas)
    z=sample.z(lo=lo,hi=hi,y=y,xmat=xmat,betas=betas,sig2=sig2)
    
    #store results
    store.betas[i,]=betas
    store.sig2[i]=sig2
  }
  list(betas=store.betas,sig2=store.sig2)
}
