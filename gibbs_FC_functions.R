tnorm <- function(n,lo,hi,mu,sig){   #generates truncated normal variates based on cumulative normal distribution
  #normal truncated lo and hi
  
  if(length(lo) == 1 & length(mu) > 1)lo <- rep(lo,length(mu))
  if(length(hi) == 1 & length(mu) > 1)hi <- rep(hi,length(mu))
  
  q1 <- pnorm(lo,mu,sig) #cumulative distribution
  q2 <- pnorm(hi,mu,sig) #cumulative distribution
  
  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  z[z == -Inf]  <- lo[z == -Inf]
  z[z == Inf]   <- hi[z == Inf]
  z
}
#------------------------------------
rmvnorm=function (n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), 
                  method = c("eigen", "svd", "chol")) 
{
  if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps), 
                   check.attributes = FALSE)) {
    stop("sigma must be a symmetric matrix")
  }
  if (length(mean) != nrow(sigma)) {
    stop("mean and sigma have non-conforming size")
  }
  sigma1 <- sigma
  dimnames(sigma1) <- NULL
  if (!isTRUE(all.equal(sigma1, t(sigma1)))) {
    warning("sigma is numerically not symmetric")
  }
  method <- match.arg(method)
  if (method == "eigen") {
    ev <- eigen(sigma, symmetric = TRUE)
    if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
      warning("sigma is numerically not positive definite")
    }
    retval <- ev$vectors %*% diag(sqrt(ev$values), length(ev$values)) %*% 
      t(ev$vectors)
  }
  else if (method == "svd") {
    sigsvd <- svd(sigma)
    if (!all(sigsvd$d >= -sqrt(.Machine$double.eps) * abs(sigsvd$d[1]))) {
      warning("sigma is numerically not positive definite")
    }
    retval <- t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
  }
  else if (method == "chol") {
    retval <- chol(sigma, pivot = TRUE)
    o <- order(attr(retval, "pivot"))
    retval <- retval[, o]
  }
  retval <- matrix(rnorm(n * ncol(sigma)), nrow = n) %*% retval
  retval <- sweep(retval, 2, mean, "+")
  colnames(retval) <- names(mean)
  retval
}
#------------------------------------
sample.betas=function(xtx,xmat,sig2,invT,z){
  prec=(1/sig2)*xtx+invT
  var1=solve(prec)
  pmedia=(1/sig2)*t(xmat)%*%z
  t(rmvnorm(1,var1%*%pmedia,var1))
}
#------------------------------------
sample.sig2=function(n,a.sig2,b.sig2,z,xmat,betas){
  media=xmat%*%betas
  sse=sum((z-media)^2)
  a1=(n/2)+a.sig2
  b1=b.sig2+(sse/2)
  1/rgamma(1,a1,b1)
}
#------------------------------------
sample.z=function(lo,hi,y,xmat,betas,sig2){
  media=xmat%*%betas
  z=y
  cond=y==hi
  z[cond]=tnorm(sum(cond),lo=hi,hi=Inf,mu=media[cond],sig=sqrt(sig2))
  cond=y==lo
  z[cond]=tnorm(sum(cond),lo=-Inf,hi=lo,mu=media[cond],sig=sqrt(sig2))
  z
}
#------------------------------------
summarize.param=function(betas,nburn,ngibbs,nomes){
  seq1=nburn:ngibbs
  betas1=apply(betas[seq1,],2,quantile,c(0.025,0.5,0.975))
  pgreater=apply(betas[seq1,]>0,2,mean)
  plower=apply(betas[seq1,]<0,2,mean)
  pval=apply(cbind(plower,pgreater),1,min)
  betas2=cbind(betas1[2,],betas1[1,],betas1[3,],pval)
  colnames(betas2)=c('Estimate','2.5% CI','97.5% CI','p-value')
  rownames(betas2)=nomes
  round(betas2,3)
}