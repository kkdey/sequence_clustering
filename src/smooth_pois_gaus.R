mu=rep(0.1,1024)
mu[400:600]=1
y=rpois(1024,mu)

plot(mu,type='l')
lines(ashsmooth.pois(y),type='l',col=2)
lines(ashsmooth.gaus(y),col=4)

vest=ashsmooth.gaus(y,v.est=TRUE)
plot(,type='l')




##debug
tt=titable(y)$difftable
    v1=(y-lshift(y))^2/2
    v2=(rshift(y)-y)^2/2
    vv=(v1+v2)/2
vvv=titable(vv)$sumtable

vvv=titable(vest)$sumtable

ww=matrix(0,nrow=10,ncol=1024)
for(j in 0:9){
ind.nnull=(tt[j+2,]!=0|vvv[j+2,]!=0)
za=ash(tt[j+2,],sqrt(vvv[j+2,]),prior="nullbiased",multiseqoutput=TRUE,pointmass=TRUE,nullcheck=TRUE,VB=FALSE,mixsd=NULL,mixcompdist="normal",gridmult=2,lambda1=1,lambda2=0,df=NULL,trace=FALSE)
ww[j+1,ind.nnull]=za$PosteriorMean/2
ww[j+1,!ind.nnull]=0
print(j)
}


plot(reverse.gwave(sum(y),ww,-ww),type='l')





#' Interleave two vectors
#' @param x,y: two vectors of the same length
#' @return a vector of length twice that of x (or y)
interleave=function(x,y){
  return(as.vector(rbind(x,y)))
}


#' Shift a vector one unit to the right
#' @param x: a vector
#' @return a vector of the same length as that of x 
rshift = function(x){L=length(x); return(c(x[L],x[-L]))}
lshift = function(x){return(c(x[-1],x[1]))}


#' Produces two TI tables. One table contains the difference between adjacent pairs of data in the same resolution, and the other table contains the sum.
#' @param sig: a signal of length a power of 2
#' @return a list of two TI tables in the form of matrices
titable=function(sig){
  n = length(sig)
  J = log2(n)

  dmat = matrix(0, nrow=J+1, ncol=n)
  ddmat = matrix(0, nrow=J+1, ncol=n)
  dmat[1,] = sig
  ddmat[1,] = sig
  
  #dmat[1,] = as.matrix(sig)  
  #dmat2 = matrix(0, nrow=J, ncol=2*n) #the parent table
  
  for(D in 0:(J-1)){
    nD = 2^(J-D); 
    nDo2 = nD/2;
    twonD = 2*nD;
    for(l in 0:(2^D-1)){
      ind = (l*nD+1):((l+1)*nD)
      #ind2 = (l*twonD+1):((l+1)*twonD)
      x = dmat[D+1,ind]
      lsumx = x[seq(from=1,to=nD-1, by=2)] + x[seq(from=2,to=nD,by=2)]
      ldiffx = x[seq(from=1,to=nD-1, by=2)] - x[seq(from=2,to=nD,by=2)]
      rx = rshift(x);
      rsumx = rx[seq(from=1,to=nD-1, by=2)] + rx[seq(from=2,to=nD,by=2)]
      rdiffx = rx[seq(from=1,to=nD-1, by=2)] - rx[seq(from=2,to=nD,by=2)]
      dmat[D+2,ind] = c(lsumx,rsumx)
      ddmat[D+2,ind] = c(ldiffx,rdiffx)
    }
  }
  return(list(sumtable=dmat,difftable=ddmat))
}


#' Turns a TItable into a matrix that one can apply glm to (for Poisson denoising).
#' @param dmat: a k by n matrix
#' return a 2 by nk/2 matrix
simplify=function(dmat){
  matrix(t(dmat),nrow=2)
}

#' Reverse wavelet transform a set of wavelet coefficients in TItable format for Gaussian data.
#' @param lp: a J by n matrix of estimated wavelet coefficients. 
#' @param lq: a J by n matrix of complementary wavelet coefficients.
#' @param est: an n-vector. Usually a constant vector with each element equal to the estimated total mean.
#' @return reconstructed signal in the original data space.
reverse.gwave=function(est,lp,lq=NULL){
  if(is.null(lq)){
    lq = -lp
  }
  if(length(est)==1){est = rep(est,ncol(lp))}
  
  J=nrow(lp)
  
  for(D in J:1){
    #print(exp(est))
    #readline("press a key")
    nD = 2^(J-D+1) 
    nDo2 = nD/2
    for(l in 0:(2^(D-1)-1)){

      ind = (l*nD+1):((l+1)*nD)
      
      estvec = est[ind]/2
      lpvec = lp[D,ind]
      lqvec = lq[D,ind]

      estl = estvec[1:nDo2]
      lpl=lpvec[1:nDo2]
      lql=lqvec[1:nDo2]
      nestl = interleave(estl+lpl,estl+lql) #interleaves the two
      

      estr = estvec[(nDo2+1):nD]
      lpr = lpvec[(nDo2+1):nD]
      lqr = lqvec[(nDo2+1):nD]
      nestr = interleave(estr+lpr,estr+lqr) #interleaves the two
      nestr = lshift(nestr)

      est[ind] = 0.5*( nestl + nestr )
    }
  }
  return(est)
}


