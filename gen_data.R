
##################### the correaltion structure of gene ####################
CorMat = function(m, sig.index = "En"){
  sig.G= matrix(0,m,m)
  if (sig.index != "En"){
      for (i in 1:m){
    for(j in i:m){
      if(sig.index=="AR1"){
        sig.G[i,j]=0.2^abs(i-j)
      }else if(sig.index=="AR2"){
        sig.G[i,j]=0.8^abs(i-j)
      }else if(sig.index=="B2"){
        if(abs(i-j) ==1){
          sig.G[i,j]= 0.6
        }else if (abs(i-j)==2){
          sig.G[i,j] = 0.33
        }
      }else if(sig.index=="B1"){
        if(abs(i-j) ==1){
          sig.G[i,j]= 0.33
        }
      }else if(sig.index=="AR5"){
        sig.G[i,j]=0.5^abs(i-j)
      }
      sig.G[j,i] = sig.G[i,j]
    }
  }
}
  diag(sig.G)= 1
  return(sig.G)
}


#############################generate the data ############################
#@n the number of subjects
#@p, the dim of gene variable
#@q. the dim of envrionmental variable
#@ Gcorr, the Gene correlation
#@ Ecorr, the envrionment correlation
#@ err.dis, the error distribution of the residual
#@ a the envrionment parameter, q*1
#@ b the  gene and gene-envrionment parameter, p*(q+1)
#@family, the model, linear, logistic
#@is.SNP, boolean, if True, the gene varibale is a SNP type
#@maf, float, maf \in (0,0.5)
#@E.cat, int, the dis number of envrionment varibles
GenGE <- function(n, p, q, Gcorr="En", Ecorr="En", err.dis,
               alpha0=0, a, b, family="gaussian", is.SNP = FALSE, maf= 0.3, E.cat = 0){
  Zsigmat <- CorMat(p, Gcorr)
  Xsigmat <- CorMat(q, Ecorr)
  Z <- mvtnorm::rmvnorm(n, mean=rep(0,p), sigma=Zsigmat) 
  X <- mvtnorm::rmvnorm(n, mean=rep(0,q), sigma=Xsigmat)

  if (is.SNP){
    MAF.vec <- rep(maf, p)
    # fre.mat <- rbind((1-MAF.vec)^2, 2*MAF.vec*(1-MAF.vec), MAF.vec^2) 
    cum.fre <- rbind((1-MAF.vec)^2, 1-MAF.vec^2)
    cutPoint <- qnorm(cum.fre)
    # Z = sapply(1:p,  function(j){cut(x = Z[,j],cutPoint[,j])})
 
    Z <- sapply(1:p,  function(j){
      sapply(Z[,j], function(x){ sum(x > cutPoint[,j])})
    })

  }

  if (E.cat > 0){
      X[, 1:E.cat] <- as.numeric(X[,1:E.cat] > 0 )
  }
  
  W <- matrix(NA, nrow=n, ncol=p*(q+1))
  X1 <- cbind(1,X)
  for (j in 1:p){
    W[,(1:(q+1))+(q+1)*(j-1)] <- Z[,j]*X1
  }

  ymean = alpha0+ X%*%a + W%*%b

  ### continue response
  if (family == "gaussian"){
     ##error distribution 
      ep = rep(0,n)
      if (err.dis == "M_NC"){
        ep =  ifelse(runif(n,0,1)< 0.8,rnorm(n),rcauchy(n))
      }else if(err.dis == "M_NL"){
        ep  = ifelse(runif(n,0,1)< 0.8,rnorm(n),exp(rnorm(n)))
      }else if(err.dis== "t"){
        ep=rt(n,2)
      }else if(err.dis == "norm"){
        ep=rnorm(n)
      }
      y = ymean+ep
      data = list(Y=y,Z=Z, X=X, W=W, err.dis =ep)
  }else if(family == "binomial"){
      prob = 1/(1+exp(-ymean))
      y <- rbinom(n,1,prob)
      data = list(Y=y,Z=Z, X=X, W=W)
    }  
  return(data)
}

