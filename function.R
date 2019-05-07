###Sig : p value based method
Sig <- function(y, X, W, family, censor = F, fdr.level = .1){
  n <- length(y)
  q <- ncol(X)
  p <- ncol(W)/(q+1)
  pmat <- matrix(NA,p,(q+1))
  for (j in 1:p){
    Wj <-  W[,(1:(q+1))+(q+1)*(j-1)]
    if (censor){
      fit <- summary(survreg(Surv(exp(y), delta) ~ X + Wj))
      pmat[j,] <- fit$table[(q+2):(2*q+2),4]
    }else{
      fit <- summary(glm(y ~ X + Wj, family=family))
      pmat[j,] <- fit$coefficients[-(1:(q+1)),4]
    }
  }
  # G
  pval <- pmat[,1]
  qval = p.adjust(pval, method ="fdr")
  s <- qval < fdr.level
  #GE
  pvalGE <- as.vector(pmat[,-1])
  qvalGE <- p.adjust(pvalGE, method="fdr")
  sGE <- qvalGE < fdr.level
  return(list(pmat = pmat, G =list(p=pval,q=qval,s=s), GE=list(p= pvalGE, q=qvalGE, s=sGE)))
}


### marginal penalization with MCP
mMCP <- function(y, X, W, lambda, family = "gaussian", penalty = "MCP", 
                 weight =NULL){
  q <- ncol(X)
  p <- ncol(W)/(q+1)
  nlambda <- length(lambda)
  if (!is.null(weight)){
    y=y[weight!=0]
    X=X[weight!=0,]
    W=W[weight!=0,]
    weight=weight[weight!=0]
    n=length(y)
    weight=weight/sum(weight)*n
    
    ymean=sum(weight*y)/sum(weight)
    Xmean=colSums(X*(weight%*%matrix(1,1,q)))/sum(weight)
    Wmean=colSums(W*(weight%*%matrix(1,1,p*(q+1))))/sum(weight)
    
    y=sqrt(weight)*(y-ymean)
    X = (sqrt(weight)%*%matrix(1,1,q))*(X-matrix(1,n,1)%*%Xmean)
    W = (sqrt(weight)%*%matrix(1,1,p*(q+1)))*(W-matrix(1,n,1)%*%Wmean)
  }
  pen.coef  <- array(NA, dim=c(p,2*q+2,nlambda))
  for (j in 1:p) {
    Wj <- W[,(1:(q+1))+(q+1)*(j-1)]
    Dj <- cbind(X, Wj)
    # ###penalization
    pfit = ncvreg::ncvreg(X = Dj, y=y,family = family,lambda= lambda,penalty=penalty,
                          penalty.factor = c(rep(0,q),rep(1,q+1)))
    nl <- ncol(pfit$beta)
    pen.coef[j,,1:nl] <- pfit$beta
    if (nl != nlambda){
      pen.coef[j, ,(nl+1): nlambda] <- pfit$beta[,nl]
    }
  }
  alpha0 <- pen.coef[,1,]
  a <- pen.coef[,(2:(q+1)),]
  beta <- pen.coef[,q+2,]
  gamma <- pen.coef[,-(1:(q+2)),]
  return(list(coef=pen.coef,alpha0 = alpha0, a = a , beta = beta, gamma = gamma, lambda = lambda) )
}

####sgMCP :: the proposed method 
#@y response vector. a vector,  n*1
#@X envrionment variable,  matrix: n*q
#@W gene and gene-envrionment variable, matrix: n*(p*(q+1)). W = [Z_1, Z_1X_1, ..., Z_1X_q, Z_2X_1, ..., Z_pX_q ]
msgMCP <- function(y, X, W, Matla, family = "gaussian",weight =NULL){
  q <- ncol(X)
  p <- ncol(W)/(q+1)
  if (!is.null(weight)){
    y=y[weight!=0]
    X=X[weight!=0,]
    W=W[weight!=0,]
    weight=weight[weight!=0]
    n=length(y)
    weight=weight/sum(weight)*n
    ymean=sum(weight*y)/sum(weight)
    Xmean=colSums(X*(weight%*%matrix(1,1,q)))/sum(weight)
    Wmean=colSums(W*(weight%*%matrix(1,1,p*(q+1))))/sum(weight)
    
    y=sqrt(weight)*(y-ymean)
    X = (sqrt(weight)%*%matrix(1,1,q))*(X-matrix(1,n,1)%*%Xmean)
    W = (sqrt(weight)%*%matrix(1,1,p*(q+1)))*(W-matrix(1,n,1)%*%Wmean)
  }
  L <- nrow(Matla)
  pen.coef  <- array(NA, dim=c(p,2*q+2, L))
  for (j in 1:p) {
    Wj <- W[,(1:(q+1))+(q+1)*(j-1)]
    Dj <- cbind(X, Wj)
    # ###penalization
    pen.coef[j,,] = msgMCP.e(y, X, Wj, xi = 3 , max_iter = 100, family, Matla, eps = 0.001)$beta
  }
  alpha0 <- pen.coef[,1,]
  a <- pen.coef[,(2:(q+1)),]
  beta <- pen.coef[,q+2,]
  gamma <- pen.coef[,-(1:(q+2)),]
  return(list(coef=pen.coef,alpha0 = alpha0, a = a , beta = beta, gamma = gamma, Matla= Matla) )
}

###For each j function
msgMCP.e  <- function(y, X, Wj, xi, max_iter = 100, family, Matla, eps = 0.001){
  
  q <- ncol(X)
  ## Set up XX, yy, lambda
  XW <- scale(cbind(X,Wj), center= T, scale = T)
  XX <- XW[,1:q]
  Wjj <- XW[,-(1:q)]
  if(family == "gaussian"){
    yy <- y - mean(y)
  }else if (family == "binomial"){
    yy <- y
  }
  #
  n <- length(yy)
  L <- nrow(Matla)
  
  if (family == "gaussian"){
    res = gaussianetscpp(yy, XX, Wjj, xi, max_iter, Matla, eps, iset = q:(2*q))
    # res = gaussian.estimate(yy, XX, Wjj, xi, max_iter, Matla, eps)
    a <- rep(mean(y),L)
  }else if (family == "binomial"){
    res = binomialetscpp(yy, XX, Wjj, xi, max_iter, Matla, eps, iset = (q+1):(2*q+1))
    a <- res$beta0
  }
  betahat <- res$beta
  iter <- res$iter
  ## Unstandardize
  beta <- matrix(0, nrow=(2*q+2), ncol=L)
  ns = 1:(2*q+1)
  bb <- betahat/attr(XW, "scaled:scale")
  beta[ns+1,] <- bb
  beta[1,] <- a - crossprod(attr(XW, "scaled:center"), bb)
  re <- list(beta=beta,iter = iter)
  return (re)
}


###############################################################
####################utility function ############################
############################################################
#### K-M weight
kmw<-function(y,delta){
  y_s=y
  delta_s=delta
  kmweight<-c()
  nw<-length(y)
  
  comb<-cbind(y,delta)
  oo=order(y)
  ocomb<-comb[oo,]
  y<-ocomb[,1]
  delta<-ocomb[,2]
  kmweight[1]<-delta[1]/nw
  for(ind in 2:nw){
    tmp<-c()
    for(ind2 in 1:(ind-1)){
      tmp[ind2]<-((nw-ind2)/(nw-ind2+1))^delta[ind2]
    }
    kmweight[ind]<-delta[ind]/(nw-ind+1)*prod(tmp)
  }
  kmweight1=matrix(0,nw,0)
  kmweight1[oo]=kmweight
  return(kmweight=kmweight1)
}


tpfp <- function(predict,true){
  p <- length(true)
  signumber <- sum(true!=0)
  TP = sum(true!=0& predict!=0)
  FP =  sum(true==0 & predict !=0)
  t.FP =FP/(p-signumber)
  t.TP = TP/signumber
  return(c(FPR=t.FP, TPR=t.TP))
}

#### tpfp, 
##mat： matrix, tpfp matrix, 2*L
#beta: vector, the true parameter, 
#@top: int, m.  the true significance variable of top m selected variable
eval.result <- function(mat, beta, top){
  
  # mat <- ores[43,"G",,,"lasso"]
  Auc <-  pAuc(fpr = mat["FPR",], tpr =mat["TPR",],partial = c(0,1)) 
  pAuc.half <-  pAuc(fpr = mat["FPR",], tpr =mat["TPR",],partial = c(0, 0.5))/0.5 
  pAuc.third <- pAuc(fpr = mat["FPR",], tpr =mat["TPR",],partial = c(0, 0.3))/0.3
  p <- length(beta)
  P <- sum(beta != 0)
  N <- p - P
  ## top 
  S <- mat[1, ]* N + mat[2,]*P
  top.ind <- which(S >= top)
  S.min <- which.min(S[top.ind])
  top.tp <-ifelse(length(S.min)==0, 0, mat[2,top.ind[S.min]] * P )
  
  ##cover.min
  # cover.min <- p
  # tp1 <- which(mat[2,] ==1)
  # fp.min <- which.min(mat[1, tp1])
  # cover.min <- ifelse(length(fp.min) ==0, p,mat[1,tp1[fp.min]]*N + P )
  
  return(c(AUC = Auc, pAUC.half = pAuc.half, pAUC.third = pAuc.third, 
           Top= top.tp))
  
}

#compute the partial auc
#@ fpr, fpr sequence, vector: L*1  
#@ tpr, tpr sequence, vector: L*1
#@partial, vector: 2*1, (r,s)
pAuc<- function(fpr,tpr, partial=c(0,1)){
  # fpr = Sfpc
  # tpr = Stpc
  pts <- data.frame(fpr = as.numeric(fpr),tpr=as.numeric(tpr))
  sortpts <- pts[order(fpr),]
  low_p = partial[1]
  high_p = partial[2]
  psortpts = dplyr::filter(sortpts,(fpr <= high_p) & (fpr>= low_p))
  pauc<- sum((psortpts[-1,2] + psortpts[-nrow(psortpts),2]) * diff(psortpts[,1]))/2
  return(pauc)
}

