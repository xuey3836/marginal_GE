rm(list = ls())
library("ncvreg")
library("Rcpp")
source("function.R")
source("gen_data.R")
sourceCpp("Rcpp_fun.cpp")

#############################################################
##################### lambda squenecing#######################
###############################################g##############
## ---------------------LASSO,MCP,SCAD---lambda----------
nlambda = 100
lambda <- exp(seq(log(0.8),log(0.1),len= nlambda))
##-------------------sgMCP,sgLASSO ------lambda----------
nlambda1  = 30
nlambda2 = 50
lambda1 = exp(seq(log(1.5),log(0.1),len= nlambda1)) 
lambda2 =exp(seq(log(1.5),log(0.03),len= nlambda2))
Matla = cbind(sort(rep(lambda1,nlambda2),decreasing = TRUE),rep(lambda2,nlambda1))
rownames(Matla) <- apply(Matla,1,function(dat) paste0(dat[1],"-",dat[2]))
L <- nrow(Matla)


#####################################################################
#######################    Generate data  #############################
###################################################################
#-----------------------------setting ------------------------------
n.seq <- c(100, 500)
p <- 500
q <- 4
Nrep <- 20 # replicate time
## the signifance G and GE
signumberG <- 16
signumberGE <- 12
signumberE <- 4
setG <- sample(1:p, size = signumberG, replace = FALSE)
setE <- sample(1:q, size = signumberE, replace = FALSE)
####---------hierichical structure-----
se = cbind(sort(rep(setG,signumberE)),rep(setE,signumberG))
setGE = matrix(se[sample(1:nrow(se),size = signumberGE,replace = F),],nrow= signumberGE,ncol=2)

family.seq <- c("gaussian", "binomial")
corr.seq <- c("AR1","B2")
err.seq <- c("Norm","M_NL","t")

# the coefficient
alpha0 = 0
alpha = rep(0,q)
beta = rep(0,p)
gamma = matrix(0,p,q)
alpha[setE] = runif(signumberE,0.2,0.8)
beta[setG] = runif(signumberG, 0.2,0.8)
gamma[setGE] = runif(signumberGE, 0.2,0.8)
a = alpha
b = as.vector(t(cbind(beta,gamma)))
vgamma = as.vector(gamma)

family <- family.seq[1]
n <- n.seq[1]
Gcorr <- corr.seq[1]
Ecorr <- corr.seq[1]
err.dis <- err.seq[1]

##G/E continuous or categorical
is.SNP <- FALSE
E.cat <- 1

##generate data
dat <- GenGE(n, p, q, Gcorr=Gcorr, Ecorr=Ecorr, err.dis, alpha0=0, a=a, b=b,
			 family=family, is.SNP=is.SNP, E.cat=E.cat)
y <- dat$Y
Z <- dat$Z
X <- dat$X
W <- dat$W


##### AFT model :construct W
# library('survival')
# kmweight=kmw(y,delta)

##############################################
################Fit the model ################
##############################################
#sig
siguniv = Sig(X = X, W=W, y = y, family = family, fdr.level = 0.1,censor = FALSE)
##MCP
family = "gaussian"
############# MCP ##################
fit.MCP <- mMCP(y=y, X=X, W=W, lambda=lambda, family=family, penalty="MCP", weight = NULL)

############ sgMCP ##############
fit.sgMCP <- msgMCP(y=y, X=X, W=W, Matla = Matla, family=family, weight = NULL)

#################result ########################
## tpfp + auc, pauc, top20
##Gene
G.tpfp <- apply(fit.sgMCP$beta, 2, tpfp, true=beta)
##Gene- envrionment
GE.tpfp <- apply(apply(fit.sgMCP$gamma,3,as.vector),2,tpfp,true=vgamma)

result <- eval.result(mat =G.tpfp, beta = beta, top=20 )
result