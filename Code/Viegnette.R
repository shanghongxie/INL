
## Note: the most updated code is on Github: http://github.com/shanghongxie/INL

### Example : target networks varied by both external modality network and covariates, p=5


library(Rcpp)
library(RcppEigen)

library(MASS)



###################################
#####  Data generation part  ######
###################################
source('~/GenData.R')


N=200

p=5

set.seed(123)
distance=matrix(0,p,p)

dist=sample(1:20,size=p*(p-1)/2,replace=T)
# 
dist=dist/10


s=1
for (j in 1:(p-1)){
  for (k in (j+1):p){
    distance[j,k]=dist[s]
    s=s+1
  }
}

distance=distance+t(distance)
isSymmetric(distance)

diag(distance)=1



### compute expectation E(Bijk|Mi), E(BijkBijk'|Mi)

## Beta: q*p(p-1), ignore beta_jj=0; prob: N*p(p-1) matrix, symmetric, obtained from previous step, using gamma(t-1),eta(t-1)
q=3
p=5

Beta=matrix(0,q,p*(p-1))

indexsr=NULL
sr=1
u=1
u2=1
for (s in 1:p){
  for (r in 1:p){
    if (r!=s){
      indexsr=rbind(indexsr,c(s,r,sr))
      sr=sr+1
    }
  }
}

######################################################
########  Setting 3  weak signal  ########
######################################################

## X~N(0.1,1), truncated at (-1,1)

## beta_12
Beta[,1]=0.5*c(1,-1,1)
## beta_14
Beta[,3]=0.5*c(1,1,-1)
## beta_23
Beta[,6]=0.5*c(1,1,1)
## beta_25
Beta[,8]=0.5*c(-1,1,1)


## beta_21
Beta[,5]=Beta[,1]
## beta_41
Beta[,13]=Beta[,3]
## beta_32
Beta[,10]=Beta[,6]
## beta_52
Beta[,18]=Beta[,8]

sigma=c(0.2,0.3,0.2,0.3,0.2)

sum(apply(Beta,2,function(x){sum(x^2)})!=0)  ## 4



######################################################
########     Setting 4 large signal           ########
######################################################

## X~N(0.1,1), truncated at (-1,1)

## beta_12
Beta[,1]=c(1,-1,1)
## beta_14
Beta[,3]=c(1,1,-1)
## beta_23
Beta[,6]=c(1,1,1)
## beta_25
Beta[,8]=c(-1,1,1)


## beta_21
Beta[,5]=Beta[,1]
## beta_41
Beta[,13]=Beta[,3]
## beta_32
Beta[,10]=Beta[,6]
## beta_52
Beta[,18]=Beta[,8]

sigma=c(0.2,0.3,0.2,0.3,0.2)

sum(apply(Beta,2,function(x){sum(x^2)})!=0)  ## 4

### generate covariates X, q-dimensional, standard norm distribution
g1=NULL
e1=NULL
Xs=list()
Ms=list()
Bis=list()
Ss=list()

### sigma^2
sigmaM=diag(sigma,nrow=p,ncol=p)

###  index of non-zero in Beta, Beta_jk^T!=0 vector
indexBeta=which(apply(Beta,2,function(x) sum(abs(x)))!=0)



gamma=2
eta=1

## row: subject and node, column: node
### common DTI network across subjects (need symmetric; undirected graph)
isActiveC1=matrix(0,p,p)

isActiveC1[1,]=c(0,0,1,1,0)
isActiveC1[2,]=c(0,0,1,1,1)
isActiveC1[3,]=c(1,1,0,0,0)
isActiveC1[4,]=c(1,1,0,0,0)
isActiveC1[5,]=c(0,1,0,0,0)


library(igraph)
par(mfrow=c(1,1))
DTIgraph=graph_from_adjacency_matrix(isActiveC1, mode = "undirected", weighted = NULL, diag = TRUE,
                                     add.colnames = NULL, add.rownames = NA)

plot(DTIgraph)

set.seed(123)
data=GenData(N, p, q, Beta, gamma, eta, distance, sigma, isActiveC1)

X=data$X; M=data$M; B=data$B; S=data$S; ix=data$ix; im=data$im; isActiveC=data$isActiveC
maxeign=c(maxeign,data$maxeigen)


##############################
#####  Estimation part  ######
##############################


gamma0s=seq(1.8,2.2,by=0.1)

eta0=1


srtable=NULL
sr=0
for (s in 1:length(gamma0s)){
  for (r in 1:length(eta0s)){
    sr=sr+1
    srtable=rbind(srtable,c(s,r,sr))
  }
}


thresh=1e-5


Beta1s=NULL
Sigma1s=NULL
SSEs=NULL

g1s=matrix(NA,nSim,p)
e1s=matrix(NA,nSim,p)

llist=list()
Objlist=list()

its=NULL


gammanode=NULL
etanode=NULL

po=p*(p-1)

EBsims=list()
EBBsims=list()


### initial 
p=5; q=3
Beta0=matrix(0,q,p*(p-1))
Sigma0=rep(1,p)

threshIn=thresh

maxit=100


iter=2
### first E-stel
times=matrix(NA,nSim,3)
MSEs=NULL
RSSs=NULL

thrlow=0.1
thrhigh=0.9

##########################################################################################

###  Use Approximation Approach to calculate posterior expectation of latent variables 

##########################################################################################

sourceCpp('~/INLApproxC.cpp')
source('~/INLApproxRcode.R')
source('~/INLApproxHardThrRcode.R')


u=seq(-1.5,1.5,by=0.0001)
lambda1=c(0.001, 0.002, 0.005, 0.02, 0.1, 0.2, 0.5, 1,  2, 3,  4,  5, 8, 10, 16, 25, 32, 64)
lambda=lambda1
lambda=lambda1
Z=matrix(NA,length(u),length(lambda))

for (i in 1:length(u)){
  Z[i,]=exp(-lambda*u[i])
}

Y=exp(-u^2)

fit=lm(Y~Z-1)
coeff1=fit$coefficients


maxitnum=matrix(0,p,100)
count=1

indexfun=indexjC(p)

indexjk=indexfun$indexjk
indexjkk=indexfun$indexjkk



###################################################
###    Step 1: Estimation no pruning            ###
###################################################

result=INLApproxRcode(X, M, S, distance, Beta0, Sigma0, gamma0s, eta0, iter, maxit, thresh,
                  lambda1, coeff1, iSim)

estBeta = result$estBeta
estSigma = result$estSigma
estgammas = result$gammanode
estgamma = mean(result$gammanode)
esteta = mean(result$etanode)
EB = result$EBsim
EBB = result$EBBsim


#### Do hard-thresholding after estimate all the betas

###################################################
###    Step 2: Hard thresholding based on EBIC  ###
###################################################

### the estimates of beta and sigma from Step 1 are used as initial estimates in Step 2
### do hard-thresholding

# install.packages("pracma")
library(pracma)

## convert to undirected graph

Betau2=matrix(0,q,po/2)
indexsr=NULL
sr=1
u=1
u2=1
for (s in 1:p){
  for (r in 1:p){
    if (r!=s){
      indexsr=rbind(indexsr,c(s,r,sr))
      sr=sr+1
    }
  }
}


u=1
Betau=matrix(0,q,po/2)
for (s in 1:(p-1)){
  for (r in (s+1):p){
    sr=which(indexsr[,1]==s & indexsr[,2]==r)
    rs=which(indexsr[,1]==r & indexsr[,2]==s)
    
    Betau[,u]=(Beta[,sr]+Beta[,rs])/2
    u=u+1
  }
}

truedge=(apply(Beta^2,2,sum)!=0)*1
truedgeu=(apply(Betau^2,2,sum)!=0)*1

### set prob cutoff of Bijk, if beta>0 and P(Bijk=1|Cij)>p, then the edge presented

probcut=seq(0,0.4,by=0.1)

cutnum=50
hyper=0.5

Alpha=matrix(0,q,po)
for (nodej in 1:p){
  Alpha[,((nodej-1)*(p-1)+1):((nodej-1)*(p-1)+(p-1))]=Beta[,((nodej-1)*(p-1)+1):((nodej-1)*(p-1)+(p-1))]*sigma[nodej]
}

indexfun=indexjC(p)

indexjk=indexfun$indexjk
indexjkk=indexfun$indexjkk

  
hardthresult=INLApproxHardThrR(X, M, Bi, S, EB, EBB, Beta, distance, gamma0s, eta0, estBeta, estSigma, estgammas, eta0, iter, maxit, thresh,
                            probcut, indexsr, indexjk, indexjkk, cutnum, hyper, lambda1, coeff1)
  
Betaht=hardthresult$HTBeta

Sigmaht=hardthresult$HTSigma

Alphaht=matrix(0,q,po)
for (nodej in 1:p){
  Alphaht[,((nodej-1)*(p-1)+1):((nodej-1)*(p-1)+(p-1))]=Betaht[,((nodej-1)*(p-1)+1):((nodej-1)*(p-1)+(p-1))]*Sigmaht[nodej]
}

MSEa=sum((Alphaht-Alpha)^2)

## true positives corresponding to each alpha value (probcut)
TP=hardthresult$TPus
FP=hardthresult$FPus
TN=hardthresult$TNus
FN=hardthresult$FNus

AUC=hardthresult$auctotal
  
 
##########################################################################################

###### Calculate posterior expection of the latent variables directly (DIRECT)    ######

##########################################################################################


sourceCpp('~/INLDirectC.cpp')
source('~/INLDirectRcode.R')
source('~/INLDirectHardThrRcode.R')



library(e1071)
B9=bincombinations(9)
B8=bincombinations(8)
B7=bincombinations(7)
B6=bincombinations(6)
B5=bincombinations(5)
B4=bincombinations(4)
B3=bincombinations(3)
B2=bincombinations(2)
B1=bincombinations(1)


maxit=100
thrlow=0.1
thrhigh=0.9

###################################################
###    Step 1: Estimation no pruning            ###
###################################################

Directresult=INLDirectRcode(X, M, S, distance, Beta0, Sigma0, gamma0s, eta0s, iter, thrlow, thrhigh, maxit, thresh,
                                               B1,B2,B3,B4,B5,B6,B7,B8,B9)

estBeta = Directresult$estBeta
estSigma = Directresult$estSigma
estgammas = Directresult$gammanode
estetas =  Directresult$etanode
EB = Directresult$EBsim
EBB = Directresult$EBBsim

###################################################
###    Step 2: Hard thresholding based on EBIC  ###
###################################################
hyper=0.5
Directhardthresult=INLDirectHardThrR(X, M, Bi, S, EB, EBB, Beta, distance, gamma0s, eta0s, estBeta, estSigma,estgammas, estetas, iter, maxit, thresh,
                                                              B1,B2,B3,B4,B5,B6,B7,B8,B9, probcut, indexsr, thrlow, thrhigh, indexjk, indexjkk, cutnum, hyper,srtable)


Betaht = Directhardthresult$HTBeta
Sigmaht = Directhardthresult$HTSigma

Alphaht=matrix(0,q,po)
for (nodej in 1:p){
  Alphaht[,((nodej-1)*(p-1)+1):((nodej-1)*(p-1)+(p-1))]=Betaht[,((nodej-1)*(p-1)+1):((nodej-1)*(p-1)+(p-1))]*Sigmaht[nodej]
}

MSEa=sum((Alphaht-Alpha)^2)

## true positives corresponding to each alpha value (probcut)
TP=Directhardthresult$TPus
FP=Directhardthresult$FPus
TN=Directhardthresult$TNus
FN=Directhardthresult$FNus

AUC=Directhardthresult$auctotal

