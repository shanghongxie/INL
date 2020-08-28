

library(pracma)
INLDirectHardThrR<-function(X, M, Bi, S, EB, EBB, Beta, distance, gamma0s, eta0s, Beta1, Sigma1, g1, e1, iter, maxit, thresh,
                      B1,B2,B3,B4,B5,B6,B7,B8,B9, probcut, indexsr, thrlow, thrhigh, indexjk, indexjkk, cutnum, hyper, srtable){
  
  N=nrow(X); q=ncol(X); p=ncol(M); po=p*(p-1)
  
 
  
  TPus=NULL
  FPus=NULL
  TNus=NULL
  FNus=NULL
  FDRus=NULL
  senus=NULL
  specus=NULL
  
 
  
  edgejs=NULL


  
  for (nodej in 1:p){
    Beta1[,((nodej-1)*(p-1)+1):((nodej-1)*(p-1)+(p-1))]=Beta1[,((nodej-1)*(p-1)+1):((nodej-1)*(p-1)+(p-1))]*Sigma1[nodej]
  }
  
  BetajSim=matrix(Beta1,q,po)
  
  ### take average of beta_jk and beta_kj
  AveBetajSim=matrix(0,q,po/2)
  u=1
  for (s in 1:(p-1)){
    for (r in (s+1):p){
      sr=which(indexsr[,1]==s & indexsr[,2]==r)
      rs=which(indexsr[,1]==r & indexsr[,2]==s)
      
      AveBetajSim[,u]=(BetajSim[,sr]+BetajSim[,rs])/2
      ### fill in average beta_jk in Beta matrix
      BetajSim[,sr]=BetajSim[,rs]=AveBetajSim[,u]
      u=u+1
    }
  }
  
  
  Beta3=matrix(0,q,po)
  Sigma3=rep(0,p)
  
  seqbetaj=sort(c(abs(BetajSim)))
  seqbetaj=seqbetaj[seqbetaj!=0]
  
  seqbetaj=seq(from=max(seqbetaj), to=min(seqbetaj),length.out=cutnum)
  
  
  ### average EB_jk and EB_kj
  AveEB=matrix(0,N,po)
  for (s in 1:(p-1)){
    for (r in (s+1):p){
      sr=which(indexsr[,1]==s & indexsr[,2]==r)
      rs=which(indexsr[,1]==r & indexsr[,2]==s)
      
      ### fill in average beta_jk in Beta matrix
      AveEB[,sr]=AveEB[,rs]=(EB[,sr]+EB[,rs])/2
    }
  }
  
  ### average over subjects
  averageEB=apply(AveEB,2,mean)
  
  ### true edges for each node j
  tedges=matrix(0,p,p-1)

  for (nodej in 1:p){
    tedges[nodej,]=(apply(Beta[,((nodej-1)*(p-1)+1):((nodej-1)*(p-1)+(p-1))]^2,2,sum)!=0)*1
  }
  
  
  probB=apply(Bi,2,sum)/N
  
  truedgeB=probB*truedge
  truedgeuB=rep(0,po/2)
  u=1
  for (s in 1:(p-1)){
    for (r in (s+1):p){
      sr=which(indexsr[,1]==s & indexsr[,2]==r)
      rs=which(indexsr[,1]==r & indexsr[,2]==s)
      
      truedgeuB[u]=(truedgeB[sr]+truedgeB[rs])/2
      u=u+1
    }
  }
  
  
  ### true edges with DTI network
  tedgeBs=tedges*matrix(probB,p,p-1)
  
  
  senjs=matrix(0,cutnum*length(probcut),p)
  specjs=matrix(0,cutnum*length(probcut),p)
  
  TPtotal=rep(0,cutnum*length(probcut))
  FPtotal=rep(0,cutnum*length(probcut))
  TNtotal=rep(0,cutnum*length(probcut))
  FNtotal=rep(0,cutnum*length(probcut))
  
  for (nodej in 1:p){
    
    # gamma1=g1[nodej]
    # eta1=e1[nodej]
    
    gamma1=mean(g1)
    
    eta1=eta0s
    
    ### extract beta_j for node j
    inibetaj=c(BetajSim[,((nodej-1)*(p-1)+1):((nodej-1)*(p-1)+(p-1))])
    
    inisigmaj=Sigma1[nodej]
    
    betaK=matrix(0,length(inibetaj),cutnum)
    
    
    ##### calculate AUC for different hard-thresholding values
    count=1
    for (pc in 1:length(probcut)){  
      for (rk in 1:cutnum){ 
        
        ### truncate 
        betaK[,rk]=inibetaj*(abs(inibetaj)>=seqbetaj[rk])
        
        betajrk=matrix(betaK[,rk],q,p-1)
        
        # estedge=(apply(betajrk^2,2,sum)>0)*1
        
        estedge=(apply(betajrk^2,2,sum)>0)*averageEB[((nodej-1)*(p-1)+1):((nodej-1)*(p-1)+(p-1))]
        
        TPj=sum((estedge>probcut[pc]) & (tedgeBs[nodej,]>probcut[pc]))
        
        ## calculate TP over all nodes
        TPtotal[count]=TPtotal[count]+TPj
        # TPjs=c(TPjs,TPj)
        
        FPj=sum((estedge>probcut[pc]) & (tedgeBs[nodej,]<=probcut[pc]))
        # FPjs=c(FPjs, FPj)
        FPtotal[count]=FPtotal[count]+FPj
        
        TNj=sum((estedge<=probcut[pc]) & (tedgeBs[nodej,]<=probcut[pc]))
        # TNjs=c(TNjs,TNj)
        TNtotal[count]=TNtotal[count]+TNj
        
        FNj=sum((estedge<=probcut[pc]) & (tedgeBs[nodej,]>probcut[pc]))
        # FNjs=c(FNjs,FNj)
        FNtotal[count]=FNtotal[count]+FNj
        
        senj=TPj/(TPj+FNj)
        senjs[count,nodej]=senj
        
        specj=TNj/(TNj+FPj)
        specjs[count,nodej]=specj
        
        count=count+1
        
      }  ## end rk (vary hard-thresholding values)
      
    }  ### end of probcut
    
    
    fprj=1-specjs[,nodej]
    fprj=sort(fprj)
    
    tprj=sort(senjs[,nodej])
    

    
    EBj=EB[,((nodej-1)*(p-1)+1):((nodej-1)*(p-1)+(p-1))]
    
    EBj1=(EBj>thrhigh)*1
    EBj0=(EBj<thrlow)*1
    EBji=(EBj<=thrhigh & EBj>thrlow)*1
    
    
    hardthresh=HardThresholdDTIjC(inibetaj, inisigmaj, seqbetaj, betaK, X, M, EB, EBB,  EBji, EBj1,
                                  indexjk, indexjkk, S, distance,  gamma1, eta1, B1, B2, B3,
                                  B4, B5, B6, B7, B8, B9, cutnum, nodej-1)
    
    
    sigmajs=hardthresh$sigmajs
    
    llj=hardthresh$Likeli
    
    
    #### calculate EBIC to choose lambda
    
    EBIC=NULL
    
    ### number of non-zero ||beta_jk||_2 (L2 norm)
    E=rep(0,cutnum)
    
    for (rk in 1:cutnum){
      Betark=matrix(betaK[,rk],q,p-1)
      ### L2 norm
      E[rk]=sum(apply(Betark^2,2,sum)!=0)
    }
    
    
    ## set hyperparameter=0.5
    # hyper=0.5
    EBIC=-2*llj+E*log(N)+2*hyper*E*log(p-1)
    
    indexmin=which.min(EBIC)
    
    ### minimize EBIC
    betaj=matrix(betaK[,indexmin],q,p-1)
    
    ### number of edges
    edgejs[nodej]=E[indexmin]
    
    sigmaj=sigmajs[indexmin]
    
    Beta3[,((nodej-1)*(p-1)+1):((nodej-1)*(p-1)+(p-1))]=betaj/sigmaj
    Sigma3[nodej]=sigmaj
    
    
  }  ### end of node j
  
  TPtotal=TPtotal/2
  FPtotal=FPtotal/2
  FNtotal=FNtotal/2
  TNtotal=TNtotal/2
  
  sentotal=TPtotal/(TPtotal+FNtotal)
  
  spectotal=TNtotal/(TNtotal+FPtotal)
  
  fprtotal=1-spectotal
  fprtotal=sort(fprtotal)
  fprtotal=c(0,fprtotal,1)
  
  tprtotal=sort(sentotal)
  tprtotal=c(0,tprtotal,1)
  
  ### AUC total varying hardthreshold and proportion
  auctotal=trapz(fprtotal,tprtotal)
  
  
  
  ### take average of beta_jk and beta_kj
  AveBeta3=matrix(0,q,po/2)
  u=1
  for (s in 1:(p-1)){
    for (r in (s+1):p){
      sr=which(indexsr[,1]==s & indexsr[,2]==r)
      rs=which(indexsr[,1]==r & indexsr[,2]==s)
      
      AveBeta3[,u]=(Beta3[,sr]+Beta3[,rs])/2
      Beta3[,sr]=Beta3[,rs]=(Beta3[,sr]+Beta3[,rs])/2
      
      u=u+1
    }
  }
  
  
  
  
  ### average EB_jk and EB_kj (undirected)
  AveEBu=matrix(0,N,po/2)
  u=1
  for (s in 1:(p-1)){
    for (r in (s+1):p){
      sr=which(indexsr[,1]==s & indexsr[,2]==r)
      rs=which(indexsr[,1]==r & indexsr[,2]==s)
      
      AveEBu[,u]=(EB[,sr]+EB[,rs])/2
      u=u+1
    }
  }
  
  averageEBu=apply(AveEBu,2,mean)
  
  ### not average
  edge=(apply(Beta3^2,2,sum)!=0)*averageEB
  
  ### average
  edgeu=(apply(AveBeta3^2,2,sum)!=0)*averageEBu

  
  for (pc in 1:length(probcut)){
  
    
    
    ##### undirected edges
    TPus[pc]=sum((edgeu>probcut[pc]) & (truedgeuB>probcut[pc]))
    FPus[pc]=sum((edgeu>probcut[pc]) & (truedgeuB<=probcut[pc]))
    TNus[pc]=sum((edgeu<=probcut[pc]) & (truedgeuB<=probcut[pc]))
    FNus[pc]=sum((edgeu<=probcut[pc]) & (truedgeuB>probcut[pc]))
    FDRus[pc]=FPus[pc]/(FPus[pc]+TPus[pc])
    
    senus[pc]=TPus[pc]/(TPus[pc]+FNus[pc])
    specus[pc]=TNus[pc]/(TNus[pc]+FPus[pc])
    
    
    
  }  ### end of probcut
  

  SSEus2=0
  MSEs2=0
  
  # 
  output=list(HTBeta=Beta3, HTSigma=Sigma3, 
              TPus=TPus, FPus=FPus, TNus=TNus, FNus=FNus, FDRus=FDRus, senus=senus, specus=specus,
              sentotal=sentotal,spectotal=spectotal,TPtotal=TPtotal,
              FPtotal=FPtotal,
              FNtotal=FNtotal,
              TNtotal=TNtotal, fprtotal=fprtotal,
              tprtotal=tprtotal,
              auctotal=auctotal)
  
  
  return (output)
  
  
  
}
