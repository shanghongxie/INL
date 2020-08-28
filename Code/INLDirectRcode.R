


INLDirectRcode<-function(X, M, S, distance, Beta0, Sigma0, gamma0s, eta0s, iter, thrlow, thrhigh, maxit, thresh,
                   B1,B2,B3,B4,B5,B6,B7,B8,B9){
  
  N=nrow(X); q=ncol(X); p=ncol(M)
  # ### store the full betas and sigmas of all the nodes
  Beta3=matrix(0,q,p*(p-1))
  Sigma3=rep(1,p)
  
  po=p*(p-1)
  EBsim=matrix(0,N,po)
  EBBsim=matrix(0,N,p*(p-1)*(p-2)/2)
  
  llj=matrix(NA, length(gamma0s)*length(eta0s),p)
  
  itsrs=matrix(NA, length(gamma0s)*length(eta0s),p)
  
  its=NULL

  
  gammanode=NULL
  etanode=NULL

  
  flag=matrix(0,p,length(gamma0s))
  
  for (nodej in 1:p){
    ### step 1
    
    itsr=NULL
    sr=0
    ll=NULL
    
    
    
    ## store Beta, sigma at each gamma/eta
    Betasr=list()
    Sigmasr=list()
    
    EBs=list()
    EBBs=list()
    
    for (s in 1:length(gamma0s)){
      for (r in 1:length(eta0s)){
        
        Beta1=Beta0
        Sigma1=Sigma0
        
        sr=sr+1
        
        it=0
        while(1){
          it=it+1
          # for (it in 1:27){  
          
          ### did not use prob here
          estep=Estep1jC(Beta1, X, M, thresh, nodej-1)
          
          indexAct=estep$indexAct
          indexjk=estep$indexjk
          indexjkk=estep$indexjkk
          # 
          XBM=estep$XBM
          UMi=estep$UMi
          UMMi=estep$UMMi
          
          range(XBM)
          range(UMi)
          range(UMMi)
          
          
          ### E step (used prob here)
          
          if (it<iter){
            
            expect=ExpectjC(Beta1, X, M,
                            indexAct, indexjk, indexjkk, 
                            XBM, UMi, UMMi, B8, B7,
                            Sigma1,  S,  distance, gamma0s[s],  eta0s[r],  nodej-1)
            
          } else{
            
            XBMj=XBM[,((nodej-1)*(p-1)+1):((nodej-1)*(p-1)+(p-1))]
            
            #### posterior
            EBj=EB[,((nodej-1)*(p-1)+1):((nodej-1)*(p-1)+(p-1))]
            
            EBj1=(EBj>thrhigh)*1
            EBj0=(EBj<thrlow)*1
            EBji=(EBj<=thrhigh & EBj>thrlow)*1
            
            expect=ExpectDTIjC(M, EBji, EBj1,
                               indexjk, indexjkk, 
                               XBMj, B1, B2, B3,
                               B4, B5, B6, B7, B8, B9,
                               Sigma1,  S, distance, gamma0s[s], eta0s[r], nodej-1)
            
          }
          
          EBB=expect$EBB
          EB=expect$EB
          
          EB=round(EB,7)
          EBB=round(EBB,7)
          
          range(EB)
          range(EBB)
          ### M step, update Beta only
          mstep2=MstepjC(X, M, EB, EBB, 
                         indexjk,indexjkk, nodej-1)
          
          Beta2=mstep2$Beta
          
          Sigma2=mstep2$Sigma
          
          if (any(Sigma2<0) || any(is.na(Sigma2)) )
          { flag[nodej,sr]=1
          break; 
          }else{
            
            if (it>maxit) {Beta1=Beta2
            Sigma1=Sigma2; break};
            
            if ((sum((Beta2-Beta1)^2)<thresh*sum(Beta2^2)) & (sum((Sigma2-Sigma1)^2)<thresh*sum(Sigma2^2)))
            {
              
              Beta1=Beta2
              Sigma1=Sigma2
              
              break;
            }   
            
            Beta1=Beta2
            Sigma1=Sigma2
            
          }
          
          
        } ### end of EM, end of while
        
        ### if estimation is wrong, set likelihood=-inf
        if (any(Sigma2<0) || any(is.na(Sigma2))){
          
          ll[sr]=-Inf
          llj[sr,nodej]=ll[sr]
          
        }else{
          
          
          itsrs[sr,nodej]=it
          
          EBs[[sr]]=EB
          EBBs[[sr]]=EBB    
          
          itsr[sr]=it
          
          Betasr[[sr]]=Beta1
          Sigmasr[[sr]]=Sigma1
          
          
          ### update XBM using the final beta
          estep=Estep1jC(Beta1, X, M, thresh, nodej-1)
          
          XBM=estep$XBM
          
          
          ### after update beta in EM
          #### calculate likelihood function
          
          XBMj=XBM[,((nodej-1)*(p-1)+1):((nodej-1)*(p-1)+(p-1))]
          
          #### posterior
          EBj=EB[,((nodej-1)*(p-1)+1):((nodej-1)*(p-1)+(p-1))]
          
          EBj1=(EBj>thrhigh)*1
          EBj0=(EBj<thrlow)*1
          EBji=(EBj<=thrhigh & EBj>thrlow)*1
          
          likeli=LikeliObsDTIjC(M, EBji, EBj1,
                                indexjk, indexjkk, 
                                XBMj, B1, B2, B3,
                                B4, B5, B6, B7, B8, B9,
                                Sigma1, S, distance, gamma0s[s], eta0s[r], nodej-1)
          
          ll[sr]=likeli$Likeli
          
          llj[sr,nodej]=likeli$Likeli
          
        }
        
      }  ### end of r
    }  ### end of s
    
    
    if (!all(ll==-Inf)) {
      
      ### find optimal gamma, eta 
      srmax=which.max(ll)
      
      gamma1=gamma0s[srtable[srmax,1]]
      eta1=eta0s[srtable[srmax,2]]
      
      gammanode[nodej]=gamma1
      etanode[nodej]=eta1
      
      Beta1=Betasr[[srmax]]
      Sigma1=Sigmasr[[srmax]]
      
      its[nodej]=itsr[srmax]
      
      EBsim[,((nodej-1)*(p-1)+1):((nodej-1)*(p-1)+(p-1))]=EBs[[srmax]][,((nodej-1)*(p-1)+1):((nodej-1)*(p-1)+(p-1))]
      EBBsim[,((nodej-1)*(p-1)*(p-2)/2+1):((nodej-1)*(p-1)*(p-2)/2+(p-1)*(p-2)/2)]=EBBs[[srmax]][,((nodej-1)*(p-1)*(p-2)/2+1):((nodej-1)*(p-1)*(p-2)/2+(p-1)*(p-2)/2)]
      
      Beta3[,((nodej-1)*(p-1)+1):((nodej-1)*(p-1)+(p-1))]=Beta1[,((nodej-1)*(p-1)+1):((nodej-1)*(p-1)+(p-1))]/Sigma1[nodej]
      Sigma3[nodej]=Sigma1[nodej]
      
    }
    
    
  } ### end of nodej
  

  ### if for any node, llj=-inf all all gamma values
  if (any(apply(llj,2,function(x){all(x==-Inf)}))){
    output=NULL
  }else{
    
    output=list(estBeta=Beta3, estSigma=Sigma3, gammanode=gammanode, etanode=etanode,EBsim=EBsim, EBBsim=EBBsim, time=time, its=its, itsrs=itsrs, llj=llj)
  }
  

  
  return (output)
  
}

