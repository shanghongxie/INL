
INLApproxRcode<-function(X, M, S, distance, Beta0, Sigma0, gamma0s, eta0, iter, maxit, thresh,
                   lambda1, coeff1, iSim){  

  N=nrow(X); q=ncol(X); p=ncol(M)
  # ### store the full betas and sigmas of all the nodes
  Beta3=matrix(0,q,p*(p-1))
  Sigma3=rep(1,p)
  
  po=p*(p-1)
  EBsim=matrix(0,N,po)
  EBBsim=matrix(0,N,p*(p-1)*(p-2)/2)
  
  llj=matrix(NA, length(gamma0s),p)
  
  itsrs=matrix(NA, length(gamma0s),p)
  
  indexfun=indexjC(p)
  
  indexjk=indexfun$indexjk
  indexjkk=indexfun$indexjkk
  
  its=NULL

  gammanode=NULL
  etanode=NULL
  
  numlike=NULL
  
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
      # for (r in 1:length(eta0s)){
      
      
      prob=probC(N, p, S, distance, gamma0s[s], eta0)
      
      Beta1=Beta0
      Sigma1=Sigma0
      
      sr=sr+1
      
      
      it=0
      while(1){
        it=it+1

        ### did not use prob here
        estep=Estep1jC(Beta1, X, M, indexjk, indexjkk, thresh, nodej-1)
        
        # 
        XBM=estep$XBM
        UMi=estep$UMi
        UMMi=estep$UMMi
        
       
        ### E step (used prob here)

        
        expect=ExpectApproxjC(Beta1, X, M,
                              indexjk, indexjkk, 
                              XBM, UMi, UMMi, 
                              Sigma1, prob, lambda1,  coeff1, nodej-1)
        
   
        PB1=expect$PB1
        PB2=expect$PB2
        
        PBB1=expect$PBB1
        PBB2=expect$PBB2
        PBB3=expect$PBB3
        PBB4=expect$PBB4
        
        EBB=expect$EBB
        EB=expect$EB
        
      
        EB=round(EB,7)
        EBB=round(EBB,7)
        
        EB[EB<0]=0
        EB[EB>1]=1
        # 
        EBB[EBB<0]=0
        EBB[EBB>1]=1
        # 
        range(EB)
        range(EBB)
        ### M step, update Beta only
        mstep2=MstepjC(X, M, EB, EBB, 
                       indexjk,indexjkk, nodej-1)
        
       
        
        Beta2=mstep2$Beta
        
        Sigma2=mstep2$Sigma
        
        Sigma2
        
        if (any(Sigma2<0) || any(is.na(Sigma2)) )
        { flag[nodej,sr]=1
        break; 
        }else{
          if (it>maxit) 
          {Beta1=Beta2
          Sigma1=Sigma2; 
          break};
          
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
        estep=Estep1jC(Beta1, X, M, indexjk, indexjkk, thresh, nodej-1)
        
        XBM=estep$XBM
        
       
        ### after update beta in EM
        #### calculate likelihood function
        likeli=LikeliObsApproxjC(M, 
                                 indexjk, 
                                 XBM,
                                 Sigma1, S, distance, gamma0s[s], eta0,lambda1, coeff1, nodej-1)
        
        
        
        likehood=likeli$Likelihoodi
        
        numlike=c(numlike,sum(likehood<0))
        
        const=likeli$constant

        
        ll[sr]=sum(log(likehood[likehood>0]))/sum(likehood>0)
        

        
        llj[sr,nodej]=ll[sr]
        
      }
      
      # }  ### end of r
    }  ### end of s
    
    if (any(is.na(ll))){
      ll=rep(NA,length(gamma0s))
      
    }else{
      if (!all(ll==-Inf)) {
        
        ### find optimal gamma, eta 
        srmax=which.max(ll)
        
        
        gamma1=gamma0s[srmax]
        
        gammanode[nodej]=gamma1
        etanode[nodej]=eta0
        
        Beta1=Betasr[[srmax]]
        Sigma1=Sigmasr[[srmax]]
        
        its[nodej]=itsr[srmax]
        
        EBsim[,((nodej-1)*(p-1)+1):((nodej-1)*(p-1)+(p-1))]=EBs[[srmax]][,((nodej-1)*(p-1)+1):((nodej-1)*(p-1)+(p-1))]
        EBBsim[,((nodej-1)*(p-1)*(p-2)/2+1):((nodej-1)*(p-1)*(p-2)/2+(p-1)*(p-2)/2)]=EBBs[[srmax]][,((nodej-1)*(p-1)*(p-2)/2+1):((nodej-1)*(p-1)*(p-2)/2+(p-1)*(p-2)/2)]
        
        Beta3[,((nodej-1)*(p-1)+1):((nodej-1)*(p-1)+(p-1))]=Beta1[,((nodej-1)*(p-1)+1):((nodej-1)*(p-1)+(p-1))]/Sigma1[nodej]
        Sigma3[nodej]=Sigma1[nodej]
        
      } ## end of if
      
    }
    
    
  } ### end of nodej
 
  
  ### if for any node, llj=-inf all all gamma values
  if (any(is.na(llj)) || any(apply(llj,2,function(x){all(x==-Inf)}))){
    output=NULL
  }else{
    output=list(estBeta=Beta3, estSigma=Sigma3, gammanode=gammanode, etanode=etanode,EBsim=EBsim, EBBsim=EBBsim,its=its, itsrs=itsrs, llj=llj,
                numlike=numlike,flag=flag)
  }
  

  
  return (output)
  
}

