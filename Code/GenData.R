library(MASS)

GenData<-function(N, p, q, Beta, gamma, eta, distance, sigma, isActiveC1){
  
  ###### generate data for one simulation
  
  i=1; X=NULL; M=NULL; S=NULL; isActiveC=NULL; it=0; B=NULL; ix=0; im=0; im2=0
  isActiveC=NULL
  
  maxeigen=NULL
  
  repeat{
    
    it=it+1
    
    ##################################
    ####        generate X       #####
    ##################################

    
    xi=mvrnorm(1,mu=0.1+numeric(q),Sigma=diag(q))
    
    if ( any(abs(xi)>1) ) {ix=ix+1; next}
    
    # if ( any(xi>1) || any(xi<0)  ) {ix=ix+1; next}
    
    #### creat DTI network
    isActiveCi=isActiveC1
    
    if (xi[1]>0){
      isActiveCi[4,5]=1
      isActiveCi[5,4]=1
    }
    
    if (xi[3]>0){
      isActiveCi[3,5]=1
      isActiveCi[5,3]=1
    }
    
    
    
    
    #### count common neighbors
    Si=matrix(0,p,p)
    
    for (j in 1:(p-1)){
      for (k in (j+1):p){
        Si[j,k]=sum(isActiveCi[j,]*isActiveCi[k,])
        Si[k,j]=Si[j,k]
      }  ## end of k
    }  ## end of j
    
    
    degree=apply(isActiveCi,1,sum)
    degree
    
    
    numerator=matrix(0,p,p)
    prob=matrix(0,p,p)
    
    for (j in 1:(p-1)){
      for (k in (j+1):p){
        
        numerator[j,k]=((1+Si[j,k])^(gamma))/(distance[j,k]^(eta))
        numerator[k,j]=numerator[j,k]
        
        prob[j,k]=numerator[j,k]/(1+numerator[j,k])
        prob[k,j]=prob[j,k]
        
      } ## end of k
    }  ## end of j
    
    
    Bi=matrix(0,p,p)
    
    for (j in 1:(p-1)){
      for (k in (j+1):p){
        
        Bi[j,k]=sample(c(0,1),1,prob=c(1-prob[j,k],prob[j,k]))
        Bi[k,j]=Bi[j,k]
        
      } ## end of k
    } ## end of j
    
    
    ##################################
    ####        generate M       #####
    ##################################
    
    ## M
    
    Bi2=diag(x=NA,nrow=p,ncol=p)
    s=1
    for (j in 1:p){
      for (k in 1:p){
        if (k != j){
          # if (B[(i-1)*p+j,k] != 0){
          Bi2[j,k]=-Bi[j,k]*sum(xi*Beta[,s])
          # }
          s=s+1
        } ## end of k!=j
      } ## end of k
    } ## end of j
    diag(Bi2)=1/sigma
    
    maxeigen=c(maxeigen, max(eigen(Bi2)$values))
    
    if (min(eigen(Bi2)$values)<0.001) {im=im+1; next}
    
    mi=mvrnorm(1, mu=numeric(p), solve(Bi2))
    
    # if ( any(abs(mi)>40) ) {im2=im2+1; next}
    
    ### Record
    X=rbind(X,xi)
    M=rbind(M, mi)  
    
    
    bbi=NULL
    s=1
    for (j in 1:p){
      for (k in 1:p){
        if(k!=j){
          bbi[s]=Bi[j,k]
          s=s+1
        }
      } ### end of k
    } ### end of j
    
    
    B=rbind(B,bbi)
    
    
    S=rbind(S,Si)
    
    isActiveC=rbind(isActiveC,isActiveCi)
    
    if (i==N) break
    i=i+1
    
  }  ### end of repeat
  
  return(list(X=X, M=M, B=B, S=S, isActiveC=isActiveC,ix=ix, im=im, maxeigen=maxeigen))
  
  # , im2=im2
  
}







