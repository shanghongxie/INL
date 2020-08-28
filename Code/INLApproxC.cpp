// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::interfaces(r,cpp)]]


#include <RcppEigen.h>
#include <Rcpp.h>
// #include <cstdio>
// #include <iostream>
// #include <iomanip>
// #include <string>
// #include <map>
#include <random>
#include <math.h>



using namespace Rcpp;

// sth wrong with this function after n>13, can not use 'int' to define gamma, since there would be overflow problem,
// instead to use 'double';

/*****  Center and standardize  *****/
// [[Rcpp::export]]
List scaleC(Eigen::MatrixXd X){
  Eigen::VectorXd mX=X.colwise().mean(), sdX(X.cols()), sdXi;
  X.rowwise()-=mX.transpose();
  sdX=X.colwise().norm()/sqrt((double)X.rows());
  sdXi=1.0/sdX.array();
  X=X*sdXi.asDiagonal();
  return List::create(Named("x")=X, Named("sd")=sdX);
}


// /***  factorial  ****/
// [[Rcpp::export]]
double tgammaC(int n){
  if (n<=1) return 1;
  else
    return n*tgammaC(n-1);
}



/***** Calculate posteR
 *********  ior mean, normalization constant  */
// [[Rcpp::export]]

Eigen::MatrixXd posteriorC(Eigen::VectorXd c, Eigen::VectorXd d,int order, int q){
  
  int N=2*order+1;
  Eigen::MatrixXd m(N, q);
  
  int j,k,l;
  // // double sum;
  //  
  m.row(0)=Eigen::VectorXd::Ones(q);
  //    
  m(0,0)=exp(d(0));
  m(0,0)=m(0,0)/(1+m(0,0));
  
  for (j=1; j<N; ++j) {
    m(j,0)=m(j-1,0)*c(0);
  }
  
  // product of (1+e^{dk})
  double dd=1;
  for (k=0; k<q; ++k){
    dd=dd*(1+exp(d(k)));
  }
  
  double ek;
  for (j=1; j<N; ++j) {
    // j=13;
    // k=1;
    for (k=1; k<q; ++k){
      ek=exp(d(k));
      m(j,k)=dd*m(j,k-1)/(1+ek);
      for (l=0; l<=j; l++){
        //  m(j,k)=m(j,k)+ek/(1+ek)*m(l,k-1)*std::pow(c(k),j-l)*std::tgamma(j+1)/std::tgamma(j-l+1)/std::tgamma(l+1);
        m(j,k)=m(j,k)+dd*ek/(1+ek)*m(l,k-1)*std::pow(c(k),j-l)*tgammaC(j+1)/tgammaC(j-l+1)/tgammaC(l+1);
      }
    }
  }
  
  return (m);
}



// Monte Carol calculate E(exp(-Y^Tc)^2)  *** /

// [[Rcpp::export]]
Eigen::VectorXd monteC(Eigen::MatrixXd di, Eigen::MatrixXd ci, int replicate){
  
  
  
  int p=di.cols(), k, it, N=di.rows(), subject;
  double ek, probk;
  
  Eigen::MatrixXd edi=1+di.array().exp();
  
  // prod_k=1^q (1+exp(dk))
  Eigen::VectorXd ed=Eigen::VectorXd::Zero(N);
  
  ed=edi.rowwise().prod();
  
  // for (subject=0; subject<N; ++subject){
  //   ed(subject)=(edi.row(subject)).prod();
  // }
  
  Eigen::VectorXd expecti=Eigen::VectorXd::Zero(N);
  
  // constant C
  Eigen::VectorXd consti=Eigen::VectorXd::Zero(N);
  
  
  it=replicate;
  
  for (subject=0; subject<N; ++subject){
    Eigen::MatrixXd Y=Eigen::MatrixXd::Zero(replicate,p);
    Eigen::VectorXd Zi;
    Eigen::VectorXd Zi2;
    
    for (k=0; k<p; ++k){
      ek=exp(di(subject,k));
      probk=ek/(1+ek);
      
      std::random_device rd{}; // use to seed the rng 
      std::mt19937 rng{rd()}; // rng
      
      std::bernoulli_distribution d(probk);
      
      // generate 5000 runs
      for(std::size_t i = 0; i < replicate; ++i){
        Y(i,k)=d(rng);
      }
    }
    
    Zi=Y*(ci.row(subject)).transpose();
    
    Zi2=(Zi.array()).square();
    
    expecti(subject)=((-Zi2.array()).exp()).sum()/(Y.rows());
    
    std::setprecision(6);
    
    consti(subject)=expecti(subject)*ed(subject);
  }
  
  return consti;
}





// [[Rcpp::export]]
Eigen::VectorXd linearApprox2C(Eigen::MatrixXd di, Eigen::MatrixXd ci,  Eigen::VectorXd lambda1, Eigen::VectorXd lambda2, 
                               Eigen::VectorXd coeff1, Eigen::VectorXd coeff2){
  
  int p=ci.cols(), N=ci.rows();
  int i,j, k;
  int m=1;
  
  // I(cj<0)
  // Eigen::MatrixXd indicatorc=Eigen::MatrixXd::Zero(N,p);
  
  // dj*sign(cj)
  Eigen::MatrixXd disignc=di;
  
  Eigen::MatrixXd absc=ci;
  
  absc=ci.cwiseAbs();
  
  // Eigen::MatrixXd disignc1=di1;
  
  // cj-=max(-cj,0)
  Eigen::MatrixXd nci=Eigen::MatrixXd::Zero(N,p);
  //cj*I(cj>0)
  Eigen::MatrixXd pci=Eigen::MatrixXd::Zero(N,p);
  
  // djI(cj<0)
  Eigen::MatrixXd dci=Eigen::MatrixXd::Zero(N,p);
  
  Eigen::VectorXd a=Eigen::VectorXd::Zero(N);
  Eigen::VectorXd ap=Eigen::VectorXd::Zero(N);
  
  Eigen::VectorXd b=Eigen::VectorXd::Zero(N);
  
  for (i=0; i<N; ++i){
    for (j=0; j<p; ++j){
      if (ci(i,j)<0){
        nci(i,j)=-ci(i,j);
      } // end of if (ci(i,j)<0)
      else{
        pci(i,j)=ci(i,j);
      }
    } // end of j
  } // end of i
  
  
  // a=sum_j=1^p cj-
  a=nci.rowwise().sum();
  
  // ap=sum_j=1^p cj*I(cj>0)
  ap=pci.rowwise().sum();
  
  
  Eigen::VectorXd constant=Eigen::VectorXd::Zero(N);
  
  Eigen::VectorXd sump=Eigen::VectorXd::Zero(N);
  
  Eigen::VectorXd coeffvalue;
  Eigen::VectorXd lambdai=lambda1;
  
  //  i=0;
  for (i=0; i<N; ++i){
    coeffvalue=coeff1;
    // if sum_j cj*I(cj>0) or \sum_j cj- >1.5, use the coefficients for range (-4, 4)
    if (a(i)>2 | ap(i)>2){
      coeffvalue=coeff2;
      lambdai=lambda2;
    }
    
    m=lambdai.size();
    
    Eigen::VectorXd prodj;
    Eigen::VectorXd prod=Eigen::VectorXd::Ones(m);
    for (k=0; k<m; ++k){
      Eigen::VectorXd value=di.row(i)-lambdai(k)*(ci.row(i));
      // 1+exp(dj-lambda(k)*cj)
      prodj=Eigen::VectorXd::Ones(p);
      prodj=prodj.array()+value.array().exp();
      double prodk=prodj.array().prod();
      prod[k]=prodk;
    } // end of k
    //  coeffvalue=coeffvalue.array()*((lambdai*a(i)).exp()).array();
    sump(i)=coeffvalue.dot(prod);
    
  }
  
  
  //Eigen::VectorXd front=(-a.array()*a.array()+(dci.array()).rowwise().sum()-sumdci1.array()).exp();
  
  // Eigen::VectorXd front=((dci.array()).rowwise().sum()).exp();
  
  // Eigen::VectorXd front=(-a.array()*a.array()+(dci.array()).rowwise().sum()+ba.array()*ba.array()/4).exp();
  //  constant=front.array()*sump.array();
  // constant=front.array();
  //  constant=front;
  
  constant=sump;
  
  return constant;
  
  //  return prod;
  
  
  // return(List::create(Named("front")=front, Named("sump")=sump,Named("a")=a,Named("sumC")=sumC,Named("ba")=ba,
  //                           Named("disignc")=disignc,Named("constant")=constant,Named("pexp")=pexp)); 
  
  
}



// [[Rcpp::export]]
Eigen::VectorXd linearApproxC(Eigen::MatrixXd di, Eigen::MatrixXd ci,  Eigen::VectorXd lambda1,
                              Eigen::VectorXd coeff1){
  
  int p=ci.cols(), N=ci.rows();
  int i,j, k;
  int m1=lambda1.size();
  //, m2=lambda2.size();
  
  // I(cj<0)
  // Eigen::MatrixXd indicatorc=Eigen::MatrixXd::Zero(N,p);
  
  // dj*sign(cj)
  // Eigen::MatrixXd disignc=di;
  // 
  // Eigen::MatrixXd absc=ci;
  // 
  // absc=ci.cwiseAbs();
  // 
  // // Eigen::MatrixXd disignc1=di1;
  // 
  // // cj-=max(-cj,0)
  // Eigen::MatrixXd nci=Eigen::MatrixXd::Zero(N,p);
  // //cj*I(cj>0)
  // Eigen::MatrixXd pci=Eigen::MatrixXd::Zero(N,p);
  // 
  // // djI(cj<0)
  // Eigen::MatrixXd dci=Eigen::MatrixXd::Zero(N,p);
  // 
  // Eigen::VectorXd a=Eigen::VectorXd::Zero(N);
  // Eigen::VectorXd ap=Eigen::VectorXd::Zero(N);
  // 
  // // Eigen::VectorXd b=Eigen::VectorXd::Zero(N);
  // 
  // // |dj/cj|
  // // Eigen::VectorXd d_c=Eigen::VectorXd::Zero(N);
  // 
  // //  Eigen::VectorXd maxdisignc=Eigen::VectorXd::Zero(N);
  // 
  // for (i=0; i<N; ++i){
  //   for (j=0; j<p; ++j){
  //     if (ci(i,j)<0){
  //       // dj*sign(cj)
  //       disignc(i,j)=-di(i,j);
  //       
  //       nci(i,j)=-ci(i,j);
  //       
  //       // dj*I(cj<0)
  //       dci(i,j)=di(i,j);
  //       
  //     } // end of if (ci(i,j)<0)
  //     else{
  //       pci(i,j)=ci(i,j);
  //     }
  //     
  //     //    d_c(i,j)=di(i,j)/ci(i,j);
  //     
  //   } // end of j
  //   
  //   // b=(K+1)max(|dj/cj|), K=1
  //   //  b(i)=2*((d_c.row(i)).cwiseAbs()).maxCoeff();
  // } // end of i
  // 
  // 
  // // a=sum_j=1^p cj-
  // a=nci.rowwise().sum();
  // 
  // // ap=sum_j=1^p cj*I(cj>0)
  // ap=pci.rowwise().sum();
  
  
  Eigen::VectorXd constant=Eigen::VectorXd::Zero(N);
  
  //  Eigen::VectorXd sump=Eigen::VectorXd::Zero(N);
  
  // Eigen::VectorXd coeffvalue=coeff1;
  // Eigen::VectorXd lambdai=lambda1;
  
  Eigen::MatrixXd value=Eigen::MatrixXd::Ones(N,p);
  Eigen::MatrixXd prodj=Eigen::MatrixXd::Ones(N,p);
  Eigen::MatrixXd prods=Eigen::MatrixXd::Zero(N,m1);
  
  
  Eigen::VectorXd prodk=Eigen::VectorXd::Zero(N);
  
  for (k=0; k<m1; ++k){
    
    value=di-lambda1(k)*ci;
    
    prodj=Eigen::MatrixXd::Ones(N,p);
    // 1+exp(djsign(cj)-lambda(k)|cj|)
    prodj=prodj.array()+value.array().exp();
    prodk=prodj.rowwise().prod();
    // exp(lambda(k)*a)*\prod_j=1^p (1+exp(djsign(cj)-lambda(k)*|cj|))
    //   prodk=prodk.array()*(lambda1(k)*a.array()).exp();
    prods.col(k)=prodk;
  } // end of k
  
  constant=prods*coeff1;
  
  // for (i=0; i<N; ++i){
  //     coeffvalue=coeff1;
  // //  
  //    // if sum_j cj*I(cj>0) or \sum_j cj- >1.5, use the coefficients for range (-4, 4)
  //     if (a(i)>1.5 | ap(i)>1.5){
  //          coeffvalue=coeff2;
  //   //   // lambdai=lambda2;
  //   //   
  //   //  // m1=lambdai.size();
  //   //   
  //   //   prodj=Eigen::VectorXd::Zero(p);
  //   //   prods=Eigen::VectorXd::Zero(m2);
  //   //   
  //   //   double prodk=0.0;
  //   //   
  //   //   for (k=0; k<m2; ++k){
  //   //     
  //   //     value=disignc.row(i)-lambda2(k)*(absc.row(i));
  //   //     
  //   //     prodj=Eigen::VectorXd::Ones(p);
  //   //     // 1+exp(djsign(cj)-lambda(k)|cj|)
  //   //     prodj=prodj.array()+value.array().exp();
  //   //     prodk=prodj.array().prod();
  //   //     // exp(lambda(k)*a)*\prod_j=1^p (1+exp(djsign(cj)-lambda(k)*|cj|))
  //   //     prodk=prodk*exp(lambda2(k)*a(i));
  //   //     prods[k]=prodk;
  //   //   } // end of k
  //   //   //  coeffvalue=coeffvalue.array()*((lambdai*a(i)).exp()).array();
  //   //   sump(i)=coeff2.dot(prods);
  //   //   
  //   // } // end of a(i)
  //   
  //   
  //   prodj=Eigen::VectorXd::Zero(p);
  //   prods=Eigen::VectorXd::Zero(m1);
  //   
  //   double prodk=0.0;
  //   
  //   for (k=0; k<m1; ++k){
  //     
  //     value=disignc.row(i)-lambda1(k)*(absc.row(i));
  //     
  //     prodj=Eigen::VectorXd::Ones(p);
  //     // 1+exp(djsign(cj)-lambda(k)|cj|)
  //     prodj=prodj.array()+value.array().exp();
  //     prodk=prodj.array().prod();
  //     // exp(lambda(k)*a)*\prod_j=1^p (1+exp(djsign(cj)-lambda(k)*|cj|))
  //     prodk=prodk*exp(lambda1(k)*a(i));
  //     prods[k]=prodk;
  //   } // end of k
  // //  coeffvalue=coeffvalue.array()*((lambdai*a(i)).exp()).array();
  //   sump(i)=coeff1.dot(prods);
  //   
  // } // end of i
  
  
  //Eigen::VectorXd front=(-a.array()*a.array()+(dci.array()).rowwise().sum()-sumdci1.array()).exp();
  
  //  Eigen::VectorXd front=((dci.array()).rowwise().sum()).exp();
  
  // Eigen::VectorXd front=(-a.array()*a.array()+(dci.array()).rowwise().sum()+ba.array()*ba.array()/4).exp();
  //  constant=front.array()*sump.array();
  // constant=front.array();
  //  constant=front;
  
  //  constant=sump;
  
  return constant;
  
  
  // return(List::create(Named("front")=front, Named("sump")=sump,Named("a")=a,Named("sumC")=sumC,Named("ba")=ba,
  //                           Named("disignc")=disignc,Named("constant")=constant,Named("pexp")=pexp)); 
  
  
}


/********* E-step, expectation ********/
// [[Rcpp::export]]
List indexjC(int p){
  
  // location of (j,k) in jk vector
  Eigen::MatrixXd indexjk=Eigen::MatrixXd::Zero(p,p);
  
  // location of (j,k,k2) 
  Eigen::MatrixXd indexjkk=Eigen::MatrixXd::Zero(p*p,p);
  
  int j,k, k2, jk, jkk;
  
  //  XB=X*Beta;
  
  // indexAct non-zero vector beta_jk
  jk=0;
  jkk=0;
  for (j=0; j<p; ++j) {
    for (k=0; k<p; ++k) {
      if (k != j) {
        // if (Beta.col(jk).squaredNorm()>0) {
        //   indexAct(j,k)=jk+1;
        // }
        indexjk(j,k)=jk+1;
        ++jk;
        
        if (k<(p-1)){
          for (k2=k+1;k2<p;++k2){
            if (k2 != j){
              indexjkk(j*p+k,k2)=jkk+1;
              ++jkk;
            }
          } // end of k2
          
        } // end of k<p-1 
      } // end of k!=j
    } // end of k
  } // end of j
  
  return(List::create(Named("indexjk")=indexjk,
                      Named("indexjkk")=indexjkk)); 
}


/********* E-step, expectation ********/
// [[Rcpp::export]]
List Estep1jC(Eigen::MatrixXd Beta, Eigen::MatrixXd X, Eigen::MatrixXd M, Eigen::MatrixXd indexjk, Eigen::MatrixXd indexjkk,
              double thresh, int nodej){
  
  // Beta: q*p(p-1), ignore beta_jj=0; prob: N*p(p-1) matrix, symmetric, obtained from previous step, using gamma(t-1),eta(t-1)
  int p=M.cols(), N=M.rows(), q=X.cols(), po=p*(p-1), pp=p*(p-1)*(p-2)/2;
  int i, j, k, k2, l, jk, jkk, jk2;
  //  Eigen::MatrixXd XB=Eigen::MatrixXd::Zero(N, p*p);
  
  // // non-zero vector beta_jk
  // Eigen::MatrixXd indexAct=Eigen::MatrixXd::Zero(p,p);
  // 
  // // location of (j,k) in jk vector
  // Eigen::MatrixXd indexjk=Eigen::MatrixXd::Zero(p,p);
  // 
  // // location of (j,k,k2) 
  // Eigen::MatrixXd indexjkk=Eigen::MatrixXd::Zero(p*p,p);
  
  //  XB=X*Beta;
  
  // // indexAct non-zero vector beta_jk
  // jk=0;
  // jkk=0;
  // for (j=0; j<p; ++j) {
  //   for (k=0; k<p; ++k) {
  //     if (k != j) {
  //       // if (Beta.col(jk).squaredNorm()>0) {
  //       //   indexAct(j,k)=jk+1;
  //       // }
  //       indexjk(j,k)=jk+1;
  //       ++jk;
  //       
  //       if (k<(p-1)){
  //         for (k2=k+1;k2<p;++k2){
  //           if (k2 != j){
  //             indexjkk(j*p+k,k2)=jkk+1;
  //             ++jkk;
  //           }
  //         } // end of k2
  //         
  //       } // end of k<p-1 
  //     } // end of k!=j
  //   } // end of k
  // } // end of j
  
  
  Eigen::MatrixXd XBM=Eigen::MatrixXd::Zero(N, po);
  // // Mij-b(beta_jk^T*Xi)Mik
  Eigen::MatrixXd UMi=Eigen::MatrixXd::Zero(N, po);
  // 
  // 
  // // Mij-b(beta_jk^T*Xi)Mik-b'(beta_jk'^T*Xi)Mik'
  Eigen::MatrixXd UMMi=Eigen::MatrixXd::Zero(N, pp);
  // 
  
  j=nodej;
  //  jkk=0;
  // for (j=0; j<p; ++j){
  for (k=0; k<p; ++k){
    if (k!=j){
      jk=indexjk(j,k)-1;
      //  if (indexAct(j,k)>0){
      XBM.col(jk)=(X*Beta.col(jk)).array()*(M.col(k).array());
      UMi.col(jk)=M.col(j).array()-XBM.col(jk).array();
      
      if (k<(p-1)){
        for (k2=k+1; k2<p; ++k2){
          if (k2!=j){
            // cl=0, then p(Bijk,Bijk',Mi)=0
            
            jk2=indexjk(j,k2)-1;
            jkk=indexjkk(j*p+k,k2)-1;
            XBM.col(jk2)=(X*Beta.col(jk2)).array()*(M.col(k2).array());
            UMMi.col(jkk)=UMi.col(jk).array()-XBM.col(jk2).array();
          } 
        } // end of k2
      } // end of k<(p-1)
      
    } // if k != j
  } // end of k
  //  } // end of j
  
  return(List::create(Named("XBM")=XBM,Named("UMi")=UMi,Named("UMMi")=UMMi)); 
}


/********* E-step, expectation ********/
// [[Rcpp::export]]
Eigen::MatrixXd probC(int N, int p, Eigen::MatrixXd S, Eigen::MatrixXd distance,  double gamma, double eta){  
  
  int po=p*(p-1);
  int i, j, k, s;
  Eigen::MatrixXd numerator=Eigen::MatrixXd::Zero(N*p,p);
  Eigen::MatrixXd prob0=Eigen::MatrixXd::Zero(N*p,p);
  Eigen::MatrixXd prob=Eigen::MatrixXd::Zero(N,po);
  
  // log of product of degrees: sijk=sij*sik / common neighbors in DTI
  Eigen::MatrixXd logS=Eigen::MatrixXd::Zero(N*p, p);
  logS=log(1+S.array());
  
  // log of distance 
  Eigen::MatrixXd logdist=Eigen::MatrixXd::Zero(p, p);
  logdist=log(distance.array());
  logdist.diagonal()=Eigen::VectorXd::Zero(p);
  
  
  for (i=0; i<N; ++i){
    
    for (j=0; j<(p-1); ++j){
      for (k=j+1; k<p; ++k){
        // for (k=0; k<p; ++k){
        //   if (k != j){
        // jk=indexjk(j,k)-1;
        numerator(i*p+j,k)=exp(gamma*logS(i*p+j,k)-eta*logdist(j,k));
        numerator(i*p+k,j)=numerator(i*p+j,k);
        
        prob0(i*p+j,k)=numerator(i*p+j,k)/(1+numerator(i*p+j,k));
        prob0(i*p+k,j)=prob0(i*p+j,k);
        
        //   } // end of k!=j
        
      } // end of k
    }  // end of j
    
  } // end of i
  
  
  for (i=0; i<N; ++i){
    
    s=0;
    for (j=0; j<p; ++j){
      for (k=0; k<p; ++k){
        if(k!=j){
          prob(i,s)=prob0(i*p+j,k);
          ++s;
        }
      }
    }
    
  }  // subject i
  
  return prob;
  
}


/********* E-step, expectation ********/
// [[Rcpp::export]]
List ExpectjC(Eigen::MatrixXd Beta, Eigen::MatrixXd X, Eigen::MatrixXd M,
              Eigen::MatrixXd indexAct, Eigen::MatrixXd indexjk, Eigen::MatrixXd indexjkk, 
              Eigen::MatrixXd XBM, Eigen::MatrixXd UMi, Eigen::MatrixXd UMMi, Eigen::MatrixXd B2, Eigen::MatrixXd B3,
              Eigen::VectorXd sigma,  Eigen::MatrixXd prob, int nodej){   
  
  //  int j, int k, int k2
  int N=X.rows(), p=M.cols(), pp=p*(p-1)*(p-2)/2, po=p*(p-1);
  int l,s,l1;
  int i,j,k,k2;
  int jl, jl1, jk, kj, jk2, jkk, kjk,countb, count, countk;  
  //  int m=lambda.size();
  int p2=p-3;
  int row2, row3;
  int rownum2=B2.rows(), rownum3=B3.rows();
  
  
  Eigen::VectorXd d1;
  Eigen::VectorXd d2;
  Eigen::VectorXd c1;
  Eigen::VectorXd c2;
  
  // b=1
  Eigen::VectorXd llk1;
  // b=0
  Eigen::VectorXd llk2;
  
  // b=1
  Eigen::VectorXd lll1;
  // b=0
  Eigen::VectorXd lll2, lll3, lll4;
  
  // Eigen::MatrixXd numerator=Eigen::MatrixXd::Zero(N*p,p);
  // Eigen::MatrixXd prob0=Eigen::MatrixXd::Zero(N*p,p);
  // Eigen::MatrixXd prob=Eigen::MatrixXd::Zero(N,po);
  // 
  // // log of product of degrees: sijk=sij*sik / common neighbors in DTI
  // Eigen::MatrixXd logS=Eigen::MatrixXd::Zero(N*p, p);
  // logS=log(1+S.array());
  // 
  // // log of distance 
  // Eigen::MatrixXd logdist=Eigen::MatrixXd::Zero(p, p);
  // logdist=log(distance.array());
  // logdist.diagonal()=Eigen::VectorXd::Zero(p);
  // 
  // 
  // for (i=0; i<N; ++i){
  //   
  //   for (j=0; j<(p-1); ++j){
  //     for (k=j+1; k<p; ++k){
  //       // for (k=0; k<p; ++k){
  //       //   if (k != j){
  //       // jk=indexjk(j,k)-1;
  //       numerator(i*p+j,k)=exp(gamma*logS(i*p+j,k)-eta*logdist(j,k));
  //       numerator(i*p+k,j)=numerator(i*p+j,k);
  //       
  //       prob0(i*p+j,k)=numerator(i*p+j,k)/(1+numerator(i*p+j,k));
  //       prob0(i*p+k,j)=prob0(i*p+j,k);
  //       
  //       //   } // end of k!=j
  //       
  //     } // end of k
  //   }  // end of j
  //   
  // } // end of i
  // 
  // 
  // for (i=0; i<N; ++i){
  //   
  //   s=0;
  //   for (j=0; j<p; ++j){
  //     for (k=0; k<p; ++k){
  //       if(k!=j){
  //         prob(i,s)=prob0(i*p+j,k);
  //         ++s;
  //       }
  //     }
  //   }
  //   
  // }  // subject i
  
  
  
  
  // P(Bijk=b, Bijk'=b', Mi) b=1, b'=1
  Eigen::MatrixXd PBB1=Eigen::MatrixXd::Zero(N, pp);
  // P(Bijk=b, Bijk'=b', Mi) b=1, b'=0
  Eigen::MatrixXd PBB2=Eigen::MatrixXd::Zero(N, pp);
  // P(Bijk=b, Bijk'=b', Mi) b=0, b'=1
  Eigen::MatrixXd PBB3=Eigen::MatrixXd::Zero(N, pp);
  // P(Bijk=b, Bijk'=b', Mi) b=0, b'=0
  Eigen::MatrixXd PBB4=Eigen::MatrixXd::Zero(N, pp);
  
  // P(Bijk=1,  Mi) 
  Eigen::MatrixXd PB1=Eigen::MatrixXd::Zero(N, po);
  // P(Bijk=0, Mi) 
  Eigen::MatrixXd PB2=Eigen::MatrixXd::Zero(N, po);
  
  // E(BijkBijk'|Mi)
  Eigen::MatrixXd EBB=Eigen::MatrixXd::Zero(N, pp);
  
  // E(Bijk|Mi)
  Eigen::MatrixXd EB=Eigen::MatrixXd::Zero(N, po);
  
  
  j=nodej;
  
  // j=0; k=2; k2=3;
  //  for (j=0; j<p; ++j){
  // prod_{k!=j}(1-pijk)
  Eigen::VectorXd prodprobj=Eigen::VectorXd::Ones(N);
  
  for (l=0; l<p; ++l){
    if (l!=j){
      jl=indexjk(j,l)-1;
      prodprobj=(prodprobj.array())*(1-prob.col(jl).array());
    } // end of l!=j
  } // end of l
  
  
  for (k=0; k<p; ++k){
    if (k!=j){
      
      jk=indexjk(j,k)-1;
      
      // calculate E(Bijk)
      
      // b=1
      Eigen::MatrixXd db1=Eigen::MatrixXd::Zero(N, p-2);
      // b=0
      //   Eigen::MatrixXd db2=Eigen::MatrixXd::Zero(N, p-2);
      
      Eigen::MatrixXd cb1=Eigen::MatrixXd::Zero(N, p-2);
      
      //   Eigen::MatrixXd cb2=Eigen::MatrixXd::Zero(N, p-2);
      
      
      d1=Eigen::VectorXd::Zero(N);
      d2=Eigen::VectorXd::Zero(N);
      c1=Eigen::VectorXd::Zero(N);
      c2=Eigen::VectorXd::Zero(N);
      
      llk1=Eigen::VectorXd::Zero(N);
      llk2=Eigen::VectorXd::Zero(N);
      
      for (row2=0; row2<rownum2; ++row2){
        
        countb=0;
        for (l1=0; l1<p; ++l1){
          if ((l1!=j) & (l1!=k)){
            jl1=indexjk(j,l1)-1;
            
            // compute P(Bijk=b,Mi)
            
            // Bijl*log(pijl/(1-pijl))
            db1.col(countb)=B2(row2,countb)*log(prob.col(jl1).array()/(1-prob.col(jl1).array()));
            
            // Bijl*betajl*Xi*Mil
            cb1.col(countb)=B2(row2,countb)*XBM.col(jl1).array();
            
            ++countb;
          } // end of (l1!=j) & (l1!=k)
        } // end of l1
        
        // b=1: log(pijk/(1-pijk))+sum_{l!=j,k}Bijl*log(pijl/(1-pijl))
        d1=db1.rowwise().sum().array()+log(prob.col(jk).array()/(1-prob.col(jk).array()));
        // b=0
        d2=db1.rowwise().sum();
        
        // b=1: betajk*Xi*Mik+sum_{l!=j,k}Bijl*betajl*Xi*Mil
        c1=cb1.rowwise().sum().array()+XBM.col(jk).array();
        
        // b=0; sum_{l!=j,k}Bijl*betajl*Xi*Mil
        c2=cb1.rowwise().sum();
        
        // sum over Bijl in {0,1}, b=1
        llk1=llk1.array()+exp(-(M.col(j).array()-c1.array())*(M.col(j).array()-c1.array())/sigma(j)*0.5)*exp(d1.array());
        
        // sum over Bijl in {0,1}, b=0
        llk2=llk2.array()+exp(-(M.col(j).array()-c2.array())*(M.col(j).array()-c2.array())/sigma(j)*0.5)*exp(d2.array());
        
        
      }  // end of row2
      
      llk1=(prodprobj.array())*llk1.array()/sqrt(2*M_PI*sigma(j));
      
      llk2=(prodprobj.array())*llk2.array()/sqrt(2*M_PI*sigma(j));
      
      // P(Bijk=1, Mi) 
      PB1.col(jk)=llk1;
      
      // P(Bijk=0, Mi) 
      PB2.col(jk)=llk2;
      
      // expectation E(Bijk)
      EB.col(jk)=PB1.col(jk).array()/(PB1.col(jk).array()+PB2.col(jk).array());
      
      
      ////// start computing EBB
      
      if (k<(p-1)){
        for (k2=k+1; k2<p; ++k2){
          if (k2!=j){
            
            // b=1, b'=1
            Eigen::MatrixXd ddb1=Eigen::MatrixXd::Zero(N, p2);
            
            Eigen::MatrixXd ccb=Eigen::MatrixXd::Zero(N, p2);
            
            
            
            Eigen::VectorXd dd1=Eigen::VectorXd::Zero(N, p2);
            // b=1,b'=0
            Eigen::VectorXd dd2=Eigen::VectorXd::Zero(N, p2);
            // b=0,b'=1
            Eigen::VectorXd dd3=Eigen::VectorXd::Zero(N, p2);
            //b=0, b'=0
            Eigen::VectorXd dd4=Eigen::VectorXd::Zero(N, p2);
            
            // b=1, b'=1
            Eigen::VectorXd cc1=Eigen::VectorXd::Zero(N);
            // b=1, b'=0
            Eigen::VectorXd cc2=Eigen::VectorXd::Zero(N);
            // b=0,b'=1
            Eigen::VectorXd cc3=Eigen::VectorXd::Zero(N);
            //b=0, b'=0
            Eigen::VectorXd cc4=Eigen::VectorXd::Zero(N);
            
            lll1=Eigen::VectorXd::Zero(N);
            lll2=Eigen::VectorXd::Zero(N);
            lll3=Eigen::VectorXd::Zero(N);
            lll4=Eigen::VectorXd::Zero(N);
            
            
            jk2=indexjk(j,k2)-1;
            
            jkk=indexjkk(j*p+k,k2)-1;
            
            for (row3=0; row3<rownum3; ++row3){
              
              count=0;
              
              for (l=0; l<p; ++l){
                if ((l!=j) & (l!=k) & (l!=k2)){
                  jl=indexjk(j,l)-1;
                  
                  ddb1.col(count)=B3(row3,count)*log(prob.col(jl).array()/(1-prob.col(jl).array()));
                  
                  ccb.col(count)=B3(row3,count)*(XBM.col(jl).array());
                  
                  
                  ++count;
                } // end of (l!=j) & (l!=k) & (l!=k2)
              } // end of l
              
              // b=1, b'=1
              dd1=ddb1.rowwise().sum().array()+log(prob.col(jk).array()/(1-prob.col(jk).array()))+log(prob.col(jk2).array()/(1-prob.col(jk2).array()));
              
              // b=1, b'=0
              dd2=ddb1.rowwise().sum().array()+log(prob.col(jk).array()/(1-prob.col(jk).array()));
              
              // b=0, b'=1
              dd3=ddb1.rowwise().sum().array()+log(prob.col(jk2).array()/(1-prob.col(jk2).array()));
              
              // b=0, b'=0
              dd4=ddb1.rowwise().sum().array();
              
              
              // b=1, b'=1
              cc1=ccb.rowwise().sum().array()+XBM.col(jk).array()+XBM.col(jk2).array();
              // b=1, b'=0
              cc2=ccb.rowwise().sum().array()+XBM.col(jk).array();
              // b=0, b'=1
              cc3=ccb.rowwise().sum().array()+XBM.col(jk2).array();
              // b=0, b'=0
              cc4=ccb.rowwise().sum().array();
              
              // sum over Bijl in {0,1}, b=1, b'=1
              lll1=lll1.array()+exp(-(M.col(j).array()-cc1.array())*(M.col(j).array()-cc1.array())/sigma(j)*0.5)*exp(dd1.array());
              
              // sum over Bijl in {0,1}, b=1, b'=0
              lll2=lll2.array()+exp(-(M.col(j).array()-cc2.array())*(M.col(j).array()-cc2.array())/sigma(j)*0.5)*exp(dd2.array());
              
              // sum over Bijl in {0,1}, b=0, b'=1
              lll3=lll3.array()+exp(-(M.col(j).array()-cc3.array())*(M.col(j).array()-cc3.array())/sigma(j)*0.5)*exp(dd3.array());
              
              // sum over Bijl in {0,1}, b=0, b'=0
              lll4=lll4.array()+exp(-(M.col(j).array()-cc4.array())*(M.col(j).array()-cc4.array())/sigma(j)*0.5)*exp(dd4.array());
              
              
            } // end of row 3
            
            lll1=(prodprobj.array())*lll1.array()/sqrt(2*M_PI*sigma(j));
            lll2=(prodprobj.array())*lll2.array()/sqrt(2*M_PI*sigma(j));
            lll3=(prodprobj.array())*lll3.array()/sqrt(2*M_PI*sigma(j));
            lll4=(prodprobj.array())*lll4.array()/sqrt(2*M_PI*sigma(j));
            
            
            // P(Bijk=b, Bijk'=b', Mi) b=1, b'=1
            PBB1.col(jkk)=lll1;
            
            // P(Bijk=b, Bijk'=b', Mi) b=1, b'=0
            PBB2.col(jkk)=lll2;
            // P(Bijk=b, Bijk'=b', Mi) b=0, b'=1
            PBB3.col(jkk)=lll3;
            
            // P(Bijk=b, Bijk'=b', Mi) b=0, b'=0 
            PBB4.col(jkk)=lll4;
            
            // expectation E(BijkBijk')
            EBB.col(jkk)=PBB1.col(jkk).array()/(PBB1.col(jkk).array()+PBB2.col(jkk).array()+PBB3.col(jkk).array()+PBB4.col(jkk).array());
            
          } // end of  k2 != j
        } // end of k2
      } // end of if k<p-1
    } // end of k != j
  } // end of k
  
  
  //  } // end of j
  
  return(List::create( Named("prodprobj")=prodprobj,  Named("PBB1")=PBB1,Named("PBB2")=PBB2,Named("PBB3")=PBB3,Named("PBB4")=PBB4,
                       Named("PB1")=PB1, Named("PB2")=PB2, Named("EBB")=EBB, Named("EB")=EB
  ));
  
}


/********* E-step, expectation ********/
// [[Rcpp::export]]
Eigen::MatrixXd probjC(int N, int p, Eigen::MatrixXd S, Eigen::MatrixXd distance,  double gamma, double eta, int nodej){   
  
  // Bj: E(B_ijk) for node j only in previous step (posterior); XBM: XBM for node j only
  //  int j, int k, int k2
  int i,j,k, s;
  
  // only need calculate for node j
  Eigen::MatrixXd numerator=Eigen::MatrixXd::Zero(N,p-1);
  Eigen::MatrixXd prob=Eigen::MatrixXd::Zero(N,p-1);
  // Eigen::MatrixXd prob=Eigen::MatrixXd::Zero(N,po);
  
  // log of product of degrees: sijk=sij*sik / common neighbors in DTI
  Eigen::MatrixXd logS=Eigen::MatrixXd::Zero(N*p, p);
  logS=log(1+S.array());
  
  // log of distance 
  Eigen::MatrixXd logdist=Eigen::MatrixXd::Zero(p, p);
  logdist=log(distance.array());
  logdist.diagonal()=Eigen::VectorXd::Zero(p);
  
  j=nodej;
  
  s=0;
  for (k=0; k<p; ++k){
    if (k != j){
      // jk=indexjk(j,k)-1;
      
      for (i=0; i<N; ++i){
        
        numerator(i,s)=exp(gamma*logS(i*p+j,k)-eta*logdist(j,k));
        //  numerator(i*p+k,j)=numerator(i*p+j,k);
        
        prob(i,s)=numerator(i,s)/(1+numerator(i,s));
        //  prob0(i*p+k,j)=prob0(i*p+j,k);
        
      } // end of subject i
      
      s=s+1;
      
    } // end of k!=j
    
  } // end of k
  
  return prob;
  
}




/********* E-step, expectation ********/
// [[Rcpp::export]]
List ExpectMarginaljC2(Eigen::MatrixXd M, Eigen::MatrixXd prob,
                       Eigen::MatrixXd indexjk, Eigen::MatrixXd indexjkk, Eigen::MatrixXd EBji, Eigen::MatrixXd EBj1,
                       Eigen::MatrixXd XBM, Eigen::MatrixXd B1, Eigen::MatrixXd B2, Eigen::MatrixXd B3,
                       Eigen::MatrixXd B4, Eigen::MatrixXd B5, Eigen::MatrixXd B6, Eigen::MatrixXd B7, Eigen::MatrixXd B8, Eigen::MatrixXd B9,
                       Eigen::VectorXd sigma,  int nodej){   
  
  // Bj: E(B_ijk) for node j only in previous step (posterior); XBM: XBM for node j only
  //  int j, int k, int k2
  int N=M.rows(), p=M.cols(), pp=p*(p-1)*(p-2)/2, po=p*(p-1);
  int l,s,l1;
  int i,j,k,k2;
  int jl, jl1, jk, kj, jk2, jkk, kjk,countb, count, countk;  
  //  int m=lambda.size();
  int p2=p-3;
  int row2, row3;
  int rownum2=B8.rows(), rownum3=B3.rows();
  int lenBi,lenBi1, lenBBi, lenBBi1;
  int ind,ind1, indb, indb1;
  int indk;
  int indk2;
  
  Eigen::VectorXd indexk=Eigen::VectorXd::Zero(p);
  
  j=nodej;
  
  s=0;
  for (k=0; k<p; ++k){
    if (k != j){
      indexk[k]=s;
      s=s+1;
    } // end of k!=j
  } // end of k
  
  Eigen::VectorXd probi;
  Eigen::VectorXd XBMi;
  
  // E(Bijk)>0.9
  Eigen::VectorXd probi1;
  Eigen::VectorXd XBMi1;
  
  Eigen::VectorXd probi2;
  Eigen::VectorXd XBMi2;
  
  // E(Bijk, Bijk')>0.9
  Eigen::VectorXd probi21;
  Eigen::VectorXd XBMi21;
  
  
  double d1;
  double d2;
  double c1;
  double c2;
  
  // b=1
  double llk1;
  // b=0
  double llk2;
  
  // b=1
  double lll1;
  // b=0
  double lll2, lll3, lll4;
  
  Eigen::MatrixXd BB;
  Eigen::MatrixXd BB2;
  
  Eigen::VectorXd db1;
  Eigen::VectorXd cb1;
  
  Eigen::VectorXd ddb1;
  Eigen::VectorXd ccb1;
  
  double dd1;
  // b=1,b'=0
  double dd2;
  // b=0,b'=1
  double dd3;
  //b=0, b'=0
  double dd4;
  
  // b=1, b'=1
  double cc1;
  // b=1, b'=0
  double cc2;
  // b=0,b'=1
  double cc3;
  //b=0, b'=0
  double cc4;
  
  
  // P(Bijk=b, Bijk'=b', Mi) b=1, b'=1
  Eigen::MatrixXd PBB1=Eigen::MatrixXd::Zero(N, pp);
  // P(Bijk=b, Bijk'=b', Mi) b=1, b'=0
  Eigen::MatrixXd PBB2=Eigen::MatrixXd::Zero(N, pp);
  // P(Bijk=b, Bijk'=b', Mi) b=0, b'=1
  Eigen::MatrixXd PBB3=Eigen::MatrixXd::Zero(N, pp);
  // P(Bijk=b, Bijk'=b', Mi) b=0, b'=0
  Eigen::MatrixXd PBB4=Eigen::MatrixXd::Zero(N, pp);
  
  // P(Bijk=1,  Mi) 
  Eigen::MatrixXd PB1=Eigen::MatrixXd::Zero(N, po);
  // P(Bijk=0, Mi) 
  Eigen::MatrixXd PB2=Eigen::MatrixXd::Zero(N, po);
  
  // E(BijkBijk'|Mi)
  Eigen::MatrixXd EBB=Eigen::MatrixXd::Zero(N, pp);
  
  // E(Bijk|Mi)
  Eigen::MatrixXd EB=Eigen::MatrixXd::Zero(N, po);
  
  
  j=nodej;
  Eigen::VectorXd prodprobj=Eigen::VectorXd::Ones(N);
  
  
  for (s=0; s<(p-1); ++s){
    prodprobj=(prodprobj.array())*(1-prob.col(s).array());
  } // end of s
  
  
  // calculate E(Bijk)
  
  
  
  // k=1;
  llk1=0.0;
  llk2=0.0;
  
  for (k=0; k<p; ++k){
    if (k!=j){
      
      indk=indexk[k];
      jk=indexjk(j,k)-1;
      
      BB=B8;
      
      for (i=0; i<N; ++i){
        
        // for each k, for summation
        probi=Eigen::VectorXd::Zero(p-2);
        XBMi=Eigen::VectorXd::Zero(p-2);
        
        
        // if E(Bijk|...)>0.9, let bijk=1; if E(Bijk|...)<0.1, let bijk=0; if 0.1<=E(Bijk|...)<=0.9, calculate in the summation
        ind=0;
        ind1=0;
        
        for (l1=0; l1<p; ++l1){
          if ((l1!=j) & (l1!=k)){
            jl1=indexk[l1];
            
            probi[ind]=prob(i,jl1);
            XBMi[ind]=XBM(i,jl1);
            ++ind;
            
          } // end of l1!=j, l1!=k, Bj(i,l1)
        } // end of l1
        
        
        db1=Eigen::VectorXd::Zero(p-2);
        cb1=Eigen::VectorXd::Zero(p-2);
        
        llk1=0.0;
        llk2=0.0;
        
        d1=0.0;
        d2=0.0;
        c1=0.0;
        c2=0.0;
        
        for (row2=0; row2<rownum2; ++row2){
          
          // countb=0;
          // for (l1=0; l1<p; ++l1){
          //   if ((l1!=j) & (l1!=k)){
          //     jl1=indexk[l1];
          
          //    db1(countb)=B8(row2,countb)*log(prob(i,jl1)/(1-prob(i,jl1)));
          
          
          
          db1=BB.row(row2).array();
          db1=db1.array()*log(probi.array()/(1-probi.array()));
          
          //     cb1(countb)=B8(row2,countb)*XBM(i,jl1);
          
          //  cb1=(BB.row(row2)).array()*(XBMi.array());
          
          cb1=BB.row(row2).array();
          cb1=cb1.array()*XBMi.array();
          
          //    ++countb;
          //  } // end of l1!=j, l1!=k
          // } // end of l1
          
          // Bijl*betajl*Xi*Mil
          //   cb1=BB.row(row2).array()*XBMi.array();
          
          // b=1: log(pijk/(1-pijk))+sum_{l!=j,k}Bijl*log(pijl/(1-pijl))
          d1=db1.sum()+log(prob(i,indk)/(1-prob(i,indk)));
          // b=1: betajk*Xi*Mik+sum_{l!=j,k}Bijl*betajl*Xi*Mil
          c1=cb1.sum()+XBM(i,indk);
          
          
          // b=0
          d2=db1.sum();
          // b=0; sum_{l!=j,k}Bijl*betajl*Xi*Mil
          c2=cb1.sum();
          
          // sum over Bijl in {0,1}, b=1
          llk1=llk1+exp(-(M(i,j)-c1)*(M(i,j)-c1)/sigma(j)*0.5)*exp(d1);
          
          // sum over Bijl in {0,1}, b=0
          llk2=llk2+exp(-(M(i,j)-c2)*(M(i,j)-c2)/sigma(j)*0.5)*exp(d2);
          
        } // end of row 2
        
        // llk1=(prodprobj[i])*llk1/sqrt(2*M_PI*sigma(j));
        // 
        // llk2=(prodprobj[i])*llk2/sqrt(2*M_PI*sigma(j));
        
        
        // P(Bijk=1, Mi) 
        PB1(i,jk)=(prodprobj[i])*llk1/sqrt(2*M_PI*sigma(j));
        
        // P(Bijk=0, Mi) 
        PB2(i,jk)=(prodprobj[i])*llk2/sqrt(2*M_PI*sigma(j));
        
        // expectation E(Bijk)
        EB(i,jk)=PB1(i,jk)/(PB1(i,jk)+PB2(i,jk));
        
        
        
      } // end of subject i
      
      
    } // end of k != j
    
  } // end of k
  
  
  
  
  
  //  } // end of j
  
  return(List::create(Named("prodprobj")=prodprobj,  Named("PBB1")=PBB1,Named("PBB2")=PBB2,Named("PBB3")=PBB3,Named("PBB4")=PBB4,
                      Named("PB1")=PB1, Named("PB2")=PB2, Named("EBB")=EBB, Named("EB")=EB, Named("indexk")=indexk, Named("probi")=probi, 
                            Named("lenBi")=lenBi, Named("lenBi1")=lenBi1, Named("BB")=BB, Named("ind")=ind, Named("db1")=db1, Named("cb1")=cb1, Named("BBi")=BB.row(row2).array(),
                                  Named("probi")=probi
  ));
  
}





/********* E-step, expectation ********/
// [[Rcpp::export]]
List ExpectMarginaljC(Eigen::MatrixXd M, Eigen::MatrixXd prob,
                      Eigen::MatrixXd indexjk, Eigen::MatrixXd indexjkk, Eigen::MatrixXd EBji, Eigen::MatrixXd EBj1,
                      Eigen::MatrixXd XBM, Eigen::MatrixXd B1, Eigen::MatrixXd B2, Eigen::MatrixXd B3,
                      Eigen::MatrixXd B4, Eigen::MatrixXd B5, Eigen::MatrixXd B6, Eigen::MatrixXd B7, Eigen::MatrixXd B8, Eigen::MatrixXd B9,
                      Eigen::VectorXd sigma,  int nodej){   
  
  // Bj: E(B_ijk) for node j only in previous step (posterior); XBM: XBM for node j only
  //  int j, int k, int k2
  int N=M.rows(), p=M.cols(), pp=p*(p-1)*(p-2)/2, po=p*(p-1);
  int l,s,l1;
  int i,j,k,k2;
  int jl, jl1, jk, kj, jk2, jkk, kjk,countb, count, countk;  
  //  int m=lambda.size();
  int p2=p-3;
  int row2, row3;
  int rownum2, rownum3;
  int lenBi,lenBi1, lenBBi, lenBBi1;
  int ind,ind1, indb, indb1;
  int indk;
  int indk2;
  
  Eigen::VectorXd indexk=Eigen::VectorXd::Zero(p);
  
  j=nodej;
  
  s=0;
  for (k=0; k<p; ++k){
    if (k != j){
      indexk[k]=s;
      s=s+1;
    } // end of k!=j
    
  } // end of k
  
  
  Eigen::VectorXd probi;
  Eigen::VectorXd XBMi;
  
  // E(Bijk)>0.9
  Eigen::VectorXd probi1;
  Eigen::VectorXd XBMi1;
  
  Eigen::VectorXd probi2;
  Eigen::VectorXd XBMi2;
  
  
  // E(Bijk, Bijk')>0.9
  Eigen::VectorXd probi21;
  Eigen::VectorXd XBMi21;
  
  
  
  double d1;
  double d2;
  double c1;
  double c2;
  
  // b=1
  double llk1;
  // b=0
  double llk2;
  
  // b=1
  double lll1;
  // b=0
  double lll2, lll3, lll4;
  
  Eigen::MatrixXd BB;
  Eigen::MatrixXd BB2;
  
  Eigen::VectorXd db1;
  Eigen::VectorXd cb1;
  
  Eigen::VectorXd ddb1;
  Eigen::VectorXd ccb1;
  
  double dd1;
  // b=1,b'=0
  double dd2;
  // b=0,b'=1
  double dd3;
  //b=0, b'=0
  double dd4;
  
  // b=1, b'=1
  double cc1;
  // b=1, b'=0
  double cc2;
  // b=0,b'=1
  double cc3;
  //b=0, b'=0
  double cc4;
  
  
  // P(Bijk=b, Bijk'=b', Mi) b=1, b'=1
  Eigen::MatrixXd PBB1=Eigen::MatrixXd::Zero(N, pp);
  // P(Bijk=b, Bijk'=b', Mi) b=1, b'=0
  Eigen::MatrixXd PBB2=Eigen::MatrixXd::Zero(N, pp);
  // P(Bijk=b, Bijk'=b', Mi) b=0, b'=1
  Eigen::MatrixXd PBB3=Eigen::MatrixXd::Zero(N, pp);
  // P(Bijk=b, Bijk'=b', Mi) b=0, b'=0
  Eigen::MatrixXd PBB4=Eigen::MatrixXd::Zero(N, pp);
  
  // P(Bijk=1,  Mi) 
  Eigen::MatrixXd PB1=Eigen::MatrixXd::Zero(N, po);
  // P(Bijk=0, Mi) 
  Eigen::MatrixXd PB2=Eigen::MatrixXd::Zero(N, po);
  
  // E(BijkBijk'|Mi)
  Eigen::MatrixXd EBB=Eigen::MatrixXd::Zero(N, pp);
  
  // E(Bijk|Mi)
  Eigen::MatrixXd EB=Eigen::MatrixXd::Zero(N, po);
  
  
  j=nodej;
  Eigen::VectorXd prodprobj=Eigen::VectorXd::Ones(N);
  
  
  for (s=0; s<(p-1); ++s){
    prodprobj=(prodprobj.array())*(1-prob.col(s).array());
  } // end of s
  
  
  // calculate E(Bijk)
  
  for (i=0; i<N; ++i){
    
    // k=1;
    
    for (k=0; k<p; ++k){
      if (k!=j){
        
        indk=indexk[k];
        jk=indexjk(j,k)-1;
        
        // 0.1<=E(Bijk)<=0.9
        lenBi=EBji.row(i).sum()-EBji(i,indk);
        
        // E(Bijk)>0.9
        lenBi1=EBj1.row(i).sum()-EBj1(i,indk);
        
        if (lenBi==1){
          BB=B1;
        }
        if (lenBi==2){
          BB=B2;
        }
        if (lenBi==3){
          BB=B3;
        }
        if (lenBi==4){
          BB=B4;
        }
        if (lenBi==5){
          BB=B5;
        }
        if (lenBi==6){
          BB=B6;
        }
        if (lenBi==7){
          BB=B7;
        }
        if (lenBi==8){
          BB=B8;
        }
        
        
        // for each k, for summation
        probi=Eigen::VectorXd::Zero(lenBi);
        XBMi=Eigen::VectorXd::Zero(lenBi);
        
        // set bijk=1
        probi1=Eigen::VectorXd::Zero(lenBi1);
        XBMi1=Eigen::VectorXd::Zero(lenBi1);
        // 
        // // set bijk=0
        // probi0=Eigen::VectorXd::Zero(lenBi);
        // XBMi0=Eigen::VectorXd::Zero(lenBi);
        
        
        // if E(Bijk|...)>0.9, let bijk=1; if E(Bijk|...)<0.1, let bijk=0; if 0.1<=E(Bijk|...)<=0.9, calculate in the summation
        ind=0;
        ind1=0;
        
        for (l1=0; l1<p; ++l1){
          if ((l1!=j) & (l1!=k)){
            jl1=indexk[l1];
            
            // extract index needed in the summation (0.1<=E(Bijk)<=0.9)
            if (EBji(i,jl1)==1){
              probi[ind]=prob(i,jl1);
              XBMi[ind]=XBM(i,jl1);
              ++ind;
            } // end of EBji
            
            // E(Bijk)>0.9
            if (EBj1(i,jl1)==1){
              probi1[ind1]=prob(i,jl1);
              XBMi1[ind1]=XBM(i,jl1);
              ++ind1;
            }
            
            // // E(Bijk)<0.1
            // if (EBj0(i,jl1)==1){
            //   probi0[ind0]=prob(i,jl1);
            //   XBMi0[ind0]=XBM(i,jl1);
            //   ++ind0;
            // }
            
          } // end of l1!=j, l1!=k, Bj(i,l1)
        } // end of l1
        
        
        db1=Eigen::VectorXd::Zero(lenBi);
        cb1=Eigen::VectorXd::Zero(lenBi);
        
        llk1=0.0;
        llk2=0.0;
        
        d1=0.0;
        d2=0.0;
        c1=0.0;
        c2=0.0;
        
        
        if (lenBi==0){
          // if E(Bijk)>0.9, always bijk=1
          
          if (lenBi1>0){
            c1=XBM(i,indk)+XBMi1.sum();
            c2=XBMi1.sum();
            d1=log(prob(i,indk)/(1-prob(i,indk)))+(log(probi1.array()/(1-probi1.array()))).sum();
            d2=(log(probi1.array()/(1-probi1.array()))).sum();
            
          } else{
            c1=XBM(i,indk);
            c2=0;
            d1=log(prob(i,indk)/(1-prob(i,indk)));
            d2=0;
          }
          
          
          llk1=exp(-(M(i,j)-c1)*(M(i,j)-c1)/sigma(j)*0.5)*exp(d1);
          
          // sum over Bijl in {0,1}, b=0
          llk2=exp(-(M(i,j)-c2)*(M(i,j)-c2)/sigma(j)*0.5)*exp(d2);
          
          llk1=(prodprobj[i])*llk1/sqrt(2*M_PI*sigma(j));
          
          llk2=(prodprobj[i])*llk2/sqrt(2*M_PI*sigma(j));
          
        } else{
          
          rownum2=BB.rows();
          
          for (row2=0; row2<rownum2; ++row2){
            
            db1=BB.row(row2).array();
            db1=db1.array()*log(probi.array()/(1-probi.array()));
            
            // Bijl*betajl*Xi*Mil
            cb1=BB.row(row2).array();
            cb1=cb1.array()*XBMi.array();
            
            if (lenBi1>0){
              // b=1: log(pijk/(1-pijk))+sum_{l!=j,k}Bijl*log(pijl/(1-pijl))
              d1=db1.sum()+log(prob(i,indk)/(1-prob(i,indk)))+(log(probi1.array()/(1-probi1.array()))).sum();
              // b=1: betajk*Xi*Mik+sum_{l!=j,k}Bijl*betajl*Xi*Mil
              c1=cb1.sum()+XBM(i,indk)+XBMi1.sum();
              
              // b=0
              d2=db1.sum()+(log(probi1.array()/(1-probi1.array()))).sum();
              // b=0; sum_{l!=j,k}Bijl*betajl*Xi*Mil
              c2=cb1.sum()+XBMi1.sum();
              
            } else{
              // b=1: log(pijk/(1-pijk))+sum_{l!=j,k}Bijl*log(pijl/(1-pijl))
              d1=db1.sum()+log(prob(i,indk)/(1-prob(i,indk)));
              // b=1: betajk*Xi*Mik+sum_{l!=j,k}Bijl*betajl*Xi*Mil
              c1=cb1.sum()+XBM(i,indk);
              
              // b=0
              d2=db1.sum();
              // b=0; sum_{l!=j,k}Bijl*betajl*Xi*Mil
              c2=cb1.sum();
              
            }
            
            
            
            // sum over Bijl in {0,1}, b=1
            llk1=llk1+exp(-(M(i,j)-c1)*(M(i,j)-c1)/sigma(j)*0.5)*exp(d1);
            
            // sum over Bijl in {0,1}, b=0
            llk2=llk2+exp(-(M(i,j)-c2)*(M(i,j)-c2)/sigma(j)*0.5)*exp(d2);
          } // end of row 2
          
          llk1=(prodprobj[i])*llk1/sqrt(2*M_PI*sigma(j));
          
          llk2=(prodprobj[i])*llk2/sqrt(2*M_PI*sigma(j));
          
        }  // end of else
        
        
        // P(Bijk=1, Mi) 
        PB1(i,jk)=llk1;
        
        // P(Bijk=0, Mi) 
        PB2(i,jk)=llk2;
        
        // expectation E(Bijk)
        EB(i,jk)=PB1(i,jk)/(PB1(i,jk)+PB2(i,jk));
        
        
        
        
        //////////////////////////////////////////
        //////    start computing EBB     ////////
        //////////////////////////////////////////
        
        if (k<(p-1)){
          for (k2=k+1; k2<p; ++k2){
            if (k2!=j){
              
              indk2=indexk[k2];
              
              lenBBi=EBji.row(i).sum()-EBji(i,indk)-EBji(i,indk2);
              
              // E(Bijk)>0.9
              lenBBi1=EBj1.row(i).sum()-EBj1(i,indk)-EBj1(i,indk2);
              
              
              if (lenBBi==1){
                BB2=B1;
              }
              
              if (lenBBi==2){
                BB2=B2;
              }
              if (lenBBi==3){
                BB2=B3;
              }
              if (lenBBi==4){
                BB2=B4;
              }
              if (lenBBi==5){
                BB2=B5;
              }
              if (lenBBi==6){
                BB2=B6;
              }
              if (lenBBi==7){
                BB2=B7;
              }
              
              jk2=indexjk(j,k2)-1;
              jkk=indexjkk(j*p+k,k2)-1;
              
              probi2=Eigen::VectorXd::Zero(lenBBi);
              XBMi2=Eigen::VectorXd::Zero(lenBBi);
              
              // E(Bijk)>0.9
              probi21=Eigen::VectorXd::Zero(lenBBi1);
              XBMi21=Eigen::VectorXd::Zero(lenBBi1);
              
              indb=0;
              indb1=0;
              for (l=0; l<p; ++l){
                if ((l!=j) & (l!=k) & (l!=k2)){
                  jl=indexk[l];
                  
                  // extract index needed in the summation
                  if (EBji(i,jl)==1){
                    
                    probi2[indb]=prob(i,jl);
                    XBMi2[indb]=XBM(i,jl);
                    ++indb;
                  }  // end of EBji
                  
                  // E(Bijk)>0.9
                  if (EBj1(i,jl)==1){
                    
                    probi21[indb1]=prob(i,jl);
                    XBMi21[indb1]=XBM(i,jl);
                    ++indb1;
                  }  // end of EBji
                  
                  
                } // end of l!=j, l!=k, l!=k2
              } // end of l
              
              
              ddb1=Eigen::VectorXd::Zero(lenBBi);
              ccb1=Eigen::VectorXd::Zero(lenBBi);
              
              dd1=0.0;
              // b=1,b'=0
              dd2=0.0;
              // b=0,b'=1
              dd3=0.0;
              //b=0, b'=0
              dd4=0.0;
              
              // b=1, b'=1
              cc1=0.0;
              // b=1, b'=0
              cc2=0.0;
              // b=0,b'=1
              cc3=0.0;
              //b=0, b'=0
              cc4=0.0;
              
              lll1=0.0;
              lll2=0.0;
              lll3=0.0;
              lll4=0.0;
              
              if (lenBBi==0){
                
                if (lenBBi1>0){
                  // // b=1, b'=1
                  dd1=log(prob(i,indk)/(1-prob(i,indk)))+log(prob(i,indk2)/(1-prob(i,indk2)))+(log(probi21.array()/(1-probi21.array()))).sum();
                  
                  // b=1, b'=0
                  dd2=log(prob(i,indk)/(1-prob(i,indk)))+(log(probi21.array()/(1-probi21.array()))).sum();
                  // 
                  // b=0, b'=1
                  dd3=log(prob(i,indk2)/(1-prob(i,indk2)))+(log(probi21.array()/(1-probi21.array()))).sum();
                  
                  // b=0, b'=0
                  dd4=(log(probi21.array()/(1-probi21.array()))).sum();
                  
                  // b=1, b'=1
                  cc1=XBM(i,indk)+XBM(i,indk2)+XBMi21.sum();
                  // b=1, b'=0
                  cc2=XBM(i,indk)+XBMi21.sum();
                  // b=0, b'=1
                  cc3=XBM(i,indk2)+XBMi21.sum();
                  // b=0, b'=0
                  cc4=XBMi21.sum();
                  
                } else {
                  
                  // b=1, b'=1
                  dd1=log(prob(i,indk)/(1-prob(i,indk)))+log(prob(i,indk2)/(1-prob(i,indk2)));
                  
                  // b=1, b'=0
                  dd2=log(prob(i,indk)/(1-prob(i,indk)));
                  
                  // b=0, b'=1
                  dd3=log(prob(i,indk2)/(1-prob(i,indk2)));
                  
                  dd4=0;
                  
                  // b=1, b'=1
                  cc1=XBM(i,indk)+XBM(i,indk2);
                  // b=1, b'=0
                  cc2=XBM(i,indk);
                  // b=0, b'=1
                  cc3=XBM(i,indk2);
                  
                  cc4=0;
                }  // end of else
                
                lll1=exp(-(M(i,j)-cc1)*(M(i,j)-cc1)/sigma(j)*0.5)*exp(dd1);
                
                // sum over Bijl in {0,1}, b=1, b'=0
                lll2=exp(-(M(i,j)-cc2)*(M(i,j)-cc2)/sigma(j)*0.5)*exp(dd2);
                
                // sum over Bijl in {0,1}, b=0, b'=1
                lll3=exp(-(M(i,j)-cc3)*(M(i,j)-cc3)/sigma(j)*0.5)*exp(dd3);
                
                // sum over Bijl in {0,1}, b=0, b'=0
                lll4=exp(-(M(i,j)-cc4)*(M(i,j)-cc4)/sigma(j)*0.5)*exp(dd4);
                
                lll1=(prodprobj[i])*lll1/sqrt(2*M_PI*sigma(j));
                lll2=(prodprobj[i])*lll2/sqrt(2*M_PI*sigma(j));
                lll3=(prodprobj[i])*lll3/sqrt(2*M_PI*sigma(j));
                lll4=(prodprobj[i])*lll4/sqrt(2*M_PI*sigma(j));
                
                
              } else {
                
                rownum3=BB2.rows();
                
                for (row3=0; row3<rownum3; ++row3){
                  
                  ddb1=BB2.row(row3).array();
                  ddb1=ddb1.array()*log(probi2.array()/(1-probi2.array()));
                  
                  ccb1=BB2.row(row3).array();
                  ccb1=ccb1.array()*(XBMi2.array());
                  
                  
                  if (lenBBi1>0){
                    
                    // b=1, b'=1
                    dd1=ddb1.sum()+log(prob(i,indk)/(1-prob(i,indk)))+log(prob(i,indk2)/(1-prob(i,indk2)))+(log(probi21.array()/(1-probi21.array()))).sum();
                    
                    // b=1, b'=0
                    dd2=ddb1.sum()+log(prob(i,indk)/(1-prob(i,indk)))+(log(probi21.array()/(1-probi21.array()))).sum();
                    
                    // b=0, b'=1
                    dd3=ddb1.sum()+log(prob(i,indk2)/(1-prob(i,indk2)))+(log(probi21.array()/(1-probi21.array()))).sum();
                    
                    // b=0, b'=0
                    dd4=ddb1.sum()+(log(probi21.array()/(1-probi21.array()))).sum();
                    
                    
                    // b=1, b'=1
                    cc1=ccb1.sum()+XBM(i,indk)+XBM(i,indk2)+XBMi21.sum();
                    // b=1, b'=0
                    cc2=ccb1.sum()+XBM(i,indk)+XBMi21.sum();
                    // b=0, b'=1
                    cc3=ccb1.sum()+XBM(i,indk2)+XBMi21.sum();
                    // b=0, b'=0
                    cc4=ccb1.sum()+XBMi21.sum();
                    
                    
                  } else {
                    
                    // b=1, b'=1
                    dd1=ddb1.sum()+log(prob(i,indk)/(1-prob(i,indk)))+log(prob(i,indk2)/(1-prob(i,indk2)));
                    
                    // b=1, b'=0
                    dd2=ddb1.sum()+log(prob(i,indk)/(1-prob(i,indk)));
                    
                    // b=0, b'=1
                    dd3=ddb1.sum()+log(prob(i,indk2)/(1-prob(i,indk2)));
                    
                    // b=0, b'=0
                    dd4=ddb1.sum();
                    
                    
                    // b=1, b'=1
                    cc1=ccb1.sum()+XBM(i,indk)+XBM(i,indk2);
                    // b=1, b'=0
                    cc2=ccb1.sum()+XBM(i,indk);
                    // b=0, b'=1
                    cc3=ccb1.sum()+XBM(i,indk2);
                    // b=0, b'=0
                    cc4=ccb1.sum();
                    
                    
                  } // end of else
                  
                  
                  // sum over Bijl in {0,1}, b=1, b'=1
                  lll1=lll1+exp(-(M(i,j)-cc1)*(M(i,j)-cc1)/sigma(j)*0.5)*exp(dd1);
                  
                  // sum over Bijl in {0,1}, b=1, b'=0
                  lll2=lll2+exp(-(M(i,j)-cc2)*(M(i,j)-cc2)/sigma(j)*0.5)*exp(dd2);
                  
                  // sum over Bijl in {0,1}, b=0, b'=1
                  lll3=lll3+exp(-(M(i,j)-cc3)*(M(i,j)-cc3)/sigma(j)*0.5)*exp(dd3);
                  
                  // sum over Bijl in {0,1}, b=0, b'=0
                  lll4=lll4+exp(-(M(i,j)-cc4)*(M(i,j)-cc4)/sigma(j)*0.5)*exp(dd4);
                  
                  
                } // end of row 3
                
                lll1=(prodprobj[i])*lll1/sqrt(2*M_PI*sigma(j));
                lll2=(prodprobj[i])*lll2/sqrt(2*M_PI*sigma(j));
                lll3=(prodprobj[i])*lll3/sqrt(2*M_PI*sigma(j));
                lll4=(prodprobj[i])*lll4/sqrt(2*M_PI*sigma(j));
                
                
              }
              
              
              // P(Bijk=b, Bijk'=b', Mi) b=1, b'=1
              PBB1(i,jkk)=lll1;
              
              // P(Bijk=b, Bijk'=b', Mi) b=1, b'=0
              PBB2(i,jkk)=lll2;
              // P(Bijk=b, Bijk'=b', Mi) b=0, b'=1
              PBB3(i,jkk)=lll3;
              
              // P(Bijk=b, Bijk'=b', Mi) b=0, b'=0 
              PBB4(i,jkk)=lll4;
              
              // expectation E(BijkBijk')
              EBB(i,jkk)=PBB1(i,jkk)/(PBB1(i,jkk)+PBB2(i,jkk)+PBB3(i,jkk)+PBB4(i,jkk));
              
            } // end of  k2 != j
            
          } // end of k2
        } // end of if k<p-1
        
      } // end of k != j
      
    } // end of k
    
  } // end of subject i
  
  
  
  //  } // end of j
  
  return(List::create(Named("prob")=prob,  Named("PBB1")=PBB1,Named("PBB2")=PBB2,Named("PBB3")=PBB3,Named("PBB4")=PBB4,
                            Named("PB1")=PB1, Named("PB2")=PB2, Named("EBB")=EBB, Named("EB")=EB, Named("indexk")=indexk, Named("probi")=probi, 
                                  Named("lenBi")=lenBi, Named("lenBi1")=lenBi1
  ));
  
}





/********* E-step, expectation ********/
// [[Rcpp::export]]
List ExpectDTIjC(Eigen::MatrixXd M, Eigen::MatrixXd EBji, Eigen::MatrixXd EBj1,
                 Eigen::MatrixXd indexjk, Eigen::MatrixXd indexjkk, 
                 Eigen::MatrixXd XBM, Eigen::MatrixXd B1, Eigen::MatrixXd B2, Eigen::MatrixXd B3,
                 Eigen::MatrixXd B4, Eigen::MatrixXd B5, Eigen::MatrixXd B6, Eigen::MatrixXd B7, Eigen::MatrixXd B8, Eigen::MatrixXd B9,
                 Eigen::VectorXd sigma,  Eigen::MatrixXd S, Eigen::MatrixXd distance,  double gamma, double eta, int nodej){   
  
  // EBj: E(B_ijk) for node j only in previous step (posterior); XBM: XBM for node j only
  //  int j, int k, int k2
  int N=M.rows(), p=M.cols(), pp=p*(p-1)*(p-2)/2, po=p*(p-1);
  int l,s,l1;
  int i,j,k,k2;
  int jl, jl1, jk, kj, jk2, jkk, kjk,countb, count, countk;  
  //  int m=lambda.size();
  int p2=p-3;
  int row2, row3;
  int rownum2, rownum3;
  int lenBi,lenBi1, lenBBi, lenBBi1;
  int ind,ind1, indb, indb1;
  int indk;
  int indk2;
  
  Eigen::VectorXd indexk=Eigen::VectorXd::Zero(p);
  
  Eigen::VectorXd probi;
  Eigen::VectorXd XBMi;
  
  // E(Bijk)>0.9
  Eigen::VectorXd probi1;
  Eigen::VectorXd XBMi1;
  
  Eigen::VectorXd probi2;
  Eigen::VectorXd XBMi2;
  
  
  // E(Bijk, Bijk')>0.9
  Eigen::VectorXd probi21;
  Eigen::VectorXd XBMi21;
  
  
  
  double d1;
  double d2;
  double c1;
  double c2;
  
  // b=1
  double llk1;
  // b=0
  double llk2;
  
  // b=1
  double lll1;
  // b=0
  double lll2, lll3, lll4;
  
  Eigen::MatrixXd BB;
  Eigen::MatrixXd BB2;
  
  Eigen::VectorXd db1;
  Eigen::VectorXd cb1;
  
  Eigen::VectorXd ddb1;
  Eigen::VectorXd ccb1;
  
  double dd1;
  // b=1,b'=0
  double dd2;
  // b=0,b'=1
  double dd3;
  //b=0, b'=0
  double dd4;
  
  // b=1, b'=1
  double cc1;
  // b=1, b'=0
  double cc2;
  // b=0,b'=1
  double cc3;
  //b=0, b'=0
  double cc4;
  
  
  
  // only need calculate for node j
  Eigen::MatrixXd numerator=Eigen::MatrixXd::Zero(N,p-1);
  Eigen::MatrixXd prob=Eigen::MatrixXd::Zero(N,p-1);
  // Eigen::MatrixXd prob=Eigen::MatrixXd::Zero(N,po);
  
  // log of product of degrees: sijk=sij*sik / common neighbors in DTI
  Eigen::MatrixXd logS=Eigen::MatrixXd::Zero(N*p, p);
  logS=log(1+S.array());
  
  // log of distance 
  Eigen::MatrixXd logdist=Eigen::MatrixXd::Zero(p, p);
  logdist=log(distance.array());
  logdist.diagonal()=Eigen::VectorXd::Zero(p);
  
  j=nodej;
  
  s=0;
  for (k=0; k<p; ++k){
    if (k != j){
      // jk=indexjk(j,k)-1;
      
      for (i=0; i<N; ++i){
        
        numerator(i,s)=exp(gamma*logS(i*p+j,k)-eta*logdist(j,k));
        //  numerator(i*p+k,j)=numerator(i*p+j,k);
        
        prob(i,s)=numerator(i,s)/(1+numerator(i,s));
        //  prob0(i*p+k,j)=prob0(i*p+j,k);
        
      } // end of subject i
      
      indexk[k]=s;
      
      s=s+1;
      
    } // end of k!=j
    
  } // end of k
  
  
  // P(Bijk=b, Bijk'=b', Mi) b=1, b'=1
  Eigen::MatrixXd PBB1=Eigen::MatrixXd::Zero(N, pp);
  // P(Bijk=b, Bijk'=b', Mi) b=1, b'=0
  Eigen::MatrixXd PBB2=Eigen::MatrixXd::Zero(N, pp);
  // P(Bijk=b, Bijk'=b', Mi) b=0, b'=1
  Eigen::MatrixXd PBB3=Eigen::MatrixXd::Zero(N, pp);
  // P(Bijk=b, Bijk'=b', Mi) b=0, b'=0
  Eigen::MatrixXd PBB4=Eigen::MatrixXd::Zero(N, pp);
  
  // P(Bijk=1,  Mi) 
  Eigen::MatrixXd PB1=Eigen::MatrixXd::Zero(N, po);
  // P(Bijk=0, Mi) 
  Eigen::MatrixXd PB2=Eigen::MatrixXd::Zero(N, po);
  
  // E(BijkBijk'|Mi)
  Eigen::MatrixXd EBB=Eigen::MatrixXd::Zero(N, pp);
  
  // E(Bijk|Mi)
  Eigen::MatrixXd EB=Eigen::MatrixXd::Zero(N, po);
  
  
  j=nodej;
  Eigen::VectorXd prodprobj=Eigen::VectorXd::Ones(N);
  
  
  for (s=0; s<(p-1); ++s){
    prodprobj=(prodprobj.array())*(1-prob.col(s).array());
  } // end of s
  
  
  // calculate E(Bijk)
  
  for (i=0; i<N; ++i){
    
    // k=1;
    
    for (k=0; k<p; ++k){
      if (k!=j){
        
        indk=indexk[k];
        jk=indexjk(j,k)-1;
        
        // 0.1<=E(Bijk)<=0.9
        lenBi=EBji.row(i).sum()-EBji(i,indk);
        
        // E(Bijk)>0.9
        lenBi1=EBj1.row(i).sum()-EBj1(i,indk);
        
        if (lenBi==1){
          BB=B1;
        }
        if (lenBi==2){
          BB=B2;
        }
        if (lenBi==3){
          BB=B3;
        }
        if (lenBi==4){
          BB=B4;
        }
        if (lenBi==5){
          BB=B5;
        }
        if (lenBi==6){
          BB=B6;
        }
        if (lenBi==7){
          BB=B7;
        }
        if (lenBi==8){
          BB=B8;
        }
        
        
        // for each k, for summation
        probi=Eigen::VectorXd::Zero(lenBi);
        XBMi=Eigen::VectorXd::Zero(lenBi);
        
        // set bijk=1
        probi1=Eigen::VectorXd::Zero(lenBi1);
        XBMi1=Eigen::VectorXd::Zero(lenBi1);
        // 
        // // set bijk=0
        // probi0=Eigen::VectorXd::Zero(lenBi);
        // XBMi0=Eigen::VectorXd::Zero(lenBi);
        
        
        // if E(Bijk|...)>0.9, let bijk=1; if E(Bijk|...)<0.1, let bijk=0; if 0.1<=E(Bijk|...)<=0.9, calculate in the summation
        ind=0;
        ind1=0;
        
        for (l1=0; l1<p; ++l1){
          if ((l1!=j) & (l1!=k)){
            jl1=indexk[l1];
            
            // extract index needed in the summation (0.1<=E(Bijk)<=0.9)
            if (EBji(i,jl1)==1){
              probi[ind]=prob(i,jl1);
              XBMi[ind]=XBM(i,jl1);
              ++ind;
            } // end of EBji
            
            // E(Bijk)>0.9
            if (EBj1(i,jl1)==1){
              probi1[ind1]=prob(i,jl1);
              XBMi1[ind1]=XBM(i,jl1);
              ++ind1;
            }
            
            // // E(Bijk)<0.1
            // if (EBj0(i,jl1)==1){
            //   probi0[ind0]=prob(i,jl1);
            //   XBMi0[ind0]=XBM(i,jl1);
            //   ++ind0;
            // }
            
          } // end of l1!=j, l1!=k, Bj(i,l1)
        } // end of l1
        
        
        db1=Eigen::VectorXd::Zero(lenBi);
        cb1=Eigen::VectorXd::Zero(lenBi);
        
        llk1=0.0;
        llk2=0.0;
        
        d1=0.0;
        d2=0.0;
        c1=0.0;
        c2=0.0;
        
        
        if (lenBi==0){
          // if E(Bijk)>0.9, always bijk=1
          
          if (lenBi1>0){
            c1=XBM(i,indk)+XBMi1.sum();
            c2=XBMi1.sum();
            d1=log(prob(i,indk)/(1-prob(i,indk)))+(log(probi1.array()/(1-probi1.array()))).sum();
            d2=(log(probi1.array()/(1-probi1.array()))).sum();
            
          } else{
            c1=XBM(i,indk);
            c2=0;
            d1=log(prob(i,indk)/(1-prob(i,indk)));
            d2=0;
          }
          
          
          llk1=exp(-(M(i,j)-c1)*(M(i,j)-c1)/sigma(j)*0.5)*exp(d1);
          
          // sum over Bijl in {0,1}, b=0
          llk2=exp(-(M(i,j)-c2)*(M(i,j)-c2)/sigma(j)*0.5)*exp(d2);
          
          llk1=(prodprobj[i])*llk1/sqrt(2*M_PI*sigma(j));
          
          llk2=(prodprobj[i])*llk2/sqrt(2*M_PI*sigma(j));
          
        } else{
          
          rownum2=BB.rows();
          
          for (row2=0; row2<rownum2; ++row2){
            
            db1=BB.row(row2).array();
            db1=db1.array()*log(probi.array()/(1-probi.array()));
            
            // Bijl*betajl*Xi*Mil
            cb1=BB.row(row2).array();
            cb1=cb1.array()*XBMi.array();
            
            if (lenBi1>0){
              // b=1: log(pijk/(1-pijk))+sum_{l!=j,k}Bijl*log(pijl/(1-pijl))
              d1=db1.sum()+log(prob(i,indk)/(1-prob(i,indk)))+(log(probi1.array()/(1-probi1.array()))).sum();
              // b=1: betajk*Xi*Mik+sum_{l!=j,k}Bijl*betajl*Xi*Mil
              c1=cb1.sum()+XBM(i,indk)+XBMi1.sum();
              
              // b=0
              d2=db1.sum()+(log(probi1.array()/(1-probi1.array()))).sum();
              // b=0; sum_{l!=j,k}Bijl*betajl*Xi*Mil
              c2=cb1.sum()+XBMi1.sum();
              
            } else{
              // b=1: log(pijk/(1-pijk))+sum_{l!=j,k}Bijl*log(pijl/(1-pijl))
              d1=db1.sum()+log(prob(i,indk)/(1-prob(i,indk)));
              // b=1: betajk*Xi*Mik+sum_{l!=j,k}Bijl*betajl*Xi*Mil
              c1=cb1.sum()+XBM(i,indk);
              
              // b=0
              d2=db1.sum();
              // b=0; sum_{l!=j,k}Bijl*betajl*Xi*Mil
              c2=cb1.sum();
            }
            
            
            
            // sum over Bijl in {0,1}, b=1
            llk1=llk1+exp(-(M(i,j)-c1)*(M(i,j)-c1)/sigma(j)*0.5)*exp(d1);
            
            // sum over Bijl in {0,1}, b=0
            llk2=llk2+exp(-(M(i,j)-c2)*(M(i,j)-c2)/sigma(j)*0.5)*exp(d2);
          } // end of row 2
          
          llk1=(prodprobj[i])*llk1/sqrt(2*M_PI*sigma(j));
          
          llk2=(prodprobj[i])*llk2/sqrt(2*M_PI*sigma(j));
          
        }  // end of else
        
        
        // P(Bijk=1, Mi) 
        PB1(i,jk)=llk1;
        
        // P(Bijk=0, Mi) 
        PB2(i,jk)=llk2;
        
        // expectation E(Bijk)
        EB(i,jk)=PB1(i,jk)/(PB1(i,jk)+PB2(i,jk));
        
        
        
        
        //////////////////////////////////////////
        //////    start computing EBB     ////////
        //////////////////////////////////////////
        
        if (k<(p-1)){
          for (k2=k+1; k2<p; ++k2){
            if (k2!=j){
              
              indk2=indexk[k2];
              
              lenBBi=EBji.row(i).sum()-EBji(i,indk)-EBji(i,indk2);
              
              // E(Bijk)>0.9
              lenBBi1=EBj1.row(i).sum()-EBj1(i,indk)-EBj1(i,indk2);
              
              
              if (lenBBi==1){
                BB2=B1;
              }
              
              if (lenBBi==2){
                BB2=B2;
              }
              if (lenBBi==3){
                BB2=B3;
              }
              if (lenBBi==4){
                BB2=B4;
              }
              if (lenBBi==5){
                BB2=B5;
              }
              if (lenBBi==6){
                BB2=B6;
              }
              if (lenBBi==7){
                BB2=B7;
              }
              
              jk2=indexjk(j,k2)-1;
              jkk=indexjkk(j*p+k,k2)-1;
              
              probi2=Eigen::VectorXd::Zero(lenBBi);
              XBMi2=Eigen::VectorXd::Zero(lenBBi);
              
              // E(Bijk)>0.9
              probi21=Eigen::VectorXd::Zero(lenBBi1);
              XBMi21=Eigen::VectorXd::Zero(lenBBi1);
              
              indb=0;
              indb1=0;
              for (l=0; l<p; ++l){
                if ((l!=j) & (l!=k) & (l!=k2)){
                  jl=indexk[l];
                  
                  // extract index needed in the summation
                  if (EBji(i,jl)==1){
                    
                    probi2[indb]=prob(i,jl);
                    XBMi2[indb]=XBM(i,jl);
                    ++indb;
                  }  // end of EBji
                  
                  // E(Bijk)>0.9
                  if (EBj1(i,jl)==1){
                    
                    probi21[indb1]=prob(i,jl);
                    XBMi21[indb1]=XBM(i,jl);
                    ++indb1;
                  }  // end of EBji
                  
                  
                } // end of l!=j, l!=k, l!=k2
              } // end of l
              
              
              ddb1=Eigen::VectorXd::Zero(lenBBi);
              ccb1=Eigen::VectorXd::Zero(lenBBi);
              
              dd1=0.0;
              // b=1,b'=0
              dd2=0.0;
              // b=0,b'=1
              dd3=0.0;
              //b=0, b'=0
              dd4=0.0;
              
              // b=1, b'=1
              cc1=0.0;
              // b=1, b'=0
              cc2=0.0;
              // b=0,b'=1
              cc3=0.0;
              //b=0, b'=0
              cc4=0.0;
              
              lll1=0.0;
              lll2=0.0;
              lll3=0.0;
              lll4=0.0;
              
              if (lenBBi==0){
                
                if (lenBBi1>0){
                  // // b=1, b'=1
                  dd1=log(prob(i,indk)/(1-prob(i,indk)))+log(prob(i,indk2)/(1-prob(i,indk2)))+(log(probi21.array()/(1-probi21.array()))).sum();
                  
                  // b=1, b'=0
                  dd2=log(prob(i,indk)/(1-prob(i,indk)))+(log(probi21.array()/(1-probi21.array()))).sum();
                  // 
                  // b=0, b'=1
                  dd3=log(prob(i,indk2)/(1-prob(i,indk2)))+(log(probi21.array()/(1-probi21.array()))).sum();
                  
                  // b=0, b'=0
                  dd4=(log(probi21.array()/(1-probi21.array()))).sum();
                  
                  // b=1, b'=1
                  cc1=XBM(i,indk)+XBM(i,indk2)+XBMi21.sum();
                  // b=1, b'=0
                  cc2=XBM(i,indk)+XBMi21.sum();
                  // b=0, b'=1
                  cc3=XBM(i,indk2)+XBMi21.sum();
                  // b=0, b'=0
                  cc4=XBMi21.sum();
                  
                } else {
                  
                  // b=1, b'=1
                  dd1=log(prob(i,indk)/(1-prob(i,indk)))+log(prob(i,indk2)/(1-prob(i,indk2)));
                  
                  // b=1, b'=0
                  dd2=log(prob(i,indk)/(1-prob(i,indk)));
                  
                  // b=0, b'=1
                  dd3=log(prob(i,indk2)/(1-prob(i,indk2)));
                  
                  dd4=0;
                  
                  // b=1, b'=1
                  cc1=XBM(i,indk)+XBM(i,indk2);
                  // b=1, b'=0
                  cc2=XBM(i,indk);
                  // b=0, b'=1
                  cc3=XBM(i,indk2);
                  
                  cc4=0;
                }  // end of else
                
                lll1=exp(-(M(i,j)-cc1)*(M(i,j)-cc1)/sigma(j)*0.5)*exp(dd1);
                
                // sum over Bijl in {0,1}, b=1, b'=0
                lll2=exp(-(M(i,j)-cc2)*(M(i,j)-cc2)/sigma(j)*0.5)*exp(dd2);
                
                // sum over Bijl in {0,1}, b=0, b'=1
                lll3=exp(-(M(i,j)-cc3)*(M(i,j)-cc3)/sigma(j)*0.5)*exp(dd3);
                
                // sum over Bijl in {0,1}, b=0, b'=0
                lll4=exp(-(M(i,j)-cc4)*(M(i,j)-cc4)/sigma(j)*0.5)*exp(dd4);
                
                lll1=(prodprobj[i])*lll1/sqrt(2*M_PI*sigma(j));
                lll2=(prodprobj[i])*lll2/sqrt(2*M_PI*sigma(j));
                lll3=(prodprobj[i])*lll3/sqrt(2*M_PI*sigma(j));
                lll4=(prodprobj[i])*lll4/sqrt(2*M_PI*sigma(j));
                
                
              } else {
                
                rownum3=BB2.rows();
                
                for (row3=0; row3<rownum3; ++row3){
                  
                  ddb1=BB2.row(row3).array();
                  ddb1=ddb1.array()*log(probi2.array()/(1-probi2.array()));
                  
                  ccb1=BB2.row(row3).array();
                  ccb1=ccb1.array()*(XBMi2.array());
                  
                  
                  if (lenBBi1>0){
                    
                    // b=1, b'=1
                    dd1=ddb1.sum()+log(prob(i,indk)/(1-prob(i,indk)))+log(prob(i,indk2)/(1-prob(i,indk2)))+(log(probi21.array()/(1-probi21.array()))).sum();
                    
                    // b=1, b'=0
                    dd2=ddb1.sum()+log(prob(i,indk)/(1-prob(i,indk)))+(log(probi21.array()/(1-probi21.array()))).sum();
                    
                    // b=0, b'=1
                    dd3=ddb1.sum()+log(prob(i,indk2)/(1-prob(i,indk2)))+(log(probi21.array()/(1-probi21.array()))).sum();
                    
                    // b=0, b'=0
                    dd4=ddb1.sum()+(log(probi21.array()/(1-probi21.array()))).sum();
                    
                    
                    // b=1, b'=1
                    cc1=ccb1.sum()+XBM(i,indk)+XBM(i,indk2)+XBMi21.sum();
                    // b=1, b'=0
                    cc2=ccb1.sum()+XBM(i,indk)+XBMi21.sum();
                    // b=0, b'=1
                    cc3=ccb1.sum()+XBM(i,indk2)+XBMi21.sum();
                    // b=0, b'=0
                    cc4=ccb1.sum()+XBMi21.sum();
                    
                    
                  } else {
                    
                    // b=1, b'=1
                    dd1=ddb1.sum()+log(prob(i,indk)/(1-prob(i,indk)))+log(prob(i,indk2)/(1-prob(i,indk2)));
                    
                    // b=1, b'=0
                    dd2=ddb1.sum()+log(prob(i,indk)/(1-prob(i,indk)));
                    
                    // b=0, b'=1
                    dd3=ddb1.sum()+log(prob(i,indk2)/(1-prob(i,indk2)));
                    
                    // b=0, b'=0
                    dd4=ddb1.sum();
                    
                    
                    // b=1, b'=1
                    cc1=ccb1.sum()+XBM(i,indk)+XBM(i,indk2);
                    // b=1, b'=0
                    cc2=ccb1.sum()+XBM(i,indk);
                    // b=0, b'=1
                    cc3=ccb1.sum()+XBM(i,indk2);
                    // b=0, b'=0
                    cc4=ccb1.sum();
                    
                    
                  } // end of else
                  
                  
                  // sum over Bijl in {0,1}, b=1, b'=1
                  lll1=lll1+exp(-(M(i,j)-cc1)*(M(i,j)-cc1)/sigma(j)*0.5)*exp(dd1);
                  
                  // sum over Bijl in {0,1}, b=1, b'=0
                  lll2=lll2+exp(-(M(i,j)-cc2)*(M(i,j)-cc2)/sigma(j)*0.5)*exp(dd2);
                  
                  // sum over Bijl in {0,1}, b=0, b'=1
                  lll3=lll3+exp(-(M(i,j)-cc3)*(M(i,j)-cc3)/sigma(j)*0.5)*exp(dd3);
                  
                  // sum over Bijl in {0,1}, b=0, b'=0
                  lll4=lll4+exp(-(M(i,j)-cc4)*(M(i,j)-cc4)/sigma(j)*0.5)*exp(dd4);
                  
                  
                } // end of row 3
                
                lll1=(prodprobj[i])*lll1/sqrt(2*M_PI*sigma(j));
                lll2=(prodprobj[i])*lll2/sqrt(2*M_PI*sigma(j));
                lll3=(prodprobj[i])*lll3/sqrt(2*M_PI*sigma(j));
                lll4=(prodprobj[i])*lll4/sqrt(2*M_PI*sigma(j));
                
                
              }
              
              
              // P(Bijk=b, Bijk'=b', Mi) b=1, b'=1
              PBB1(i,jkk)=lll1;
              
              // P(Bijk=b, Bijk'=b', Mi) b=1, b'=0
              PBB2(i,jkk)=lll2;
              // P(Bijk=b, Bijk'=b', Mi) b=0, b'=1
              PBB3(i,jkk)=lll3;
              
              // P(Bijk=b, Bijk'=b', Mi) b=0, b'=0 
              PBB4(i,jkk)=lll4;
              
              // expectation E(BijkBijk')
              EBB(i,jkk)=PBB1(i,jkk)/(PBB1(i,jkk)+PBB2(i,jkk)+PBB3(i,jkk)+PBB4(i,jkk));
              
            } // end of  k2 != j
            
          } // end of k2
        } // end of if k<p-1
        
      } // end of k != j
      
    } // end of k
    
  } // end of subject i
  
  
  
  //  } // end of j
  
  return(List::create(Named("prob")=prob, Named("prodprobj")=prodprobj, Named("PBB1")=PBB1,Named("PBB2")=PBB2,Named("PBB3")=PBB3,Named("PBB4")=PBB4,
                            Named("PB1")=PB1, Named("PB2")=PB2, Named("EBB")=EBB, Named("EB")=EB, Named("indexk")=indexk, Named("probi")=probi, Named("indk")=indk
  ));
  
}



/********* E-step, expectation ********/
// [[Rcpp::export]]
List ExpectApproxjC(Eigen::MatrixXd Beta, Eigen::MatrixXd X, Eigen::MatrixXd M,
                    Eigen::MatrixXd indexjk, Eigen::MatrixXd indexjkk, 
                    Eigen::MatrixXd XBM, Eigen::MatrixXd UMi, Eigen::MatrixXd UMMi,
                    Eigen::VectorXd sigma,  Eigen::MatrixXd prob,  Eigen::VectorXd lambda1, 
                    Eigen::VectorXd coeff1, int nodej){   
  
  //  int j, int k, int k2
  int N=X.rows(), p=M.cols(), pp=p*(p-1)*(p-2)/2, po=p*(p-1);
  int l,s,l1;
  int i,j,k,k2;
  int jl, jl1, jk, kj, jk2, jkk, kjk,countb, count, countk;  
  // int m=lambda.size();
  int p2=p-3;
  
  // Eigen::MatrixXd numerator=Eigen::MatrixXd::Zero(N*p,p);
  // Eigen::MatrixXd prob0=Eigen::MatrixXd::Zero(N*p,p);
  // Eigen::MatrixXd prob=Eigen::MatrixXd::Zero(N,po);
  // 
  // // log of product of degrees: sijk=sij*sik / common neighbors in DTI
  // Eigen::MatrixXd logS=Eigen::MatrixXd::Zero(N*p, p);
  // logS=log(1+S.array());
  // 
  // // log of distance 
  // Eigen::MatrixXd logdist=Eigen::MatrixXd::Zero(p, p);
  // logdist=log(distance.array());
  // logdist.diagonal()=Eigen::VectorXd::Zero(p);
  // 
  // for (i=0; i<N; ++i){
  //   
  //   for (j=0; j<(p-1); ++j){
  //     for (k=j+1; k<p; ++k){
  // 
  //       numerator(i*p+j,k)=exp(gamma*logS(i*p+j,k)-eta*logdist(j,k));
  //       numerator(i*p+k,j)=numerator(i*p+j,k);
  //       
  //       prob0(i*p+j,k)=numerator(i*p+j,k)/(1+numerator(i*p+j,k));
  //       prob0(i*p+k,j)=prob0(i*p+j,k);
  //       
  //     } // end of k
  //   }  // end of j
  //   
  // } // end of i
  // 
  // 
  // 
  // for (i=0; i<N; ++i){
  //   s=0;
  //   for (j=0; j<p; ++j){
  //     for (k=0; k<p; ++k){
  //       if(k!=j){
  //         prob(i,s)=prob0(i*p+j,k);
  //         ++s;
  //       }
  //     }
  //   }
  // }
  
  
  // P(Bijk=b, Bijk'=b', Mi) b=1, b'=1
  Eigen::MatrixXd PBB1=Eigen::MatrixXd::Zero(N, pp);
  // P(Bijk=b, Bijk'=b', Mi) b=1, b'=0
  Eigen::MatrixXd PBB2=Eigen::MatrixXd::Zero(N, pp);
  // P(Bijk=b, Bijk'=b', Mi) b=0, b'=1
  Eigen::MatrixXd PBB3=Eigen::MatrixXd::Zero(N, pp);
  // P(Bijk=b, Bijk'=b', Mi) b=0, b'=0
  Eigen::MatrixXd PBB4=Eigen::MatrixXd::Zero(N, pp);
  
  // P(Bijk=1,  Mi) 
  Eigen::MatrixXd PB1=Eigen::MatrixXd::Zero(N, po);
  // P(Bijk=0, Mi) 
  Eigen::MatrixXd PB2=Eigen::MatrixXd::Zero(N, po);
  
  // E(BijkBijk'|Mi)
  Eigen::MatrixXd EBB=Eigen::MatrixXd::Zero(N, pp);
  
  // E(Bijk|Mi)
  Eigen::MatrixXd EB=Eigen::MatrixXd::Zero(N, po);
  
  
  j=nodej;
  
  Eigen::VectorXd prodprobj=Eigen::VectorXd::Ones(N);
  
  for (l=0; l<p; ++l){
    if (l!=j){
      jl=indexjk(j,l)-1;
      prodprobj=(prodprobj.array())*(1-prob.col(jl).array());
    } // end of l!=j
  } // end of l
  
  
  //  k=1;
  for (k=0; k<p; ++k){
    if (k!=j){
      jk=indexjk(j,k)-1;
      // calculate E(Bijk)
      
      // b=1, 
      Eigen::MatrixXd db1=Eigen::MatrixXd::Zero(N, p-2);
      // b=0
      Eigen::MatrixXd db2=Eigen::MatrixXd::Zero(N, p-2);
      
      Eigen::MatrixXd cb=Eigen::MatrixXd::Zero(N, p-2);
      
      countb=0;
      
      for (l1=0; l1<p; ++l1){
        if ((l1!=j) & (l1!=k)){
          jl1=indexjk(j,l1)-1;
          
          // compute P(Bijk=b,Mi)
          db1.col(countb)=log(prob.col(jl1).array()/(1-prob.col(jl1).array()));
          db2.col(countb)=db1.col(countb);
          
          // if (indexAct(j,l1)>0){
          cb.col(countb)=(XBM.col(jl1).array())/sqrt(2*sigma(j));
          
          // b=1
          db1.col(countb)=db1.col(countb).array()+(UMi.col(jk).array())*(XBM.col(jl1).array())/sigma(j);
          // b=0
          db2.col(countb)=db2.col(countb).array()+(M.col(j).array())*(XBM.col(jl1).array())/sigma(j);
          
          ++countb;
        } // end of (l1!=j) & (l1!=k)
      } // end of l1
      
      Eigen::VectorXd constantb1=Eigen::VectorXd::Zero(N);
      Eigen::VectorXd constantb2=Eigen::VectorXd::Zero(N);
      
      constantb1=linearApproxC(db1, cb,  lambda1, coeff1);
      
      constantb2=linearApproxC(db2, cb,  lambda1, coeff1);
      //      
      //     
      //     // P(Bijk=1, Mi) 
      PB1.col(jk)=(prodprobj.array())*exp(-(UMi.col(jk).array())*(UMi.col(jk).array())/sigma(j)*0.5)*(prob.col(jk).array())/(1-prob.col(jk).array())/sqrt(2*M_PI*sigma(j));
      //     
      //     
      //     //   PB1.col(jk)=exp(-(UMi.col(jk).array())*(UMi.col(jk).array())/sigma(j)*0.5);
      //     
      //     
      PB1.col(jk)=(PB1.col(jk).array())*(constantb1.array());
      //     
      //     
      //     // P(Bijk=0, Mi) 
      PB2.col(jk)=(prodprobj.array())*exp(-(M.col(j).array())*(M.col(j).array())/sigma(j)*0.5)/sqrt(2*M_PI*sigma(j));
      //     
      //     
      //     //    PB2.col(jk)=(prodprobj.array())/sqrt(2*M_PI*sigma(j));
      //     
      PB2.col(jk)=(PB2.col(jk).array())*(constantb2.array());
      //     
      for (i=0; i<N; ++i){
        if (PB1(i,jk)==0 & PB2(i,jk)==0){
          EB(i,jk)=1;
        }else{
          // expectation E(Bijk)
          EB(i,jk)=PB1(i,jk)/(PB1(i,jk)+PB2(i,jk));
        }
      } // end of i
      
      
      //// start computing E(BijkBijk')
      if (k<(p-1)){
        for (k2=k+1; k2<p; ++k2){
          if (k2!=j){
            //           
            //           // b=1, b'=1
            Eigen::MatrixXd d1=Eigen::MatrixXd::Zero(N, p2);
            // b=1,b'=0
            Eigen::MatrixXd d2=Eigen::MatrixXd::Zero(N, p2);
            // b=0,b'=1
            Eigen::MatrixXd d3=Eigen::MatrixXd::Zero(N, p2);
            //b=0, b'=0
            Eigen::MatrixXd d4=Eigen::MatrixXd::Zero(N, p2);
            
            Eigen::MatrixXd c=Eigen::MatrixXd::Zero(N, p2);
            
            jk2=indexjk(j,k2)-1;
            
            jkk=indexjkk(j*p+k,k2)-1;
            
            count=0;
            //           
            for (l=0; l<p; ++l){
              if ((l!=j) & (l!=k) & (l!=k2)){
                jl=indexjk(j,l)-1;
                //               
                d1.col(count)=log(prob.col(jl).array()/(1-prob.col(jl).array()));
                d2.col(count)=d1.col(count);
                d3.col(count)=d1.col(count);
                d4.col(count)=d1.col(count);
                //               
                //               // if (indexAct(j,l)>0){
                c.col(count)=(XBM.col(jl).array())/sqrt(2*sigma(j));
                
                //               // b=1, b'=1
                d1.col(count)=d1.col(count).array()+(UMMi.col(jkk).array())*(XBM.col(jl).array())/sigma(j);
                // b=1, b'=0
                d2.col(count)=d2.col(count).array()+(UMi.col(jk).array())*(XBM.col(jl).array())/sigma(j);
                // b=0, b'=1
                d3.col(count)=d3.col(count).array()+(UMi.col(jk2).array())*(XBM.col(jl).array())/sigma(j);
                // b=0, b'=0
                d4.col(count)=d4.col(count).array()+(M.col(j).array())*(XBM.col(jl).array())/sigma(j);
                
                ++count;
              } // end of (l!=j) & (l!=k) & (l!=k2)
            } // end of l
            
            Eigen::VectorXd constant1=Eigen::VectorXd::Zero(N);
            Eigen::VectorXd constant2=Eigen::VectorXd::Zero(N);
            Eigen::VectorXd constant3=Eigen::VectorXd::Zero(N);
            Eigen::VectorXd constant4=Eigen::VectorXd::Zero(N);
            
            constant1=linearApproxC(d1,  c, lambda1, coeff1);
            //          
            constant2=linearApproxC(d2,  c, lambda1, coeff1);
            
            constant3=linearApproxC(d3,  c, lambda1, coeff1);
            //       
            constant4=linearApproxC(d4,  c, lambda1, coeff1);
            // 
            //           
            // P(Bijk=b, Bijk'=b', Mi) b=1, b'=1
            PBB1.col(jkk)=(prodprobj.array())*exp(-(UMMi.col(jkk).array())*(UMMi.col(jkk).array())/sigma(j)*0.5)*(prob.col(jk).array()/(1-prob.col(jk).array()))*
              (prob.col(jk2).array()/(1-prob.col(jk2).array()))/sqrt(2*M_PI*sigma(j));
            //           // 
            //           
            //           // PBB1.col(jkk)=(prodprobj.array())*(prob.col(jk).array()/(1-prob.col(jk).array()))*
            //           //   (prob.col(jk2).array()/(1-prob.col(jk2).array()))/sqrt(2*M_PI*sigma(j));
            //           
            //           
            PBB1.col(jkk)=(PBB1.col(jkk).array())*(constant1.array());
            //           
            //           // P(Bijk=b, Bijk'=b', Mi) b=1, b'=0
            //           //   PBB2.col(jkk)=(prodprobj.array())*(prob.col(jk).array())/(1-prob.col(jk).array())/sqrt(2*M_PI*sigma(j));
            //           
            //           
            PBB2.col(jkk)=(prodprobj.array())*exp(-(UMi.col(jk).array())*(UMi.col(jk).array())/sigma(j)*0.5)*(prob.col(jk).array())/(1-prob.col(jk).array())/sqrt(2*M_PI*sigma(j));
            //           
            PBB2.col(jkk)=(PBB2.col(jkk).array())*(constant2.array());
            //           
            //           // P(Bijk=b, Bijk'=b', Mi) b=0, b'=1
            PBB3.col(jkk)=(prodprobj.array())*exp(-(UMi.col(jk2).array())*(UMi.col(jk2).array())/sigma(j)*0.5)*(prob.col(jk2).array())/(1-prob.col(jk2).array())/sqrt(2*M_PI*sigma(j));
            //           
            //           
            //           //    PBB3.col(jkk)=(prodprobj.array())*(prob.col(jk2).array())/(1-prob.col(jk2).array())/sqrt(2*M_PI*sigma(j));
            PBB3.col(jkk)=(PBB3.col(jkk).array())*(constant3.array());
            //           
            //           // P(Bijk=b, Bijk'=b', Mi) b=0, b'=0 
            PBB4.col(jkk)=(prodprobj.array())*exp(-(M.col(j).array())*(M.col(j).array())/sigma(j)*0.5)/sqrt(2*M_PI*sigma(j));
            //           
            //           
            //           //  PBB4.col(jkk)=(prodprobj.array())/sqrt(2*M_PI*sigma(j));
            PBB4.col(jkk)=(PBB4.col(jkk).array())*(constant4.array());
            //           
            //           
            //           // expectation E(BijkBijk')
            for (i=0; i<N; ++i){
              if (PBB1(i,jkk)==0 & PBB2(i,jkk)==0 & PBB3(i,jkk)==0 & PBB4(i,jkk)==0){
                EBB(i,jkk)=1;
              }else{
                //               // expectation E(Bijk)
                EBB(i,jkk)=PBB1(i,jkk)/(PBB1(i,jkk)+PBB2(i,jkk)+PBB3(i,jkk)+PBB4(i,jkk));
              }
            }
            //           
            //           // EBB.col(jkk)=PBB1.col(jkk).array()/(PBB1.col(jkk).array()+PBB2.col(jkk).array()+PBB3.col(jkk).array()+PBB4.col(jkk).array());
            //           
          } // end of  k2 != j
        } // end of k2
      } // end of if k<p-1
      
    } // end of k != j
  } // end of k
  
  
  return(List::create(Named("prodprobj")=prodprobj,
                      Named("PBB1")=PBB1,Named("PBB2")=PBB2,Named("PBB3")=PBB3,Named("PBB4")=PBB4,
                            Named("PB1")=PB1, Named("PB2")=PB2, Named("EBB")=EBB, Named("EB")=EB
                        // , Named("cb")=cb,
                        //     Named("db1")=db1, Named("db2")=db2, Named("constantb1")=constantb1, Named("constantb2")=constantb2
  ));
  //        Named("eMbeta1")=Mbeta1, Named("eMbeta2")=Mbeta2
  
  
}





/********* Compute obs likelihood function *******/
// [[Rcpp::export]]
List LikeliObsjC(Eigen::MatrixXd Beta, Eigen::MatrixXd X, Eigen::MatrixXd M,
                 Eigen::MatrixXd indexAct, Eigen::MatrixXd indexjk, Eigen::MatrixXd indexjkk, 
                 Eigen::MatrixXd XBM, Eigen::MatrixXd B,
                 Eigen::VectorXd sigma, Eigen::MatrixXd S, Eigen::MatrixXd distance,  double gamma, double eta,  Eigen::VectorXd lambda,
                 Eigen::VectorXd coeff, int nodej){   
  
  //  int j, int k, int k2
  int N=X.rows(), p=M.cols(), pp=p*(p-1)*(p-2)/2, po=p*(p-1);
  int l,s,l1;
  int i,j,k,k2;
  int jl, jl1, jk, kj, jk2, jkk, kjk,  count, countk;  
  int row;
  int rownum=B.rows();
  int m=lambda.size();
  int p2=p-3;
  
  
  Eigen::MatrixXd numerator=Eigen::MatrixXd::Zero(N*p,p);
  Eigen::MatrixXd prob0=Eigen::MatrixXd::Zero(N*p,p);
  
  Eigen::MatrixXd prob=Eigen::MatrixXd::Zero(N,po);
  
  Eigen::MatrixXd db, cb;
  
  Eigen::VectorXd d1, c1;
  
  Eigen::VectorXd llk;
  
  // log of product of degrees: sijk=sij*sik / common neighbors in DTI
  Eigen::MatrixXd logS=Eigen::MatrixXd::Zero(N*p, p);
  logS=log(1+S.array());
  
  // log of distance 
  Eigen::MatrixXd logdist=Eigen::MatrixXd::Zero(p, p);
  logdist=log(distance.array());
  logdist.diagonal()=Eigen::VectorXd::Zero(p);
  
  Eigen::VectorXd constant;
  
  Eigen::VectorXd Likelihoodi;
  double Likeli;
  double sumLikeli=0.0;  
  
  Eigen::VectorXd prodprobj;
  
  
  for (i=0; i<N; ++i){
    
    for (j=0; j<(p-1); ++j){
      for (k=j+1; k<p; ++k){
        // for (k=0; k<p; ++k){
        //   if (k != j){
        // jk=indexjk(j,k)-1;
        numerator(i*p+j,k)=exp(gamma*logS(i*p+j,k)-eta*logdist(j,k));
        numerator(i*p+k,j)=numerator(i*p+j,k);
        
        prob0(i*p+j,k)=numerator(i*p+j,k)/(1+numerator(i*p+j,k));
        prob0(i*p+k,j)=prob0(i*p+j,k);
        
        //    } // end of k!=j
      } // end of k
    }  // end of j
    
    
  } // end of i
  
  
  for (i=0; i<N; ++i){
    s=0;
    for (j=0; j<p; ++j){
      for (k=0; k<p; ++k){
        if(k!=j){
          prob(i,s)=prob0(i*p+j,k);
          ++s;
        }
      }
    }
  }  // end of i
  
  
  
  // compute likelihood function for observed data: P(Mij, Bijk, k!=j | Mik, k!=j) 
  Likelihoodi=Eigen::VectorXd::Zero(N);
  //  Likeli=Eigen::VectorXd::Zero(p);
  
  Likeli=0.0;
  
  //  for (j=0; j<p; ++j){
  j=nodej;
  // prod_{k!=j}(1-pijk)
  prodprobj=Eigen::VectorXd::Ones(N);
  
  for (l=0; l<p; ++l){
    if (l!=j){
      jl=indexjk(j,l)-1;
      prodprobj=(prodprobj.array())*(1-prob.col(jl).array());
    } // end of l!=j
  } // end of l
  
  
  db=Eigen::MatrixXd::Zero(N, p-1);
  cb=Eigen::MatrixXd::Zero(N, p-1);
  
  d1=Eigen::VectorXd::Zero(N);
  c1=Eigen::VectorXd::Zero(N);
  
  llk=Eigen::VectorXd::Zero(N);
  
  for (row=0; row<rownum; ++row){
    
    count=0;
    for (k=0; k<p; ++k){
      if (k!=j){
        
        jk=indexjk(j,k)-1;
        
        // Bijk*log(pijk/(1-pijk))
        db.col(count)=B(row, count)*log(prob.col(jk).array()/(1-prob.col(jk).array()));
        
        // Bijk*betajk*Xi*Mik
        cb.col(count)=B(row,count)*XBM.col(jk).array();
        
        //  db.col(count)=db.col(count).array()+(M.col(j).array())*(XBM.col(jk).array())/sigma(j);
        
        // /sqrt(2*sigma(j))
        
        ++count;
        
      } // end of k != j
    } // end of k
    
    // sum_{k!=j} Bijk*log(pijk/(1-pijk))
    d1=db.rowwise().sum();
    // sum_{k!=j} Bijk*betajk*Xi*Mik
    c1=cb.rowwise().sum();
    
    // sum over Bijk in {0,1}
    llk=llk.array()+exp(-(M.col(j).array()-c1.array())*(M.col(j).array()-c1.array())/sigma(j)*0.5)*exp(d1.array());
    
  } // end of row
  
  Likelihoodi=(prodprobj.array())*llk.array()/sqrt(2*M_PI*sigma(j));
  
  // } // end of j
  
  // take log, loglikelihood
  Likelihoodi=(Likelihoodi.array()).log();
  
  // sum_{i:n}
  Likeli=Likelihoodi.sum();
  
  // add likelihood of p nodes together
  // sumLikeli=Likeli.sum();
  
  
  return(List::create(Named("constant")=constant, Named("Likelihoodi")=Likelihoodi, Named("Likeli")=Likeli
  ));
  
}



/********* Compute obs likelihood function *******/
// [[Rcpp::export]]
List LikeliObsDTIjC(Eigen::MatrixXd M, Eigen::MatrixXd EBji, Eigen::MatrixXd EBj1,
                    Eigen::MatrixXd indexjk, Eigen::MatrixXd indexjkk, 
                    Eigen::MatrixXd XBM, Eigen::MatrixXd B1, Eigen::MatrixXd B2, Eigen::MatrixXd B3,
                    Eigen::MatrixXd B4, Eigen::MatrixXd B5, Eigen::MatrixXd B6, Eigen::MatrixXd B7, Eigen::MatrixXd B8, Eigen::MatrixXd B9,
                    Eigen::VectorXd sigma, Eigen::MatrixXd S, Eigen::MatrixXd distance,  double gamma, double eta, int nodej){   
  
  // Bj: previous iteration posterior E(B_ijk) for node j only; XBM: XBM for node j only
  //  int j, int k, int k2
  int N=M.rows(), p=M.cols(), pp=p*(p-1)*(p-2)/2, po=p*(p-1);
  int l,s,l1;
  int i,j,k,k2;
  int jl, jl1, jk, kj, jk2, jkk, kjk,  count, countk;  
  int row;
  int rownum;
  // int m=lambda.size();
  int p2=p-3;
  int lenB, lenB1;
  int ind, ind1;
  int indk;
  
  Eigen::MatrixXd BB;
  
  Eigen::VectorXd indexk=Eigen::VectorXd::Zero(p);
  Eigen::VectorXd probi;
  Eigen::VectorXd XBMi;
  
  Eigen::VectorXd probi1;
  Eigen::VectorXd XBMi1;
  
  // only need calculate for node j
  Eigen::MatrixXd numerator=Eigen::MatrixXd::Zero(N,p-1);
  Eigen::MatrixXd prob=Eigen::MatrixXd::Zero(N,p-1);
  
  Eigen::VectorXd db, cb;
  
  double d1, c1;
  
  double llk;
  
  // log of product of degrees: sijk=sij*sik / common neighbors in DTI
  Eigen::MatrixXd logS=Eigen::MatrixXd::Zero(N*p, p);
  logS=log(1+S.array());
  
  // log of distance 
  Eigen::MatrixXd logdist=Eigen::MatrixXd::Zero(p, p);
  logdist=log(distance.array());
  logdist.diagonal()=Eigen::VectorXd::Zero(p);
  
  Eigen::VectorXd constant;
  
  Eigen::VectorXd Likelihoodi;
  double Likeli;
  double sumLikeli=0.0;  
  
  Eigen::VectorXd prodprobj;
  
  
  j=nodej;
  
  s=0;
  for (k=0; k<p; ++k){
    if (k != j){
      // jk=indexjk(j,k)-1;
      
      for (i=0; i<N; ++i){
        
        numerator(i,s)=exp(gamma*logS(i*p+j,k)-eta*logdist(j,k));
        //  numerator(i*p+k,j)=numerator(i*p+j,k);
        
        prob(i,s)=numerator(i,s)/(1+numerator(i,s));
        //  prob0(i*p+k,j)=prob0(i*p+j,k);
        
      } // end of subject i
      
      indexk[k]=s;
      
      s=s+1;
      
    } // end of k!=j
    
  } // end of k
  
  
  // compute likelihood function for observed data: P(Mij, Bijk, k!=j | Mik, k!=j) 
  Likelihoodi=Eigen::VectorXd::Zero(N);
  //  Likeli=Eigen::VectorXd::Zero(p);
  
  Likeli=0.0;
  
  
  prodprobj=Eigen::VectorXd::Ones(N);
  
  for (s=0; s<(p-1); ++s){
    prodprobj=(prodprobj.array())*(1-prob.col(s).array());
  } // end of s
  
  
  for (i=0; i<N; ++i){
    
    // 0.1<=E(Bijk)<=0.9
    lenB=EBji.row(i).sum();
    
    // E(Bijk)>0.9
    lenB1=EBj1.row(i).sum();
    
    
    if (lenB==1){
      BB=B1;
    }
    if (lenB==2){
      BB=B2;
    }
    if (lenB==3){
      BB=B3;
    }
    if (lenB==4){
      BB=B4;
    }
    if (lenB==5){
      BB=B5;
    }
    if (lenB==6){
      BB=B6;
    }
    if (lenB==7){
      BB=B7;
    }
    if (lenB==8){
      BB=B8;
    }
    if (lenB==9){
      BB=B9;
    }
    
    // for each k, for summation
    probi=Eigen::VectorXd::Zero(lenB);
    XBMi=Eigen::VectorXd::Zero(lenB);
    
    // set bijk=1
    probi1=Eigen::VectorXd::Zero(lenB1);
    XBMi1=Eigen::VectorXd::Zero(lenB1);
    
    // if E(Bijk|...)>0.9, let bijk=1; if E(Bijk|...)<0.1, let bijk=0; if 0.1<=E(Bijk|...)<=0.9, calculate in the summation
    ind=0;
    ind1=0;
    for (k=0; k<p; ++k){
      if (k!=j){
        indk=indexk[k];
        jk=indexjk(j,k)-1;
        
        // extract index needed in the calculation
        if (EBji(i,indk)==1) {
          
          probi[ind]=prob(i,indk);
          XBMi[ind]=XBM(i,indk);
          
          ++ind;
          
        } // end of EBji
        
        // E(Bijk)>0.9
        if (EBj1(i,indk)==1) {
          probi1[ind1]=prob(i,indk);
          XBMi1[ind1]=XBM(i,indk);
          
          ++ind1;
          
        } // end of Bj
        
        
      } // end of k!=j
    } // end of k
    
    
    db=Eigen::VectorXd::Zero(lenB);
    cb=Eigen::VectorXd::Zero(lenB);
    
    d1=0.0;
    c1=0.0;
    
    llk=0.0;
    
    if (lenB==0){
      
      if (lenB1>0){
        c1=XBMi1.sum();
        d1=(log(probi1.array()/(1-probi1.array()))).sum();
      } else{
        
        c1=0;
        d1=0;
      }
      
      llk=exp(-(M(i,j)-c1)*(M(i,j)-c1)/sigma(j)*0.5)*exp(d1);
      
    } else{
      
      rownum=BB.rows();
      
      for (row=0; row<rownum; ++row){
        
        db=BB.row(row).array();
        
        db=db.array()*log(probi.array()/(1-probi.array()));
        
        // Bijk*betajk*Xi*Mik
        cb=BB.row(row).array();
        cb=cb.array()*XBMi.array();
        
        if (lenB1>0){
          d1=db.sum()+(log(probi1.array()/(1-probi1.array()))).sum();
          c1=cb.sum()+XBMi1.sum();
        } else {
          // sum_{k!=j} Bijk*log(pijk/(1-pijk))
          d1=db.sum();
          // sum_{k!=j} Bijk*betajk*Xi*Mik
          c1=cb.sum(); 
          
        }
        
        // sum over Bijk in {0,1}
        llk=llk+exp(-(M(i,j)-c1)*(M(i,j)-c1)/sigma(j)*0.5)*exp(d1);
      }  // end of row
      
    } // end of else
    
    Likelihoodi[i]=(prodprobj[i])*llk/sqrt(2*M_PI*sigma(j));
    
  }  // end of subject i
  
  // take log, loglikelihood
  Likelihoodi=(Likelihoodi.array()).log();
  
  // sum_{i:n}
  Likeli=Likelihoodi.sum();
  
  // add likelihood of p nodes together
  // sumLikeli=Likeli.sum();
  
  
  return(List::create(Named("constant")=constant, Named("Likelihoodi")=Likelihoodi, Named("Likeli")=Likeli
  ));
  
}




/********* Compute obs likelihood function appox *******/
// [[Rcpp::export]]
List LikeliObsApproxjC(Eigen::MatrixXd M,
                       Eigen::MatrixXd indexjk,
                       Eigen::MatrixXd XBM, 
                       Eigen::VectorXd sigma, Eigen::MatrixXd S, Eigen::MatrixXd distance,  double gamma, double eta,  Eigen::VectorXd lambda1, 
                       Eigen::VectorXd coeff1, int nodej){   
  
  //  int j, int k, int k2
  int N=M.rows(), p=M.cols(), pp=p*(p-1)*(p-2)/2, po=p*(p-1);
  int l,s,l1;
  int i,j,k,k2;
  int jl, jl1, jk, kj, jk2, jkk, kjk,  count, countk;  
  // int m=lambda.size();
  int p2=p-3;
  
  
  // // only need calculate for node j
  // Eigen::MatrixXd numerator=Eigen::MatrixXd::Zero(N,p-1);
  // Eigen::MatrixXd prob=Eigen::MatrixXd::Zero(N,p-1);
  
  
  Eigen::MatrixXd numerator=Eigen::MatrixXd::Zero(N*p,p);
  Eigen::MatrixXd prob0=Eigen::MatrixXd::Zero(N*p,p);
  // 
  Eigen::MatrixXd prob=Eigen::MatrixXd::Zero(N,po);
  
  // log of product of degrees: sijk=sij*sik / common neighbors in DTI
  Eigen::MatrixXd logS=Eigen::MatrixXd::Zero(N*p, p);
  logS=log(1+S.array());
  
  // log of distance 
  Eigen::MatrixXd logdist=Eigen::MatrixXd::Zero(p, p);
  logdist=log(distance.array());
  logdist.diagonal()=Eigen::VectorXd::Zero(p);
  
  Eigen::VectorXd constant;
  
  Eigen::VectorXd add;
  
  Eigen::VectorXd Likelihoodi;
  double Likeli;
  
  double sumLikeli=0.0;  
  
  Eigen::VectorXd prodprobj;
  
  Eigen::MatrixXd db, cb;
  
  
  for (i=0; i<N; ++i){
    
    for (j=0; j<(p-1); ++j){
      for (k=j+1; k<p; ++k){
        // for (k=0; k<p; ++k){
        //   if (k != j){
        // jk=indexjk(j,k)-1;
        numerator(i*p+j,k)=exp(gamma*logS(i*p+j,k)-eta*logdist(j,k));
        numerator(i*p+k,j)=numerator(i*p+j,k);
        
        prob0(i*p+j,k)=numerator(i*p+j,k)/(1+numerator(i*p+j,k));
        prob0(i*p+k,j)=prob0(i*p+j,k);
        
        //       //    } // end of k!=j
      } // end of k
    }  // end of j
    
  } // end of i
  
  for (i=0; i<N; ++i){
    s=0;
    for (j=0; j<p; ++j){
      for (k=0; k<p; ++k){
        if(k!=j){
          prob(i,s)=prob0(i*p+j,k);
          ++s;
        } // end of k!=j
      } // end of k
    } // end of j
  }  // end of i
  
  // Eigen::VectorXd indexk=Eigen::VectorXd::Zero(p);
  // 
  // j=nodej;
  // 
  // s=0;
  // for (k=0; k<p; ++k){
  //   if (k != j){
  //     // jk=indexjk(j,k)-1;
  //     
  //     for (i=0; i<N; ++i){
  //       
  //       numerator(i,s)=exp(gamma*logS(i*p+j,k)-eta*logdist(j,k));
  //       //  numerator(i*p+k,j)=numerator(i*p+j,k);
  //       
  //       prob(i,s)=numerator(i,s)/(1+numerator(i,s));
  //       //  prob0(i*p+k,j)=prob0(i*p+j,k);
  //       
  //     } // end of subject i
  //     
  //     indexk[k]=s;
  //     
  //     s=s+1;
  //     
  //   } // end of k!=j
  //   
  // } // end of k
  
  
  
  // compute likelihood function for observed data: P(Mij, Bijk, k!=j | Mik, k!=j) 
  // Likelihoodi=Eigen::MatrixXd::Zero(N, p);
  // Likeli=Eigen::VectorXd::Zero(p);
  
  Likelihoodi=Eigen::VectorXd::Zero(N);
  Likeli=0.0;
  
  j=nodej;
  // for (j=0; j<p; ++j){
  // prod_{k!=j}(1-pijk)
  prodprobj=Eigen::VectorXd::Ones(N);
  
  for (l=0; l<p; ++l){
    if (l!=j){
      jl=indexjk(j,l)-1;
      prodprobj=(prodprobj.array())*(1-prob.col(jl).array());
    } // end of l!=j
  } // end of l
  
  // for (s=0; s<(p-1); ++s){
  //   prodprobj=(prodprobj.array())*(1-prob.col(s).array());
  // } // end of s
  // 
  
  db=Eigen::MatrixXd::Zero(N, p-1);
  cb=Eigen::MatrixXd::Zero(N, p-1);
  
  count=0;
  // s=0;
  for (k=0; k<p; ++k){
    if (k!=j){
      
      jk=indexjk(j,k)-1;
      
      db.col(count)=log(prob.col(jk).array()/(1-prob.col(jk).array()));
      
      cb.col(count)=(XBM.col(jk).array())/sqrt(2*sigma(j));
      
      db.col(count)=db.col(count).array()+(M.col(j).array())*(XBM.col(jk).array())/sigma(j);
      
      ++count;
      
      // s=s+1;
      
    } // end of k != j
  } // end of k
  
  constant=Eigen::VectorXd::Zero(N);
  
  add=Eigen::VectorXd::Ones(N);
  
  //  Eigen::VectorXd eMbeta=exp(-(M.col(j).array())*(M.col(j).array())/sigma(j)*0.5);
  
  constant=linearApproxC(db, cb, lambda1, coeff1);
  
  // for (i=0; i<N; ++i){
  //    if (constant(i)<0) constant(i)=0;
  // }
  
  Likelihoodi=(prodprobj.array())*exp(-(M.col(j).array())*(M.col(j).array())/sigma(j)*0.5)/sqrt(2*M_PI*sigma(j));
  
  // add 1 in order to avoid taking log 
  //  Likelihoodi=(Likelihoodi.array())*(constant.array())+add.array();
  Likelihoodi=(Likelihoodi.array())*(constant.array());
  
  //  } // end of j
  
  // take log, loglikelihood, probably not take log
  // Likelihoodi=(Likelihoodi.array()).log();
  
  // Likeli=Likelihoodi.sum();
  
  
  // use likelihood directly
  Likeli=Likelihoodi.prod();
  
  // // sum_{i:n}log(f(Mij | Mik, Bijk, k!=j))
  // Likeli=Likelihoodi.colwise().sum();
  // 
  // // add likelihood of p nodes together
  // sumLikeli=Likeli.sum();
  
  
  return(List::create(Named("constant")=constant, Named("Likelihoodi")=Likelihoodi, Named("Likeli")=Likeli,
                            Named("db")=db, Named("cb")=cb
  ));
  
}



/**** Update gamma and eta by optimize objective function  ****/
// [[Rcpp::export]]
List dream3DC(Eigen::MatrixXd EB, Eigen::MatrixXd S, Eigen::MatrixXd distance, Eigen::MatrixXd indexjk,
              Eigen::VectorXd gammas, Eigen::VectorXd etas, int N, int p){
  
  
  Eigen::MatrixXd numerator=Eigen::MatrixXd::Zero(N*p,p);
  Eigen::MatrixXd prob=Eigen::MatrixXd::Zero(N*p,p);
  
  // number of gamma eta combination
  int len=gammas.size()*etas.size();
  
  // sum_{i=1:N}sum_{k!=j}log(1-pijk)
  Eigen::MatrixXd Q1=Eigen::MatrixXd::Zero(len,p);
  
  // sum_{i=1:N}sum_{k!=j}log(pijk/(1-pijk))E(Bijk|Mi)
  Eigen::MatrixXd Q2=Eigen::MatrixXd::Zero(len,p);
  
  // Q1+Q2
  Eigen::MatrixXd Q=Eigen::MatrixXd::Zero(len,p);
  
  
  int i,j,k, jk, s,r, sr;
  
  // log of product of degrees: sijk=sij*sik / common neighbors in DTI
  Eigen::MatrixXd logS=Eigen::MatrixXd::Zero(N*p, p);
  logS=log(1+S.array());
  
  // log of distance 
  Eigen::MatrixXd logdist=Eigen::MatrixXd::Zero(p, p);
  logdist=log(distance.array());
  logdist.diagonal()=Eigen::VectorXd::Zero(p);
  
  Eigen::VectorXd Qsr;
  
  Eigen::MatrixXd probi;
  
  // log(1-pijk)
  Eigen::MatrixXd lognprobi;
  
  //log(pijk/(1-pijk))
  Eigen::MatrixXd logprobi;
  
  // log(pijk/(1-pijk))*E(Bijk|Mi)
  Eigen::MatrixXd Eprobi;
  
  
  
  sr=0;
  for (s=0; s<gammas.size(); ++s){
    for (r=0; r<etas.size(); ++r){
      
      Qsr=Eigen::VectorXd::Zero(p);
      
      for (i=0; i<N; ++i){
        
        // EB for subject i
        Eigen::MatrixXd EBi=Eigen::MatrixXd::Zero(p, p);
        
        
        for (j=0; j<p; ++j){
          for (k=0; k<p; ++k){
            if (k != j){
              
              jk=indexjk(j,k)-1;
              EBi(j,k)=EB(i,jk);
              
              numerator(i*p+j,k)=exp(gammas[s]*logS(i*p+j,k)-etas[r]*logdist(j,k));
            } // end of k!=j
            
          } // end of k
        }  // end of j
        
        
        // prob for subject i
        probi=Eigen::MatrixXd::Zero(p, p);
        
        // log(1-pijk)
        lognprobi=Eigen::MatrixXd::Zero(p, p);
        
        //log(pijk/(1-pijk))
        logprobi=Eigen::MatrixXd::Zero(p, p);
        
        // log(pijk/(1-pijk))*E(Bijk|Mi)
        Eprobi=Eigen::MatrixXd::Zero(p, p);
        
        
        
        for (j=0; j<p; ++j){
          
          for (k=j+1; k<p; ++k){
            prob(i*p+j,k)=numerator(i*p+j,k)/(1+(numerator.row(i*p+j)).sum()+(numerator.row(i*p+k)).sum());
            prob(i*p+k,j)=prob(i*p+j,k);
            
            
            probi(j,k)=prob(i*p+j,k);
            probi(k,j)=probi(j,k);
            
            
          } // end of k
          
        } // end of j
        
        lognprobi=log(1-probi.array());
        lognprobi.diagonal()=Eigen::VectorXd::Zero(p);
        
        logprobi=log(probi.array()/(1-probi.array()));
        logprobi.diagonal()=Eigen::VectorXd::Zero(p);
        
        Eprobi=logprobi.array()*EBi.array();
        
        // .rowwise().sum(): sum_k!=j
        Qsr=Qsr+lognprobi.rowwise().sum()+Eprobi.rowwise().sum();
        
        //  Q1.row(sr)=Q1.row(sr)+lognprobi.rowwise().sum();
        //  Q2.row(sr)=Q2.row(sr)+Eprobi.rowwise().sum();
        
      } // end of i
      
      // want to maximize Q
      Q.row(sr)=Qsr;
      
      
      ++sr;
      
    }  // end of r
  } // end of s
  
  // Eigen::VectorXd a=lognprobi.rowwise().sum();
  // Eigen::VectorXd b=Eprobi.rowwise().sum();
  
  
  
  return(List::create(Named("Q")=Q, Named("Q1")=Q1, Named("Q2")=Q2,
                      Named("numerator")=numerator, Named("Qsr")=Qsr, Named("Eprobi")=Eprobi, Named("lognprobi")=lognprobi,
                      Named("prob")=prob));
}



/***** newton-R
 aphson estimate eta and gamma  ****/

// [[Rcpp::export]]
List dream1stepC(Eigen::MatrixXd EB, Eigen::MatrixXd S, Eigen::MatrixXd distance, Eigen::MatrixXd indexjk,
                 double inigamma, double inieta, int N, int p,  int nodej){
  
  // EB=matrix(N, po)
  
  // B: true connectivity in sMRI, (N*p,p) matrix
  // isActiveC: connectivity in DTI (0/1), (N*p,p) matrix
  // distance: each subject has the same distance between nodes, (p,p) matrix
  // degree: degree of each node for each subject (N,p) matrix; isActiveC: whether the node l belongs to Cij ;
  // N: number of subjects; p: number of nodes
  
  // product of degrees: sijk=sij*sik / common neighbors in DTI
  // Eigen::MatrixXd S=Eigen::MatrixXd::Zero(N*p, p);
  // // degree for the ith subject
  // Eigen::MatrixXd Si=Eigen::MatrixXd::Zero(p, p);
  
  int i,j,k,l, jk, jl;
  
  
  // log of product of degrees: sijk=sij*sik / common neighbors in DTI
  Eigen::MatrixXd logS=Eigen::MatrixXd::Zero(N*p, p);
  logS=log(1+S.array());
  
  // log of distance 
  Eigen::MatrixXd logdist=Eigen::MatrixXd::Zero(p, p);
  logdist=log(distance.array());
  logdist.diagonal()=Eigen::VectorXd::Zero(p);
  
  
  // initial:
  
  double gamma=inigamma;
  double eta=inieta;
  double d1=0.0, d2=0.0, d3=0.0, d4=0.0, d5=0.0, d6=0.0;
  // Eigen::VectorXd sumActiveC=Eigen::VectorXd::Zero(N*q);
  
  // |Cij|: number of links for node j of subject i;
  // sumActiveC=isActiveC.rowwise().sum();
  
  double gamma1=0.0, eta1=0.0;
  
  // for gamma stop criteria
  double diffstop=0.0;
  double denominator=0.0;
  
  // // for eta stop criteria
  // double diffstope=0.0;
  // double denominatore=0.0;
  
  
  
  // first derivatives and second derivatives;
  double dgamma=0.0, deta=0.0, dgamma2=0.0, deta2=0.0, dgammaeta=0.0;
  
  
  //   int it;
  // 1-step iteration
  //  for (it=0; it<5; it++) {  
  double sumBs=0.0;
  double sumBd=0.0;
  
  double sumC2=0.0;
  double sumC3=0.0;
  double sumC4=0.0;
  double sumC5=0.0;
  double sumC6=0.0;
  
  
  // Eigen::VectorXd sumBjs=Eigen::VectorXd::Zero(p);
  // Eigen::VectorXd sumBjd=Eigen::VectorXd::Zero(p);
  // 
  // 
  // Eigen::VectorXd sumCj2=Eigen::VectorXd::Zero(p);
  // Eigen::VectorXd sumCj3=Eigen::VectorXd::Zero(p);
  // Eigen::VectorXd sumCj4=Eigen::VectorXd::Zero(p);
  // Eigen::VectorXd sumCj5=Eigen::VectorXd::Zero(p);
  // Eigen::VectorXd sumCj6=Eigen::VectorXd::Zero(p);
  
  
  //   i=5;
  for (i=0; i<N; ++i){
    
    // for each subject i
    //   Eigen::VectorXd sumD1=Eigen::VectorXd::Zero(p), sumD2=Eigen::VectorXd::Zero(p), sumD3=Eigen::VectorXd::Zero(p), sumD4=Eigen::VectorXd::Zero(p), sumD5=Eigen::VectorXd::Zero(p), sumD6=Eigen::VectorXd::Zero(p);
    
    Eigen::MatrixXd Bs=Eigen::MatrixXd::Zero(p, p), Bd=Eigen::MatrixXd::Zero(p, p), B2=Eigen::MatrixXd::Zero(p, p),B3=Eigen::MatrixXd::Zero(p, p),B4=Eigen::MatrixXd::Zero(p, p),B5=Eigen::MatrixXd::Zero(p, p),B6=Eigen::MatrixXd::Zero(p, p);
    
    //     Eigen::MatrixXd C1=Eigen::MatrixXd::Zero(p, p);
    Eigen::MatrixXd C2=Eigen::MatrixXd::Zero(p, p);
    Eigen::MatrixXd C3=Eigen::MatrixXd::Zero(p, p);
    Eigen::MatrixXd C4=Eigen::MatrixXd::Zero(p, p);
    Eigen::MatrixXd C5=Eigen::MatrixXd::Zero(p, p);
    Eigen::MatrixXd C6=Eigen::MatrixXd::Zero(p, p);
    
    //    Eigen::MatrixXd sumActiveCi=Eigen::MatrixXd::Zero(q, q);
    
    j=nodej;
    //  for (j=0; j<p; ++j) {
    for (l=0; l<p; ++l){
      // isActiveC: DTI connected set; only calculate the upper diagonal part, symmetric
      if (l!=j){
        
        d1=exp(gamma*logS(i*p+j,l)-eta*logdist(j,l));
        d2=logS(i*p+j,l)*d1;
        d3=logdist(j,l)*d1;
        d4=logS(i*p+j,l)*d2;
        d5=logdist(j,l)*d3;
        d6=logdist(j,l)*d2;
        
        jl=indexjk(j,l)-1;
        Bs(j,l)=EB(i,jl)*logS(i*p+j,l);
        
        Bd(j,l)=EB(i,jl)*logdist(j,l);
        
        C2(j,l)=d2/(d1+1);
        
        C3(j,l)=d3/(d1+1);
        
        C4(j,l)=d4/(d1+1)/(d1+1);
        
        C5(j,l)=d5/(d1+1)/(d1+1);
        
        C6(j,l)=d6/(d1+1)/(d1+1);
        
      } // end of l != j
    } // end of l
    //    } // end of j
    
    
    
    // for each subject, calculate row sum of B1, sum_{l!=j}, since diagonal=0; sum over subject i=1:N
    //   sumBjs=sumBjs+Bs.rowwise().sum();
    
    sumBs=sumBs+(Bs.row(j)).sum();  
    
    //   sumBjd=sumBjd+Bd.rowwise().sum();
    
    sumBd=sumBd+(Bd.row(j)).sum();
    
    // C2.rowwise().sum(): p-dimensional vector, sum over l!=j
    // sumCj2: p-dimensional vector
    // sum over subject i=1:N
    // sumCj2=sumCj2+C2.rowwise().sum();
    // sumCj3=sumCj3+C3.rowwise().sum();
    // sumCj4=sumCj4+C4.rowwise().sum();
    // sumCj5=sumCj5+C5.rowwise().sum();
    // sumCj6=sumCj6+C6.rowwise().sum();
    
    
    sumC2=sumC2+(C2.row(j)).sum();
    sumC3=sumC3+(C3.row(j)).sum();
    sumC4=sumC4+(C4.row(j)).sum();
    sumC5=sumC5+(C5.row(j)).sum();
    sumC6=sumC6+(C6.row(j)).sum();
    
    
    
  } // end of subject i
  
  
  // // sum over j=1:p
  // sumBs=sumBjs.sum();
  // sumBd=sumBjd.sum();
  // 
  // sumC2=sumCj2.sum();
  // sumC3=sumCj3.sum();
  // sumC4=sumCj4.sum();
  // sumC5=sumCj5.sum();
  // sumC6=sumCj6.sum();
  
  // first derivative respect to gamma;
  dgamma=(-sumBs+sumC2)/N;
  
  deta=(sumBd-sumC3)/N;
  
  dgamma2=sumC4/N;
  
  deta2=sumC5/N;
  
  dgammaeta=-sumC6/N;
  
  double dd=0.0;
  
  dd=dgamma2*deta2-dgammaeta*dgammaeta;
  
  
  // use newton-ralphson update gamma and eta;
  
  // difference between gamma_j(n+1) and gamma_j(n);
  double difgamma=0.0;
  difgamma=(dgamma*deta2-deta*dgammaeta)/dd;
  
  // difference between eta_j(n+1) and eta_j(n); 
  double difeta=0.0;
  difeta=(-dgamma*dgammaeta+dgamma2*deta)/dd;
  
  //  difeta=0.0;
  
  // update gamma and eta
  
  
  gamma1=gamma-difgamma;
  eta1=eta-difeta;
  
  gamma=gamma1;
  eta=eta1;
  
  //  } // end of it
  
  
  
  return(List::create(Named("gamma")=gamma, Named("eta")=eta, Named("DegreeS")=S,
                      Named("gamma1")=gamma1, Named("eta1")=eta1,
                      Named("logS")=logS, Named("logdist")=logdist,
                      //      Named("difgamma")=difgamma,
                      //      Named("difeta")=difeta ,
                      Named("denominator")=denominator, 
                      //Named("denominator")=denominator,
                      Named("dgamma")=dgamma, Named("deta")=deta,
                      Named("dgamma2")=dgamma2, Named("deta2")=deta2,
                      Named("dgammaeta")=dgammaeta
                        //  ,Named("gamma1s")=gamma1s, Named("eta1s")=eta1s
                        // Named("sumBs")=sumBs, Named("sumB2")=sumB2,
                        //          ,Named("sumC2")=sumC2
  ));
  
  // , Named("iterition")=it
}




/***** M-step: update beta and sigma^2 ****/
// [[Rcpp::export]]
List MstepjC(Eigen::MatrixXd X, Eigen::MatrixXd M,  Eigen::MatrixXd EB, Eigen::MatrixXd EBB, 
             Eigen::MatrixXd indexjk, Eigen::MatrixXd indexjkk, int nodej){ 
  // solve Beta directly
  int  N=X.rows(), p=M.cols(), q=X.cols(), po=p*(p-1);
  
  int  i, j, k, l, jk, s, r, sr, jl, jkl, jkk, it=0, il, iadd, ia=0; 
  
  int countk, countl;
  
  Eigen::MatrixXd beta=Eigen::MatrixXd::Zero(q*(p-1),p);
  Eigen::MatrixXd Beta=Eigen::MatrixXd::Zero(q,p*(p-1));
  
  
  Eigen::VectorXd Sigma=Eigen::VectorXd::Ones(p);
  
  Eigen::VectorXd MEB;
  Eigen::VectorXd MEBXv;
  
  Eigen::MatrixXd MEBX;
  Eigen::MatrixXd MEBXX;
  
  Eigen::VectorXd MEBBXXdiag;
  
  Eigen::MatrixXd MMEB;
  // right hand: MMEBX
  Eigen::VectorXd RightHand;
  // right hand matrix
  Eigen::MatrixXd RightHandM;
  
  Eigen::VectorXd betaj;
  
  Eigen::MatrixXd betajX;
  Eigen::MatrixXd EMMj;
  
  double EMMBetaxj; 
  double EBBMMBetaxj;
  
  Eigen::VectorXd betajXk; 
  Eigen::VectorXd betajXl;
  Eigen::VectorXd EMMjk;
  
  
  Eigen::VectorXd mmebj;
  
  
  // M^2
  Eigen::MatrixXd M2=M.array().square();
  
  // XX
  Eigen::MatrixXd XX=Eigen::MatrixXd::Zero(N,q*(q+1)/2);
  
  sr=0;
  for (s=0; s<q; ++s){
    for (r=s; r<q; ++r){
      XX.col(sr)=X.col(s).array()*X.col(r).array();
      ++sr;
    }
  }
  
  // update Beta, update for each beta_{jk} separately, k varies
  
  j=nodej;
  //  for (j=0; j<p; ++j){
  // j=1;
  //  MEB=Eigen::MatrixXd::Zero(N,(p-1)*p/2);
  
  
  MEBXX=Eigen::MatrixXd::Zero(q*(p-1),q*(p-1));
  MEBBXXdiag=Eigen::VectorXd::Zero(q*(p-1));
  
  // right hand: MMEBX
  RightHand=Eigen::VectorXd::Zero(q*(p-1));
  // right hand matrix
  RightHandM=Eigen::MatrixXd::Zero(p-1,q);
  MMEB=Eigen::MatrixXd::Zero(N,p-1);
  
  countk=0;
  for (k=0; k<p; ++k){
    if (k != j){
      jk=indexjk(j,k)-1;
      
      MMEB.col(countk)=M.col(j).array()*M.col(k).array()*EB.col(jk).array();
      
      countl=countk;
      for (l=k; l<p; ++l){
        if (l != j){
          MEB=Eigen::VectorXd::Zero(N);
          MEBXv=Eigen::VectorXd::Zero(q*(q+1)*0.5);
          // turn vector to q*q matrix
          
          MEBX=Eigen::MatrixXd::Zero(q,q);
          
          if (l==k){
            MEB=M.col(k).array()*M.col(k).array()*EB.col(jk).array();
          }
          if (l>k) {
            jkl=indexjkk(j*p+k,l)-1;
            MEB=M.col(k).array()*M.col(l).array()*EBB.col(jkl).array();
          }
          
          // a q(q+1)/2 dimensional vector
          MEBXv=MEB.transpose()*XX;
          
          // turn the vector to a q*q matrix
          sr=0;
          for (s=0; s<q; ++s){
            for (r=s; r<q; ++r){
              MEBX(s,r)=MEBXv[sr];
              MEBX(r,s)=MEBXv[sr];
              ++sr;
            }
          }
          
          // block of size (p,q), starting at (i,j) matrix.block(i,j,p,q);
          MEBXX.block(countk*q,countl*q,q,q)=MEBX;
          MEBXX.block(countl*q,countk*q,q,q)=MEBX;
          
          ++countl;
        } // end of l!=j
        
      } // end of l
      ++countk;
    } // end of k!=j
  } // end of k
  
  // (p-1)*q matrix
  RightHandM=MMEB.transpose()*X; // ?? check
  
  countk=0;
  for (k=0; k<p; ++k){
    if (k != j){
      RightHand.segment(countk*q,q)=RightHandM.row(countk);
      ++countk;
    }
  }
  
  
  
  // left hand
  // MEBBXXdiag=MEBXX.diagonal();
  // MEBXX.diagonal()=Eigen::VectorXd::Zero(q*(p-1));
  // MEBXX=MEBXX+MEBXX.transpose();
  // MEBXX.diagonal()=MEBBXXdiag;
  
  // solve Beta.col(j), a q*(p-1) dimensional vector
  betaj=MEBXX.colPivHouseholderQr().solve(RightHand);
  beta.col(j)=betaj;
  
  
  
  countk=0;
  for (k=0; k<p; ++k){
    if (k != j){
      // vector.segment(i,n): Block containing n elements, starting at position i
      Beta.col((p-1)*j+countk)=betaj.segment(countk*q,q);
      ++countk;
    }
  }
  
  
  ///////////////////////////////////////////////  
  ////////////     solve sigma(j)    ///////////
  //////////////////////////////////////////////
  
  betajX=Eigen::MatrixXd::Zero(N,p-1);
  
  betajXk=Eigen::VectorXd::Zero(p-1);
  
  EMMj=Eigen::MatrixXd::Zero(N,p-1);
  EMMjk=Eigen::VectorXd::Zero(p-1);
  
  betajXl=Eigen::VectorXd::Zero(p-1);
  
  EMMBetaxj=0.0;  
  EBBMMBetaxj=0.0;
  
  
  // sumprobj=Eigen::VectorXd::Zero(N);
  // sumBprobj=Eigen::VectorXd::Zero(N);
  
  countk=0;  
  for (k=0; k<p; ++k){
    if (k != j){
      jk=indexjk(j,k)-1;
      
      
      // // sum over k!=j; sum log(1-pijk)
      // sumprobj=sumprobj.array()+log(1-prob.col(jk).array());
      // 
      // // sum over k!=j; sum E(Bijk)*log(pijk/1-pijk)
      // sumBprobj=sumBprobj.array()+log(prob.col(jk).array()/(1-prob.col(jk).array()))*EB.col(jk).array();
      
      // betajX.col(countk)=X*betaj(countk*q,q);
      betajXk=X*betaj.segment(countk*q,q);
      EMMjk=EB.col(jk).array()*M.col(j).array()*M.col(k).array();
      
      // first part
      EMMBetaxj=EMMBetaxj+EMMjk.dot(betajXk);
      
      mmebj=Eigen::VectorXd::Zero(N);
      
      countl=0;
      for (l=0; l<p; ++l){
        if (l!=j){
          jl=indexjk(j,l)-1;
          
          if (l==k){
            betajXl=betajXk;
            betajXl=betajXl.array()*betajXk.array();
            // second part
            mmebj=M.col(k).array()*M.col(k).array()*EB.col(jk).array();
            EBBMMBetaxj=EBBMMBetaxj+betajXl.dot(mmebj);
          }
          
          if (l>k){
            betajXl=X*betaj.segment(countl*q,q);
            jkl=indexjkk(j*p+k,l)-1;
            betajXl=betajXl.array()*betajXk.array();
            
            mmebj=M.col(k).array()*M.col(l).array()*EBB.col(jkl).array();
            EBBMMBetaxj=EBBMMBetaxj+(betajXl).dot(mmebj);
          }
          
          if (l<k){
            betajXl=X*betaj.segment(countl*q,q);
            jkl=indexjkk(j*p+l,k)-1;
            betajXl=betajXl.array()*betajXk.array();
            
            mmebj=M.col(k).array()*M.col(l).array()*EBB.col(jkl).array();
            EBBMMBetaxj=EBBMMBetaxj+(betajXl).dot(mmebj);
          }
          ++countl;
        } // end of l!=j
        
      } // end of l
      ++countk;
    } // end of k!=j
  } // end of k
  
  // sum over subjects 
  // sumprob(j)=sumprobj.sum();
  // sumBprob(j)=sumBprobj.sum();
  
  
  //  Sigma[j]=M.col(j).squaredNorm()/N; 
  Sigma[j]=M.col(j).squaredNorm()/N-(2*EMMBetaxj-EBBMMBetaxj)/N;
  
  //   MB(j)=M.col(j).squaredNorm()/N;
  
  // logsigma(j)=0.5*(log(2*M_PI*Sigma(j)));  
  // 
  // Ellj(j)=-logsigma(j)+sumprob(j)/N+sumBprob(j)/N-0.5;
  
  
  // } // end of j
  
  // objective function
  //   double Obj=0.0;
  //   
  //   // // sum_{j=1}^p log(2*pi*sigma_j)
  //   // double logsigma=0.0;
  //   // logsigma=0.5*(log(2*M_PI*Sigma.array())).sum();
  //   
  // //  Obj=(log(prob.array()/(1-prob.array()))*EB.array()).sum()/N;
  // //
  //   Obj=logsigma-(log(1-prob.array())).sum()/N-(log(prob.array()/(1-prob.array()))*EB.array()).sum()/N+p*0.5;
  //   
  
  return(List::create(Named("Beta")=Beta,Named("beta")=beta, Named("Sigma")=Sigma
                        //, 
                        //Named("RightHandM")=RightHandM, Named("RightHand")=RightHand, Named("EMMBetaxj")=EMMBetaxj, Named("EBBMMBetaxj")=EBBMMBetaxj
                        // , Named("EMMBetaxj")=EMMBetaxj, Named("EBBMMBetaxj")=EBBMMBetaxj ,Named("MEBX")=MEBX, Named("MEBXv")=MEBXv, Named("XX")=XX, Named("countk")=countk, Named("countl")=countl,
                        //   Named("j")=j,Named("k")=k,Named("l")=l, Named("betajXk")=betajXk, Named("MEBXX")=MEBXX, Named("MMEB")=MMEB
                        //    Named("Ellj")=Ellj,Named("logsigma")=logsigma
  ));
}

//,Named("countk")=countk,Named("countl")=countl



/***** M-step: update beta and sigma^2 ****/
// [[Rcpp::export]]
List MstepRefitjC(Eigen::MatrixXd X, Eigen::MatrixXd M,  Eigen::MatrixXd EB, Eigen::MatrixXd EBB, 
                  Eigen::MatrixXd indexjk, Eigen::MatrixXd indexjkk, Eigen::MatrixXd isActive, int pA,int nodej){ 
  // solve Beta directly
  // pA: number of selected edges for nodej
  int  N=X.rows(), p=M.cols(), q=X.cols(), po=p*(p-1);
  
  int  i, j, k, l, jk, s, r, sr, jl, jkl, jkk, it=0, il, iadd, ia=0; 
  
  int countk, countl;
  
  double sigmaj=1;
  
  Eigen::VectorXd MEB;
  Eigen::VectorXd MEBXv;
  
  Eigen::MatrixXd MEBX;
  Eigen::MatrixXd MEBXX;
  
  Eigen::VectorXd MEBBXXdiag;
  
  Eigen::MatrixXd MMEB;
  // right hand: MMEBX
  Eigen::VectorXd RightHand;
  // right hand matrix
  Eigen::MatrixXd RightHandM;
  
  Eigen::VectorXd betaj;
  
  Eigen::MatrixXd betajX;
  Eigen::MatrixXd EMMj;
  
  double EMMBetaxj; 
  double EBBMMBetaxj;
  
  Eigen::VectorXd betajXk; 
  Eigen::VectorXd betajXl;
  Eigen::VectorXd EMMjk;
  
  
  Eigen::VectorXd mmebj;
  
  
  // M^2
  Eigen::MatrixXd M2=M.array().square();
  
  // XX
  Eigen::MatrixXd XX=Eigen::MatrixXd::Zero(N,q*(q+1)/2);
  
  sr=0;
  for (s=0; s<q; ++s){
    for (r=s; r<q; ++r){
      XX.col(sr)=X.col(s).array()*X.col(r).array();
      ++sr;
    }
  }
  
  // update Beta, update for each beta_{jk} separately, k varies
  
  j=nodej;
  //  for (j=0; j<p; ++j){
  // j=1;
  //  MEB=Eigen::MatrixXd::Zero(N,(p-1)*p/2);
  
  
  MEBXX=Eigen::MatrixXd::Zero(q*pA,q*pA);
  MEBBXXdiag=Eigen::VectorXd::Zero(q*pA);
  
  // right hand: MMEBX
  RightHand=Eigen::VectorXd::Zero(q*pA);
  // right hand matrix
  RightHandM=Eigen::MatrixXd::Zero(pA,q);
  
  MMEB=Eigen::MatrixXd::Zero(N,pA);
  
  countk=0;
  for (k=0; k<p; ++k){
    if (k != j){
      jk=indexjk(j,k)-1;
      
      // only calculate for the ones selected 
      if (isActive(j,k) == 1 ){
        
        MMEB.col(countk)=M.col(j).array()*M.col(k).array()*EB.col(jk).array();
        
        countl=countk;
        for (l=k; l<p; ++l){
          if (l != j){
            // jl=indexjk(j,l)-1;
            
            // only calculate for the ones selected , node l should be connected to node j
            if (isActive(j,l) == 1){
              MEB=Eigen::VectorXd::Zero(N);
              MEBXv=Eigen::VectorXd::Zero(q*(q+1)*0.5);
              // turn vector to q*q matrix
              
              MEBX=Eigen::MatrixXd::Zero(q,q);
              
              if (l==k){
                MEB=M.col(k).array()*M.col(k).array()*EB.col(jk).array();
              }
              if (l>k) {
                jkl=indexjkk(j*p+k,l)-1;
                MEB=M.col(k).array()*M.col(l).array()*EBB.col(jkl).array();
              }
              
              // a q(q+1)/2 dimensional vector
              MEBXv=MEB.transpose()*XX;
              
              // turn the vector to a q*q matrix
              sr=0;
              for (s=0; s<q; ++s){
                for (r=s; r<q; ++r){
                  MEBX(s,r)=MEBXv[sr];
                  MEBX(r,s)=MEBXv[sr];
                  ++sr;
                }
              }
              
              // block of size (p,q), starting at (i,j) matrix.block(i,j,p,q);
              MEBXX.block(countk*q,countl*q,q,q)=MEBX;
              MEBXX.block(countl*q,countk*q,q,q)=MEBX;
              
              ++countl;
              
            } // end of isActive_jl
            
          } // end of l!=j
          
        } // end of l
        ++countk;
        
      }  // end of isActive_jk
      
    }  // end of k!=j
    
  } // end of k
  
  // (p-1)*q matrix
  RightHandM=MMEB.transpose()*X; // ?? check
  
  countk=0;
  for (k=0; k<p; ++k){
    if (k != j){
      jk=indexjk(j,k)-1;
      
      // only calculate for the ones selected 
      if (isActive(j,k) == 1 ){
        
        RightHand.segment(countk*q,q)=RightHandM.row(countk);
        ++countk;
        
      } // end of isActive_jk  
      
    } // end of k!=j
    
    
  } // end of k
  
  
  
  // left hand
  // MEBBXXdiag=MEBXX.diagonal();
  // MEBXX.diagonal()=Eigen::VectorXd::Zero(q*(p-1));
  // MEBXX=MEBXX+MEBXX.transpose();
  // MEBXX.diagonal()=MEBBXXdiag;
  
  // solve Beta.col(j), a q*(p-1) dimensional vector
  betaj=MEBXX.colPivHouseholderQr().solve(RightHand);
  
  
  ///////////////////////////////////////////////  
  ////////////     solve sigma(j)    ///////////
  //////////////////////////////////////////////
  
  betajX=Eigen::MatrixXd::Zero(N,pA);
  betajXk=Eigen::VectorXd::Zero(pA);
  EMMj=Eigen::MatrixXd::Zero(N,pA);
  EMMjk=Eigen::VectorXd::Zero(pA);
  
  betajXl=Eigen::VectorXd::Zero(pA);
  
  EMMBetaxj=0.0;  
  EBBMMBetaxj=0.0;
  
  
  // sumprobj=Eigen::VectorXd::Zero(N);
  // sumBprobj=Eigen::VectorXd::Zero(N);
  
  countk=0;  
  for (k=0; k<p; ++k){
    if (k != j){
      jk=indexjk(j,k)-1;
      
      if (isActive(j,k) !=0 ){
        betajXk=X*betaj.segment(countk*q,q);
        EMMjk=EB.col(jk).array()*M.col(j).array()*M.col(k).array();
        
        // first part
        EMMBetaxj=EMMBetaxj+EMMjk.dot(betajXk);
        
        mmebj=Eigen::VectorXd::Zero(N);
        
        countl=0;
        for (l=0; l<p; ++l){
          if (l!=j){
            jl=indexjk(j,l)-1;
            
            if (isActive(j,l) !=0 ){
              if (l==k){
                betajXl=betajXk;
                betajXl=betajXl.array()*betajXk.array();
                // second part
                mmebj=M.col(k).array()*M.col(k).array()*EB.col(jk).array();
                EBBMMBetaxj=EBBMMBetaxj+betajXl.dot(mmebj);
              }
              
              if (l>k){
                betajXl=X*betaj.segment(countl*q,q);
                jkl=indexjkk(j*p+k,l)-1;
                betajXl=betajXl.array()*betajXk.array();
                
                mmebj=M.col(k).array()*M.col(l).array()*EBB.col(jkl).array();
                EBBMMBetaxj=EBBMMBetaxj+(betajXl).dot(mmebj);
              }
              
              if (l<k){
                betajXl=X*betaj.segment(countl*q,q);
                jkl=indexjkk(j*p+l,k)-1;
                betajXl=betajXl.array()*betajXk.array();
                
                mmebj=M.col(k).array()*M.col(l).array()*EBB.col(jkl).array();
                EBBMMBetaxj=EBBMMBetaxj+(betajXl).dot(mmebj);
              }
              ++countl;
              
            } // end of isActive_jl
            
          } // end of l!=j
          
        } // end of l
        ++countk;
        
      } // end of isActive_jk
      
    } // end of k!=j
  } // end of k
  
  //  Sigma[j]=M.col(j).squaredNorm()/N; 
  sigmaj=M.col(j).squaredNorm()/N-(2*EMMBetaxj-EBBMMBetaxj)/N;
  
  
  
  
  return(List::create(Named("betaj")=betaj, Named("sigmaj")=sigmaj
                        //, 
                        //Named("RightHandM")=RightHandM, Named("RightHand")=RightHand, Named("EMMBetaxj")=EMMBetaxj, Named("EBBMMBetaxj")=EBBMMBetaxj
                        // , Named("EMMBetaxj")=EMMBetaxj, Named("EBBMMBetaxj")=EBBMMBetaxj ,Named("MEBX")=MEBX, Named("MEBXv")=MEBXv, Named("XX")=XX, Named("countk")=countk, Named("countl")=countl,
                        //   Named("j")=j,Named("k")=k,Named("l")=l, Named("betajXk")=betajXk, Named("MEBXX")=MEBXX, Named("MMEB")=MMEB
                        //    Named("Ellj")=Ellj,Named("logsigma")=logsigma
  ));
}

//,Named("countk")=countk,Named("countl")=countl


/***** sigma^2 ****/
// [[Rcpp::export]]
List SigmaRefitjC(Eigen::MatrixXd X, Eigen::MatrixXd M,  Eigen::MatrixXd EB, Eigen::MatrixXd EBB, 
                  Eigen::MatrixXd indexjk, Eigen::MatrixXd indexjkk, Eigen::MatrixXd isActive, int pA, Eigen::VectorXd beta, int nodej){ 
  // solve Beta directly
  // pA: number of selected edges for nodej
  int  N=X.rows(), p=M.cols(), q=X.cols(), po=p*(p-1);
  
  int  i, j, k, l, jk, s, r, sr, jl, jkl, jkk, it=0, il, iadd, ia=0; 
  
  int countk, countl;
  
  double sigmaj=1;
  
  Eigen::MatrixXd betajX;
  Eigen::MatrixXd EMMj;
  
  double EMMBetaxj; 
  double EBBMMBetaxj;
  
  Eigen::VectorXd betajXk; 
  Eigen::VectorXd betajXl;
  Eigen::VectorXd EMMjk;
  
  
  Eigen::VectorXd mmebj;
  
  
  // M^2
  Eigen::MatrixXd M2=M.array().square();
  
  // XX
  Eigen::MatrixXd XX=Eigen::MatrixXd::Zero(N,q*(q+1)/2);
  
  sr=0;
  for (s=0; s<q; ++s){
    for (r=s; r<q; ++r){
      XX.col(sr)=X.col(s).array()*X.col(r).array();
      ++sr;
    }
  }
  
  // update Beta, update for each beta_{jk} separately, k varies
  
  j=nodej;
  
  
  
  ///////////////////////////////////////////////  
  ////////////     solve sigma(j)    ///////////
  //////////////////////////////////////////////
  
  betajX=Eigen::MatrixXd::Zero(N,pA);
  betajXk=Eigen::VectorXd::Zero(pA);
  EMMj=Eigen::MatrixXd::Zero(N,pA);
  EMMjk=Eigen::VectorXd::Zero(pA);
  
  betajXl=Eigen::VectorXd::Zero(pA);
  
  EMMBetaxj=0.0;  
  EBBMMBetaxj=0.0;
  
  
  // sumprobj=Eigen::VectorXd::Zero(N);
  // sumBprobj=Eigen::VectorXd::Zero(N);
  
  countk=0;  
  for (k=0; k<p; ++k){
    if (k != j){
      jk=indexjk(j,k)-1;
      
      if (isActive(j,k) !=0 ){
        betajXk=X*beta.segment(countk*q,q);
        EMMjk=EB.col(jk).array()*M.col(j).array()*M.col(k).array();
        
        // first part
        EMMBetaxj=EMMBetaxj+EMMjk.dot(betajXk);
        
        mmebj=Eigen::VectorXd::Zero(N);
        
        countl=0;
        for (l=0; l<p; ++l){
          if (l!=j){
            jl=indexjk(j,l)-1;
            
            if (isActive(j,l) !=0 ){
              if (l==k){
                betajXl=betajXk;
                betajXl=betajXl.array()*betajXk.array();
                // second part
                mmebj=M.col(k).array()*M.col(k).array()*EB.col(jk).array();
                EBBMMBetaxj=EBBMMBetaxj+betajXl.dot(mmebj);
              }
              
              if (l>k){
                betajXl=X*beta.segment(countl*q,q);
                jkl=indexjkk(j*p+k,l)-1;
                betajXl=betajXl.array()*betajXk.array();
                
                mmebj=M.col(k).array()*M.col(l).array()*EBB.col(jkl).array();
                EBBMMBetaxj=EBBMMBetaxj+(betajXl).dot(mmebj);
              }
              
              if (l<k){
                betajXl=X*beta.segment(countl*q,q);
                jkl=indexjkk(j*p+l,k)-1;
                betajXl=betajXl.array()*betajXk.array();
                
                mmebj=M.col(k).array()*M.col(l).array()*EBB.col(jkl).array();
                EBBMMBetaxj=EBBMMBetaxj+(betajXl).dot(mmebj);
              }
              ++countl;
              
            } // end of isActive_jl
            
          } // end of l!=j
          
        } // end of l
        ++countk;
        
      } // end of isActive_jk
      
    } // end of k!=j
  } // end of k
  
  //  Sigma[j]=M.col(j).squaredNorm()/N; 
  sigmaj=M.col(j).squaredNorm()/N-(2*EMMBetaxj-EBBMMBetaxj)/N;
  
  
  
  
  return(List::create(Named("sigmaj")=sigmaj
                        //, 
                        //Named("RightHandM")=RightHandM, Named("RightHand")=RightHand, Named("EMMBetaxj")=EMMBetaxj, Named("EBBMMBetaxj")=EBBMMBetaxj
                        // , Named("EMMBetaxj")=EMMBetaxj, Named("EBBMMBetaxj")=EBBMMBetaxj ,Named("MEBX")=MEBX, Named("MEBXv")=MEBXv, Named("XX")=XX, Named("countk")=countk, Named("countl")=countl,
                        //   Named("j")=j,Named("k")=k,Named("l")=l, Named("betajXk")=betajXk, Named("MEBXX")=MEBXX, Named("MMEB")=MMEB
                        //    Named("Ellj")=Ellj,Named("logsigma")=logsigma
  ));
}




/***** Calculate likelihood ****/
// [[Rcpp::export]]
List LikelihoodC(Eigen::MatrixXd X, Eigen::MatrixXd M,  Eigen::MatrixXd Bi, Eigen::MatrixXd prob, Eigen::MatrixXd Beta, Eigen::VectorXd Sigma, 
                 Eigen::MatrixXd indexjk, Eigen::MatrixXd indexjkk){ 
  
  // look at likelihood for each node j instead of all nodes
  
  int  N=X.rows(), p=M.cols(), q=X.cols(), po=p*(p-1);
  
  int  i, j, k, l, jk, s, r, sr, jl; 
  
  int countk, countl;
  
  Eigen::MatrixXd betajX;
  
  double BMMBetaxj; 
  double BBMMBetaxj;
  
  Eigen::VectorXd betajXk; 
  Eigen::VectorXd betajXl;
  Eigen::VectorXd BMMjk;
  
  Eigen::VectorXd mmbbj;
  
  
  // loglikelihood function : each node j separately
  Eigen::VectorXd llj=Eigen::VectorXd::Zero(p);
  
  //  log(2*pi*sigma_j)
  Eigen::VectorXd logsigma=Eigen::VectorXd::Zero(p);
  logsigma=0.5*(log(2*M_PI*Sigma.array()));
  
  Eigen::VectorXd MB=Eigen::VectorXd::Zero(p);
  
  Eigen::VectorXd sumprob=Eigen::VectorXd::Zero(p);
  
  Eigen::VectorXd sumBprob=Eigen::VectorXd::Zero(p);
  
  // N dimensional vector for jth probj
  Eigen::VectorXd sumprobj;
  
  // N dimensional vector for jth probj
  Eigen::VectorXd sumBprobj;
  
  for (j=0; j<p; ++j){
    
    //    j=0;
    
    betajX=Eigen::MatrixXd::Zero(N,p-1);
    //   
    betajXk=Eigen::VectorXd::Zero(p-1);
    
    BMMjk=Eigen::VectorXd::Zero(p-1);
    //   
    betajXl=Eigen::VectorXd::Zero(p-1);
    
    BMMBetaxj=0.0;  
    BBMMBetaxj=0.0;
    
    
    sumprobj=Eigen::VectorXd::Zero(N);
    sumBprobj=Eigen::VectorXd::Zero(N);
    
    for (k=0; k<p; ++k){
      if (k != j){
        jk=indexjk(j,k)-1;
        
        // sum over k!=j; sum log(1-pijk)
        sumprobj=sumprobj.array()+log(1-prob.col(jk).array());
        //   sumprob(j)=sumprobj.sum();
        
        // sum over k!=j; sum Bijk*log(pijk/1-pijk)
        sumBprobj=sumBprobj.array()+log(prob.col(jk).array()/(1-prob.col(jk).array()))*Bi.col(jk).array();
        //   sumBprob(j)=sumBprobj.sum();
        
        betajXk=X*Beta.col(jk);
        BMMjk=Bi.col(jk).array()*M.col(j).array()*M.col(k).array();
        
        // first part: (Bijk*Mij*Mik*(beta_jk*Xi))
        BMMBetaxj=BMMBetaxj+BMMjk.dot(betajXk);
        
        // second part: BijkBijk'MikMik'(beta_jk*Xi)(beta_jk'*Xi)       
        mmbbj=Eigen::VectorXd::Zero(N);
        
        for (l=0; l<p; ++l){
          if (l!=j){
            jl=indexjk(j,l)-1;
            
            if (l==k){
              betajXl=betajXk;
              betajXl=betajXl.array()*betajXk.array();
              // second part
              mmbbj=M.col(k).array()*M.col(k).array()*Bi.col(jk).array();
              BBMMBetaxj=BBMMBetaxj+betajXl.dot(mmbbj);
            }else{
              
              betajXl=X*Beta.col(jl);
              betajXl=betajXl.array()*betajXk.array();
              
              mmbbj=M.col(k).array()*M.col(l).array()*Bi.col(jk).array()*Bi.col(jl).array();
              BBMMBetaxj=BBMMBetaxj+(betajXl).dot(mmbbj);
              
            }
          } // end of l!=j
          
        } // end of l
        
      } // end of k!=j
    } // end of k
    
    sumprob(j)=sumprobj.sum();
    sumBprob(j)=sumBprobj.sum();
    
    MB(j)=(M.col(j).squaredNorm()/N-(2*BMMBetaxj-BBMMBetaxj)/N)/Sigma(j);
    
    //   MB(j)=M.col(j).squaredNorm()/N;
    
    llj(j)=-logsigma(j)+sumprob(j)/N+sumBprob(j)/N-MB(j)*0.5;
    
  }  // end of j
  
  
  return(List::create(Named("Beta")=Beta,Named("Sigma")=Sigma, Named("BMMBetaxj")=BMMBetaxj, Named("BBMMBetaxj")=BBMMBetaxj,
                            Named("sumprob")=sumprob, Named("sumBprob")=sumBprob,
                            Named("llj")=llj,Named("MB")=MB));
}




/********* Compute EBIC *******/
// [[Rcpp::export]]
List HardThresholdApproxjC(Eigen::VectorXd inibetaj, double inisigmaj, Eigen::VectorXd seqbetaj, Eigen::MatrixXd betaK, Eigen::MatrixXd X, Eigen::MatrixXd M, Eigen::MatrixXd EB, Eigen::MatrixXd EBB, 
                           Eigen::MatrixXd indexjk, Eigen::MatrixXd indexjkk, Eigen::MatrixXd S, Eigen::MatrixXd distance,  double gamma, double eta,  Eigen::VectorXd lambda1, 
                           Eigen::VectorXd coeff1, int cutnum, int nodej){   
  
  //  int j, int k, int k2
  int N=X.rows(), p=M.cols(), q=X.cols(), pp=p*(p-1)*(p-2)/2, po=p*(p-1);
  int l,s,l1;
  int i,j,k,k2;
  int jl, jl1, jk, kj, jkl, jk2, jkk, kjk,  count, countk, rk, countl;  
  //  int m=lambda.size();
  // length of cutoff values
  int cut=cutnum;
  int p2=p-3;
  
  double sigmaj=1;
  
  Eigen::VectorXd sigmajs=Eigen::VectorXd::Zero(cut);
  
  Eigen::MatrixXd Beta;
  
  Eigen::MatrixXd betajs=Eigen::MatrixXd::Zero(inibetaj.size(), cut);
  
  Eigen::MatrixXd XBM=Eigen::MatrixXd::Zero(N, po);
  
  Eigen::VectorXd betark;
  
  Eigen::VectorXd constant;
  
  //  Eigen::VectorXd Likelihoodi;
  
  Eigen::MatrixXd Likelihoodi;
  Eigen::VectorXd Likeli;
  
  //  double sumLikeli=0.0;  
  
  
  
  Eigen::MatrixXd numerator=Eigen::MatrixXd::Zero(N*p,p);
  Eigen::MatrixXd prob0=Eigen::MatrixXd::Zero(N*p,p);
  
  Eigen::MatrixXd prob=Eigen::MatrixXd::Zero(N,po);
  
  // log of product of degrees: sijk=sij*sik / common neighbors in DTI
  Eigen::MatrixXd logS=Eigen::MatrixXd::Zero(N*p, p);
  logS=log(1+S.array());
  
  // log of distance 
  Eigen::MatrixXd logdist=Eigen::MatrixXd::Zero(p, p);
  logdist=log(distance.array());
  logdist.diagonal()=Eigen::VectorXd::Zero(p);
  
  Eigen::VectorXd prodprobj;
  
  
  for (i=0; i<N; ++i){
    
    for (j=0; j<(p-1); ++j){
      for (k=j+1; k<p; ++k){
        // for (k=0; k<p; ++k){
        //   if (k != j){
        // jk=indexjk(j,k)-1;
        numerator(i*p+j,k)=exp(gamma*logS(i*p+j,k)-eta*logdist(j,k));
        numerator(i*p+k,j)=numerator(i*p+j,k);
        
        prob0(i*p+j,k)=numerator(i*p+j,k)/(1+numerator(i*p+j,k));
        prob0(i*p+k,j)=prob0(i*p+j,k);
        
        //    } // end of k!=j
      } // end of k
    }  // end of j
    
    
  } // end of i
  
  
  for (i=0; i<N; ++i){
    s=0;
    for (j=0; j<p; ++j){
      for (k=0; k<p; ++k){
        if(k!=j){
          prob(i,s)=prob0(i*p+j,k);
          ++s;
        } // end of k!=j
      } // end of k
    } // end of j
  }  // end of i
  
  
  
  Eigen::VectorXd MEB;
  Eigen::VectorXd MEBXv;
  
  Eigen::MatrixXd MEBX;
  Eigen::MatrixXd MEBXX;
  
  Eigen::VectorXd MEBBXXdiag;
  
  Eigen::MatrixXd MMEB;
  // right hand: MMEBX
  Eigen::VectorXd RightHand;
  // right hand matrix
  Eigen::MatrixXd RightHandM;
  
  Eigen::VectorXd betaj;
  
  Eigen::MatrixXd betajX;
  Eigen::MatrixXd EMMj;
  
  double EMMBetaxj; 
  double EBBMMBetaxj;
  
  Eigen::VectorXd betajXk; 
  Eigen::VectorXd betajXl;
  Eigen::VectorXd EMMjk;
  
  Eigen::VectorXd mmebj;
  
  
  //  Eigen::VectorXd prodprobj;
  
  Eigen::MatrixXd db, cb;
  
  
  //  ones=Eigen::VectorXd::Ones(inibetaj.size());
  // compute likelihood function for observed data: P(Mij, Bijk, k!=j | Mik, k!=j) 
  j=nodej;
  
  prodprobj=Eigen::VectorXd::Ones(N);
  
  for (l=0; l<p; ++l){
    if (l!=j){
      jl=indexjk(j,l)-1;
      prodprobj=(prodprobj.array())*(1-prob.col(jl).array());
    } // end of l!=j
  } // end of l
  
  
  Likeli=Eigen::VectorXd::Zero(cut);
  Likelihoodi=Eigen::MatrixXd::Zero(N,cut);
  
  // rk=1;
  for (rk=0; rk<cut; ++rk){
    
    Beta=Eigen::MatrixXd::Zero(q,p*(p-1));
    
    //  Likelihoodi=Eigen::VectorXd::Zero(N);
    
    
    // hard-thresholded betaj
    betark=betaK.col(rk);
    
    betajs.col(rk)=betark;
    
    
    // update the corresponding sigma_j
    betajX=Eigen::MatrixXd::Zero(N,p-1);
    
    betajXk=Eigen::VectorXd::Zero(p-1);
    
    EMMj=Eigen::MatrixXd::Zero(N,p-1);
    EMMjk=Eigen::VectorXd::Zero(p-1);
    
    betajXl=Eigen::VectorXd::Zero(p-1);
    
    EMMBetaxj=0.0;  
    EBBMMBetaxj=0.0;
    
    countk=0;  
    for (k=0; k<p; ++k){
      if (k != j){
        jk=indexjk(j,k)-1;
        
        // vector.segment(i,n): Block containing n elements, starting at position i
        Beta.col((p-1)*j+countk)=betark.segment(countk*q,q);
        
        // betajX.col(countk)=X*betaj(countk*q,q);
        betajXk=X*betark.segment(countk*q,q);
        EMMjk=EB.col(jk).array()*M.col(j).array()*M.col(k).array();
        
        // first part
        EMMBetaxj=EMMBetaxj+EMMjk.dot(betajXk);
        
        mmebj=Eigen::VectorXd::Zero(N);
        
        countl=0;
        for (l=0; l<p; ++l){
          if (l!=j){
            jl=indexjk(j,l)-1;
            
            if (l==k){
              betajXl=betajXk;
              betajXl=betajXl.array()*betajXk.array();
              // second part
              mmebj=M.col(k).array()*M.col(k).array()*EB.col(jk).array();
              EBBMMBetaxj=EBBMMBetaxj+betajXl.dot(mmebj);
            }
            
            if (l>k){
              betajXl=X*betark.segment(countl*q,q);
              jkl=indexjkk(j*p+k,l)-1;
              betajXl=betajXl.array()*betajXk.array();
              
              mmebj=M.col(k).array()*M.col(l).array()*EBB.col(jkl).array();
              EBBMMBetaxj=EBBMMBetaxj+(betajXl).dot(mmebj);
            }
            
            if (l<k){
              betajXl=X*betark.segment(countl*q,q);
              jkl=indexjkk(j*p+l,k)-1;
              betajXl=betajXl.array()*betajXk.array();
              
              mmebj=M.col(k).array()*M.col(l).array()*EBB.col(jkl).array();
              EBBMMBetaxj=EBBMMBetaxj+(betajXl).dot(mmebj);
            }
            ++countl;
          } // end of l!=j
          
        } // end of l
        ++countk;
      } // end of k!=j
    } // end of k
    
    //  Sigma[j]=M.col(j).squaredNorm()/N; 
    sigmaj=M.col(j).squaredNorm()/N-(2*EMMBetaxj-EBBMMBetaxj)/N;
    
    sigmajs[rk]=sigmaj;
    
    
    // compute likelihood
    db=Eigen::MatrixXd::Zero(N, p-1);
    cb=Eigen::MatrixXd::Zero(N, p-1);
    
    count=0;
    for (k=0; k<p; ++k){
      if (k!=j){
        
        jk=indexjk(j,k)-1;
        
        XBM.col(jk)=(X*Beta.col(jk)).array()*(M.col(k).array());
        
        //   XBM.col(jk)=Beta.col(jk).array()*(M.col(k).array());
        
        db.col(count)=log(prob.col(jk).array()/(1-prob.col(jk).array()));
        
        cb.col(count)=(XBM.col(jk).array())/sqrt(2*sigmaj);
        
        db.col(count)=db.col(count).array()+(M.col(j).array())*(XBM.col(jk).array())/sigmaj;
        
        ++count;
        
      } // end of k != j
    } // end of k
    
    constant=Eigen::VectorXd::Zero(N);
    
    constant=linearApproxC(db, cb,  lambda1, coeff1);
    
    Likelihoodi.col(rk)=(prodprobj.array())*exp(-(M.col(j).array())*(M.col(j).array())/sigmaj*0.5)/sqrt(2*M_PI*sigmaj);
    
    Likelihoodi.col(rk)=(Likelihoodi.col(rk).array())*(constant.array());
    
    //  } // end of j
    
    // take log, loglikelihood, probably not take log
    //  Likelihoodi=(Likelihoodi.array()).log();
    
    // sum_{i:n}
    //  Likeli[rk]=Likelihoodi.sum();
    
    // prod_{1:n}, use likelihood directly instead of log-likelihood
    Likeli[rk]=Likelihoodi.col(rk).prod();
    
    
    
  }  // end of rk, hard-thresholding
  
  
  
  return(List::create(Named("constant")=constant,Named("cb")=cb,Named("db")=db, Named("XBM")=XBM, Named("Beta")=Beta, Named("Likelihoodi")=Likelihoodi, Named("Likeli")=Likeli, Named("betajs")=betajs, Named("sigmajs")=sigmajs
  ));
  
}



/********* Compute EBIC *******/
// [[Rcpp::export]]
List HardThresholdjC(Eigen::VectorXd inibetaj, double inisigmaj, Eigen::VectorXd seqbetaj, Eigen::MatrixXd betaK, Eigen::MatrixXd X, Eigen::MatrixXd M, Eigen::MatrixXd EB, Eigen::MatrixXd EBB, 
                     Eigen::MatrixXd indexjk, Eigen::MatrixXd indexjkk, Eigen::MatrixXd S, Eigen::MatrixXd distance,  double gamma, double eta, Eigen::MatrixXd B, int cutnum, int nodej){   
  
  //  int j, int k, int k2
  int N=X.rows(), p=M.cols(), q=X.cols(), pp=p*(p-1)*(p-2)/2, po=p*(p-1);
  int l,s,l1;
  int i,j,k,k2;
  int jl, jl1, jk, kj, jkl, jk2, jkk, kjk,  count, countk, rk, countl;  
  //  int m=lambda.size();
  // length of cutoff values
  int cut=cutnum;
  int p2=p-3;
  int row;
  int rownum=B.rows();
  
  double sigmaj=1;
  
  Eigen::VectorXd sigmajs=Eigen::VectorXd::Zero(cut);
  
  Eigen::MatrixXd Beta;
  
  Eigen::MatrixXd betajs=Eigen::MatrixXd::Zero(inibetaj.size(), cut);
  
  Eigen::MatrixXd XBM=Eigen::MatrixXd::Zero(N, po);
  
  Eigen::VectorXd betark;
  
  Eigen::VectorXd constant;
  
  Eigen::VectorXd Likelihoodi;
  Eigen::VectorXd Likeli;
  
  //  double sumLikeli=0.0;  
  
  Eigen::MatrixXd db, cb;
  
  Eigen::VectorXd d1, c1;
  
  Eigen::VectorXd llk;
  
  
  
  
  Eigen::MatrixXd numerator=Eigen::MatrixXd::Zero(N*p,p);
  Eigen::MatrixXd prob0=Eigen::MatrixXd::Zero(N*p,p);
  
  Eigen::MatrixXd prob=Eigen::MatrixXd::Zero(N,po);
  
  // log of product of degrees: sijk=sij*sik / common neighbors in DTI
  Eigen::MatrixXd logS=Eigen::MatrixXd::Zero(N*p, p);
  logS=log(1+S.array());
  
  // log of distance 
  Eigen::MatrixXd logdist=Eigen::MatrixXd::Zero(p, p);
  logdist=log(distance.array());
  logdist.diagonal()=Eigen::VectorXd::Zero(p);
  
  Eigen::VectorXd prodprobj;
  
  
  for (i=0; i<N; ++i){
    
    for (j=0; j<(p-1); ++j){
      for (k=j+1; k<p; ++k){
        // for (k=0; k<p; ++k){
        //   if (k != j){
        // jk=indexjk(j,k)-1;
        numerator(i*p+j,k)=exp(gamma*logS(i*p+j,k)-eta*logdist(j,k));
        numerator(i*p+k,j)=numerator(i*p+j,k);
        
        prob0(i*p+j,k)=numerator(i*p+j,k)/(1+numerator(i*p+j,k));
        prob0(i*p+k,j)=prob0(i*p+j,k);
        
        //    } // end of k!=j
      } // end of k
    }  // end of j
    
    
  } // end of i
  
  
  for (i=0; i<N; ++i){
    s=0;
    for (j=0; j<p; ++j){
      for (k=0; k<p; ++k){
        if(k!=j){
          prob(i,s)=prob0(i*p+j,k);
          ++s;
        }
      }
    }
  }  // end of i
  
  
  
  Eigen::VectorXd MEB;
  Eigen::VectorXd MEBXv;
  
  Eigen::MatrixXd MEBX;
  Eigen::MatrixXd MEBXX;
  
  Eigen::VectorXd MEBBXXdiag;
  
  Eigen::MatrixXd MMEB;
  // right hand: MMEBX
  Eigen::VectorXd RightHand;
  // right hand matrix
  Eigen::MatrixXd RightHandM;
  
  Eigen::VectorXd betaj;
  
  Eigen::MatrixXd betajX;
  Eigen::MatrixXd EMMj;
  
  double EMMBetaxj; 
  double EBBMMBetaxj;
  
  Eigen::VectorXd betajXk; 
  Eigen::VectorXd betajXl;
  Eigen::VectorXd EMMjk;
  
  Eigen::VectorXd mmebj;
  
  
  //  Eigen::VectorXd prodprobj
  
  
  //  ones=Eigen::VectorXd::Ones(inibetaj.size());
  // compute likelihood function for observed data: P(Mij, Bijk, k!=j | Mik, k!=j) 
  j=nodej;
  
  prodprobj=Eigen::VectorXd::Ones(N);
  
  for (l=0; l<p; ++l){
    if (l!=j){
      jl=indexjk(j,l)-1;
      prodprobj=(prodprobj.array())*(1-prob.col(jl).array());
    } // end of l!=j
  } // end of l
  
  
  Likeli=Eigen::VectorXd::Zero(cut);
  
  // rk=1;
  for (rk=0; rk<cut; ++rk){
    
    Beta=Eigen::MatrixXd::Zero(q,p*(p-1));
    
    Likelihoodi=Eigen::VectorXd::Zero(N);
    
    
    // hard-thresholded betaj
    betark=betaK.col(rk);
    
    betajs.col(rk)=betark;
    
    
    // update the corresponding sigma_j
    betajX=Eigen::MatrixXd::Zero(N,p-1);
    
    betajXk=Eigen::VectorXd::Zero(p-1);
    
    EMMj=Eigen::MatrixXd::Zero(N,p-1);
    EMMjk=Eigen::VectorXd::Zero(p-1);
    
    betajXl=Eigen::VectorXd::Zero(p-1);
    
    EMMBetaxj=0.0;  
    EBBMMBetaxj=0.0;
    
    countk=0;  
    for (k=0; k<p; ++k){
      if (k != j){
        jk=indexjk(j,k)-1;
        
        // vector.segment(i,n): Block containing n elements, starting at position i
        Beta.col((p-1)*j+countk)=betark.segment(countk*q,q);
        
        // betajX.col(countk)=X*betaj(countk*q,q);
        betajXk=X*betark.segment(countk*q,q);
        EMMjk=EB.col(jk).array()*M.col(j).array()*M.col(k).array();
        
        // first part
        EMMBetaxj=EMMBetaxj+EMMjk.dot(betajXk);
        
        mmebj=Eigen::VectorXd::Zero(N);
        
        countl=0;
        for (l=0; l<p; ++l){
          if (l!=j){
            jl=indexjk(j,l)-1;
            
            if (l==k){
              betajXl=betajXk;
              betajXl=betajXl.array()*betajXk.array();
              // second part
              mmebj=M.col(k).array()*M.col(k).array()*EB.col(jk).array();
              EBBMMBetaxj=EBBMMBetaxj+betajXl.dot(mmebj);
            }
            
            if (l>k){
              betajXl=X*betark.segment(countl*q,q);
              jkl=indexjkk(j*p+k,l)-1;
              betajXl=betajXl.array()*betajXk.array();
              
              mmebj=M.col(k).array()*M.col(l).array()*EBB.col(jkl).array();
              EBBMMBetaxj=EBBMMBetaxj+(betajXl).dot(mmebj);
            }
            
            if (l<k){
              betajXl=X*betark.segment(countl*q,q);
              jkl=indexjkk(j*p+l,k)-1;
              betajXl=betajXl.array()*betajXk.array();
              
              mmebj=M.col(k).array()*M.col(l).array()*EBB.col(jkl).array();
              EBBMMBetaxj=EBBMMBetaxj+(betajXl).dot(mmebj);
            }
            ++countl;
          } // end of l!=j
          
        } // end of l
        ++countk;
      } // end of k!=j
    } // end of k
    
    //  Sigma[j]=M.col(j).squaredNorm()/N; 
    sigmaj=M.col(j).squaredNorm()/N-(2*EMMBetaxj-EBBMMBetaxj)/N;
    
    sigmajs[rk]=sigmaj;
    
    
    // compute likelihood
    
    db=Eigen::MatrixXd::Zero(N, p-1);
    cb=Eigen::MatrixXd::Zero(N, p-1);
    
    d1=Eigen::VectorXd::Zero(N);
    c1=Eigen::VectorXd::Zero(N);
    
    llk=Eigen::VectorXd::Zero(N);
    
    for (row=0; row<rownum; ++row){
      
      count=0;
      for (k=0; k<p; ++k){
        if (k!=j){
          
          jk=indexjk(j,k)-1;
          
          XBM.col(jk)=(X*Beta.col(jk)).array()*(M.col(k).array());
          
          // Bijk*log(pijk/(1-pijk))
          db.col(count)=B(row, count)*log(prob.col(jk).array()/(1-prob.col(jk).array()));
          
          // Bijk*betajk*Xi*Mik
          cb.col(count)=B(row,count)*XBM.col(jk).array();
          
          //  db.col(count)=db.col(count).array()+(M.col(j).array())*(XBM.col(jk).array())/sigma(j);
          
          // /sqrt(2*sigma(j))
          
          ++count;
          
        } // end of k != j
      } // end of k
      
      // sum_{k!=j} Bijk*log(pijk/(1-pijk))
      d1=db.rowwise().sum();
      // sum_{k!=j} Bijk*betajk*Xi*Mik
      c1=cb.rowwise().sum();
      
      // sum over Bijk in {0,1}
      llk=llk.array()+exp(-(M.col(j).array()-c1.array())*(M.col(j).array()-c1.array())/sigmaj*0.5)*exp(d1.array());
      
    } // end of row
    
    Likelihoodi=(prodprobj.array())*llk.array()/sqrt(2*M_PI*sigmaj);
    
    //  } // end of j
    
    // take log, loglikelihood
    Likelihoodi=(Likelihoodi.array()).log();
    
    // sum_{i:n}
    Likeli[rk]=Likelihoodi.sum();
    
    
    
  }  // end of rk, hard-thresholding
  
  
  
  return(List::create(Named("constant")=constant, Named("Likelihoodi")=Likelihoodi, Named("Likeli")=Likeli, Named("betajs")=betajs, Named("sigmajs")=sigmajs,
                      Named("prodprobj")=prodprobj, Named("llk")=llk
  ));
  
}



/********* Compute EBIC *******/
// [[Rcpp::export]]
List HardThresholdDTIjC(Eigen::VectorXd inibetaj, double inisigmaj, Eigen::VectorXd seqbetaj, Eigen::MatrixXd betaK, Eigen::MatrixXd X, Eigen::MatrixXd M, Eigen::MatrixXd EB, Eigen::MatrixXd EBB, Eigen::MatrixXd EBji, Eigen::MatrixXd EBj1,
                        Eigen::MatrixXd indexjk, Eigen::MatrixXd indexjkk, Eigen::MatrixXd S, Eigen::MatrixXd distance,  double gamma, double eta, Eigen::MatrixXd B1, Eigen::MatrixXd B2, Eigen::MatrixXd B3,
                        Eigen::MatrixXd B4, Eigen::MatrixXd B5, Eigen::MatrixXd B6, Eigen::MatrixXd B7, Eigen::MatrixXd B8, Eigen::MatrixXd B9, int cutnum, int nodej){
  
  //  int j, int k, int k2
  int N=X.rows(), p=M.cols(), q=X.cols(), pp=p*(p-1)*(p-2)/2, po=p*(p-1);
  int l,s,l1;
  int i,j,k,k2;
  int jl, jl1, jk, kj, jkl, jk2, jkk, kjk,  count, countk, rk, countl;
  //  int m=lambda.size();
  // length of cutoff values
  int cut=cutnum;
  int p2=p-3;
  int row;
  int rownum;
  int lenB, lenB1;
  int ind, ind1;
  int indk;
  
  double sigmaj=1;
  
  Eigen::MatrixXd BB;
  
  Eigen::VectorXd indexk=Eigen::VectorXd::Zero(p);
  Eigen::VectorXd probi;
  Eigen::VectorXd XBMi;
  
  Eigen::VectorXd probi1;
  Eigen::VectorXd XBMi1;
  
  Eigen::VectorXd sigmajs=Eigen::VectorXd::Zero(cut);
  
  Eigen::MatrixXd Beta;
  
  Eigen::MatrixXd betajs=Eigen::MatrixXd::Zero(inibetaj.size(), cut);
  
  Eigen::MatrixXd XBM=Eigen::MatrixXd::Zero(N, po);
  
  Eigen::VectorXd betark;
  
  Eigen::VectorXd constant;
  
  Eigen::VectorXd Likelihoodi;
  Eigen::VectorXd Likeli;
  
  //  double sumLikeli=0.0;
  
  Eigen::MatrixXd db, cb;
  
  double d1, c1;
  
  double llk;
  
  
  Eigen::MatrixXd numerator=Eigen::MatrixXd::Zero(N,p-1);
  Eigen::MatrixXd prob0=Eigen::MatrixXd::Zero(N,p-1);
  
  Eigen::MatrixXd prob=Eigen::MatrixXd::Zero(N,po);
  
  // log of product of degrees: sijk=sij*sik / common neighbors in DTI
  Eigen::MatrixXd logS=Eigen::MatrixXd::Zero(N*p, p);
  logS=log(1+S.array());
  
  // log of distance
  Eigen::MatrixXd logdist=Eigen::MatrixXd::Zero(p, p);
  logdist=log(distance.array());
  logdist.diagonal()=Eigen::VectorXd::Zero(p);
  
  Eigen::VectorXd prodprobj;
  
  j=nodej;
  
  s=0;
  for (k=0; k<p; ++k){
    if (k != j){
      // jk=indexjk(j,k)-1;
      
      for (i=0; i<N; ++i){
        
        numerator(i,s)=exp(gamma*logS(i*p+j,k)-eta*logdist(j,k));
        //  numerator(i*p+k,j)=numerator(i*p+j,k);
        
        prob(i,s)=numerator(i,s)/(1+numerator(i,s));
        //  prob0(i*p+k,j)=prob0(i*p+j,k);
        
      } // end of subject i
      
      indexk[k]=s;
      
      s=s+1;
      
    } // end of k!=j
    
  } // end of k
  
  
  Eigen::VectorXd MEB;
  Eigen::VectorXd MEBXv;
  
  Eigen::MatrixXd MEBX;
  Eigen::MatrixXd MEBXX;
  
  Eigen::VectorXd MEBBXXdiag;
  
  Eigen::MatrixXd MMEB;
  // right hand: MMEBX
  Eigen::VectorXd RightHand;
  // right hand matrix
  Eigen::MatrixXd RightHandM;
  
  Eigen::VectorXd betaj;
  
  Eigen::MatrixXd betajX;
  Eigen::MatrixXd EMMj;
  
  double EMMBetaxj;
  double EBBMMBetaxj;
  
  Eigen::VectorXd betajXk;
  Eigen::VectorXd betajXl;
  Eigen::VectorXd EMMjk;
  
  Eigen::VectorXd mmebj;
  
  
  //  Eigen::VectorXd prodprobj
  
  
  //  ones=Eigen::VectorXd::Ones(inibetaj.size());
  // compute likelihood function for observed data: P(Mij, Bijk, k!=j | Mik, k!=j)
  //   j=nodej;
  
  prodprobj=Eigen::VectorXd::Ones(N);
  
  for (s=0; s<(p-1); ++s){
    prodprobj=(prodprobj.array())*(1-prob.col(s).array());
  } // end of s
  
  
  
  Likeli=Eigen::VectorXd::Zero(cut);
  
  // rk=1;
  for (rk=0; rk<cut; ++rk){
    
    Beta=Eigen::MatrixXd::Zero(q,p*(p-1));
    
    Likelihoodi=Eigen::VectorXd::Zero(N);
    
    
    // hard-thresholded betaj
    betark=betaK.col(rk);
    
    betajs.col(rk)=betark;
    
    
    // update the corresponding sigma_j
    betajX=Eigen::MatrixXd::Zero(N,p-1);
    
    betajXk=Eigen::VectorXd::Zero(p-1);
    
    EMMj=Eigen::MatrixXd::Zero(N,p-1);
    EMMjk=Eigen::VectorXd::Zero(p-1);
    
    betajXl=Eigen::VectorXd::Zero(p-1);
    
    EMMBetaxj=0.0;
    EBBMMBetaxj=0.0;
    
    countk=0;
    for (k=0; k<p; ++k){
      if (k != j){
        jk=indexjk(j,k)-1;
        
        // vector.segment(i,n): Block containing n elements, starting at position i
        Beta.col((p-1)*j+countk)=betark.segment(countk*q,q);
        
        // betajX.col(countk)=X*betaj(countk*q,q);
        betajXk=X*betark.segment(countk*q,q);
        EMMjk=EB.col(jk).array()*M.col(j).array()*M.col(k).array();
        
        // first part
        EMMBetaxj=EMMBetaxj+EMMjk.dot(betajXk);
        
        mmebj=Eigen::VectorXd::Zero(N);
        
        countl=0;
        for (l=0; l<p; ++l){
          if (l!=j){
            jl=indexjk(j,l)-1;
            
            if (l==k){
              betajXl=betajXk;
              betajXl=betajXl.array()*betajXk.array();
              // second part
              mmebj=M.col(k).array()*M.col(k).array()*EB.col(jk).array();
              EBBMMBetaxj=EBBMMBetaxj+betajXl.dot(mmebj);
            }
            
            if (l>k){
              betajXl=X*betark.segment(countl*q,q);
              jkl=indexjkk(j*p+k,l)-1;
              betajXl=betajXl.array()*betajXk.array();
              
              mmebj=M.col(k).array()*M.col(l).array()*EBB.col(jkl).array();
              EBBMMBetaxj=EBBMMBetaxj+(betajXl).dot(mmebj);
            }
            
            if (l<k){
              betajXl=X*betark.segment(countl*q,q);
              jkl=indexjkk(j*p+l,k)-1;
              betajXl=betajXl.array()*betajXk.array();
              
              mmebj=M.col(k).array()*M.col(l).array()*EBB.col(jkl).array();
              EBBMMBetaxj=EBBMMBetaxj+(betajXl).dot(mmebj);
            }
            ++countl;
          } // end of l!=j
          
        } // end of l
        ++countk;
      } // end of k!=j
    } // end of k
    
    //  Sigma[j]=M.col(j).squaredNorm()/N;
    sigmaj=M.col(j).squaredNorm()/N-(2*EMMBetaxj-EBBMMBetaxj)/N;
    
    sigmajs[rk]=sigmaj;
    
    
    // compute likelihood
    
    for (k=0; k<p; ++k){
      if (k!=j){
        
        jk=indexjk(j,k)-1;
        
        XBM.col(jk)=(X*Beta.col(jk)).array()*(M.col(k).array());
        
      } // end of k!=j
    } // end of k
    
    
    for (i=0; i<N; ++i){
      
      // 0.1<=E(Bijk)<=0.9
      lenB=EBji.row(i).sum();
      
      // E(Bijk)>0.9
      lenB1=EBj1.row(i).sum();
      
      
      if (lenB==1){
        BB=B1;
      }
      if (lenB==2){
        BB=B2;
      }
      if (lenB==3){
        BB=B3;
      }
      if (lenB==4){
        BB=B4;
      }
      if (lenB==5){
        BB=B5;
      }
      if (lenB==6){
        BB=B6;
      }
      if (lenB==7){
        BB=B7;
      }
      if (lenB==8){
        BB=B8;
      }
      if (lenB==9){
        BB=B9;
      }
      
      // for each k, for summation
      probi=Eigen::VectorXd::Zero(lenB);
      XBMi=Eigen::VectorXd::Zero(lenB);
      
      // set bijk=1
      probi1=Eigen::VectorXd::Zero(lenB1);
      XBMi1=Eigen::VectorXd::Zero(lenB1);
      
      // if E(Bijk|...)>0.9, let bijk=1; if E(Bijk|...)<0.1, let bijk=0; if 0.1<=E(Bijk|...)<=0.9, calculate in the summation
      ind=0;
      ind1=0;
      for (k=0; k<p; ++k){
        if (k!=j){
          indk=indexk[k];
          jk=indexjk(j,k)-1;
          
          // extract index needed in the calculation
          if (EBji(i,indk)==1) {
            
            probi[ind]=prob(i,indk);
            XBMi[ind]=XBM(i,jk);
            
            ++ind;
            
          } // end of EBji
          
          // E(Bijk)>0.9
          if (EBj1(i,indk)==1) {
            probi1[ind1]=prob(i,indk);
            XBMi1[ind1]=XBM(i,jk);
            
            ++ind1;
            
          } // end of Bj
          
          
        } // end of k!=j
      } // end of k
      
      db=Eigen::VectorXd::Zero(lenB);
      cb=Eigen::VectorXd::Zero(lenB);
      
      d1=0.0;
      c1=0.0;
      
      llk=0.0;
      
      if (lenB==0){
        
        if (lenB1>0){
          c1=XBMi1.sum();
          d1=(log(probi1.array()/(1-probi1.array()))).sum();
        } else{
          
          c1=0;
          d1=0;
        }
        
        llk=exp(-(M(i,j)-c1)*(M(i,j)-c1)/sigmaj*0.5)*exp(d1);
        
      } else{
        
        rownum=BB.rows();
        
        for (row=0; row<rownum; ++row){
          
          db=BB.row(row).array();
          
          db=db.array()*log(probi.array()/(1-probi.array()));
          
          // Bijk*betajk*Xi*Mik
          cb=BB.row(row).array();
          cb=cb.array()*XBMi.array();
          
          if (lenB1>0){
            d1=db.sum()+(log(probi1.array()/(1-probi1.array()))).sum();
            c1=cb.sum()+XBMi1.sum();
          } else {
            // sum_{k!=j} Bijk*log(pijk/(1-pijk))
            d1=db.sum();
            // sum_{k!=j} Bijk*betajk*Xi*Mik
            c1=cb.sum(); 
            
          }
          
          // sum over Bijk in {0,1}
          llk=llk+exp(-(M(i,j)-c1)*(M(i,j)-c1)/sigmaj*0.5)*exp(d1);
        }  // end of row
        
      } // end of else
      
      Likelihoodi[i]=(prodprobj[i])*llk/sqrt(2*M_PI*sigmaj);
      
    }  // end of subject i
    
    
    // take log, loglikelihood
    Likelihoodi=(Likelihoodi.array()).log();
    
    // sum_{i:n}
    Likeli[rk]=Likelihoodi.sum();
    
    
    
  }  // end of rk, hard-thresholding
  
  
  
  return(List::create(Named("constant")=constant, Named("Likelihoodi")=Likelihoodi, Named("Likeli")=Likeli, Named("betajs")=betajs, Named("sigmajs")=sigmajs,
                      Named("prodprobj")=prodprobj, Named("llk")=llk
  ));
  
}
