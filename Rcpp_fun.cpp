#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]

double Sfun1c(double z, double t){
  double a = 1 - t/std::abs(z);
  if(a> 0){
    return z*a;
  }else {
    return 0;
  }
}
// [[Rcpp::export]]
arma::vec Sfunc(arma::vec z, double t){
  double a = 1- t/norm(z, 2);
  arma::vec b = zeros(size(z));
  if ( a> 0){
    return z*a;
  }else{
    return b;
  }
}


// [[Rcpp::export]]
arma::vec GMCPbcc(int n, int q, arma::mat  Wjt, arma::vec bj, arma::vec r, double kappa1, double kappa2, double xi){
  double eps = 0.05;
  arma::vec muj = (Wjt*r)/n + bj;
  double normbj = norm(bj,2);
  double ghat;
  arma::vec bjnew;
  if (kappa1 ==0){
    ghat = 1 + eps;
  }else{
    if (normbj == 0){
      ghat = 2^31;
    }else if(normbj > xi*sqrt(q+1)*kappa1){
      ghat = 1+eps;
    }else{
      ghat = 1+eps + (1/normbj)*sqrt(q+1)*kappa1-normbj/xi;
    }
  }
  arma::vec vj = muj;
  for (int k=1; k< q+1; k ++){
    if (std::abs(muj(k))<=xi*kappa2*ghat){
      vj(k) = Sfun1c(muj(k),kappa2)/(1-1/(xi*ghat));
    }
  }
  double normvj = norm(vj,2);
  if (normvj > (xi*sqrt(q+1)*kappa1)) {
    bjnew = vj;
  } else {
    bjnew = (xi/((1+eps)*xi-1))*Sfunc(vj,sqrt(q+1)*kappa1);
  }
  
  return bjnew;
}



// [[Rcpp::export]]

List gaussianetscpp(arma::vec yy,  arma::mat XX, arma::mat Wjj, 
                          double xi, int max_iter, arma::mat Matla, 
                          double eps, arma::uvec iset){
  arma::mat L_mat, L_tran_inverse, Wjortht;
  arma::mat Wjorth;
  int n = yy.size();
  int L = Matla.n_rows;
  int q = XX.n_cols;
  double lambda1, lambda2, diff;
  arma::mat beta(2*q+1,L);
  beta.zeros();
  arma::vec iter(L);
  iter.zeros();
  arma::vec beta0(2*q+1); 
  arma::vec anew, bnew, ya(n), r(n), beta1;
  arma::mat Wjjt = Wjj.t();
  L_mat = chol(Wjjt*Wjj/n, "upper");
  L_tran_inverse = inv(L_mat);
  Wjorth = Wjj*L_tran_inverse;
  Wjortht = Wjorth.t();
  arma::mat XX_t = XX.t();
  arma::mat Mat = inv(XX_t*XX)*XX_t;

  for(int l=0; l < L; l++){
    lambda1 = Matla(l,0);
    lambda2 = Matla(l,1);
    if(l==0){
      beta0 = beta.col(0);
    }else{
      beta0 = beta.col(l-1);
      if(max(beta0) > 3){
        beta0.zeros();
      }
    }
    beta0.zeros();
    while(iter(l) < max_iter){
      bnew = beta0(iset);
      iter(l) = iter(l) + 1;
      //update a
      ya = yy - Wjj*bnew;
      anew = Mat*ya;
      //update b
      r = ya - XX*anew;
      bnew = GMCPbcc(n,q, Wjortht, bnew,r,lambda1,lambda2,xi);
      bnew = L_tran_inverse*bnew;
      // beta.col(l) <- arma::join_cols(anew, bnew);
      for(int i = 0; i<q; i++){
        beta(i,l) = anew(i);
      }
      for(int i = 0; i<q+1; i++){
        beta(i+q,l) = bnew(i);
      }
      diff = norm(beta.col(l) - beta0 ,2);
      if (diff < eps or max(beta.col(l)) > 3){
        break;
      }
      beta0 = beta.col(l);
    }
  }
  List res = List::create(Named("beta")=  beta, _["iter"] = iter);
  return(res);
}


// [[Rcpp::export]]

List binomialetscpp(arma::vec y,  arma::mat XX, arma::mat Wjj, 
                     double xi, int max_iter, arma::mat Matla, 
                     double eps, arma::uvec iset){
  arma::mat L_mat, L_tran_inverse, Wjortht;
  arma::mat Wjorth;
  int n = y.size();
  int L = Matla.n_rows;
  int q = XX.n_cols;
  double lambda1, lambda2, diff;
  arma::mat beta(2*q+2,L), Mat;
  beta.zeros();
  arma::vec iter(L);
  iter.zeros();
  arma::vec beta0(2*q+2); 
  arma::vec anew, bnew, ya(n), r(n), beta1;
  arma::mat I1 = arma::ones<arma::mat>(n,1);
  XX.insert_cols(0,I1);
  arma::mat XW = arma::join_rows(XX, Wjj);
  arma::vec eta(n), mu(n), w(n), yy(n);
  arma::mat Weight(n,n), XWn, Wjjt;
  Weight.zeros();
  
  for(int l=0; l < L; l++){
    lambda1 = Matla(l,0);
    lambda2 = Matla(l,1);
    if(l==0){
      beta0 = beta.col(0);
    }else{
      beta0 = beta.col(l-1);
      if(max(beta0) > 3){
        beta0.zeros();
      }
    }
    while(iter(l) < max_iter){
      //change

      eta = XW * beta0;
      mu = 1/(1+arma::exp(-eta));

      for (int i=0; i<n; i++){
        // mu(i) = 1/(1+ exp(-eta(i)));
        w(i) = std::max(mu(i)*(1-mu(i)),0.0001);
        yy(i)= sqrt(w(i))* (eta(i) + (y(i) - mu(i))/ w(i));
        Weight(i,i) = sqrt(w(i));
      }

      XWn = Weight * XW;
      XX = XW.cols(0,q);
      Wjj = XW.cols(q+1,2*q+1);
      
      Wjjt = Wjj.t();
      L_mat = chol(Wjjt*Wjj/n, "upper");
      L_tran_inverse = inv(L_mat);
      Wjorth = Wjj*L_tran_inverse;
      Wjortht = Wjorth.t();
      arma::mat XX_t = XX.t();
      Mat = inv(XX_t*XX)*XX_t;
      
      bnew = beta0(iset);
      iter(l) = iter(l) + 1;
      //update a
      ya = yy - Wjj*bnew;
      anew = Mat*ya;
      //update b
      r = ya - XX*anew;
      bnew = GMCPbcc(n,q, Wjortht, bnew,r,lambda1,lambda2,xi);
      bnew = L_tran_inverse*bnew;
      // beta.col(l) <- arma::join_cols(anew, bnew);
      for(int i = 0; i<q+1; i++){
        beta(i,l) = anew(i);
      }
      for(int i = 0; i<q+1; i++){
        beta(i+q+1,l) = bnew(i);
      }
      diff = norm(beta.col(l) - beta0 ,2);
      if (diff < eps or max(beta.col(l)) > 3){
        break;
      }
      beta0 = beta.col(l);
   }
  }
  List res = List::create(Named("beta")=  beta.rows(1,2*q+1), 
                          _["beta0"]= beta.row(0), _["iter"] = iter);
  return(res);
}



