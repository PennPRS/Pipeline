#include <arma-strict-R-headers.h>
#include <utils.h>
#include <SFBM.h>
// --------------------------------//--------------------------------
  //#include <omp.h>
  // //[[Rcpp::plugins(openmp)]]
// --------------------------------//--------------------------------
  
  // Important: without // [[Rcpp::depends(RcppTN)]] gives error

// [[Rcpp::depends(RcppTN)]]
using namespace arma;
using namespace Rcpp;
using namespace std;
using namespace stats;


const double log2pi = std::log(2.0 * M_PI);

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
vec Mahalanobis(mat x, rowvec center, mat cov) {
  int n = x.n_rows;
  mat x_cen;
  x_cen.copy_size(x);
  for (int i=0; i < n; i++) {
    x_cen.row(i) = x.row(i) - center;
  }
  return sum((x_cen * cov.i()) % x_cen, 1);
}

// [[Rcpp::export]]
mat rmvncpp(int n, vec mu, mat V){
  int p = mu.n_elem;
  //if(any(is.na(match(dim(V),p))))
    //stop("Dimension problem!")
  mat D = chol(V); //, tol = 1e-50
  mat T1(n,p);
  mat T2(n,p);
  for(int co = 0; co < p; co ++)
  {
    T1.col(co) = (vec)rnorm(n,0,1);
    vec tem(n);
    tem.fill(mu(co));
    T2.col(co) = tem;
  }
  mat T = (T1 * D + T2).t();
  return(T);
}

// [[Rcpp::export]]
vec dmvnormcpp(mat x, rowvec mean, mat sigma) {
  vec distval = Mahalanobis(x,  mean, sigma);
  double logdet = sum(log(eig_sym(sigma)));
  vec logretval = -( (x.n_cols * log2pi + logdet + distval)/2  ) ;
  return(exp(logretval));
}

// sample bh and update its posterior mean for overlapping SNPs; 
// return: 
  
  //update_beta_ind_cpp(bh1[x], bjh1[s-1,TaggingSNPinx1[[x]]], bjh1[s-1,x], Di1[[x]], sigmasq1, hsqh1[s-1], p10)
//update.bind = function(bh, bjtag, bj0, dtag, sigmasq, hsq, p0)
  
  // [[Rcpp::export]]
mat update_beta1(vec snp_index, Environment cormat, const NumericVector& Bh, vec Bjh0, const NumericVector& Sigmasq,
                 double hsq, double p0, int sparse){
  XPtr<SFBM> Cmat = cormat["address"];
  //uvec snp_indext = snp_indext0 - 1;
  double len_snp_index = snp_index.n_elem;
  //arma::vec Bjh0(Bjh00.begin(), M);
  mat outputs(2,len_snp_index);
  for (int i = 0; i < len_snp_index; i++){
    int x = snp_index(i) - 1;
    //int x0 = indmat(x) - 1; // snp index
    double bh = Bh(x);
    double sigmasq = Sigmasq(x);
    //uvec tagstem = TaggingSNPinx[x];
    //uvec tags = tagstem - 1;
    //vec bjtag = Bjh0.elem(tags);
    
    double bj0 = Bjh0(x); //double bj0 = Bjh0(x);
    //vec dtag = Di[x];
    
    //double mu = dot(dtag, bjtag) - bj0;
    //vec Bjh0t =  Bjh0.elem(snp_indext);
    //double mu = Cmat->dot_col(x0, Bjh0t) - bj0;
    double mu = Cmat->dot_col(x, Bjh0) - bj0;
    double v = hsq + sigmasq;
    double fnormal = p0 * R::dnorm(bh, mu, sqrt(v), 0);
    double f0 = (1-p0) * R::dnorm(bh, mu, sqrt(sigmasq), 0);
    double p_n = fnormal;
    double p_d = fnormal + f0;
    double ptilde = 0;
    double postvar = 0;
    double bjnew = 0;
    double postmean = 0;
    double pm = 0;
    //double dens = 0;
    if (p_n > 0){
      ptilde = p_n/p_d;
      postvar = 1/(1/sigmasq + 1/hsq);
      postmean = ((bh - mu)/sigmasq) * postvar;
      //NumericVector g = runif(1);
      //double gamma = g(0);
      double gamma = ::unif_rand();
      if (gamma < ptilde){
        // mat tem = rnorm(1, postmean, sqrt(postvar));
        // bjnew = tem(0);
        bjnew = postmean + ::norm_rand() * ::sqrt(postvar);
      }
      pm = ptilde * postmean;
      // sparse
      if (sparse == 1){
        if (ptilde < p0){
          bjnew = 0;
          pm = 0;
        }
      }
    }
    //outputs(0,i) = ptilde; outputs(1,i) = postmean;
    outputs(0,i) = bjnew;
    outputs(1,i) = pm;
    Bjh0(x) = bjnew;
  }
  return(outputs);
}


// [[Rcpp::export]]
List update_beta2(vec snp_index, mat indmat, Environment cormat1, Environment cormat2,
                  vec Bh1, vec Bh2, vec Bjh1, vec Bjh2,
                  const NumericVector& Sigmasq1, const NumericVector& Sigmasq2,
                  double hsq1, double hsq2, const NumericVector& rhos,
                  double p_11, double p_10, double p_01, vec sparse){
  XPtr<SFBM> Cmat1 = cormat1["address"];
  XPtr<SFBM> Cmat2 = cormat2["address"];
  
  int m1 = Bh1.n_elem; int m2 = Bh2.n_elem;
  myassert_size(Cmat1->nrow(), m1); myassert_size(Cmat1->ncol(), m1);
  myassert_size(Cmat2->nrow(), m2); myassert_size(Cmat2->ncol(), m2);
  
  double len_snp_index = snp_index.n_elem;
  mat ptilde(4,len_snp_index);
  ptilde.fill(0);
  mat bjnew(2,len_snp_index);
  bjnew.fill(0);
  mat pm(2,len_snp_index);
  pm.fill(0);
  for (int i = 0; i < len_snp_index; i++){
    int x = snp_index(i) - 1; // snp index
    int x1 = indmat(x,0) - 1; // original index1
    int x2 = indmat(x,1) - 1; // original index2
    
    double bh1 = Bh1(x1); double bh2 = Bh2(x2);
    double sigmasq1 = Sigmasq1(x1); double sigmasq2 = Sigmasq2(x2);
    double rho = rhos(x);
    double bj1 = Bjh1(x1);
    double bj2 = Bjh2(x2);
    
    double mu1 = Cmat1->dot_col(x1, Bjh1) - bj1;
    double mu2 = Cmat2->dot_col(x2, Bjh2) - bj2;
    
    rowvec mu(2);
    mu.col(0) = mu1; mu.col(1) = mu2;
    double v11 = hsq1 + sigmasq1;
    double v22 = hsq2 + sigmasq2;
    double v12 = sqrt(hsq1*hsq2) * rho;
    double s1 = hsq1 + sigmasq1;
    double s2 = hsq2 + sigmasq2;
    mat v(2,2);
    v(0,0) = v11;
    v(1,0) = v(0,1) = v12;
    v(1,1) = v22;
    mat bh(1,2);
    bh(0,0) = bh1;
    bh(0,1) = bh2;
    vec f11tem = p_11 * dmvnormcpp(bh, mu, v);
    double f11 = f11tem(0);
    double f01 = p_01 * R::dnorm(bh1,mu1,sqrt(sigmasq1),0) * R::dnorm(bh2,mu2,sqrt(s2),0);
    double f10 = p_10 * R::dnorm(bh1,mu1,sqrt(s1),0) * R::dnorm(bh2,mu2,sqrt(sigmasq2),0);
    double f00 = (1-p_11-p_01-p_10) * R::dnorm(bh1,mu1,sqrt(sigmasq1),0) * R::dnorm(bh2,mu2,sqrt(sigmasq2),0);
    double p_d = f11 + f01 + f10 + f00;
    
    if (p_d>0){
      vec fs(4);
      fs(0) = f11/p_d; fs(1) = f01/p_d; fs(2) = f10/p_d; fs(3) = f00/p_d;
      ptilde.col(i) = fs;
      
      mat v0(2,2);
      v0(0,0) = hsq1; v0(0,1) = v0(1,0) = sqrt(hsq1*hsq2)*rho; v0(1,1) = hsq2;
      mat v01 = inv(v0);
      v01(0,0) = v01(0,0) + 1/sigmasq1; v01(1,1) = v01(1,1) + 1/sigmasq2;
      mat postvar = inv(v01);
      mat meantem(2,1);
      meantem = bh.t() - mu.t();
      meantem(0,0) = meantem(0,0)/sigmasq1;
      meantem(1,0) = meantem(1,0)/sigmasq2;
      mat postmean = postvar * meantem;
      double postvar1 = 1/(1/sigmasq1 + 1/hsq1);
      double postmean1 = ((bh1 - mu1)/sigmasq1) * postvar1;
      double postvar2 = 1/(1/sigmasq2 + 1/hsq2);
      double postmean2 = ((bh2 - mu2)/sigmasq2) * postvar2;
      
      NumericVector g = runif(1);
      double gamma = g(0);
      //
        if (gamma < fs(0)){
          bjnew.col(i) = rmvncpp(1, postmean, postvar);
        }
      if ((gamma >= fs(0))&(gamma < (fs(0)+fs(1)))){
        NumericVector tem = rnorm(1, postmean2, sqrt(postvar2));
        bjnew(1,i) = tem(0);
      }
      if ((gamma >= (fs(0)+fs(1)))&(gamma < (fs(0)+fs(1)+fs(2)))){
        NumericVector tem = rnorm(1, postmean1, sqrt(postvar1));
        bjnew(0,i) = tem(0);
      }
      double pm1tem = fs(0)*postmean(0,0) + fs(2)*postmean1;
      double pm2tem = fs(0)*postmean(1,0) + fs(1)*postmean2;
      pm(0,i) = pm1tem;
      pm(1,i) = pm2tem;
      
      // sparse
      if (sparse(0) == 1){
        double fs1 = fs(0) + fs(2);
        double fs1_0 = p_11 + p_10;
        if (fs1 < fs1_0){
          bjnew(0,i) = 0;
          pm(0,i) = 0;
          postmean(0,0) = 0;
          postmean1 = 0;
        }
      }
      if (sparse(1) == 1){
        double fs2 = fs(0) + fs(1);
        double fs2_0 = p_11 + p_01;
        if (fs2 < fs2_0){
          bjnew(1,i) = 0;
          pm(1,i) = 0;
          postmean(1,0) = 0;
          postmean2 = 0;
        }
      }
    }
    Bjh1(x1) = bjnew(0,i);
    Bjh2(x2) = bjnew(1,i);
  }
  return List::create(_["bjnew"] = bjnew, _["pm"] = pm); //_["ptilde"] = ptilde,
}



// [[Rcpp::export]]
List update_beta3(vec snp_index, mat indmat, Environment cormat1, Environment cormat2, Environment cormat3,
                  vec Bh1, vec Bh2, vec Bh3, vec Bjh1, vec Bjh2, vec Bjh3,
                  const NumericVector& Sigmasq1, const NumericVector& Sigmasq2, const NumericVector& Sigmasq3,
                  double hsq1, double hsq2, double hsq3,
                  const NumericVector& rhos12, const NumericVector& rhos13, const NumericVector& rhos23,
                  double p111, double p110, double p101, double p011,
                  double p100, double p010, double p001, double p000, vec sparse){
  XPtr<SFBM> Cmat1 = cormat1["address"];
  XPtr<SFBM> Cmat2 = cormat2["address"];
  XPtr<SFBM> Cmat3 = cormat3["address"];
  
  int m1 = Bh1.n_elem; int m2 = Bh2.n_elem; int m3 = Bh3.n_elem;
  myassert_size(Cmat1->nrow(), m1); myassert_size(Cmat1->ncol(), m1);
  myassert_size(Cmat2->nrow(), m2); myassert_size(Cmat2->ncol(), m2);
  myassert_size(Cmat3->nrow(), m3); myassert_size(Cmat3->ncol(), m3);
  
  double len_snp_index = snp_index.n_elem;
  mat ptilde(8,len_snp_index);
  ptilde.fill(0);
  mat bjnew(3,len_snp_index);
  bjnew.fill(0);
  mat pm(3,len_snp_index);
  pm.fill(0);
  for (int i = 0; i < len_snp_index; i++){
    int x = snp_index(i) - 1; // snp index
    int x1 = indmat(x,0) - 1; // original index1
    int x2 = indmat(x,1) - 1; // original index2
    int x3 = indmat(x,2) - 1; // original index3
    
    double bh1 = Bh1(x1); double bh2 = Bh2(x2); double bh3 = Bh3(x3);
    double sigmasq1 = Sigmasq1(x1); double sigmasq2 = Sigmasq2(x2); double sigmasq3 = Sigmasq3(x3);
    double rho12 = rhos12(x); double rho13 = rhos13(x); double rho23 = rhos23(x);
    double bj1 = Bjh1(x1); double bj2 = Bjh2(x2); double bj3 = Bjh3(x3);
    
    double mu1 = Cmat1->dot_col(x1, Bjh1) - bj1;
    double mu2 = Cmat2->dot_col(x2, Bjh2) - bj2;
    double mu3 = Cmat3->dot_col(x3, Bjh3) - bj3;
    
    rowvec mu(3);
    mu.col(0) = mu1; mu.col(1) = mu2; mu.col(2) = mu3;
    double v11 = hsq1 + sigmasq1; double v22 = hsq2 + sigmasq2; double v33 = hsq3 + sigmasq3;
    double v12 = sqrt(hsq1*hsq2) * rho12; double v13 = sqrt(hsq1*hsq3) * rho13; double v23 = sqrt(hsq2*hsq3) * rho23;
    double s1 = hsq1 + sigmasq1; double s2 = hsq2 + sigmasq2; double s3 = hsq3 + sigmasq3;
    mat v(3,3);
    v(0,0) = v11; v(1,1) = v22; v(2,2) = v33;
    v(1,0) = v(0,1) = v12; v(2,0) = v(0,2) = v13; v(2,1) = v(1,2) = v23;
    
    mat bh(1,3);
    bh(0,0) = bh1; bh(0,1) = bh2; bh(0,2) = bh3;
    
    vec f111tem = p111 * dmvnormcpp(bh, mu, v);
    double f111 = f111tem(0);
    uvec f110ind(2); f110ind(0) = 0; f110ind(1) = 1;
    vec f110tem = dmvnormcpp(bh.cols(f110ind), mu.cols(f110ind), v.submat(f110ind,f110ind));
    double f110 = p110 * f110tem(0) * R::dnorm(bh3,mu3,sqrt(sigmasq3),0);
    uvec f101ind(2); f101ind(0) = 0; f101ind(1) = 2;
    vec f101tem = dmvnormcpp(bh.cols(f101ind), mu.cols(f101ind), v.submat(f101ind,f101ind));
    double f101 = p101 * f101tem(0) * R::dnorm(bh2,mu2,sqrt(sigmasq2),0);
    uvec f011ind(2); f011ind(0) = 1; f011ind(1) = 2;
    vec f011tem = dmvnormcpp(bh.cols(f011ind), mu.cols(f011ind), v.submat(f011ind,f011ind));
    double f011 = p011 * f011tem(0) * R::dnorm(bh1,mu1,sqrt(sigmasq1),0);
    double f100 = p100 * R::dnorm(bh1,mu1,sqrt(s1),0) * R::dnorm(bh2,mu2,sqrt(sigmasq2),0) * R::dnorm(bh3,mu3,sqrt(sigmasq3),0);
    double f010 = p010 * R::dnorm(bh2,mu2,sqrt(s2),0) * R::dnorm(bh1,mu1,sqrt(sigmasq1),0) * R::dnorm(bh3,mu3,sqrt(sigmasq3),0);
    double f001 = p001 * R::dnorm(bh3,mu3,sqrt(s3),0) * R::dnorm(bh1,mu1,sqrt(sigmasq1),0) * R::dnorm(bh2,mu2,sqrt(sigmasq2),0);
    double f000 = p000 * R::dnorm(bh1,mu1,sqrt(sigmasq1),0) * R::dnorm(bh2,mu2,sqrt(sigmasq2),0) * R::dnorm(bh3,mu3,sqrt(sigmasq3),0);
    double p_d = f111 + f110 + f101 + f011 + f100 + f010 + f001 + f000;
    
    
    if (p_d>0){
      vec fs(8);
      fs(0) = f111/p_d; fs(1) = f110/p_d; fs(2) = f101/p_d; fs(3) = f011/p_d;
      fs(4) = f100/p_d; fs(5) = f010/p_d; fs(6) = f001/p_d; fs(7) = f000/p_d;
      ptilde.col(i) = fs;
      
      mat v0(3,3);
      v0(0,0) = hsq1; v0(1,1) = hsq2; v0(2,2) = hsq3;
      v0(0,1) = v0(1,0) = sqrt(hsq1*hsq2)*rho12;
      v0(0,2) = v0(2,0) = sqrt(hsq1*hsq3)*rho13;
      v0(1,2) = v0(2,1) = sqrt(hsq2*hsq3)*rho23;
      mat v0inv = inv(v0);
      v0inv(0,0) = v0inv(0,0) + 1/sigmasq1; v0inv(1,1) = v0inv(1,1) + 1/sigmasq2; v0inv(2,2) = v0inv(2,2) + 1/sigmasq3;
      mat postvar = inv(v0inv);
      mat meantem(3,1);
      meantem = bh.t() - mu.t();
      meantem(0,0) = meantem(0,0)/sigmasq1;
      meantem(1,0) = meantem(1,0)/sigmasq2;
      meantem(2,0) = meantem(2,0)/sigmasq3;
      mat postmean = postvar * meantem;
      
      mat v0_12 = v0.submat(f110ind,f110ind); mat v0_12inv = inv(v0_12);
      v0_12inv(0,0) = v0_12inv(0,0) + 1/sigmasq1; v0_12inv(1,1) = v0_12inv(1,1) + 1/sigmasq2;
      mat postvar12 = inv(v0_12inv); mat meantem12(2,1);
      meantem12 = meantem.rows(f110ind); mat postmean12 = postvar12 * meantem12;
      
      mat v0_13 = v0.submat(f101ind,f101ind); mat v0_13inv = inv(v0_13);
      v0_13inv(0,0) = v0_13inv(0,0) + 1/sigmasq1; v0_13inv(1,1) = v0_13inv(1,1) + 1/sigmasq3;
      mat postvar13 = inv(v0_13inv); mat meantem13(2,1);
      meantem13 = meantem.rows(f101ind); mat postmean13 = postvar13 * meantem13;
      
      mat v0_23 = v0.submat(f011ind,f011ind); mat v0_23inv = inv(v0_23);
      v0_23inv(0,0) = v0_23inv(0,0) + 1/sigmasq2; v0_23inv(1,1) = v0_23inv(1,1) + 1/sigmasq3;
      mat postvar23 = inv(v0_23inv); mat meantem23(2,1);
      meantem23 = meantem.rows(f011ind); mat postmean23 = postvar23 * meantem23;
      
      
      double postvar1 = 1/(1/sigmasq1 + 1/hsq1); double postmean1 = ((bh1 - mu1)/sigmasq1) * postvar1;
      double postvar2 = 1/(1/sigmasq2 + 1/hsq2); double postmean2 = ((bh2 - mu2)/sigmasq2) * postvar2;
      double postvar3 = 1/(1/sigmasq3 + 1/hsq3); double postmean3 = ((bh3 - mu3)/sigmasq3) * postvar3;
      
      NumericVector g = runif(1); double gamma = g(0);
      if (gamma < fs(0)){
        bjnew.col(i) = rmvncpp(1, postmean, postvar);
      }
      
      uvec indi(1); indi(0) = i;
      // 1 1 0
      if ((gamma >= fs(0))&(gamma < (fs(0)+fs(1)))){
        bjnew.submat(f110ind,indi) = rmvncpp(1, postmean12, postvar12);
      }
      // 1 0 1
      if ((gamma >= (fs(0)+fs(1)))&(gamma < (fs(0)+fs(1)+fs(2)))){
        bjnew.submat(f101ind,indi) = rmvncpp(1, postmean13, postvar13);
      }
      // 0 1 1
      if ((gamma >= (fs(0)+fs(1)+fs(2)))&(gamma < (fs(0)+fs(1)+fs(2)+fs(3)))){
        bjnew.submat(f011ind,indi) = rmvncpp(1, postmean23, postvar23);
      }
      // 1 0 0
      if ((gamma >= (fs(0)+fs(1)+fs(2)+fs(3)))&(gamma < (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)))){
        NumericVector tem = rnorm(1, postmean1, sqrt(postvar1));
        bjnew(0,i) = tem(0);
      }
      // 0 1 0
      if ((gamma >= (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)))&(gamma < (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)))){
        NumericVector tem = rnorm(1, postmean2, sqrt(postvar2));
        bjnew(1,i) = tem(0);
      }
      // 0 0 1
      if ((gamma >= (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)))&(gamma < (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)))){
        NumericVector tem = rnorm(1, postmean3, sqrt(postvar3));
        bjnew(2,i) = tem(0);
      }
      
      double pm1tem = fs(0)*postmean(0,0) + fs(1)*postmean12(0,0) + fs(2)*postmean13(0,0) + fs(4)*postmean1;
      double pm2tem = fs(0)*postmean(1,0) + fs(1)*postmean12(1,0) + fs(3)*postmean23(0,0) + fs(5)*postmean2;
      double pm3tem = fs(0)*postmean(2,0) + fs(2)*postmean13(1,0) + fs(3)*postmean23(1,0) + fs(6)*postmean3;
      pm(0,i) = pm1tem;
      pm(1,i) = pm2tem;
      pm(2,i) = pm3tem;
      
      // sparse
      if (sparse(0) == 1){
        double fs1 = fs(0) + fs(1) + fs(2) + fs(4);
        double fs1_0 = p111 + p110 + p101 + p100;
        if (fs1 < fs1_0){
          bjnew(0,i) = 0;
          pm(0,i) = 0;
          postmean(0,0) = 0;
          postmean1 = 0;
        }
      }
      if (sparse(1) == 1){
        double fs2 = fs(0) + fs(1) + fs(3) + fs(5);
        double fs2_0 = p111 + p110 + p011 + p010;
        if (fs2 < fs2_0){
          bjnew(1,i) = 0;
          pm(1,i) = 0;
          postmean(1,0) = 0;
          postmean2 = 0;
        }
      }
      if (sparse(2) == 1){
        double fs3 = fs(0) + fs(2) + fs(3) + fs(6);
        double fs3_0 = p111 + p101 + p011 + p001;
        if (fs3 < fs3_0){
          bjnew(2,i) = 0;
          pm(2,i) = 0;
          postmean(2,0) = 0;
          postmean3 = 0;
        }
      }
    }
    Bjh1(x1) = bjnew(0,i);
    Bjh2(x2) = bjnew(1,i);
    Bjh3(x3) = bjnew(2,i);
  }
  return List::create(_["bjnew"] = bjnew, _["pm"] = pm); //_["ptilde"] = ptilde,
}


// [[Rcpp::export]]
List update_beta4(vec snp_index, mat indmat,
                  Environment cormat1, Environment cormat2, Environment cormat3, Environment cormat4,
                  vec Bh1, vec Bh2, vec Bh3, vec Bh4, vec Bjh1, vec Bjh2, vec Bjh3, vec Bjh4,
                  const NumericVector& Sigmasq1, const NumericVector& Sigmasq2, const NumericVector& Sigmasq3, const NumericVector& Sigmasq4,
                  double hsq1, double hsq2, double hsq3, double hsq4,
                  const NumericVector& rhos12, const NumericVector& rhos13, const NumericVector& rhos14,
                  const NumericVector& rhos23, const NumericVector& rhos24, const NumericVector& rhos34,
                  double p1111, double p1110, double p1101, double p1011, double p0111,
                  double p1100, double p1010, double p1001, double p0110, double p0101, double p0011,
                  double p1000, double p0100, double p0010, double p0001, double p0000, vec sparse){
  XPtr<SFBM> Cmat1 = cormat1["address"];
  XPtr<SFBM> Cmat2 = cormat2["address"];
  XPtr<SFBM> Cmat3 = cormat3["address"];
  XPtr<SFBM> Cmat4 = cormat4["address"];
  
  int m1 = Bh1.n_elem; int m2 = Bh2.n_elem; int m3 = Bh3.n_elem; int m4 = Bh4.n_elem;
  myassert_size(Cmat1->nrow(), m1); myassert_size(Cmat1->ncol(), m1);
  myassert_size(Cmat2->nrow(), m2); myassert_size(Cmat2->ncol(), m2);
  myassert_size(Cmat3->nrow(), m3); myassert_size(Cmat3->ncol(), m3);
  myassert_size(Cmat4->nrow(), m4); myassert_size(Cmat4->ncol(), m4);
  
  double len_snp_index = snp_index.n_elem;
  mat ptilde(16,len_snp_index);
  ptilde.fill(0);
  mat bjnew(4,len_snp_index);
  bjnew.fill(0);
  mat pm(4,len_snp_index);
  pm.fill(0);
  
  uvec ind12(2); ind12(0) = 0; ind12(1) = 1;
  uvec ind13(2); ind13(0) = 0; ind13(1) = 2;
  uvec ind14(2); ind14(0) = 0; ind14(1) = 3;
  uvec ind23(2); ind23(0) = 1; ind23(1) = 2;
  uvec ind24(2); ind24(0) = 1; ind24(1) = 3;
  uvec ind34(2); ind34(0) = 2; ind34(1) = 3;
  uvec ind123(3); ind123(0) = 0; ind123(1) = 1; ind123(2) = 2;
  uvec ind124(3); ind124(0) = 0; ind124(1) = 1; ind124(2) = 3;
  uvec ind134(3); ind134(0) = 0; ind134(1) = 2; ind134(2) = 3;
  uvec ind234(3); ind234(0) = 1; ind234(1) = 2; ind234(2) = 3;
  uvec f1110ind(3); f1110ind(0) = 0; f1110ind(1) = 1; f1110ind(2) = 2;
  uvec f1101ind(3); f1101ind(0) = 0; f1101ind(1) = 1; f1101ind(2) = 3;
  uvec f1011ind(3); f1011ind(0) = 0; f1011ind(1) = 2; f1011ind(2) = 3;
  uvec f0111ind(3); f0111ind(0) = 1; f0111ind(1) = 2; f0111ind(2) = 3;
  uvec f1100ind(2); f1100ind(0) = 0; f1100ind(1) = 1;
  uvec f1010ind(2); f1010ind(0) = 0; f1010ind(1) = 2;
  uvec f1001ind(2); f1001ind(0) = 0; f1001ind(1) = 3;
  uvec f0110ind(2); f0110ind(0) = 1; f0110ind(1) = 2;
  uvec f0101ind(2); f0101ind(0) = 1; f0101ind(1) = 3;
  uvec f0011ind(2); f0011ind(0) = 2; f0011ind(1) = 3;
  
  for (int i = 0; i < len_snp_index; i++){
    int x = snp_index(i) - 1; // snp index
    int x1 = indmat(x,0) - 1; // original index1
    int x2 = indmat(x,1) - 1; // original index2
    int x3 = indmat(x,2) - 1; // original index3
    int x4 = indmat(x,3) - 1; // original index4
    
    double bh1 = Bh1(x1); double bh2 = Bh2(x2); double bh3 = Bh3(x3); double bh4 = Bh4(x4);
    double sigmasq1 = Sigmasq1(x1); double sigmasq2 = Sigmasq2(x2); double sigmasq3 = Sigmasq3(x3); double sigmasq4 = Sigmasq4(x4);
    double rho12 = rhos12(x); double rho13 = rhos13(x); double rho14 = rhos14(x);
    double rho23 = rhos23(x); double rho24 = rhos24(x); double rho34 = rhos34(x);
    double bj1 = Bjh1(x1); double bj2 = Bjh2(x2); double bj3 = Bjh3(x3); double bj4 = Bjh4(x4);
    
    double mu1 = Cmat1->dot_col(x1, Bjh1) - bj1; double mu2 = Cmat2->dot_col(x2, Bjh2) - bj2;
    double mu3 = Cmat3->dot_col(x3, Bjh3) - bj3; double mu4 = Cmat4->dot_col(x4, Bjh4) - bj4;
    
    // need to restart editing from here:
      rowvec mu(4);
    mu.col(0) = mu1; mu.col(1) = mu2; mu.col(2) = mu3; mu.col(3) = mu4;
    double v11 = hsq1 + sigmasq1; double v22 = hsq2 + sigmasq2; double v33 = hsq3 + sigmasq3; double v44 = hsq4 + sigmasq4;
    double v12 = sqrt(hsq1*hsq2) * rho12; double v13 = sqrt(hsq1*hsq3) * rho13; double v14 = sqrt(hsq1*hsq4) * rho14;
    double v23 = sqrt(hsq2*hsq3) * rho23; double v24 = sqrt(hsq2*hsq4) * rho24; double v34 = sqrt(hsq3*hsq4) * rho34;
    double s1 = hsq1 + sigmasq1; double s2 = hsq2 + sigmasq2; double s3 = hsq3 + sigmasq3; double s4 = hsq4 + sigmasq4;
    mat v(4,4);
    v(0,0) = v11; v(1,1) = v22; v(2,2) = v33; v(3,3) = v44;
    v(1,0) = v(0,1) = v12; v(2,0) = v(0,2) = v13; v(3,0) = v(0,3) = v14;
    v(2,1) = v(1,2) = v23; v(3,1) = v(1,3) = v24; v(3,2) = v(2,3) = v34;
    
    mat bh(1,4);
    bh(0,0) = bh1; bh(0,1) = bh2; bh(0,2) = bh3; bh(0,3) = bh4;
    
    vec f1111tem = p1111 * dmvnormcpp(bh, mu, v);
    double f1111 = f1111tem(0);
    vec f1110tem = dmvnormcpp(bh.cols(f1110ind), mu.cols(f1110ind), v.submat(f1110ind,f1110ind));
    double f1110 = p1110 * f1110tem(0) * R::dnorm(bh4,mu4,sqrt(sigmasq4),0);
    vec f1101tem = dmvnormcpp(bh.cols(f1101ind), mu.cols(f1101ind), v.submat(f1101ind,f1101ind));
    double f1101 = p1101 * f1101tem(0) * R::dnorm(bh3,mu3,sqrt(sigmasq3),0);
    vec f1011tem = dmvnormcpp(bh.cols(f1011ind), mu.cols(f1011ind), v.submat(f1011ind,f1011ind));
    double f1011 = p1011 * f1011tem(0) * R::dnorm(bh2,mu2,sqrt(sigmasq2),0);
    vec f0111tem = dmvnormcpp(bh.cols(f0111ind), mu.cols(f0111ind), v.submat(f0111ind,f0111ind));
    double f0111 = p0111 * f0111tem(0) * R::dnorm(bh1,mu1,sqrt(sigmasq1),0);
    
    vec f1100tem = dmvnormcpp(bh.cols(f1100ind), mu.cols(f1100ind), v.submat(f1100ind,f1100ind));
    double f1100 = p1100 * f1100tem(0) * R::dnorm(bh3,mu3,sqrt(sigmasq3),0) * R::dnorm(bh4,mu4,sqrt(sigmasq4),0);
    vec f1010tem = dmvnormcpp(bh.cols(f1010ind), mu.cols(f1010ind), v.submat(f1010ind,f1010ind));
    double f1010 = p1010 * f1010tem(0) * R::dnorm(bh2,mu2,sqrt(sigmasq2),0) * R::dnorm(bh4,mu4,sqrt(sigmasq4),0);
    vec f1001tem = dmvnormcpp(bh.cols(f1001ind), mu.cols(f1001ind), v.submat(f1001ind,f1001ind));
    double f1001 = p1001 * f1001tem(0) * R::dnorm(bh2,mu2,sqrt(sigmasq2),0) * R::dnorm(bh3,mu3,sqrt(sigmasq3),0);
    vec f0110tem = dmvnormcpp(bh.cols(f0110ind), mu.cols(f0110ind), v.submat(f0110ind,f0110ind));
    double f0110 = p0110 * f0110tem(0) * R::dnorm(bh1,mu1,sqrt(sigmasq1),0) * R::dnorm(bh4,mu4,sqrt(sigmasq4),0);
    vec f0101tem = dmvnormcpp(bh.cols(f0101ind), mu.cols(f0101ind), v.submat(f0101ind,f0101ind));
    double f0101 = p0101 * f0101tem(0) * R::dnorm(bh1,mu1,sqrt(sigmasq1),0) * R::dnorm(bh3,mu3,sqrt(sigmasq3),0);
    vec f0011tem = dmvnormcpp(bh.cols(f0011ind), mu.cols(f0011ind), v.submat(f0011ind,f0011ind));
    double f0011 = p0011 * f0011tem(0) * R::dnorm(bh1,mu1,sqrt(sigmasq1),0) * R::dnorm(bh2,mu2,sqrt(sigmasq2),0);
    
    double f1000 = p1000 * R::dnorm(bh1,mu1,sqrt(s1),0) * R::dnorm(bh2,mu2,sqrt(sigmasq2),0) * R::dnorm(bh3,mu3,sqrt(sigmasq3),0) * R::dnorm(bh4,mu4,sqrt(sigmasq4),0);
    double f0100 = p0100 * R::dnorm(bh2,mu2,sqrt(s2),0) * R::dnorm(bh1,mu1,sqrt(sigmasq1),0) * R::dnorm(bh3,mu3,sqrt(sigmasq3),0) * R::dnorm(bh4,mu4,sqrt(sigmasq4),0);
    double f0010 = p0010 * R::dnorm(bh3,mu3,sqrt(s3),0) * R::dnorm(bh1,mu1,sqrt(sigmasq1),0) * R::dnorm(bh2,mu2,sqrt(sigmasq2),0) * R::dnorm(bh4,mu4,sqrt(sigmasq4),0);
    double f0001 = p0001 * R::dnorm(bh4,mu4,sqrt(s4),0) * R::dnorm(bh1,mu1,sqrt(sigmasq1),0) * R::dnorm(bh2,mu2,sqrt(sigmasq2),0) * R::dnorm(bh3,mu3,sqrt(sigmasq3),0);
    double f0000 = p0000 * R::dnorm(bh1,mu1,sqrt(sigmasq1),0) * R::dnorm(bh2,mu2,sqrt(sigmasq2),0) * R::dnorm(bh3,mu3,sqrt(sigmasq3),0) * R::dnorm(bh4,mu4,sqrt(sigmasq4),0);
    double p_d = f1111 + f1110 + f1101 + f1011 + f1110 + f1100 + f1010 + f1001 + f0110 + f0101 + f0011 + f1000 + f0100 + f0010 + f0001 + f0000;
    
    if (p_d>0){
      vec fs(16);
      fs(0) = f1111/p_d; fs(1) = f1110/p_d; fs(2) = f1101/p_d; fs(3) = f1011/p_d; fs(4) = f0111/p_d;
      fs(5) = f1100/p_d; fs(6) = f1010/p_d; fs(7) = f1001/p_d; fs(8) = f0110/p_d; fs(9) = f0101/p_d; fs(10) = f0011/p_d;
      fs(11) = f1000/p_d; fs(12) = f0100/p_d; fs(13) = f0010/p_d; fs(14) = f0001/p_d; fs(15) = f0000/p_d;
      ptilde.col(i) = fs; //
        
        mat v0(4,4);
      v0(0,0) = hsq1; v0(1,1) = hsq2; v0(2,2) = hsq3; v0(3,3) = hsq4;
      v0(0,1) = v0(1,0) = sqrt(hsq1*hsq2)*rho12; v0(0,2) = v0(2,0) = sqrt(hsq1*hsq3)*rho13; v0(0,3) = v0(3,0) = sqrt(hsq1*hsq4)*rho14;
      v0(1,2) = v0(2,1) = sqrt(hsq2*hsq3)*rho23; v0(1,3) = v0(3,1) = sqrt(hsq2*hsq4)*rho24; v0(2,3) = v0(3,2) = sqrt(hsq3*hsq4)*rho34;
      mat v0inv = inv(v0);
      v0inv(0,0) = v0inv(0,0) + 1/sigmasq1; v0inv(1,1) = v0inv(1,1) + 1/sigmasq2;
      v0inv(2,2) = v0inv(2,2) + 1/sigmasq3; v0inv(3,3) = v0inv(3,3) + 1/sigmasq4;
      mat postvar = inv(v0inv);
      mat meantem(4,1);
      meantem = bh.t() - mu.t();
      meantem(0,0) = meantem(0,0)/sigmasq1; meantem(1,0) = meantem(1,0)/sigmasq2;
      meantem(2,0) = meantem(2,0)/sigmasq3; meantem(3,0) = meantem(3,0)/sigmasq4;
      mat postmean = postvar * meantem;
      
      mat v0_12 = v0.submat(ind12,ind12); mat v0_12inv = inv(v0_12);
      v0_12inv(0,0) = v0_12inv(0,0) + 1/sigmasq1; v0_12inv(1,1) = v0_12inv(1,1) + 1/sigmasq2;
      mat postvar12 = inv(v0_12inv); mat meantem12(2,1);
      meantem12 = meantem.rows(ind12); mat postmean12 = postvar12 * meantem12;
      
      mat v0_13 = v0.submat(ind13,ind13); mat v0_13inv = inv(v0_13);
      v0_13inv(0,0) = v0_13inv(0,0) + 1/sigmasq1; v0_13inv(1,1) = v0_13inv(1,1) + 1/sigmasq3;
      mat postvar13 = inv(v0_13inv); mat meantem13(2,1);
      meantem13 = meantem.rows(ind13); mat postmean13 = postvar13 * meantem13;
      
      mat v0_14 = v0.submat(ind14,ind14); mat v0_14inv = inv(v0_14);
      v0_14inv(0,0) = v0_14inv(0,0) + 1/sigmasq1; v0_14inv(1,1) = v0_14inv(1,1) + 1/sigmasq4;
      mat postvar14 = inv(v0_14inv); mat meantem14(2,1);
      meantem14 = meantem.rows(ind14); mat postmean14 = postvar14 * meantem14;
      
      mat v0_23 = v0.submat(ind23,ind23); mat v0_23inv = inv(v0_23);
      v0_23inv(0,0) = v0_23inv(0,0) + 1/sigmasq2; v0_23inv(1,1) = v0_23inv(1,1) + 1/sigmasq3;
      mat postvar23 = inv(v0_23inv); mat meantem23(2,1);
      meantem23 = meantem.rows(ind23); mat postmean23 = postvar23 * meantem23;
      
      mat v0_24 = v0.submat(ind24,ind24); mat v0_24inv = inv(v0_24);
      v0_24inv(0,0) = v0_24inv(0,0) + 1/sigmasq2; v0_24inv(1,1) = v0_24inv(1,1) + 1/sigmasq4;
      mat postvar24 = inv(v0_24inv); mat meantem24(2,1);
      meantem24 = meantem.rows(ind24); mat postmean24 = postvar24 * meantem24;
      
      mat v0_34 = v0.submat(ind34,ind34); mat v0_34inv = inv(v0_34);
      v0_34inv(0,0) = v0_34inv(0,0) + 1/sigmasq3; v0_34inv(1,1) = v0_34inv(1,1) + 1/sigmasq4;
      mat postvar34 = inv(v0_34inv); mat meantem34(2,1);
      meantem34 = meantem.rows(ind34); mat postmean34 = postvar34 * meantem34;
      
      mat v0_123 = v0.submat(ind123,ind123); mat v0_123inv = inv(v0_123);
      v0_123inv(0,0) = v0_123inv(0,0) + 1/sigmasq1; v0_123inv(1,1) = v0_123inv(1,1) + 1/sigmasq2; v0_123inv(2,2) = v0_123inv(2,2) + 1/sigmasq3;
      mat postvar123 = inv(v0_123inv); mat meantem123(3,1);
      meantem123 = meantem.rows(ind123); mat postmean123 = postvar123 * meantem123;
      
      mat v0_124 = v0.submat(ind124,ind124); mat v0_124inv = inv(v0_124);
      v0_124inv(0,0) = v0_124inv(0,0) + 1/sigmasq1; v0_124inv(1,1) = v0_124inv(1,1) + 1/sigmasq2; v0_124inv(2,2) = v0_124inv(2,2) + 1/sigmasq4;
      mat postvar124 = inv(v0_124inv); mat meantem124(3,1);
      meantem124 = meantem.rows(ind124); mat postmean124 = postvar124 * meantem124;
      
      mat v0_134 = v0.submat(ind134,ind134); mat v0_134inv = inv(v0_134);
      v0_134inv(0,0) = v0_134inv(0,0) + 1/sigmasq1; v0_134inv(1,1) = v0_134inv(1,1) + 1/sigmasq3; v0_134inv(2,2) = v0_134inv(2,2) + 1/sigmasq4;
      mat postvar134 = inv(v0_134inv); mat meantem134(3,1);
      meantem134 = meantem.rows(ind134); mat postmean134 = postvar134 * meantem134;
      
      mat v0_234 = v0.submat(ind234,ind234); mat v0_234inv = inv(v0_234);
      v0_234inv(0,0) = v0_234inv(0,0) + 1/sigmasq2; v0_234inv(1,1) = v0_234inv(1,1) + 1/sigmasq3; v0_234inv(2,2) = v0_234inv(2,2) + 1/sigmasq4;
      mat postvar234 = inv(v0_234inv); mat meantem234(3,1);
      meantem234 = meantem.rows(ind234); mat postmean234 = postvar234 * meantem234;
      
      double postvar1 = 1/(1/sigmasq1 + 1/hsq1); double postmean1 = ((bh1 - mu1)/sigmasq1) * postvar1;
      double postvar2 = 1/(1/sigmasq2 + 1/hsq2); double postmean2 = ((bh2 - mu2)/sigmasq2) * postvar2;
      double postvar3 = 1/(1/sigmasq3 + 1/hsq3); double postmean3 = ((bh3 - mu3)/sigmasq3) * postvar3;
      double postvar4 = 1/(1/sigmasq4 + 1/hsq4); double postmean4 = ((bh4 - mu4)/sigmasq4) * postvar4;
      
      NumericVector g = runif(1); double gamma = g(0);
      
      // error: chol(): decomposition failed
      // 1 1 1 1
      if (gamma < fs(0)){
        bjnew.col(i) = rmvncpp(1, postmean, postvar);
      }
      uvec indi(1); indi(0) = i;
      // 1 1 1 0
      if ((gamma >= fs(0))&(gamma < (fs(0)+fs(1)))){
        bjnew.submat(ind123,indi) = rmvncpp(1, postmean123, postvar123);
      }
      // 1 1 0 1
      if ((gamma >= (fs(0)+fs(1)))&(gamma < (fs(0)+fs(1)+fs(2)))){
        bjnew.submat(ind124,indi) = rmvncpp(1, postmean124, postvar124);
      }
      // 1 0 1 1
      if ((gamma >= (fs(0)+fs(1)+fs(2)))&(gamma < (fs(0)+fs(1)+fs(2)+fs(3)))){
        bjnew.submat(ind134,indi) = rmvncpp(1, postmean134, postvar134);
      }
      // 0 1 1 1
      if ((gamma >= (fs(0)+fs(1)+fs(2)+fs(3)))&(gamma < (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)))){
        bjnew.submat(ind234,indi) = rmvncpp(1, postmean234, postvar234);
      }
      // 1 1 0 0
      if ((gamma >= (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)))&(gamma < (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)))){
        bjnew.submat(ind12,indi) = rmvncpp(1, postmean12, postvar12);
      }
      // 1 0 1 0
      if ((gamma >= (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)))&(gamma < (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)))){
        bjnew.submat(ind13,indi) = rmvncpp(1, postmean13, postvar13);
      }
      // 1 0 0 1
      if ((gamma >= (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)))&(gamma < (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)))){
        bjnew.submat(ind14,indi) = rmvncpp(1, postmean14, postvar14);
      }
      // 0 1 1 0
      if ((gamma >= (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)))&(gamma < (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)+fs(8)))){
        bjnew.submat(ind23,indi) = rmvncpp(1, postmean23, postvar23);
      }
      // 0 1 0 1
      if ((gamma >= (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)+fs(8)))&(gamma < (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)+fs(8)+fs(9)))){
        bjnew.submat(ind24,indi) = rmvncpp(1, postmean24, postvar24);
      }
      // 0 0 1 1
      if ((gamma >= (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)+fs(8)+fs(9)))&(gamma < (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)+fs(8)+fs(9)+fs(10)))){
        bjnew.submat(ind34,indi) = rmvncpp(1, postmean34, postvar34);
      }
      // 1 0 0 0
      if ((gamma >= (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)+fs(8)+fs(9)+fs(10)))&(gamma < (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)+fs(8)+fs(9)+fs(10)+fs(11)))){
        NumericVector tem = rnorm(1, postmean1, sqrt(postvar1));
        bjnew(0,i) = tem(0);
      }
      // 0 1 0 0
      if ((gamma >= (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)+fs(8)+fs(9)+fs(10)+fs(11)))&(gamma < (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)+fs(8)+fs(9)+fs(10)+fs(11)+fs(12)))){
        NumericVector tem = rnorm(1, postmean2, sqrt(postvar2));
        bjnew(1,i) = tem(0);
      }
      // 0 0 1 0
      if ((gamma >= (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)+fs(8)+fs(9)+fs(10)+fs(11)+fs(12)))&(gamma < (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)+fs(8)+fs(9)+fs(10)+fs(11)+fs(12)+fs(13)))){
        NumericVector tem = rnorm(1, postmean3, sqrt(postvar3));
        bjnew(2,i) = tem(0);
      }
      // 0 0 0 1
      if ((gamma >= (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)+fs(8)+fs(9)+fs(10)+fs(11)+fs(12)+fs(13)))&(gamma < (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)+fs(8)+fs(9)+fs(10)+fs(11)+fs(12)+fs(13)+fs(14)))){
        NumericVector tem = rnorm(1, postmean4, sqrt(postvar4));
        bjnew(3,i) = tem(0);
      }
      double pm1tem = fs(0)*postmean(0,0) + fs(1)*postmean123(0,0) + fs(2)*postmean124(0,0) + fs(3)*postmean134(0,0) +
        fs(5)*postmean12(0,0) + fs(6)*postmean13(0,0) + fs(7)*postmean14(0,0) + fs(11)*postmean1;
      double pm2tem = fs(0)*postmean(1,0) + fs(1)*postmean123(1,0) + fs(2)*postmean124(1,0) + fs(4)*postmean234(0,0) +
        fs(5)*postmean12(1,0) + fs(8)*postmean23(0,0) + fs(9)*postmean24(0,0) + fs(12)*postmean2;
      double pm3tem = fs(0)*postmean(2,0) + fs(1)*postmean123(2,0) + fs(3)*postmean134(1,0) + fs(4)*postmean234(1,0) +
        fs(6)*postmean13(1,0) + fs(8)*postmean23(1,0) + fs(10)*postmean34(0,0) + fs(13)*postmean3;
      double pm4tem = fs(0)*postmean(3,0) + fs(2)*postmean124(2,0) + fs(3)*postmean134(2,0) + fs(4)*postmean234(2,0) +
        fs(7)*postmean14(1,0) + fs(9)*postmean24(1,0) + fs(10)*postmean34(1,0) + fs(14)*postmean4;
      pm(0,i) = pm1tem;
      pm(1,i) = pm2tem;
      pm(2,i) = pm3tem;
      pm(3,i) = pm4tem;
      // sparse
      if (sparse(0) == 1){
        double fs1 = fs(0) + fs(1) + fs(2) + fs(3) + fs(5) + fs(6) + fs(7) + fs(11);
        double fs1_0 = p1111 + p1110 + p1101 + p1011 + p1100 + p1010 + p1001 + p1000;
        if (fs1 < fs1_0){
          bjnew(0,i) = 0;
          pm(0,i) = 0;
          postmean(0,0) = 0;
          postmean1 = 0;
        }
      }
      if (sparse(1) == 1){
        double fs2 = fs(0) + fs(1) + fs(2) + fs(4) + fs(5) + fs(8) + fs(9) + fs(12);
        double fs2_0 = p1111 + p1110 + p1101 + p0111 + p1100 + p0110 + p0101 + p0100;
        if (fs2 < fs2_0){
          bjnew(1,i) = 0;
          pm(1,i) = 0;
          postmean(1,0) = 0;
          postmean2 = 0;
        }
      }
      if (sparse(2) == 1){
        double fs3 = fs(0) + fs(1) + fs(3) + fs(4) + fs(6) + fs(8) + fs(10) + fs(13);
        double fs3_0 = p1111 + p1110 + p1011 + p0111 + p1010 + p0110 + p0011 + p0010;
        if (fs3 < fs3_0){
          bjnew(2,i) = 0;
          pm(2,i) = 0;
          postmean(2,0) = 0;
          postmean3 = 0;
        }
      }
      if (sparse(3) == 1){
        double fs4 = fs(0) + fs(2) + fs(3) + fs(4) + fs(7) + fs(9) + fs(10) + fs(14);
        double fs4_0 = p1111 + p1101 + p1011 + p0111 + p1001 + p0101 + p0011 + p0001;
        if (fs4 < fs4_0){
          bjnew(3,i) = 0;
          pm(3,i) = 0;
          postmean(3,0) = 0;
          postmean4 = 0;
        }
      }
    }
    Bjh1(x1) = bjnew(0,i);
    Bjh2(x2) = bjnew(1,i);
    Bjh3(x3) = bjnew(2,i);
    Bjh4(x4) = bjnew(3,i);
  }
  //return List::create(_["ptilde"] = ptilde, _["bjnew"] = bjnew, _["pm"] = pm);
  return List::create(_["bjnew"] = bjnew, _["pm"] = pm); //_["ptilde"] = ptilde,
}



// [[Rcpp::export]]
List update_beta5(vec snp_index, mat indmat,
                  Environment cormat1, Environment cormat2, Environment cormat3, Environment cormat4, Environment cormat5,
                  vec Bh1, vec Bh2, vec Bh3, vec Bh4, vec Bh5, 
                  vec Bjh1, vec Bjh2, vec Bjh3, vec Bjh4, vec Bjh5,
                  const NumericVector& Sigmasq1, const NumericVector& Sigmasq2, const NumericVector& Sigmasq3, const NumericVector& Sigmasq4, const NumericVector& Sigmasq5,
                  double hsq1, double hsq2, double hsq3, double hsq4, double hsq5,
                  const NumericVector& rhos12, const NumericVector& rhos13, const NumericVector& rhos14, const NumericVector& rhos15,
                  const NumericVector& rhos23, const NumericVector& rhos24, const NumericVector& rhos25, 
                  const NumericVector& rhos34, const NumericVector& rhos35, const NumericVector& rhos45,
                  vec pall, vec sparse){
  XPtr<SFBM> Cmat1 = cormat1["address"]; XPtr<SFBM> Cmat2 = cormat2["address"];
  XPtr<SFBM> Cmat3 = cormat3["address"]; XPtr<SFBM> Cmat4 = cormat4["address"]; XPtr<SFBM> Cmat5 = cormat5["address"];
  
  int m1 = Bh1.n_elem; int m2 = Bh2.n_elem; int m3 = Bh3.n_elem; int m4 = Bh4.n_elem; int m5 = Bh5.n_elem;
  myassert_size(Cmat1->nrow(), m1); myassert_size(Cmat1->ncol(), m1);
  myassert_size(Cmat2->nrow(), m2); myassert_size(Cmat2->ncol(), m2);
  myassert_size(Cmat3->nrow(), m3); myassert_size(Cmat3->ncol(), m3);
  myassert_size(Cmat4->nrow(), m4); myassert_size(Cmat4->ncol(), m4);
  myassert_size(Cmat5->nrow(), m5); myassert_size(Cmat5->ncol(), m5);
  
  double len_snp_index = snp_index.n_elem;
  mat ptilde(32,len_snp_index);
  ptilde.fill(0);
  mat bjnew(5,len_snp_index);
  bjnew.fill(0);
  mat pm(5,len_snp_index);
  pm.fill(0);
  
  uvec ind12(2); ind12(0) = 0; ind12(1) = 1;
  uvec ind13(2); ind13(0) = 0; ind13(1) = 2;
  uvec ind14(2); ind14(0) = 0; ind14(1) = 3;
  uvec ind15(2); ind15(0) = 0; ind15(1) = 4;
  uvec ind23(2); ind23(0) = 1; ind23(1) = 2;
  uvec ind24(2); ind24(0) = 1; ind24(1) = 3;
  uvec ind25(2); ind25(0) = 1; ind25(1) = 4;
  uvec ind34(2); ind34(0) = 2; ind34(1) = 3;
  uvec ind35(2); ind35(0) = 2; ind35(1) = 4;
  uvec ind45(2); ind45(0) = 3; ind45(1) = 4;
  uvec ind123(3); ind123(0) = 0; ind123(1) = 1; ind123(2) = 2;
  uvec ind124(3); ind124(0) = 0; ind124(1) = 1; ind124(2) = 3;
  uvec ind125(3); ind125(0) = 0; ind125(1) = 1; ind125(2) = 4;
  uvec ind134(3); ind134(0) = 0; ind134(1) = 2; ind134(2) = 3;
  uvec ind135(3); ind135(0) = 0; ind135(1) = 2; ind135(2) = 4;
  uvec ind145(3); ind145(0) = 0; ind145(1) = 3; ind145(2) = 4;
  uvec ind234(3); ind234(0) = 1; ind234(1) = 2; ind234(2) = 3;
  uvec ind235(3); ind235(0) = 1; ind235(1) = 2; ind235(2) = 4;
  uvec ind245(3); ind245(0) = 1; ind245(1) = 3; ind245(2) = 4;
  uvec ind345(3); ind345(0) = 2; ind345(1) = 3; ind345(2) = 4;
  
  uvec ind1234(4); ind1234(0) = 0; ind1234(1) = 1; ind1234(2) = 2; ind1234(3) = 3;
  uvec ind1235(4); ind1235(0) = 0; ind1235(1) = 1; ind1235(2) = 2; ind1235(3) = 4;
  uvec ind1245(4); ind1245(0) = 0; ind1245(1) = 1; ind1245(2) = 3; ind1245(3) = 4;
  uvec ind1345(4); ind1345(0) = 0; ind1345(1) = 2; ind1345(2) = 3; ind1345(3) = 4;
  uvec ind2345(4); ind2345(0) = 1; ind2345(1) = 2; ind2345(2) = 3; ind2345(3) = 4;
  
  uvec ind12345(5); ind12345(0) = 0; ind12345(1) = 1; ind12345(2) = 2; ind12345(3) = 3; ind12345(4) = 4;
  
  uvec f11110ind(4); f11110ind(0) = 0; f11110ind(1) = 1; f11110ind(2) = 2; f11110ind(3) = 3;
  uvec f11101ind(4); f11101ind(0) = 0; f11101ind(1) = 1; f11101ind(2) = 2; f11101ind(3) = 4;
  uvec f11011ind(4); f11011ind(0) = 0; f11011ind(1) = 1; f11011ind(2) = 3; f11011ind(3) = 4;
  uvec f10111ind(4); f10111ind(0) = 0; f10111ind(1) = 2; f10111ind(2) = 3; f10111ind(3) = 4;
  uvec f01111ind(4); f01111ind(0) = 1; f01111ind(1) = 2; f01111ind(2) = 3; f01111ind(3) = 4;
  
  uvec f01110ind(3); f01110ind(0) = 1; f01110ind(1) = 2; f01110ind(2) = 3;
  uvec f01101ind(3); f01101ind(0) = 1; f01101ind(1) = 2; f01101ind(2) = 4;
  uvec f01011ind(3); f01011ind(0) = 1; f01011ind(1) = 3; f01011ind(2) = 4;
  uvec f00111ind(3); f00111ind(0) = 2; f00111ind(1) = 3; f00111ind(2) = 4;
  
  uvec f11100ind(3); f11100ind(0) = 0; f11100ind(1) = 1; f11100ind(2) = 2;
  uvec f11010ind(3); f11010ind(0) = 0; f11010ind(1) = 1; f11010ind(2) = 3;
  uvec f11001ind(3); f11001ind(0) = 0; f11001ind(1) = 1; f11001ind(2) = 4;
  uvec f10110ind(3); f10110ind(0) = 0; f10110ind(1) = 2; f10110ind(2) = 3;
  uvec f10101ind(3); f10101ind(0) = 0; f10101ind(1) = 2; f10101ind(2) = 4;
  uvec f10011ind(3); f10011ind(0) = 0; f10011ind(1) = 3; f10011ind(2) = 4;
  
  uvec f00011ind(2); f00011ind(0) = 3; f00011ind(1) = 4;
  uvec f00101ind(2); f00101ind(0) = 2; f00101ind(1) = 4;
  uvec f01001ind(2); f01001ind(0) = 1; f01001ind(1) = 4;
  uvec f10001ind(2); f10001ind(0) = 0; f10001ind(1) = 4;
  uvec f00110ind(2); f00110ind(0) = 2; f00110ind(1) = 3;
  uvec f01010ind(2); f01010ind(0) = 1; f01010ind(1) = 3;
  uvec f10010ind(2); f10010ind(0) = 0; f10010ind(1) = 3;
  uvec f01100ind(2); f01100ind(0) = 1; f01100ind(1) = 2;
  uvec f10100ind(2); f10100ind(0) = 0; f10100ind(1) = 2;
  uvec f11000ind(2); f11000ind(0) = 0; f11000ind(1) = 1;
  
  double p11111 = pall(0); double p11110 = pall(1); double p11101 = pall(2); 
  double p11011 = pall(3); double p10111 = pall(4);
  double p11100 = pall(5); double p11010 = pall(6); double p11001 = pall(7); 
  double p10110 = pall(8); double p10101 = pall(9); double p10011 = pall(10); 
  double p11000 = pall(11); double p10100 = pall(12); double p10010 = pall(13); 
  double p10001 = pall(14); double p10000 = pall(15); double p01111 = pall(16); 
  double p01110 = pall(17); double p01101 = pall(18); double p01011 = pall(19); 
  double p00111 = pall(20); double p01100 = pall(21); double p01010 = pall(22); 
  double p01001 = pall(23); double p00110 = pall(24); double p00101 = pall(25);
  double p00011 = pall(26); double p01000 = pall(27); double p00100 = pall(28); 
  double p00010 = pall(29); double p00001 = pall(30); double p00000 = pall(31); 
  for (int i = 0; i < len_snp_index; i++){
    int x = snp_index(i) - 1; // snp index
    int x1 = indmat(x,0) - 1; // original index1
    int x2 = indmat(x,1) - 1; // original index2
    int x3 = indmat(x,2) - 1; // original index3
    int x4 = indmat(x,3) - 1; // original index4
    int x5 = indmat(x,4) - 1; // original index5
    
    double bh1 = Bh1(x1); double bh2 = Bh2(x2); double bh3 = Bh3(x3); double bh4 = Bh4(x4); double bh5 = Bh5(x5);
    double sigmasq1 = Sigmasq1(x1); double sigmasq2 = Sigmasq2(x2);
    double sigmasq3 = Sigmasq3(x3); double sigmasq4 = Sigmasq4(x4); double sigmasq5 = Sigmasq5(x5);
    double rho12 = rhos12(x); double rho13 = rhos13(x); double rho14 = rhos14(x); double rho15 = rhos15(x);
    double rho23 = rhos23(x); double rho24 = rhos24(x); double rho25 = rhos25(x); 
    double rho34 = rhos34(x); double rho35 = rhos35(x); double rho45 = rhos45(x);
    double bj1 = Bjh1(x1); double bj2 = Bjh2(x2); double bj3 = Bjh3(x3); double bj4 = Bjh4(x4); double bj5 = Bjh5(x5);
    
    double mu1 = Cmat1->dot_col(x1, Bjh1) - bj1; double mu2 = Cmat2->dot_col(x2, Bjh2) - bj2;
    double mu3 = Cmat3->dot_col(x3, Bjh3) - bj3; double mu4 = Cmat4->dot_col(x4, Bjh4) - bj4;
    double mu5 = Cmat5->dot_col(x5, Bjh5) - bj5;
    
    // need to restart editing from here:
      rowvec mu(5);
    mu.col(0) = mu1; mu.col(1) = mu2; mu.col(2) = mu3; mu.col(3) = mu4; mu.col(4) = mu5;
    double v11 = hsq1 + sigmasq1; double v22 = hsq2 + sigmasq2; double v33 = hsq3 + sigmasq3; double v44 = hsq4 + sigmasq4; double v55 = hsq5 + sigmasq5;
    double v12 = sqrt(hsq1*hsq2) * rho12; double v13 = sqrt(hsq1*hsq3) * rho13; double v14 = sqrt(hsq1*hsq4) * rho14; double v15 = sqrt(hsq1*hsq5) * rho15;
    double v23 = sqrt(hsq2*hsq3) * rho23; double v24 = sqrt(hsq2*hsq4) * rho24; double v25 = sqrt(hsq2*hsq5) * rho25;  
    double v34 = sqrt(hsq3*hsq4) * rho34; double v35 = sqrt(hsq3*hsq5) * rho35; double v45 = sqrt(hsq4*hsq5) * rho45;
    double s1 = hsq1 + sigmasq1; double s2 = hsq2 + sigmasq2; double s3 = hsq3 + sigmasq3; double s4 = hsq4 + sigmasq4; double s5 = hsq5 + sigmasq5;
    mat v(5,5);
    v(0,0) = v11; v(1,1) = v22; v(2,2) = v33; v(3,3) = v44; v(4,4) = v55;
    v(1,0) = v(0,1) = v12; v(2,0) = v(0,2) = v13; v(3,0) = v(0,3) = v14; v(4,0) = v(0,4) = v15;
    v(2,1) = v(1,2) = v23; v(3,1) = v(1,3) = v24; v(4,1) = v(1,4) = v25; 
    v(3,2) = v(2,3) = v34; v(4,2) = v(2,4) = v35; v(4,3) = v(3,4) = v45;
    
    mat bh(1,5);
    bh(0,0) = bh1; bh(0,1) = bh2; bh(0,2) = bh3; bh(0,3) = bh4; bh(0,4) = bh5;
    
    vec f11111tem = p11111 * dmvnormcpp(bh, mu, v);
    double f11111 = f11111tem(0);
    vec f11110tem = dmvnormcpp(bh.cols(f11110ind), mu.cols(f11110ind), v.submat(f11110ind,f11110ind));
    double f11110 = p11110 * f11110tem(0) * R::dnorm(bh5,mu5,sqrt(sigmasq5),0);
    vec f11101tem = dmvnormcpp(bh.cols(f11101ind), mu.cols(f11101ind), v.submat(f11101ind,f11101ind));
    double f11101 = p11101 * f11101tem(0) * R::dnorm(bh4,mu4,sqrt(sigmasq4),0);
    vec f11011tem = dmvnormcpp(bh.cols(f11011ind), mu.cols(f11011ind), v.submat(f11011ind,f11011ind));
    double f11011 = p11011 * f11011tem(0) * R::dnorm(bh3,mu3,sqrt(sigmasq3),0);
    vec f10111tem = dmvnormcpp(bh.cols(f10111ind), mu.cols(f10111ind), v.submat(f10111ind,f10111ind));
    double f10111 = p10111 * f10111tem(0) * R::dnorm(bh2,mu2,sqrt(sigmasq2),0);
    vec f01111tem = dmvnormcpp(bh.cols(f01111ind), mu.cols(f01111ind), v.submat(f01111ind,f01111ind));
    double f01111 = p01111 * f01111tem(0) * R::dnorm(bh1,mu1,sqrt(sigmasq1),0);
    
    vec f11100tem = dmvnormcpp(bh.cols(f11100ind), mu.cols(f11100ind), v.submat(f11100ind,f11100ind));
    double f11100 = p11100 * f11100tem(0) * R::dnorm(bh4,mu4,sqrt(sigmasq4),0) * R::dnorm(bh5,mu5,sqrt(sigmasq5),0);
    vec f11010tem = dmvnormcpp(bh.cols(f11010ind), mu.cols(f11010ind), v.submat(f11010ind,f11010ind));
    double f11010 = p11010 * f11010tem(0) * R::dnorm(bh3,mu3,sqrt(sigmasq3),0) * R::dnorm(bh5,mu5,sqrt(sigmasq5),0);
    vec f11001tem = dmvnormcpp(bh.cols(f11001ind), mu.cols(f11001ind), v.submat(f11001ind,f11001ind));
    double f11001 = p11001 * f11001tem(0) * R::dnorm(bh3,mu3,sqrt(sigmasq3),0) * R::dnorm(bh4,mu4,sqrt(sigmasq4),0);
    vec f10110tem = dmvnormcpp(bh.cols(f10110ind), mu.cols(f10110ind), v.submat(f10110ind,f10110ind));
    double f10110 = p10110 * f10110tem(0) * R::dnorm(bh2,mu2,sqrt(sigmasq2),0) * R::dnorm(bh5,mu5,sqrt(sigmasq5),0);
    vec f10101tem = dmvnormcpp(bh.cols(f10101ind), mu.cols(f10101ind), v.submat(f10101ind,f10101ind));
    double f10101 = p10101 * f10101tem(0) * R::dnorm(bh2,mu2,sqrt(sigmasq2),0) * R::dnorm(bh4,mu4,sqrt(sigmasq4),0);
    vec f10011tem = dmvnormcpp(bh.cols(f10011ind), mu.cols(f10011ind), v.submat(f10011ind,f10011ind));
    double f10011 = p10011 * f10011tem(0) * R::dnorm(bh2,mu2,sqrt(sigmasq2),0) * R::dnorm(bh3,mu3,sqrt(sigmasq3),0);
    vec f01110tem = dmvnormcpp(bh.cols(f01110ind), mu.cols(f01110ind), v.submat(f01110ind,f01110ind));
    double f01110 = p01110 * f01110tem(0) * R::dnorm(bh1,mu1,sqrt(sigmasq1),0) * R::dnorm(bh5,mu5,sqrt(sigmasq5),0);
    vec f01101tem = dmvnormcpp(bh.cols(f01101ind), mu.cols(f01101ind), v.submat(f01101ind,f01101ind));
    double f01101 = p01101 * f01101tem(0) * R::dnorm(bh1,mu1,sqrt(sigmasq1),0) * R::dnorm(bh4,mu4,sqrt(sigmasq4),0);
    vec f01011tem = dmvnormcpp(bh.cols(f01011ind), mu.cols(f01011ind), v.submat(f01011ind,f01011ind));
    double f01011 = p01011 * f01011tem(0) * R::dnorm(bh1,mu1,sqrt(sigmasq1),0) * R::dnorm(bh3,mu3,sqrt(sigmasq3),0);
    vec f00111tem = dmvnormcpp(bh.cols(f00111ind), mu.cols(f00111ind), v.submat(f00111ind,f00111ind));
    double f00111 = p00111 * f00111tem(0) * R::dnorm(bh1,mu1,sqrt(sigmasq1),0) * R::dnorm(bh2,mu2,sqrt(sigmasq2),0);
    
    vec f11000tem = dmvnormcpp(bh.cols(f11000ind), mu.cols(f11000ind), v.submat(f11000ind,f11000ind));
    double f11000 = p11000 * f11000tem(0) * R::dnorm(bh3,mu3,sqrt(sigmasq3),0) * R::dnorm(bh4,mu4,sqrt(sigmasq4),0) * R::dnorm(bh5,mu5,sqrt(sigmasq5),0);
    vec f10100tem = dmvnormcpp(bh.cols(f10100ind), mu.cols(f10100ind), v.submat(f10100ind,f10100ind));
    double f10100 = p10100 * f10100tem(0) * R::dnorm(bh2,mu2,sqrt(sigmasq2),0) * R::dnorm(bh4,mu4,sqrt(sigmasq4),0) * R::dnorm(bh5,mu5,sqrt(sigmasq5),0);
    vec f10010tem = dmvnormcpp(bh.cols(f10010ind), mu.cols(f10010ind), v.submat(f10010ind,f10010ind));
    double f10010 = p10010 * f10010tem(0) * R::dnorm(bh2,mu2,sqrt(sigmasq2),0) * R::dnorm(bh3,mu3,sqrt(sigmasq3),0) * R::dnorm(bh5,mu5,sqrt(sigmasq5),0);
    vec f10001tem = dmvnormcpp(bh.cols(f10001ind), mu.cols(f10001ind), v.submat(f10001ind,f10001ind));
    double f10001 = p10001 * f10001tem(0) * R::dnorm(bh2,mu2,sqrt(sigmasq2),0) * R::dnorm(bh3,mu3,sqrt(sigmasq3),0) * R::dnorm(bh4,mu4,sqrt(sigmasq4),0);
    vec f01100tem = dmvnormcpp(bh.cols(f01100ind), mu.cols(f01100ind), v.submat(f01100ind,f01100ind));
    double f01100 = p01100 * f01100tem(0) * R::dnorm(bh1,mu1,sqrt(sigmasq1),0) * R::dnorm(bh4,mu4,sqrt(sigmasq4),0) * R::dnorm(bh5,mu5,sqrt(sigmasq5),0);
    vec f01010tem = dmvnormcpp(bh.cols(f01010ind), mu.cols(f01010ind), v.submat(f01010ind,f01010ind));
    double f01010 = p01010 * f01010tem(0) * R::dnorm(bh1,mu1,sqrt(sigmasq1),0) * R::dnorm(bh3,mu3,sqrt(sigmasq3),0) * R::dnorm(bh5,mu5,sqrt(sigmasq5),0);
    vec f01001tem = dmvnormcpp(bh.cols(f01001ind), mu.cols(f01001ind), v.submat(f01001ind,f01001ind));
    double f01001 = p01001 * f01001tem(0) * R::dnorm(bh1,mu1,sqrt(sigmasq1),0) * R::dnorm(bh3,mu3,sqrt(sigmasq3),0) * R::dnorm(bh4,mu4,sqrt(sigmasq4),0);
    vec f00110tem = dmvnormcpp(bh.cols(f00110ind), mu.cols(f00110ind), v.submat(f00110ind,f00110ind));
    double f00110 = p00110 * f00110tem(0) * R::dnorm(bh1,mu1,sqrt(sigmasq1),0) * R::dnorm(bh2,mu2,sqrt(sigmasq3),0) * R::dnorm(bh5,mu5,sqrt(sigmasq5),0);
    vec f00101tem = dmvnormcpp(bh.cols(f00101ind), mu.cols(f00101ind), v.submat(f00101ind,f00101ind));
    double f00101 = p00101 * f00101tem(0) * R::dnorm(bh1,mu1,sqrt(sigmasq1),0) * R::dnorm(bh2,mu2,sqrt(sigmasq3),0) * R::dnorm(bh4,mu4,sqrt(sigmasq4),0);
    vec f00011tem = dmvnormcpp(bh.cols(f00011ind), mu.cols(f00011ind), v.submat(f00011ind,f00011ind));
    double f00011 = p00011 * f00011tem(0) * R::dnorm(bh1,mu1,sqrt(sigmasq1),0) * R::dnorm(bh2,mu2,sqrt(sigmasq3),0) * R::dnorm(bh3,mu3,sqrt(sigmasq3),0);
    
    double f10000 = p10000 * R::dnorm(bh1,mu1,sqrt(s1),0) * R::dnorm(bh2,mu2,sqrt(sigmasq2),0) * R::dnorm(bh3,mu3,sqrt(sigmasq3),0) * R::dnorm(bh4,mu4,sqrt(sigmasq4),0) * R::dnorm(bh5,mu5,sqrt(sigmasq5),0);
    double f01000 = p01000 * R::dnorm(bh2,mu2,sqrt(s2),0) * R::dnorm(bh1,mu1,sqrt(sigmasq1),0) * R::dnorm(bh3,mu3,sqrt(sigmasq3),0) * R::dnorm(bh4,mu4,sqrt(sigmasq4),0) * R::dnorm(bh5,mu5,sqrt(sigmasq5),0);
    double f00100 = p00100 * R::dnorm(bh3,mu3,sqrt(s3),0) * R::dnorm(bh1,mu1,sqrt(sigmasq1),0) * R::dnorm(bh2,mu2,sqrt(sigmasq2),0) * R::dnorm(bh4,mu4,sqrt(sigmasq4),0) * R::dnorm(bh5,mu5,sqrt(sigmasq5),0);
    double f00010 = p00010 * R::dnorm(bh4,mu4,sqrt(s4),0) * R::dnorm(bh1,mu1,sqrt(sigmasq1),0) * R::dnorm(bh2,mu2,sqrt(sigmasq2),0) * R::dnorm(bh3,mu3,sqrt(sigmasq3),0) * R::dnorm(bh5,mu5,sqrt(sigmasq5),0);
    double f00001 = p00001 * R::dnorm(bh5,mu5,sqrt(s5),0) * R::dnorm(bh1,mu1,sqrt(sigmasq1),0) * R::dnorm(bh2,mu2,sqrt(sigmasq2),0) * R::dnorm(bh3,mu3,sqrt(sigmasq3),0) * R::dnorm(bh4,mu4,sqrt(sigmasq4),0);
    double f00000 = p00000 * R::dnorm(bh1,mu1,sqrt(sigmasq1),0) * R::dnorm(bh2,mu2,sqrt(sigmasq2),0) * R::dnorm(bh3,mu3,sqrt(sigmasq3),0) * R::dnorm(bh4,mu4,sqrt(sigmasq4),0) * R::dnorm(bh5,mu5,sqrt(sigmasq5),0); 
    double p_d = f11111 + f11110 + f11101 + f11011 + f10111 + f01111 + 
      f11100 + f11010 + f11001 + f10110 + f10101 + f10011 + f01110 + f01101 + f01011 + f00111 +
      f11000 + f10100 + f10010 + f10001 + f01100 + f01010 + f01001 + f00110 + f00101 + f00011 +
      f10000 + f01000 + f00100 + f00010 + f00001 + f00000;
    
    if (p_d>0){
      vec fs(32);
      fs(0) = f11111/p_d; 
      fs(1) = f11110/p_d; fs(2) = f11101/p_d; fs(3) = f11011/p_d; fs(4) = f10111/p_d; fs(5) = f01111/p_d;
      fs(6) = f11100/p_d; fs(7) = f11010/p_d; fs(8) = f11001/p_d; fs(9) = f10110/p_d; fs(10) = f10101/p_d; fs(11) = f10011/p_d;
      fs(12) = f01110/p_d; fs(13) = f01101/p_d; fs(14) = f01011/p_d; fs(15) = f00111/p_d;
      fs(16) = f11000/p_d; fs(17) = f10100/p_d; fs(18) = f10010/p_d; fs(19) = f10001/p_d; 
      fs(20) = f01100/p_d; fs(21) = f01010/p_d; fs(22) = f01001/p_d;
      fs(23) = f00110/p_d; fs(24) = f00101/p_d; fs(25) = f00011/p_d;
      fs(26) = f10000/p_d; fs(27) = f01000/p_d; fs(28) = f00100/p_d; fs(29) = f00010/p_d; fs(30) = f00001/p_d;
      fs(31) = f00000/p_d;
      ptilde.col(i) = fs; //
        double fs_0_20 = fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)+fs(8)+fs(9)+fs(10)+fs(11)+fs(12)+fs(13)+fs(14)+fs(15)+fs(16)+fs(17)+fs(18)+fs(19)+fs(20);
        
        mat v0(5,5);
        v0(0,0) = hsq1; v0(1,1) = hsq2; v0(2,2) = hsq3; v0(3,3) = hsq4; v0(4,4) = hsq5;
        v0(0,1) = v0(1,0) = sqrt(hsq1*hsq2)*rho12; v0(0,2) = v0(2,0) = sqrt(hsq1*hsq3)*rho13; v0(0,3) = v0(3,0) = sqrt(hsq1*hsq4)*rho14; v0(0,4) = v0(4,0) = sqrt(hsq1*hsq5)*rho15;
        v0(1,2) = v0(2,1) = sqrt(hsq2*hsq3)*rho23; v0(1,3) = v0(3,1) = sqrt(hsq2*hsq4)*rho24; v0(1,4) = v0(4,1) = sqrt(hsq2*hsq5)*rho25; 
        v0(2,3) = v0(3,2) = sqrt(hsq3*hsq4)*rho34; v0(2,4) = v0(4,2) = sqrt(hsq3*hsq5)*rho35;
        v0(3,4) = v0(4,3) = sqrt(hsq4*hsq5)*rho45;
        mat v0inv = inv(v0);
        v0inv(0,0) = v0inv(0,0) + 1/sigmasq1; v0inv(1,1) = v0inv(1,1) + 1/sigmasq2;
        v0inv(2,2) = v0inv(2,2) + 1/sigmasq3; v0inv(3,3) = v0inv(3,3) + 1/sigmasq4; v0inv(4,4) = v0inv(4,4) + 1/sigmasq5;
        mat postvar = inv(v0inv);
        mat meantem(5,1);
        meantem = bh.t() - mu.t();
        meantem(0,0) = meantem(0,0)/sigmasq1; meantem(1,0) = meantem(1,0)/sigmasq2;
        meantem(2,0) = meantem(2,0)/sigmasq3; meantem(3,0) = meantem(3,0)/sigmasq4; meantem(4,0) = meantem(4,0)/sigmasq5;
        mat postmean = postvar * meantem;
        
        mat v0_12 = v0.submat(ind12,ind12); mat v0_12inv = inv(v0_12);
        v0_12inv(0,0) = v0_12inv(0,0) + 1/sigmasq1; v0_12inv(1,1) = v0_12inv(1,1) + 1/sigmasq2;
        mat postvar12 = inv(v0_12inv); mat meantem12(2,1);
        meantem12 = meantem.rows(ind12); mat postmean12 = postvar12 * meantem12;
        
        mat v0_13 = v0.submat(ind13,ind13); mat v0_13inv = inv(v0_13);
        v0_13inv(0,0) = v0_13inv(0,0) + 1/sigmasq1; v0_13inv(1,1) = v0_13inv(1,1) + 1/sigmasq3;
        mat postvar13 = inv(v0_13inv); mat meantem13(2,1);
        meantem13 = meantem.rows(ind13); mat postmean13 = postvar13 * meantem13;
        
        mat v0_14 = v0.submat(ind14,ind14); mat v0_14inv = inv(v0_14);
        v0_14inv(0,0) = v0_14inv(0,0) + 1/sigmasq1; v0_14inv(1,1) = v0_14inv(1,1) + 1/sigmasq4;
        mat postvar14 = inv(v0_14inv); mat meantem14(2,1);
        meantem14 = meantem.rows(ind14); mat postmean14 = postvar14 * meantem14;
        
        mat v0_15 = v0.submat(ind15,ind15); mat v0_15inv = inv(v0_15);
        v0_15inv(0,0) = v0_15inv(0,0) + 1/sigmasq1; v0_15inv(1,1) = v0_15inv(1,1) + 1/sigmasq5;
        mat postvar15 = inv(v0_15inv); mat meantem15(2,1);
        meantem15 = meantem.rows(ind15); mat postmean15 = postvar15 * meantem15;
        
        mat v0_23 = v0.submat(ind23,ind23); mat v0_23inv = inv(v0_23);
        v0_23inv(0,0) = v0_23inv(0,0) + 1/sigmasq2; v0_23inv(1,1) = v0_23inv(1,1) + 1/sigmasq3;
        mat postvar23 = inv(v0_23inv); mat meantem23(2,1);
        meantem23 = meantem.rows(ind23); mat postmean23 = postvar23 * meantem23;
        
        mat v0_24 = v0.submat(ind24,ind24); mat v0_24inv = inv(v0_24);
        v0_24inv(0,0) = v0_24inv(0,0) + 1/sigmasq2; v0_24inv(1,1) = v0_24inv(1,1) + 1/sigmasq4;
        mat postvar24 = inv(v0_24inv); mat meantem24(2,1);
        meantem24 = meantem.rows(ind24); mat postmean24 = postvar24 * meantem24;
        
        mat v0_25 = v0.submat(ind25,ind25); mat v0_25inv = inv(v0_25);
        v0_25inv(0,0) = v0_25inv(0,0) + 1/sigmasq2; v0_25inv(1,1) = v0_25inv(1,1) + 1/sigmasq5;
        mat postvar25 = inv(v0_25inv); mat meantem25(2,1);
        meantem25 = meantem.rows(ind25); mat postmean25 = postvar25 * meantem25;
        
        mat v0_34 = v0.submat(ind34,ind34); mat v0_34inv = inv(v0_34);
        v0_34inv(0,0) = v0_34inv(0,0) + 1/sigmasq3; v0_34inv(1,1) = v0_34inv(1,1) + 1/sigmasq4;
        mat postvar34 = inv(v0_34inv); mat meantem34(2,1);
        meantem34 = meantem.rows(ind34); mat postmean34 = postvar34 * meantem34;
        
        mat v0_35 = v0.submat(ind35,ind35); mat v0_35inv = inv(v0_35);
        v0_35inv(0,0) = v0_35inv(0,0) + 1/sigmasq3; v0_35inv(1,1) = v0_35inv(1,1) + 1/sigmasq5;
        mat postvar35 = inv(v0_35inv); mat meantem35(2,1);
        meantem35 = meantem.rows(ind35); mat postmean35 = postvar35 * meantem35;
        
        mat v0_45 = v0.submat(ind45,ind45); mat v0_45inv = inv(v0_45);
        v0_45inv(0,0) = v0_45inv(0,0) + 1/sigmasq4; v0_45inv(1,1) = v0_45inv(1,1) + 1/sigmasq5;
        mat postvar45 = inv(v0_45inv); mat meantem45(2,1);
        meantem45 = meantem.rows(ind45); mat postmean45 = postvar45 * meantem45;
        
        mat v0_123 = v0.submat(ind123,ind123); mat v0_123inv = inv(v0_123);
        v0_123inv(0,0) = v0_123inv(0,0) + 1/sigmasq1; v0_123inv(1,1) = v0_123inv(1,1) + 1/sigmasq2; v0_123inv(2,2) = v0_123inv(2,2) + 1/sigmasq3;
        mat postvar123 = inv(v0_123inv); mat meantem123(3,1);
        meantem123 = meantem.rows(ind123); mat postmean123 = postvar123 * meantem123;
        
        mat v0_124 = v0.submat(ind124,ind124); mat v0_124inv = inv(v0_124);
        v0_124inv(0,0) = v0_124inv(0,0) + 1/sigmasq1; v0_124inv(1,1) = v0_124inv(1,1) + 1/sigmasq2; v0_124inv(2,2) = v0_124inv(2,2) + 1/sigmasq4;
        mat postvar124 = inv(v0_124inv); mat meantem124(3,1);
        meantem124 = meantem.rows(ind124); mat postmean124 = postvar124 * meantem124;
        
        mat v0_125 = v0.submat(ind125,ind125); mat v0_125inv = inv(v0_125);
        v0_125inv(0,0) = v0_125inv(0,0) + 1/sigmasq1; v0_125inv(1,1) = v0_125inv(1,1) + 1/sigmasq2; v0_125inv(2,2) = v0_125inv(2,2) + 1/sigmasq5;
        mat postvar125 = inv(v0_125inv); mat meantem125(3,1);
        meantem125 = meantem.rows(ind125); mat postmean125 = postvar125 * meantem125;
        
        mat v0_134 = v0.submat(ind134,ind134); mat v0_134inv = inv(v0_134);
        v0_134inv(0,0) = v0_134inv(0,0) + 1/sigmasq1; v0_134inv(1,1) = v0_134inv(1,1) + 1/sigmasq3; v0_134inv(2,2) = v0_134inv(2,2) + 1/sigmasq4;
        mat postvar134 = inv(v0_134inv); mat meantem134(3,1);
        meantem134 = meantem.rows(ind134); mat postmean134 = postvar134 * meantem134;
        
        mat v0_135 = v0.submat(ind135,ind135); mat v0_135inv = inv(v0_135);
        v0_135inv(0,0) = v0_135inv(0,0) + 1/sigmasq1; v0_135inv(1,1) = v0_135inv(1,1) + 1/sigmasq3; v0_135inv(2,2) = v0_135inv(2,2) + 1/sigmasq5;
        mat postvar135 = inv(v0_135inv); mat meantem135(3,1);
        meantem135 = meantem.rows(ind135); mat postmean135 = postvar135 * meantem135;
        
        mat v0_145 = v0.submat(ind145,ind145); mat v0_145inv = inv(v0_145);
        v0_145inv(0,0) = v0_145inv(0,0) + 1/sigmasq1; v0_145inv(1,1) = v0_145inv(1,1) + 1/sigmasq4; v0_145inv(2,2) = v0_145inv(2,2) + 1/sigmasq5;
        mat postvar145 = inv(v0_145inv); mat meantem145(3,1);
        meantem145 = meantem.rows(ind145); mat postmean145 = postvar145 * meantem145;
        
        mat v0_234 = v0.submat(ind234,ind234); mat v0_234inv = inv(v0_234);
        v0_234inv(0,0) = v0_234inv(0,0) + 1/sigmasq2; v0_234inv(1,1) = v0_234inv(1,1) + 1/sigmasq3; v0_234inv(2,2) = v0_234inv(2,2) + 1/sigmasq4;
        mat postvar234 = inv(v0_234inv); mat meantem234(3,1);
        meantem234 = meantem.rows(ind234); mat postmean234 = postvar234 * meantem234;
        
        mat v0_235 = v0.submat(ind235,ind235); mat v0_235inv = inv(v0_235);
        v0_235inv(0,0) = v0_235inv(0,0) + 1/sigmasq2; v0_235inv(1,1) = v0_235inv(1,1) + 1/sigmasq3; v0_235inv(2,2) = v0_235inv(2,2) + 1/sigmasq5;
        mat postvar235 = inv(v0_235inv); mat meantem235(3,1);
        meantem235 = meantem.rows(ind235); mat postmean235 = postvar235 * meantem235;
        
        mat v0_245 = v0.submat(ind245,ind245); mat v0_245inv = inv(v0_245);
        v0_245inv(0,0) = v0_245inv(0,0) + 1/sigmasq2; v0_245inv(1,1) = v0_245inv(1,1) + 1/sigmasq4; v0_245inv(2,2) = v0_245inv(2,2) + 1/sigmasq5;
        mat postvar245 = inv(v0_245inv); mat meantem245(3,1);
        meantem245 = meantem.rows(ind245); mat postmean245 = postvar245 * meantem245;
        
        mat v0_345 = v0.submat(ind345,ind345); mat v0_345inv = inv(v0_345);
        v0_345inv(0,0) = v0_345inv(0,0) + 1/sigmasq3; v0_345inv(1,1) = v0_345inv(1,1) + 1/sigmasq4; v0_345inv(2,2) = v0_345inv(2,2) + 1/sigmasq5;
        mat postvar345 = inv(v0_345inv); mat meantem345(3,1);
        meantem345 = meantem.rows(ind345); mat postmean345 = postvar345 * meantem345;
        
        mat v0_1234 = v0.submat(ind1234,ind1234); mat v0_1234inv = inv(v0_1234);
        v0_1234inv(0,0) = v0_1234inv(0,0) + 1/sigmasq1; v0_1234inv(1,1) = v0_1234inv(1,1) + 1/sigmasq2; v0_1234inv(2,2) = v0_1234inv(2,2) + 1/sigmasq3; v0_1234inv(3,3) = v0_1234inv(3,3) + 1/sigmasq4;
        mat postvar1234 = inv(v0_1234inv); mat meantem1234(4,1);
        meantem1234 = meantem.rows(ind1234); mat postmean1234 = postvar1234 * meantem1234;
        
        mat v0_1235 = v0.submat(ind1235,ind1235); mat v0_1235inv = inv(v0_1235);
        v0_1235inv(0,0) = v0_1235inv(0,0) + 1/sigmasq1; v0_1235inv(1,1) = v0_1235inv(1,1) + 1/sigmasq2; v0_1235inv(2,2) = v0_1235inv(2,2) + 1/sigmasq3; v0_1235inv(3,3) = v0_1235inv(3,3) + 1/sigmasq5;
        mat postvar1235 = inv(v0_1235inv); mat meantem1235(4,1);
        meantem1235 = meantem.rows(ind1235); mat postmean1235 = postvar1235 * meantem1235;
        
        mat v0_1245 = v0.submat(ind1245,ind1245); mat v0_1245inv = inv(v0_1245);
        v0_1245inv(0,0) = v0_1245inv(0,0) + 1/sigmasq1; v0_1245inv(1,1) = v0_1245inv(1,1) + 1/sigmasq2; v0_1245inv(2,2) = v0_1245inv(2,2) + 1/sigmasq4; v0_1245inv(3,3) = v0_1245inv(3,3) + 1/sigmasq5;
        mat postvar1245 = inv(v0_1245inv); mat meantem1245(4,1);
        meantem1245 = meantem.rows(ind1245); mat postmean1245 = postvar1245 * meantem1245;
        
        mat v0_1345 = v0.submat(ind1345,ind1345); mat v0_1345inv = inv(v0_1345);
        v0_1345inv(0,0) = v0_1345inv(0,0) + 1/sigmasq1; v0_1345inv(1,1) = v0_1345inv(1,1) + 1/sigmasq3; v0_1345inv(2,2) = v0_1345inv(2,2) + 1/sigmasq4; v0_1345inv(3,3) = v0_1345inv(3,3) + 1/sigmasq5;
        mat postvar1345 = inv(v0_1345inv); mat meantem1345(4,1);
        meantem1345 = meantem.rows(ind1345); mat postmean1345 = postvar1345 * meantem1345;
        
        mat v0_2345 = v0.submat(ind2345,ind2345); mat v0_2345inv = inv(v0_2345);
        v0_2345inv(0,0) = v0_2345inv(0,0) + 1/sigmasq2; v0_2345inv(1,1) = v0_2345inv(1,1) + 1/sigmasq3; v0_2345inv(2,2) = v0_2345inv(2,2) + 1/sigmasq4; v0_2345inv(3,3) = v0_2345inv(3,3) + 1/sigmasq5;
        mat postvar2345 = inv(v0_2345inv); mat meantem2345(4,1);
        meantem2345 = meantem.rows(ind2345); mat postmean2345 = postvar2345 * meantem2345;
        
        double postvar1 = 1/(1/sigmasq1 + 1/hsq1); double postmean1 = ((bh1 - mu1)/sigmasq1) * postvar1;
        double postvar2 = 1/(1/sigmasq2 + 1/hsq2); double postmean2 = ((bh2 - mu2)/sigmasq2) * postvar2;
        double postvar3 = 1/(1/sigmasq3 + 1/hsq3); double postmean3 = ((bh3 - mu3)/sigmasq3) * postvar3;
        double postvar4 = 1/(1/sigmasq4 + 1/hsq4); double postmean4 = ((bh4 - mu4)/sigmasq4) * postvar4;
        double postvar5 = 1/(1/sigmasq5 + 1/hsq5); double postmean5 = ((bh5 - mu5)/sigmasq5) * postvar5;
        
        NumericVector g = runif(1); double gamma = g(0);
        
        // 1 1 1 1 1
        if (gamma < fs(0)){
          bjnew.col(i) = rmvncpp(1, postmean, postvar);
        }
        uvec indi(1); indi(0) = i;
        // 1 1 1 1 0
        if ((gamma >= fs(0))&(gamma < (fs(0)+fs(1)))){
          bjnew.submat(ind1234,indi) = rmvncpp(1, postmean1234, postvar1234);
        }
        // 1 1 1 0 1
        if ((gamma >= (fs(0)+fs(1)))&(gamma < (fs(0)+fs(1)+fs(2)))){
          bjnew.submat(ind1235,indi) = rmvncpp(1, postmean1235, postvar1235);
        }
        // 1 1 0 1 1
        if ((gamma >= (fs(0)+fs(1)+fs(2)))&(gamma < (fs(0)+fs(1)+fs(2)+fs(3)))){
          bjnew.submat(ind1245,indi) = rmvncpp(1, postmean1245, postvar1245);
        }
        // 1 0 1 1 1
        if ((gamma >= (fs(0)+fs(1)+fs(2)+fs(3)))&(gamma < (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)))){
          bjnew.submat(ind1345,indi) = rmvncpp(1, postmean1345, postvar1345);
        }
        // 0 1 1 1 1
        if ((gamma >= (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)))&(gamma < (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)))){
          bjnew.submat(ind2345,indi) = rmvncpp(1, postmean2345, postvar2345);
        }
        // 1 1 1 0 0
        if ((gamma >= (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)))&(gamma < (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)))){
          bjnew.submat(ind123,indi) = rmvncpp(1, postmean123, postvar123);
        }
        // 1 1 0 1 0
        if ((gamma >= (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)))&(gamma < (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)))){
          bjnew.submat(ind124,indi) = rmvncpp(1, postmean124, postvar124);
        }
        // 1 1 0 0 1
        if ((gamma >= (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)))&(gamma < (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)+fs(8)))){
          bjnew.submat(ind125,indi) = rmvncpp(1, postmean125, postvar125);
        }
        // 1 0 1 1 0
        if ((gamma >= (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)+fs(8)))&(gamma < (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)+fs(8)+fs(9)))){
          bjnew.submat(ind134,indi) = rmvncpp(1, postmean134, postvar134);
        }
        // 1 0 1 0 1
        if ((gamma >= (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)+fs(8)+fs(9)))&(gamma < (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)+fs(8)+fs(9)+fs(10)))){
          bjnew.submat(ind135,indi) = rmvncpp(1, postmean135, postvar135);
        }
        // 1 0 0 1 1
        if ((gamma >= (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)+fs(8)+fs(9)+fs(10)))&(gamma < (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)+fs(8)+fs(9)+fs(10)+fs(11)))){
          bjnew.submat(ind145,indi) = rmvncpp(1, postmean145, postvar145);
        }
        // 0 1 1 1 0
        if ((gamma >= (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)+fs(8)+fs(9)+fs(10)+fs(11)))&(gamma < (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)+fs(8)+fs(9)+fs(10)+fs(11)+fs(12)))){
          bjnew.submat(ind234,indi) = rmvncpp(1, postmean234, postvar234);
        }
        // 0 1 1 0 1
        if ((gamma >= (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)+fs(8)+fs(9)+fs(10)+fs(11)+fs(12)))&(gamma < (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)+fs(8)+fs(9)+fs(10)+fs(11)+fs(12)+fs(13)))){
          bjnew.submat(ind235,indi) = rmvncpp(1, postmean235, postvar235);
        }
        // 0 1 0 1 1
        if ((gamma >= (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)+fs(8)+fs(9)+fs(10)+fs(11)+fs(12)+fs(13)))&(gamma < (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)+fs(8)+fs(9)+fs(10)+fs(11)+fs(12)+fs(13)+fs(14)))){
          bjnew.submat(ind245,indi) = rmvncpp(1, postmean245, postvar245);
        }
        // 0 0 1 1 1
        if ((gamma >= (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)+fs(8)+fs(9)+fs(10)+fs(11)+fs(12)+fs(13)+fs(14)))&(gamma < (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)+fs(8)+fs(9)+fs(10)+fs(11)+fs(12)+fs(13)+fs(14)+fs(15)))){
          bjnew.submat(ind345,indi) = rmvncpp(1, postmean345, postvar345);
        }
        // 1 1 0 0 0
        if ((gamma >= (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)+fs(8)+fs(9)+fs(10)+fs(11)+fs(12)+fs(13)+fs(14)+fs(15)))&(gamma < (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)+fs(8)+fs(9)+fs(10)+fs(11)+fs(12)+fs(13)+fs(14)+fs(15)+fs(16)))){
          bjnew.submat(ind12,indi) = rmvncpp(1, postmean12, postvar12);
        }
        // 1 0 1 0 0
        if ((gamma >= (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)+fs(8)+fs(9)+fs(10)+fs(11)+fs(12)+fs(13)+fs(14)+fs(15)+fs(16)))&(gamma < (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)+fs(8)+fs(9)+fs(10)+fs(11)+fs(12)+fs(13)+fs(14)+fs(15)+fs(16)+fs(17)))){
          bjnew.submat(ind13,indi) = rmvncpp(1, postmean13, postvar13);
        }
        // 1 0 0 1 0
        if ((gamma >= (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)+fs(8)+fs(9)+fs(10)+fs(11)+fs(12)+fs(13)+fs(14)+fs(15)+fs(16)+fs(17)))&(gamma < (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)+fs(8)+fs(9)+fs(10)+fs(11)+fs(12)+fs(13)+fs(14)+fs(15)+fs(16)+fs(17)+fs(18)))){
          bjnew.submat(ind14,indi) = rmvncpp(1, postmean14, postvar14);
        }
        // 1 0 0 0 1
        if ((gamma >= (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)+fs(8)+fs(9)+fs(10)+fs(11)+fs(12)+fs(13)+fs(14)+fs(15)+fs(16)+fs(17)+fs(18)))&(gamma < (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)+fs(8)+fs(9)+fs(10)+fs(11)+fs(12)+fs(13)+fs(14)+fs(15)+fs(16)+fs(17)+fs(18)+fs(19)))){
          bjnew.submat(ind15,indi) = rmvncpp(1, postmean15, postvar15);
        }
        // 0 1 1 0 0
        if ((gamma >= (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)+fs(8)+fs(9)+fs(10)+fs(11)+fs(12)+fs(13)+fs(14)+fs(15)+fs(16)+fs(17)+fs(18)+fs(19)))&(gamma < (fs(0)+fs(1)+fs(2)+fs(3)+fs(4)+fs(5)+fs(6)+fs(7)+fs(8)+fs(9)+fs(10)+fs(11)+fs(12)+fs(13)+fs(14)+fs(15)+fs(16)+fs(17)+fs(18)+fs(19)+fs(20)))){
          bjnew.submat(ind23,indi) = rmvncpp(1, postmean23, postvar23);
        }
        // 0 1 0 1 0
        if ((gamma >= (fs_0_20))&(gamma < (fs_0_20+fs(21)))){
          bjnew.submat(ind24,indi) = rmvncpp(1, postmean24, postvar24);
        }
        // 0 1 0 0 1
        if ((gamma >= (fs_0_20+fs(21)))&(gamma < (fs_0_20+fs(21)+fs(22)))){
          bjnew.submat(ind25,indi) = rmvncpp(1, postmean25, postvar25);
        }
        // 0 0 1 1 0
        if ((gamma >= (fs_0_20+fs(21)+fs(22)))&(gamma < (fs_0_20+fs(21)+fs(22)+fs(23)))){
          bjnew.submat(ind34,indi) = rmvncpp(1, postmean34, postvar34);
        }
        // 0 0 1 0 1
        if ((gamma >= (fs_0_20+fs(21)+fs(22)+fs(23)))&(gamma < (fs_0_20+fs(21)+fs(22)+fs(23)+fs(24)))){
          bjnew.submat(ind35,indi) = rmvncpp(1, postmean35, postvar35);
        }
        // 0 0 0 1 1
        if ((gamma >= (fs_0_20+fs(21)+fs(22)+fs(23)+fs(24)))&(gamma < (fs_0_20+fs(21)+fs(22)+fs(23)+fs(24)+fs(25)))){
          bjnew.submat(ind45,indi) = rmvncpp(1, postmean45, postvar45);
        }
        // 1 0 0 0 0
        if ((gamma >= (fs_0_20+fs(21)+fs(22)+fs(23)+fs(24)+fs(25)))&(gamma < (fs_0_20+fs(21)+fs(22)+fs(23)+fs(24)+fs(25)+fs(26)))){
          NumericVector tem = rnorm(1, postmean1, sqrt(postvar1));
          bjnew(0,i) = tem(0);
        }
        // 0 1 0 0 0
        if ((gamma >= (fs_0_20+fs(21)+fs(22)+fs(23)+fs(24)+fs(25)+fs(26)))&(gamma < (fs_0_20+fs(21)+fs(22)+fs(23)+fs(24)+fs(25)+fs(26)+fs(27)))){
          NumericVector tem = rnorm(1, postmean2, sqrt(postvar2));
          bjnew(1,i) = tem(0);
        }
        // 0 0 1 0 0
        if ((gamma >= (fs_0_20+fs(21)+fs(22)+fs(23)+fs(24)+fs(25)+fs(26)+fs(27)))&(gamma < (fs_0_20+fs(21)+fs(22)+fs(23)+fs(24)+fs(25)+fs(26)+fs(27)+fs(28)))){
          NumericVector tem = rnorm(1, postmean3, sqrt(postvar3));
          bjnew(2,i) = tem(0);
        }
        // 0 0 0 1 0
        if ((gamma >= (fs_0_20+fs(21)+fs(22)+fs(23)+fs(24)+fs(25)+fs(26)+fs(27)+fs(28)))&(gamma < (fs_0_20+fs(21)+fs(22)+fs(23)+fs(24)+fs(25)+fs(26)+fs(27)+fs(28)+fs(29)))){
          NumericVector tem = rnorm(1, postmean4, sqrt(postvar4));
          bjnew(3,i) = tem(0);
        }
        // 0 0 0 0 1
        if ((gamma >= (fs_0_20+fs(21)+fs(22)+fs(23)+fs(24)+fs(25)+fs(26)+fs(27)+fs(28)+fs(29)))&(gamma < (fs_0_20+fs(21)+fs(22)+fs(23)+fs(24)+fs(25)+fs(26)+fs(27)+fs(28)+fs(29)+fs(30)))){
          NumericVector tem = rnorm(1, postmean5, sqrt(postvar5));
          bjnew(4,i) = tem(0);
        }
        double pm1tem = fs(0)*postmean(0,0) + fs(1)*postmean1234(0,0) + fs(2)*postmean1235(0,0) + fs(3)*postmean1245(0,0) +
          fs(4)*postmean1345(0,0) + fs(6)*postmean123(0,0) + fs(7)*postmean124(0,0) + fs(8)*postmean125(0,0) +
          fs(9)*postmean134(0,0) +fs(10)*postmean135(0,0) + fs(11)*postmean145(0,0) +
          fs(16)*postmean12(0,0) + fs(17)*postmean13(0,0) + fs(18)*postmean14(0,0) + fs(19)*postmean15(0,0) +
          fs(26)*postmean1;
        double pm2tem = fs(0)*postmean(1,0) + fs(1)*postmean1234(1,0) + fs(2)*postmean1235(1,0) + fs(3)*postmean1245(1,0) +
          fs(5)*postmean2345(0,0) + fs(6)*postmean123(1,0) + fs(7)*postmean124(1,0) + fs(8)*postmean125(1,0) +
          fs(12)*postmean234(0,0) +fs(13)*postmean235(0,0) + fs(14)*postmean245(0,0) +
          fs(16)*postmean12(1,0) + fs(20)*postmean23(0,0) + fs(21)*postmean24(0,0) + fs(22)*postmean25(0,0) +
          fs(27)*postmean2;
        double pm3tem = fs(0)*postmean(2,0) + fs(1)*postmean1234(2,0) + fs(2)*postmean1235(2,0) + fs(4)*postmean1345(1,0) +
          fs(5)*postmean2345(1,0) + fs(6)*postmean123(2,0) + fs(9)*postmean134(1,0) + fs(10)*postmean135(1,0) +
          fs(12)*postmean234(1,0) +fs(13)*postmean235(1,0) + fs(15)*postmean345(0,0) +
          fs(17)*postmean13(1,0) + fs(20)*postmean23(1,0) + fs(23)*postmean34(0,0) + fs(24)*postmean35(0,0) +
          fs(28)*postmean3;
        double pm4tem = fs(0)*postmean(3,0) + fs(1)*postmean1234(3,0) + fs(3)*postmean1245(2,0) + fs(4)*postmean1345(2,0) +
          fs(5)*postmean2345(2,0) + fs(7)*postmean124(2,0) + fs(9)*postmean134(2,0) + fs(11)*postmean145(1,0) +
          fs(12)*postmean234(2,0) +fs(14)*postmean245(1,0) + fs(15)*postmean345(1,0) +
          fs(18)*postmean14(1,0) + fs(21)*postmean24(1,0) + fs(23)*postmean34(1,0) + fs(25)*postmean45(0,0) +
          fs(29)*postmean4;
        double pm5tem = fs(0)*postmean(4,0) + fs(2)*postmean1235(3,0) + fs(3)*postmean1245(3,0) + fs(4)*postmean1345(3,0) +
          fs(5)*postmean2345(3,0) + fs(8)*postmean125(2,0) + fs(10)*postmean135(2,0) + fs(11)*postmean145(2,0) +
          fs(13)*postmean235(2,0) +fs(14)*postmean245(2,0) + fs(15)*postmean345(2,0) +
          fs(19)*postmean15(1,0) + fs(22)*postmean25(1,0) + fs(24)*postmean35(1,0) + fs(25)*postmean45(1,0) +
          fs(30)*postmean5;
        pm(0,i) = pm1tem; pm(1,i) = pm2tem; pm(2,i) = pm3tem; pm(3,i) = pm4tem; pm(4,i) = pm5tem;
        // sparse
        if (sparse(0) == 1){
          double fs1 = fs(0) + fs(1) + fs(2) + fs(3) + fs(4) + 
            fs(6) + fs(7) + fs(8) + fs(9) + fs(10) + fs(11) +  
            fs(16) + fs(17) + fs(18) + fs(19) + fs(26);
          double fs1_0 = p11111 + p11110 + p11101 + p11011 + p10111 + 
            p11100 + p11010 + p11001 + p10110 + p10101 + p10011 + 
            p11000 + p10100 + p10010 + p10001 + p10000;
          if (fs1 < fs1_0){
            bjnew(0,i) = 0;
            pm(0,i) = 0;
            postmean(0,0) = 0;
            postmean1 = 0;
          }
        }
        if (sparse(1) == 1){
          double fs2 = fs(0) + fs(1) + fs(2) + fs(3) + fs(5) + 
            fs(6) + fs(7) + fs(8) + fs(12) + fs(13) + fs(14) +  
            fs(16) + fs(20) + fs(21) + fs(22) + fs(27);
          double fs2_0 = p11111 + p11110 + p11101 + p11011 + p01111 + 
            p11100 + p11010 + p11001 + p01110 + p01101 + p01011 + 
            p11000 + p01100 + p01010 + p01001 + p01000;
          if (fs2 < fs2_0){
            bjnew(1,i) = 0;
            pm(1,i) = 0;
            postmean(1,0) = 0;
            postmean2 = 0;
          }
        }
        if (sparse(2) == 1){
          double fs3 = fs(0) + fs(1) + fs(2) + fs(4) + fs(5) + 
            fs(6) + fs(9) + fs(10) + fs(12) + fs(13) + fs(15) +  
            fs(17) + fs(20) + fs(23) + fs(24) + fs(28);
          double fs3_0 = p11111 + p11110 + p11101 + p10111 + p01111 + 
            p11100 + p10110 + p10101 + p01110 + p01101 + p00111 + 
            p10100 + p01100 + p00110 + p00101 + p00100;
          if (fs3 < fs3_0){
            bjnew(2,i) = 0;
            pm(2,i) = 0;
            postmean(2,0) = 0;
            postmean3 = 0;
          }
        }
        if (sparse(3) == 1){
          double fs4 = fs(0) + fs(1) + fs(3) + fs(4) + fs(5) + 
            fs(7) + fs(9) + fs(11) + fs(12) + fs(14) + fs(15) +  
            fs(18) + fs(21) + fs(23) + fs(25) + fs(29);
          double fs4_0 = p11111 + p11110 + p11011 + p10111 + p01111 + 
            p11010 + p10110 + p10011 + p01110 + p01011 + p00111 + 
            p10010 + p01010 + p00110 + p00011 + p00010;
          if (fs4 < fs4_0){
            bjnew(3,i) = 0;
            pm(3,i) = 0;
            postmean(3,0) = 0;
            postmean4 = 0;
          }
        }
        if (sparse(4) == 1){
          double fs5 = fs(0) + fs(2) + fs(3) + fs(4) + fs(5) + 
            fs(8) + fs(10) + fs(11) + fs(13) + fs(14) + fs(15) +  
            fs(19) + fs(22) + fs(24) + fs(25) + fs(30);
          double fs5_0 = p11111 + p11101 + p11011 + p10111 + p01111 + 
            p11001 + p10101 + p10011 + p01101 + p01011 + p00111 + 
            p10001 + p01001 + p00101 + p00011 + p00001;
          if (fs5 < fs5_0){
            bjnew(4,i) = 0;
            pm(4,i) = 0;
            postmean(4,0) = 0;
            postmean5 = 0;
          }
        }
    }
    Bjh1(x1) = bjnew(0,i);
    Bjh2(x2) = bjnew(1,i);
    Bjh3(x3) = bjnew(2,i);
    Bjh4(x4) = bjnew(3,i);
    Bjh5(x5) = bjnew(4,i);
  }
  //return List::create(_["ptilde"] = ptilde, _["bjnew"] = bjnew, _["pm"] = pm);
  return List::create(_["bjnew"] = bjnew, _["pm"] = pm); //_["ptilde"] = ptilde,
}




// Update rho

// [[Rcpp::export]]
List densitycpp(vec indx1, vec indx2, vec overlapping_indx, vec bjh10, vec bjh20,
                double hsq1, double hsq2, double rho)
{
  double len_overlap = overlapping_indx.n_elem;
  double len_nonoverlap1 = indx1.n_elem;
  double len_nonoverlap2 = indx2.n_elem;
  int n_snp = len_overlap+len_nonoverlap1+len_nonoverlap2;
  vec dens1(n_snp); // density use for hsq1
  vec dens2(n_snp); // density use for hsq2
  vec dens12(n_snp); // density use for both hsq1 and hsq2
  dens1.fill(0);
  dens2.fill(0);
  dens12.fill(0);
  // s-1 density
  // 1. population 1 nonoverlapping density
  for (int i = 1; i < len_nonoverlap1; i++){
    int x = indx1(i) - 1; // snp index
    double beta1 = bjh10(x);
    if (beta1 != 0){
      dens1(x) = R::dnorm(beta1,0,sqrt(hsq1),0);
      dens1(x) = log(dens1(x));
    }
  }
  // 2. population 2 nonoverlapping density
  for (int i = 1; i < len_nonoverlap2; i++){
    int x = indx2(i) - 1; // snp index
    // population 2
    double beta2 = bjh20(x);
    if (beta2 != 0){
      dens2(x) = R::dnorm(beta2,0,sqrt(hsq2),0);
      dens2(x) = log(dens2(x));
    }
  }
  // 3. population 1&2 overlapping density
  double v11 = hsq1;
  double v22 = hsq2;
  double v12 = sqrt(hsq1*hsq2) * rho;
  mat v(2,2);
  v(0,0) = v11;
  v(1,0) = v(0,1) = v12;
  v(1,1) = v22;
  double rho_density = 0;
  for (int i = 1; i < len_overlap; i++){
    int x = overlapping_indx(i) - 1; // snp index
    double beta1 = bjh10(x);
    double beta2 = bjh20(x);
    if ((beta1 != 0) & (beta2 != 0)){
      rowvec mu0(2);
      mu0.fill(0);
      mat betas(1,2);
      betas(0,0) = beta1;
      betas(0,1) = beta2;
      vec tem = dmvnormcpp(betas, mu0, v);
      dens12(x) = tem(0);
      // work with log-densities for computational reasons
      dens12(x) = log(dens12(x));
      rho_density = rho_density + dens12(x);
    }
    if ((beta1 == 0) & (beta2 != 0)){
      dens12(x) = R::dnorm(beta2,0,sqrt(hsq2),0);
      dens12(x) = log(dens12(x));
    }
    if ((beta1 != 0) & (beta2 == 0)){
      dens12(x) = R::dnorm(beta1,0,sqrt(hsq1),0);
      dens12(x) = log(dens12(x));
    }
  }
  mat dense(n_snp,3);
  dense.col(0) = dens1;
  dense.col(1) = dens2;
  dense.col(2) = dens12;
  return List::create(_["dense"] = dense,_["rho_density"] = rho_density);
}
//update covariance matrix of beta1j,beta2j:
  
  // [[Rcpp::export]]
List Supdate(vec overlap_indx, vec bjh10, vec bjh20, double nu, mat Omega0){
  double len_co = overlap_indx.n_elem;
  mat BBt(2,2);
  double df = nu;
  BBt.fill(0);
  for (int i = 1; i < len_co; i++){
    int x = overlap_indx(i) - 1; // snp index
    mat B(2,1);
    B(0,0) = bjh10(x);
    B(1,0) = bjh20(x);
    if ((B(0,0) != 0)&(B(1,0) != 0)){
      mat tem = B * B.t();
      BBt = BBt + tem;
      df = df + 1;
    }
  }
  mat S = Omega0 + BBt;
  return List::create(_["df"] = df, _["BBt"] = BBt);
}


// [[Rcpp::export]]
List h2update(vec indx, vec bjh1, vec bjh2, double shape, double scale){
  double len_indx = indx.n_elem;
  double alpha1 = shape;
  double beta1 = scale;
  double alpha2 = shape;
  double beta2 = scale;
  for (int i = 1; i < len_indx; i++){
    int x = indx(i) - 1; // snp index
    if ((bjh1(x) != 0)&(bjh2(x) == 0)){
      alpha1 = alpha1 + 0.5;
      beta1 = beta1 + bjh1(x)*bjh1(x)*0.5;
    }
    if ((bjh1(x) == 0)&(bjh2(x) != 0)){
      alpha2 = alpha2 + 0.5;
      beta2 = beta2 + bjh2(x)*bjh2(x)*0.5;
    }
  }
  return List::create(_["alpha1"] = alpha1, _["beta1"] = beta1, _["alpha2"] = alpha2, _["beta2"] = beta2);
}

