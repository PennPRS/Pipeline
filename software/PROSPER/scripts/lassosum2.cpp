#include <cstdlib>
#include <time.h>
#include <string>
#include <iostream>// std::cout
#include <cmath>
#include <algorithm>// std::min
#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;
using namespace std;
using namespace stats;

//--------------------------------
  //--------------------------------
  // Construct path

//' Construct lambda path for single ethnic group
//'
//' @param summ summary statistics
//' @param delta tuning parameter for L1/L2
//' @param nlambda number of lambda in the path
//' @param lambda_min_ratio the ratio of the maximum and minimum of lambdas
//' @return lambda path
//' @keywords internal
// [[Rcpp::export]]
vec l_path(double summmax, int nlambda=20, double lambda_min_ratio = 0.1){
  nlambda += 1;
  double tmp = abs(summmax);
  double loglambda_max = log(tmp);
  double loglambda_min = std::min(log(0.001), log( tmp * lambda_min_ratio));
  vec lambdapath(nlambda-1);
  for (int i=1; i < nlambda; i++) {
    lambdapath(i-1) = exp(loglambda_max - (loglambda_max-loglambda_min)/(nlambda-1)  * i);
  }
  return(lambdapath);
}

// [[Rcpp::export]]
vec d_path(double max=100, double min=0.5, int ndelta = 10){
  double sqrt_max =  pow(max, 1.0/3.0);
  double sqrt_min =  pow(min, 1.0/3.0);
  vec deltapath(ndelta);
  for (int i=0; i < ndelta; ++i) {
    deltapath[ndelta-1-i] = pow( sqrt_max - (sqrt_max-sqrt_min)/(ndelta-1)  * i, 3);
  }
  return(deltapath);
}


//--------------------------------
  //--------------------------------
  // Multi-ethic analysis


//' Compute beta joint
//'
//' @ NOTE: input data should be cleaned to standard format (SNP order should be same with reference SNP list)
//' @param summ: summary statistics  (set to 0 for ethnic-specific SNPs)
//' @param R: correlation matrix of all SNPs (set to 0 for ethnic-specific SNPs)
//' @param indx: indicator of whether a SNP in a certain ethnic group
//' @param M: number of ethnic groups
//' @param delta: best tuning parameter delta in single-ethnic lassosum2 (len=M)
//' @param lambda: lambda penalty (len=M)
//' @param c: correlation penalty (dim=M*M)
//' @return beta vector and computation details
//' @keywords internal
// [[Rcpp::export]]
List enet_singlethnic_bl(vec summ,
                         mat R,
                         double lambda,
                         double delta,
                         double thresh=1e-04, int maxiter=1000){

  int p = summ.n_elem; // number of SNP

  if(p==0){
    return List::create(Named("conv") = 0,
                        Named("niter") = 0,
                        Named("lambda") = lambda,
                        Named("delta") = delta,
                        Named("b") = summ);
  }
    
  vec b(p); b.fill(0.0); // store coefficient

  double dlx, del, tmp, ui;

  int conv=0;
  int niter;
  double denom=1.0+delta;

  for(int k=0; k<maxiter ;++k) {
    dlx=0.0;

    for (int i=0; i<p; ++i) {
      tmp = b(i);
      ui = summ(i) - (dot(R.col(i), b)-b(i));
      b(i) = sign(ui) * std::max(0.0, abs(ui)-lambda) / denom;
      if(b(i)==tmp){
        continue;
      }else{
        del=b(i)-tmp;
        dlx=std::max(dlx,abs(del));
      }
    }
    if(dlx < thresh) {
      conv=1;
      niter=k;
      break;
    }
  }

  if(conv==0){
    niter=maxiter;
    b.fill(0.0);
  }

  return List::create(Named("conv") = conv,
                      Named("niter") = niter,
                      Named("lambda") = lambda,
                      Named("delta") = delta,
                      Named("b") = b);
}

//' Repeatedly compute beta by blocks for fixed tuning parameters
//'
//' @ NOTE: input data should be cleaned to standard format (SNP order should be same with reference SNP list)
//' @param summ: list of summary statistics by blocks
//' @param R: LD correlation matrix by block
//' @param M: number of ethnic groups
//' @param indx: indicator of whether a SNP in a certain ethnic group
//' @param indx_block: indicator of whether an ethnic group has data in a certain block
//' @param delta: best tuning parameter delta in single-ethnic lassosum2 (len=M)
//' @param lambda: lambda penalty (len=M)
//' @param c: correlation penalty (dim=M*M)
//' @keywords internal
// [[Rcpp::export]]
List enet_singlethnic_fixtuning(List summ, List R,
                                double lambda, double delta,
                                double thresh=1e-04, int maxiter=1000){
  // List summ, List R

  int nblock = summ.size();

  vec conv(nblock);
  vec niter(nblock);
  List b(nblock);

  for(int bl = 0; bl < nblock; ++bl) {

    if(R[bl]==R_NilValue) {
      conv(bl) = NA_INTEGER;
      niter(bl) = NA_INTEGER;
      b[bl] =  R_NilValue;
    }else{
      vec summ_tmp = summ[bl];
      mat R_tmp = R[bl];
      List beta_bl = enet_singlethnic_bl(summ_tmp,R_tmp,
                                         lambda, delta, thresh, maxiter);
      conv(bl) = beta_bl["conv"];
      niter(bl) = beta_bl["niter"];
      b[bl] =  beta_bl["b"];
    }
  }

  return List::create(Named("conv") = conv,
                      Named("niter") = niter,
                      Named("lambda") = lambda,
                      Named("delta") = delta,
                      Named("b") = b);

}


//' Repeatedly compute beta by blocks for all tuning parameters
//'
//' @ NOTE: input data should be cleaned to standard format (SNP order should be same with reference SNP list)
//' @param summ: list of summary statistics by block
//' @param R: LD correlation matrix by block
//' @param M: number of ethnic groups
//' @param indx: indicator of whether a SNP in a certain ethnic group
//' @param indx_block: indicator of whether an ethnic group has data in a certain block
//' @param delta: best tuning parameter delta in single-ethnic lassosum2 (len=M)
//' @param lambdapath: lambda penalty (dim=M*Ntuning)
//' @param cpath: list of correlation penalty (size=Ntuning)
//' @keywords internal
// [[Rcpp::export]]
List enet_singlethnic(List summ, List R,
                      vec deltapath, vec lambdapath, 
                      double thresh=1e-04, int maxiter=1000, int verbose=2){

  double li, di;
  int Ll = lambdapath.n_elem;
  int Ld = deltapath.n_elem;
  int nblock = summ.size();
  int alltuning = Ll*Ld;

  mat conv(nblock, alltuning);
  mat niter(nblock, alltuning);
  List b(alltuning);
  vec lambda(alltuning);
  vec delta(alltuning);

//  if (verbose == 2){
//    Rcout << "* Starting model fitting --" << std::endl;
//  }

  int ii = 0;

  for (int ll = 0; ll < Ll; ++ll) {

    li = lambdapath(ll);

    for (int ld = 0; ld < Ld; ++ld) {
      di = deltapath(ld);
      List res_fixtuning = enet_singlethnic_fixtuning(summ, R, li, di, thresh, maxiter);
      lambda(ii) = li;
      delta(ii) = di;
      
      conv.col(ii) = as<vec>(res_fixtuning["conv"]);
      niter.col(ii) = as<vec>(res_fixtuning["niter"]);
      b[ii] = res_fixtuning["b"];

      if (verbose == 2){
        if (ii % 10 == 0){
          Rcout << " " << ii << ".." ;
        }
      }
      ii += 1;
    }
  }

  if (verbose == 2){
    Rcout << " Model fitting has been completed" << std::endl;
  }
  
  return List::create(Named("conv") = conv,
                      Named("niter") = niter,
                      Named("b") = b,
                      Named("lambda") = lambda,
                      Named("delta") = delta);
  
}

