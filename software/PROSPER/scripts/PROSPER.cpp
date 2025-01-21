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
vec r_path(double rmax, double rmin, int nr=20){
  double tmp;
  tmp = abs(rmax);
  double logr_max = log(tmp);
  tmp = abs(rmin);
  double logr_min = log(tmp);
  vec rpath(nr);
  for (int i=0; i < nr; i++) {
    rpath(i) = exp(logr_max - (logr_max-logr_min)/(nr-1)  * i);
  }
  return(rpath);
}

//' Construct c path for multi ethnic group
//'
//' @param lambda_max1 maximum lambda1
//' @param lambda_max2 maximum lambda2
//' @param nlambda number of lambda3 in the path
//' @return lambda path for lambda3 (smalles is 0)
//' @keywords internal
//'
// [[Rcpp::export]]
vec c_path(double maxc=100, double minc=0.5, int nc=10){
  double sqrtc_max = pow(maxc, 1.0/4);
  double sqrtc_min = pow(minc, 1.0/4);
  vec cpath(nc);
  for (int i=0; i < nc; i++) {
    cpath(i) = pow(sqrtc_max - (sqrtc_max-sqrtc_min)/(nc-1)  * i, 4.0);
  }
  return(cpath);
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
List enet_multiethnic_bl(mat summ, List R, mat indx, int M,
                         vec delta, vec lambda, mat c,
                         double thresh=1e-04, int maxiter=1000){
  int p = indx.n_rows; // number of SNP
  
  if(p==0){
    return List::create(Named("conv") = 0,
                        Named("niter") = 0,
                        Named("b") = summ);
  }

  mat b(p, M); b.fill(0.0); // store coefficient

  double dlx, del, tmp, ui, denom;

  int conv=0;
  int niter;
  for(int k=0; k<maxiter ;++k) {
    dlx=0.0;

    for(int l=0; l<M; ++l){
      vec summl = summ.col(l);
      mat Rl = as<mat>(R[l]);

      for (int i=0; i<p; ++i) {

        tmp = b(i,l);
        if(indx(i,l)!=0){
          ui = summl(i) - (dot(Rl.col(i), b.col(l))-b(i,l));
          denom=1.0+delta(l);
          for(int ll=0; ll<M; ++ll){
            if( (ll!=l) & (indx(i,ll)!=0) ){
              ui += c(ll,l)*b(i,ll);
              denom += c(ll,l);
            }
          }
          b(i,l) = sign(ui) * std::max(0.0, abs(ui)-lambda(l)) / denom;
        }else{
          continue;
        }

        if(b(i,l)==tmp){
          continue;
        }else{
          del=b(i,l)-tmp;
          dlx=std::max(dlx,abs(del));
        }
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
List enet_multiethnic_fixtuning(List summ, List R,
                                int M, List indx,
                                vec indx_block,
                                vec delta, vec lambda, mat c,
                                double thresh=1e-04, int maxiter=1000){
  // mat summ, List R, mat indx

  int nblock = indx_block.n_elem;

  vec conv(nblock);
  vec niter(nblock);
  List b(nblock);

  for(int bl = 0; bl < nblock; ++bl) {

    if(indx_block(bl) == 0){
      conv(bl) = NA_INTEGER;
      niter(bl) = NA_INTEGER;
      b[bl] =  R_NilValue;
    }else{
      mat summ_tmp = as<mat>(summ[bl]);
      List R_tmp = as<List>(R[bl]);
      mat indx_tmp = as<mat>(indx[bl]);

      List beta_bl = enet_multiethnic_bl(summ_tmp,R_tmp,
                                         indx_tmp, M,
                                         delta, lambda, c,
                                         thresh, maxiter);

      conv(bl) = beta_bl["conv"];
      niter(bl) = beta_bl["niter"];
      b[bl] = as<mat>(beta_bl["b"]);

    }

  }

  return List::create(Named("conv") = conv,
                      Named("niter") = niter,
                      Named("lambda") = lambda,
                      Named("c") = c,
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
List enet_multiethnic(List summ, List R,
                      int M, List indx,
                      vec indx_block,
                      vec delta, mat lambdapath, List cpath,
                      double thresh=1e-04, int maxiter=1000, int verbose=2){

  int Ll = lambdapath.n_cols;
  int Lc = cpath.size();

  int nblock = indx_block.n_elem;

  int alltuning = Ll*Lc;

  mat conv(nblock, alltuning);
  mat niter(nblock, alltuning);
  List b(alltuning);
  mat lambda(M,alltuning);
  List c(alltuning);

  if (verbose == 2){
    Rcout << "* Starting model fitting --" << std::endl;
  }

  int ii = 0;

  for (int ll = 0; ll < Ll; ++ll) {

    vec li = lambdapath.col(ll);

    for (int lc = 0; lc < Lc; ++lc) {
      mat ci = as<mat>(cpath[lc]);

      List res_fixtuning = enet_multiethnic_fixtuning(summ, R,
                                                      M, indx,
                                                      indx_block,
                                                      delta, li, ci,
                                                      thresh, maxiter);
      
      lambda.col(ii) = li;
      c[ii] = ci;
      conv.col(ii) = as<vec>(res_fixtuning["conv"]);
      niter.col(ii) = as<vec>(res_fixtuning["niter"]);
      b[ii] = res_fixtuning["b"];
      if (verbose == 2){
        if (ii % 10 == 0){
          Rcout << "Fitting " << ii << "/" << alltuning << " "  << std::endl;
        }
      }
      ii += 1;
      
    }
  }
  if (verbose == 2){
    Rcout << "* Model fitting for all tuning parameters has been completed!" << std::endl;
  }
  
  return List::create(Named("conv") = conv,
                      Named("niter") = niter,
                      Named("b") = b,
                      Named("lambda") = lambda,
                      Named("c") = c);
  
}



