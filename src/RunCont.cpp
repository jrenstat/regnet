
#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include"ContCD.h"
#include"RobustCD.h"
#include"QR.h"
#include"RunCont.h"
#include"Utilities.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List RunCont(arma::mat const &xc, arma::mat const &xg, arma::vec const &y, double lamb1, double lamb2, arma::vec bc, arma::vec bg, double r, arma::mat const &a, arma::vec const &triRowAbsSums, int p, int pc, char method, int maxit, double tol)
{
  int count = 0, n = xc.n_rows;
  int niter = 0;
  arma::vec bold(p, fill::none), yc, yg, u; // bc = bc0, bg = bg0;
  arma::vec const inp = arma::sum(arma::square(xg),0).t()/n;
  bool converged = false;
  double diff = NA_REAL;
  if(method == 'n') u = inp + lamb2 * triRowAbsSums;

  while(count < maxit){
    Rcpp::checkUserInterrupt();
    niter++;
		yc = y - xg * bg;
		bc = fastLm(yc, xc);
		yg = y - xc * bc;
	bold = bg;
    if(method == 'n'){
      ContNet(xg, yg, lamb1, lamb2, bg, r, a, u, n, p);
    }else if(method == 'm'){
      ContMCP(xg, yg, lamb1, bg, r, inp, n, p);
    }else{
      ContLasso(xg, yg, lamb1, bg, inp, n, p);
	    }
	    diff = arma::accu(arma::abs(bg - bold))/(arma::accu(bg != 0)+0.1);
	    if(diff < tol){
	      converged = true;
	      break;
	    }
	    else{
	      count++;
	    }
	  }
	  arma::vec b1(pc+p, fill::none);
	  b1.subvec(0, pc-1) = bc;
	  b1.subvec(pc, pc+p-1) = bg;
  return Rcpp::List::create(Rcpp::Named("b") = b1,
                            Rcpp::Named("converged") = converged,
                            Rcpp::Named("niter") = niter,
                            Rcpp::Named("diff") = diff);
}


// [[Rcpp::export]]
Rcpp::List RunCont_robust(arma::mat const &xc, arma::mat const &xg, arma::vec const &y, double lamb1, double lamb2, arma::vec bc, arma::vec bg, double r, arma::mat const &a, int p, int pc, char method, bool debugging, int maxit, double tol)
{
  int count = 0, n = xc.n_rows;
  int niter = 0;
  arma::vec bold(p, fill::none), yc, yg; // bc = bc0, bg = bg0;
  arma::mat const wc = arma::abs(xc)/n;
  arma::vec const totalWeights = arma::sum(wc, 0).t();
  bool converged = false;
  double diff = NA_REAL;

  while(count < maxit){
    Rcpp::checkUserInterrupt();
    niter++;
		yc = y - xg * bg;
		// bc = QRWMR(xc, yc, bc);
		QRWMR(xc, yc, bc, wc, totalWeights, debugging);
	yg = y - xc * bc;
	bold = bg;
    if(method == 'n'){
      LadNet(xg, yg, lamb1, lamb2, bg, r, a, n, p, debugging);
    }else if(method == 'm'){
      LadMCP(xg, yg, lamb1, bg, r, n, p, debugging);
    }else{
      LadLasso(xg, yg, lamb1, bg, n, p, debugging);
	    }
	    diff = arma::accu(arma::abs(bg - bold))/(arma::accu(bg != 0)+0.1);
		//Rcpp::Rcout << "diff: " << diff <<std::endl;
	    if(diff < tol){
	      converged = true;
	      break;
	    }
	    else{
	      // bg = bnew;
	      count++;
	    }
	  }
	  arma::vec b1(pc+p, fill::none);
	  b1.subvec(0, pc-1) = bc;
	  b1.subvec(pc, pc+p-1) = bg;
  return Rcpp::List::create(Rcpp::Named("b") = b1,
                            Rcpp::Named("converged") = converged,
                            Rcpp::Named("niter") = niter,
                            Rcpp::Named("diff") = diff);
}






