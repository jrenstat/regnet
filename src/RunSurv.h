#ifndef RUNSURV_h
#define RUNSURV_h

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

Rcpp::List RunSurv_robust(arma::mat const &xc, arma::mat const &xg, arma::vec const &y, double lamb1, double lamb2, arma::vec bc, arma::vec bg, double r, arma::mat const &a, int p, int pc, char method, bool debugging, int maxit, double tol);

void RunSurv_robust_warm(arma::mat const &xc, arma::mat const &xg, arma::vec const &y, double lamb1, double lamb2, arma::vec &bc, arma::vec &bg, double r, arma::mat const &a, int p, int pc, char method, int maxit, double tol);

Rcpp::List RunSurv(arma::mat const &xc, arma::mat const &xg, arma::vec const &y, double lamb1, double lamb2, arma::vec bc, arma::vec bg, double r, arma::mat const &a, arma::vec const &triRowAbsSums, int p, int pc, char method, int maxit, double tol);

#endif
