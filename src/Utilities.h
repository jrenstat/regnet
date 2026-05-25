#ifndef UTILITIES_h
#define UTILITIES_h

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

double Soft(double z, double lambda);
arma::vec sigmoid_stable(arma::vec const &eta);
arma::vec TriRowAbsSums(arma::mat const &a);
double validation_LS(arma::mat  const &x, arma::vec  const &y, arma::vec  const &b);
double validation_LS(arma::mat const &xc, arma::mat const &xg, arma::vec const &y, arma::vec const &b, int p, int pc);
double validation_LAD(arma::mat  const &x, arma::vec  const &y, arma::vec  const &b);
double validation_LAD(arma::mat const &xc, arma::mat const &xg, arma::vec const &y, arma::vec const &bc, arma::vec const &bg);
double validation_logit(arma::mat const &x0, arma::vec const &y0, arma::vec const &b);
arma::vec validation_logit2(arma::mat const &x0, arma::vec const &y0, arma::vec const &b);
arma::vec fastLm(arma::vec const &y, arma::mat const &X);

#endif
