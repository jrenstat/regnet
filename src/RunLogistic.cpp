#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include"LogitCD.h"
#include"RunLogistic.h"

using namespace Rcpp;
using namespace arma;

arma::vec RunLogit_fit(arma::mat const &x, arma::vec const &y, double lamb1, double lamb2, arma::vec b, double r, arma::mat const &a, arma::vec const &triRowAbsSums, int p, double alpha, char method, int maxit, double tol, bool &converged, int &niter, double &diff)
{
  int count = 0, n = x.n_rows;
  arma::vec bnew = b, u;
  converged = false;
  niter = 0;
  diff = NA_REAL;
  if(method == 'n') u = 0.25 + lamb2 * triRowAbsSums;
  
  while(count < maxit){
    Rcpp::checkUserInterrupt();
    niter++;
    if(method == 'n'){
        bnew = Network(x, y, lamb1, lamb2, b, r, a, u, n, p);
    }else if(method == 'm'){
        bnew = MCP(x, y, lamb1, b, r, n, p);
    }else{
        bnew = Elastic(x, y, lamb1, b, alpha, n, p);
    }
    diff = arma::accu(arma::abs(b - bnew))/(arma::accu(b != 0)+0.1);
    if(diff < tol){
      converged = true;
      break;
    }
    else{
      b = bnew;
      count++;
    }
  }
  return bnew;
}

// [[Rcpp::export]]
Rcpp::List RunLogit(arma::mat const &x, arma::vec const &y, double lamb1, double lamb2, arma::vec b, double r, arma::mat const &a, arma::vec const &triRowAbsSums, int p, double alpha, char method, int maxit, double tol)
{
  bool converged;
  int niter;
  double diff;
  arma::vec bnew = RunLogit_fit(x, y, lamb1, lamb2, b, r, a, triRowAbsSums, p, alpha, method, maxit, tol, converged, niter, diff);

  return Rcpp::List::create(Rcpp::Named("b") = bnew,
                            Rcpp::Named("converged") = converged,
                            Rcpp::Named("niter") = niter,
                            Rcpp::Named("diff") = diff);
}

// [[Rcpp::export]]
Rcpp::List RunNet(arma::mat& x, arma::vec& y, double lamb1, double lamb2, arma::vec b, double r, arma::mat& a, int p, int maxit, double tol)
{
  int count = 0, n = x.n_rows;
  int niter = 0;
  arma::vec bnew = b;
  bool converged = false;
  double diff = NA_REAL;
  while(count < maxit){
    Rcpp::checkUserInterrupt();
    niter++;
    bnew = Network(x, y, lamb1, lamb2, b, r, a, n, p);
    diff = arma::accu(arma::abs(b - bnew))/p;
    if(diff < tol){
      converged = true;
      break;
    }
    else{
      b = bnew;
      count++;
    }
  }
  return Rcpp::List::create(Rcpp::Named("b") = bnew,
                            Rcpp::Named("converged") = converged,
                            Rcpp::Named("niter") = niter,
                            Rcpp::Named("diff") = diff);
}

// [[Rcpp::export]]
Rcpp::List RunMCP(arma::mat& x, arma::vec& y, double lambda, arma::vec b, double r, int p, int maxit, double tol)
{
  int count = 0, n = x.n_rows;
  int niter = 0;
  arma::vec bnew = b;
  bool converged = false;
  double diff = NA_REAL;
  while(count < maxit){
    Rcpp::checkUserInterrupt();
    niter++;
    bnew = MCP(x, y, lambda, b, r, n, p);
    diff = arma::accu(arma::abs(b - bnew))/p;
    if(diff < tol){
      converged = true;
      break;
    }
    else{
      b = bnew;
      count++;
    }
  }
  return Rcpp::List::create(Rcpp::Named("b") = bnew,
                            Rcpp::Named("converged") = converged,
                            Rcpp::Named("niter") = niter,
                            Rcpp::Named("diff") = diff);
}

// [[Rcpp::export]]
Rcpp::List RunElastic(arma::mat& x, arma::vec& y, double lambda, arma::vec b, double alpha, int p, int maxit, double tol)
{
  int count = 0, n = x.n_rows;
  int niter = 0;
  arma::vec bnew = b;
  bool converged = false;
  double diff = NA_REAL;
  while(count < maxit){
    Rcpp::checkUserInterrupt();
    niter++;
    bnew = Elastic(x, y, lambda, b, alpha, n, p);
    diff = arma::accu(arma::abs(b - bnew))/p;
    if(diff < tol){
      converged = true;
      break;
    }
    else{
      b = bnew;
      count++;
    }
  }
  return Rcpp::List::create(Rcpp::Named("b") = bnew,
                            Rcpp::Named("converged") = converged,
                            Rcpp::Named("niter") = niter,
                            Rcpp::Named("diff") = diff);
}
