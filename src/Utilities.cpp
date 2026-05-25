#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include"Utilities.h"

using namespace Rcpp;
using namespace arma;


double Soft(double z, double lambda){
  if(z > lambda) return(z - lambda);
  else if(z < -lambda) return(z + lambda);
  else return(0);
}

arma::vec sigmoid_stable(arma::vec const &eta){
  arma::vec out(eta.n_elem, fill::none);
  out.fill(datum::nan);

  arma::uvec pos = arma::find(eta >= 0);
  arma::uvec neg = arma::find(eta < 0);

  if(pos.n_elem > 0){
    out.elem(pos) = 1/(1 + arma::exp(-eta.elem(pos)));
  }
  if(neg.n_elem > 0){
    arma::vec e = arma::exp(eta.elem(neg));
    out.elem(neg) = e/(1 + e);
  }

  return(out);
}

arma::vec TriRowAbsSums(arma::mat const &a){
	int p = a.n_rows;
	arma::vec rowSums(p, fill::zeros);
	for(int m = 0; m < (p-1); m++){
		// rowSums(m) = arma::accu(arma::abs(a.row(m).subvec(m+1, p-1)));
		rowSums(m) = arma::accu(arma::abs(a.submat(m, m+1, m, p-1)));
	}
	return(rowSums);
}

//Non-robust MSE
double validation_LS(arma::mat const &x, arma::vec const &y, arma::vec const &b){
	//double mse = arma::accu(pow(y - x*b, 2))/y.n_elem;
    double mse = arma::accu(arma::square(y - x*b));
	return(mse);
}

double validation_LS(arma::mat const &xc, arma::mat const &xg, arma::vec const &y, arma::vec const &b, int p, int pc){
    double mse = arma::accu(arma::square(y - xc*b.subvec(0, pc-1)-xg*b.subvec(pc, pc+p-1)));
	return(mse);
}

double validation_LAD(arma::mat const &x, arma::vec const &y, arma::vec const &b){
	//double lad = arma::accu(arma::abs(y - x*b))/y.n_elem;
    double lad = arma::accu(arma::abs(y - x*b));
	return(lad);
}

double validation_LAD(arma::mat const &xc, arma::mat const &xg, arma::vec const &y, arma::vec const &bc, arma::vec const &bg){
    double lad = arma::accu(arma::abs(y - xc*bc - xg*bg));
	return(lad);
}

double validation_logit(arma::mat const &x0, arma::vec const &y0, arma::vec const &b){
    arma::vec yi = sigmoid_stable(x0*b);
    arma::uvec y = (yi > 0.5);
    double mc = arma::accu(arma::abs(y0 - y));
	return(mc);
}

arma::vec validation_logit2(arma::mat const &x0, arma::vec const &y0, arma::vec const &b){
    arma::vec yi = sigmoid_stable(x0*b), mc(2);
    mc(1) = arma::accu(arma::abs(y0 - yi));
	arma::uvec y = (yi > 0.5);
    mc(0) = arma::accu(arma::abs(y0 - y));
	return(mc);
}

arma::vec fastLm(arma::vec const &y, arma::mat const &X){
    arma::colvec coef = arma::solve(X, y);    
    return coef;
}
