
#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include"RunLogistic.h"
#include"Utilities.h"


using namespace Rcpp;
using namespace arma;
//using namespace R;

// [[Rcpp::export()]]
Rcpp::List LogitGrid(arma::mat const &x, arma::vec const &y, arma::mat const &x2, arma::vec const &y2, arma::vec const &lamb1, arma::vec const &lamb2, arma::vec b, double r, arma::mat const &a, int p, double alpha, char method, int maxit, double tol)
{
	arma::vec btmp, triRowAbsSums, mc(2);
	Rcpp::List fit;
    arma::mat CVM(lamb1.n_elem, lamb2.n_elem, fill::zeros), CVM2(lamb1.n_elem, lamb2.n_elem, fill::zeros);
	if(method == 'n') triRowAbsSums = TriRowAbsSums(a);
	
		for(unsigned int j=0; j<lamb2.n_elem ; j++){
			for(unsigned int i=0; i<lamb1.n_elem ; i++){
				fit = RunLogit(x, y, lamb1(i), lamb2(j), b, r, a, triRowAbsSums, p, alpha, method, maxit, tol);
				btmp = Rcpp::as<arma::vec>(fit["b"]);
	            CVM(i, j) = validation_logit(x2, y2, btmp);
				mc = validation_logit2(x2, y2, btmp);
			CVM(i, j) = mc(0);
			CVM2(i, j) = mc(1);
			Rcpp::checkUserInterrupt();
		}
	}
    // return CVM;
	return Rcpp::List::create(Rcpp::Named("CVM") = CVM,
							  Rcpp::Named("CVM2") = CVM2);
}

/* 
// [[Rcpp::export()]]
arma::mat LogitGrid_MC(arma::mat const &x, arma::vec const &y, arma::mat const &x2, arma::vec const &y2, arma::vec const &lamb1, arma::vec const &lamb2, arma::vec b, double r, arma::mat const &a, int p, double alpha, char method, int ncores)
{
	arma::vec btmp, triRowAbsSums;
    arma::mat CVM(lamb1.n_elem, lamb2.n_elem, fill::zeros);
	if(method == 'n') triRowAbsSums = TriRowAbsSums(a);
	
	omp_set_num_threads(ncores);
	#pragma omp parallel for collapse(2) private(btmp)
	for(unsigned int j=0; j<lamb2.n_elem ; j++){
		for(unsigned int i=0; i<lamb1.n_elem ; i++){
			btmp = RunLogit(x, y, lamb1(i), lamb2(j), b, r, a, triRowAbsSums, p, alpha, method);
            CVM(i, j) = validation_logit(x2, y2, btmp);
		}
	}
    return CVM;
}
 */
 
 
 
