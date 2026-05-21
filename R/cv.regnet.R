#' k-folds cross-validation for regnet
#'
#' This function does k-fold cross-validation for regnet and returns the optimal value(s) of lambda.
#'
#' @keywords models
#' @param X X matrix without intercept (see \code{\link{regnet}}).
#' @param Y the response variable Y (see \code{\link{regnet}}).
#' @param response the response type. regnet supports three types of response: "binary", "continuous" and "survival".
#' @param penalty the penalty type. regnet provides three choices for the penalty function: "network", "mcp" and "lasso".
#' @param lamb.1 a user-supplied sequence of \eqn{\lambda_{1}} values, which serves as a tuning parameter to impose sparsity.
#' If it is left as NULL, regnet will compute its own sequence.
#' @param lamb.2 a user-supplied sequence of \eqn{\lambda_{2}} values for network method. \eqn{\lambda_{2}} controls the smoothness
#' among coefficient profiles. If it is left as NULL, a default sequence will be used.
#' @param folds the number of folds for cross-validation; the default is 5.
#' @param foldid an optional integer vector of fold assignments with length equal to the number of observations.
#' If supplied, \code{folds} is ignored and the number of folds is inferred from \code{foldid}.
#' @param r the regularization parameter in MCP; default is 5. For binary response, r should be larger than 4.
#' @param clv a value or a vector, indexing variables that are not subject to penalty. clv only works for continuous and survival responses in the current version of regnet,
#' and will be ignored for other types of responses.
#' @param initiation the method for initiating the coefficient vector. The default method is elastic-net.
#' @param alpha.i the elastic-net mixing parameter. The program can use the elastic-net for choosing initial values of
#' the coefficient vector. alpha.i is the elastic-net mixing parameter, with 0 \eqn{\le} alpha.i \eqn{\le} 1. alpha.i=1 is the
#' lasso penalty, and alpha.i=0 is the ridge penalty. If the user chooses a method other than elastic-net for initializing
#' coefficients, alpha.i will be ignored.
#' @param robust a logical flag. Whether or not to use robust methods. Robust methods are available for survival and continuous response.
#' @param adjacency the adjacency construction rule. \code{"thresholded"} uses
#' \eqn{A = r^\alpha I(|r| > cutoff)}, while \code{"full"} uses
#' \eqn{A = r^\alpha}.
#' @param adjacency.alpha a positive integer power used in constructing the adjacency matrix. The default is 5.
#' @param maxit a positive integer specifying the maximum number of coordinate descent iterations.
#' @param tol a positive convergence tolerance for coordinate descent.
#' @param keep.data a logical flag. If \code{TRUE}, \code{X} and \code{Y} are stored in the returned object so \code{\link{regnet}} can refit from the \code{cv.regnet} object without resupplying the data.
#' @param verbo output progress to the console.
#' @param debugging a logical flag. If TRUE, extra information will be returned.
#'
#' @details When lamb.1 is left as NULL, regnet computes its own sequence. You can find the lamb.1 sequence used by the program in
#' the returned CVM matrix (see the 'Value' section). If you find the default sequence does not work well, you can try (1) standardizing
#' the response vector Y; or (2) providing a customized lamb.1 sequence for your data.
#'
#' Sometimes multiple optimal values(pairs) of lambda(s) can be found (see 'Value'). This is usually normal when the response is binary.
#' However, if the response is survival or continuous, you may want to check (1) if the sequence of lambda is too large
#' (i.e. all coefficients are shrunken to zero under all values of lambda) ; or (2) if the sequence is too small
#' (i.e. all coefficients are non-zero under all values of lambda). If neither, simply choose the value(pair) of lambda based on your preference.
#'
#' @return an object of class "cv.regnet" is returned, which is a list with components:
#' \item{lambda}{the optimal value(s) of \eqn{\lambda}. More than one value will be returned, if multiple lambdas have the cross-validated error =
#' min(cross-validated errors). If the network penalty is used, lambda contains optimal pair(s) of \eqn{\lambda_{1}} and \eqn{\lambda_{2}}.}
#' \item{mcvm}{the cross-validated error of the optimal \eqn{\lambda}. For binary response, the error is the misclassification rate.
#' For continuous response, mean squared error (MSE) is used. For survival response,
#' the MSE is used for non-robust methods, and the criterion for robust methods is the least absolute deviation (LAD).}
#' \item{CVM}{a matrix of the mean cross-validated errors of all lambdas used in the fits. The row names of CVM are the values of \eqn{\lambda_{1}}.
#' If the network penalty was used, the column names are the values of \eqn{\lambda_{2}}.}
#' \item{foldid}{the fold assignments used for cross-validation.}
#'
#' @references
#' Ren, J., Du, Y., Li, S., Ma, S., Jiang,Y. and Wu, C. (2019). Robust network-based regularization
#' and variable selection for high dimensional genomics data in cancer prognosis.
#' {\emph{Genet. Epidemiol.}, 43:276-291} \doi{10.1002/gepi.22194}
#'
#' Ren, J., He, T., Li, Y., Liu, S., Du, Y., Jiang, Y., and Wu, C. (2017).
#' Network-based regularization for high dimensional SNP data in the case-control study of
#' Type 2 diabetes. {\emph{BMC Genetics}, 18(1):44} \doi{10.1186/s12863-017-0495-5}
#'
#' @seealso \code{\link{regnet}}
#'
#' @examples
#' \donttest{
#' ## Binary response using network method
#' data(LogisticExample)
#' X = rgn.logi$X
#' Y = rgn.logi$Y
#' out = cv.regnet(X, Y, response="binary", penalty="network", folds=5, r = 4.5)
#' out$lambda
#' fit = regnet(X, Y, "binary", "network", out$lambda[1,1], out$lambda[1,2], r = 4.5)
#' index = which(rgn.logi$beta != 0)
#' pos = which(fit$coeff != 0)
#' tp = length(intersect(index, pos))
#' fp = length(pos) - tp
#' list(tp=tp, fp=fp)
#'
#' ## Binary response using MCP method
#' out = cv.regnet(X, Y, response="binary", penalty="mcp", folds=5, r = 4.5)
#' out$lambda
#' fit = regnet(X, Y, "binary", "mcp", out$lambda[1], r = 4.5)
#' index = which(rgn.logi$beta != 0)
#' pos = which(fit$coeff != 0)
#' tp = length(intersect(index, pos))
#' fp = length(pos) - tp
#' list(tp=tp, fp=fp)
#' }
#'
#' @export

cv.regnet <- function(X, Y, response=c("binary", "continuous", "survival"), penalty=c("network", "mcp", "lasso"), lamb.1=NULL, lamb.2=NULL,
                      folds=5, foldid=NULL, r=NULL, clv=NULL, initiation=NULL, alpha.i=1, robust=FALSE,
                      adjacency=c("thresholded", "full"), adjacency.alpha=5, maxit=20, tol=1e-3, keep.data=TRUE,
                      verbo = FALSE, debugging = FALSE)
{
  standardize=TRUE
  response = match.arg(response)
  penalty = match.arg(penalty)
  adjacency = match.arg(adjacency)
  conv = validate_convergence_control(maxit, tol)
  this.call = match.call()
  Y.input = Y
  X = as.matrix(X)

  if(response == "survival"){
    if(ncol(Y) != 2) stop("Y should be a two-column matrix")
    if(!setequal(colnames(Y), c("time", "status"))) stop("Y should be a two-column matrix with columns named 'time' and 'status'")
    Y0 = Y[,"time"]
    status = as.numeric(Y[,"status"])
    if(any(Y0<=0)) stop("survival times need to be positive")
    if(length(Y0) != nrow(X))  stop("the number of rows of Y does not match the number of rows of X");
    if(!all(status%in% c(0,1))) stop("status has to be a binary variable of 1 and 0")
  }else{
    if(length(Y) != nrow(X))  stop("length of Y does not match the number of rows of X");
  }

  if(response=="binary"){
    if(!all(Y%in% c(0,1))) stop("Y must be a binary variable of 1 and 0")
    if(robust){
      message("robust methods are not available for ", response, " response.")
      robust = FALSE
    }
  }

  fold.info = make_foldid(nrow(X), folds, foldid)
  folds = fold.info$folds
  foldid = fold.info$foldid

  if(alpha.i>1 | alpha.i<0) stop("alpha.i should be between 0 and 1")
  # if(is.null(initiation)){
  #   if(response == "survival") initiation = "zero" else initiation = "elnet"
  # }

  if(is.null(r)) r = 5
  if(penalty != "network"){
    lamb.2 = 0
  }else{
    if(ncol(X)<3) stop("too less variables for network penalty.")
  }
  alpha = alpha.i # temporarily
  # ncores = as.integer(ncores)
  # if(ncores<1) stop("incorrect value of ncores")
  ncores=1
  
  fit=switch (response,
    "binary" = CV.Logit(X, Y, penalty, lamb.1, lamb.2, folds, foldid, r, alpha, init=initiation, alpha.i, standardize, ncores, verbo, debugging, adjacency, adjacency.alpha, conv$maxit, conv$tol),
    "continuous" = CV.Cont(X, Y, penalty, lamb.1, lamb.2, folds, foldid, clv=clv, r, alpha, init=initiation, alpha.i, robust, standardize, ncores, verbo, debugging, adjacency, adjacency.alpha, conv$maxit, conv$tol),
    "survival" = CV.Surv(X, Y0, status, penalty, lamb.1, lamb.2, folds, foldid, clv=clv, r, init=initiation, alpha.i, robust, standardize, ncores, verbo, debugging, adjacency, adjacency.alpha, conv$maxit, conv$tol)
  )
  fit$call = this.call
  fit$foldid = foldid
  fit$para = list(response=response, penalty=penalty, robust=robust, r=r, clv=clv,
                  initiation=initiation, alpha.i=alpha.i, folds=folds,
                  adjacency=adjacency, adjacency.alpha=adjacency.alpha,
                  maxit=conv$maxit, tol=conv$tol,
                  lamb.1=lamb.1, lamb.2=lamb.2)
  if(keep.data) fit$data = list(X=X, Y=Y.input)
  class(fit) = "cv.regnet"
  fit
}
