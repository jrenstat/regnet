#' fit a regression for given lambda with network-based regularization
#'
#' Network-based penalization regression for given values of \eqn{\lambda_{1}} and \eqn{\lambda_{2}}.
#' Typical usage is to have the cv.regnet function compute the optimal lambdas, then provide the resulting
#' \code{cv.regnet} object or selected lambdas to the regnet function. Users could also use MCP or Lasso.
#'
#' @keywords models
#' @param X a matrix of predictors without intercept. Each row should be an observation vector. A column of 1 will be added to the X matrix
#' by the program as the intercept.
#' @param ... other arguments passed to S3 methods.
#' @param Y the response variable. For response="binary", Y should be a numeric vector with zeros and ones. For response="survival", Y should be a
#' two-column matrix with columns named 'time' and 'status'. The latter is a binary variable, with '1' indicating an event, and '0'
#' indicating censoring.
#' @param response the response type. regnet supports three types of response: "binary", "continuous" and "survival".
#' @param penalty the penalty type. regnet provides three choices for the penalty function: "network", "mcp" and "lasso".
#' @param lamb.1 the tuning parameter \eqn{\lambda_{1}} that imposes sparsity.
#' @param lamb.2 the tuning parameter \eqn{\lambda_{2}} that controls the smoothness among coefficient profiles. \eqn{\lambda_{2}} is  needed
#' for network penalty.
#' @param r the regularization parameter in MCP. For binary response, r should be larger than 4.
#' @param clv a value or a vector, indexing variables that are not subject to penalty. clv only works for continuous and survival responses
#' for now, and will be ignored for other types of responses.
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
#' @param debugging a logical flag. If TRUE, extra information will be returned.
#'
#' @details The current version of regnet supports three types of responses: "binary", "continuous" and "survival".
#' \itemize{
#' \item regnet(..., response="binary", penalty="network") fits a network-based penalized logistic regression.
#' \item regnet(..., response="continuous", penalty="network") fits a network-based least square regression.
#' \item regnet(..., response="survival", penalty="network", robust=TRUE) fits a robust regularized AFT model using network penalty.
#' }
#' By default, regnet uses non-robust methods for all types of responses. To use robust methods, simply set robust=TRUE. It is recommended to use robust methods for survival response. Please see the references for more details about the models. Users could also use MCP or Lasso penalty.
#'
#' The coefficients are always estimated on a standardized X matrix. regnet standardizes each column of X to have unit variance
#' (using 1/n rather than 1/(n-1) formula). If the coefficients on the original scale are needed, the user can refit a standard model
#' using the subset of variables that have non-zero coefficients.
#'
#' @return an object of class "regnet" is returned, which is a list with components:
#' \itemize{
#' \item coeff: a vector of estimated coefficients. Please note that, if there are variables not subject to penalty (indicated by clv),
#' the order of returned vector is c(Intercept, unpenalized coefficients of clv variables, penalized coefficients of other variables).
#' \item Adj: a matrix of adjacency measures of the identified genetic variants. Identified genetic variants are those that have non-zero estimated coefficients. The adjacency matrix is returned for all penalties, but it affects estimation only when \code{penalty="network"}.
#' \item convergence: if \code{debugging=TRUE}, a list containing convergence diagnostics: \code{converged}, \code{niter}, and \code{diff}.
#' }
#'
#' @references Ren, J., He, T., Li, Y., Liu, S., Du, Y., Jiang, Y., and Wu, C. (2017).
#' Network-based regularization for high dimensional SNP data in the case-control study of
#' Type 2 diabetes. {\emph{BMC Genetics}, 18(1):44} \doi{10.1186/s12863-017-0495-5}
#'
#' Ren, J., Du, Y., Li, S., Ma, S., Jiang,Y. and Wu, C. (2019). Robust network-based regularization
#' and variable selection for high dimensional genomics data in cancer prognosis.
#' {\emph{Genet. Epidemiol.}, 43:276-291} \doi{10.1002/gepi.22194}
#'
#' @seealso \code{\link{cv.regnet}}
#'
#' @examples
#' ## Survival response
#' data(SurvExample)
#' X = rgn.surv$X
#' Y = rgn.surv$Y
#' clv = c(1:5) # variables 1 to 5 are clinical variables which we choose not to penalize.
#' penalty = "network"
#' fit = regnet(X, Y, "survival", penalty, rgn.surv$lamb1, rgn.surv$lamb2, clv=clv, robust=TRUE)
#' index = which(rgn.surv$beta[-(1:6)] != 0)
#' pos = which(fit$coeff[-(1:6)] != 0)
#' tp = length(intersect(index, pos))
#' fp = length(pos) - tp
#' list(tp=tp, fp=fp)
#'
#' @export
regnet <- function(X, ...)
{
  UseMethod("regnet")
}

#' @rdname regnet
#' @export
regnet.default <- function(X, Y, response=c("binary", "continuous", "survival"), penalty=c("network", "mcp", "lasso"), lamb.1=NULL, lamb.2=NULL,
                   r=NULL, clv=NULL, initiation=NULL, alpha.i=1, robust=FALSE,
                   adjacency=c("thresholded", "full"), adjacency.alpha=5,
                   maxit=20, tol=1e-3, debugging = FALSE, ...)
{
  # intercept = TRUE
  standardize=TRUE
  response = match.arg(response)
  penalty = match.arg(penalty)
  adjacency = match.arg(adjacency)
  conv = validate_convergence_control(maxit, tol)
  X = as.matrix(X)
  if(penalty != "network"){
    lamb.2 = 0
  }else{
    if(ncol(X)<3) stop("too less variables for network penalty.")
  }
  if(is.null(lamb.1)) stop("Both lambda1 and lambda2 need to be provided")
  if(penalty == "network" && is.null(lamb.2)) stop("Lambda2 needs to be provided for network method")

  this.call = match.call()
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
    if(!all(Y%in% c(0,1))) stop("Y has to be a binary variable of 1 and 0.")
    if(robust){
      message("robust methods are not available for ", response, " response.")
      robust = FALSE
    }
  }
  if(alpha.i>1 | alpha.i<0) stop("alpha.i should be between 0 and 1")
  # if(is.null(initiation)){
  #   if(response == "survival") initiation = "zero" else initiation = "elnet"
  # }

  if(is.null(r)) r = 5
  alpha = alpha.i # temporary

  out=switch (response,
              "binary" = LogitCD(X, Y, penalty, lamb.1, lamb.2, r, alpha, init=initiation, alpha.i,standardize, adjacency, adjacency.alpha, conv$maxit, conv$tol),
              "continuous" = ContCD(X, Y, penalty, lamb.1, lamb.2, clv, r, alpha, init=initiation, alpha.i, robust, standardize, debugging, adjacency, adjacency.alpha, conv$maxit, conv$tol),
              "survival" = SurvCD(X, Y0, status, penalty, lamb.1, lamb.2, clv, r, init=initiation, alpha.i, robust, standardize, debugging, adjacency, adjacency.alpha, conv$maxit, conv$tol)
  )
  para = list(response=response, penalty=penalty, robust=robust, lamb.1=lamb.1, lamb.2=lamb.2,
              r=r, clv=clv, initiation=initiation, alpha.i=alpha.i,
              adjacency=adjacency, adjacency.alpha=adjacency.alpha,
              maxit=conv$maxit, tol=conv$tol)
  fit = list(call = this.call, coeff = out$b, Adj=out$Adj, para=para)
  if(debugging) fit$convergence = out$convergence
  class(fit) = "regnet"
  fit
}

#' Fit a regnet model from cross-validation results
#'
#' This method refits \code{\link{regnet}} using the data and model settings
#' stored in a \code{\link{cv.regnet}} object. The tuning parameter values must
#' be supplied explicitly through \code{lamb.1} and, for network models,
#' \code{lamb.2}.
#'
#' @param X a \code{cv.regnet} object returned by \code{\link{cv.regnet}}.
#' @param lamb.1 the tuning parameter \eqn{\lambda_{1}} that imposes sparsity.
#' @param lamb.2 the tuning parameter \eqn{\lambda_{2}} that controls smoothness among coefficient profiles. \eqn{\lambda_{2}} is needed for network penalty.
#' @param ... additional named arguments passed to \code{regnet.default}. These can be used to override settings stored in \code{X}.
#'
#' @return an object of class \code{regnet}.
#'
#' @seealso \code{\link{cv.regnet}}, \code{\link{regnet}}
#'
#' @examples
#' \donttest{
#' data(LogisticExample)
#' X = rgn.logi$X
#' Y = rgn.logi$Y
#' cvfit = cv.regnet(X, Y, response="binary", penalty="network", folds=5, r = 4.5)
#' fit = regnet(cvfit, cvfit$lambda[1,1], cvfit$lambda[1,2])
#' }
#'
#' @export
regnet.cv.regnet <- function(X, lamb.1=NULL, lamb.2=NULL, ...)
{
  cvfit = X
  if(is.null(cvfit$data) || is.null(cvfit$data$X) || is.null(cvfit$data$Y)){
    stop("The cv.regnet object does not contain X and Y. Re-run cv.regnet() with keep.data=TRUE.")
  }

  para = cvfit$para
  if(is.null(para$response) || is.null(para$penalty)){
    stop("The cv.regnet object does not contain enough model information to refit.")
  }

  args = list(X=cvfit$data$X, Y=cvfit$data$Y, response=para$response, penalty=para$penalty,
              lamb.1=lamb.1, lamb.2=lamb.2, r=para$r, clv=para$clv,
              initiation=para$initiation, alpha.i=para$alpha.i, robust=para$robust,
              adjacency=para$adjacency, adjacency.alpha=para$adjacency.alpha,
              maxit=if(is.null(para$maxit)) 20 else para$maxit,
              tol=if(is.null(para$tol)) 1e-3 else para$tol)
  extra = list(...)
  if(length(extra)){
    if(is.null(names(extra)) || any(names(extra) == "")) stop("Additional arguments in ... must be named")
    args[names(extra)] = extra
  }
  fit = do.call(regnet.default, args)
  fit$call = match.call()
  fit$cv.call = cvfit$call
  fit
}
