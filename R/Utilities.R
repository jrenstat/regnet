lambda.n = rev(exp(seq(1,45,1)/5 -7))
lambda.m = rev(exp(seq(1,45,1)/5 -7))
#lambda.n = rev(exp(seq(1,45,1)/4 -9))
#lambda.m = rev(exp(seq(1,45,1)/4 -9))
lambda.e = rev(exp(seq(1,45,1)/4 -9))
lambda.l = rev(exp(seq(1,45,1)/4 -9))

initiation <- function(x, y, alpha, family="gaussian"){
  lasso.cv <- suppressWarnings(glmnet::cv.glmnet(x,y, family=family, alpha=alpha, nfolds=5, nlambda=50))
  lambda <- lasso.cv$lambda.min
  lasso.fit <- glmnet::glmnet(x, y, family, alpha=alpha, nlambda=50)
  coef0 <- as.vector(stats::predict(lasso.fit, s=lambda, type="coefficients"))[-1]
}

initiation_cox <- function(x, y0, d){
  y = cbind(time = y0, status = d)
  lasso.cv = suppressWarnings(glmnet::cv.glmnet(x, y, alpha=1, family="cox", nfolds=5, nlambda=30, standardize=FALSE))
  alpha = 2*(lasso.cv$lambda.min)
  lasso.fit = glmnet::glmnet(x,y,family="cox", alpha=1, nlambda=30, standardize=FALSE)
  coef0 = as.numeric(stats::predict(lasso.fit, s=alpha, type="coefficients"))
}

TruePos <- function(b, b.true){
  index = which(b.true != 0)
  pos = which(b != 0)
  tp = length(intersect(index, pos))
  fp = length(pos) - tp
  list(tp=tp, fp=fp)
}

make_foldid <- function(n, folds=5, foldid=NULL){
  if(is.null(foldid)){
    folds = as.integer(folds)
    if(n < folds) stop("sample size too small for ", folds, "-fold cross-validation.")
    if(folds < 2) stop("incorrect value of folds")
    rs = sample(seq_len(n))
    foldid = integer(n)
    fold.size = floor(n/folds)
    remainder = n%%folds
    start = 1L
    for(f in seq_len(folds)){
      size = fold.size + as.integer(f > folds - remainder)
      test = rs[start:(start + size - 1L)]
      foldid[test] = f
      start = start + size
    }
  }else{
    if(length(foldid) != n) stop("foldid length must equal the number of observations")
    if(any(is.na(foldid))) stop("foldid cannot contain missing values")
    if(!is.numeric(foldid)) stop("foldid should be an integer vector")
    if(any(foldid != as.integer(foldid))) stop("foldid should be an integer vector")
    foldid = match(as.integer(foldid), sort(unique(as.integer(foldid))))
    folds = length(unique(foldid))
    if(folds < 2) stop("foldid should contain at least two folds")
  }
  list(foldid = as.integer(foldid), folds = as.integer(folds))
}

setup_clv <- function(clv, p){
  if(is.null(clv)){
    clv0 = integer(0)
  }else{
    if(!is.numeric(clv)) stop("clv should contain valid column indices of X")
    if(any(is.na(clv)) || any(clv != as.integer(clv))) stop("clv should contain valid column indices of X")
    clv0 = unique(as.integer(clv))
    if(any(clv0 < 1) || any(clv0 > p)) stop("clv should contain valid column indices of X")
  }
  clv = c(1L, clv0 + 1L)
  list(original = clv0, internal = clv, penalized = setdiff(seq_len(p + 1L), clv))
}

check_clv_rank <- function(x, label="clv"){
  x = as.matrix(x)
  if(ncol(x) > nrow(x) || qr(x)$rank < ncol(x)){
    stop("Unpenalized variables in ", label, " are rank deficient")
  }
  invisible(TRUE)
}

validate_convergence_control <- function(maxit, tol){
  if(length(maxit) != 1 || !is.numeric(maxit) || is.na(maxit) || maxit < 1 || maxit != as.integer(maxit)){
    stop("maxit should be a positive integer")
  }
  if(length(tol) != 1 || !is.numeric(tol) || is.na(tol) || tol <= 0){
    stop("tol should be a positive number")
  }
  invisible(list(maxit=as.integer(maxit), tol=as.numeric(tol)))
}

Adjacency = function(x, alpha=5, type=c("thresholded", "full"))
{
  type = match.arg(type)
  if(length(alpha) != 1 || !is.numeric(alpha) || is.na(alpha) || alpha <= 0) stop("adjacency.alpha should be a positive number")
  if(abs(alpha - round(alpha)) > sqrt(.Machine$double.eps)) stop("adjacency.alpha should be a positive integer")
  n = nrow(x)
  p = ncol(x)
  if(p < 2){
    A = matrix(0, p, p)
    colnames(A) = rownames(A) = colnames(x)
    return(A)
  }
  r0 = stats::cor(x)
  r0[is.na(r0)] = 0
  diag(r0) = 1
  if(type == "full"){
    A = (r0)^alpha
    diag(A) = 0
    return(A)
  }
  r = r0; r[which(r==1)] = 1 - 0.01
  z = 0.5*log((1+r[upper.tri(r)])/(1-r[upper.tri(r)]))
  z = z[is.finite(z)]
  if(length(z) == 0 || n <= 3){
    cutoff = 1
  }else{
    c0 = mean(sqrt(n-3)*z) + 2*stats::sd(sqrt(n-3)*z)
    cutoff = (exp(2*c0/sqrt(n-3))-1)/(exp(2*c0/sqrt(n-3))+1)
    if(!is.finite(cutoff)) cutoff = 1
  }
  r = r0
  A = (r)^alpha*(abs(r)>cutoff)
  diag(A) = 0
  A
}

.onUnload <- function (libpath) {
  library.dynam.unload("regnet", libpath)
}
