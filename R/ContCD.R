
ContCD <- function(X, y, penalty=c("network", "mcp", "lasso"), lamb.1=NULL, lamb.2=NULL, clv=NULL, r=5, alpha=1,
                    init=NULL, alpha.i=1, robust=FALSE, standardize=TRUE, debugging=FALSE,
                    adjacency=c("thresholded", "full"), adjacency.alpha=5,
                    maxit=20, tol=1e-3)
{
  intercept = TRUE
  adjacency = match.arg(adjacency)
  conv = validate_convergence_control(maxit, tol)
  clv.info = setup_clv(clv, ncol(X))
  clv = clv.info$internal

  n = nrow(X); p.c = length(clv); p = ncol(X)-p.c+intercept;
  vname = colnames(X)
  x = as.matrix(X); y = as.matrix(y)
  b0 = rep(0, p+p.c)
  method = substr(penalty, 1, 1)
  #---------------------------------------------- Main Loop -----------------------------------------
  V0 = apply(X, 2, function(t) stats::sd(t)*sqrt((n-1)/n));
  if(any(V0==0) & (penalty == "network")) stop("X columns have standard deviation equal zero");
  if(standardize){
    V0[V0==0|is.na(V0)]=1
    X = scale(X, center = TRUE, scale = V0)
  }
  X = cbind(1, X)
  init = match.arg(init, choices = c("elnet","zero"))
  if(init == "elnet") b0 = initiation(X, y, alpha.i, "gaussian")

  x.c=X[, clv, drop = FALSE]; x.g = X[, -clv, drop = FALSE]
  check_clv_rank(x.c)
  # if(penalty == "network") a = Adjacency(x.g) else a = as.matrix(0)
  a = Adjacency(x.g, alpha=adjacency.alpha, type=adjacency)

  if(robust){
    fit = RunCont_robust(x.c, x.g, y, lamb.1, lamb.2, b0[clv], b0[-clv], r, a, p, p.c, method, debugging, conv$maxit, conv$tol)
  }else{
    triRowAbsSums = rowSums(abs(a*upper.tri(a, diag = FALSE)))
    fit = RunCont(x.c, x.g, y, lamb.1, lamb.2, b0[clv], b0[-clv], r, a, triRowAbsSums, p, p.c, method, conv$maxit, conv$tol)
  }
  b = as.numeric(fit$b)
  convergence = list(converged=fit$converged, niter=fit$niter, diff=fit$diff)

  if(!is.null(vname)){
    names(b) = c("Intercept", vname[clv.info$original], vname[-clv.info$original])
  }else if(p.c==1){
    names(b) = c("Intercept", paste("g", seq = (1:p), sep=""))
  }else{
    names(b) = c("Intercept", paste("clv", seq = (1:(p.c-1)), sep=""), paste("g", seq = (1:p), sep=""))
  }

  sub = which(utils::tail(b,p)!=0)
  out = list(b=drop(b), Adj=a[sub,sub,drop=FALSE], convergence=convergence)
}
