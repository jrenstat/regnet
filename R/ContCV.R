
CV.Cont <- function(X, Y, penalty=c("network", "mcp", "lasso"), lamb.1=NULL, lamb.2=NULL, folds=5, foldid=NULL, clv=NULL,
                    r=5, alpha=1, init=NULL, alpha.i=1, robust=FALSE, standardize=TRUE, ncores=1, verbo = FALSE, debugging = FALSE,
                    adjacency=c("thresholded", "full"), adjacency.alpha=5,
                    maxit=20, tol=1e-3)
{
  intercept = TRUE
  adjacency = match.arg(adjacency)
  conv = validate_convergence_control(maxit, tol)
  clv.info = setup_clv(clv, ncol(X))
  clv = clv.info$internal
  n = nrow(X); p.c = length(clv); p = ncol(X)-p.c+intercept;
  X = as.matrix(X); Y = as.matrix(Y)
  fold.info = make_foldid(n, folds, foldid)
  foldid = fold.info$foldid
  folds = fold.info$folds

  V0 = apply(X, 2, function(t) stats::sd(t)*sqrt((n-1)/n)); 
  if(any(V0==0) & (penalty == "network")) stop("X columns have standard deviation equal zero");
  V0[V0==0|is.na(V0)]=1
  X = scale(X, center = TRUE, scale = V0)

  # if(is.null(lamb.2)) lamb.2 = c(0.1, 1, 10)
  if(is.null(lamb.2)){
    if(robust){
      lamb.2 = c(0.01, 0.1, 1)
    }else{
      lamb.2 = c(0.1, 1, 10)
    }
  }

  if(is.null(lamb.1)){
    lasso.fit = glmnet::glmnet(X, Y, family="gaussian", nlambda=10)
    LL = log(min(lasso.fit$lambda)); UL = log(max(lasso.fit$lambda)+log10(max(1,lamb.2)))
    lamb.1 = rev(exp(seq(LL,UL,length.out = 50)))
  }

  init = match.arg(init, choices = c("elnet","zero"))
  b0 = rep(0, (p+p.c))
  CVM = matrix(0, length(lamb.1), length(lamb.2))
  method = substr(penalty, 1, 1)
  if(penalty == "network") a = Adjacency(cbind(1, X)[,-clv,drop=FALSE], alpha=adjacency.alpha, type=adjacency) else a = as.matrix(0)

  # # OpenMP
  # Sys.setenv("PKG_CXXFLAGS" = "-fopenmp")
  # Sys.setenv("PKG_LIBS" = "-fopenmp")
  #----------------------------------------- Main Loop ---------------------------------------
  for(f in 1:folds){
    if(verbo) cat("CrossValidation: ",f, "/", folds, "\n")
    test = which(foldid == f)

    x = X[-test,,drop=FALSE]; y = Y[-test]
    x2 = X[test,,drop=FALSE]; y2 = Y[test]
    if(!robust){
      V1 = apply(x, 2, function(t) stats::sd(t)*sqrt((n-1)/n)); V1[V1==0|is.na(V1)]=1
      V2 = apply(x2, 2, function(t) stats::sd(t)*sqrt((n-1)/n)); V2[V2==0|is.na(V2)]=1
      x = scale(x, center = FALSE, scale = V1 )
      x2 = scale(x2, center = FALSE, scale = V2)
    }
    x = cbind(1, x); x2 = cbind(1, x2)
    if(init == "elnet") b0 = initiation(x, y, alpha.i, "gaussian") # which(is.na(x), arr.ind = TRUE)

    x.c=x[, clv, drop = FALSE]; x.g = x[, -clv, drop = FALSE];
    check_clv_rank(x.c)
    x2 = cbind(x2[,clv, drop = FALSE], x2[,-clv, drop = FALSE])
    # if(penalty == "network") a = Adjacency(x.g) else a = as.matrix(0)
	
    # if(ncores>1){
      # CVM = CVM + ContGrid_MC(x.c, x.g, y, x2, y2, lamb.1, lamb.2, b0[clv], b0[-clv], r, a, p, p.c, robust, method, ncores, debugging)
    # }else{
      CVM = CVM + ContGrid(x.c, x.g, y, x2, y2, lamb.1, lamb.2, b0[clv], b0[-clv], r, a, p, p.c, robust, method, debugging, conv$maxit, conv$tol)
    # }
  }
  CVM = CVM/n
  mcvm = min(CVM)
  inds = which(CVM == mcvm, arr.ind=TRUE)
  lambda = lambda1 = lamb.1[inds[,1]]
  lambda2 = lamb.2[inds[,2]]
  if(length(lambda)>1) message("multiple optimal values(pairs) of lambda(s) are found.")
  rownames(CVM) = signif(lamb.1, digits = 3)
  if(penalty == "network"){
    lambda = cbind(lambda1, lambda2)
    colnames(CVM) = lamb.2
  }
  outlist = list(lambda=lambda, mcvm=mcvm, CVM=CVM, penalty=penalty)
  class(outlist) = "cv.cont"
  outlist
}
