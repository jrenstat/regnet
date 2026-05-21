
CV.Surv <- function(X0, Y0, status, penalty=c("network", "mcp", "lasso"), lamb.1=NULL, lamb.2=NULL, folds=5, foldid=NULL, clv=NULL, r=5,
                    init=NULL, alpha.i=1, robust=TRUE, standardize=TRUE, ncores, verbo = FALSE, debugging = FALSE,
                    adjacency=c("thresholded", "full"), adjacency.alpha=5,
                    maxit=20, tol=1e-3)
{
  intercept = TRUE
  adjacency = match.arg(adjacency)
  conv = validate_convergence_control(maxit, tol)
  clv.info = setup_clv(clv, ncol(X0))
  clv = clv.info$internal

  n = nrow(X0); p.c = length(clv); p = ncol(X0)-p.c+intercept;
  fold.info = make_foldid(n, folds, foldid)
  foldid = fold.info$foldid
  folds = fold.info$folds
  V0 = apply(X0, 2, function(t) stats::sd(t)*sqrt((n-1)/n)); 
  if(any(V0==0) & (penalty == "network")) stop("X columns have standard deviation equal zero");
  if(standardize){
    V0[V0==0|is.na(V0)]=1
    X1 = scale(X0, center = TRUE, scale = V0)
  }
  if(intercept) X1 = cbind(Intercept = rep(1, n), X1)
  Y1 = Y0

  out = KMweight(X1, Y1, status, robust)
  X = out$X + 10^-9
  Y = out$Y
  init = match.arg(init, choices = c("zero","elnet","cox"))

  # cat("var(Y): ", stats::var(Y), "\n")
  # if(is.null(lamb.1)){
  #   u=abs(t(X) %*% Y)
  #   if(!robust && penalty!="network") u=u/sqrt(n)
  #   LL = log(stats::quantile(u, 0.1)); UL = log(max(u))
  #   lamb.1 = rev(exp(seq(LL,UL,length.out = 35)))
  #   cat("XtY: ", range(lamb.1), "\n")
  #   # lamb.1=NULL
  # }

  if(is.null(lamb.1)){
    lasso.fit = glmnet::glmnet(X, Y, family="gaussian", nlambda=10)
    LL = log(min(lasso.fit$lambda)); UL = log(max(lasso.fit$lambda))
    lamb.1 = rev(exp(seq(LL,UL,length.out = 35))) * max(1,log10(1/stats::var(Y)))
    # cat("glm: ", range(rev(exp(seq(LL,UL,length.out = 35)))), "\n")
    # cat("glm: ", range(lamb.1), "\n")
  }

  if(is.null(lamb.2)){
    if(robust){
      lamb.2 = c(0.001, 0.01, 0.1, 1)
    }else{
      lamb.2 = c(0.01, 0.1, 1, 10)
    }
  }
  CVM = matrix(0, length(lamb.1), length(lamb.2));

  if(init == "zero"){
    b0 = rep(0, (p+p.c))
  } else if(init == "elnet"){
    b0 = initiation(X, Y, alpha.i)
  } else{
    b0 = initiation_cox(out$Xo, out$Yo, out$So)
  }
  a = Adjacency(X[,-clv,drop=FALSE], alpha=adjacency.alpha, type=adjacency)
  method = substr(penalty, 1, 1)
  #---------------------------------------------- Main Loop -----------------------------------------
  for(f in 1:folds){
    if(verbo) cat("CrossValidation: ",f, "/", folds, "\n")
    test = which(foldid == f)

    x = X[-test,,drop=FALSE]; y = Y[-test]
    x2 = X[test,,drop=FALSE]; y2 = Y[test]

    # if(!robust){
    #   x = scale(x, center = TRUE, scale = FALSE)
    #   x2 = scale(x2, center = TRUE, scale = FALSE)
    # }
    # if(init == "elnet") b0 = initiation(x, y, alpha.i)

    x.c=x[, clv, drop = FALSE]; x.g = x[, -clv, drop = FALSE];
    check_clv_rank(x.c)
    x2 = cbind(x2[,clv, drop = FALSE], x2[,-clv, drop = FALSE])

    # if(ncores>1){
      # CVM = CVM + SurvGrid_MC(x.c, x.g, y, x2, y2, lamb.1, lamb.2, b0[clv], b0[-clv], r, a, p, p.c, robust, method, ncores, debugging)
    # }else{
      CVM = CVM + SurvGrid(x.c, x.g, y, x2, y2, lamb.1, lamb.2, b0[clv], b0[-clv], r, a, p, p.c, robust, method, debugging, conv$maxit, conv$tol)
    # }


  }
  CVM = CVM/n
  mcvm = min(CVM)
  inds = which(CVM == mcvm, arr.ind=TRUE)
  lambda = lambda1 = lamb.1[inds[,1]]
  lambda2 = lamb.2[inds[,2]]
  if(length(lambda)>1){
    message("multiple optimal values(pairs) of lambda(s) are found. Lambda sequence (lamb.1/lamb.2) may need to be adjusted.")
  }else if(lambda==lamb.1[1]){
    message("Values of the lambda sequence (lamb.1) maybe too small.")
  }else if(lambda==lamb.1[length(lamb.1)]){
    message("Values of the lambda sequence (lamb.1) maybe too large.")
  }
  rownames(CVM) = signif(lamb.1, digits = 3)
  if(penalty == "network"){
    lambda = cbind(lambda1, lambda2)
    colnames(CVM) = lamb.2
  }
  outlist = list(lambda=lambda, mcvm=mcvm, CVM=CVM, penalty=penalty)
  class(outlist) = "cv.surv"
  outlist
}
