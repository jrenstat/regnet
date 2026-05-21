
CV.Logit <- function(X, Y, penalty=c("network", "mcp", "lasso"), lamb.1=NULL, lamb.2=NULL, folds=5, foldid=NULL, r=5, alpha=1,
                     init=NULL, alpha.i=1, standardize=TRUE, ncores=1, verbo = FALSE, debugging = FALSE,
                    adjacency=c("thresholded", "full"), adjacency.alpha=5,
                    maxit=20, tol=1e-3)
{
  adjacency = match.arg(adjacency)
  conv = validate_convergence_control(maxit, tol)
  # if(is.null(lamb.1)){
  #   lamb.1 = switch (penalty,
  #                    "network" = lambda.n,
  #                    "mcp" = lambda.m,
  #                    "lasso" = lambda.l)
  # }
  if(is.null(lamb.2)) lamb.2 = c(0.1, 1, 10)
  init = match.arg(init, choices = c("elnet","zero"))

  n = nrow(X); p = ncol(X);
  X = as.matrix(X); Y = as.matrix(Y)
  fold.info = make_foldid(n, folds, foldid)
  foldid = fold.info$foldid
  folds = fold.info$folds
  # X = scale(X, center = TRUE, scale = FALSE)
  
  V0 = apply(X, 2, function(t) stats::sd(t)*sqrt((n-1)/n));
  if(any(V0==0) & (penalty == "network")) stop("X columns have standard deviation equal zero");
  if(standardize){
    V0[V0==0|is.na(V0)]=1
    X = scale(X, center = TRUE, scale = V0)
  }

  if(is.null(lamb.1)){
  # if(TRUE){
    lasso.fit = glmnet::glmnet(X, Y, family="binomial", nlambda=10); # lasso.fit$lambda
    LL = log(min(lasso.fit$lambda)); UL = log(max(lasso.fit$lambda)+log10(max(1,lamb.2))); # exp(UL)
    # LL = log(min(lasso.fit$lambda)); UL = log(max(lasso.fit$lambda)*max(1,log(lamb.2))); # exp(UL)
    # UL = log(max(lasso.fit$lambda)); exp(UL)
    lamb.1 = rev(exp(seq(LL,UL,length.out = 45)))
  }

  b0 = rep(0, p+1)
  CVM = CVM2 = matrix(0, length(lamb.1), length(lamb.2))
  method = substr(penalty, 1, 1)
  if(penalty == "network") a = Adjacency(X, alpha=adjacency.alpha, type=adjacency) else a = as.matrix(0)
  #---------------------------------------------- Main Loop -----------------------------------------
  for(f in 1:folds){
    if(verbo) cat("CrossValidation: ",f, "/", folds, "\n")
    test = which(foldid == f)

    x = X[-test,,drop=FALSE]; y = Y[-test]
    x2 = X[test,,drop=FALSE]; y2 = Y[test]
    if(standardize){
      V1 = apply(x, 2, function(t) stats::sd(t)*sqrt((n-1)/n)); V1[V1==0|is.na(V1)]=1
      V2 = apply(x2, 2, function(t) stats::sd(t)*sqrt((n-1)/n)); V2[V2==0|is.na(V2)]=1
      x = scale(x, center = FALSE, scale = V1 )
      x2 = scale(x2, center = FALSE, scale = V2)
    }
    # if(penalty == "network") a = Adjacency(x) else a = as.matrix(0)

    x = cbind(1, x); x2 = cbind(1, x2)
    if(init == "elnet") b0 = initiation(x, y, alpha.i, "binomial")

    # if(ncores>1){
      # CVM = CVM + LogitGrid_MC(x, y, x2, y2, lamb.1, lamb.2, b0, r, a, p, alpha, method, ncores)
    # }else{
      CVMs = LogitGrid(x, y, x2, y2, lamb.1, lamb.2, b0, r, a, p, alpha, method, conv$maxit, conv$tol)
      CVM = CVM + CVMs$CVM
      CVM2 = CVM2 + CVMs$CVM2
    # }

  }

  CVM = CVM/n
  mcvm = min(CVM)
  inds = which(CVM == mcvm, arr.ind=TRUE)
  inds0 = NULL
  if(penalty == "network"){ #*&%*@!
    inds = inds[!duplicated(inds[,2]),,drop=FALSE]; inds0=inds
    inds = inds[which(CVM2[inds] == min(CVM2[inds])),,drop=FALSE]
  }

  lambda = lambda1 = lamb.1[inds[,1]]
  lambda2 = lamb.2[inds[,2]]
  if(length(lambda)>1){
    message("multiple optimal values(pairs) of lambda(s) are found.")
  }else if(lambda==lamb.1[1]){
    message("Values of the lambda sequence (lamb.1) maybe too small.")
  }else if(lambda==lamb.1[length(lamb.1)]){
    message("Values of the lambda sequence (lamb.1) maybe too large.")
  }
  rownames(CVM) = rownames(CVM2) = signif(lamb.1, digits = 3)
  if(penalty == "network"){
    lambda = cbind(lambda1, lambda2)
    colnames(CVM) = colnames(CVM2) = lamb.2
  }
  outlist = list(lambda=lambda, mcvm=mcvm, CVM=CVM, penalty=penalty)
  if(debugging){
    para = list(lamb.1=lamb.1, lamb.2=lamb.2, inds=inds, inds0=inds0, CVM2=CVM2)
    outlist$para = para
  }
  class(outlist) = "cv.logit"
  outlist
}
