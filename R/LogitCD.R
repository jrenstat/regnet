
LogitCD <- function(X, Y, penalty=c("network", "mcp", "lasso"), lamb.1=NULL, lamb.2=NULL, r=5, alpha=1,
                     init=NULL, alpha.i=1, standardize=TRUE,
                     adjacency=c("thresholded", "full"), adjacency.alpha=5)
{
  adjacency = match.arg(adjacency)
  n = nrow(X); p = ncol(X);
  x = as.matrix(X); y = as.matrix(Y)
  vname = colnames(x)
  b0 = rep(0, p+1)
  method = substr(penalty, 1, 1)
  #---------------------------------------------- Main Loop -----------------------------------------
  V0 = apply(x, 2, function(t) stats::sd(t)*sqrt((n-1)/n));
  if(any(V0==0) & (penalty == "network")) stop("X columns have standard deviation equal zero");
  if(standardize){
    V0[V0==0|is.na(V0)]=1
    x = scale(x, center = TRUE, scale = V0)
  }

  a = Adjacency(x, alpha=adjacency.alpha, type=adjacency)
  x = cbind(rep(1,n), x)
  init = match.arg(init, choices = c("elnet","zero"))
  if(init == "elnet") b0 = initiation(x, y, alpha.i, "binomial")

  triRowAbsSums = rowSums(abs(a*upper.tri(a, diag = FALSE)))
  b = RunLogit(x, y, lamb.1, lamb.2, b0, r, a, triRowAbsSums, p, alpha, method)
  b = as.numeric(b)

  if(!is.null(vname)){
    names(b) = c("Intercept", vname)
  }else{
    names(b) = c("Intercept", paste("v", seq = (1:p), sep=""))
  }

  sub = which(b[-1]!=0)
  out = list(b=drop(b), Adj=a[sub,sub,drop=FALSE])
}
