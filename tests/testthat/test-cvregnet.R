test_that("cv.regnet_data_format_cont", {
  n = 27; p = 3;
  X = scale(matrix(rnorm(n*p,0,5), n, p), scale=TRUE)
  Y= 5 + X%*%c(0,3,0) + rnorm(n)
  expect_error(cv.regnet(X[,1], Y, "c", "n"), "too less variables for network penalty")
  expect_error(cv.regnet(X[1:4,], Y[1:4], "c", "n"), "sample size too small")
  expect_error(cv.regnet(X, Y[1:4], "c", "n"), "length of Y does not match")
  expect_error(cv.regnet(X, Y, "c", "n", folds = -5), "incorrect value of folds")
  expect_error(cv.regnet(X, Y, "c", "n", alpha.i = 2), "alpha.i should be between 0 and 1")
  # expect_error(cv.regnet(X, Y, "c", "n", ncore = 0), "incorrect value of ncores")

  out = cv.regnet(X, Y, "c", "n", lamb.2 = 1)
  expect_equal(ncol(out$lambda), 2)
  expect_equal(colnames(out$CVM), "1")
  expect_equal(out$penalty, "network")

  out = cv.regnet(X, Y, "c", "n")
  expect_equal(colnames(out$CVM), c("0.1", "1", "10"))

  out = cv.regnet(X, Y, "c", "m")
  expect_equal(ncol(out$CVM), 1)
  expect_equal(out$penalty, "mcp")

  out = cv.regnet(X, Y, "c", "l")
  expect_equal(ncol(out$CVM), 1)
  expect_equal(out$penalty, "lasso")
  # fit = regnet(X, Y, "c", "n", out$lambda[1], lamb.2=1); fit$coeff
})

test_that("cv.regnet_data_format_surv", {
  n = 27; p = 3;
  X = scale(matrix(rnorm(n*p,0,1), n, p), scale=TRUE)
  Y0= exp(2 + X%*%c(0,3,0) + rnorm(n)); Y1 = sample(rep(c(0,1,1,1),n/2),n)
  Y = data.frame(time=(Y0+Y0*(Y1-1)*0.2), status=Y1)

  expect_error(cv.regnet(X, Y0, "s", "n"), "Y should be a two-column matrix")
  expect_error(cv.regnet(X, cbind(Y[,1], Y[,2]), "s", "n"), "columns named 'time' and 'status'")
  expect_error(cv.regnet(X, data.frame(time=Y[,1]-10, status=Y[,2]), "s", "n"), "survival times need to be positive")
  expect_error(cv.regnet(X, data.frame(time=Y[,1], status=Y[,2]*2), "s", "n"), "binary variable of 1 and 0")
  expect_error(cv.regnet(X, Y[1:4,], "s", "n"), "the number of rows of Y")

  out = cv.regnet(X, Y, "s", "n", lamb.2 = 1)
  expect_equal(ncol(out$lambda), 2)
  expect_equal(colnames(out$CVM), "1")
  expect_equal(out$penalty, "network")

  out = cv.regnet(X, Y, "s", "m")
  expect_equal(ncol(out$CVM), 1)
  expect_equal(out$penalty, "mcp")
  # fit = regnet(X, Y, "s", "n", out$lambda[1], lamb.2=1); fit$coeff
})

test_that("cv.regnet_data_format_logit", {
  n = 52; p = 3; Y = rep(0,n)
  X = scale(matrix(rnorm(n*p,0,5), n, p), scale=TRUE)
  Y0 = 1/(1+exp(-X%*%c(0,5,0)))
  Y = rbinom(n,1,Y0);

  expect_error(cv.regnet(X, Y0, NULL, "n"), "Y must be a binary variable")
  expect_error(cv.regnet(X, Y[1:4], "b", "n"), "length of Y does not match")
  expect_message(cv.regnet(X, Y, "b", "n", robust=TRUE), "robust methods are not available")
  out.robust = suppressMessages(cv.regnet(X, Y, "b", "n", robust=TRUE))
  expect_false(out.robust$para$robust)

  out = cv.regnet(X, Y, "b", "n", lamb.2 = 1)
  expect_equal(ncol(out$lambda), 2)
  expect_equal(colnames(out$CVM), "1")
  expect_equal(out$penalty, "network")

  out = cv.regnet(X, Y, "b", "m")
  expect_equal(ncol(out$CVM), 1)
  expect_equal(out$penalty, "mcp")
  # fit = regnet(X, Y, "b", "m", out$lambda[1]); fit$coeff

})

test_that("foldid and regnet.cv.regnet refit work", {
  n = 30; p = 5
  X = scale(matrix(rnorm(n*p,0,1), n, p), scale=TRUE)
  colnames(X) = paste0("x", 1:p)
  Y = 2 + X[,2] * 3 + rnorm(n)
  foldid = rep(1:5, length.out=n)

  out = cv.regnet(X, Y, "c", "l", lamb.1=c(0.5, 0.1), foldid=foldid)
  expect_equal(out$foldid, foldid)
  expect_equal(out$para$folds, 5)
  expect_equal(out$para$penalty, "lasso")

  fit = regnet(out, out$lambda[1])
  expect_s3_class(fit, "regnet")
  expect_equal(fit$para$penalty, "lasso")
  expect_equal(fit$para$lamb.1, out$lambda[1])
})

test_that("default fold construction uses all folds", {
  set.seed(103)
  fold.info = make_foldid(6, folds=5)
  expect_equal(sort(unique(fold.info$foldid)), 1:5)
  expect_false(any(tabulate(fold.info$foldid, nbins=5) == 0))
})

test_that("adjacency mode and non-contiguous clv are handled", {
  n = 35; p = 6
  x = matrix(rnorm(n*p,0,1), n, p)
  X = scale(data.frame(x1=x[,1], x2=x[,2], x3=x[,3], x4=x[,4], x5=x[,5], x6=x[,6]), scale=TRUE)
  Y = 1 + X[,2] * 2 - X[,5] + rnorm(n)

  A.thresholded = Adjacency(X, type="thresholded")
  A.full = Adjacency(X, type="full")
  expect_equal(unname(diag(A.full)), rep(0, ncol(X)))
  expect_gte(sum(A.full != 0), sum(A.thresholded != 0))

  fit = regnet(X, Y, "c", "l", lamb.1=0.1, clv=c(2,4), initiation="zero")
  expect_named(fit$coeff, c("Intercept", "x2", "x4", "x1", "x3", "x5", "x6"))
})

test_that("rank-deficient unpenalized clv is rejected", {
  set.seed(101)
  n = 24; p = 4
  X = matrix(rnorm(n*p), n, p)
  X[,2] = X[,1]
  Y = rnorm(n)
  foldid = rep(1:3, length.out=n)

  expect_error(
    regnet(X, Y, "c", "l", lamb.1=0.1, clv=c(1,2)),
    "rank deficient"
  )
  expect_error(
    cv.regnet(X, Y, "c", "l", lamb.1=0.1, clv=c(1,2), foldid=foldid),
    "rank deficient"
  )

  Y.surv = data.frame(time=exp(rnorm(n)), status=rep(c(1, 1, 0), length.out=n))
  expect_error(
    regnet(X, Y.surv, "s", "l", lamb.1=0.1, clv=c(1,2), robust=FALSE),
    "rank deficient"
  )
  expect_error(
    cv.regnet(X, Y.surv, "s", "l", lamb.1=0.1, clv=c(1,2), foldid=foldid, robust=FALSE),
    "rank deficient"
  )
})

test_that("convergence controls are validated and passed through", {
  set.seed(102)
  n = 30; p = 5
  X = scale(matrix(rnorm(n*p), n, p), scale=TRUE)
  Y = 1 + X[,2] * 2 + rnorm(n)
  foldid = rep(1:3, length.out=n)

  expect_error(
    regnet(X, Y, "c", "l", lamb.1=0.01, maxit=0),
    "maxit"
  )
  expect_error(
    cv.regnet(X, Y, "c", "l", lamb.1=0.01, foldid=foldid, tol=0),
    "tol"
  )
  fit = expect_warning(
    regnet(X, Y, "c", "l", lamb.1=0.01, maxit=1, tol=.Machine$double.eps),
    NA
  )
  expect_null(fit$convergence)

  fit.debug = expect_warning(
    regnet(X, Y, "c", "l", lamb.1=0.01, maxit=1, tol=.Machine$double.eps, debugging=TRUE),
    NA
  )
  expect_named(fit.debug$convergence, c("converged", "niter", "diff"))
  expect_false(fit.debug$convergence$converged)
  expect_equal(fit.debug$convergence$niter, 1L)
  expect_true(is.numeric(fit.debug$convergence$diff))

  out = expect_warning(
    cv.regnet(X, Y, "c", "l", lamb.1=0.01, foldid=foldid, maxit=1, tol=.Machine$double.eps),
    NA
  )
  expect_equal(out$para$maxit, 1L)
  expect_equal(out$para$tol, .Machine$double.eps)
  fit.from.cv = expect_warning(
    regnet(out, out$lambda[1], debugging=TRUE),
    NA
  )
  expect_named(fit.from.cv$convergence, c("converged", "niter", "diff"))
})
