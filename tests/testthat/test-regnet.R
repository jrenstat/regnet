test_that("regnet_surv", {
  n = 25; p = 5;
  x = matrix(rnorm(n*p,0,5), n, p)
  X = scale(data.frame(x,X6=x[,4]+x[,5]*0.5,X7=x[,4]*0.2-x[,5]), scale=TRUE); #Adjacency(X)
  Y0= exp(2 + X[,4:7]%*%c(3,2,3,-2) + rnorm(n)); Y1 = sample(rep(c(0,1,1,1),n/2),n)
  Y = data.frame(time=(Y0+Y0*(Y1-1)*runif(n,0.2,0.5)), status=Y1)

  out = cv.regnet(X, Y, "s", "n")
  fit = regnet(X, Y, "s", "n", out$lambda[1,1], out$lambda[1,2]); fit$coeff
  expect_length(fit$coeff, ncol(X)+1)
  expect_named(fit$coeff)
  expect_equal(ncol(fit$Adj), sum(fit$coeff[-1]!=0))
  expect_equal(fit$para$penalty, "network")
  expect_false(fit$para$robust)

  out = cv.regnet(X, Y, "s", "n",robust = TRUE)
  fit = regnet(X, Y, "s", "n", out$lambda[1,1], out$lambda[1,2], robust = TRUE); #fit$coeff
  expect_equal(ncol(fit$Adj), sum(fit$coeff[-1]!=0))
  expect_equal(fit$para$penalty, "network")
  expect_true(fit$para$robust)

  out = cv.regnet(X, Y, "s", "m", folds = n) #LOOCV
  fit = regnet(X, Y, "s", "m", out$lambda[1]); fit$coeff
  expect_equal(ncol(fit$Adj), sum(fit$coeff[-1]!=0))
  expect_equal(fit$para$penalty, "mcp")
  expect_false(fit$para$robust)

  out = cv.regnet(X, Y, "s", "m",robust = TRUE)
  fit = regnet(X, Y, "s", "m", out$lambda[1], robust = TRUE); #fit$coeff
  expect_equal(ncol(fit$Adj), sum(fit$coeff[-1]!=0))
  expect_equal(fit$para$penalty, "mcp")
  expect_true(fit$para$robust)

  out = cv.regnet(X, Y, "s", "l")
  fit = regnet(X, Y, "s", "l", out$lambda[1]); fit$coeff
  expect_equal(ncol(fit$Adj), sum(fit$coeff[-1]!=0))
  expect_equal(fit$para$penalty, "lasso")
  expect_false(fit$para$robust)

  out = cv.regnet(X, Y, "s", "l",robust = TRUE)
  fit = regnet(X, Y, "s", "l", out$lambda[1], robust = TRUE); #fit$coeff
  expect_equal(ncol(fit$Adj), sum(fit$coeff[-1]!=0))
  expect_equal(fit$para$penalty, "lasso")
  expect_true(fit$para$robust)
})


test_that("regnet_data_format_cont", {
  n = 25; p = 5;
  x = matrix(rnorm(n*p,0,5), n, p)
  X = scale(data.frame(x,X6=x[,4]+x[,5]*0.5,X7=x[,4]*0.2-x[,5]), scale=TRUE); #Adjacency(X)
  Y= 3 + X[,4:7]%*%c(3,2,3,-2) + rnorm(n)
  out = cv.regnet(X, Y, "c", "n")
  expect_error(regnet(X[,1], Y, "c", "n",out$lambda[1,1], out$lambda[1,2]), "too less variables for network penalty")
  expect_error(regnet(X, Y[1:4], "c", "n",out$lambda[1,1], out$lambda[1,2]), "length of Y does not match")
  expect_error(regnet(X, Y, "c", "n",out$lambda[1,1], out$lambda[1,2], alpha.i = 2), "alpha.i should be between 0 and 1")
})

test_that("regnet_data_format_surv", {
  n = 25; p = 5;
  x = matrix(rnorm(n*p,0,5), n, p)
  X = scale(data.frame(x,X6=x[,4]+x[,5]*0.5,X7=x[,4]*0.2-x[,5]), scale=TRUE); #Adjacency(X)
  Y0= exp(2 + X[,4:7]%*%c(3,2,3,-2) + rnorm(n)); Y1 = sample(rep(c(0,1,1,1),n/2),n)
  Y = data.frame(time=(Y0+Y0*(Y1-1)*runif(n,0.2,0.5)), status=Y1)
  out = cv.regnet(X, Y, "s", "n")

  expect_error(regnet(X, Y, "s", "n"), "Both lambda1 and lambda2 need to be provided")
  expect_error(regnet(X, Y, "s", "n", out$lambda[1,1]), "Lambda2 needs to be provided")
  expect_error(regnet(X, Y0, "s", "n",out$lambda[1,1], out$lambda[1,2]), "Y should be a two-column matrix")
  expect_error(regnet(X, cbind(Y[,1], Y[,2]), "s", "n",out$lambda[1,1], out$lambda[1,2]), "columns named 'time' and 'status'")
  expect_error(regnet(X, data.frame(time=Y[,1]-10, status=Y[,2]), "s", "n",out$lambda[1,1], out$lambda[1,2]), "survival times need to be positive")
  expect_error(regnet(X, data.frame(time=Y[,1], status=Y[,2]*2), "s", "n",out$lambda[1,1], out$lambda[1,2]), "binary variable of 1 and 0")
  expect_error(regnet(X, Y[1:4,], "s", "n",out$lambda[1,1], out$lambda[1,2]), "the number of rows of Y does not match")
})

test_that("regnet_cont", {
  n = 15; p = 5;
  x = matrix(rnorm(n*p,0,5), n, p)
  X = scale(data.frame(x,X6=x[,4]+x[,5]*0.5,X7=x[,4]*0.2-x[,5]), scale=TRUE); #Adjacency(X)
  Y= 3 + X[,4:7]%*%c(3,2,3,-2) + rnorm(n)

  out = cv.regnet(X, Y, "c", "n", folds = n) #LOOCV
  fit = regnet(X, Y, "c", "n", out$lambda[1,1], out$lambda[1,2]); fit$coeff
  expect_length(fit$coeff, ncol(X)+1)
  expect_named(fit$coeff)
  expect_equal(ncol(fit$Adj), sum(fit$coeff[-1]!=0))
  expect_equal(fit$para$penalty, "network")
  expect_false(fit$para$robust)

  out = cv.regnet(X, Y, "c", "n",robust = TRUE)
  fit = regnet(X, Y, "c", "n", out$lambda[1,1], out$lambda[1,2], robust = TRUE)
  expect_equal(ncol(fit$Adj), sum(fit$coeff[-1]!=0))
  expect_equal(fit$para$penalty, "network")
  expect_true(fit$para$robust)

  out = cv.regnet(X, Y, "c", "m")
  fit = regnet(X, Y, "c", "m", out$lambda[1]);
  expect_equal(ncol(fit$Adj), sum(fit$coeff[-1]!=0))
  expect_equal(fit$para$penalty, "mcp")
  expect_false(fit$para$robust)

  out = cv.regnet(X, Y, "c", "m",robust = TRUE, folds = n) #LOOCV
  fit = regnet(X, Y, "c", "m", out$lambda[1], robust = TRUE); #fit$coeff
  expect_equal(ncol(fit$Adj), sum(fit$coeff[-1]!=0))
  expect_equal(fit$para$penalty, "mcp")
  expect_true(fit$para$robust)

  out = cv.regnet(X, Y, "c", "l")
  fit = regnet(X, Y, "c", "l", out$lambda[1]); fit$coeff
  expect_equal(ncol(fit$Adj), sum(fit$coeff[-1]!=0))
  expect_equal(fit$para$penalty, "lasso")
  expect_false(fit$para$robust)

  out = cv.regnet(X, Y, "c", "l",robust = TRUE)
  fit = regnet(X, Y, "c", "l", out$lambda[1], robust = TRUE); #fit$coeff
  expect_equal(ncol(fit$Adj), sum(fit$coeff[-1]!=0))
  expect_equal(fit$para$penalty, "lasso")
  expect_true(fit$para$robust)
})


test_that("regnet_logit", {
  n = 40; p = 5; Y = rep(0,n)
  x = matrix(rnorm(n*p,0,5), n, p)
  X = scale(data.frame(x,X6=x[,4]+x[,5]*0.5,X7=x[,4]*0.2-x[,5]), scale=TRUE); #Adjacency(X)
  Y0 = 1/(1+exp(-X[,4:7]%*%c(3,2,3,-2)))
  Y = rbinom(n,1,Y0); #hist(Y)

  out = cv.regnet(X, Y, "b", "n")
  fit = regnet(X, Y, "b", "n", out$lambda[1,1], out$lambda[1,2]); fit$coeff
  expect_length(fit$coeff, ncol(X)+1)
  expect_named(fit$coeff)
  expect_equal(ncol(fit$Adj), sum(fit$coeff[-1]!=0))
  expect_equal(fit$para$penalty, "network")
  expect_false(fit$para$robust)
  expect_message(regnet(X, Y, "b", "n", out$lambda[1,1], out$lambda[1,2], robust = TRUE),"robust methods are not available")
  fit.robust = suppressMessages(regnet(X, Y, "b", "n", out$lambda[1,1], out$lambda[1,2], robust = TRUE))
  expect_false(fit.robust$para$robust)


  out = cv.regnet(X, Y, "b", "m")
  fit = regnet(X, Y, "b", "m", out$lambda[1]); fit$coeff
  expect_equal(ncol(fit$Adj), sum(fit$coeff[-1]!=0))
  expect_equal(fit$para$penalty, "mcp")
  expect_false(fit$para$robust)


  out = cv.regnet(X, Y, "b", "l", folds = n) #LOOCV
  fit = regnet(X, Y, "b", "l", out$lambda[1]); fit$coeff
  expect_equal(ncol(fit$Adj), sum(fit$coeff[-1]!=0))
  expect_equal(fit$para$penalty, "lasso")
  expect_false(fit$para$robust)

})
