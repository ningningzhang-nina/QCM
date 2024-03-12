partial_linear_regression<-function(X,Y,Z,h1,h2){
  source("kernel_regression_nw.R")
  gy_hat<-mNW(Z, Z, Y, h1, K = dnorm)
  gx_hat<-mNW(Z, Z, X, h2, K = dnorm)
  error1<-Y-gy_hat
  error2<-X-gx_hat
  beta_hat<-(sum(error2^2))^(-1)*sum(error2*error1)
  return(beta_hat)
}