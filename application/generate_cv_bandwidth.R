# `https://www.stat.cmu.edu/~cshalizi/350/hw/06/np/html/np.plregression.html`1
library(R.matlab)
library(ggplot2)
library(latex2exp)
library(np)
setwd("/Users/gina/Downloads/projects/Project_implied_conditional_moments/ningning/codes/!final_application/")
indexnames<-list('USD_CAD','NZD_USD','AUD_USD')
coef <- matrix(nrow = 3, ncol = 3)
b1 <- matrix(nrow = 3, ncol = 3)
b2 <- matrix(nrow = 3, ncol = 3)
b3 <- matrix(nrow = 3, ncol = 3)
for (i in 1:3) {
  start_time <- Sys.time()
  filename1 <- file.path(paste0("final", "_","final", "_", indexnames[i], "_moments.mat"))
  moments <- readMat(filename1)
  sigma0 <- moments$moments[1,]
  h_t <- sigma0^2
  skew0 <- moments$moments[2,]
  kurt0 <- moments$moments[3,]
  
  epsilon <- moments$res;
  eta <- epsilon/sigma0
  n <- length(h_t)
  ###
  Y1 <- h_t[-1]
  X1 <- h_t[-n]
  Z1 <- epsilon[-n]
  
  b11 <- bw.cv.grid(X = Z1, Y = Y1, plot.cv = TRUE)
  b21 <- bw.cv.grid(X = Z1, Y = X1, plot.cv = TRUE)
  start_time <- Sys.time()
  
  end_time <- Sys.time()
  end_time - start_time
  data<-cbind(Y1,X1,Z1)
  #
  bw1 <- npplregbw(xdat=X1, zdat=Z1, ydat=Y1, tol=.1, ftol=.1)
  pl1 <- npplreg(bws=bw1)
  summary(pl1)
  
  
  Y11 <- Y1-pl1$xcoef[[1]]*X1
  b31 <- bw.cv.grid(X = Z1, Y = Y11, plot.cv = TRUE)
  
  bw1 <- npregbw(xdat=Z1,ydat=Y11)
  model1 <- npreg(bws = bw1, gradients = TRUE)
  summary(model1)
  # 
  # simdata<-data.frame(x=Z1,y=Y11)
  # h <- seq(from=0,to=1,length=10)
  # find_hcv(simdata, nw, h)
  # lscv(h,simdata, nw)
  # np.gcv(data = simdata, h.seq=NULL, num.h = 50, estimator = "NW", 
  #        kernel = "gaussian")
  
  
  # data$Y11=Y11
  # data$X1=X1
  # mdl<-nls(Y11~a0+a1*X1^2,data,start = list(a0 = 0.1, a1 = 0.2))
  # theta_hat=coef(mdl)
  # test_nonpara_para(X1,Y11,length(X1),model1$bws[[1]],1,theta_hat,100)
  
  
  ### skewness
  Y2 <- skew0[-1]
  X2 <- skew0[-n]
  Z2 <- eta[-n]
  
  b12 <- bw.cv.grid(X = Z2, Y = Y2, plot.cv = TRUE)
  b22 <- bw.cv.grid(X = Z2, Y = X2, plot.cv = TRUE)
  data<-cbind(Y2,X2,Z2)
  
  
  bw2 <- npplregbw(xdat=X2, zdat=Z2, ydat=Y2, tol=.1, ftol=.1)
  pl2 <- npplreg(bws=bw2)
  summary(pl2)
  
  Y22 <- Y2-pl2$xcoef[[1]]*X2
  b32 <- bw.cv.grid(X = Z2, Y = Y22, plot.cv = TRUE)
  
  bw2 <- npregbw(xdat=Z2,ydat=Y22)
  model2 <- npreg(bws = bw2, gradients = TRUE)
  summary(model2)
  
  ###
  Y3 <- kurt0[-1]
  X3 <- kurt0[-n]
  Z3 <- eta[-n]
  
  b13 <- bw.cv.grid(X = Z3, Y = Y3, plot.cv = TRUE)
  b23 <- bw.cv.grid(X = Z3, Y = X3, plot.cv = TRUE)
  
  data<-cbind(Y3,X3,Z3)
  
  
  bw3 <- npplregbw(xdat=X3, zdat=Z3, ydat=Y3, tol=.1, ftol=.1)
  pl3 <- npplreg(bws=bw3)
  summary(pl3)
  
  Y33 <- Y3-pl3$xcoef[[1]]*X3
  b33 <- bw.cv.grid(X = Z3, Y = Y33, plot.cv = TRUE)
  bw3 <- npregbw(xdat=Z3,ydat=Y33)
  model3 <- npreg(bws = bw3, gradients = TRUE)
  summary(model3)
  
  
  coef[i,]<-cbind(pl1$xcoef[[1]],pl2$xcoef[[1]],pl3$xcoef[[1]])
  # 
  b1[i,]<-cbind(pl1$bw$bandwidth[[1]],pl2$bw$bandwidth[[1]],pl3$bw$bandwidth[[1]])
  
  b2[i,]<-cbind(pl1$bw$bandwidth[[2]],pl2$bw$bandwidth[[2]],pl3$bw$bandwidth[[2]])
  
  b3[i,]<-cbind(model1$bws[[1]],model2$bws[[1]],model3$bws[[1]])
  b1[i,]<-cbind(b11,b12,b13)
  b2[i,]<-cbind(b21,b22,b23)
  b3[i,]<-cbind(b31,b32,b33)
  end_time <- Sys.time()
  end_time - start_time
}

writeMat("nic_bandwidth.mat",coef=coef,b1=b1,b2=b2,b3=b3)



test_nonpara_para<-function (X,Y,n,h,model,theta_hat,nboot) {
  # reference: W. Hardle and E.Mammen 'comparing nonparametric versus parametric regression fits'
  yhat=function_theta(X,theta_hat,1)
  T_n <- n*(h^(0.5))*integrate(f_hat,min(X),max(X),X=X,h=h,Y=Y,yhat=yhat)
  Y_fit<-m_h(X,X,h,Y)
  Y_par <- function_theta(X,theta_hat,type)
  # start to do wild boostrap
  epsilons = Y - m_h(X,X,h,Y)
  WB = Wild_Boot(epsilons, nboot)
  T_n_star = matrix(0,nrow=1,ncol=nboot)
  for(j in 1:nboot){
    y_star = Y_par + WB[,j]
    # reconstruct parametric regression using x and y_star
    # parametric model
    # para_model22 = @(b,x)b(1)+b(2)^2*x.^2+b(3)^2*(x<0).*x.^2;
    beta0 <- theta_hat
    data$Y=Y
    data$X=X
    if (model==1){
      mdl2 = nls(y_star~a0+a1*X^2,start=beta0,data,start=list(a0=theta_hat[1],a1=theta_hat[2]))
    }else if (model==2){
      mdl2 = nls(y_star~a0+a1*X^3,start=beta0,data,start=list(a0=theta_hat[1],a1=theta_hat[2]))
    }else if (model==3){
      mdl2 = nls(y_star~a0+a1*X^4,start=beta0,data,start=list(a0=theta_hat[1],a1=theta_hat[2]))
    }else if (model==4){
      mdl2 = nls(y_star~a0+a1*X^2+a2*X^2*ifelse(X<0,1,0),data,start=list(a0=theta_hat[1],a1=theta_hat[2],a2=theta_hat[3]))
    }else if (model==5){
      mdl2 = nls(y_star~a0+a1*X+a2*X*ifelse(X<0,1,0),data,start=list(a0=theta_hat[1],a1=theta_hat[2],a2=theta_hat[3]))
    }else {
      mdl2 = nls(y_star~a0+a1*(X+a3)^2+a2*(X+a3)^2*ifelse(X+a3<0,1,0),data,start=list(a0=theta_hat[1],a1=theta_hat[2],a2=theta_hat[3],a4=theta_hat[4]))
    }
  }
  theta_hat<-coef(mdl2)
  T_n_star[1,j] = n*h^(1/2)*integrate(f_hat,min(X),max(X),X=X,h=h,Y=y_star,theta_hat=theta_hat,model=1)
  T_n
  T_n_star[1,j]
}

kerf_gau <- function(z,h){
  kerf<-exp(-z*z/2)/sqrt(2*pi)/h; # gaussian kernel function
  return(kerf)
}

f_hat <- function(x,X,h,Y,theta_hat,model){
  f.hat<-(m_h(x,X,h,Y)-m_h(x,X,h,function_theta(X,theta_hat,model)))^2
  return(f.hat)
}

f_hat <- function(x,X,h,Y,yhat){
  numerator1=0
  numerator2=0
  denominator=0
  for (i in 1:length(X)){
    numerator1=numerator1+m.hat(x,X[i],h,Y[i])
    numerator2=numerator2+m.hat(x,X[i],h,yhat[i])
    denominator=denominator+m.hat(x,X[i],h,1)
  }
  f.hat<-(numerator1/denominator-numerator2/denominator)^2
  return(f.hat)
}
m.hat<-function(x,xx,h,y){
  m_hat=kerf_gau((x-X)/h,h)*y
  return(m_hat)
}

m_h<-function(x,X,h,Y){
  if (length(x)==1) {
    m.h = sum(kerf_gau((x-X)/h,h)*Y)/sum(kerf_gau((x-X)/h,h))
  } else {
    di<-(t(replicate(length(X), X1))-replicate(length(X), X))/h
    m.h = rowSums(kerf_gau(di,h)*Y)/rowSums(kerf_gau(di,h))
  }
  return(m.h)
}


Wild_Boot<-function(x,nboot){
  n<-length(x)
  p<-matrix(rbinom(n*nboot,1,(sqrt(5)+1)/(2*sqrt(5))),nrow=n)
  z<- (1 - sqrt(5))/2 * p + (1 + sqrt(5))/2 * (1 - p)
  WB<- z*replicate(nboot, x)
  return(WB)
}

function_theta<-function(X,beta0,model){
  if (model==1){
    function.theta<-beta0[1]+beta0[2]*X^2
  }else if(model==2){
    function.theta<-beta0[1]+beta0[2]*X^3
  }else if (model==3){
    function.theta<-beta0[1]+beta0[2]*X^4
  }else if (model==4){
    function.theta<-beta0[1]+beta0[2]*X^2+beta0[3]*X^2*ifelse(X<0,1,0)
  }else if (model==5){
    function.theta<-beta0[1]+beta0[2]*X^3+beta0[3]*X^3*ifelse(X<0,1,0)
  }else {
    function.theta<-beta0[1]+beta0[2]*(X+beta0[4])^2+beta0[3]*(X+beta0[4])^2*ifelse(X+beta0[4]<0,1,0)
  }
  return(function.theta)
}
my<-function(x){
  res<-sum((rep(x,3)-c(1,2,3))^2)+2^2
  return(res)
}

my<-function(x,y){
  
  
  res<-((x-y[1])^2+(x-y[2])^2)/((x+y[1])^2)
  #res<-sum((x-y)^2)
  return(res)
}

testFn0 <- function(x) {
  prod(cos(x))
}

adaptIntegrate(testFn0, rep(0,2), rep(1,2), tol=1e-4)
# 
# 
# 
# 
# 
mNW <- function(x, X, Y, h, K = dnorm) {
  
  # Arguments
  # x: evaluation points
  # X: vector (size n) with the predictors
  # Y: vector (size n) with the response variable
  # h: bandwidth
  # K: kernel
  
  # Matrix of size n x length(x)
  Kx <- sapply(X, function(Xi) K((x - Xi) / h) / h)
  
  # Weights
  W <- Kx / rowSums(Kx) # Column recycling!
  
  # Means at x ("drop" to drop the matrix attributes)
  drop(W %*% Y)
  
}
cvNW <- function(X, Y, h, K = dnorm) {
  sum(((Y - mNW(x = X, X = X, Y = Y, h = h, K = K)) /
         (1 - K(0) / colSums(K(outer(X, X, "-") / h))))^2)
}

bw.cv.grid <- function(X, Y,
                       h.grid = diff(range(X)) * (seq(0.1, 1, l = 200))^2,
                       K = dnorm, plot.cv = FALSE) {
  obj <- sapply(h.grid, function(h) cvNW(X = X, Y = Y, h = h, K = K))
  h <- h.grid[which.min(obj)]
  if (plot.cv) {
    plot(h.grid, obj, type = "o")
    rug(h.grid)
    abline(v = h, col = 2, lwd = 2)
  }
  h
}