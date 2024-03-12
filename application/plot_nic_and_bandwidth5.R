library(R.matlab)
library(pracma)
#library(ggplot2)
library(latex2exp)
#library(np)
setwd("/Users/gina/Downloads/Research/2.QCM/Code files/final_final_application/")
source("kernel_regression_nw.R")
source("partial_linear_regression.R")
source("threshold_quadratic_regression.R")
source("threshold_reverse_quadratic_regression.R")
indexnames<-list('AUD_USD','NZD_USD','CAD_USD')
coef <- matrix(nrow = 3, ncol = 3)
adjusted_rsquare <- matrix(nrow = 3, ncol = 6)
b=readMat("nic_bandwidth.mat")
b1=b$b1
b2=b$b2
b3=b$b3
# b1=matrix(,nrow=8,ncol=3)
# b2=matrix(,nrow=8,ncol=3)
# b3=matrix(,nrow=8,ncol=3)
setEPS()
postscript("NIC.eps",width=6,height=6)
par(mfrow=c(3,3),mar=c(1.5,2.5,2,0.1),mgp = c(1.5, 0.4, 0),font.lab=1.5, font.main=1.5)
new_indexnames<-list('AUD/USD','NZD/USD','CAD/USD')
#moments_names <- list(TeX(r'($g_h(\cdot)$)'), TeX(r'($g_s(\cdot)$)'), TeX(r'($g_k(\cdot)$)'))
moments_names <- list(TeX("NICs for $h_t$"), TeX("NICs for $s_t$"), TeX("NICs for $k_t$"))

for (i in 1:3) {
  filename1 <- file.path(paste0("final", "_", indexnames[i], "_moments.mat"))
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
  
  Z1_hat<-seq(min(Z1),max(Z2),by=0.01)
  
  
  
  # b11 <- bw.cv.grid(X = Z1, Y = Y1, plot.cv = TRUE)
  # b21 <- bw.cv.grid(X = Z1, Y = X1, plot.cv = TRUE)
  # b1[i,1]=b11
  # b2[i,1]=b21
  
  beta_hat<-partial_linear_regression(X1,Y1,Z1,b1[i,1],b2[i,1])
  beta_hat
  Y11 <- Y1-beta_hat*X1
  # b31 <- bw.cv.grid(X = Z1, Y = Y11, plot.cv = TRUE)
  # b3[i,1]=b31
  Y11_hat<-mNW(Z1, Z1, Y11, b3[i,1], K = dnorm)
  # gjr
  Z1_square<-Z1^2
  Z1_square_negative<-Z1^2*ifelse(Z1<0,1,0)
  model1=lm(Y11~Z1_square+Z1_square_negative)
  print(summary(model1))
  #coeff1=model11$coefficients[,1]
  para.results1=model1$fitted.values
  adjusted_rsquare[i,1]=1-var(Y11-para.results1)/var(Y11)*(length(Y11)-1)/(length(Y11)-4)
  para1=cbind(Z1,para.results1)
  para1=sortrows(para1, k = 1)
  # garch
  model11=lm(Y11~Z1_square)
  print(summary(model11))
  #coeff1=model11$coefficients[,1]
  para.results11=model11$fitted.values
  adjusted_rsquare[i,2]=1-var(Y11-para.results11)/var(Y11)*(length(Y11)-1)/(length(Y11)-3)
  para11=cbind(Z1,para.results11)
  para11=sortrows(para11, k = 1)
  # kernel
  kernel_results1=cbind(Z1,Y11_hat)
  kernel_results1=sortrows(kernel_results1)
  
  if (i==1) {
    plot(kernel_results1[,1],kernel_results1[,2],type="l",lty=2,ylab=new_indexnames[i],main=moments_names[1])
  } else {
    plot(kernel_results1[,1],kernel_results1[,2],type="l",lty=2,ylab=new_indexnames[i],main="")
  } 
  
  lines(para1[,1],para1[,2], lty=1)
  lines(para11[,1],para11[,2], lty=3)
  
  
  
  ### skewness
  Y2 <- skew0[-1]
  X2 <- skew0[-n]
  Z2 <- eta[-n]
  
  # b12 <- bw.cv.grid(X = Z2, Y = Y2, plot.cv = TRUE)
  # b22 <- bw.cv.grid(X = Z2, Y = X2, plot.cv = TRUE)
  # b1[i,2]=b12
  # b2[i,2]=b22
  beta_hat<-partial_linear_regression(X2,Y2,Z2,b1[i,2],b2[i,2])
  
  Y22 <- Y2-beta_hat*X2
  # b32 <- bw.cv.grid(X = Z2, Y = Y22, plot.cv = TRUE)
  # b3[i,2]=b32
  Y22_hat<-mNW(Z2, Z2, Y22, b3[i,2], K = dnorm)
  # our skew nic
  Z3<-Z2*ifelse(Z2<0,1,0)
  #model2=lm(Y22~Z2+Z3)
  # if (i==1){
  #   model2=lm(Y22~poly(Z2,9))
  #   print(summary(model2))
  #   #coeff1=model11$coefficients[,1]
  #   para.results2=model2$fitted.values
  #   adjusted_rsquare[i,3]=1-var(Y22-para.results2)/var(Y22)*(length(Y22)-1)/(length(Y22)-11)
  # } else if (i==2) {
  #   model2=lm(Y22~poly(Z2,10))
  #   print(summary(model2))
  #   #coeff1=model11$coefficients[,1]
  #   para.results2=model2$fitted.values
  #   adjusted_rsquare[i,3]=1-var(Y22-para.results2)/var(Y22)*(length(Y22)-1)/(length(Y22)-12)
  # } else {
  #   model2=lm(Y22~poly(Z2,7))
  #   print(summary(model2))
  #   #coeff1=model11$coefficients[,1]
  #   para.results2=model2$fitted.values
  #   adjusted_rsquare[i,3]=1-var(Y22-para.results2)/var(Y22)*(length(Y22)-1)/(length(Y22)-9)
  # }
  #   
  # 
  # 
  # para2=cbind(Z2,para.results2)
  # para2=sortrows(para2, k = 1)
  # previous skew nic
  Z2_triple<-Z2^3
  model22=lm(Y22~Z2_triple)
  print(summary(model22))
  #coeff1=model11$coefficients[,1]
  para.results22=model22$fitted.values
  adjusted_rsquare[i,4]=1-var(Y22-para.results22)/var(Y22)*(length(Y22)-1)/(length(Y22)-3)
  para22=cbind(Z2,para.results22)
  para22=sortrows(para22, k = 1)
  # kernel nic
  kernel_results2=cbind(Z2,Y22_hat)
  kernel_results2=sortrows(kernel_results2)
  if (i==1){
    plot(kernel_results2[,1],kernel_results2[,2],type="l",lty=2,ylab="",main=moments_names[2])
  }else {
    plot(kernel_results2[,1],kernel_results2[,2],type="l",lty=2,ylab="",main="")
  }
  
  #lines(para2[,1],para2[,2], lty=1)
  lines(para22[,1],para22[,2], lty=3)
  
  
  ###
  # if (i==1){
  #   kurt0=kurt0+7.0748-3.7057
  # } else if (i==2){
  #   kurt0=kurt0+7.1818-3.5129
  # }else{
  #   kurt0=kurt0+5.3841-3.6003
  # }
  Y3 <- kurt0[-1]
  X3 <- kurt0[-n]
  Z3 <- eta[-n]
  
  # b13 <- bw.cv.grid(X = Z3, Y = Y3, plot.cv = TRUE)
  # b23 <- bw.cv.grid(X = Z3, Y = X3, plot.cv = TRUE)
  # b1[i,3]=b13
  # b2[i,3]=b23
  beta_hat<-partial_linear_regression(X3,Y3,Z3,b1[i,3],b2[i,3])
  
  
  Y33 <- Y3-beta_hat*X3
  # b33 <- bw.cv.grid(X = Z3, Y = Y33, plot.cv = TRUE)
  # b3[i,3]=b33
  Y33_hat<-mNW(Z3, Z3, Y33, b3[i,3], K = dnorm)
  # our kurt nic
  #fit<-threshold_reverse_quadratic_regression(NULL,Y33,Z3)
  #para.results3=fit[[1]]
  # if (i==1){
  #   model3=lm(Y33~poly(Z3,9))
  # } else if (i==2){
  #   model3=lm(Y33~poly(Z3,8))
  # } else {
  #   model3=lm(Y33~poly(Z3,9))
  # }
  # 
  # print(summary(model3))
  # #coeff1=model11$coefficients[,1]
  # para.results3=model3$fitted.values
  # 
  # 
  # para3=cbind(Z3,para.results3)
  # para3=sortrows(para3, k = 1)
  # if (i==1){
  #   model3=lm(Y33~poly(Z3,4))
  # } else if (i==2) {
  #   model3=lm(Y33~poly(Z3,4))
  # } else {
  #   model3=lm(Y33~poly(Z3,4))
  # }
  # 
  # 
  # summary(model3)
  #coeff1=model11$coefficients[,1]
  #para.results3=model3$fitted.values
  #para3=cbind(Z3,para.results3)
  #para3=sortrows(para3, k = 1)
  #adjusted_rsquare[i,5]=1-var(Y33-para.results3)/var(Y33)*(length(Y33)-1)/(length(Y33)-12)
  # previous nic
  Z3_four<-Z3^4
  
  model33=lm(Y33~Z3_four)
  print(summary(model33))
  #coeff1=model11$coefficients[,1]
  para.results33=model33$fitted.values
  adjusted_rsquare[i,6]=1-var(Y33-para.results33)/var(Y33)*(length(Y33)-1)/(length(Y33)-3)
  para33=cbind(Z3,para.results33)
  para33=sortrows(para33, k = 1)
  # kernel 
  kernel_results3=cbind(Z3,Y33_hat)
  kernel_results3=sortrows(kernel_results3)
  if (i==1){
    plot(kernel_results3[,1],kernel_results3[,2],type="l",lty=2,ylab="",main=moments_names[3])
  } else {
    plot(kernel_results3[,1],kernel_results3[,2],type="l",lty=2,ylab="",main="")
  }
  
  #lines(para3[,1],para3[,2], lty=1)
  lines(para33[,1],para33[,2], lty=3)
}
#writeMat("nic_bandwidth.mat",b1=b1,b2=b2,b3=b3)
dev.off()