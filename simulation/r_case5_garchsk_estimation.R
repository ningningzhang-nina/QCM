library("R.matlab")
library("GARCHSK")
library("MASS")
setwd("/Users/gina/Downloads/Research/2.QCM/Code files/final_simulation/results/")
#######CASE 1
#xdata<-readMat("new_simulation_garch.mat")
xdata<-readMat("new_DGP1_garch.mat")
muhat<-matrix(0,100,1000)
sigma2hat<-matrix(0,100,1000)
skewhat<-matrix(0,100,1000)
kurthat<-matrix(0,100,1000)
for (rep in 1:100){
  Est<-garchsk_est(xdata$yy[rep,])
  results<-garchsk_construct(Est$params,xdata$yy[rep,])
  muhat[rep,]<-results$mu
  sigma2hat[rep,]<-results$h
  skewhat[rep,]<-results$sk
  kurthat[rep,]<-results$ku
}
write.matrix(muhat,file="new_dgp1_case5_muhat.csv")
write.matrix(sigma2hat,file="new_dgp1_case5_sigma2hat.csv")
write.matrix(skewhat,file="new_dgp1_case5_skewhat.csv")
write.matrix(kurthat,file="new_dgp1_case5_kurthat.csv")

#######CASE 2
#xdata<-readMat("final_simulation_studentt_timevaring_original_noar.mat")
xdata<-readMat("new_DGP2_studentt.mat")
muhat<-matrix(0,100,1000)
sigma2hat<-matrix(0,100,1000)
skewhat<-matrix(0,100,1000)
kurthat<-matrix(0,100,1000)
for (rep in 1:100){
  Est<-garchsk_est(xdata$y[rep,])
  results<-garchsk_construct(Est$params,xdata$y[rep,])
  muhat[rep,]<-results$mu
  sigma2hat[rep,]<-results$h
  skewhat[rep,]<-results$sk
  kurthat[rep,]<-results$ku
}
write.matrix(muhat,file="new_dgp2_case5_muhat.csv")
write.matrix(sigma2hat,file="new_dgp2_case5_sigma2hat.csv")
write.matrix(skewhat,file="new_dgp2_case5_skewhat.csv")
write.matrix(kurthat,file="new_dgp2_case5_kurthat.csv")
#######CASE 3
#xdata<-readMat("simulation_garch_snp_caviar.mat")
xdata<-readMat("new_DGP3_garch_snp.mat")
muhat<-matrix(0,100,1000)
sigma2hat<-matrix(0,100,1000)
skewhat<-matrix(0,100,1000)
kurthat<-matrix(0,100,1000)
for (rep in 1:100){
  Est<-garchsk_est(xdata$yy[rep,])
  results<-garchsk_construct(Est$params,xdata$yy[rep,])
  muhat[rep,]<-results$mu
  sigma2hat[rep,]<-results$h
  skewhat[rep,]<-results$sk
  kurthat[rep,]<-results$ku
}
write.matrix(muhat,file="new_dgp3_case5_muhat.csv")
write.matrix(sigma2hat,file="new_dgp3_case5_sigma2hat.csv")
write.matrix(skewhat,file="new_dgp3_case5_skewhat.csv")
write.matrix(kurthat,file="new_dgp3_case5_kurthat.csv")
#######CASE 4
xdata<-readMat("simulation_gjr_skew_studentt_caviar.mat")
muhat<-matrix(0,100,1000)
sigma2hat<-matrix(0,100,1000)
skewhat<-matrix(0,100,1000)
kurthat<-matrix(0,100,1000)
for (rep in 1:100){
  Est<-garchsk_est(xdata$yyy[rep,])
  results<-garchsk_construct(Est$params,xdata$yyy[rep,])
  muhat[rep,]<-results$mu
  sigma2hat[rep,]<-results$h
  skewhat[rep,]<-results$sk
  kurthat[rep,]<-results$ku
}
write.matrix(muhat,file="dgp4_case5_muhat.csv")
write.matrix(sigma2hat,file="dgp4_case5_sigma2hat.csv")
write.matrix(skewhat,file="dgp4_case5_skewhat.csv")
write.matrix(kurthat,file="dgp4_case5_kurthat.csv")
#######CASE 4
xdata<-readMat("new_dpg4_gir_skew_studentt.mat")
muhat<-matrix(0,100,1000)
sigma2hat<-matrix(0,100,1000)
skewhat<-matrix(0,100,1000)
kurthat<-matrix(0,100,1000)
for (rep in 1:100){
  Est<-garchsk_est(xdata$yyy[rep,])
  results<-garchsk_construct(Est$params,xdata$yyy[rep,])
  muhat[rep,]<-results$mu
  sigma2hat[rep,]<-results$h
  skewhat[rep,]<-results$sk
  kurthat[rep,]<-results$ku
}
write.matrix(muhat,file="new_dgp4_case5_muhat.csv")
write.matrix(sigma2hat,file="new_dgp4_case5_sigma2hat.csv")
write.matrix(skewhat,file="new_dgp4_case5_skewhat.csv")
write.matrix(kurthat,file="new_dgp4_case5_kurthat.csv")
#######CASE 5
xdata<-readMat("new_simulation_engle.mat")
muhat<-matrix(0,100,1000)
sigma2hat<-matrix(0,100,1000)
skewhat<-matrix(0,100,1000)
kurthat<-matrix(0,100,1000)
for (rep in 1:100){
  Est<-garchsk_est(xdata$yy[rep,])
  results<-garchsk_construct(Est$params,xdata$yy[rep,])
  muhat[rep,]<-results$mu
  sigma2hat[rep,]<-results$h
  skewhat[rep,]<-results$sk
  kurthat[rep,]<-results$ku
}
write.matrix(muhat,file="dgp5_case5_muhat.csv")
write.matrix(sigma2hat,file="dgp5_case5_sigma2hat.csv")
write.matrix(skewhat,file="dgp5_case5_skewhat.csv")
write.matrix(kurthat,file="dgp5_case5_kurthat.csv")
