threshold_reverse_quadratic_regression<-function(X,Y,Z){
  lower_limit=quantile(Z,0.05)
  upper_limit=quantile(Z,0.95)
  change_points=seq(0,upper_limit,by=0.001)
  print(change_points)
  mse=c()
  
  for (i in (1:length(change_points))){
    #newZ=Z+change_points[i]
    ZZ=Z
    Z_square_negative<-ZZ^2*ifelse((ZZ<change_points[i])&(ZZ>-change_points[i]),-1,1)
    #Z_square1<-ZZ^2*ifelse((ZZ<-change_points[i])|(ZZ>change_points[i]),1,0)
    #ZZ=Z
    #Z_square2<-ZZ*ifelse((ZZ>change_points[i]),1,0)
    if (is.null(X)){
      
      model=lm(Y~Z_square_negative)
    } else{
      model=lm(Y~X+Z_square_negative)
    }
    mse[i]=mean((model$residuals)^2)
  }
  #plot(change_points,mse)
  ind=which.min(mse)
  #points(change_points[ind],mse[ind],col=2)
  
  opt_change=change_points[ind]
  #print(opt_change)
  #newZ=(Z+opt_change)
  ZZ=Z
  
  Z_square_negative<-ZZ^2*ifelse((ZZ<opt_change)&(ZZ>-opt_change),-1,1)
  #Z_square1<-ZZ^2*ifelse((ZZ<-opt_change)|(ZZ>opt_change),1,0)
  #ZZ=Z
  #Z_square2<-ZZ*ifelse((ZZ>opt_change),1,0)
 
  if (is.null(X)) {
    model=lm(Y~Z_square_negative)
  } else {
    model=lm(Y~X+Z_square_negative)
  }
  model_summary=summary(model)
  coeff=model_summary$coefficients[,c(1,2,4)]
  adjusted_rsquare=summary(model)$r.squared
  lbqtest1=Box.test(model$resid, lag = 1, type = "Ljung")$p.value
  lbqtest2=Box.test((model$resid)^2, lag = 1, type = "Ljung")$p.value
  fittedvalues=model$fitted.values
  results=list(fittedvalues,coeff,adjusted_rsquare,lbqtest1,lbqtest2,opt_change)
  return(results)
}