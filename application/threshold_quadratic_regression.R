threshold_quadratic_regression<-function(X,Y,Z){
  lower_limit=quantile(Z,0.05)
  upper_limit=quantile(Z,0.95)
  change_points=seq(lower_limit,upper_limit,by=0.005)
  mse=c()
  
  for (i in (1:length(change_points))){
    newZ=Z+change_points[i]
    Z_square<-newZ^2
    Z_square_negative<-newZ^2*ifelse(newZ<0,1,0)
    if (is.null(X)){
      
      model=lm(Y~Z_square+Z_square_negative)
    } else{
    model=lm(Y~X+Z_square+Z_square_negative)
    }
    mse[i]=mean((model$residuals)^2)
  }
  #plot(change_points,mse)
  ind=which.min(mse)
  #points(change_points[ind],mse[ind],col=2)
  
  opt_change=change_points[ind]
  #print(opt_change)
  newZ=(Z+opt_change)
  Z_square<-newZ^2
  Z_square_negative<-newZ^2*ifelse(newZ<0,1,0)
  if (is.null(X)) {
    model=lm(Y~Z_square+Z_square_negative)
  } else {
    model=lm(Y~X+Z_square+Z_square_negative)
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