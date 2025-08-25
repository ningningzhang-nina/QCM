library(R.matlab)
library(latex2exp)
library("GARCHSK")
setwd("/Users/gina/Downloads/Research/2.QCM/Code files/final_final_application/")

setEPS()
postscript("new_qcm.eps",width=18,height=4)
par(mfrow=c(3,3),mar=c(1.5,2.7,2,1.2),mgp = c(1.5, 0.4, 0),font.lab=1.5, font.main=1.5, cex.main=1.5,cex.lab=1.5,cex.axis=1.2)
indexnames<-list('AUD/USD','NZD/USD','CAD/USD')
moments <- list(TeX(r'(Estimators of $h_t$)'), TeX(r'(Estimators of $s_t$)'), TeX(r'(Estimators of $k_t$)'))
AUDUSD <- readMat("final_AUD_USD_moments.mat")
yyyy <- readxl::read_excel("./data/AUD_USD.xlsx")
yyy <- yyyy[[2]]
yy <- rev(yyy)
y <- (log(yy[-1]) - log(yy[-length(yy)])) * 100
res <- read.csv("AUD_USD_res.csv",header=FALSE)  # 默认UTF-8编码

Est<-garchsk_est(res$V1)
results<-garchsk_construct(Est$params,res$V1)
muhat<-results$mu
sigma2hat<-results$h
skewhat<-results$sk
kurthat<-results$ku
results$y <- y
write.csv(results, 'garchsk_AUD_USD_results.csv')

n<-length(AUDUSD$moments[1,])
plot(AUDUSD$moments[1,(n-861):(n-774)]^2,ylab=indexnames[1],xlab="",main=moments[1],type="l",lwd=1.5, xaxt = "n", ylim = range(c(0, 2)))
par(new = TRUE)
plot(sigma2hat[(n-861):(n-774)],type="l",lty=4,lwd=2,axes = FALSE, xlab = "", ylab = "")
axis(4)  # 右侧y轴
#mtext("y2", side = 4, line = 3)
abline(v=c(49,57),lty=2)
axis(side=1,at=c(1,49,57,88),labels =c("Jan. 1","Mar. 9","19","May 1"))
plot(AUDUSD$moments[2,(n-861):(n-774)],main=moments[2],xlab="",ylab="",type="l",lwd=1.5, xaxt = "n", ylim = range(c(-0.6, 0.1)))
par(new = TRUE)
plot(skewhat[(n-861):(n-774)],type="l",lty=4,lwd=2,axes = FALSE, xlab = "", ylab = "")
axis(4)  # 右侧y轴
#mtext("y2", side = 4, line = 3)
abline(v=c(49,57),lty=2)
axis(side=1,at=c(1,49,57,88),labels =c("Jan. 1","Mar. 9","19","May 1"))
plot(AUDUSD$moments[3,(n-861):(n-774)],main=moments[3],xlab="",ylab="",type="l",lwd=1.5, xaxt = "n")
par(new = TRUE)
plot(kurthat[(n-861):(n-774)],type="l",lty=4,lwd=2,axes = FALSE, xlab = "", ylab = "")
axis(4)  # 右侧y轴
mtext("y2", side = 4, line = 3)
abline(v=c(49,57),lty=2)
axis(side=1,at=c(1,49,57,88),labels =c("Jan. 1","Mar. 9","19","May 1"))




NZDUSD <- readMat("final_NZD_USD_moments.mat")
yyyy <- readxl::read_excel("./data/NZD_USD.xlsx")
yyy <- yyyy[[2]]
yy <- rev(yyy)
y <- (log(yy[-1]) - log(yy[-length(yy)])) * 100
res <- read.csv("NZD_USD_res.csv",header=FALSE)  # 默认UTF-8编码

Est<-garchsk_est(res$V1)
results<-garchsk_construct(Est$params,res$V1)
muhat<-results$mu
sigma2hat<-results$h
skewhat<-results$sk
kurthat<-results$ku
results$y <- y
write.csv(results, 'garchsk_NZD_USD_results.csv')

plot(NZDUSD$moments[1,(n-861):(n-774)]^2, ylim = c(0, 1.2),xlab="",ylab=indexnames[2],type="l",lwd=1.5, xaxt = "n")
lines(sigma2hat[(n-861):(n-774)],lty=4,lwd=2)
abline(v=c(49,58),lty=2)
axis(side=1,at=c(1,49,58,88),labels =c("Jan. 1","Mar. 9","20","May 1"))
plot(NZDUSD$moments[2,(n-861):(n-774)],xlab="",ylab="",type="l",lwd=1.5, xaxt = "n")
lines(skewhat[(n-861):(n-774)],lty=4,lwd=2)
abline(v=c(49,58),lty=2)
axis(side=1,at=c(1,49,58,88),labels =c("Jan. 1","Mar. 9","20","May 1"))
plot(NZDUSD$moments[3,(n-861):(n-774)],xlab="",ylab="",type="l",lwd=1.5, xaxt = "n")
par(new = TRUE)
plot(kurthat[(n-861):(n-774)],type="l",lty=4,lwd=2,axes = FALSE, xlab = "", ylab = "")
axis(4)  # 右侧y轴
mtext("y2", side = 4, line = 3)
abline(v=c(49,58),lty=2)
axis(side=1,at=c(1,49,58,88),labels =c("Jan. 1","Mar. 9","20","May 1"))





USDCAD <- readMat("final_CAD_USD_moments.mat")
yyyy <- readxl::read_excel("./data/USD_CAD.xlsx")
yyy <- yyyy[[2]]
yy <- rev(yyy)
y <- (log(yy[-1]) - log(yy[-length(yy)])) * 100
res <- read.csv("CAD_USD_res.csv",header=FALSE)  # 默认UTF-8编码

Est<-garchsk_est(res$V1)
results<-garchsk_construct(Est$params,res$V1)
muhat<-results$mu
sigma2hat<-results$h
skewhat<-results$sk
kurthat<-results$ku
results$y <- y
write.csv(results, 'garchsk_CAD_USD_results.csv')


m<-length(USDCAD$moments[1,])
plot(USDCAD$moments[1,(n-861):(n-774)]^2,xlab="",ylab=indexnames[3],type="l",lwd=1.5, xaxt = "n", ylim = range(c(0, 1)))
lines(sigma2hat[(n-861):(n-774)],lty = 4,lwd=2)
abline(v=c(49,57),lty=2)
axis(side=1,at=c(1,49,57,88),labels =c("Jan. 1","Mar. 9","19","May 1"))


plot(USDCAD$moments[2,(n-861):(n-774)], ylim = c(-0.4, 0.2),xlab="",ylab="",type="l",lwd=1.5, xaxt = "n")
lines(skewhat[(n-861):(n-774)],lty = 4,lwd=2)
abline(v=c(49,57),lty=2)
axis(side=1,at=c(1,49,57,88),labels =c("Jan. 1","Mar. 9","19","May 1"))


plot(USDCAD$moments[3,(n-861):(n-774)],xlab="",ylab="",type="l",lwd=1.5, xaxt = "n")
par(new = TRUE)
plot(kurthat[(n-861):(n-774)],type="l",lty=4,lwd=2,axes = FALSE, xlab = "", ylab = "")
axis(4)  # 右侧y轴
mtext("y2", side = 4, line = 3)
abline(v=c(49,57),lty=2)
axis(side=1,at=c(1,49,57,88),labels =c("Jan. 1","Mar. 9","19","May 1"))

dev.off()

