rm(list=ls())
library(spdep)
library(starma)
library(plm)
library(tseries)
library(forecast)
library(FinTS)
library(psych)
library(ggplot2)
library(igraph)
w1<- read.gal("C:/Users/thinkpad/Desktop/starma/starma1.gal",override.id=TRUE)
w1[[17]]<- 26L
w2<- read.gal("C:/Users/thinkpad/Desktop/starma/starma2.gal",override.id=TRUE)
w2[[17]]<- as.integer(c("13","14","19","21"))
tlist <- list(order0=diag(27),
              order1=nb2mat(w1),
              order2=nb2mat(w2))

dsp<-read.csv('C:/Users/thinkpad/Desktop/starma/starma.csv',stringsAsFactors = FALSE)
date<-dsp[,1]
date<-as.Date(date)
dsp2<-dsp[,2:28]
describe(dsp2)
purtest(object=dsp2,exo="intercept",test="ips",lags="AIC",pmax=4)
par(mfrow=c(2,2))
for (i in 1:4){
  plot(date,dsp2[,i],ylab=i,col=2,type="l")}

#
b <- function(a) return(diff(log(a)))
psd <- apply(dsp2, 2, b)
purtest(object=psd,exo="none",test="levinlin",lags="AIC",pmax=4)
purtest(object=psd,exo="trend",test="levinlin",lags="AIC",pmax=4)
purtest(object=psd,exo="intercept",test="levinlin",lags="AIC",pmax=4)
psd1=psd[1:51,]
psd2=psd[52:62,]
plot(date[1:62],psd[,1],type="l")



#step1
par(mfrow=c(1,1))
stacf(psd1, tlist)
stpacf(psd1, tlist)
plot(date[1:51],psd1[,1],type="l")

# stcor.test can be applied to a space-time series to test the null hypothesis of non correlation
stcor.test(psd1, tlist, tlag=NULL, slag=NULL, fitdf=0)
# Compute covariance between 2-nd and 1-st order neighbours, at time lag 1
stcov(psd1, tlist,slag1=2, slag2=1, tlag=1)

max_ar <- 3
max_ma <- 3
model_df <- expand.grid(ar = 0:max_ar, ma = 0:max_ma,iterate = 1)[-1, ]
model_df$bic <- NA
for (j in 1:nrow(model_df)) {
  print (j)
  model_df$bic[j] <- starma(psd1, tlist, model_df$ar[j], model_df$ma[j])$bic
  #model_df$bic[j] <- fit_model$bic
}

opt_ar <- model_df$ar[which.min(model_df$bic)]
opt_ma <- model_df$ma[which.min(model_df$bic)]
starma_model <- starma(psd1, tlist, opt_ar, opt_ma)
starma_model 
summary(starma_model)

est<- psd1-starma_model$residuals
plot(date[1:51],psd[1:51,1],col=2,ylab = '浙江省',type="l")
lines(date[1:51],est[1:51,1],col=3,type="l")

ar1 <- matrix(c(1, 1, 1,  
               1, 1, 1,
               1, 1, 1), 3, 3)
ma1 <- matrix(c(0, 0, 0,
               0, 0, 0,
               0, 0, 0), 3, 3)
model1 <- starma(psd1, tlist, ar1, ma1,iterate = 10)
model1
summary(model1)
est1<- psd1-model1$residuals
 


r211<-array()
r212<-array()
mse1 <- array()
mse11 <- array()
armapre <- matrix(nr=62,nc=27)
for (i in 1:27){
                armafit1 <- auto.arima(psd1[,i])
                armaest1 <- psd1[,i]-armafit1$residuals
                armapre[1:51,i] <- armaest1
                armapre[52:62,i]<-forecast(armafit1,h=11)$mean
                r212[i]=1-sum(armafit1$residuals*armafit1$residuals)/sum((psd1[,i]-mean(psd1[,i]))*(psd1[,i]-mean(psd1[,i])))
                r211[i]=1-sum(model1$residuals[,i]*model1$residuals[,i])/sum((psd1[,i]-mean(psd1[,i]))*(psd1[,i]-mean(psd1[,i])))
                mse1[i]=sum(model1$residuals[,i]*model1$residuals[,i])/32
                mse11[i]=sum(armafit1$residuals*armafit1$residuals)/32
                if(mse1[i]<mse11[i]){print (i)}
                plot(date[1:51],psd1[,i],ylab=i,col=2,type="l")
                lines(date[1:51],est1[,i],col=3,type="l")
                lines(date[1:51],armaest1,col=4,type="l")
                legend("topright",pch=c(15,15,15),legend=c("real","starma","arma"),col=c(2,3,4),bty="n")
}

summse2=sum(model1$residuals*model1$residuals)/30/27
summse21=sum(mse11)/27

stcor.test(model1$residuals, tlist,tlag=3, slag=2, fitdf=-4)
stacf(model1$residuals, tlist)
stpacf(model1$residuals, tlist)
ArchTest(model1$residuals[,1], lags=1, demean = FALSE)

#dynamic forecast[52:62]
psd11 <-rbind(psd1,matrix(nr=11,nc=27))
for(i in 52:62)
{psd11[i,]<- model1$phi[1,] %*% t(tlist[["order0"]] %*% t(rbind(psd11[i-1,],psd11[i-2,],psd11[i-3,])))
             +model1$phi[2,] %*% t(tlist[["order1"]] %*% t(rbind(psd11[i-1,],psd11[i-2,],psd11[i-3,])))
             +model1$phi[3,] %*% t(tlist[["order2"]] %*% t(rbind(psd11[i-1,],psd11[i-2,],psd11[i-3,]))) }


msed1 <- array()
msed11 <- array()
for (i in 1:27){
  plot(date[52:62],psd[52:62,i],ylab=i,col=2,type="l")
  lines(date[52:62],psd11[52:62,i],col=3,type="l")
  msed1[i]=sum((psd11[52:62,i]-psd[52:62,i])*(psd11[52:62,i]-psd[52:62,i]))/11
  lines(date[52:62],armapre[52:62,i],col=4,type="l")
  msed11[i]=sum((armapre[52:62,i]-psd[52:62,i])*(armapre[52:62,i]-psd[52:62,i]))/11
  if(msed1[i]<msed11[i]){print (i)}
  legend("topright",pch=c(15,15,15),legend=c("real","starma","arma"),col=c(2,3,4),bty="n")}

psdx1 <- rbind(est1[1:51,],psd11[52:62,])
for (i in 1:27){
  plot(date[1:62],psd[1:62,i],ylab=i,col=2,type="l")
  lines(date[1:62],psdx1[1:62,i],col=3,type="l")
  lines(date[1:62],armapre[1:62,i],col=4,type="l")
  legend("topright",pch=c(15,15,15),legend=c("real","starma","arma"),col=c(2,3,4),bty="n")}


#static forecast[52:62]
psd21 <- matrix(nr=62,nc=27)
for(i in 52:62)
{psd21[i,]<- model1$phi[1,] %*% t(tlist[["order0"]] %*% t(rbind(psd[i-1,],psd[i-2,],psd[i-3,])))
            +model1$phi[2,] %*% t(tlist[["order1"]] %*% t(rbind(psd[i-1,],psd[i-2,],psd[i-3,])))
            +model1$phi[3,] %*% t(tlist[["order2"]] %*% t(rbind(psd[i-1,],psd[i-2,],psd[i-3,]))) }

for (i in 1:27){
  plot(date[52:62],psd[52:62,i],ylab=i,col=2,type="l")
  lines(date[52:62],psd21[52:62,i],col=3,type="l")
  lines(date[52:62],armapre[52:62,i],col=4,type="l")
  legend("topright",pch=c(15,15,15),legend=c("real","starma","arma"),col=c(2,3,4),bty="n")}


psdx2 <- rbind(est1[1:51,],psd21[52:62,])
for (i in 1:27){
plot(date[1:62],psd[1:62,i],ylab=i,col=2,type="l")
lines(date[1:62],psdx2[1:62,i],col=3,type="l")
lines(date[1:62],armapre[1:62,i],col=4,type="l")
legend("topright",pch=c(15,15,15),legend=c("real","starma","arma"),col=c(2,3,4),bty="n")}


#starmapre end
starmapre1 <- rbind(psd[1,],psdx2[2:62,],matrix(nr=1,nc=27))
starmapre <- matrix(nr=63,nc=27)
for(i in 2:63) {
  for (j in 1:27)
  starmapre[i,j] <- exp(starmapre1[i-1,j]+log(dsp2[i-1,j]))
}
for(i in 1:27){starmapre[1,i] <- dsp2[1,i]}

#armapre end
armapre1 <-rbind(psd[1,],armapre[1:62,])
armapree <- matrix(nr=63,nc=27)
for(i in 2:63) {
  for (j in 1:27)
    armapree[i,j] <- exp(armapre1[i-1,j]+log(dsp2[i-1,j]))
}
for(i in 1:27){armapree[1,i] <- dsp2[1,i]}

for(i in 1:27){
plot(date[1:63],dsp2[1:63,i],ylab=i,col=2,type="l")
lines(date[1:63],starmapre[1:63,i],col=3,type="l")
lines(date[1:63],armapree[1:63,i],col=4,type="l")
legend("topright",pch=c(15,15,15),legend=c("real","starma","arma"),col=c(2,3,4),bty="n")}












































#
par(mfrow=c(1,1))
b <- function(a) return(diff(diff(log(a),12)))
psd <- apply(dsp2, 2, b)
purtest(object=psd,exo="none",test="levinlin",lags="AIC",pmax=4)
purtest(object=psd,exo="trend",test="levinlin",lags="AIC",pmax=4)
purtest(object=psd,exo="intercept",test="levinlin",lags="AIC",pmax=4)
psd1=psd[1:45,]
psd2=psd[46:50,]
plot(date[1:50],psd[,1],type="l")



#step1
stacf(psd1, tlist,10)
stpacf(psd1, tlist,10)
plot(date[1:50],psd1[,1],type="l")

# stcor.test can be applied to a space-time series to test the null hypothesis of non correlation
stcor.test(psd1, tlist, tlag=NULL, slag=NULL, fitdf=0)
# Compute covariance between 2-nd and 1-st order neighbours, at time lag 1
stcov(psd1, tlist,slag1=2, slag2=1, tlag=1)

max_ar <- 3
max_ma <- 3
model_df <- expand.grid(ar = 0:max_ar, ma = 0:max_ma,iterate = 1)[-1, ]
model_df$bic <- NA
for (j in 1:nrow(model_df)) {
  print (j)
  model_df$bic[j] <- starma(psd1, tlist, model_df$ar[j], model_df$ma[j])$bic
  #model_df$bic[j] <- fit_model$bic
}

opt_ar <- model_df$ar[which.min(model_df$bic)]
opt_ma <- model_df$ma[which.min(model_df$bic)]
starma_model <- starma(psd1, tlist, opt_ar, opt_ma)
starma_model 
summary(starma_model)

est<- psd1-starma_model$residuals
plot(date[1:45],psd[1:45,1],col=2,ylab = '浙江省',type="l")
lines(date[1:45],est[1:45,1],col=3,type="l")

ar1 <- matrix(c(1, 0, 0,  
                1, 0, 0,
                1, 0, 0), 3, 3)
ma1 <- matrix(c(0, 0, 0,
                0, 0, 0,
                0, 0, 0), 3, 3)
model1 <- starma(psd1, tlist, ar1, ma1,iterate = 10)
model1
summary(model1)
est1<- psd1-model1$residuals

stcor.test(model1$residuals, tlist, tlag=1, slag=2, fitdf=0)
stacf(model1$residuals, tlist,10)
stpacf(model1$residuals, tlist,10)

r211<-array()
r212<-array()
mse1 <- array()
mse11 <- array()
armapre <- matrix(nr=50,nc=27)
for (i in 1:27){
  armafit1 <- auto.arima(psd1[,i],trace = TRUE)
  armaest1 <- psd1[,i]-armafit1$residuals
  armapre[1:45,i] <- armaest1
  armapre[46:50,i]<-forecast(armafit1,h=5)$mean
  r212[i]=1-sum(armafit1$residuals*armafit1$residuals)/sum((psd1[,i]-mean(psd1[,i]))*(psd1[,i]-mean(psd1[,i])))
  r211[i]=1-sum(model1$residuals[,i]*model1$residuals[,i])/sum((psd1[,i]-mean(psd1[,i]))*(psd1[,i]-mean(psd1[,i])))
  mse1[i]=sum(model1$residuals[,i]*model1$residuals[,i])/32
  mse11[i]=sum(armafit1$residuals*armafit1$residuals)/32
  if(mse1[i]<mse11[i]){print (i)}
  plot(date[1:45],psd1[,i],ylab=i,col=2,type="l")
  lines(date[1:45],est1[,i],col=3,type="l")
  lines(date[1:45],armaest1,col=4,type="l")
  legend("topright",pch=c(15,15,15),legend=c("real","starma","arma"),col=c(2,3,4),bty="n")
}


summse2=sum(model1$residuals*model1$residuals)/30/27
summse21=sum(mse11)/27

stcor.test(model1$residuals, tlist,tlag=1, slag=2, fitdf=0)
stacf(model1$residuals, tlist,10)
stpacf(model1$residuals, tlist,10)
ArchTest(model1$residuals[,1], lags=1, demean = FALSE)

#dynamic forecast[52:62]

psd11 <-rbind(psd[1:45,],matrix(nr=5,nc=27))
for(i in 46:50)
{psd11[i,]<- model1$phi[,1] * tlist[["order0"]] %*% psd11[i-1,]
+model1$phi[,2] * tlist[["order1"]] %*% psd11[i-1,]
+model1$phi[,3] * tlist[["order2"]] %*% psd11[i-1,] }

msed1 <- array()
msed11 <- array()
maed1 <- array()
maed11 <- array()
r2d1 <- array()
r2d11 <- array()
for (i in 1:27){
  plot(date[46:50],psd[46:50,i],ylab=i,col=2,type="l")
  lines(date[46:50],psd11[46:50,i],col=3,type="l")
  msed1[i]=sum((psd11[46:50,i]-psd[46:50,i])*(psd11[46:50,i]-psd[46:50,i]))/5
  maed1[i]=sum(abs(psd11[46:50,i]-psd[46:50,i]))/5
  r2d1[i]=1-sum((psd11[46:50,i]-psd[46:50,i])*(psd11[46:50,i]-psd[46:50,i]))/sum((psd[46:50,i]-mean(psd[46:50,i]))*(psd[46:50,i]-mean(psd[46:50,i])))
  lines(date[46:50],armapre[46:50,i],col=4,type="l")
  msed11[i]=sum((armapre[46:50,i]-psd[46:50,i])*(armapre[46:50,i]-psd[46:50,i]))/5
  maed11[i]=sum(abs(armapre[46:50,i]-psd[46:50,i]))/5
  r2d11[i]=1-sum((armapre[46:50,i]-psd[46:50,i])*(armapre[46:50,i]-psd[46:50,i]))/sum((psd[46:50,i]-mean(psd[46:50,i]))*(psd[46:50,i]-mean(psd[46:50,i])))
  if(msed1[i]<msed11[i]){print (i)}
  legend("topright",pch=c(15,15,15),legend=c("real","starma","arma"),col=c(2,3,4),bty="n")}

sum(msed1)/27
sum(msed11)/27
sum(maed1)/27
sum(maed11)/27


psdx1 <- rbind(est1[1:45,],psd11[46:50,])
for (i in 1:27){
  plot(date[1:50],psd[1:50,i],ylab=i,col=2,type="l")
  lines(date[1:50],psdx1[1:50,i],col=3,type="l")
  lines(date[1:50],armapre[1:50,i],col=4,type="l")
  legend("topright",pch=c(15,15,15),legend=c("real","starma","arma"),col=c(2,3,4),bty="n")}


#static forecast[45:50]
psd21 <- -rbind(psd[1:45,],matrix(nr=5,nc=27))
for(i in 46:50)
  {psd21[i,]<- model1$phi[,1] * tlist[["order0"]] %*% psd[i-1,]
  +model1$phi[,2] * tlist[["order1"]] %*% psd[i-1,]
  +model1$phi[,3] * tlist[["order2"]] %*% psd[i-1,] }

  

for (i in 1:27){
  plot(date[46:50],psd[46:50,i],ylab=i,col=2,type="l")
  lines(date[46:50],psd21[46:50,i],col=3,type="l")
  lines(date[46:50],armapre[46:50,i],col=4,type="l")
  legend("topright",pch=c(15,15,15),legend=c("real","starma","arma"),col=c(2,3,4),bty="n")}


psdx2 <- rbind(est1[1:45,],psd21[46:50,])
for (i in 1:27){
  plot(date[1:50],psd[1:50,i],ylab=i,col=2,type="l")
  lines(date[1:50],psdx2[1:50,i],col=3,type="l")
  lines(date[1:50],armapre[1:50,i],col=4,type="l")
  legend("topright",pch=c(15,15,15),legend=c("real","starma","arma"),col=c(2,3,4),bty="n")}


#starmapre end
pb <- function(a) return(diff(log(a),12))
psdd <- apply(dsp2, 2, pb)
starmapre1 <- rbind(psdd[1:2,],psdx2[2:50,])
tt <- matrix(nr=51,nc=27)
colnames(tt) <- colnames(dsp2)
starmapre2 <- rbind(dsp2[1:12,],tt)
#first difference
#train
for(i in 3:46) {
  for (j in 1:27){
    starmapre1[i,j] <- psdd[i-1,j]+psdx1[i-1,j]
  }}
#test
for(i in 47:51) {
  for (j in 1:27){
    starmapre1[i,j] <- starmapre1[i-1,j]+psdx1[i-1,j]
  }}

#seceond difference
#train
for(i in 13:58){
  for(j in 1:27){
    starmapre2[i,j] <- exp(log(dsp2[i-12,j])+starmapre1[i-12,j])
  }}
#test
for(i in 59:63){
  for(j in 1:27){
    starmapre2[i,j] <- exp(log(starmapre2[i-12,j])+starmapre1[i-12,j])
  }}
plot(date[1:63],dsp2[1:63,5],ylab=i,col=2,type="l")
lines(date[1:63],starmapre2[1:63,5],col=3,type="l")


#armapre end
armapre1 <-rbind(psdd[1,],armapre[1:50,])
armapre2 <- rbind(dsp2[1:12,],tt)
#first difference
#train
for(i in 2:46) {
  for (j in 1:27){
    armapre1[i,j] <- psdd[i-1,j]+armapre[i-1,j]
  }}
#test
for(i in 47:51) {
  for (j in 1:27){
    armapre1[i,j] <- armapre1[i-1,j]+armapre[i-1,j]
  }}

#seceond difference
#train
for(i in 13:58){
  for(j in 1:27){
    armapre2[i,j] <- exp(log(dsp2[i-12,j])+armapre1[i-12,j])
  }}
#test
for(i in 59:63){
  for(j in 1:27){
    armapre2[i,j] <- exp(log(armapre2[i-12,j])+armapre1[i-12,j])
  }}


for(i in 1:27){
  plot(date[1:63],dsp2[1:63,i],ylab=i,col=2,type="l")
  lines(date[1:63],starmapre2[1:63,i],col=3,type="l")
  lines(date[1:63],armapre2[1:63,i],col=4,type="l")
  legend("topright",pch=c(15,15,15),legend=c("real","starma","arma"),col=c(2,3,4),bty="n")}


  par(mfrow=c(2,2))
  plot(date[1:63],dsp2[1:63,2],ylab="云南",col=2,type="l")
  lines(date[1:63],starmapre2[1:63,2],col=3,type="l")
  lines(date[1:63],armapre2[1:63,2],col=4,type="l")
  legend("topright",pch=c(15,15,15),legend=c("real","starma","arma"),col=c(2,3,4),bty="n")

  plot(date[1:63],dsp2[1:63,3],ylab="新疆",col=2,type="l")
  lines(date[1:63],starmapre2[1:63,3],col=3,type="l")
  lines(date[1:63],armapre2[1:63,3],col=4,type="l")
  legend("topright",pch=c(15,15,15),legend=c("real","starma","arma"),col=c(2,3,4),bty="n")
  
  plot(date[1:63],dsp2[1:63,6],ylab="陕西",col=2,type="l")
  lines(date[1:63],starmapre2[1:63,6],col=3,type="l")
  lines(date[1:63],armapre2[1:63,6],col=4,type="l")
  legend("topright",pch=c(15,15,15),legend=c("real","starma","arma"),col=c(2,3,4),bty="n")
  
  
  plot(date[1:63],dsp2[1:63,7],ylab="山西",col=2,type="l")
  lines(date[1:63],starmapre2[1:63,7],col=3,type="l")
  lines(date[1:63],armapre2[1:63,7],col=4,type="l")
  legend("topright",pch=c(15,15,15),legend=c("real","starma","arma"),col=c(2,3,4),bty="n")
  
  plot(date[1:63],dsp2[1:63,8],ylab="山东",col=2,type="l")
  lines(date[1:63],starmapre2[1:63,8],col=3,type="l")
  lines(date[1:63],armapre2[1:63,8],col=4,type="l")
  legend("topright",pch=c(15,15,15),legend=c("real","starma","arma"),col=c(2,3,4),bty="n")
  
  plot(date[1:63],dsp2[1:63,9],ylab="青海",col=2,type="l")
  lines(date[1:63],starmapre2[1:63,9],col=3,type="l")
  lines(date[1:63],armapre2[1:63,9],col=4,type="l")
  legend("topright",pch=c(15,15,15),legend=c("real","starma","arma"),col=c(2,3,4),bty="n")
  
  plot(date[1:63],dsp2[1:63,12],ylab="辽宁",col=2,type="l")
  lines(date[1:63],starmapre2[1:63,12],col=3,type="l")
  lines(date[1:63],armapre2[1:63,12],col=4,type="l")
  legend("topright",pch=c(15,15,15),legend=c("real","starma","arma"),col=c(2,3,4),bty="n")
  
  plot(date[1:63],dsp2[1:63,13],ylab="江西",col=2,type="l")
  lines(date[1:63],starmapre2[1:63,13],col=3,type="l")
  lines(date[1:63],armapre2[1:63,13],col=4,type="l")
  legend("topright",pch=c(15,15,15),legend=c("real","starma","arma"),col=c(2,3,4),bty="n")
  
  
  plot(date[1:63],dsp2[1:63,14],ylab="湖南",col=2,type="l")
  lines(date[1:63],starmapre2[1:63,14],col=3,type="l")
  lines(date[1:63],armapre2[1:63,14],col=4,type="l")
  legend("topright",pch=c(15,15,15),legend=c("real","starma","arma"),col=c(2,3,4),bty="n")
  
  plot(date[1:63],dsp2[1:63,15],ylab="湖北",col=2,type="l")
  lines(date[1:63],starmapre2[1:63,15],col=3,type="l")
  lines(date[1:63],armapre2[1:63,15],col=4,type="l")
  legend("topright",pch=c(15,15,15),legend=c("real","starma","arma"),col=c(2,3,4),bty="n")
  
  plot(date[1:63],dsp2[1:63,16],ylab="北京",col=2,type="l")
  lines(date[1:63],starmapre2[1:63,16],col=3,type="l")
  lines(date[1:63],armapre2[1:63,16],col=4,type="l")
  legend("topright",pch=c(15,15,15),legend=c("real","starma","arma"),col=c(2,3,4),bty="n")
  
  plot(date[1:63],dsp2[1:63,17],ylab="海南",col=2,type="l")
  lines(date[1:63],starmapre2[1:63,17],col=3,type="l")
  lines(date[1:63],armapre2[1:63,17],col=4,type="l")
  legend("topright",pch=c(15,15,15),legend=c("real","starma","arma"),col=c(2,3,4),bty="n")
  
  plot(date[1:63],dsp2[1:63,18],ylab="贵州",col=2,type="l")
  lines(date[1:63],starmapre2[1:63,18],col=3,type="l")
  lines(date[1:63],armapre2[1:63,18],col=4,type="l")
  legend("topright",pch=c(15,15,15),legend=c("real","starma","arma"),col=c(2,3,4),bty="n")
  
  plot(date[1:63],dsp2[1:63,19],ylab="广西",col=2,type="l")
  lines(date[1:63],starmapre2[1:63,19],col=3,type="l")
  lines(date[1:63],armapre2[1:63,19],col=4,type="l")
  legend("topright",pch=c(15,15,15),legend=c("real","starma","arma"),col=c(2,3,4),bty="n")
  
  plot(date[1:63],dsp2[1:63,22],ylab="安徽",col=2,type="l")
  lines(date[1:63],starmapre2[1:63,22],col=3,type="l")
  lines(date[1:63],armapre2[1:63,22],col=4,type="l")
  legend("topright",pch=c(15,15,15),legend=c("real","starma","arma"),col=c(2,3,4),bty="n")
  
  plot(date[1:63],dsp2[1:63,24],ylab="重庆",col=2,type="l")
  lines(date[1:63],starmapre2[1:63,24],col=3,type="l")
  lines(date[1:63],armapre2[1:63,24],col=4,type="l")
  legend("topright",pch=c(15,15,15),legend=c("real","starma","arma"),col=c(2,3,4),bty="n")
  


