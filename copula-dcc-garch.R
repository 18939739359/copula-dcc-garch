##
library(fGarch)
library(rmgarch)
library(rugarch)
library(tseries)
library(zoo)
library(copula)
library(VineCopula)
library(rugarch)
library(xts)
library(psych)
library(ggplot2)
library(FinTS)
library(forecast)
library(swfscMisc)
library(mapdata)
library(dplyr)
setwd("F:/r/tgarchcovar")
options(max.print=1000000)

price<-read.csv('debugdata.csv',stringsAsFactors = FALSE)
date<-price[,1]
date<-as.Date(date)
priceret<-price[,2:18]
LogReturn <- function(vx) return(100*diff(log(vx)))
logret <- apply(priceret, 2, LogReturn)
logret <- logret[1:1240,]
rets<-as.xts(logret,date[1:1240])
date<-date[1:1240]



garch11.spec <- ugarchspec(mean.model = list(armaOrder = c(0,0)), 
                           variance.model = list(garchOrder = c(1,1), 
                                                 model = "sGARCH"), 
                           distribution.model = "sstd")


dcc.garch11.spec = dccspec(uspec = multispec(replicate(2, garch11.spec)), 
                           dccOrder = c(1,1),  distribution = "mvt")


cgarchspe=cgarchspec(uspec= multispec(replicate(2, garch11.spec)),lag.max = NULL,lag.criterion = c("AIC"),
                     dccOrder = c(1, 1), asymmetric = TRUE,
                     distribution.model = list(copula = c("mvt"),
                                               method = c("Kendall"), time.varying = TRUE,
                                               transformation = c("spd")))


##garch
mean <- matrix(nr=1,nc=17)
for(i in 2:27)
{
  garch.fit.i <- ugarchfit(garch11.spec,data=rets[,i],solver="solnp")
  mean[1,i]<- mean(rets[,i]-residuals(garch.fit.i))
}

##dcc garch
dcccovar <- rets
dccVaRidown <-rets
dccVaRiupside<-rets
dccVaRi50<-rets
dcccovari<-rets
dcc <- rets
deltadcccovari <- rets
dcc.res<-rets
for (i in 2:17)
{
  data= cbind(rets[,1],rets[,i])
  names(data) <- c("1","i")
  dcc.fit.i <- dccfit(dcc.garch11.spec, data = data, solver="solnp")
  dcc.fit.i
  dcc[,i] <- rcor(dcc.fit.i)[1,2,]
  dcc.res[,i] <- dcc.fit.i@model[["residuals"]][,2]
  dcc.res[,1] <- dcc.fit.i@model[["residuals"]][,1]
  dccVaRidown[,i] = mean[1,i]+ dcc.fit.i@model[["sigma"]][,2]*qdist("sstd", p=0.05, mu = 0, sigma = 1,
                                                                    skew = dcc.fit.i@mfit[["coef"]][["[i].skew"]], 
                                                                    shape= dcc.fit.i@mfit[["coef"]][["[i].shape"]])
  
  dccVaRiupside[,i] = mean[1,i]+ dcc.fit.i@model[["sigma"]][,2]*qdist("sstd", p=0.95, mu = 0, sigma = 1,
                                                                      skew = dcc.fit.i@mfit[["coef"]][["[i].skew"]], 
                                                                      shape= dcc.fit.i@mfit[["coef"]][["[i].shape"]])
  
  dccVaRi50[,i] = mean[1,i]+ dcc.fit.i@model[["sigma"]][,2]*qdist("sstd", p=0.50, mu = 0, sigma = 1,
                                                                  skew = dcc.fit.i@mfit[["coef"]][["[i].skew"]], 
                                                                  shape= dcc.fit.i@mfit[["coef"]][["[i].shape"]])
  deltadcccovari[,i]=dcc[,i]*dcc.fit.i@model[["sigma"]][,1]/dcc.fit.i@model[["sigma"]][,2]*(dccVaRidown[,i]-dccVaRi50[,i])
  dcccovari[,i]=dcc[,i]*dcc.fit.i@model[["sigma"]][,1]/dcc.fit.i@model[["sigma"]][,2]*dccVaRidown[,i]
  #plot (date,dcc[,i],col=2,type="l")
  dcc.resdd <- pobs(pnorm(dcc.res))
  print(BiCopSelect(dcc.resdd[,1] ,dcc.resdd[,i], familyset = NA,selectioncrit = "AIC"))
  y <- BiCopSelect(dcc.resdd[,1] ,dcc.resdd[,i], familyset = NA,selectioncrit = "AIC")
  print(y[["AIC"]])
}

dcc.resdd <- pobs(pnorm(dcc.res))
write.csv(dcc.resdd,file = "d:/dcc.csv")

##copula dcc garch

copuladcc <- rets
copulacovar <- rets
VaRidown <-rets
VaRiupside<-rets
VaRi50<-rets
covari<-rets
deltacovari<-rets
dcc.sigma <-rets
for (i in 2 :17)
{data= cbind(rets[,1],rets[,i])
names(data) <- c("1","i")
dcccopula.fit.i<- cgarchfit(cgarchspe, data=data, parallel = parallel, parallel.control = parallel.control, 
                            fit.control = list(eval.se=TRUE))
copuladcc[,i] <- rcor(dcccopula.fit.i)[1,2,]
dcc.sigma[,i] <- dcccopula.fit.i@model[["sigma"]][,2]
VaRidown[,i] = mean[1,i]+ dcccopula.fit.i@model[["sigma"]][,2]*qdist("sstd", p=0.05, mu = 0, sigma = 1,
                                                                     skew = dcccopula.fit.i@mfit[["coef"]][["[i].skew"]], 
                                                                     shape= dcccopula.fit.i@mfit[["coef"]][["[i].shape"]])

VaRiupside[,i] = mean[1,i]+ dcccopula.fit.i@model[["sigma"]][,2]*qdist("sstd", p=0.95, mu = 0, sigma = 1,
                                                                       skew = dcccopula.fit.i@mfit[["coef"]][["[i].skew"]], 
                                                                       shape= dcccopula.fit.i@mfit[["coef"]][["[i].shape"]])

VaRi50[,i] = mean[1,i]+ dcccopula.fit.i@model[["sigma"]][,2]*qdist("sstd", p=0.50, mu = 0, sigma = 1,
                                                                   skew = dcccopula.fit.i@mfit[["coef"]][["[i].skew"]], 
                                                                   shape= dcccopula.fit.i@mfit[["coef"]][["[i].shape"]])
deltacovari[,i]=copuladcc[,i]*dcccopula.fit.i@model[["sigma"]][,1]/dcccopula.fit.i@model[["sigma"]][,2]*(VaRidown[,i]-VaRi50[,i])
covari[,i]=copuladcc[,i]*dcccopula.fit.i@model[["sigma"]][,1]/dcccopula.fit.i@model[["sigma"]][,2]*VaRidown[,i]
print(i)
#show(dcccopula.fit.i)
}


par(mfrow=c(2,2))
for (i in 2 :17)
{
  plot (date,dcc[,i],ylab=i,col=2,type="l",ylim=c(-0.1,1.1))
  lines(date,copuladcc[,i],col=3,type="l")
  legend("topright",pch=c(15,15),legend=c("dcc","copuladcc"),col=c(2,3),bty="n",cex=0.6)
  
}

for (i in 2 :17){
  print(mean(dcc[,i]))
  print(mean(copuladcc[,i]))
}

VaR.back.test=function(xdata,VaRdata,p,con.level){
  p=p;
  con.level=con.level;
  prob=1-p;
  N=sum(xdata>VaRdata);
  T=length(xdata);
  LR_POF=-2*log((prob^N)*(1-prob)^(T-N))+2*log(((N/T)^N)*(1-N/T)^(T-N));
  critical=qchisq(con.level,df=1);
  P_value=pchisq(LR_POF,df=1,lower.tail=F);
  list(LR_POF,P_value)
}

for (i in 2:17){
  dccvartesti=VaR.back.test(rets[,i],dccVaRidown[,i],0.05,0.95)
  vartesti=VaR.back.test(rets[,i],VaRidown[,i],0.05,0.95)
  print(i)
  options(digits=22)
  print(vartesti)
  print(dccvartesti)
}
for (i in 2:17){
  plot (date,deltacovari[,i],ylab=i,col=2,type="l",ylim=c(-7,7))
  lines(date,deltadcccovari[,i],col=3,type="l")
  legend("topright",pch=c(15,15),legend=c("copuladcc","dcc"),col=c(2,3),bty="n",cex=0.8)
}

for (i in 2:17){
  dccvartesti=VaR.back.test(rets[,i],dcccovari[,i],0.05,0.95)
  vartesti=VaR.back.test(rets[,i],covari[,i],0.05,0.95)
  print(i)
  options(digits=4)
  print(vartesti)
  print(dccvartesti)
}

for (i in 2:17){
  plot (date,dcccovari[,i],ylab=2,col=2,type="l")
  lines(date,covari[,i],col=3,type="l")
  lines(date,rets[,i],col=4,type="l")
  #lines(date,rets[,i],col=4,type="l")
  legend("topright",pch=c(15,15),legend=c("copuladcc","dcc"),col=c(2,3),bty="n",cex=0.8)
}


x=dcc.sigma[,3]

p= density(x)
plot(p)
describe(x)


for (i in 2:17)
{
  data= cbind(rets[,1],rets[,i])
  names(data) <- c("1","i")
  dcc.fit.i <- dccfit(dcc.garch11.spec, data = data, solver="solnp")
  dcc.fit.i
  dcc[,i] <- rcor(dcc.fit.i)[1,2,]
  dcc.res[,i] <- dcc.fit.i@model[["residuals"]][,2]
  dcc.res[,1] <- dcc.fit.i@model[["residuals"]][,1]
  dcc.resdd <- pobs(pnorm(dcc.res))
  print(i)
  print(BiCopSelect(dcc.resdd[,1] ,dcc.resdd[,i], familyset = 2,selectioncrit = "AIC"))
  y <- BiCopSelect(dcc.resdd[,1] ,dcc.resdd[,i], familyset = 2,selectioncrit = "AIC")
  print(y[["AIC"]])
}
write.csv(dcc.sigma,file = "d:/dccsigma.csv")
for (i in 2:17)
{
  plot(density(dcc.sigma[,i]))
  
}
describe(dcc.sigma)
ks.test(dcc.sigma[,2],pt,df=4)