#Pris en main codes de David pour faire une analyse fonctionelle des contours de SIE - Feb 2021
setwd("~/Documents/Research/David_Nerini/FDA_SIE")
library(fda)
require(shapes)
library(stR)
library(GENEAread)
library(fftwtools)
library(swdft)
library(Matrix)
library(rootSolve)
library(princurve)
library(fields)
library(oce)

source("OPAfd.r")
source("GPAfd.r")
source("utilGPAfd.r")
source("fct.pour.poly-root.r")

load("ice_sheet_GEBCO.RData")
load("Sea_Ice_Extent_1978-11_2020-12.RData")

### EXTRACTION DATA
n=length(Extent)
vec=cbind(Extent[[35]]$x,Extent[[35]]$y)
nvec=nrow(vec)
vec=rbind(vec[620:nvec,],vec[1:619,])
Extent[[35]]$x=vec[,1]
Extent[[35]]$y=vec[,2]

vec=cbind(Extent[[36]]$x,Extent[[36]]$y)
nvec=nrow(vec)
vec=rbind(vec[620:nvec,],vec[1:619,])
Extent[[36]]$x=vec[,1]
Extent[[36]]$y=vec[,2]

### QQS GRAPHES
limit=c(-41e05,41e05)
X11()
par(pty="s")
isheet=rbind(isheet,isheet[1,])

for(k in 1:n){
  plot(Extent[[k]],xlim=limit,ylim=limit,type="l",main=date[k])
  polygon(isheet,col = "dark blue")
  text(locator(1))
}

##### ON PLONGE DANS UNE BASE DE FOURIER COMPLEXE
X=complex(re=isheet[,1],im=isheet[,2])

m=30
period=2*pi
time<-seq(0,period,length=length(X))
Cbasecomplex=eval.complex.basis(time,m,period)
Bmat=t(Conj(Cbasecomplex))%*%Cbasecomplex
Dmat=t(Conj(Cbasecomplex))%*%X
coefs=solve(Bmat)%*%Dmat
Cbasecomplex=eval.complex.basis(seq(0,period,length=10000),m,period)
Xfd=Cbasecomplex%*%coefs

#PLOT
x11()
par(pty="s")
rge=c(-5e06,5e06)
plot(X,col="grey",xlim=rge,ylim=rge,type="l")
abline(h=0,v=0,lty=3)
lines(Xfd,col=c(1,2))

################ FIT DE TS LE MONDE
#### 
m=100
period=2*pi
nobs=length(Extent)

Xfd=matrix(NA,2*m+1,nobs)

for(k in 1:nobs){
  cat("k=",k,fill=TRUE)
  npts=nrow(Extent[[k]]$x)
  Xk<-complex(re=Extent[[k]]$x,im=Extent[[k]]$y)
  time<-seq(0,period,length=length(Xk))
  Cbasecomplex=eval.complex.basis(time,m,period)
  Bmat=t(Conj(Cbasecomplex))%*%Cbasecomplex
  Dmat=t(Conj(Cbasecomplex))%*%Xk
  coefs=solve(Bmat)%*%Dmat
  Xfd[,k]<-coefs
}
##### PLOT DATA + FIT
tabs=seq(0,period,length=1000)
Cbasecomplex=eval.complex.basis(tabs,m,period)
Xline=Cbasecomplex%*%Xfd

x11()
for(k in 1:nobs){
  plot(Extent[[k]],xlim=limit,ylim=limit,type="l",main=date[k],col="grey",lwd=2)
  polygon(isheet,col = "dark blue")
  lines(Xline[,k],col=2)
  points(Xfd[101,k],col="white",pch=19)
  text(locator(1))
}

#### CALCUL SURFACE
Deriv=complex(re=0,im=2*pi/period*seq(-m,m,1))
#tabs=seq(0,period,length=100000)
#Cbasecomplex=eval.complex.basis(tabs,m,period)
#dt=diff(seq(0,period,length=100000))[1]
#Snum=dt*sum(Im(Cbasecomplex%*%coefs)*Re(Cbasecomplex%*%(Deriv*coefs)))
surf1=NULL
for(k in 1:nobs){
  Sfd=Im(-period/2*t(Conj(Deriv*Xfd[,k]))%*%(Xfd[,k]))
  surf1=c(surf1,abs(Sfd))
}
### PLOT SURFACE DE GLACE
tdate=paste(rep(1978:2020,each=12),c(paste(0,1:9,sep=""),"10","11","12"),"15",sep="-")
tdate=tdate[11:516]
tdate=as.Date(tdate,origin="1899-12-30")
plot(tdate,surf1/1e12,type="l")

## TRAITEMENT DE LA SERIE TEMPORELLE
##### AUTRE METHODE DE DECOMPOSITION : STR, HYNDMAN
surf=surf1/1e12
vec.ts=ts(surf,freq=12,start=c(1978,11),end=c(2020,12))
vec.str=AutoSTR(vec.ts,conf=c(0.025,0.975), trace=TRUE)
#plot(vec.str,main="ICE AREA")
plot(components(vec.str))
STR = components(vec.str)
#STR = cbind(STR,STR[,1]-STR[,3]) #To run only once
#Remove the december 1987 month
ind = which(tdate=="1987-12-15")
STR[ind,] = apply(STR[c(ind-1,ind+1),],2,mean)
#
yyl = c(expression(paste("SIE raw (10"^"6"," km"^"2",")")), "SIE trend", "SIE season", "SIE residual", "SIE raw-season")
letter = c("a)","b)","c)","d)","e)")
range = c(2000+((10/12)*0)*10^-1,2014+((10/12)*3)*10^-1,2017+((10/12)*5)*10^-1)  #decimal years (the month year is the starting month...)
ttime = seq(1978+((10/12)*10)*10^-1,2020+((10/12)*11)*10^-1,((10/12)*1)*10^-1) + ((10/12)*1)*10^-1/2  #Nov 1979 to Dec 2020 in decimal year plotted in the middle of the month

#
#Trends of subplo #1
n1 = which(tdate=="2000-01-15")
n2 = which(tdate=="2014-04-15")
n3 = which(tdate=="2017-06-15")
n4 = which(tdate=="2020-12-15")
nn = list(); nn[[1]] = 1:n1; nn[[2]] = n1:n2; nn[[3]] = n2:n3; nn[[4]] = n3:n4
trend = NULL

png("STR.png",width=7, height=7, units="in", res=600)
layout(matrix(c(1,2,3,4,5,5), nrow = 6, ncol = 1, byrow = TRUE))
par(mar = c(0,5,1,1))
for(t in 1:5){
  if(t ==5){par(mar = c(3,5,1,1))
    plot(STR[,t],ylab = yyl[t],axes = F,ylim = c(21.5,28),yaxs = "i")
  }else{plot(STR[,t],ylab = yyl[t],axes = F)}
  abline(v = range,lty = 3,col = 2,lwd = 2)
  abline(v = seq(1980,2020,10),lty = 3,col = 'grey',lwd = 1)
  abline(h = 0,lty = 3,col = 'grey',lwd = 1)
  axis(2,las = 1); axis(3,labels = F)
  if(t ==5){axis(1)}else{axis(1,labels = F)}
  box()
  mtext(letter[t], side = 3, line = 0, outer = FALSE, cex = NA,adj = 0) 
}
#Plot the trends
for(i in 1:4){
  p = nn[[i]]
  reg = lm(STR[p,t]~tdate[p])
  Y <- predict(reg, newdata=data.frame(tdate[p]))
  lines(ttime[p], y=Y,col = 2,lwd = 2)
  ye1 = substr(tdate[nn[[i]]][1],1,4)
  mo1 = month.abb[as.numeric(substr(tdate[nn[[i]]][1],6,7))]
  ye2 = substr(tail(tdate[nn[[i]]],n = 1),1,4)
  mo2 = month.abb[as.numeric(substr(tail(tdate[nn[[i]]],n = 1),6,7))]
  trend[i] = paste(round(reg$coefficients[2]*12*10^3,2)," (",mo1," ",ye1," - ",mo2," ",ye2,")",sep = "")
}
legend("bottomleft",ncol = 2,x.intersp = 0,col = 2,lwd = c(2,"","","","","")
  ,legend = c(expression(paste("Linear trends (10"^"3"," km"^"2"," year"^"-1",") :")),trend[1:2],"",trend[3:4]))

dev.off()


##############Study the evolution of the trend STR[,2]
boss = c(1992.5,1997.3,2002.1,2006.4,2011.2,2016)
plot(STR[,2],las = 1)
abline(v = boss,lty = 3)
grid()
#
mean(diff(boss)) #period of 4.7 years?

#Try to compute the period exactly
vec.ts=ts(STR[,2],freq=12,start=c(1978,11),end=c(2020,12))
vec.str=AutoSTR(vec.ts,conf=c(0.025,0.975), trace=TRUE)
plot(components(vec.str))



### TRANSFORMEE DE FOURIER LOCALE (what is going on here?)
vec.stft=stft(vec.ts)
vec.swdft=swdft(surf,n=12,m=2,smooth='daniell')
plot(vec.swdft)


###########################################################
#### PLOT DES DATAS
###########################################################

xyrange=c(-5e06,5e06)

tabs=seq(0,period,length=1000)
Cbasecomplex=eval.complex.basis(tabs,m,period)
Xline=Cbasecomplex%*%Xfd

#x11(width=15,height=5)
png("Data.png",width=15, height=5, units="in", res=200)
par(mfrow=c(1,3))
plot(Xline[,1],type="l",xlab="x",ylab="y",xlim=xyrange,ylim=xyrange)
points(Xline[1,1],pch=19)
ind=1:nobs
for(k in ind){
  points(Xline[,k],type="l")
  points(Xline[1,k])
}
#
plot(tabs,Re(Xline[,1]),type="l")
for(k in ind){
  points(tabs,Re(Xline[,k]),type="l")
}
#
plot(tabs,Im(Xline[,1]),type="l")
for(k in ind){
  points(tabs,Im(Xline[,k]),type="l")
}
dev.off()

#save(Xfd,Xline,Cbasecomplex,m,period,nobs,file = "Contour15.RData")
###########################################################
#### PLOT DES RESULTATS - STAT DE BASE SUR DATA BRUTE FITTEE FOURIER
###########################################################
graphics.off()

xyrange=c(-5e06,5e06)

tabs=seq(0,period,length=1000)
Cbasecomplex=eval.complex.basis(tabs,m,period)
Xline=Cbasecomplex%*%Xfd
Xm=apply(Xline,1,mean)

#ind=sample(1:402,5)
png("Data2.png",width=10, height=10, units="in", res=200)
plot(Xline[,1],type="l",xlab="x",ylab="y",xlim=xyrange,ylim=xyrange,col=2)
ind=1:nobs
for(k in ind){
  points(Xline[,k],type="l",col="red")
  #Sys.sleep(2)
}
points(Xline[1,],pch=19)
lines(Xm,lwd=2)
dev.off()

#Data relative to mean
x11()
par(mfrow=c(1,3))
ind=1:nobs
for(k in ind){
  plot(Xm,type="l",xlab="x",ylab="y",xlim=xyrange,ylim=xyrange,col="light grey",main=paste("nobs=",k))
  points(Xline[,k],type="l",col="red")
  points(Xline[1,k],pch=19)
  
  plot(tabs,Re(Xm),type="l",xlab="time",ylab="X move")
  points(tabs,Re(Xline[,k]),type="l",col=2)
  plot(tabs,Im(Xm),type="l",xlab="time",ylab="Y move")
  points(tabs,Im(Xline[,k]),type="l",col=2)
  text(locator(1))
}

# 
############################### RECALAGE DE PHASE BASIQUE SANS NORMALISATION
### ON DERIVE SSE EN FONCTION DE delta
### ON CHERCHE LES ZEROS DE LA DERIVEE
### ON ISOLE LE deltaelu QUI MIIMISE L'ERREUR QUADRATIQUE SSE
### ON FAIT TOURNER L'OBJET ET ON STOCKE DANS Xfdopt ET ROULE MA POULE
source('pour_optim_rot_phase.r')
##### ALGO GPA - PREND UN PEU DE TEMPS....
##########################################################################

delta=seq(0,2*pi,length=1000)
inter=range(delta)

Xfdopt=Xfd
nloop=5

numfreq=seq(-m,m,1)
puls=2*pi/period

for(k in 1:nloop){
  
  for(i in 1:nobs){
    cat("k=",k,"  i=",i,fill=TRUE)
    ind=i
    X1fd=meanfd(Xfdopt[,-i],norm=TRUE)
    X2fd=Xfdopt[,i]
    
    mindelta=uniroot.all(fderiv_phase, interval=inter, lower=min(inter), upper=max(inter), 
      tol=.Machine$double.eps^0.2, maxiter=1000, n=500,X1fd=X1fd,X2fd=X2fd,m=m,period=period)
    
    resid=SSE_phase(mindelta,X1fd,X2fd,m,period)
    deltaelu=mindelta[which.min(resid)]
    
    R=diag(complex(mod=1,arg=puls*numfreq*deltaelu))
    Xfdopt[,i]=R%*%Xfdopt[,i]
  }
  
}
#save(Xfdopt,file = "SIE_opt_rot.RData")

