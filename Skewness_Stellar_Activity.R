rm(list=ls(all=TRUE))  # clear all variables

################################################
###import dataset for Alpha Centauri B
which.dataset_ccf = "data/HD128621_CCF.rdb"
ccf=read.table(which.dataset_ccf, head=T)

which.dataset_time = "data/HD128621_data.rdb"
time=read.table(which.dataset_time, head=T)
################################################

################################################
###import dataset for Gliese 785
which.dataset_ccf = "data/HD192310_CCF.rdb"
ccf=read.table(which.dataset_ccf, head=T)

which.dataset_time = "data/HD192310_data.rdb"
time=read.table(which.dataset_time, head=T)
################################################

################################################
###preliminary operations: the data are factors, we need to adjust them
head(ccf)

ccf=ccf[-1,]
head(ccf)
dim(ccf)

#attach(ccf)

#rv=as.numeric(as.character(rv))

matrix_ccf=matrix(NA,dim(ccf)[1],dim(ccf)[2])
dim(matrix_ccf)

for(i in 1:dim(ccf)[2])
matrix_ccf[,i]=as.numeric(as.character(ccf[,i]))
#matrix_ccf[,-1]=matrix_ccf[,-1]/10000000

################################################
####cut the queue?
#matrix_ccf=matrix(NA,100,dim(ccf)[2])
#dim(matrix_ccf)

#for(i in 1:dim(ccf)[2])
#matrix_ccf[,i]=as.numeric(as.character(ccf[-c(1:30,131:dim(ccf)[1]),i]))
################################################

#colfunc = colorRampPalette(c("black","yellow"))
#col=colfunc(2000)
#plot(1-matrix_ccf[,2]~matrix_ccf[,1], type="l", col=col[1],ylim=c(10000000,255000000))
#for(i in 3:2000)
#lines(matrix_ccf[,i]~matrix_ccf[,1],type="l", col=col[i])

head(matrix_ccf)

dim(time)
head(time)
time=time[-1,]
head(time)
dim(time)

################################################
##############
####CCF Examples: Least Squares fit Skew Normal with DP:
install.packages("sn")
library(sn)
#library(csn)
rv = matrix_ccf[,1]
x=rv

###one particulat CCF is selected
y=1-(matrix_ccf[,451]/max(matrix_ccf[,451]))

#plot(x,y)

###start to the mle of the normal
fitG =
function(x,y,mu,sig,scale){

  f = function(p){
  	#d = p[3]*(1-(dnorm(x,mean=p[1],sd=p[2])))
    d = p[3]*(dnorm(x,mean=p[1],sd=p[2]))
    sum((d-y)^2)
  }

  optim(c(mu,sig,scale),f,method="L-BFGS-B", lower=c(-Inf,0,-Inf), upper=c(Inf,Inf,Inf), hessian=T)
 }
ml=fitG(x,y,mean(x),1,1)
p=ml$par

###fit from a skew normal: direct rep (DP)
fitSN =
function(x,y,mu,sig,a,scale){

  f = function(p){
    d = p[4]*(dsn(x,xi=p[1],omega=p[2],alpha=p[3]))
    sum((d-y)^2)
  }

  optim(c(mu,sig,a,scale),f,method="BFGS", lower=c(-Inf,0.1,-Inf,-Inf), upper=c(Inf,Inf,Inf,Inf), hessian=T)
 }

###central rep from a package: slow and not sure about this
#fitSNc =
#function(x,y,mean,sigma,gamma,scale){
#
#  f = function(p){
#    d = p[4]*(dcsn(x,mu=p[1],sigma=p[2],gamma=p[3],delta=1,nu=0))
#
#    sum((d-y)^2)
#  }
#
#  optim(c(mean,sigma,gamma,scale),f,method="BFGS", lower=c(-Inf,0.01,-Inf,-Inf), upper=c(Inf,Inf,Inf,Inf), hessian=T, control=list(trace=TRUE))
# }

####fit from a skew normal: central rep (CP)
fitSN_cp =
function(x,y,mu,sig,a,scale){

  f = function(z){
  	
  	dp=cp2dp(z[-4], family="SN")
    d = z[4]*(dsn(x,xi=dp[1],omega=dp[2],alpha=dp[3]))
    sum((d-y)^2)
  }

  optim(c(mu,sig,a,scale),f,method="BFGS", lower=c(-Inf,0.1,-Inf,-Inf), upper=c(Inf,Inf,Inf,Inf), hessian=T)
 }
 
####fit from a skew t: central rep (CP)
fitST_cp =
function(x,y,mu,sig,a,scale){

  f = function(z){
  	
  	dp=cp2dp(z[-4], family="SN")
    d = z[4]*(dst(x,xi=dp[1],omega=dp[2],alpha=dp[3]))
    sum((d-y)^2)
  }

  optim(c(mu,sig,a,scale),f,method="BFGS", lower=c(-Inf,0.1,-Inf,-Inf), upper=c(Inf,Inf,Inf,Inf), hessian=T)
 }

###fit with the direct rep (DP)
skew.est=fitSN(x,y,p[1],p[2],5,p[3])
skew.est

###fit with the central rep (CP)
skew.est_cp=fitSN_cp(x,y,p[1],p[2],.1,p[3])
skew.est_cp

###fit with the skew t central rep (CP)
skew.t.est_cp=fitST_cp(x,y,p[1],p[2],.1,p[3])
skew.t.est_cp

###come back to the direct from the plot: note that they're equal to skew.est$par, obviously !!!
###[-4] because it is the norm. constant.
dp.est=cp2dp(skew.est_cp$par[-4], family="SN")

#dev.new()
plot(x,y)
#lines(x, p[3]*(1-(dnorm(x, p[1], p[2]))),col="blue")
lines(x, p[3]*(dnorm(x, p[1], p[2])),col="blue")
lines(x, skew.est$par[4]*(dsn(x, skew.est$par[1], skew.est$par[2], skew.est$par[3])),col="red")
lines(x, skew.est_cp$par[4]*(dsn(x, dp.est[1], dp.est[2], dp.est[3])),col="yellow")
legend("topright", c("Normal Fit","Skew Normal Fit", "Skew Normal Fit CP" ), col = c("blue","red","yellow"), lty = c(1,1,1),lwd = 3)
#dev.off()

####confront between mean by using a skew normal
#loc.skewnormal.umberto=c(NA,dim(matrix_ccf)[2])
mean.skewnormal.umberto=c(NA,dim(matrix_ccf)[2])
#alpha.skewnormal.umberto=c(NA,dim(matrix_ccf)[2])
gamma.skewnormal.umberto=c(NA,dim(matrix_ccf)[2])

#mean.skew.t.umberto=c(NA,dim(matrix_ccf)[2])
#gamma.skew.t.umberto=c(NA,dim(matrix_ccf)[2])

chisq.norm=c(NA,dim(matrix_ccf)[2])
chisq.skew=c(NA,dim(matrix_ccf)[2])

for(i in 2:dim(matrix_ccf)[2])
{ 
	 x=rv
     y=1-(matrix_ccf[,i]/max(matrix_ccf[,i]))
     
     ###start to the mle of the normal
     ml=fitG(x,y,mean(x),1,1)
     p=ml$par
     chisq.norm[i-1]=ml$value
     
     #skew.est=fitSN(x,y,p[1],p[2],5,p[3])
     skew.est=fitSN_cp(x,y,p[1],p[2],.5,p[3])
     chisq.skew[i-1]=skew.est$value
     
     mean.skewnormal.umberto[i-1]=skew.est$par[1]
    
     #loc.skewnormal.umberto[i-1]=skew.est$par[1]
     #mean.skewnormal.umberto[i-1]=dp2cp(skew.est$par[-4], family="SN")[1]	
     
     #alpha.skewnormal.umberto[i-1]=skew.est$par[3]
     #gamma.skewnormal.umberto[i-1]=dp2cp(skew.est$par[-4], family="SN")[3]
     
     gamma.skewnormal.umberto[i-1]=skew.est$par[3]
     
     #skew.t.est=fitST_cp(x,y,p[1],p[2],.5,p[3])
     #mean.skew.t.umberto[i-1]=skew.t.est$par[1]
     #gamma.skew.t.umberto[i-1]=skew.t.est$par[3]
}

###store the results in a file
write.table(gamma.skewnormal.umberto, file="gamma.skewnormal.umberto.dat")
write.table(mean.skewnormal.umberto, file="mean.skewnormal.umberto.dat")
####################################################
head(time)
julian.day=as.numeric(as.character(time[,1]))
mean.bis.span=as.numeric(as.character(time[,3]))
bis.span=as.numeric(as.character(time[,4]))

###for Alpha Centauri only
mean.skewnormal.xavier=as.numeric(as.character(time[,5]))
gamma.skewnormal.xavier=as.numeric(as.character(time[,6]))

###set of Plots Alpha Centauri only
#dev.new()
boxplot(mean.bis.span,mean.skewnormal.xavier,mean.skewnormal.umberto, main="mean_bis vs. mean_Xavier vs. mean_Umberto of the SN", border=c("black","blue","red"))
#dev.off()

#dev.new()
boxplot(bis.span,gamma.skewnormal.xavier,gamma.skewnormal.umberto, main="mean_bis vs. mean_Xavier vs. mean_Umberto of the SN", border=c("black","blue","red"))
#dev.off()

#dev.new()
plot(julian.day, mean.bis.span, col="black", ylim=c(-22.75,-22.69),lwd=.5)
points(julian.day, mean.skewnormal.xavier, lwd=.5, col="blue")
points(julian.day, mean.skewnormal.umberto, lwd=.5, col="red")
legend("topright", c("Normal Bisector","Skew Normal Xavier", "Skew Normal Umberto"), col = c("black","blue","red"), lty = c(1,1,2),lwd = 3)
#dev.off()

#dev.new()
#par(new = TRUE)
plot(julian.day, bis.span, col="black",lwd=0.3, ylim=c(-0.004,0.02))
points(julian.day, gamma.skewnormal.xavier, lwd=0.4, col="blue")
points(julian.day, gamma.skewnormal.umberto, lwd=0.5,col="red")
legend("topright", c("Normal Bisector","Skew Normal Xavier", "Skew Normal Umberto"), col = c("black","blue","red"), lty = c(1,1,2),lwd = 3)
#dev.off()

###Set of plots for Gliese 785 only
#dev.new()
boxplot(mean.bis.span,mean.skewnormal.umberto, main="mean_bis vs. mean_Umberto of the SN", border=c("black","red"))
#dev.off()

#dev.new()
boxplot(bis.span,gamma.skewnormal.umberto, main="mean_bis vs. mean_Umberto of the SN", border=c("black","red"))
#dev.off()

#dev.new()
plot(julian.day, mean.bis.span, col="black", ylim=c(min(mean.bis.span)-0.002,max(mean.bis.span)+0.02),lwd=0.3)
points(julian.day, mean.skewnormal.umberto, lwd=0.3, col="red")
legend("topright", c("Normal Bisector","Skew Normal Umberto"), col = c("black","red"), lty = c(1,1),lwd = 3)
#dev.off()

#dev.new()
#par(new = TRUE)
plot(julian.day, bis.span, col="black",lwd=0.3, ylim=c(-0.04,max(bis.span)+0.02))
points(julian.day, gamma.skewnormal.umberto, lwd=0.3,col="red")
legend("topright", c("Normal Bisector","Skew Normal Umberto"), col = c("black","red"), lty = c(1,1),lwd = 3)
#dev.off()

###Residuals: chisq norm vs chisq skew norm
#dev.new()
boxplot(chisq.norm,chisq.skew, main="Residuals: Norm vs. Skew", border=c("black","red"),ylim=c(0,0.01))
#dev.off()

#dev.new()
#par(new = TRUE)
plot(julian.day, chisq.norm, col="black",lwd=0.3, ylim=c(0.001,0.002))
points(julian.day, chisq.skew, lwd=0.4,col="red")
legend("topright", c("Normal Bisector","Skew Normal Umberto"), col = c("black","red"), lty = c(1,1),lwd = 3)
#dev.off()

length(which(chisq.norm>chisq.skew))/dim(ccf)[2]

which(chisq.norm<chisq.skew)
###Interesting: cases in which the LS with the normal is smaller than the one obtained with the skew

#dev.new()
plot(julian.day[which(chisq.norm<chisq.skew)], mean.bis.span[which(chisq.norm<chisq.skew)], col="black", ylim=c(min(mean.bis.span)-0.002,max(mean.bis.span)+0.02),lwd=0.3)
points(julian.day[which(chisq.norm<chisq.skew)], mean.skewnormal.umberto[which(chisq.norm<chisq.skew)], lwd=0.3, col="red")
legend("topright", c("Normal Bisector","Skew Normal Umberto"), col = c("black","red"), lty = c(1,1),lwd = 3)
#dev.off()

#dev.new()
#par(new = TRUE)
plot(julian.day[which(chisq.norm<chisq.skew)], bis.span[which(chisq.norm<chisq.skew)], col="black",lwd=0.3, ylim=c(-0.04,max(bis.span)+0.02))
points(julian.day[which(chisq.norm<chisq.skew)], gamma.skewnormal.umberto[which(chisq.norm<chisq.skew)], lwd=0.3,col="red")
legend("topright", c("Normal Bisector","Skew Normal Umberto"), col = c("black","red"), lty = c(1,1),lwd = 3)
#dev.off()










