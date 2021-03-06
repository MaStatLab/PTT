x=seq(from=0.0001,to=0.9999,by=0.002)
set.seed(100)
nobs = 5000
max.resol=11
size.max = nobs

p1=0.1;p2=0.3;p3=0.4;p4=0.2
sample.ind = runif(size.max)

X = as.matrix((sample.ind<=p1)*runif(size.max)+(p1<sample.ind & sample.ind<=p1+p2)*(0.25+rbeta(size.max,1,1)*0.25)+(p1+p2<sample.ind & sample.ind<=p1+p2+p3)*(0.25+rbeta(size.max,2,2)*0.25) +(sample.ind>p1+p2+p3)*rbeta(size.max,5000,2000))

true.den = p1*dunif(x) + p2*dunif(x,0.25,0.5) + p3 *4*dbeta((x-0.25)*4,2,2) + p4*dbeta(x,5000,2000)

## Fit a Markov APT
ans = apt(X=X,max.resol=max.resol,rho0=0.2)
pred.den = apt(X=X,Xpred=x,max.resol=max.resol,rho0=0.2)$pred

xlim = c(0,1)
ylim = c(0,15)
plot(x,pred.den,type='l',xlim=xlim,ylim=ylim,xlab="x",ylab="Density",main="PPD for Markov-APT")
par(new=TRUE)
plot(x,true.den,type='l',lty=2,col='red',xlab="",ylab="",xlim=xlim,ylim=ylim)


## Fit OPT instead
ans.opt = opt(X=X,max.resol=max.resol,rho0=0.2)
pred.den.opt = opt(X=X,Xpred=x,max.resol=max.resol,rho0=0.2)$pred


xlim = c(0,1)
ylim = c(0,15)
plot(x,pred.den.opt,type='l',xlim=xlim,ylim=ylim,xlab="x",ylab="Density",main="PPD for OPT")
par(new=TRUE)
plot(x,true.den,type='l',lty=2,col='red',xlab="",ylab="",xlim=xlim,ylim=ylim)




