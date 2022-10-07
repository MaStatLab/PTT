library(mvtnorm)


### Generate a data set
nobs = 2000;
sample.ind=runif(nobs)

mean=c(0.5,0.5);sigma=diag(c(0.01,0.0064))
mean.local=c(0.8,0.2); sigma.local = sigma/25


norm.obs.mat = rmvnorm(n=nobs,mean=mean, sigma=sigma)
norm.obs.mat.local = rmvnorm(n=nobs, mean=mean.local, sigma=sigma.local)
p1 = 0.85
X = (sample.ind < p1) * cbind(rbeta(nobs,20,10),rbeta(nobs,1,1)) + (sample.ind >= p1) * cbind(rbeta(nobs,10,20),rbeta(nobs,1,1))
Y = (sample.ind < p1) * norm.obs.mat + (sample.ind >= p1 ) * norm.obs.mat.local

### The values at which predictive density is evaluated
y1.grid = seq(0.0000001,0.9999999,by=0.01)
y2.grid = seq(0.0000001,0.9999999,by=0.01)
y.grid = expand.grid(y1.grid,y2.grid)
xpred = c(0.8,0.5)
x.grid = matrix(rep(xpred,nrow(y.grid)),byrow=TRUE,ncol=2)

## Fit the cond-OPT
cond.opt.2D.fit = cond.opt(X=X,Y=Y,Xpred=x.grid,Ypred=y.grid,rho0.X = 0.5,rho0.Y=0.5,max.resX=7,max.resY=7)
## Plot the predictive density for X=(0.2,0.5)
image(x=y1.grid,y=y2.grid,z=matrix(cond.opt.2D.fit$predictive_densities,nrow=length(y1.grid)),col=topo.colors(100),
      main=paste("Predictive conditional density at X=(",xpred[1],",",xpred[2],")",sep=""),xlab="Y1",ylab="Y2")


## Plot the HMAP partition on the predictor space
part.hmap = cond.opt.2D.fit$part_points_hmap

terminal.part = part.hmap[which(part.hmap[,"state"] == Inf),]
part.to.plot = cbind(terminal.part,den = 1)
part.to.plot[,2] = part.to.plot[,2]+1
part.to.plot[,4] = part.to.plot[,4]+1
part.to.plot[,1:4] = part.to.plot[,1:4] / 2^10
colnames(part.to.plot)[1:4] = c("xmin","xmax","ymin","ymax")

## Plot the entire predictor space partition
plot(x=c(0,1),y=c(0,1),type="n",axes=FALSE,xlab="X1",ylab="X2",main="HMAP partition on the predictor space")
rect(xleft=part.to.plot[,1],xright=part.to.plot[,2],ybottom=part.to.plot[,3],ytop=part.to.plot[,4])
