
library(mvtnorm)

transform=FALSE
x.grid = seq(0.0000001,0.9999999,by=0.01)
y.grid = seq(0.0000001,0.9999999,by=0.01)
xy.grid = expand.grid(x.grid,y.grid)

x.grid.thin = x.grid[(1:length(x.grid))%% 4 == 0]
y.grid.thin = x.grid[(1:length(y.grid))%% 4 == 0]


nobs = 3000;
sample.ind=runif(nobs)

mean=c(0.5,0.5);sigma=diag(c(0.01,0.0064))
mean.local=c(0.8,0.2); sigma.local = sigma/25


norm.obs.mat =rmvnorm(n=nobs,mean=mean, sigma=sigma)
norm.obs.mat.local = rmvnorm(n=nobs, mean=mean.local, sigma=sigma.local)

p1 = 0.85
obs.mat = (sample.ind < p1) * norm.obs.mat + (sample.ind >= p1 ) * norm.obs.mat.local
true.den = function(xy.grid) { p1*dmvnorm(xy.grid,mean=mean,sigma=sigma) + (1-p1)*dmvnorm(xy.grid,mean=mean.local,sigma=sigma.local) }
true.den.grid = true.den(xy.grid)




n.post.sample=500
max.dim = 11

library(grid)

markov.apt.2D.fit = markov.apt.density(obs.mat,max.dim=max.dim,rho0=0.1,rho0.mode=0,n.s=5,tran.mode=2,beta=0,
                                       lognu.lb=-1,lognu.ub=4,n.grid=5,x.mat.pred=xy.grid,n.post.samples = n.post.sample)

for (i in 1:n.post.sample) {

  part.post.sample = markov.apt.2D.fit$part_points_post_samples[[i]]

terminal.part = part.post.sample[which(part.post.sample[,"nu"] == Inf),]
part.to.plot = cbind(terminal.part,den = exp(terminal.part[,"logp"])*2^(terminal.part[,"level"]))
part.to.plot[,2] = part.to.plot[,2]+1
part.to.plot[,4] = part.to.plot[,4]+1
part.to.plot[,1:4] = part.to.plot[,1:4] / 2^max.dim
colnames(part.to.plot) = c("xmin","xmax","ymin","ymax","level","nu","lopp","den")

xlim=c(0,1)
ylim=c(0,1)
zlim=c(0,max(part.to.plot[,"den"]))
border=FALSE
plot.part(part.to.plot,xlim=xlim,ylim=ylim,zlim=zlim,border=border,plot.scale=FALSE,color.fun=topo.colors)

Sys.sleep(0.5)
}
