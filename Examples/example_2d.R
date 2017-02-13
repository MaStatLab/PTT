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
max.resol = 11

markov.apt.2D.fit = apt(X=obs.mat,Xpred=xy.grid,max.resol=max.resol,rho0=0.2,
                        tran.mode=2,beta=0,n.post.samples = n.post.sample)

for (i in 1:n.post.sample) {

  part.post.sample = markov.apt.2D.fit$part_points_post_samples[[i]]

  terminal.part = part.post.sample[which(part.post.sample[,"nu"] == Inf),]
  part.to.plot = cbind(terminal.part,den = exp(terminal.part[,"logp"])*2^(terminal.part[,"level"]))
  part.to.plot[,2] = part.to.plot[,2]+1
  part.to.plot[,4] = part.to.plot[,4]+1
  part.to.plot[,1:4] = part.to.plot[,1:4] / 2^max.resol
  colnames(part.to.plot) = c("xmin","xmax","ymin","ymax","level","nu","lopp","den")

  xlim=c(0,1)
  ylim=c(0,1)
  zlim=c(0,max(part.to.plot[,"den"]))
  border=FALSE
  plot.part(part.to.plot,xlim=xlim,ylim=ylim,zlim=zlim,border=border,plot.scale=FALSE,color.fun=topo.colors)

  Sys.sleep(0.5)
}
