\name{opt}
\alias{opt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Optional Polya tree}

\description{
This function fits the optional Polya tree model introduced in Wong and Ma (2010).
}
\usage{

opt(X, Xpred = NULL, Omega.type = "unit", max.resol = 10, rho0=0.5, rho0.mode = 0, n.post.samples=0)
}

\arguments{
  \item{X}{Matrix of the data. Each row is an observation. Each column is a dimension.}

  \item{Xpred}{Matrix of locations in the sample space where posterior predictive density is to be computed.}

  \item{Omega.type}{Specifies the type of samples space. Under the default \code{"unit"} option, the sample space is the unit rectangle [0,1]^k where k is the number of dimensions (i.e., the number of columns in \code{X}. Under the \code{"standardized"} option, the sample space is a hypercube defined by the empirical range of each dimension of \code{X}.}


  \item{max.resol}{The maximum resolution level of the partition tree. The default is \code{10}. The maximum is \code{14}.}

  \item{rho0.mode}{The way prior probability for complete shrinkage (i.e., the stopping state) is specified. The default is \code{rho0.mode=0} under which the prior probability for complete shrinkage is the same constant \code{rho0} for all nodes in the partition. If \code{rho0.mode=1} then the probability grows to 1 with the level of the node quadratically: \code{1-(1-rho0)^k}, where k is the level of the node. If \code{rho0.mode=2} then the probability grows to 1 with the level of the node exponentially.}


  \item{n.post.samples}{The number of posterior samples to be drawn.}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{logrho}{Posterior probability of complete shrinkage on the entire sample space.}
\item{logphi}{Log marginal likelihood.}
\item{part_points_hmap}{A representative (the hMAP) partition on the sample space.}
\item{predictive_densities}{The posterior predictive densities at the points specified in \code{Xpred}.}
\item{part_points_post_samples}{The partition for each posterior sample.}
\item{nu_and_prob_post_samples}{The shrinkage and probability assignment of each posterior sample.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Li Ma
}

\examples{
### 1D Example
x=seq(from=0.0001,to=0.9999,by=0.002)
set.seed(100)
nobs = 5000
max.resol=11
size.max = nobs

p1=0.1;p2=0.3;p3=0.4;p4=0.2
sample.ind = runif(size.max)

X = as.matrix((sample.ind<=p1)*runif(size.max)+(p1<sample.ind & sample.ind<=p1+p2)*(0.25+rbeta(size.max,1,1)*0.25)+(p1+p2<sample.ind & sample.ind<=p1+p2+p3)*(0.25+rbeta(size.max,2,2)*0.25) +(sample.ind>p1+p2+p3)*rbeta(size.max,5000,2000))

true.den = p1*dunif(x) + p2*dunif(x,0.25,0.5) + p3 *4*dbeta((x-0.25)*4,2,2) + p4*dbeta(x,5000,2000)

## Fit an OPT
ans = opt(X=X)
pred.den = opt(X=X,Xpred=x)$pred

xlim = c(0,1)
ylim = c(0,15)
plot(x,pred.den,type='l',xlim=xlim,ylim=ylim,xlab="x",ylab="Density",main="PPD for OPT")
par(new=TRUE)
plot(x,true.den,type='l',lty=2,col='red',xlab="",ylab="",xlim=xlim,ylim=ylim)

x.grid = seq(0.0000001,0.9999999,by=0.01)
y.grid = seq(0.0000001,0.9999999,by=0.01)
xy.grid = expand.grid(x.grid,y.grid)

x.grid.thin = x.grid[(1:length(x.grid))%% 4 == 0]
y.grid.thin = x.grid[(1:length(y.grid))%% 4 == 0]


### 2D Example

nobs = 3000;
sample.ind=runif(nobs)

mean=c(0.5,0.5);sigma=diag(c(0.01,0.0064))
mean.local=c(0.8,0.2); sigma.local = sigma/25


norm.obs.mat = mvtnorm::rmvnorm(n=nobs,mean=mean, sigma=sigma)
norm.obs.mat.local = mvtnorm::rmvnorm(n=nobs, mean=mean.local, sigma=sigma.local)

p1 = 0.85
obs.mat = (sample.ind < p1) * norm.obs.mat + (sample.ind >= p1 ) * norm.obs.mat.local
true.den = function(xy.grid) { p1*dmvnorm(xy.grid,mean=mean,sigma=sigma) + (1-p1)*dmvnorm(xy.grid,mean=mean.local,sigma=sigma.local) }
true.den.grid = true.den(xy.grid)

opt.2D.fit = opt(X=obs.mat,Xpred=xy.grid,n.post.samples = 1000)

## Plot the posterior sample densities
for (i in 1:n.post.sample) {

  part.post.sample = opt.2D.fit$part_points_post_samples[[i]]

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

}
