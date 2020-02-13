\name{apt}
\alias{apt}
\title{Adaptive Polya tree}

\description{
This function fits the adaptive Polya tree model with independent or stochastically increasing (Markov) shrinkage.
}
\usage{
apt( X, Xpred = NULL, Omega.type = "unit", max.resol = 10, rho0=0.2, rho0.mode = 0,tran.mode=1,
                            lognu.lb=-1, lognu.ub=4, n.grid=5, n.s=5,beta=0.1,n.post.samples=0)
}

\arguments{
  \item{X}{Matrix of the data. Each row is an observation. Each column is a dimension.}

  \item{Xpred}{Matrix of locations in the sample space where posterior predictive density is to be computed.}

  \item{Omega.type}{Specifies the type of samples space. Under the default \code{"unit"} option, the sample space is the unit rectangle [0,1]^k where k is the number of dimensions (i.e., the number of columns in \code{X}. Under the \code{"standardized"} option, the sample space is a hypercube defined by the empirical range of each dimension of \code{X}.}

  \item{max.resol}{The maximum resolution level of the partition tree. The default is \code{10}. The maximum is \code{14}.}

  \item{rho0.mode}{The way prior probability for complete shrinkage (i.e., the stopping state) is specified. The default is \code{rho0.mode=0} under which the prior probability for complete shrinkage is the same constant \code{rho0} for all nodes in the partition. If \code{rho0.mode=1} then the probability grows to 1 with the level of the node quadratically: \code{1-(1-rho0)^k}, where k is the level of the node. If \code{rho0.mode=2} then the probability grows to 1 with the level of the node exponentially.}

  \item{rho0}{Prior probability for complete shrinkage.}

  \item{tran.mode}{Type of shrinkage. The default \code{tran.mode=1} corresponds to stochastically increasing shrinkage with uniform prior probability over higher shrinkage states. If \code{tran.mode=2} then shrinakage is stochastically increasing with kernel-based exponential probability over higher shrinkage states. If \code{tran.mode=0} then shrinakage is independent over nodes and there is no stochastically increasing shrinkage.}

  \item{lognu.lb}{The lowerbound of the support for the shrinkage parameter on the log 10 scale.}

  \item{lognu.ub}{The upperbound of the support for the shrinakge parameter on the log 10 scale.}

  \item{n.grid}{The number of grid points used for numerical integration on the shrinkage parameter in computing the marginal likelihood.}

  \item{n.s}{The number of shrinkage states (excluding the complete shrinkage state).}

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

## Fit a Markov APT
ans = apt(X=X,max.resol=max.resol)
pred.den = apt(X=X,Xpred=x,max.resol=max.resol)$pred

xlim = c(0,1)
ylim = c(0,15)
plot(x,pred.den,type='l',xlim=xlim,ylim=ylim,xlab="x",ylab="Density",main="PPD for Markov-APT")
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

markov.apt.2D.fit = apt(X=obs.mat,Xpred=xy.grid,n.post.samples = 1000)

## Plot the posterior sample densities
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

}

