markov.apt.density.kD= function( x.mat, max.dim = 10, rho0=0.2, rho0.mode = 0,tran.mode=1,lognu.lb=-1, lognu.ub=4, n.grid=5, n.s=2,beta=0.1,x.mat.pred=NULL,n.post.samples=0){
  
  if (!is.loaded("../src/markov_apt_kD_density_C.so"))
    dyn.load('../src/markov_apt_kD_density_C.so')

  if (is.vector(x.mat)) {
    x.mat = rbind(NULL,x.mat)
  }

  
  n.pred = as.integer(ncol(x.mat))
  if (is.null(x.mat.pred)) {
    x.mat.pred=rbind(NULL,rep(0.5,n.pred))
  } else {
    if (is.vector(x.mat.pred)) {
      x.mat.pred = rbind(NULL,x.mat.pred)
    }
  }

  
  ans = .Call('markov_apt_kD_C', as.numeric(t(x.mat)), as.integer(ncol(x.mat)), as.integer(max.dim), as.numeric(rho0), as.integer(rho0.mode),as.integer(tran.mode),as.numeric(lognu.lb),as.numeric(lognu.ub),as.integer(n.grid),as.integer(n.s),as.numeric(beta),as.integer(n.post.samples),as.numeric(t(x.mat.pred)))

  if (n.post.samples > 0) {

    names(ans)=c("logrho","logphi","partition","ppd","part.post.samples")

    for (i in 1:n.post.samples) {
      nvar = ncol(x.mat)
      colnames(ans$part.post.samples[[i]])[(2*nvar+1):(2*nvar+3)] = c("level","nu","logp")
    }
    
  } else {
    names(ans)=c("logrho","logphi","partition","ppd")
  }
  
  return(ans)
}



