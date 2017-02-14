apt = function( X, Xpred = NULL, Omega.type = "unit", max.resol = 10, rho0=0.2, rho0.mode = 0,tran.mode=1,
                                 lognu.lb=-1, lognu.ub=4, n.grid=5, n.s=5,beta=0.1,n.post.samples=0){

  X = as.matrix(X)
  p = ncol(X)

  if (Omega.type == "unit") {

    Omega = matrix(rep(c(0,1),2*p),nrow=p,ncol=2,byrow=TRUE)

  } else if(Omega.type == "standardized") {

    Omega = t(apply(rbind(X,Xpred),2,range))
    Omega[,2] = Omega[,2]*1.00001
  }
  else
  {
    print("ERROR: Sample space 'Omega' incorrectly specified")
    return(0);
  }


  if (is.null(Xpred)) {
    Xpred = matrix(rep(0.5,p),ncol=p,nrow=1)
  } else {
    Xpred = as.matrix(Xpred)
  }

  ans = fitPTTcpp(X,Xpred,Omega,max.resol,rho0,rho0.mode,tran.mode,lognu.lb,lognu.ub,n.grid,n.s,beta,n.post.samples)

  part.hmap = matrix(unlist(ans$part_points_hmap),byrow=TRUE,ncol=2*p+2)
  colnames(part.hmap)[2*p+(1:2)] = c("level","state")
  part.hmap[part.hmap[,"state"] == n.s+1, "state"] = Inf

  ans[which(names(ans) == "part_points_hmap")][[1]] = part.hmap

  if (n.post.samples > 0) {
    for (i in 1:n.post.samples) {
      part.post.sample = matrix(unlist(ans$part_points_post_samples[[i]]),byrow=TRUE,ncol=2*p+1)
      colnames(part.post.sample)[2*p+1] = "level"
      nu.post.sample = matrix(unlist(ans$nu_and_prob_post_samples[[i]]),byrow=TRUE,ncol=2)
      colnames(nu.post.sample) = c("nu","logp")

      ans$part_points_post_samples[[i]] = cbind(part.post.sample,nu.post.sample)
    }

  }
  else {
    ans[which(names(ans) == "part_points_post_samples")] = NULL
  }
  ans[which(names(ans) == "nu_and_prob_post_samples")] = NULL

  return(ans)
}

## Fitting an optional Polya tree model
opt = function( X, Xpred = NULL, Omega.type = "unit", max.resol = 10, rho0=0.2, rho0.mode = 0, n.post.samples=0) {

  ans = apt(X=X, Xpred=Xpred, Omega.type=Omega.type, max.resol=max.resol, rho0=rho0, rho0.mode = rho0.mode,tran.mode=1,
                           lognu.lb=0, lognu.ub=0, n.grid=1, n.s=1,beta=0,n.post.samples=n.post.samples)

  return(ans)
}
