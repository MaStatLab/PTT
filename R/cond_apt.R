## Fitting a conditional APT
cond.apt = function( X, Y, Xpred = NULL, Ypred=NULL, OmegaX.type = "unit", OmegaY.type = "unit", max.resX = 5, max.resY=8,
                     rho0.X=0.2, rho0.mode.X = 0, rho0.Y=0.2, rho0.mode.Y = 0,
                     tran.mode=1,lognu.lb=-1, lognu.ub=4, n.grid=5, n.s=5,beta=0.1,n.post.samples=0)
{

  X = as.matrix(X)
  Y = as.matrix(Y)
  p.X = ncol(X)
  p.Y = ncol(Y)
  if (OmegaX.type == "unit") {

    OmegaX = matrix(rep(c(0,1),2*p.X),nrow=p.X,ncol=2,byrow=TRUE)

  } else if(OmegaX.type == "standardized") {

    OmegaX = t(apply(rbind(X,Xpred),2,range))
    OmegaX[,2] = OmegaX[,2]*1.00001
  }
  else
  {
    print("ERROR: Sample space for X 'Omega X' incorrectly specified")
    return(0);
  }

  if (OmegaY.type == "unit") {

    OmegaY = matrix(rep(c(0,1),2*p.Y),nrow=p.Y,ncol=2,byrow=TRUE)

  } else if(OmegaY.type == "standardized") {

    OmegaY = t(apply(rbind(Y,Ypred),2,range))
    OmegaY[,2] = OmegaY[,2]*1.00001
  }
  else
  {
    print("ERROR: Sample space for Y 'Omega Y' incorrectly specified")
    return(0);
  }

  if (is.null(Xpred)) {
    Xpred = matrix(rep(0.5,p.X),ncol=p.X,nrow=1)
    Ypred = matrix(rep(0.5,p.Y),ncol=p.Y,nrow=1)
  } else {
    Xpred = as.matrix(Xpred)
    Ypred = as.matrix(Ypred)
  }

  ans = fitCondPTTcpp(X,Y,Xpred,Ypred,OmegaX,OmegaY,max.resX,max.resY,rho0.X,rho0.mode.X,rho0.Y,rho0.mode.Y,
                      tran.mode,lognu.lb,lognu.ub,n.grid,n.s,beta,n.post.samples)
  part.hmap = matrix(unlist(ans$part_points_hmap),byrow=TRUE,ncol=2*p.X+2)
  colnames(part.hmap) = c(paste(rep(paste("X",1:p.X,sep=""),rep(2,p.X)),c(".l",".u"),sep=""),"level","state")
  part.hmap[part.hmap[,"state"] == 2,"state"] = Inf
  # print(part.hmap)

  rhos.hmap = matrix(unlist(ans$rhos_hmap),byrow=TRUE,ncol=1)
  colnames(rhos.hmap)[1] = c("rho")
  ans[which(names(ans) == "part_points_hmap")][[1]] = cbind(part.hmap,rhos.hmap)
  ans[which(names(ans) == "rhos_hmap")] = NULL

  if (n.post.samples > 0) {
    for (i in 1:n.post.samples) {
      part.post.sample = matrix(unlist(ans$part_points_post_samples[[i]]),byrow=TRUE,ncol=2*p.X+1)
      colnames(part.post.sample) = c(paste(rep(paste("X",1:p.X,sep=""),rep(2,p.X)),c(".l",".u"),sep=""),"level")

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


## Fitting a conditional OPT
cond.opt = function( X, Y, Xpred = NULL, Ypred=NULL, OmegaX.type = "unit", OmegaY.type = "unit", max.resX = 7, max.resY=7,
                     rho0.X=0.2, rho0.mode.X = 0, rho0.Y=0.2, rho0.mode.Y = 0, tran.mode=1,
                     lognu.lb=-1, lognu.ub=4, n.grid=5, n.s=5,beta=0.1,n.post.samples=0){

  ans = cond.apt(X,Y,Xpred,Ypred,OmegaX.type,OmegaY.type,max.resX,max.resY,rho0.X,rho0.mode.X,rho0.Y,rho0.mode.Y,
                 tran.mode=1,lognu.lb=0, lognu.ub=0, n.grid=1, n.s=1,beta=0,n.post.samples=n.post.samples)

  return(ans)
}

