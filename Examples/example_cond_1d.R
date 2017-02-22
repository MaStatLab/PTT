## An example for conditional optional Polya tree for univariate response and predictor
size=500
size.pred = 100 ## Size of a testing set

x=rbeta(size+size.pred, 2,2)
y=rbeta(size+size.pred, 10,30)
y[x<0.25] = rbeta(sum(x<0.25),30,20)
y[x>0.5] = rbeta(sum(x>0.5),0.5,0.5)



X = matrix(x[1:size])
Y = matrix(y[1:size])
Xpred = matrix(x[(size+1):(size+size.pred)])
Ypred = matrix(y[(size+1):(size+size.pred)])

## Fit cond-OPT
ans.cond.opt=cond.opt(X=X,Y=Y,Xpred=Xpred,Ypred=Ypred,max.resX=10,max.resY=10,rho0.X=0.5,rho0.Y=0.5)
print(ans.cond.opt$part_points_hmap) # the hMAP partition
print(sum(log(ans.cond.opt$predictive_densities))) # log predictive score for the testing set
