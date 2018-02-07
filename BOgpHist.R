####
# Bayesian optimisation using gpHist...
##
options(error = recover)
options(warnings = recover)

## Used for illustration. Only works in 2D and 3D, set to FALSE if not needed
PLOT=TRUE

#used for 3d illustration
library(rgl)
require(akima)

##requires the gpHist package from the gitHub repro
require(gpHist)
## used for optimisation of the surrogate function
require(rgenoud)

#getY function to be optimized
#npoints number of points used in the initial grid
#ndim number of dimensions of the optimisation function
#lower optimisation function lower bound
#upper optimisation function upper bound
#paramLower array of lower bound for parameter required for the transformation function
#paramUpper array of upper bound for parameter required for the transformation function
#datatransform function that transforms the data. 
#ntransform number of parameters required for the transformation function
#it number of iterations
#acq name of the acquisition function ('ei' || 'ucb' || 'poi')
#kappa kappa parameter for acquisition function (ucb)
#eps  eps parameter for the acquisition function
#nEigenVals percentage of eigenvalues used to approximate the variance of the gpHist used eigenvalues = (number of examples/nEigenVals)



BOgpHist =function(getY,npoints,ndim,lower,upper,paramLower,paramUpper,datatransform=NULL, ntransform=0, it=100,acq = "ucb",kappa= 2.576, eps=0.0,nEigenVals=1){
  evaluations = 0;
  improvement = c()
  
  if(acq!= 'ei'&& acq != 'ucb' && acq != 'poi' ){
    print('Not valid aquisition function!')
    return (0)
  }
  ##create 2d grid
  #npoints_c = ndim^2
  if (ndim == 2){
    len = npoints
    print('create grid')
    X = expand.grid(seq(lower,upper,length=len), seq(lower,upper,length=len)) 
    X = as.matrix(X)
  }else if(ndim==1){
    X =as.matrix( seq(lower,upper,length=npoints))
  }else if(ndim>2){
    expand.grids <- function(x,d){
      expand.grid(replicate(d, x, simplify=FALSE))
    }
    X = expand.grids(  seq(lower,upper,length=npoints),ndim)
    X = as.matrix(X)
  }
  cat (c('init grid with :',nrow(X),'points\n') )
  Y = matrix(,ncol=1,nrow=nrow(X)) ##create Y
  
  for(i in 1:nrow(X)){
    Y[i] = getY(X[i,])
  }
  evaluations = evaluations + nrow(X)
  ymin = min(Y)
   

  for( i in 1:it){
    ##estimates hyper parameters. the number of eigenvectors is set to 1. this is the fastest but also most inaccurate setting
    # consider to chose a different parameter for this. potentially the parameter chosen for the variance approximation
    start = proc.time() ##used for time measure
    
    resmat = estimateHyperParameters(X, Y, paramLower=paramLower,paramUpper=paramUpper,datatransform=datatransform,nParams=ntransform,it = 200,tol=0.0001,k=1);
  
    end = proc.time()
    print('time for estimating hyperparameters:')
    print(end-start)
    
    start = proc.time() ##used for time measure
    
    px = resmat[2:length(resmat)]
    if(is.null(datatransform)){
      GP = gpHist(X= X,Y=Y,sigma = resmat[1],k=ceiling(nrow(X)/nEigenVals) )
    }else{
      X_trans = datatransform( X,px)
      GP = gpHist(X= X_trans,Y=Y,sigma = resmat[1],k=ceiling(nrow(X)/nEigenVals) )
    }
    end = proc.time()
    print('time to fit GP:')
    print(end-start)
    
    if(is.na(GP$logmarginal)){ ## this can happen but should not happen
      print('log marginal contains NA')
      stop()
    }
    
    
    start = proc.time()
    utmat = Utility_Max(X_trans,GP,ymin,lower,upper,datatransform = datatransform,px=px,kappa = kappa, eps = eps,acq=acq) ##maybe more than 10
    end = proc.time()
    print('time for optimizing acq function:')
    print(end-start)
    
    ###this is time consuming and only used for illustration
    if(PLOT & ndim==1){
      
      ##plot true function
      truePointsX =seq(lower,upper,0.05)
      truePointsY =c() 
      for(k in 1:length(truePointsX)){
        truePointsY = c(truePointsY,getY(truePointsX[k]))
      }
      
      plot(truePointsX,truePointsY,type='l',col='black',xlab='x', ylab='y',lwd=2)#,ylim=c(-100 ,100) )
      points(X,Y,type='o',lty=0,cex=2)
      
      title(i)
      predPoints = matrix(seq(lower,upper,0.1))
      
      if(is.null(datatransform)){
        x_pred_trans = predPoints
        X_trans = X
      }else{
        x_pred_trans = datatransform(predPoints ,px)
      }
      ##plot stuff  
      preds = gpHistPredict(GP=GP,X=X_trans,x_pred = x_pred_trans)
      lines(predPoints,preds,col='red',lwd=2)
      vars =  gpHistVariance(GP=GP,X=X_trans,x_pred = x_pred_trans)
      vars = abs(vars)
      
      lines(predPoints,preds-sqrt(vars),col='red',lty=2 ,lwd=2)
      lines(predPoints,preds+sqrt(vars),col='red',lty=2 ,lwd=2)
      
      
      eivals = c()
      for(w in 1:length(x_pred_trans)){
        ret = Utility(x_pred_trans[w], X_trans,GP,ymin, acq = acq , kappa = kappa, eps=eps,datatransform=NULL);
        eivals =c(eivals,ret)
      }
      if(acq=='ei'){
        lines(predPoints,-eivals,col='blue')
      }else{
        lines(predPoints,eivals,col='blue')
      }
     
      ###plot next point to evaluate
      for( w in 1:nrow(utmat)){
        points(utmat[w,1],utmat[w,2],pch='+',col='black')
      }
      
     legend("topleft", legend = c("True function ", "Approximation","Standard deviation",'Upper confidence bound','Sample points'),
             text.width = strwidth("1,000,000"),
            lty = c(1,1,2,1,NA), pch= c(NA,NA,NA,NA,'o'),col=c('black','red','red' ,'blue'), xjust = 1, yjust = 1,cex=1.0)
      
      #####
    }else if(PLOT & ndim==2){
     # plot3d(X[,1],X[,2], Y, col="red", size=5,xlab='X',ylab='Y',zlab='f(X,Y)') 
      grid = 0.1
      pd = (upper-lower)/ grid +1
      points = matrix(,nrow=pd^2, ncol=2)
      rowidx = 1
      for(xi in seq(lower,upper, grid)){
        for(yi in seq(lower,upper, grid)){
          points[rowidx,] = c(xi,yi)
          rowidx = rowidx + 1
        }
      }
      
      
      if(is.null(datatransform)){
        x_pred_trans = points
        X_trans = X
      }else{
        x_pred_trans = datatransform(points ,px)
      }
      
      ##plot stuff  
      preds = gpHistPredict(GP=GP,X=X_trans,x_pred = x_pred_trans)
      
      tp=c()
      for(i in 1:nrow(points)){
        tp[i] = getY(points[i,])
      }
      plot3d(points[,1],points[,2], preds, col="red", size=5,xlab='X',ylab='Y',zlab='f(X,Y)') 
      plot3d(points[,1],points[,2], tp, col="blue", size=5,add=T) 
      
      s=interp(points[,1], points[,2],preds)
      k=interp(points[,1], points[,2],tp)
      
      persp3d(s$x,s$y,s$z, aspect=c(10, 1, 0.5),  col = "red", add=TRUE )
      persp3d(k$x,k$y,k$z, aspect=c(10, 1, 0.5),  col = "lightblue", add=TRUE )
      plot3d(X[,1],X[,2], Y, col="green", size=5,add=T) 
      
      vars =  gpHistVariance(GP=GP,X=X_trans,x_pred = x_pred_trans)
      if(any(vars<0) ){
        print('negative variances')
        vars[vars<0]=0
      }
      v=interp(points[,1], points[,2],preds-sqrt(vars) )
      persp3d(v$x,v$y,v$z, aspect=c(10, 1, 0.5),  col = "yellow", add=TRUE )
      
      legend3d("topright", legend = c('Optimisation function','HIK approximation','Standard deviation'), pch = 16, col = c('blue','red','yellow'), cex=2, inset=c(0.02))
      
      ##add 2d
      n.grid <- 20
      x.grid <- y.grid <- seq(0,1,length=n.grid)
      design.grid <- expand.grid(x.grid, y.grid)
      response.grid <- apply(design.grid, 1, getY)
      z.grid <- matrix(response.grid, n.grid, n.grid)
      contour(x.grid, y.grid, z.grid, 40)
      title("objective function")
      points(X[,1],X[,2], pch=17, col="blue") ## known points
      
    }
    
    
    ##select the minimum of all estimated points of the acquisition function
    # this is tedious, because we make sure that this point is not already known.
    # this is debatable if this can or should still happen, it is certainly time consuming

    start = proc.time()
    
      ##querry points
    newX = matrix(utmat[1,1:ndim],nrow=1)
      
    rowcheck  <- function(df1, df2){
        xx <- apply(df1, 1, paste, collapse = "")
        yy <- apply(df2, 1, paste, collapse = "")
        zz <- xx %in% yy
        return(zz)
    }
    xidx = rowcheck(X,newX)
      
    if(any(xidx) ){
      print('no new point') #we can update now or leave
      # it is debatable if the hyper-parameters should be estimated again.
      next;
    }
    
    end = proc.time()
    print('time to select next sample point:')
    print(end-start)
    
    newX = matrix(newX,nrow=1)
    newY = getY(newX ) ###get real Y
   
    evaluations = evaluations +1 # add 1
  
    ##add new points  
    if(PLOT& ndim==1){
      points(newX,newY,col='blue')
    }else if (PLOT & ndim ==2){
      points(newX[1],newX[2], pch=17, col="red") ## new point
    }
    
     X= rbind(X, newX) ## rbind is slow, but we add the new sample
     Y= rbind(Y,newY )
     
    improvement = c(improvement, newY )
    if(newY < ymin){ ##improvement
      ymin = newY
    }
    ##here one could add checks if improvement subsides or other criteria to stop the process
     
  }
  
  yminpos = which.min(Y)
  xminpos = X[yminpos,]
  
  list('X'=X,'Y'=Y,'Xpos'=xminpos,'Ymin'=Y[yminpos],'GP'=GP,'resmat'=resmat,'X_trans'=X_trans,'evaluations'=evaluations,'improvement'=improvement)
}

##########UTILITY FUNCTIONS

Utility <- function(x_vec, X,GP,y_max, acq = "ucb", kappa, eps,datatransform=NULL,px=NULL) {
  
  if(!is.null(datatransform)){
    testpoint = matrix(datatransform(matrix(x_vec,nrow=1),px),nrow=1)  
  }else{
    testpoint = matrix(x_vec,nrow=1)  
  }
  
  GP_Mean <- gpHistPredict(GP ,X,testpoint ) # 
  GP_MSE <- gpHistVariance(GP,X,testpoint)#
  if(is.na(GP_Mean)){
    print('MEAN is NA' )
    GP_Mean =0
  }
  if(GP_MSE<0){ ##to avoid errors!
    print('force lower approximation')
   ##use a slightly worse approximation
    N_eigenvals = length(GP$lambda)
    while(GP_MSE<0 & N_eigenvals>1 ){
      N_eigenvals = N_eigenvals -1
      GP$lambda = GP$lambda[1:N_eigenvals]
      GP$vectors = GP$vectors[,1:N_eigenvals]
      GP_MSE <- gpHistVariance(GP,X,testpoint)
    }

  }
  if(is.na(GP_MSE)){
    print('GP_MSE IS NA !!')
    GP_MSE = 0
  }

  # Utility Function Type
  if (acq == "ucb") {
    squ = sqrt(GP_MSE)
    if(is.na(squ)){
      print('SQRT(GP_MSE) IS NA !!')
      squ = 0
    }
    value <- GP_Mean - kappa * squ

  } else if (acq == "ei") {
    z <- (y_max- GP_Mean  - eps) / sqrt(GP_MSE)
    value <-  (y_max- GP_Mean- eps) * pnorm(z) + sqrt(GP_MSE) * dnorm(z)
    value = value *-1;
    if(is.na(value)){
      print('calulated value is NA !?')
      value = 0
    }
  } else if (acq == "poi") {
    z <- (y_max-GP_Mean - eps) / sqrt(GP_MSE)
    value <- -pnorm(z)
  }
  
  if(is.na(value) ){
    print('value == NA')
    value = 0
  }
  return(value)
}

###try to find next value...
Utility_Max <- function(X,GP,  y_max,lower=0,upper=1,acq=  "ucb", datatransform=NULL,px=NULL,kappa= 2.576, eps=0.0) {
  ndim = ncol(X)
  if(length(lower)==1){
    if(ndim>1){
      lower = rep(lower,ndim )
      upper = rep(upper,ndim )
    }
  }

  #only limited amount of generations is used as well as population size these parameters are debatable  
  res = genoud(Utility, ndim,Domains= cbind(lower,upper), max.generations=50 ,hard.generation.limit=T,pop.size=40, BFGSburnin=2,boundary.enforcement=2,X= X, GP = GP, acq = acq, y_max = y_max,datatransform=datatransform,px=px, kappa = kappa, eps = eps)
  cbind(t(res$par),res$value)
}



## Rastrigin test function
fctest = function(p){
  A=10
  ret = A + p^2-A*cos(2*pi*p)
}

###transformation
datatrans = function(X,p){
  if(!is.matrix(X)){
    X = matrix(X,nrow=1)  ###one sample
  }
  X = (X+10+p[2]) ## aplly a shift into positive range and estimate an additional shift from the data
  X = (( (exp( p[1]*abs(X) )-1 ) / (exp(p[1])-1 ))  ) # use exponential transformation
  X
}

### utilzed parameters:
# testfunction, function to be optimized 
#6 initial points,
#1 dimension
# range between -10 and 10
# first parameter is for sigma - sigma between 0.0001-0.1 , exponential parameter, and additional shift
#upper range 
#number of parameters for the data transformation function
# number of iterations.
# kappa parameter not used for ei
# use ei acquisition function
# number of eigenvalues used for variance approximation
# eps parameter set to 0
res = BOgpHist(fctest,6,1,-10,10,c(0.0001,0.0001,0),c(0.1,1,100),datatransform=datatrans, ntransform=2,it=10,kappa=2.5,acq='ei',nEigenVals=4,eps=0.0)



####
datatrans = function(X,p){
  if(!is.matrix(X)){
    X = matrix(X,nrow=1)  ###one sample
  }
  X = (X+p[2]) ## estimate an additional shift from the data
  X = (( (exp( p[1]*abs(X) )-1 ) / (exp(p[1])-1 ))  ) # use exponential transformation
  X
}


branin = function (x){
  x1 <- x[1] * 15 - 5
  x2 <- x[2] * 15
  (x2 - 5/(4 * pi^2) * (x1^2) + 5/pi * x1 - 6)^2 + 10 * (1 - 
                                                           1/(8 * pi)) * cos(x1) + 10
}


### in 2D
res = BOgpHist(getY=branin,npoints=3,ndim=2,lower=0,upper=1,c(0.0001,0.0001,0,0.0001,0),c(0.01,1,10,1,10),datatransform=datatrans, ntransform=4,it=19,kappa=5.576,acq='ucb',nEigenVals=4,eps=0.0)

