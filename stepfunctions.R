

####################################### Bunch of random functions ###############################################################

# Scale values between 0 and 1
scaler <- function(x){(x-min(x))/(max(x)-min(x))}

# Returns the closest value in supplied vector to a given one: f(given_vect, value). 
closest <- function(x, your){
  res <-c()
  for (num in 1:length(your)){ 
    res <- c(res, x[which(abs(x-your[num])==min(abs(x-your[num])))])
  }
  return(res)
}

# Calculate binomial coefficient for large numbers. Note the function does not return exp()
ramanujan <- function(n){
  n*log(n) - n + log(n*(1 + 4*n*(1+2*n)))/6 + log(pi)/2
}

bignchoosek <- function(n,k){
  #exp(ramanujan(n) - ramanujan(k) - ramanujan(n-k))
  if (k == 0 | n == k){
    return (log(choose(n,k)))
  }
  else{
    return (ramanujan(n) - ramanujan(k) - ramanujan(n-k))
  }
}

createTopo <- function(la1, la2, lo1, lo2){
  
  #Import topography data
  data(ETOPO5)
  
  # Need to have 2 ranges because of negative coordinates (-10 = 350)
  glat1 = c(la1, la2)
  glon1 = c(0, lo2)
  
  glat2 = c(la1, la2)
  glon2 = c(360+lo1, 359.99)
  
  # Extract topology data for lon and lat ranges specified above
  topo1 = getETOPO(ETOPO5, glat1, glon1)
  topo2 = getETOPO(ETOPO5, glat2, glon2)
  
  resolution <- 11#(1 - highest, 11 - lowest)
  
  # Extract coordinates for the points corresponding to topology data and reduce resolution (default 12pts per degree)
  topoLon1 <- attr(topo1, "lon")[c(TRUE,rep(FALSE,resolution))]
  topoLon2 <- attr(topo2, "lon")[c(TRUE,rep(FALSE,resolution))]-360 # Change from for example 340 to -10
  topoLat <- attr(topo2, "lat")[c(TRUE,rep(FALSE,resolution))]
  topoLon <- c(topoLon2, topoLon1)
  
  topo <- rbind(topo2, topo1)
  topo <- t(topo[, ncol(topo):1])
  
  topo <- topo[c(TRUE,rep(FALSE,resolution)), c(TRUE,rep(FALSE,resolution))] # Extract every 12th point
  returns <- list(topo=topo, topoLon=topoLon, topoLat=topoLat)
  return(returns)
}

FischerSolver <- function(Nx, Ny, s, sigmax, sigmay, D, x0, y0, xtran, ytran, timePoints, oldnew=NULL, initMat=NULL){
  
  # Setup a grid
  x.grid <- setup.grid.1D(x.up = la1, x.down = la2, N = Nx)
  y.grid <- setup.grid.1D(x.up = lo1, x.down = lo2, N = Ny)
  
  grid2D <- setup.grid.2D(x.grid, y.grid)
  
  distx <- x.grid$dx[1] # Distance between two points in x direction on a grid in degrees
  disty <- y.grid$dx[1] # Distance between two points in y direction on a grid in degrees
  
  distlon <- sapply(x.grid$x.mid,function(lat){
    return(distHaversine(c(0,lat),c(disty,lat))/1000)
  })
  distlon <- matrix(rep(distlon[length(distlon):1],Ny+1),nrow=Nx,ncol=Ny+1)
  
  distlat <- distHaversine(c(0,0),c(0,distx))/1000 # Distance between two latitudes in km
  distlat <- matrix(distlat, nrow=Nx+1, ncol=Ny)
  
  # Initial frequency at a point
  p0 = 1/(2*D*mean(distlat)*distlon[x0])
  
  D.grid <- list()
  D.grid$x.int <- 0.5*(distx*(sigmay/distlat))^2
  D.grid$x.int <- D.grid$x.int[,ncol(D.grid$x.int):1]
  D.grid$y.int <- 0.5*(disty*(sigmax/distlon))^2
  
  # Add a land bridge to GB
  gb_lon1 <- which(topoLon == closest(topoLon, -6))
  gb_lon2 <- which(topoLon == closest(topoLon, 1))
  gb_lat1 <- which(rev(topoLat) == closest(rev(topoLat), 49))
  gb_lat2 <- which(rev(topoLat) == closest(rev(topoLat), 52))
  topo[gb_lat1:gb_lat2,gb_lon1:gb_lon2] <- 1
  topo[gb_lat1:gb_lat2,gb_lon1:gb_lon2] <- 1
  
  # Add a land bridge to Sardinia
  sard_lon1 <- which(topoLon == closest(topoLon, 9))
  sard_lon2 <- which(topoLon == closest(topoLon, 11))
  sard_lat1 <- which(rev(topoLat) == closest(rev(topoLat), 41))
  sard_lat2 <- which(rev(topoLat) == closest(rev(topoLat), 42))
  topo[sard_lat1:sard_lat2,sard_lon1:sard_lon2] <- 1
  topo[sard_lat1:sard_lat2,sard_lon1:sard_lon2] <- 1
  
  # Assign 0 deffusion where topo is negative
  D.grid$x.int[1:nrow(topo),1:ncol(topo)][which(topo < 0)] <- 0
  D.grid$y.int[1:nrow(topo),1:ncol(topo)][which(topo < 0)] <- 0
  
  # Define advection grid
  v.grid <- list()
  v.grid$x.int <- (distx*(ytran/distlat))
  v.grid$x.int <- v.grid$x.int[,ncol(v.grid$x.int):1]
  v.grid$y.int <- (disty*(xtran/distlon))
  
  # Change in allele frequency due to selection (reaction term)
  deltaP <- function(p){
    res <- p*(1-p)*(p*2*s+s*(1-2*p))*10
    return(res)
  }
  
  # Define Reaction-diffusion function
  ReacDiff <- function(t, y, parms) {
    
    X1 <- matrix(nrow = Nx, ncol = Ny, data = y[1:(Nx*Ny)])
    
    dX1 <- deltaP(X1)+10*
      tran.2D (C = X1, D.grid = D.grid, dx = distx, dy = disty, grid=grid2D, v.grid=v.grid)$dC
    list(c(dX1)) }
  
  # Initial conditions specify allele freq=p0 at a single point, rest are 0
  # Different initial conditions depending time period before and after 5000 years BP
  if (oldnew == "old") {
    yini <- matrix(rep(0,Nx*Ny), ncol=Ny, nrow=Nx)#+1e-5
    yini[x0,y0] <- p0#-1e-5
  } else if (oldnew == "new") {yini <- initMat}
  
  # This is where the numeric calculations actually happen
  times <- 0:timePoints
  out <- ode.2D(y = yini, parms = NULL, func = ReacDiff, nspec = 1, dimens = c(Nx, Ny), times = times,
                lrw = 10000000)
  return(out)
}

getCI <- function(bestRes, oldnew, LL){
  
  flag <-  TRUE
  maxit <- 0
  estim <- c(coef(bestRes)[["s"]], coef(bestRes)[["sigmax"]], coef(bestRes)[["sigmay"]], coef(bestRes)[["xtran"]], coef(bestRes)[["ytran"]])
  
  while (flag & maxit != 100) {
    cat("Iteration:", maxit, "\n")
    
    lower <- c(0.001, 1, 1, -2.5, -2.5)
    upper <- c(0.1, 100, 100, 2.5, 2.5)
    
    set.seed(66)
    start_list <- list(s=estim[1], sigmax=estim[2], sigmay=estim[3], xtran=estim[4], ytran=estim[5])
    start_list_new <- mapply(x = start_list, y = rnorm(5,0,c(0.002, 1, 1, 0.05, 0.05)), 
                             function(x, y) x<-x+y)
    
    bool_v <- mapply(x = start_list_new,u = upper,l = lower,FUN = function(x, u, l) x <= u & x >= l)
    
    start_list_new[which(bool_v == FALSE)] <- start_list[which(bool_v == FALSE)] 
    start_list <- start_list_new
    
    if (oldnew == "old") {
      resCI <- mle2(minuslogl=LL, start=start_list,
                    fixed=list(Nx=Nx, Ny=Ny, D=D, x0=coef(bestRes)[["x0"]], y0=coef(bestRes)[["y0"]], timeP=timeP, timeSplit=timeSplit, oldnew=0),
                    method = "L-BFGS-B")
      
      estim <- c(coef(resCI)[["s"]], coef(resCI)[["sigmax"]], coef(resCI)[["sigmay"]], coef(resCI)[["xtran"]], coef(resCI)[["ytran"]])
      
      hes <- numDeriv::hessian(function(x) LL(1, Nx, Ny, x[1], x[2], x[3], D, coef(bestRes)[["x0"]], coef(bestRes)[["y0"]], x[4], x[5], timeP, timeSplit, 0), estim)
      hes2 <- nlme::fdHess(estim, function(x) LL(1, Nx, Ny, x[1], x[2], x[3], D, coef(bestRes)[["x0"]], coef(bestRes)[["y0"]], x[4], x[5], timeP, timeSplit, 0))
      
    } else if (oldnew == "new"){
      resCI <- mle2(minuslogl=LL, start=start_list,
                    fixed=list(Nx=Nx, Ny=Ny, D=D, x0=coef(bestRes)[["x0"]], y0=coef(bestRes)[["y0"]], timeP=timeP, timeSplit=timeSplit, oldnew=1),
                    method = "L-BFGS-B")
      
      estim <- c(coef(resCI)[["s"]], coef(resCI)[["sigmax"]], coef(resCI)[["sigmay"]], coef(resCI)[["xtran"]], coef(resCI)[["ytran"]])
      
      hes <- numDeriv::hessian(function(x) LL(1, Nx, Ny, x[1], x[2], x[3], D, coef(bestRes)[["x0"]], coef(bestRes)[["y0"]], x[4], x[5], timeP, timeSplit, 1), estim)
      hes2 <- nlme::fdHess(estim, function(x) LL(1, Nx, Ny, x[1], x[2], x[3], D, coef(bestRes)[["x0"]], coef(bestRes)[["y0"]], x[4], x[5], timeP, timeSplit, 1))
      
    }
    inv_hes <- solve(hes)
    inv_hes2 <- solve(hes2$Hessian)
    
    maxit <- maxit + 1
    if (all(diag(inv_hes)>=0) | all(diag(inv_hes2)>=0)){
      CIvect <- c()
      flag <- FALSE
      ifelse(all(diag(inv_hes)>=0), inv_hes <- inv_hes, inv_hes <- inv_hes2)
      for (i in 1:length(estim)){
        CIvect <- c(CIvect, round(estim[i]-1.96*sqrt(inv_hes[i,i]), 4), round(estim[i]+1.96*sqrt(inv_hes[i,i]), 4))
      }}
  }
  return(c(resCI, CIvect))
}

param_optim_old <- function(init_points = 2, seed = 5, num_cores = 1, LL){

  # Randomly initialise grid points for optimization
  set.seed(seed)
  pars <- data.frame(min=c(0,1,1,-2.5,-2.5,1), max=c(0.1,100,100,2.5,2.5,28))
  rownames(pars) <- c("s", "sigmax", "sigmay", "xtran", "ytran", "place")
  parSets <- Latinhyper(pars, init_points)
  parSets <- cbind(parSets, sapply(parSets[,"place"], function(x) Nx+1-which(round(topoLat) == Latpoints[round(x)])[1]),
                   sapply(parSets[,"place"], function(x) which(round(topoLon) == Lonpoints[x])[1]))
  parSets <- parSets[, !colnames(parSets) %in% c("place")] 
  colnames(parSets) <- c("s", "sigmax", "sigmay", "xtran", "ytran", "x0", "y0")
  
  samps = 1:init_points
  set.seed(seed, kind = "L'Ecuyer-CMRG")
  #Do the optimization and return best set of parameters
  res <- mclapply(samps, function(x) mle2(minuslogl=LL, start=as.list(parSets[x,]),
                                          fixed=list(Nx=Nx, Ny=Ny, D=D, timeP=timeP, timeSplit=timeSplit, oldnew=0),
                                          method = "SANN", control = list(temp = 100)), mc.cores=num_cores)

  tres <- unlist(lapply(1:length(res), function(i) logLik(res[[i]])))
  bestRes <- res[[which(tres == max(tres, na.rm=TRUE))]]
  
  return(bestRes)
}

param_optim_new <- function(init_points = 2, seed = 5, num_cores = 1, LL){
  # Randomly initialise grid points for optimization
  set.seed(seed)
  pars <- data.frame(min=c(0,1,1,-2.5,-2.5), max=c(0.1,100,100,2.5,2.5))
  rownames(pars) <- c("s", "sigmax", "sigmay", "xtran", "ytran")
  parSets <- Latinhyper(pars, init_points)
  
  samps = 1:init_points
  
  set.seed(seed, kind = "L'Ecuyer-CMRG")
  # Do the optimization and return best set of parameters
  res <- mclapply(samps, function(x) mle2(minuslogl=LL, start=as.list(parSets[x,]),
                                          fixed=list(Nx=Nx, Ny=Ny, D=D, x0=coef(bestResOld)[["x0"]], y0=coef(bestResOld)[["y0"]],
                                                     timeP=timeP, timeSplit=timeSplit, oldnew=1),
                                          method = "SANN", control = list(temp = 100)), mc.cores=num_cores)

  tres <- unlist(lapply(1:length(res), function(i) logLik(res[[i]])))
  bestRes <- res[[which(tres == max(tres, na.rm=TRUE))]]
  
  return(bestRes)
}

param_optim_fixedInit <- function(init_points = 2, seed = 5, num_cores = 1, LL, oldnew){
  
  # Randomly initialise grid points for optimization
  set.seed(seed)
  pars <- data.frame(min=c(0,1,1,-2.5,-2.5), max=c(0.1,100,100,2.5,2.5))
  rownames(pars) <- c("s", "sigmax", "sigmay", "xtran", "ytran")
  parSets <- Latinhyper(pars, init_points)
  
  colnames(parSets) <- c("s", "sigmax", "sigmay", "xtran", "ytran")
  
  samps = 1:init_points
  
  # Do the optimization and return best set of parameters
  if (oldnew == "old") {
    set.seed(seed, kind = "L'Ecuyer-CMRG")
    res <- mclapply(samps, function(x) mle2(minuslogl=LL, start=as.list(parSets[x,]),
                                            fixed=list(Nx=Nx, Ny=Ny, D=D,  
                                                       x0=x0_fixed, y0=y0_fixed,
                                                       timeP=timeP, timeSplit=timeSplit, oldnew=0),
                                            method = "SANN", control = list(temp = 100)), mc.cores=num_cores)
  } else if (oldnew == "new"){ 
    set.seed(seed, kind = "L'Ecuyer-CMRG")
    res <- mclapply(samps, function(x) mle2(minuslogl=LL, start=as.list(parSets[x,]),
                                            fixed=list(Nx=Nx, Ny=Ny, D=D,  
                                                       x0=coef(bestResOld)[["x0"]], y0=coef(bestResOld)[["y0"]],
                                                       timeP=timeP, timeSplit=timeSplit, oldnew=1),
                                            method = "SANN", control = list(temp = 100)), mc.cores=num_cores)
  }
  tres <- unlist(lapply(1:length(res), function(i) logLik(res[[i]])))
  bestRes <- res[[which(tres == max(tres, na.rm=TRUE))]]
  
  return(bestRes)
}


getCI_B <- function(bestRes, advection=TRUE, oldnew, LL){
  
  flag <-  TRUE
  maxit <- 0
  estim <- c(coef(bestRes)[["s"]], coef(bestRes)[["sigmax"]], coef(bestRes)[["sigmay"]], coef(bestRes)[["xtran"]], coef(bestRes)[["ytran"]])
  
  while (flag & maxit != 100) {
    cat("Iteration:", maxit, "\n")
    
    if (advection == TRUE){
      lower <- c(0.001, 1, 1, -2.5, -2.5)
      upper <- c(0.1, 100, 100, 2.5, 2.5)
    
      set.seed(66)
      start_list <- list(s=estim[1], sigmax=estim[2], sigmay=estim[3], xtran=estim[4], ytran=estim[5])
      start_list_new <- mapply(x = start_list, y = rnorm(5,0,c(0.002, 1, 1, 0.05, 0.05)), 
                               function(x, y) x<-x+y)
      fixed = list(Nx=Nx, Ny=Ny, D=D, x0=coef(bestRes)[["x0"]], y0=coef(bestRes)[["y0"]], timeP=timeP, timeSplit=timeSplit, oldnew=oldnew)
    } else{      
      lower <- c(0.001, 1, 1)
      upper <- c(0.1, 100, 100)
      set.seed(66)
      start_list <- list(s=estim[1], sigmax=estim[2], sigmay=estim[3])
      start_list_new <- mapply(x = start_list, y = rnorm(3,0,c(0.002, 1, 1)), 
                               function(x, y) x<-x+y)
      fixed = list(Nx=Nx, Ny=Ny, D=D, x0=coef(bestRes)[["x0"]], y0=coef(bestRes)[["y0"]], xtran = 0, ytran = 0, timeP=timeP, timeSplit=timeSplit, oldnew=oldnew)
    }      
    
    bool_v <- mapply(x = start_list_new,u = upper,l = lower,FUN = function(x, u, l) x <= u & x >= l)
    
    start_list_new[which(bool_v == FALSE)] <- start_list[which(bool_v == FALSE)] 
    start_list <- start_list_new
    
    resCI <- mle2(minuslogl=LL, start=start_list, lower=lower, upper=upper,
                  fixed=fixed,
                  method = "L-BFGS-B")
    if (advection == TRUE){
      estim <- c(coef(resCI)[["s"]], coef(resCI)[["sigmax"]], coef(resCI)[["sigmay"]], coef(resCI)[["xtran"]], coef(resCI)[["ytran"]])
      hes <- numDeriv::hessian(function(x) LL(1, Nx, Ny, x[1], x[2], x[3], D, coef(bestRes)[["x0"]], coef(bestRes)[["y0"]], x[4], x[5], timeP, timeSplit, oldnew), estim)
      hes2 <- nlme::fdHess(estim, function(x) LL(1, Nx, Ny, x[1], x[2], x[3], D, coef(bestRes)[["x0"]], coef(bestRes)[["y0"]], x[4], x[5], timeP, timeSplit, oldnew))
    } else{
      estim <- c(coef(resCI)[["s"]], coef(resCI)[["sigmax"]], coef(resCI)[["sigmay"]])
      hes <- numDeriv::hessian(function(x) LL(1, Nx, Ny, x[1], x[2], x[3], D, coef(bestRes)[["x0"]], coef(bestRes)[["y0"]], 0, 0, timeP, timeSplit, oldnew), estim)
      hes2 <- nlme::fdHess(estim, function(x) LL(1, Nx, Ny, x[1], x[2], x[3], D, coef(bestRes)[["x0"]], coef(bestRes)[["y0"]], 0, 0, timeP, timeSplit, oldnew))}
    
    inv_hes <- solve(hes)
    inv_hes2 <- solve(hes2$Hessian)
    
    maxit <- maxit + 1
    if (all(diag(inv_hes)>=0) | all(diag(inv_hes2)>=0)){
      CIvect <- c()
      flag <- FALSE
      ifelse(all(diag(inv_hes)>=0), inv_hes <- inv_hes, inv_hes <- inv_hes2)
      for (i in 1:length(estim)){
        CIvect <- c(CIvect, round(estim[i]-1.96*sqrt(inv_hes[i,i]), 4), round(estim[i]+1.96*sqrt(inv_hes[i,i]), 4))
      }}
  }
  return(c(resCI, CIvect))
}

param_optim_old_B <- function(init_points = 2, seed = 5, num_cores = 1, advection = TRUE, LL){
  
  # Randomly initialise grid points for optimization
  if (advection == TRUE){
    set.seed(seed)
    pars <- data.frame(min=c(0,1,1,-2.5,-2.5,1), max=c(0.1,100,100,2.5,2.5,28))
    rownames(pars) <- c("s", "sigmax", "sigmay", "xtran", "ytran", "place")
    parSets <- Latinhyper(pars, init_points)
    parSets <- cbind(parSets, sapply(parSets[,"place"], function(x) Nx+1-which(round(topoLat) == Latpoints[round(x)])[1]),
                     sapply(parSets[,"place"], function(x) which(round(topoLon) == Lonpoints[x])[1]))
    parSets <- parSets[, !colnames(parSets) %in% c("place")] 
    colnames(parSets) <- c("s", "sigmax", "sigmay", "xtran", "ytran", "x0", "y0")
    
    samps = 1:init_points
    
    set.seed(seed, kind = "L'Ecuyer-CMRG")
    #Do the optimization and return best set of parameters
    res <- mclapply(samps, function(x) mle2(minuslogl=LL, start=as.list(parSets[x,]),
                                            fixed=list(Nx=Nx, Ny=Ny, D=D, timeP=timeP, timeSplit=timeSplit, oldnew=0),
                                            method = "SANN", control = list(temp = 100)), mc.cores=num_cores)
  } else if (advection == FALSE){ 
    set.seed(seed)
    pars <- data.frame(min=c(0,1,1,1), max=c(0.1,100,100,28))
    rownames(pars) <- c("s", "sigmax", "sigmay", "place")
    parSets <- Latinhyper(pars, init_points)
    parSets <- cbind(parSets, sapply(parSets[,"place"], function(x) Nx+1-which(round(topoLat) == Latpoints[round(x)])[1]),
                     sapply(parSets[,"place"], function(x) which(round(topoLon) == Lonpoints[x])[1]))
    parSets <- parSets[, !colnames(parSets) %in% c("place")] 
    colnames(parSets) <- c("s", "sigmax", "sigmay", "x0", "y0")
    samps = 1:init_points
    set.seed(seed, kind = "L'Ecuyer-CMRG")
    res <- mclapply(samps, function(x) mle2(minuslogl=LL, start=as.list(parSets[x,]),
                                            fixed=list(Nx=Nx, Ny=Ny, D=D, xtran=0, ytran=0, timeP=timeP, timeSplit=timeSplit, oldnew=0),
                                            method = "SANN", control = list(temp = 100)), mc.cores=num_cores)}
  tres <- unlist(lapply(1:length(res), function(i) logLik(res[[i]])))
  bestRes <- res[[which(tres == max(tres, na.rm=TRUE))]]
  
  return(bestRes)
}

param_optim_new_B <- function(init_points = 2, seed = 5, num_cores = 1, advection = TRUE, LL){
  if (advection == TRUE){
  # Randomly initialise grid points for optimization
    set.seed(seed)
    pars <- data.frame(min=c(0,1,1,-2.5,-2.5), max=c(0.1,100,100,2.5,2.5))
    rownames(pars) <- c("s", "sigmax", "sigmay", "xtran", "ytran")
    parSets <- Latinhyper(pars, init_points)
    samps = 1:init_points
    set.seed(seed, kind = "L'Ecuyer-CMRG")
    # Do the optimization and return best set of parameters
    res <- mclapply(samps, function(x) mle2(minuslogl=LL, start=as.list(parSets[x,]),
                                            fixed=list(Nx=Nx, Ny=Ny, D=D, x0=coef(bestResOld)[["x0"]], y0=coef(bestResOld)[["y0"]],
                                                       timeP=timeP, timeSplit=timeSplit, oldnew=1),
                                            method = "SANN", control = list(temp = 100)), mc.cores=num_cores)
  } else if (advection == FALSE){    
    set.seed(seed)
    pars <- data.frame(min=c(0,1,1), max=c(0.1,100,100))
    rownames(pars) <- c("s", "sigmax", "sigmay")
    parSets <- Latinhyper(pars, init_points)
    samps = 1:init_points
    set.seed(seed, kind = "L'Ecuyer-CMRG")
    res <- mclapply(samps, function(x) mle2(minuslogl=LL, start=as.list(parSets[x,]),
                                            fixed=list(Nx=Nx, Ny=Ny, D=D, x0=coef(bestResOld)[["x0"]], y0=coef(bestResOld)[["y0"]],
                                                       xtran=0, ytran=0, timeP=timeP, timeSplit=timeSplit, oldnew=1),
                                            method = "SANN", control = list(temp = 100)), mc.cores=num_cores)}   
  tres <- unlist(lapply(1:length(res), function(i) logLik(res[[i]])))
  bestRes <- res[[which(tres == max(tres, na.rm=TRUE))]]
  
  return(bestRes)
}

################################# Functions for plotting ###############################
plot_discrete_cbar <- function(
  breaks, # Vector of breaks. If +-Inf are used, triangles will be added to the sides of the color bar
  palette = "Greys", # RColorBrewer palette to use
  colors = RColorBrewer::brewer.pal(length(breaks) - 1, palette), # Alternatively, manually set colors
  direction = 1, # Flip colors? Can be 1 or -1
  spacing = "natural", # Spacing between labels. Can be "natural" or "constant"
  border_color = NA, # NA = no border color
  legend_title = NULL,
  legend_direction = "horizontal", # Can be "horizontal" or "vertical"
  font_size = 5,
  expand_size = 1, # Controls spacing around legend plot
  spacing_scaling = 1, # Multiplicative factor for label and legend title spacing
  width = 0.1, # Thickness of color bar
  triangle_size = 0.1 # Relative width of +-Inf triangles
) {
  require(ggplot2)
  if (!(spacing %in% c("natural", "constant"))) stop("spacing must be either 'natural' or 'constant'")
  if (!(direction %in% c(1, -1))) stop("direction must be either 1 or -1")
  if (!(legend_direction %in% c("horizontal", "vertical"))) stop("legend_direction must be either 'horizontal' or 'vertical'")
  breaks = as.numeric(breaks)
  new_breaks = sort(unique(breaks))
  if (any(new_breaks != breaks)) warning("Wrong order or duplicated breaks")
  breaks = new_breaks
  if (class(colors) == "function") colors = colors(length(breaks) - 1)
  if (length(colors) != length(breaks) - 1) stop("Number of colors (", length(colors), ") must be equal to number of breaks (", length(breaks), ") minus 1")
  if (!missing(colors)) warning("Ignoring RColorBrewer palette '", palette, "', since colors were passed manually")
  
  if (direction == -1) colors = rev(colors)
  
  inf_breaks = which(is.infinite(breaks))
  if (length(inf_breaks) != 0) breaks = breaks[-inf_breaks]
  plotcolors = colors
  
  n_breaks = length(breaks)
  
  labels = breaks
  
  if (spacing == "constant") {
    breaks = 1:n_breaks
  }
  
  r_breaks = range(breaks)
  
  cbar_df = data.frame(stringsAsFactors = FALSE,
                       y = breaks,
                       yend = c(breaks[-1], NA),
                       color = as.character(1:n_breaks)
  )[-n_breaks,]
  
  xmin = 1 - width/2
  xmax = 1 + width/2
  
  cbar_plot = ggplot(cbar_df, aes(xmin=xmin, xmax = xmax, ymin = y, ymax = yend, fill = factor(color, levels = 1:length(colors)))) +
    geom_rect(show.legend = FALSE,
              color=border_color)+ 
    
    if (any(inf_breaks == 1)) { # Add < arrow for -Inf
      firstv = breaks[1]
      polystart = data.frame(
        x = c(xmin, xmax, 1),
        y = c(rep(firstv, 2), firstv - diff(r_breaks) * triangle_size)
      )
      plotcolors = plotcolors[-1]
      cbar_plot = cbar_plot +
        geom_polygon(data=polystart, aes(x=x, y=y),
                     show.legend = FALSE,
                     inherit.aes = FALSE,
                     fill = colors[1],
                     color=border_color)
    }
  if (any(inf_breaks > 1)) { # Add > arrow for +Inf
    lastv = breaks[n_breaks]
    polyend = data.frame(
      x = c(xmin, xmax, 1),
      y = c(rep(lastv, 2), lastv + diff(r_breaks) * triangle_size)
    )
    plotcolors = plotcolors[-length(plotcolors)]
    cbar_plot = cbar_plot +
      geom_polygon(data=polyend, aes(x=x, y=y),
                   show.legend = FALSE,
                   inherit.aes = FALSE,
                   fill = colors[length(colors)],
                   color=border_color)
  }
  
  if (legend_direction == "horizontal") { #horizontal legend
    mul = 1
    x = xmin
    xend = xmax
    cbar_plot = cbar_plot + coord_flip()
    angle = 0
    legend_position = xmax + 0.1 * spacing_scaling
  } else { # vertical legend
    mul = -1
    x = xmax
    xend = xmin
    angle = -90
    legend_position = xmax + 0.2 * spacing_scaling
  }
  
  cbar_plot = cbar_plot +
    geom_segment(data=data.frame(y = breaks, yend = breaks),
                 aes(y=y, yend=yend),
                 x = x - 0.05 * mul * spacing_scaling, xend = xend,
                 inherit.aes = FALSE) +
    annotate(geom = 'text', x = x - 0.1 * mul * spacing_scaling, y = breaks,
             label = labels,
             size = font_size, angle=-30) +
    scale_x_continuous(expand = c(expand_size,expand_size)) +
    scale_fill_manual(values=plotcolors) +
    #theme(plot.background=element_rect(fill="darkseagreen"))
    theme_void()
  
  if (!is.null(legend_title)) { # Add legend title
    cbar_plot = cbar_plot +
      annotate(geom = 'text', x = legend_position, y = mean(r_breaks),
               label = legend_title,
               angle = angle,
               size = font_size)
  }
  
  cbar_plot
}

# Make GIF out of diffusion results
saveGif <- function(col_num, image_num, outName, diffusion_mat, timeP, Nx, Ny, lo1, lo2, la1, la2, x0, y0){
  
  # Remove country boarders in the map
  mycrs <- "+proj=longlat +datum=WGS84 +no_defs"
  world <- maps::map("world", fill=TRUE) %$%
    maptools::map2SpatialPolygons(., IDs=names,proj4string=CRS(mycrs))
  while(rgeos::gIsValid(world)==FALSE){
    world <- rgeos::gBuffer(world, byid = TRUE, width = 0, quadsegs = 5, capStyle = "ROUND")
  }
  world <- raster::aggregate(world)
  
  # Map grid
  x <- seq(lo1, lo2, length.out = Ny)
  y <- seq(la1, la2, length.out = Nx)
  diff <- expand.grid(X=x, Y=y)
  
  time_periods <- seq(timeP,0,-1)
  timevect <- timeP+1 - time_periods

  # Colors and time vectors
  colorBreaks = c(1e-5, 3^(seq(-6, log10(1), length.out = col_num)))
  cols <- hcl.colors((length(colorBreaks) - 1L), "YlOrRd", rev=TRUE)
  breakLevels <- levels(cut(c(1), breaks = colorBreaks, right = FALSE)) 
  colorNames <- cols
  names(colorNames) <- breakLevels
  
  # Create a GIF with allele frequencies advancing in time
  saveGIF({
    displayed_periods = round(seq(1, length(time_periods), length.out=image_num))
    
    for (tp in displayed_periods){

      mat_out <- matrix(unlist(diffusion_mat[timevect[tp],]),nrow=Nx,ncol=Ny)
      mat_out[which(topo < 0, arr.ind = TRUE)] <- 0
      diff$Z <- c(t(mat_out[nrow(mat_out):1,]))
      diff$Z <- cut(diff$Z, breaks = colorBreaks, right = FALSE)
      
      mplot <- ggplot() +
        geom_sf() +
        coord_sf(xlim = c(-10, 80), ylim = c(30, 71), expand = FALSE)+
        geom_raster() +
        theme_classic()+
        scale_fill_manual(values = colorNames)+
        geom_tile(diff, mapping = aes(X, Y, fill= Z))+
        geom_polygon(data=world, aes(x=long, y=lat, group=group), fill='NA', color='black')+
        geom_point(aes(x=x0, y=y0), shape=21, size=5, colour="black", fill="darkslategray3")+
        theme(plot.title = element_text(hjust = 0.5, face = "bold"), text=element_text(size=20))+
        ggtitle(paste(round(time_periods[tp])*290, " Years BP"))+
        labs(y="Latitude (°)", x = "Longitude (°)")+
        guides(fill='none')
      
      cbar <- plot_discrete_cbar(round(colorBreaks, 3), colors=cols, spacing="constant", font_size = 4, legend_direction = "horizontal", expand_size = 0.1)
      plot(plot_grid(mplot, cbar, nrow=2, rel_heights = c(10, 1)))
    }},movie.name = outName)
}

######################## Functions for log-likelihood calculations ####################
