library(parallel)
library(MASS)
library(StatMatch)

#no_cores <- detectCores() - 6
#cl <- makeCluster(no_cores)

WORKSTATION = 3

if(WORKSTATION == 1){
  root <- '/run/user/156071/gvfs/sftp:host=ilogin.dragon.kaust.edu.sa,user=salvanmo/scratch/dragon/intel/salvanmo/'
}else if(WORKSTATION == 2){
  root <- '/scratch/dragon/intel/salvanmo/'
}else{
  root <- '/Users/salvanmo/Documents/GitHub/scaling-robot/'
}

source(file=paste(root,"NSconvo_fit_multi.R",sep=''))
source(file=paste(root,"matern_cov.R",sep=''))

n<-400
N<-sqrt(n)
t=1
grid_x = seq(from=0.5, to=1, length.out = N)
grid_y = seq(from=0.5, to=1, length.out = N)
sim_grid_locations = expand.grid(grid_x, grid_y)
sim_grid_locations <- matrix(cbind(sim_grid_locations[,1],sim_grid_locations[,2]),ncol=2)

##generate locations
matern_theta_strong=c(0.5,2,4,0.8)
thets <- matern_theta_strong
materncov <- matern_cov_regular_grid_v4(thets,c(.05001,.05001),time=t)

set.seed(12345678) #12, 1234, 112

conso_matern_params <- list()
conso_nonstationary_params <- list()

for(simulation in 1:100){
  A <- mvrnorm(1000, rep(0,dim(materncov)[1]),materncov)
  
  N.mc <- 4
  fit.radius <- 0.25 #0.15 works
  
  locations <- coords <- as.matrix(sim_grid_locations)
  N <- dim(coords)[1]
  data <- matrix(A, nrow=2*N)
  
  lon_min <- min(coords[,1])
  lon_max <- max(coords[,1])
  lat_min <- min(coords[,2])
  lat_max <- max(coords[,2])
  
  #=======================================
  # mixture component knot locations
  #=======================================
  mc_x <- seq(from = lon_min + 0.5*(lon_max - lon_min)/floor(sqrt(N.mc)),
              to = lon_max - 0.5*(lon_max - lon_min)/floor(sqrt(N.mc)),
              length = floor(sqrt(N.mc)) )
  mc_y <- seq(from = lat_min + 0.5*(lat_max - lat_min)/floor(sqrt(N.mc)),
              to = lat_max - 0.5*(lat_max - lat_min)/floor(sqrt(N.mc)),
              length = floor(sqrt(N.mc)) )
  mc.locations <- expand.grid( mc_x, mc_y )
  mc.locations <- matrix(c(mc.locations[,1], mc.locations[,2]), ncol=2, byrow=F)
  
  K <- dim(mc.locations)[1]
  
  #===========================================================================
  # Check the mixture component locations
  #===========================================================================
  check.mc.locs <- mc_N( coords, mc.locations, fit.radius )
  lambda.w <- ( 0.5*min(dist(mc.locations)) )^2
  
  mean.model = data ~ 1
  data <- A[1:400]
  OLS.model1 <- lm(mean.model, x=TRUE)
  data <- A[401:800]
  OLS.model2 <- lm(mean.model, x=TRUE)
  Xmat <- matrix( unname( OLS.model1$x ), nrow=2*N )
  
  #===========================================================================
  # Specify lower, upper, and initial parameter values for optim()
  #===========================================================================
  
  max.distance <- sqrt(sum((c(lon_min,lat_min) - c(lon_max,lat_max))^2))
  
  resid.var <- (max(summary(OLS.model1)$sigma,summary(OLS.model2)$sigma))^2
  
  lam1.LB <- 1e-05
  lam2.LB <- 1e-05
  tausq.local.LB <- 1e-05
  sigmasq.local.LB <- 1e-05
  kappa.local.LB <- 1e-05
  
  lam1.UB <- max.distance/4
  lam2.UB <- max.distance/4
  tausq.local.UB <- 4*resid.var
  sigmasq.local.UB <- 4*resid.var
  kappa.local.UB <- 30
  
  lam1.init <- max.distance/10
  lam2.init <- max.distance/10
  tausq.local.init <- 0.1*resid.var
  sigmasq.local.init <- 0.9*resid.var
  kappa.local.init <- 1
  
  # Storage for the mixture component kernels
  #mc.kernels <- array(NA, dim=c(2, 2, K))
  MLEs.save <- matrix(NA, K, 14)
  coords_mc_dist <- StatMatch::mahalanobis.dist(data.x = coords, data.y = mc.locations, vc = diag(2))
  
  #data <- matrix(A, nrow=2*N)
  
  for( k in 1:K ){
    
    ind_local <- (coords_mc_dist[,k] <= fit.radius)
    
    # Subset
    temp.locations <- coords[ind_local,]
    n.fit <- dim(temp.locations)[1]
    
    if(!is.matrix(A)){
      temp.data <- A[rep(ind_local,2*t)]
    }else{
      temp.data <- A[,rep(ind_local,2*t)]
    }
    
    Xtemp <- as.matrix(rep(1,length(which(ind_local==TRUE))),ncol=1)
    
    #####################################################
    # Local estimation
    f_loglik <- make_local_lik(locations = temp.locations, Xmat = Xtemp, data = temp.data, fixed = rep(FALSE, 14), time=t)
    
    MLEs <- optim( par = c(lam1.init, lam2.init, pi/4,lam1.init, lam2.init, pi/4,thets,1,1,1,1),
                   fn = f_loglik, control=list(maxit=20000,parscale=c(lam1.init, lam2.init, pi/4,lam1.init, lam2.init, pi/4, thets,1,1,1,1),trace=5))
    
    MLEs.save[k,] <- MLEs$par
    
    # Save the kernel matrix
    #mc.kernels[,,k] <- kernel_cov(MLEs$par[1:3])

  }
  
  beta0_lam1.1 <- MLEs.save[,1]
  beta0_lam1.2 <- MLEs.save[,2]
  beta0_phi.1 <- MLEs.save[,3]
  beta0_lam2.1 <- MLEs.save[,4]
  beta0_lam2.2 <- MLEs.save[,5]
  beta0_phi.2 <- MLEs.save[,6]
  
  #h <- rslt$lambda.w
  p1 <- matrix(0,N.mc,14)
  for (i in 1:N.mc) {
    p1[i,] <- c(log(beta0_lam1.1[i]),log(beta0_lam1.2[i]),log(beta0_phi.1[i]/(pi/2 - beta0_phi.1[i])),
                log(beta0_lam2.1[i]),log(beta0_lam2.2[i]),log(beta0_phi.2[i]/(pi/2 - beta0_phi.1[i])),
                colMeans(MLEs.save[,7:12]),MLEs.save[i,13:14])
  }
  
  nk <- N.mc
  Nk <- sqrt(nk)
  nkr <- n/nk
  NKR <- sqrt(nkr)
  local.index <- rep(0, nkr)
  for(j in 1:NKR){
    local.index[(1 + (j - 1) * NKR):(j * NKR)] <- (1 + (j - 1) * Nk * NKR):(((j - 1) * Nk * NKR) + NKR)
  }
  LL.est <- matrix(0,nk,12)
  for(l in 1 : Nk){
    for (k in 1 : Nk) {
      cat(l,k,'\n')
      index <- (k-1)*NKR+ local.index + Nk*nkr*(l-1)
      locations <- coords
      temp.locations<- locations[index,]
      
      if(!is.matrix(A)){
        temp.data <- A[rep(index,2*t)]
        temp.data <- A[c(index,index + 400)]
      }else{
        temp.data <- A[,rep(index,2*t)]
        temp.data <- A[,c(index,index + 400)]
      }
      
      #temp.data <- as.matrix(c(A[1:(n*t)][c(index,index+n,index+2*n,index+3*n,index+4*n)],A[(n*t+1):(n*t*2)][c(index,index+n,index+2*n,index+3*n,index+4*n)]), nrow=2*length(index))
      
      make.local.loglik <- local_loglik(locations = temp.locations, 
                                        data = temp.data, p = p1[l+(k-1)*Nk,], time=t) 
      #test_theta <- c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1)
      #test_theta <- rep(0.001,12)
      #test_theta <- rep(c(239.4240956,0.4565653,222.6135517,0.4972399,-4.6949158,2.3519542),2)
      #test_theta <- rep(c(-20.4240956,0.4565653,12.6135517,0.4972399,4.6949158,2.3519542),2)
      #test_theta <- rep(0,12)
      
      test_theta <- c(-1.50380504,  -0.07558557,   4.82244196,  -5.66094136,   2.99578278,  43.99618927,
                      19.24018098, -14.63719610,  -0.51861398,  15.23767899,   2.68615214,  -3.00711322)
      MLEs.local <- optim(test_theta, make.local.loglik,control=list(maxit=20000, parscale=test_theta, trace=5))
      #test_theta <- MLEs.local$par
      #MLEs.local <- optim(test_theta, make.local.loglik,control=list(maxit=10000, parscale=test_theta, trace=5))
      LL.est[l+(k-1)*Nk,]<-MLEs.local$par
    }
  }
  conso_matern_params[[simulation]] <- p1
  conso_nonstationary_params[[simulation]] <- LL.est
}
save(conso_matern_params,conso_nonstationary_params, file="/scratch/dragon/intel/salvanmo/Results/nonstationary2.Rdata")

#  -1.50380504  -0.07558557   4.82244196  -5.66094136   2.99578278  43.99618927
#  19.24018098 -14.63719610  -0.51861398  15.23767899   2.68615214  -3.00711322

#  -1.505190   1.099794   8.066450 -10.759477  -1.664848  56.165726   4.911672
#  -5.670777   3.312877   5.553745   2.910514  -7.774712