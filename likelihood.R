matern_likelihood <-function(theta){
  
  w <- theta[5:6]
  
  t=time
  q=2
  
  # create a spatial autocorrelation signature
  # coordinate list
  
  loc <- coords <- sim_grid_locations
  
  n <- nrow(sim_grid_locations)
  
  if (t==1){
    loc <- loc
  } else {
    for (tt in 1:(t-1)){
      temploc <- matrix(,ncol=2,nrow=nrow(coords))
      for(rr in 1:nrow(coords)){
        temploc[rr,] <- c(coords[rr,1]-tt*w[1],coords[rr,2]-tt*w[2])
      }
      loc <- rbind(loc, temploc)
    }
  }
  
  locations <- loc
  
  theta2 <- function (n,beta0,beta1,beta2,beta3,beta4) {
    theta3 <- beta0 + beta1*(locations[,1] - .5) + beta2*(locations[,2]-.5) + 
      beta3*(locations[,1] - .5)^2 + beta4*(locations[,2] - .5)^2
    theta3 <- matrix(theta3,nrow=nrow(locations),ncol=1)
    return(theta3)
  }
  
  #log.lam1.1<-theta2(n,-3,1,1,-6,-7)
  #log.lam1.2<-theta2(n,-5,1,1,6,-4)
  #logit.phi.1<-theta2(n,0,1,-2,0,1)
  
  #log.lam2.1<-theta2(n,-1.65,0.5,0.5,0,0)
  #log.lam2.2<-theta2(n,-2.8,-1,2,0,-7)
  #logit.phi.2<-theta2(n,-3,-1,2,0,-1)
  
  #log.lam2.1<-theta2(n,-3,-1,-1,-6,-7)
  #log.lam2.2<-theta2(n,-5,-1,-1,6,-4)
  #logit.phi.2<-theta2(n,0,-1,-2,0,1)
  
  log.lam1.1<-theta2(n,-3,1,1,-6,-7)
  log.lam1.2<-theta2(n,-5,1,1,2,-12)
  logit.phi.1<-theta2(n,0,1,-2,0,1)
  
  log.lam2.1<-theta2(n,-3,-1,-1,-6,-7)
  log.lam2.2<-theta2(n,-5,-1,-1,2,-12)
  logit.phi.2<-theta2(n,0,-1,-2,0,1)
  
  KERNEL_LIST <- list()
  
  kernel.local <- array(0, dim = c(2, 2, nrow(locations)))
  for(i in 1:nrow(locations)){
    lam1 <- exp(log.lam1.1[i,])
    lam2 <- exp(log.lam1.2[i,])
    phi <- (pi/2)*exp(logit.phi.1[i,])/(1+exp(logit.phi.1[i,]))
    Pmat <- matrix(c(cos(phi), -sin(phi), sin(phi), cos(phi)), nrow = 2, byrow = T)
    Dmat <- diag(c(lam1, lam2))
    Sigma <- Pmat %*% Dmat %*% t(Pmat)
    kernel.local[, ,i] <-  Sigma
  }
  
  KERNEL_LIST[[1]] <- kernel.local
  
  for(i in 1:nrow(locations)){
    lam1 <- exp(log.lam2.1[i,])
    lam2 <- exp(log.lam2.2[i,])
    phi <- (pi/2)*exp(logit.phi.2[i,])/(1+exp(logit.phi.2[i,]))
    Pmat <- matrix(c(cos(phi), -sin(phi), sin(phi), cos(phi)), nrow = 2, byrow = T)
    Dmat <- diag(c(lam1, lam2))
    Sigma <- Pmat %*% Dmat %*% t(Pmat)
    kernel.local[, ,i] <-  Sigma
  }
  
  KERNEL_LIST[[2]] <- kernel.local
  
  ##Calculate Matern form Nonstationary Covariance function 
  FIN_Sigma.mat <- list()
  dist0 <- list()
  
  for(KK in 1:2){
    Sigma.mat <- matrix(rep(NA, (n*t)^2), nrow = n*t)
    Q.mat <- matrix(rep(NA, (n*t)^2), nrow = n*t)
    Inv_ij <- matrix(rep(NA,4),2,2)
    
    for (i in 1:nrow(locations)) {
      #Sigma.mat[i, i] <- 1
      #Q.mat[i, i] <- 0
      Kernel_i <- KERNEL_LIST[[KK]][, , i]
      det_i <- Kernel_i[1,1] * Kernel_i[2,2] - Kernel_i[1,2] * Kernel_i[2,1]
      for (j in 1:nrow(locations)) {
        Kernel_j <- KERNEL_LIST[[KK]][, , j]
        det_j <- Kernel_j[1,1] * Kernel_j[2,2] - Kernel_j[1,2] * Kernel_j[2,1]
        Kernel_ij <- 0.5 * (Kernel_i + Kernel_j)
        Inv_ij[1,1] <- Kernel_ij[2,2] 
        Inv_ij[2,2] <- Kernel_ij[1,1] 
        Inv_ij[2,1] <- - Kernel_ij[2,1] 
        Inv_ij[1,2] <- - Kernel_ij[1,2] 
        det_ij <- Kernel_ij[1,1] * Kernel_ij[2,2] - Kernel_ij[1,2] * Kernel_ij[2,1]
        x <- c(locations[i,1] - locations[j,1], locations[i,2] - locations[j,2])
        Sigma.mat[i, j] <- sqrt(sqrt(det_i * det_j)/det_ij)
        Q.mat[i, j] <- sqrt(t(x) %*% Inv_ij %*% x/det_ij)
        #Sigma.mat[j, i] <- Sigma.mat[i, j]
        #Q.mat[j, i] <- Q.mat[i, j]
      }
    }
    FIN_Sigma.mat[[KK]] <- Sigma.mat
    dist0[[KK]] <- Q.mat
  }
  
  for (i in 1:nrow(locations)) {
    #Sigma.mat[i, i] <- 1
    #Q.mat[i, i] <- 0
    Kernel_i <- KERNEL_LIST[[1]][, , i]
    det_i <- Kernel_i[1,1] * Kernel_i[2,2] - Kernel_i[1,2] * Kernel_i[2,1]
    for (j in 1:nrow(locations)) {
      Kernel_j <- KERNEL_LIST[[2]][, , j]
      det_j <- Kernel_j[1,1] * Kernel_j[2,2] - Kernel_j[1,2] * Kernel_j[2,1]
      Kernel_ij <- 0.5 * (Kernel_i + Kernel_j)
      Inv_ij[1,1] <- Kernel_ij[2,2] 
      Inv_ij[2,2] <- Kernel_ij[1,1] 
      Inv_ij[2,1] <- - Kernel_ij[2,1] 
      Inv_ij[1,2] <- - Kernel_ij[1,2] 
      det_ij <- Kernel_ij[1,1] * Kernel_ij[2,2] - Kernel_ij[1,2] * Kernel_ij[2,1]
      x <- c(locations[i,1] - locations[j,1], locations[i,2] - locations[j,2])
      Sigma.mat[i, j] <- sqrt(sqrt(det_i * det_j)/det_ij)
      Q.mat[i, j] <- sqrt(t(x) %*% Inv_ij %*% x/det_ij)
      #Sigma.mat[j, i] <- Sigma.mat[i, j]
      #Q.mat[j, i] <- Q.mat[i, j]
    }
  }
  FIN_Sigma.mat[[3]] <- Sigma.mat
  dist0[[3]] <- Q.mat
  
  for (i in 1:nrow(locations)) {
    #Sigma.mat[i, i] <- 1
    #Q.mat[i, i] <- 0
    Kernel_i <- KERNEL_LIST[[2]][, , i]
    det_i <- Kernel_i[1,1] * Kernel_i[2,2] - Kernel_i[1,2] * Kernel_i[2,1]
    for (j in 1:nrow(locations)) {
      Kernel_j <- KERNEL_LIST[[1]][, , j]
      det_j <- Kernel_j[1,1] * Kernel_j[2,2] - Kernel_j[1,2] * Kernel_j[2,1]
      Kernel_ij <- 0.5 * (Kernel_i + Kernel_j)
      Inv_ij[1,1] <- Kernel_ij[2,2] 
      Inv_ij[2,2] <- Kernel_ij[1,1] 
      Inv_ij[2,1] <- - Kernel_ij[2,1] 
      Inv_ij[1,2] <- - Kernel_ij[1,2] 
      det_ij <- Kernel_ij[1,1] * Kernel_ij[2,2] - Kernel_ij[1,2] * Kernel_ij[2,1]
      x <- c(locations[i,1] - locations[j,1], locations[i,2] - locations[j,2])
      Sigma.mat[i, j] <- sqrt(sqrt(det_i * det_j)/det_ij)
      Q.mat[i, j] <- sqrt(t(x) %*% Inv_ij %*% x/det_ij)
      #Sigma.mat[j, i] <- Sigma.mat[i, j]
      #Q.mat[j, i] <- Q.mat[i, j]
    }
  }
  FIN_Sigma.mat[[4]] <- Sigma.mat
  dist0[[4]] <- Q.mat
  
  nu=theta[1:2]
  beta=theta[3]
  rot=theta[4]
  Beta=matrix(0,q,q)
  diag(Beta)=1
  Beta[2,1]=rot
  Beta[1,2]=rot
  
  S=matrix(NA,  q*dim(dist0[[1]])[1], q*dim(dist0[[1]])[1])
  
  for(i in 1:q){
    for(j in 1:q){
      
      temp=(i-1)*dim(dist0[[1]])[1]+1:dim(dist0[[1]])[1]
      temp1=(j-1)*dim(dist0[[1]])[1]+1:dim(dist0[[1]])[1]
      
      if(i==j){
        
        temp2=ifelse(dist0[[i]]!=0,FIN_Sigma.mat[[i]]*(dist0[[i]]/beta)^nu[i] * besselK(dist0[[i]]/beta,nu[i])/(2^(nu[i]-1)*gamma(nu[i])),FIN_Sigma.mat[[i]])
        S[temp,temp1]=temp2
        
      }
      
      if(i!=j & i<j){
        
        nu1=nu[i]
        nu2=nu[j]
        nu3=(nu[i]+nu[j])/2
        
        rho=Beta[i,j]*(gamma(nu1+3/2)/gamma(nu1))^(1/2) * (gamma(nu2+3/2)/gamma(nu2))^(1/2)*gamma(nu3)/(gamma(nu3+3/2))
        #rho <- Beta[i,j]*gamma(nu3)/sqrt(gamma(nu1)*gamma(nu2))
        
        lai=ifelse(dist0[[3]]!=0 ,(dist0[[3]]/beta)^nu3 * besselK(dist0[[3]]/beta,nu3)/(2^(nu3-1)*gamma(nu3))*sqrt(FIN_Sigma.mat[[i]] * FIN_Sigma.mat[[j]])*rho,sqrt(FIN_Sigma.mat[[i]] * FIN_Sigma.mat[[j]])*rho)
        S[temp,temp1]=lai
        #S[temp1,temp]=t(lai)
      }
      
      if(i!=j & i>j){
        
        nu1=nu[i]
        nu2=nu[j]
        nu3=(nu[i]+nu[j])/2
        
        rho=Beta[i,j]*(gamma(nu1+3/2)/gamma(nu1))^(1/2) * (gamma(nu2+3/2)/gamma(nu2))^(1/2)*gamma(nu3)/(gamma(nu3+3/2))
        #rho <- Beta[i,j]*gamma(nu3)/sqrt(gamma(nu1)*gamma(nu2))
        
        lai=ifelse(dist0[[4]]!=0 ,(dist0[[4]]/beta)^nu3 * besselK(dist0[[4]]/beta,nu3)/(2^(nu3-1)*gamma(nu3))*sqrt(FIN_Sigma.mat[[i]] * FIN_Sigma.mat[[j]])*rho,sqrt(FIN_Sigma.mat[[i]] * FIN_Sigma.mat[[j]])*rho)
        S[temp,temp1]=lai
        #S[temp1,temp]=t(lai)
      }
    }
  }
  
  return(S)
}

local_loglik<- function (params,locations, data, p) 
{   
  fixed <- rep(FALSE,6 )
  params <- fixed
  function(s) {
    params[!fixed] <- s
    beta11 <- params[1]
    beta12 <- params[2]
    beta21 <- params[3]
    beta22 <- params[4]
    beta31  <- params[5]
    beta32  <- params[6]
    
    beta10 <- p[1]
    beta20 <- p[2]
    beta30 <- p[3]
    nu <- p[4]
    sigma2<-p[5]
    tau2<-p[6]
    N <- dim(locations)[1]
    m <- dim(data)[2]
    beta1 <- c(beta10,beta11,beta12)
    beta2 <- c(beta20,beta21,beta22)
    beta3 <- c(beta30,beta31,beta32)
    Xmat <- cbind(rep(1,N), locations - colMeans(locations))
    kernel.local <- array(0, dim = c(2, 2, N))
    for ( i in 1:N){
      lam1 <- exp(sum(Xmat[i,] * beta1))
      lam2 <- exp(sum(Xmat[i,] * beta2))
      phi <- (pi/2) * exp(sum(Xmat[i,] * beta3))/(1 + exp(sum(Xmat[i,] * beta3)))
      Pmat <- matrix(c(cos(phi), -sin(phi), sin(phi), cos(phi)), nrow = 2, byrow = T)
      Dmat <- diag(c(lam1, lam2))
      Sigma <- Pmat %*% Dmat %*% t(Pmat)
      kernel.local[, , i] <-  Sigma
    }
    Sigma.mat <- matrix(rep(NA, N^2), nrow = N)
    Q.mat <- matrix(rep(NA, N^2), nrow = N)
    Inv_ij <- matrix(rep(NA,4),2,2)
    for (i in 1:N) {
      Sigma.mat[i, i] <- 1
      Q.mat[i, i] <- 0
      Kernel_i <- kernel.local[, , i]
      det_i <- Kernel_i[1,1] * Kernel_i[2,2] - Kernel_i[1,2] * Kernel_i[2,1]
      if (i < N) {
        for (j in (i + 1):N) {
          Kernel_j <- kernel.local[, , j]
          det_j <- Kernel_j[1,1] * Kernel_j[2,2] - Kernel_j[1,2] * Kernel_j[2,1]
          Kernel_ij <- 0.5 * (Kernel_i + Kernel_j)
          Inv_ij[1,1] <- Kernel_ij[2,2] 
          Inv_ij[2,2] <- Kernel_ij[1,1] 
          Inv_ij[2,1] <- -Kernel_ij[2,1] 
          Inv_ij[1,2] <- -Kernel_ij[1,2] 
          det_ij <- Kernel_ij[1,1] * Kernel_ij[2,2] - Kernel_ij[1,2] * Kernel_ij[2,1]
          x <- locations[i, ] - locations[j, ]
          Sigma.mat[i, j] <- sqrt(sqrt(det_i * det_j)/det_ij)
          Q.mat[i, j] <- sqrt(t(x) %*% Inv_ij %*% x/det_ij)
          Sigma.mat[j, i] <- Sigma.mat[i, j]
          Q.mat[j, i] <- Q.mat[i, j]
        }
      }
    }
    
    cov <- geoR::cov.spatial(Q.mat, cov.model = "matern", 
                             cov.pars = c(sigma2, 1), kappa = nu)
    NS.cov <- Sigma.mat * cov + diag(rep(tau2,N))
    
    Eigen <- eigen(NS.cov)
    Eigen_value <- diag(Eigen$values)
    Eigen_vec <- Eigen$vectors
    Cov.inv <- Eigen_vec %*% diag(1/diag(Eigen_value)) %*% t(Eigen_vec)
    
    loglikelihood <-  m * sum(log(2*pi*diag(Eigen_value))) + 
      sum(diag(t(data) %*% Cov.inv %*% data))
    if (abs(loglikelihood) == Inf) loglikelihood <- 1e+06
    return(loglikelihood)
  }
}