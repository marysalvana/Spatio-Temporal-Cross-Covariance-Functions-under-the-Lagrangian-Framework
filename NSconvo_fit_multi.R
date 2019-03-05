
mc_N <- function( coords, mc.locations, fit.radius ){
  
  
  K <- dim(mc.locations)[1]
  mc.N.fit <- rep(NA,K)
  
  for( k in 1:K ){
    
    temp.locs <- coords[ abs(coords[,1]-mc.locations[k,1]) <= fit.radius
                         & (abs(coords[,2] - mc.locations[k,2]) <= fit.radius), ]
    
    # Isolate the data/locations to be used for calculating the local kernel
    distances <- rep(NA,dim(temp.locs)[1])
    
    for(i in 1:dim(temp.locs)[1]){
      distances[i] <- sqrt(sum((temp.locs[i,] - mc.locations[k,])^2))
    }
    
    temp.locations <- temp.locs[distances <= fit.radius,]
    n.fit <- dim(temp.locations)[1]
    
    mc.N.fit[k] <- n.fit
    
  }
  
  return(mc.N.fit)
}

make_local_lik <- function( locations, data, Xmat,params,fixed = rep(FALSE, 14), time){
  params <- fixed
  function(p) {
    params[!fixed] <- p
    
    lam1.1 <- params[1]
    lam1.2 <- params[2]
    eta.1 <- params[3]
    
    lam2.1 <- params[4]
    lam2.2 <- params[5]
    eta.2 <- params[6]
    
    b0 <- params[13:14]
    XX=rbind(b0[1]*Xmat,b0[2]*Xmat)
    # Setup and covariance calculation
    N <- dim(locations)[1]
    q=2
    t=time
    
    Pmat1 <- matrix(c(cos(eta.1), -sin(eta.1), sin(eta.1), cos(eta.1)), nrow = 2, byrow = T)
    Dmat1 <- diag(c(lam1.1, lam1.2))
    Sigma1 <- Pmat1 %*% Dmat1 %*% t(Pmat1)
    distances1 <- StatMatch::mahalanobis.dist(data.x = locations, vc = Sigma1)
    
    Pmat2 <- matrix(c(cos(eta.2), -sin(eta.2), sin(eta.2), cos(eta.2)), nrow = 2, byrow = T)
    Dmat2 <- diag(c(lam2.1, lam2.2))
    Sigma2 <- Pmat2 %*% Dmat2 %*% t(Pmat2)
    distances2 <- StatMatch::mahalanobis.dist(data.x = locations, vc = Sigma2)
    
    distances3 <- StatMatch::mahalanobis.dist(data.x = locations, vc = 0.5*(Sigma2+Sigma2))
    
    NS.cov <-  matern_cov_regular_grid_v2_for_estimation_sim_step1(params[7:12],distances1,distances2,distances3)
    
    loglikelihood <-0
    
    possibleError <- tryCatch(
      cov.chol <- chol(NS.cov),
      error=function(e) e
    )
    # Check for error before calculating the likelihood
    if(inherits(possibleError, "error")){
      loglikelihood <- Inf
    } else if(params[1] < lam1.LB | params[1] > lam1.UB | params[2] < lam2.LB | params[2] > lam2.UB | params[3] < 0 | params[3] > pi/2 |
              params[4] < lam1.LB | params[4] > lam1.UB | params[5] < lam2.LB | params[5] > lam2.UB | params[6] < 0 | params[6] > pi/2 | 
              params[7] < 0.001 | params[7] > 10 | params[8] < 0.001 | params[8] > 10 | params[9] < 0.001 | params[9] > 10 |
              params[10] < -1 | params[10] > 1 | params[11] < 0.001 | params[11] > 10 | params[12] < 0.001 | params[12] > 10 |
              params[13] < -10 | params[13] > 10 | params[14] < -10 | params[14] > 10 ){
      
      loglikelihood <- Inf
    }else{
      
      if(!is.matrix(data)){
        tmp1 <- backsolve(cov.chol,  matrix(c(data[1:(length(data)/(2*t))],data[length(data)/2+1:(length(data)/(2*t))]),ncol=1) - XX, transpose = TRUE)
        #tmp1 <- backsolve(cov.chol, c(data)-c(rep(b0[1],dim(data)[1]/2),rep(b0[2],dim(data)[1]/2)), transpose = TRUE)
        
        ResCinvRes <- t(tmp1) %*% tmp1
        loglikelihood <- loglikelihood + q * sum(log(diag(cov.chol))) + 0.5 * sum(diag(ResCinvRes))
      }else{
        for(tt in 1:nrow(data)){
          tmp1 <- backsolve(cov.chol,  matrix(c(data[tt,1:(dim(data)[2]/(2*t))],data[tt,dim(data)[2]/2+1:(dim(data)[2]/(2*t))]),ncol=1) - XX, transpose = TRUE)
          #tmp1 <- backsolve(cov.chol, c(data)-c(rep(b0[1],dim(data)[1]/2),rep(b0[2],dim(data)[1]/2)), transpose = TRUE)
          
          ResCinvRes <- t(tmp1) %*% tmp1
          loglikelihood <- loglikelihood + q * sum(log(diag(cov.chol))) + 0.5 * sum(diag(ResCinvRes))
          
        }
      }
    }
    if (abs(loglikelihood) == Inf) { loglikelihood <- 1e+06 } # Make sure not Inf
    return(loglikelihood)
  }
}

make_local_lik_v2 <- function( locations, data, Xmat,params,fixed = rep(FALSE, 14), time){
  params <- fixed
  function(p) {
    params[!fixed] <- p
    
    lam1.1 <- params[1]
    lam1.2 <- params[2]
    eta.1 <- params[3]
    
    lam2.1 <- params[4]
    lam2.2 <- params[5]
    eta.2 <- params[6]
    
    b0 <- params[13:14]
    XX=rbind(b0[1]*Xmat,b0[2]*Xmat)
    # Setup and covariance calculation
    N <- dim(locations)[1]
    q=2
    t=time
    
    Pmat1 <- matrix(c(cos(eta.1), -sin(eta.1), sin(eta.1), cos(eta.1)), nrow = 2, byrow = T)
    Dmat1 <- diag(c(lam1.1, lam1.2))
    Sigma1 <- Pmat1 %*% Dmat1 %*% t(Pmat1)
    distances1 <- StatMatch::mahalanobis.dist(data.x = locations, vc = Sigma1)
    
    Pmat2 <- matrix(c(cos(eta.2), -sin(eta.2), sin(eta.2), cos(eta.2)), nrow = 2, byrow = T)
    Dmat2 <- diag(c(lam2.1, lam2.2))
    Sigma2 <- Pmat2 %*% Dmat2 %*% t(Pmat2)
    distances2 <- StatMatch::mahalanobis.dist(data.x = locations, vc = Sigma2)
    
    distances3 <- StatMatch::mahalanobis.dist(data.x = locations, vc = 0.5*(Sigma2+Sigma2))
    
    NS.cov <-  matern_cov_regular_grid_v2_for_estimation_sim_step1(params[7:12],distances1,distances2,distances3)
    
    loglikelihood <-0
    
    possibleError <- tryCatch(
      cov.chol <- chol(NS.cov),
      error=function(e) e
    )
    # Check for error before calculating the likelihood
    if(inherits(possibleError, "error")){
      loglikelihood <- 1e+020
    } else{
      
      for(tt in 1:1){
        tmp1 <- backsolve(cov.chol,  matrix(c(data[tt,1:(dim(data)[2]/(2*t))],data[tt,dim(data)[2]/2+1:(dim(data)[2]/(2*t))]),ncol=1) - XX, transpose = TRUE)
        #tmp1 <- backsolve(cov.chol, c(data)-c(rep(b0[1],dim(data)[1]/2),rep(b0[2],dim(data)[1]/2)), transpose = TRUE)
        
        ResCinvRes <- t(tmp1) %*% tmp1
        loglikelihood <- loglikelihood + q * sum(log(diag(cov.chol))) + 0.5 * sum(diag(ResCinvRes))
      }
      # Likelihood calculation
      
    }
    if (abs(loglikelihood) == Inf) { loglikelihood <- 1e+020 } # Make sure not Inf
    return(loglikelihood)
  }
}

local_loglik<- function (params,locations, data, p, time) 
{   
  fixed <- rep(FALSE,12)
  params <- fixed
  function(s) {
    params[!fixed] <- s
    beta1.11 <- params[1]
    beta1.12 <- params[2]
    beta1.21 <- params[3]
    beta1.22 <- params[4]
    beta1.31  <- params[5]
    beta1.32  <- params[6]
    
    beta2.11 <- params[7]
    beta2.12 <- params[8]
    beta2.21 <- params[9]
    beta2.22 <- params[10]
    beta2.31  <- params[11]
    beta2.32  <- params[12]
    
    beta1.10 <- p[1]
    beta1.20 <- p[2]
    beta1.30 <- p[3]
    beta2.10 <- p[4]
    beta2.20 <- p[5]
    beta2.30 <- p[6]
    
    
    NN <- dim(locations)[1]
    beta1.1 <- c(beta1.10,beta1.11,beta1.12)
    beta1.2 <- c(beta1.20,beta1.21,beta1.22)
    beta1.3 <- c(beta1.30,beta1.31,beta1.32)
    
    beta2.1 <- c(beta2.10,beta2.11,beta2.12)
    beta2.2 <- c(beta2.20,beta2.21,beta2.22)
    beta2.3 <- c(beta2.30,beta2.31,beta2.32)
    
    t <- time
    
    Xmat <- cbind(rep(1,NN), locations - colMeans(locations))
    
    KERNEL_LIST <- list()
    
    kernel.local <- array(0, dim = c(2, 2, nrow(locations)))
    for(i in 1:nrow(locations)){
      lam1 <- exp(sum(Xmat[i,] * beta1.1))
      lam2 <- exp(sum(Xmat[i,] * beta1.2))
      phi <- (pi/2) * exp(sum(Xmat[i,] * beta1.3))/(1 + exp(sum(Xmat[i,] * beta1.3)))
      Pmat <- matrix(c(cos(phi), -sin(phi), sin(phi), cos(phi)), nrow = 2, byrow = T)
      Dmat <- diag(c(lam1, lam2))
      Sigma <- Pmat %*% Dmat %*% t(Pmat)
      kernel.local[, ,i] <-  Sigma
    }
    
    KERNEL_LIST[[1]] <- kernel.local
    
    for(i in 1:nrow(locations)){
      lam1 <- exp(sum(Xmat[i,] * beta2.1))
      lam2 <- exp(sum(Xmat[i,] * beta2.2))
      phi <- (pi/2) * exp(sum(Xmat[i,] * beta2.3))/(1 + exp(sum(Xmat[i,] * beta2.3)))
      Pmat <- matrix(c(cos(phi), -sin(phi), sin(phi), cos(phi)), nrow = 2, byrow = T)
      Dmat <- diag(c(lam1, lam2))
      Sigma <- Pmat %*% Dmat %*% t(Pmat)
      kernel.local[, ,i] <-  Sigma
    }
    
    KERNEL_LIST[[2]] <- kernel.local
    
    FIN_Sigma.mat <- list()
    dist0 <- list()
    
    for(KK in 1:2){
      Sigma.mat <- matrix(rep(NA, NN^2), nrow = NN)
      Q.mat <- matrix(rep(NA, NN^2), nrow = NN)
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
    
    q=2
    
    nu = c(0.5,2)
    beta = 4
    #nu=p[7:8]
    #beta=p[9]
    #rot=p[10]
    Beta=matrix(0,q,q)
    diag(Beta)=1
    Beta[2,1]=rot
    Beta[1,2]=rot
    #var =p[11:12]
    
    b0 <- p[13:14]
    XX=matrix(c(b0[1]*Xmat[,1],b0[2]*Xmat[,1]),ncol=1)
    
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
    
    possibleError <- tryCatch(
      cov.chol <- chol(S),
      error=function(e) e
    )
    loglikelihood <-0
    
    # Check for error before calculating the likelihood
    if(inherits(possibleError, "error")){
      loglikelihood <- 1e+06
    } else{
      
      if(!is.matrix(data)){
        tmp1 <- backsolve(cov.chol,  matrix(c(data[1:(length(data)/(2*t))],data[length(data)/2+1:(length(data)/(2*t))]),ncol=1)-XX, transpose = TRUE)
        #tmp1 <- backsolve(cov.chol, c(data)-c(rep(b0[1],dim(data)[1]/2),rep(b0[2],dim(data)[1]/2)), transpose = TRUE)
        
        ResCinvRes <- t(tmp1) %*% tmp1
        loglikelihood <- loglikelihood + q * sum(log(diag(cov.chol))) + 0.5 * sum(diag(ResCinvRes))
      }else{
        for(tt in 1:nrow(data)){
          tmp1 <- backsolve(cov.chol,  matrix(c(data[tt,1:(dim(data)[2]/(2*t))],data[tt,dim(data)[2]/2+1:(dim(data)[2]/(2*t))]),ncol=1)-XX, transpose = TRUE)
          #tmp1 <- backsolve(cov.chol, c(data)-c(rep(b0[1],dim(data)[1]/2),rep(b0[2],dim(data)[1]/2)), transpose = TRUE)
          
          ResCinvRes <- t(tmp1) %*% tmp1
          loglikelihood <- loglikelihood + q * sum(log(diag(cov.chol))) + 0.5 * sum(diag(ResCinvRes))
          
        }
      }
    }
    if (abs(loglikelihood) == Inf) { loglikelihood <- 1e+06 } # Make sure not Inf
    return(loglikelihood)
  }
}

local_loglik_v3<- function (params,locations, data, p, time) 
{   
  fixed <- rep(FALSE,18)
  params <- fixed
  function(s) {
    params[!fixed] <- s
    beta1.11 <- params[1]
    beta1.12 <- params[2]
    beta1.21 <- params[3]
    beta1.22 <- params[4]
    beta1.31  <- params[5]
    beta1.32  <- params[6]
    
    beta2.11 <- params[7]
    beta2.12 <- params[8]
    beta2.21 <- params[9]
    beta2.22 <- params[10]
    beta2.31  <- params[11]
    beta2.32  <- params[12]
    
    beta1.10 <- params[13]
    beta1.20 <- params[14]
    beta1.30 <- params[15]
    beta2.10 <- params[16]
    beta2.20 <- params[17]
    beta2.30 <- params[18]
    
    #beta1.10 <- p[1]
    #beta1.20 <- p[2]
    #beta1.30 <- p[3]
    #beta2.10 <- p[4]
    #beta2.20 <- p[5]
    #beta2.30 <- p[6]
    
    
    NN <- dim(locations)[1]
    beta1.1 <- c(beta1.10,beta1.11,beta1.12)
    beta1.2 <- c(beta1.20,beta1.21,beta1.22)
    beta1.3 <- c(beta1.30,beta1.31,beta1.32)
    
    beta2.1 <- c(beta2.10,beta2.11,beta2.12)
    beta2.2 <- c(beta2.20,beta2.21,beta2.22)
    beta2.3 <- c(beta2.30,beta2.31,beta2.32)
    
    t <- time
    
    Xmat <- cbind(rep(1,NN), locations - colMeans(locations))
    
    KERNEL_LIST <- list()
    
    kernel.local <- array(0, dim = c(2, 2, nrow(locations)))
    for(i in 1:nrow(locations)){
      lam1 <- exp(sum(Xmat[i,] * beta1.1))
      lam2 <- exp(sum(Xmat[i,] * beta1.2))
      phi <- (pi/2) * exp(sum(Xmat[i,] * beta1.3))/(1 + exp(sum(Xmat[i,] * beta1.3)))
      Pmat <- matrix(c(cos(phi), -sin(phi), sin(phi), cos(phi)), nrow = 2, byrow = T)
      Dmat <- diag(c(lam1, lam2))
      Sigma <- Pmat %*% Dmat %*% t(Pmat)
      kernel.local[, ,i] <-  Sigma
    }
    
    KERNEL_LIST[[1]] <- kernel.local
    
    for(i in 1:nrow(locations)){
      lam1 <- exp(sum(Xmat[i,] * beta2.1))
      lam2 <- exp(sum(Xmat[i,] * beta2.2))
      phi <- (pi/2) * exp(sum(Xmat[i,] * beta2.3))/(1 + exp(sum(Xmat[i,] * beta2.3)))
      Pmat <- matrix(c(cos(phi), -sin(phi), sin(phi), cos(phi)), nrow = 2, byrow = T)
      Dmat <- diag(c(lam1, lam2))
      Sigma <- Pmat %*% Dmat %*% t(Pmat)
      kernel.local[, ,i] <-  Sigma
    }
    
    KERNEL_LIST[[2]] <- kernel.local
    
    FIN_Sigma.mat <- list()
    dist0 <- list()
    
    for(KK in 1:2){
      Sigma.mat <- matrix(rep(NA, NN^2), nrow = NN)
      Q.mat <- matrix(rep(NA, NN^2), nrow = NN)
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
    
    q=2
    
    nu = c(0.5,2)
    beta = 4
    #nu=p[7:8]
    #beta=p[9]
    #rot=p[10]
    Beta=matrix(0,q,q)
    diag(Beta)=1
    Beta[2,1]=rot
    Beta[1,2]=rot
    #var =p[11:12]
    
    b0 <- p[13:14]
    XX=matrix(c(b0[1]*Xmat[,1],b0[2]*Xmat[,1]),ncol=1)
    
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
    
    possibleError <- tryCatch(
      cov.chol <- chol(S),
      error=function(e) e
    )
    loglikelihood <-0
    
    # Check for error before calculating the likelihood
    if(inherits(possibleError, "error")){
      loglikelihood <- 1e+06
    } else{
      loglikelihood <- sum((S-new_cov_target)^2)
    }
    if (abs(loglikelihood) == Inf) { loglikelihood <- 1e+06 } # Make sure not Inf
    return(loglikelihood)
  }
}

local_loglik_v4<- function (params,locations, data, p, time) 
{   
  fixed <- rep(FALSE,18)
  params <- fixed
  function(s) {
    params[!fixed] <- s
    beta1.11 <- params[1]
    beta1.12 <- params[2]
    beta1.21 <- params[3]
    beta1.22 <- params[4]
    beta1.31  <- params[5]
    beta1.32  <- params[6]
    
    beta2.11 <- params[7]
    beta2.12 <- params[8]
    beta2.21 <- params[9]
    beta2.22 <- params[10]
    beta2.31  <- params[11]
    beta2.32  <- params[12]
    
    beta1.10 <- params[13]
    beta1.20 <- params[14]
    beta1.30 <- params[15]
    beta2.10 <- params[16]
    beta2.20 <- params[17]
    beta2.30 <- params[18]
    
    #beta1.10 <- p[1]
    #beta1.20 <- p[2]
    #beta1.30 <- p[3]
    #beta2.10 <- p[4]
    #beta2.20 <- p[5]
    #beta2.30 <- p[6]
    
    
    NN <- dim(locations)[1]
    beta1.1 <- c(beta1.10,beta1.11,beta1.12)
    beta1.2 <- c(beta1.20,beta1.21,beta1.22)
    beta1.3 <- c(beta1.30,beta1.31,beta1.32)
    
    beta2.1 <- c(beta2.10,beta2.11,beta2.12)
    beta2.2 <- c(beta2.20,beta2.21,beta2.22)
    beta2.3 <- c(beta2.30,beta2.31,beta2.32)
    
    t <- time
    
    Xmat <- cbind(rep(1,NN), locations - colMeans(locations))
    
    KERNEL_LIST <- list()
    
    kernel.local <- array(0, dim = c(2, 2, nrow(locations)))
    for(i in 1:nrow(locations)){
      lam1 <- exp(sum(Xmat[i,] * beta1.1))
      lam2 <- exp(sum(Xmat[i,] * beta1.2))
      phi <- (pi/2) * exp(sum(Xmat[i,] * beta1.3))/(1 + exp(sum(Xmat[i,] * beta1.3)))
      Pmat <- matrix(c(cos(phi), -sin(phi), sin(phi), cos(phi)), nrow = 2, byrow = T)
      Dmat <- diag(c(lam1, lam2))
      Sigma <- Pmat %*% Dmat %*% t(Pmat)
      kernel.local[, ,i] <-  Sigma
    }
    
    KERNEL_LIST[[1]] <- kernel.local
    
    for(i in 1:nrow(locations)){
      lam1 <- exp(sum(Xmat[i,] * beta2.1))
      lam2 <- exp(sum(Xmat[i,] * beta2.2))
      phi <- (pi/2) * exp(sum(Xmat[i,] * beta2.3))/(1 + exp(sum(Xmat[i,] * beta2.3)))
      Pmat <- matrix(c(cos(phi), -sin(phi), sin(phi), cos(phi)), nrow = 2, byrow = T)
      Dmat <- diag(c(lam1, lam2))
      Sigma <- Pmat %*% Dmat %*% t(Pmat)
      kernel.local[, ,i] <-  Sigma
    }
    
    KERNEL_LIST[[2]] <- kernel.local
    
    FIN_Sigma.mat <- list()
    dist0 <- list()
    
    for(KK in 1:2){
      Sigma.mat <- matrix(rep(NA, NN^2), nrow = NN)
      Q.mat <- matrix(rep(NA, NN^2), nrow = NN)
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
    
    q=2
    
    nu = c(0.5,2)
    beta = 4
    #nu=p[7:8]
    #beta=p[9]
    #rot=p[10]
    Beta=matrix(0,q,q)
    diag(Beta)=1
    Beta[2,1]=rot
    Beta[1,2]=rot
    #var =p[11:12]
    
    b0 <- p[13:14]
    XX=matrix(c(b0[1]*Xmat[,1],b0[2]*Xmat[,1]),ncol=1)
    
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
    
    possibleError <- tryCatch(
      cov.chol <- chol(S),
      error=function(e) e
    )
    loglikelihood <-0
    
    # Check for error before calculating the likelihood
    if(inherits(possibleError, "error")){
      loglikelihood <- 1e+06
    } else{
      
      if(!is.matrix(data)){
        tmp1 <- backsolve(cov.chol,  matrix(c(data[1:(length(data)/(2*t))],data[length(data)/2+1:(length(data)/(2*t))]),ncol=1)-XX, transpose = TRUE)
        #tmp1 <- backsolve(cov.chol, c(data)-c(rep(b0[1],dim(data)[1]/2),rep(b0[2],dim(data)[1]/2)), transpose = TRUE)
        
        ResCinvRes <- t(tmp1) %*% tmp1
        loglikelihood <- loglikelihood + q * sum(log(diag(cov.chol))) + 0.5 * sum(diag(ResCinvRes))
      }else{
        for(tt in 1:nrow(data)){
          tmp1 <- backsolve(cov.chol,  matrix(c(data[tt,1:(dim(data)[2]/(2*t))],data[tt,dim(data)[2]/2+1:(dim(data)[2]/(2*t))]),ncol=1)-XX, transpose = TRUE)
          #tmp1 <- backsolve(cov.chol, c(data)-c(rep(b0[1],dim(data)[1]/2),rep(b0[2],dim(data)[1]/2)), transpose = TRUE)
          
          ResCinvRes <- t(tmp1) %*% tmp1
          loglikelihood <- loglikelihood + q * sum(log(diag(cov.chol))) + 0.5 * sum(diag(ResCinvRes))
          
        }
      }
    }
    if (abs(loglikelihood) == Inf) { loglikelihood <- 1e+06 } # Make sure not Inf
    return(loglikelihood)
  }
}


local_loglik_v2<- function (params,locations, data, p) 
{   
  fixed <- rep(FALSE,12)
  params <- fixed
  function(s) {
    params[!fixed] <- s
    beta1.11 <- params[1]
    beta1.12 <- params[2]
    beta1.21 <- params[3]
    beta1.22 <- params[4]
    beta1.31  <- params[5]
    beta1.32  <- params[6]
    
    beta2.11 <- params[7]
    beta2.12 <- params[8]
    beta2.21 <- params[9]
    beta2.22 <- params[10]
    beta2.31  <- params[11]
    beta2.32  <- params[12]
    
    beta1.10 <- p[1]
    beta1.20 <- p[2]
    beta1.30 <- p[3]
    beta2.10 <- p[4]
    beta2.20 <- p[5]
    beta2.30 <- p[6]
    
    N <- dim(locations)[1]
    beta1.1 <- c(beta1.10,beta1.11,beta1.12)
    beta1.2 <- c(beta1.20,beta1.21,beta1.22)
    beta1.3 <- c(beta1.30,beta1.31,beta1.32)
    
    beta2.1 <- c(beta2.10,beta2.11,beta2.12)
    beta2.2 <- c(beta2.20,beta2.21,beta2.22)
    beta2.3 <- c(beta2.30,beta2.31,beta2.32)
    
    Xmat <- cbind(rep(1,N), locations - colMeans(locations))
    
    KERNEL_LIST <- list()
    
    kernel.local <- array(0, dim = c(2, 2, nrow(locations)))
    for(i in 1:nrow(locations)){
      lam1 <- exp(sum(Xmat[i,] * beta1.1))
      lam2 <- exp(sum(Xmat[i,] * beta1.2))
      phi <- (pi/2) * exp(sum(Xmat[i,] * beta1.3))/(1 + exp(sum(Xmat[i,] * beta1.3)))
      Pmat <- matrix(c(cos(phi), -sin(phi), sin(phi), cos(phi)), nrow = 2, byrow = T)
      Dmat <- diag(c(lam1, lam2))
      Sigma <- Pmat %*% Dmat %*% t(Pmat)
      kernel.local[, ,i] <-  Sigma
    }
    
    KERNEL_LIST[[1]] <- kernel.local
    
    for(i in 1:nrow(locations)){
      lam1 <- exp(sum(Xmat[i,] * beta2.1))
      lam2 <- exp(sum(Xmat[i,] * beta2.2))
      phi <- (pi/2) * exp(sum(Xmat[i,] * beta2.3))/(1 + exp(sum(Xmat[i,] * beta2.3)))
      Pmat <- matrix(c(cos(phi), -sin(phi), sin(phi), cos(phi)), nrow = 2, byrow = T)
      Dmat <- diag(c(lam1, lam2))
      Sigma <- Pmat %*% Dmat %*% t(Pmat)
      kernel.local[, ,i] <-  Sigma
    }
    
    KERNEL_LIST[[2]] <- kernel.local
    
    FIN_Sigma.mat <- list()
    dist0 <- list()
    
    for(KK in 1:2){
      Sigma.mat <- matrix(rep(NA, (N*t)^2), nrow = N*t)
      Q.mat <- matrix(rep(NA, (N*t)^2), nrow = N*t)
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
    
    nu=p[7:8]
    beta=p[9]
    rot=p[10]
    Beta=matrix(0,q,q)
    diag(Beta)=1
    Beta[2,1]=rot
    Beta[1,2]=rot
    var =p[11:12]
    
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
    
    possibleError <- tryCatch(
      cov.chol <- chol(S),
      error=function(e) e
    )
    loglikelihood <-0
    
    # Check for error before calculating the likelihood
    if(inherits(possibleError, "error")){
      loglikelihood <- 1e+06
    } else{
      
      for(tt in 1:t){
        tmp1 <- backsolve(cov.chol,  matrix(c(data[((tt-1)*nrow(data)/(q*t)+1):(tt*nrow(data)/(q*t)),],data[(nrow(data)/q+(tt-1)*nrow(data)/(q*t)+1):(nrow(data)/q+tt*nrow(data)/(q*t)),]),ncol=1), transpose = TRUE)
        #tmp1 <- backsolve(cov.chol, c(data)-c(rep(b0[1],dim(data)[1]/2),rep(b0[2],dim(data)[1]/2)), transpose = TRUE)
        
        ResCinvRes <- t(tmp1) %*% tmp1
        loglikelihood <- loglikelihood + q * sum(log(diag(cov.chol))) + 0.5 * sum(diag(ResCinvRes))
      }
    }
    if (abs(loglikelihood) == Inf) { loglikelihood <- 1e+06 } # Make sure not Inf
    return(loglikelihood)
  }
}
