empirical_st_cov <- function( data1, data2 = NULL, cross, locations, max_time_lag ){
  
  # input: data matrix with temporal replicates as rows and each column represents the locations
  # cross: TRUE or FALSE
  
  # output: correlation
  
  mean = colMeans(data1)
  scov = matrix(c(0), ncol(data1), ncol(data1))
  n_stations = ncol(data1)
  
  emp <- list()
  
  if(!cross == TRUE){
    
    for(tau in 0:max_time_lag){
      
      for(i in (tau+1):(nrow(data1)-tau)){
        scov = scov + (data1[i,] - mean)%*%t(data1[i+tau,] - mean)
      }
      
      scov = scov/(nrow(data1) - tau)
      
      sd = apply(data1, 2, sd)
      
      Dz = sd*diag(n_stations)
      
      scor = solve(Dz)%*%scov%*%solve(Dz)
      
      ## Arrange the empirical covariance values in one dataframe/matrix. 
      
      loc <- 1
      xlag <- locations[loc,1]-locations[,1]
      ylag <- locations[loc,2]-locations[,2]
      emp_vals <- scor[,loc]
      
      for(loc in 2:ncol(data1)){
        xlag <- c(xlag, locations[loc,1] - locations[,1])
        ylag <- c(ylag, locations[loc,2] - locations[,2])
        emp_vals <- c(emp_vals,scor[,loc])
      }
      emp_temp <- data.frame(xlag,ylag,emp_vals)
      colnames(emp_temp) <- c('xlag', 'ylag', 'correlation')
      emp[[tau+1]] <- emp_temp
    }
    
  }else{
    
    for(tau in 0:max_time_lag){
      
      mean2 = colMeans(data2)
      
      for(i in (tau+1):(nrow(data1)-tau)){
        scov = scov + (data1[i,] - mean)%*%t(data2[i+tau,] - mean2)
        
      }
      scov = scov/(nrow(data1) - tau)
      
      sd = apply(data1, 2, sd)
      sd2 = apply(data2, 2, sd)
      
      Dz = sd*diag(n_stations)
      Dz2 = sd2*diag(n_stations)
      
      scor = solve(Dz)%*%scov%*%solve(Dz2)
      
      ## Arrange the empirical covariance values in one dataframe/matrix. 
      
      loc <- 1
      xlag <- locations[loc,1]-locations[,1]
      ylag <- locations[loc,2]-locations[,2]
      emp_vals <- scor[,loc]
      
      for(loc in 2:ncol(data1)){
        xlag <- c(xlag, locations[loc,1] - locations[,1])
        ylag <- c(ylag, locations[loc,2] - locations[,2])
        emp_vals <- c(emp_vals,scor[,loc])
      }
      emp_temp <- data.frame(xlag,ylag,emp_vals)
      colnames(emp_temp) <- c('xlag', 'ylag', 'crosscorrelation')
      emp[[tau+1]] <- emp_temp
    }
  }
  
  return(emp)
}

empirical_covariance_dataframe <- function(data1_cov, data2_cov, cross_cov){
  
  empirical_var1 <- data1_cov[[1]]
  empirical_var2 <- data2_cov[[1]]
  empirical_var3 <- cross_cov[[1]]
  
  binned.1 <- empirical_var1 %>% group_by(xlag,ylag) %>% summarize(avg1=mean(correlation))
  binned.2 <- empirical_var2 %>% group_by(xlag,ylag) %>% summarize(avg1=mean(correlation))
  binned.3 <- empirical_var3 %>% group_by(xlag,ylag) %>% summarize(avg1=mean(crosscorrelation))
  
  binned_orig <- cbind(binned.1$xlag,binned.1$ylag,rep(0,nrow(binned.1)),binned.1$avg1,
                       binned.2$avg1,binned.3$avg1)
  
  if(length(cross) > 1){
    for (i in 2:length(cross)){
      empirical_var1 <- data1_cov[[i]]
      empirical_var2 <- data2_cov[[i]]
      empirical_var3 <- cross_cov[[i]]
      
      binned.1 <- empirical_var1 %>% group_by(xlag,ylag) %>% summarize(avg1=mean(correlation))
      binned.2 <- empirical_var2 %>% group_by(xlag,ylag) %>% summarize(avg1=mean(correlation))
      binned.3 <- empirical_var3 %>% group_by(xlag,ylag) %>% summarize(avg1=mean(crosscorrelation))
      
      binned_orig <- rbind(binned_orig,cbind(binned.1$xlag,binned.1$ylag,rep(i-1,nrow(binned.1)),binned.1$avg1,
                                             binned.2$avg1,binned.3$avg1))
    }
  }
  colnames(binned_orig) <- c('xlag', 'ylag', 'timelag', 'var1_cor', 'var2_cor', 'cross_cor')
  
  return(binned_orig)
}

empirical_st_cov_nonstationary <- function( data1, data2 = NULL, cross, locations, max_time_lag){
  
  # input: data matrix with temporal replicates as rows and each column represents the locations
  # cross: TRUE or FALSE
  
  # output: correlation
  
  mean = colMeans(data1)
  scov = matrix(c(0), ncol(data1), ncol(data1))
  n_stations = ncol(data1)
  
  emp <- list()
  
  if(!cross == TRUE){
    
    for(tau in 0:max_time_lag){
      
      for(i in (tau+1):(nrow(data1)-tau)){
        scov = scov + (data1[i,] - mean)%*%t(data1[i+tau,] - mean)
      }
      
      scov = scov/(nrow(data1) - tau)
      
      sd = apply(data1, 2, sd)
      
      Dz = sd*diag(n_stations)
      
      scor = solve(Dz)%*%scov%*%solve(Dz)
      
      ## Arrange the empirical covariance values in one dataframe/matrix. 
      
      loc <- 1
      xlag <- locations[loc,1]-locations[,1]
      ylag <- locations[loc,2]-locations[,2]
      emp_vals <- scor[,loc]
      
      for(loc in 2:ncol(data1)){
        xlag <- c(xlag, locations[loc,1] - locations[,1])
        ylag <- c(ylag, locations[loc,2] - locations[,2])
        emp_vals <- c(emp_vals,scor[,loc])
      }
      emp_temp <- data.frame(xlag,ylag,emp_vals)
      colnames(emp_temp) <- c('xlag', 'ylag', 'correlation')
      emp[[tau+1]] <- emp_temp
    }
    
  }else{
    
    for(tau in 0:max_time_lag){
      
      mean2 = colMeans(data2)
      
      for(i in (tau+1):(nrow(data1)-tau)){
        scov = scov + (data1[i,] - mean)%*%t(data2[i+tau,] - mean2)
        
      }
      scov = scov/(nrow(data1) - tau)
      
      sd = apply(data1, 2, sd)
      sd2 = apply(data2, 2, sd)
      
      Dz = sd*diag(n_stations)
      Dz2 = sd2*diag(n_stations)
      
      scor = solve(Dz)%*%scov%*%solve(Dz2)
      
      ## Arrange the empirical covariance values in one dataframe/matrix.
      
      loc <- 1
      xlag <- locations[loc,1]-locations[,1]
      ylag <- locations[loc,2]-locations[,2]
      emp_vals <- scor[,loc]
      
      for(loc in 2:ncol(data1)){
        xlag <- c(xlag, locations[loc,1] - locations[,1])
        ylag <- c(ylag, locations[loc,2] - locations[,2])
        emp_vals <- c(emp_vals,scor[,loc])
      }
      emp_temp <- data.frame(xlag,ylag,emp_vals)
      colnames(emp_temp) <- c('xlag', 'ylag', 'crosscorrelation')
      emp[[tau+1]] <- emp_temp
    }
  }
  
  return(emp)
}

empirical_covariance_dataframe_nonstationary <- function(data1_cov, data2_cov, cross_cov, params){
  
  empirical_var1 <- data1_cov[[1]]
  empirical_var2 <- data2_cov[[1]]
  empirical_var3 <- cross_cov[[1]]
  
  binned.1 <- empirical_var1 %>% group_by(xlag,ylag) %>% summarize(avg1=mean(correlation))
  binned.2 <- empirical_var2 %>% group_by(xlag,ylag) %>% summarize(avg1=mean(correlation))
  binned.3 <- empirical_var3 %>% group_by(xlag,ylag) %>% summarize(avg1=mean(crosscorrelation))
  
  binned_orig <- cbind(binned.1$xlag,binned.1$ylag,rep(0,nrow(binned.1)),binned.1$avg1,
                       binned.2$avg1,binned.3$avg1)
  
  if(length(cross) > 1){
    for (i in 2:length(cross)){
      empirical_var1 <- data1_cov[[i]]
      empirical_var2 <- data2_cov[[i]]
      empirical_var3 <- cross_cov[[i]]
      
      binned.1 <- empirical_var1 %>% group_by(xlag,ylag) %>% summarize(avg1=mean(correlation))
      binned.2 <- empirical_var2 %>% group_by(xlag,ylag) %>% summarize(avg1=mean(correlation))
      binned.3 <- empirical_var3 %>% group_by(xlag,ylag) %>% summarize(avg1=mean(crosscorrelation))
      
      binned_orig <- rbind(binned_orig,cbind(binned.1$xlag,binned.1$ylag,rep(i-1,nrow(binned.1)),binned.1$avg1,
                                             binned.2$avg1,binned.3$avg1))
    }
  }
  
  locations <- cbind(binned.1$xlag, binned.1$ylag)
  
  lam1.1 <- params[1]
  lam1.2 <- params[2]
  eta.1 <- params[3]
  
  lam2.1 <- params[4]
  lam2.2 <- params[5]
  eta.2 <- params[6]
  
  Pmat1 <- matrix(c(cos(eta.1), -sin(eta.1), sin(eta.1), cos(eta.1)), nrow = 2, byrow = T)
  Dmat1 <- diag(c(lam1.1, lam1.2))
  Sigma1 <- Pmat1 %*% Dmat1 %*% t(Pmat1)
  distances1 <- mahalanobis(locations, c(0,0), Sigma1)
  
  Pmat2 <- matrix(c(cos(eta.2), -sin(eta.2), sin(eta.2), cos(eta.2)), nrow = 2, byrow = T)
  Dmat2 <- diag(c(lam2.1, lam2.2))
  Sigma2 <- Pmat2 %*% Dmat2 %*% t(Pmat2)
  distances2 <- mahalanobis(locations, c(0,0), Sigma2)
  
  distances3 <- mahalanobis(locations, c(0,0), 0.5*(Sigma2+Sigma2))
  
  binned_orig <- cbind(binned_orig, distances1, distances2, distances3)
  
  colnames(binned_orig) <- c('xlag', 'ylag', 'timelag', 'var1_cor', 'var2_cor', 'cross_cor', 'dist1', 'dist2', 'dist3')
  
  return(binned_orig)
}
