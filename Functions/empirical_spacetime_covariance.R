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
  
  dist0 <- sqrt(binned_orig[,1]^2+binned_orig[,2]^2)/1000
  
  for (i in 2:6){
    empirical_var1 <- data1_cov[[i]]
    empirical_var2 <- data2_cov[[i]]
    empirical_var3 <- cross_cov[[i]]
    
    binned.1 <- empirical_var1 %>% group_by(xlag,ylag) %>% summarize(avg1=mean(correlation))
    binned.2 <- empirical_var2 %>% group_by(xlag,ylag) %>% summarize(avg1=mean(correlation))
    binned.3 <- empirical_var3 %>% group_by(xlag,ylag) %>% summarize(avg1=mean(crosscorrelation))
    
    binned_orig <- rbind(binned_orig,cbind(binned.1$xlag,binned.1$ylag,rep(i-1,nrow(binned.1)),binned.1$avg1,
                                           binned.2$avg1,binned.3$avg1))
  }
  colnames(binned_orig) <- c('xlag', 'ylag', 'timelag', 'var1_cor', 'var2_cor', 'cross_cor')
  
  return(binned_orig)
}
