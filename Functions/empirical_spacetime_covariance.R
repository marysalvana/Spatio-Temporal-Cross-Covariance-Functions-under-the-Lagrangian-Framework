#---------STATIONARY---------#

fitting <- function(coordinates, obs1, obs2 = NULL){
  # FORMAT OF INPUT VALUES
  # coordinates: two-column matrix with longitude in column 1 and latitude in column 2
  # obs1: T-column matrix of observations where a row represents observations of one location and the columns are hourly/monthly/yearly observations
  # obs2: another covariate
  
  data_matrix <- data_format_into_matrix(data1 = obs1, data2 = obs2, temporal_replicates = ncol(obs1), simulated = F)
  conso_cor <- empirical_st_cov(data1 = t(data_matrix[insample_loc_index, 1:(ncol(data_matrix)/2)]), data2 = t(data_matrix[insample_loc_index, (ncol(data_matrix)/2 + 1):ncol(data_matrix)]), locations = coordinates[insample_loc_index,], max_time_lag = 5, simulated = F)
  binned <- empirical_covariance_dataframe(data1_cov = var1_cov, data2_cov = var2_cov, cross_cov = cross, simulated = F)
  hlag <- sqrt(binned[which(binned[,3] == 0), 1]^2 + binned[which(binned[,3] == 0), 2]^2)
  #display plots
  par(mfrow = c(1,3))
  plot(hlag/1000, binned[which(binned[,3]==0), 4], pch=3, ylab='', col=1,xlab='Spatial Lag (km)', main='', col.main= "#4EC1DE", ylim=c(0,1))
  plot(hlag/1000, binned[which(binned[,3]==0), 5], pch=3, ylab='', col=1,xlab='Spatial Lag (km)', main='', col.main= "#4EC1DE", ylim=c(0,1))
  plot(hlag/1000, binned[which(binned[,3]==0), 6], pch=3, ylab='', col=1,xlab='Spatial Lag (km)', main='', col.main= "#4EC1DE", ylim=c(0,1))
  
}

empirical_st_cov <- function(data1, data2, locations, max_time_lag, simulated, p = 2){
  
  # input: data matrix with temporal replicates as rows and each column represents the locations
  # cross: TRUE or FALSE
  # p = 2: bivariate
  
  # output: correlation
  
  emp <- list()
  
  if(!simulated == TRUE){
    
    mean <- colMeans(data1)
    mean2 <- colMeans(data2)
    scov <- scov2 <- scov3 <- matrix(c(0), ncol(data1), ncol(data1))
    n_stations = ncol(data1)
    
    for(tau in 0:max_time_lag){
      
      for(i in (tau+1):(nrow(data1)-tau)){
        scov = scov + (data1[i,] - mean)%*%t(data1[i + tau,] - mean)
        scov2 = scov2 + (data2[i,] - mean2)%*%t(data2[i + tau,] - mean2)
        scov3 = scov3 + (data1[i,] - mean)%*%t(data2[i + tau,] - mean2)
      }
      
      scov <- scov/(nrow(data1) - tau)
      scov2 <- scov2/(nrow(data1) - tau)
      scov3 <- scov3/(nrow(data1) - tau)
      
      sd = apply(data1, 2, sd)
      sd2 = apply(data2, 2, sd)
      
      Dz = sd*diag(n_stations)
      Dz2 = sd2*diag(n_stations)
      
      scor1 = solve(Dz)%*%scov%*%solve(Dz)
      scor2 = solve(Dz2)%*%scov2%*%solve(Dz2)
      scor3 = solve(Dz)%*%scov3%*%solve(Dz2)
      
      ## Arrange the empirical covariance values in one dataframe/matrix. 
      
      loc <- 1
      xlag <- locations[loc,1]-locations[,1]
      ylag <- locations[loc,2]-locations[,2]
      emp_vals <- scor[,loc]
      emp_vals <- cbind(scor1[,loc], scor2[,loc], scor3[,loc])
      
      for(loc in 1:nrow(locations)){
        xlag <- c(xlag, locations[loc,1] - locations[,1])
        ylag <- c(ylag, locations[loc,2] - locations[,2])
        emp_vals <- rbind(emp_vals, cbind(scor1[,loc], scor2[,loc], scor3[,loc]))
      }
      
      empirical <- data.frame(xlag, ylag, emp_vals)
      colnames(empirical) <- c('xlag', 'ylag', 'var1_cor', 'var2_cor', 'cross_cor')
      
      emp[[tau + 1]] <- empirical
    }
    
  }else{
    
    empcor <- cor(data1)
    
    for(tau in 0:max_time_lag){
    
      empcov <- list()
    
      empcov[[1]] <- empcor[1:nrow(locations),(tau*nrow(locations) + 1):((tau + 1)*nrow(locations))]
      empcov[[2]] <- empcor[(max_time_lag + 1)*nrow(locations) + 1:nrow(locations),(max_time_lag + 1)*nrow(locations) + (tau*nrow(locations) + 1):((tau + 1)*nrow(locations))]
      empcov[[3]] <- empcor[1:nrow(locations),(max_time_lag + 1)*nrow(locations) + (tau*nrow(locations) + 1):((tau + 1)*nrow(locations))]
      
      loc <- 1
      xlag <- locations[loc,1]-locations[,1]
      ylag <- locations[loc,2]-locations[,2]
      emp_vals <- cbind(empcov[[1]][, loc], empcov[[2]][, loc], empcov[[3]][, loc])
      
      for(loc in 1:nrow(locations)){
        xlag <- c(xlag, locations[loc,1] - locations[,1])
        ylag <- c(ylag, locations[loc,2] - locations[,2])
        emp_vals <- rbind(emp_vals, cbind(empcov[[1]][, loc], empcov[[2]][, loc], empcov[[3]][, loc]))
      }
      
      empirical <- data.frame(xlag, ylag, emp_vals)
      colnames(empirical) <- c('xlag', 'ylag', 'var1_cor', 'var2_cor', 'cross_cor')
      
      emp[[tau + 1]] <- empirical
    }
  }
  
  return(emp)
}

empirical_covariance_dataframe <- function(data1_cov, data2_cov = NULL, cross_cov = NULL, simulated = F){
  
  if(!simulated == TRUE){
    empirical_var1 <- data1_cov[[1]]
    empirical_var2 <- data2_cov[[1]]
    empirical_var3 <- cross_cov[[1]]
    
    binned.1 <- empirical_var1 %>% group_by(xlag, ylag) %>% summarize(avg1=mean(correlation))
    binned.2 <- empirical_var2 %>% group_by(xlag, ylag) %>% summarize(avg1=mean(correlation))
    binned.3 <- empirical_var3 %>% group_by(xlag, ylag) %>% summarize(avg1=mean(crosscorrelation))
    
    binned_orig <- cbind(binned.1$xlag, binned.1$ylag, rep(0,nrow(binned.1)), binned.1$avg1,
                         binned.2$avg1, binned.3$avg1)
    
    if(length(data1_cov) > 1){
      for (i in 2:length(data1_cov)){
        empirical_var1 <- data1_cov[[i]]
        empirical_var2 <- data2_cov[[i]]
        empirical_var3 <- cross_cov[[i]]
        
        binned.1 <- empirical_var1 %>% group_by(xlag, ylag) %>% summarize(avg1=mean(correlation))
        binned.2 <- empirical_var2 %>% group_by(xlag, ylag) %>% summarize(avg1=mean(correlation))
        binned.3 <- empirical_var3 %>% group_by(xlag, ylag) %>% summarize(avg1=mean(crosscorrelation))
        
        binned_orig <- rbind(binned_orig,cbind(binned.1$xlag,binned.1$ylag,rep(i-1,nrow(binned.1)),binned.1$avg1,
                                               binned.2$avg1,binned.3$avg1))
      }
    }
    colnames(binned_orig) <- c('xlag', 'ylag', 'tlag', 'var1_cor', 'var2_cor', 'cross_cor')
  }else{
    empirical_var1 <- data1_cov[[1]]
    
    binned.1 <- empirical_var1 %>% group_by(xlag, ylag) %>% summarize(avg1=mean(var1_cor))
    binned.2 <- empirical_var1 %>% group_by(xlag, ylag) %>% summarize(avg1=mean(var2_cor))
    binned.3 <- empirical_var1 %>% group_by(xlag, ylag) %>% summarize(avg1=mean(cross_cor))
    
    binned_orig <- cbind(binned.1$xlag, binned.1$ylag, rep(0,nrow(binned.1)), binned.1$avg1,
                         binned.2$avg1, binned.3$avg1)
    
    if(length(data1_cov) > 1){
      for (i in 2:length(data1_cov)){
        empirical_var1 <- data1_cov[[i]]
        
        binned.1 <- empirical_var1 %>% group_by(xlag, ylag) %>% summarize(avg1=mean(var1_cor))
        binned.2 <- empirical_var1 %>% group_by(xlag, ylag) %>% summarize(avg1=mean(var2_cor))
        binned.3 <- empirical_var1 %>% group_by(xlag, ylag) %>% summarize(avg1=mean(cross_cor))
        
        binned_orig <- rbind(binned_orig,cbind(binned.1$xlag,binned.1$ylag,rep(i-1,nrow(binned.1)),binned.1$avg1,
                                               binned.2$avg1,binned.3$avg1))
      }
    }
    colnames(binned_orig) <- c('xlag', 'ylag', 'tlag', 'var1_cor', 'var2_cor', 'cross_cor')
  }
  return(binned_orig)
}


#---------NONSTATIONARY---------#

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
