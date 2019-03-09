#---------STATIONARY---------#

fitting <- function(coordinates, obs1, obs2 = NULL){
  # FORMAT OF INPUT VALUES
  # coordinates: two-column matrix with longitude in column 1 and latitude in column 2
  # obs1: T-column matrix of observations where a row represents observations of one location and the columns are hourly/monthly/yearly observations
  # obs2: another covariate
  
  data_matrix <- data_format_into_matrix(data1 = obs1, data2 = obs2, temporal_replicates = ncol(obs1), simulated = F)
  conso_cor <- empirical_st_cov(data1 = t(data_matrix[insample_loc_index, 1:(ncol(data_matrix)/2)]), data2 = t(data_matrix[insample_loc_index, (ncol(data_matrix)/2 + 1):ncol(data_matrix)]), locations = coordinates[insample_loc_index,], max_time_lag = 5, simulated = F)
  binned <- empirical_covariance_dataframe(data1_cov = var1_cov, data2_cov = var2_cov, cross_cov = cross, simulated = F)
  
  theta_init <- c(3.218, 3.736, 794.8, 1, 1, max(binned[which(binned[,3] == 0), 6]) - 0.001)
  emp_vel <- which.max(binned[which(binned[,3] == 1), 4])
  w_init <- binned[which(binned[,3] == 1)[emp_vel], 1:2]
  mod <- fit_model(init = theta_init, wind_init = w_init, mod = 'matern', randvel = F, weight = 3, empcov_spatial = binned[which(binned[,3]==0),], empcov_st = binned[which(binned[,3] > 0),], nug_eff = F, meters = T, num_iter = 10)
  
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

simulation_study <- function(true_param_spatial, true_param_velocity, sim_model, rand_vel, num_sim, max_u, num_variables, location, plot = T, nugget = F){
  
  mod_params <- matrix(, ncol = length(true_param_spatial) + length(true_param_velocity), nrow = num_sim)
  
  if(sim_model == 'matern' & rand_vel == T){
    sim.cov <- simulate_model(mod = sim_model, randvel = rand_vel, theta = true_param_spatial, wind = true_param_velocity[1:2], wind_var = matrix(c(true_param_velocity[3:4], true_param_velocity[4:5]), ncol=2), maxtimelag = max_u, p = num_variables, locations = location, meters = T, nugeff = nugget)
  }else if(sim_model == 'matern' & rand_vel == F){
    sim.cov <- simulate_model(mod = sim_model, randvel = rand_vel, theta = true_param_spatial, wind = true_param_velocity, maxtimelag = max_u, p = num_variables, locations = location, meters = T, nugeff = nugget)
  }else if(sim_model == 'lmc' & rand_vel == T){
    sim.cov <- simulate_model(mod = sim_model, randvel = rand_vel, theta = true_param_spatial, wind = true_param_velocity[1:4], wind_var = list(matrix(c(true_param_velocity[5:6], true_param_velocity[6:7]), ncol=2), matrix(c(true_param_velocity[8:9], true_param_velocity[9:10]), ncol=2)), maxtimelag = max_u, p = num_variables, locations = location, meters = T, nugeff = nugget)
  }else if(sim_model == 'lmc' & rand_vel == F){
    sim.cov <- simulate_model(mod = sim_model, randvel = rand_vel, theta = true_param_spatial, wind = true_param_velocity, maxtimelag = max_u, p = num_variables, locations = location, meters = T, nugeff = nugget)
  }
  
  for(iter in 1:num_sim){
    
    A <- mvrnorm(n = 1000, mu = rep(0, dim(sim.cov)[1]), Sigma = sim.cov)
    A1 <- A[, -to_remove]
    
    conso_cor <- empirical_st_cov(data1 = A1, locations = location[insample_loc_index, ], max_time_lag = max_u, simulated = T)  
    binned <- empirical_covariance_dataframe(data1_cov = conso_cor, simulated = T)
    
    if(plot == T){
      if(iter == 1){
        hlag <- sqrt(binned[which(binned[,3]==0), 1]^2 + binned[which(binned[,3]==0), 2]^2)
        #display plots
        par(mfrow = c(1,3))
        plot(hlag/1000, binned[which(binned[,3]==0), 4], pch=3, ylab='', col=1,xlab='Spatial Lag (km)', main='', col.main= "#4EC1DE", ylim=c(0,1))
        plot(hlag/1000, binned[which(binned[,3]==0), 5], pch=3, ylab='', col=1,xlab='Spatial Lag (km)', main='', col.main= "#4EC1DE", ylim=c(0,1))
        plot(hlag/1000, binned[which(binned[,3]==0), 6], pch=3, ylab='', col=1,xlab='Spatial Lag (km)', main='', col.main= "#4EC1DE", ylim=c(0,1))
      }
    }
    
    if(sim_model == 'matern'){
      theta_init <- c(3.218, 3.736, 794.8, 1, 1, max(binned[which(binned[,3] == 0), 6])-0.001) #change this values to empirical
      emp_vel <- which.max(binned[which(binned[,3] == 1), 4])
      
      if(rand_vel == T){
        w_init <- c(1 - binned[which(binned[,3] == 1)[emp_vel], 4], (1 - binned[which(binned[,3] == 1)[emp_vel], 4])/10, 1 - binned[which(binned[,3] == 1)[emp_vel], 4], binned[which(binned[,3] == 1)[emp_vel], 1:2])
        mod <- fit_model(init = theta_init, wind_init = w_init, mod = sim_model, randvel = rand_vel, weight = 3, empcov_spatial = binned[which(binned[,3]==0),], empcov_st = binned[which(binned[,3] == 1),], nug_eff = nugget, meters = T, num_iter = 0)
      }else if(rand_vel == F){
        w_init <- binned[which(binned[,3] == 1)[emp_vel], 1:2]
        mod <- fit_model(init = theta_init, wind_init = w_init, mod = sim_model, randvel = rand_vel, weight = 3, empcov_spatial = binned[which(binned[,3]==0),], empcov_st = binned[which(binned[,3] > 0),], nug_eff = nugget, meters = T, num_iter = 10)
      }
    }else if(sim_model == 'lmc'){
      theta_init <- c(3, 4, 300, 705, 0.99, 0.99, 0.838, 0.545, 0.0001, 0.999)
      if(rand_vel == T){
        w_init <- c(-2600000, -661200, 3604000, 1947000, 100, 0.00009, 100, 100, 0.00009, 100)
        mod <- fit_model(init = theta_init, wind_init = w_init, mod = sim_model, randvel = rand_vel, weight = 3, empcov_spatial = binned[which(binned[,3]==0),], empcov_st = binned[which(binned[,3] == 1),], nug_eff = nugget, meters = T, num_iter = 0)
      }else if(rand_vel == F){
        w_init <- c(-2600000, -661200, 3604000, 1947000)
        mod <- fit_model(init = theta_init, wind_init = w_init, mod = sim_model, randvel = rand_vel, weight = 3, empcov_spatial = binned[which(binned[,3]==0),], empcov_st = binned[which(binned[,3] > 0),], nug_eff = nugget, meters = T, num_iter = 10)
      }
    }
    mod_params[iter, ] <- mod$parameters
  }
  return(mod_params)
}

simulate_model <- function(mod, theta, wind, wind_var = NULL, maxtimelag, p = 2, locations, meters = T, nugeff, randvel){
  if(mod == 'matern' & randvel == T){
    cov.mod <- matern_random_cov(theta, wind, wind_var, max_time_lag = maxtimelag, q = p, new_locations = locations, nug_eff = nugeff )
  }else if(mod == 'matern' & randvel == F){
    cov.mod <- matern_cov(theta, wind, max_time_lag = maxtimelag, q = p, new_locations = locations, nug_eff = nugeff )
  }else if(mod == 'lmc' & randvel == T){
    cov.mod <- lmc_random_cov(theta, wind, wind_var, max_time_lag = maxtimelag, q = p, new_locations = locations, nug_eff = nugeff )
  }else if(mod == 'lmc' & randvel == F){
    cov.mod <- lmc_cov(theta, wind, max_time_lag = maxtimelag, q = p, new_locations = locations, nug_eff = nugeff )
  }
  return(cov.mod)
}

fit_model <- function(init = NULL, wind_init, mod, randvel, weight, empcov_spatial = NULL, empcov_st, nug_eff = F, meters = T, num_iter, 
                      est_param.temp = NULL, aniso = F, rotation_matrix = NULL){
  
  # num_iter : number of loops to run optim
  
  if(mod == 'matern' | mod == 'gneitingmatern'){
    fit1.mod <- optim(par = init[-length(init)], wls, emp_cov1 = empcov_spatial, nug_eff = F, meters = T, weights = weight, step = 1, aniso = F, model = mod, control=list(maxit = 10000, parscale = init[-length(init)], trace = 5))
    fit2.mod <- optim(par = init[length(init)], wls, emp_cov1 = empcov_spatial, nug_eff = F, meters = T, weights = weight, step = 2, est_param = fit1.mod$par, aniso = F, model = mod, method='SANN', control = list(maxit = 3000, parscale = init[length(init)], trace = 5))
    
    if(mod == 'matern'){
      if(randvel == T){
        
        fit3.mod <- optim(par = wind_init, wls, emp_cov1 = empcov_st, nug_eff = F, meters = T, weights = weight, step = 3, est_param = c(fit1.mod$par, fit2.mod$par), aniso = F, rand.vel = T, model = mod, control = list(maxit = 10000, parscale = wind_init, trace = 5))
        
        lst <- list(parameters = c(est_param.temp, fit3.mod$par), fn_value = fit1.mod$value + fit2.mod$value + fit3.mod$value)
        
        return(lst)
      }else if(randvel == F){
        
        fit3.mod <- optim(par = wind_init, wls, emp_cov1 = empcov_st, nug_eff = F, meters = T, weights = weight, step = 3, est_param = c(fit1.mod$par, fit2.mod$par), aniso = F, rand.vel = F, model = mod, control = list(maxit = 10000, parscale = wind_init, trace = 5))
        
        if(num_iter > 0){
          for(iter3 in 1:num_iter){
            new_wind_init <- fit3.mod$par
            fit3.mod <- optim(par = new_wind_init, wls, emp_cov1 = empcov_st, nug_eff = F, meters = T, weights = weight, step = 3, est_param = c(fit1.mod$par, fit2.mod$par), aniso = F, rand.vel = F, model = mod, control = list(maxit = 10000, parscale = new_wind_init, trace = 5))
          }
        }
        
        lst <- list(parameters = c(fit1.mod$par, fit2.mod$par, fit3.mod$par), fn_value = fit1.mod$value + fit2.mod$value + fit3.mod$value)
        
        return(lst)
        
      }
    }else if(mod == 'gneitingmatern'){
      fit3.mod <- optim(par = gneitingmatern_init, wls, emp_cov1 = empcov_st, nug_eff = F, meters = T, weights = weight, step = 3, est_param = c(fit1.mod$par, fit2.mod$par), aniso = F, model = mod, control = list(maxit = 10000, parscale = gneitingmatern_init, trace = 5))
      
      lst <- list(parameters = c(est_param.temp, fit3.mod$par), fn_value = fit1.mod$value + fit2.mod$value + fit3.mod$value)
      
      return(lst)
    }
  }else if(mod == 'lmc'){
    if(randvel == T){
      fit1.mod <- optim(par = init, wls, emp_cov1 = empcov_spatial, nug_eff = F, meters = T, weights = weight, step = 1, aniso = F, model = mod, control=list(maxit = 10000, parscale = init, trace = 5))
      fit2.mod <- optim(par = wind_init, wls, emp_cov1 = empcov_st, nug_eff = F, meters = T, weights = weight, step = 2, est_param = fit1.mod$par, aniso = F, rand.vel = randvel, model = mod, control = list(maxit = 10000,parscale = wind_init, trace = 5))
      
      lst <- list(parameters = c(fit1.mod$par, fit2.mod$par), fn_value = fit1.mod$value + fit2.mod$value)
      
      return(lst)
    }else if(randvel == F){
      fit1.mod <- optim(par = init, wls, emp_cov1 = empcov_spatial, nug_eff = F, meters = T, weights = weight, step = 1, aniso = F, model = mod, control=list(maxit = 10000, parscale = init, trace = 5))
      fit2.mod <- optim(par = wind_init, wls, emp_cov1 = empcov_st, nug_eff = F, meters = T, weights = weight, step = 2, est_param = fit1.mod$par, aniso = F, rand.vel = randvel, model = mod, control = list(maxit = 10000,parscale = wind_init, trace = 5))
      
      lst <- list(parameters = c(fit1.mod$par, fit2.mod$par), fn_value = fit1.mod$value + fit2.mod$value)
      
      return(lst)
    }
  }
}

wls <- function(theta, emp_cov1, weights, nug_eff, step, est_param = NULL, meters, aniso = F, model, rand.vel) {
  
  # nug_eff: T means we also need to estimate nugget effect
  
  loss <- 0
  
  if(!aniso == T){
    R <- diag(2)
  }else{
    R <- rotation_matrix
  }
  
  #R <- matrix(c(-1.1339125 *cos(116.9811232),-1.1339125 *sin(116.9811232),0.8617141*sin(116.9811232),
  #              -0.8617141*cos(116.9811232)),ncol=2,byrow=T)
  
  if(meters == T){
    h <- sqrt(emp_cov1[,1]^2 + emp_cov1[,2]^2)/1000
  }else{
    h <- sqrt(emp_cov1[,1]^2 + emp_cov1[,2]^2)
  }
  
  np <- ifelse(h != 0, 1/h, 1)
  
  if( model == 'matern'){
    if(step == 1){
      if(nug_eff == T){
        nug <- theta[6:7]
      }else{
        nug <- c(0, 0)
      }
      nu <- theta[1:2]
      beta <- theta[3]
      var <- theta[4:5]
      
      if( theta[1] < 0.0001 | theta[2] < 0.0001 | theta[3] < 0.0001 | theta[4] < 0.0001 | theta[5] < 0.0001 ){
        return(Inf)
      }else{
        for(i in 1:2){
          
          theo <- ifelse(h != 0, var[i]*(h/beta)^nu[i]*besselK(h/beta, nu[i])/(2^(nu[i]-1)*gamma(nu[i])), var[i]+nug[i])
          if (weights == 1) 
            tloss <- sum((emp_cov1[,i+3] - theo)^2)
          if (weights == 2) 
            tloss <- sum(np* (emp_cov1[,i+3] - theo)^2)
          if (weights == 3) 
            tloss <- sum(((emp_cov1[,i+3] - theo)/(1.001 - theo))^2)
          
          loss <- loss + tloss
        }
        return(loss)
      }
    }else if(step == 2){
      
      if(nug_eff == T){
        nug <- est_param[6:7]
      }else{
        nug <- c(0, 0)
      }
      nu <- est_param[1:2]
      beta <- est_param[3]
      var <- est_param[4:5]
      rho <- theta
      
      nu3 <- (nu[1]+nu[2])/2
      
      if( abs(rho) > 1 ){
        return(Inf)
      }else{
        
        theo <- ifelse(h != 0, rho*sqrt(var[1]*var[2])*(h/beta)^nu3*besselK(h/beta, nu3)/(2^(nu3-1)*gamma(nu3)), rho*sqrt(var[1]*var[2]))
        if (weights == 1) 
          tloss <- sum((theo - emp_cov1[,6])^2)
        if (weights == 2) 
          tloss <- sum(np* (emp_cov1[,6] - theo)^2)
        if (weights == 3) 
          tloss <- sum(((emp_cov1[,6] - theo)/(max(emp_cov1[,6]) - theo))^2)
        loss <- loss + tloss
        
        return(loss)
      }
    }else{
      if(nug_eff == T){
        nug <- est_param[6:7]
        rho <- est_param[8]
      }else{
        nug <- c(0, 0)
        rho <- est_param[6]
      }
      nu <- est_param[1:2]
      beta <- est_param[3]
      var <- est_param[4:5]
      
      nu1 <- nu[1]
      nu2 <- nu[2]
      nu3 <- (nu[1]+nu[2])/2
      
      if(!rand.vel == T){
        w <- theta
        
        if(meters == T){
          h <- sqrt((emp_cov1[,1] - emp_cov1[,3]*w[1])^2 + (emp_cov1[,2] - emp_cov1[,3]*w[2])^2)/1000
        }else{
          h <- sqrt((emp_cov1[,1] - emp_cov1[,3]*w[1])^2 + (emp_cov1[,2] - emp_cov1[,3]*w[2])^2)
        }
        
        np <- ifelse(h != 0, 1/h, 1)
        
        for(i in 1:2){
          
          theo <- ifelse(h!=0, var[i]*(h/beta)^nu[i]*besselK(h/beta,nu[i])/(2^(nu[i]-1)*gamma(nu[i])), var[i]+nug[i])
          if (weights == 1) 
            tloss <- sum((theo - emp_cov1[,i+3])^2)
          if (weights == 2) 
            tloss <- sum(np* (emp_cov1[,i+3] - theo)^2)
          if (weights == 3) 
            tloss <- sum(((emp_cov1[,i+3] - theo)/(1.001 - theo))^2)
          
          loss <- loss + tloss
        }
        
        for(i in 3:3){
          
          theo <- ifelse(h!=0, rho*sqrt(var[1]*var[2])*(h/beta)^nu3*besselK(h/beta,nu3)/(2^(nu3-1)*gamma(nu3)), rho*sqrt(var[1]*var[2]))
          if (weights == 1) 
            tloss <- sum((theo - emp_cov1[,i+3])^2)
          if (weights == 2) 
            tloss <- sum(np * (emp_cov1[,i+3] - theo)^2)
          if (weights == 3) 
            tloss<-sum(((emp_cov1[,i+3] - theo)/(max(emp_cov1[, 3+j]) - theo))^2)
          loss <- loss + tloss
        }
        return(loss)
      }else{
        
        Sigma <- matrix(c(theta[1], theta[2], theta[2], theta[3]), ncol=2)
        
        if(meters == T){
          #w <- emp_cov1[which.max(emp_cov1[,4]),1:2]/1000
          w <- theta[4:5]/1000
          #w2 <- theta[6:7]/1000
          hh <- emp_cov1[,1:2]%*%R/1000
          hh2 <- cbind(hh[,1] - kappa[1], hh[,2] - kappa[2])
        }else{
          #w <- emp_cov1[which.max(emp_cov1[,4]),1:2]
          w <- theta[4:5]
          #w2 <- theta[6:7]
          hh <- emp_cov1[,1:2]%*%R
          hh2 <- cbind(hh[,1] - kappa[1], hh[,2] - kappa[2])
        }
        
        if(min(eigen(Sigma)$val) < 0 | theta[1] < 0 | theta[3] < 0){
          return(Inf)
        }else{
          for(i in 1:2){
            
            Int.func <- function(c,hvec){   
              y.fun  <- function(y) y^(nu[i])*exp(-y)*dmvn(X=hvec[1:2], mu=hvec[3]*w, sigma=(hvec[3]^2*Sigma + beta^2*2*y*diag(2)))
              sapply(c, y.fun)
            }
            lai <- function(xxxx) integrate(Int.func, lower = 0, upper = Inf, hvec = xxxx)$val
            
            theo <- apply(cbind(hh, 1), 1, lai)
            if(weight == 3){
              tloss <- sum(((emp_cov1[,i+3] - 4*pi*beta^2/gamma(nu[i])*theo)/(1.001 - 4*pi*beta^2/gamma(nu[i])*theo))^2)
            }else{
              tloss <- sum(((emp_cov1[,i+3] - 4*pi*beta^2/gamma(nu[i])*theo))^2)
            }
            loss <- loss + tloss
          }
          
          for(i in 3:3){
            
            Int.func <- function(c, hvec){   
              y.fun  <- function(y) y^(nu3)*exp(-y)*dmvn(X = hvec[1:2], mu = hvec[3]*w, sigma = (hvec[3]^2*Sigma + beta^2*2*y*diag(2)))
              sapply(c, y.fun)
            }
            lai <- function(xxxx) integrate(Int.func, lower = 0, upper=Inf, hvec = xxxx)$val
            theo <- apply(cbind(hh2, 1), 1, lai)
            if(weight == 3){
              tloss <- sum(((emp_cov1[,i+3] - sqrt(var[1] * var[2])*rho*4*pi*beta^2/gamma(nu3)*theo)/(max(emp_cov1[,i+3]) - sqrt(var[1] * var[2])*rho*4*pi*beta^2/gamma(nu3)*theo))^2)
            }else{
              tloss <- sum(((emp_cov1[,i+3] - sqrt(var[1] * var[2])*rho*4*pi*beta^2/gamma(nu3)*theo))^2)
            }
            loss <- loss + tloss
          }
          return(loss)
        }
      }
    }
  } else if(model == 'lmc'){
    if(step == 1){
      if(nug_eff == T){
        nug <- theta[7:8]
        alpha <-  matrix(theta[9:12], ncol=2, byrow=T)
      }else{
        nug <- c(0, 0)
        alpha <-  matrix(theta[7:10], ncol=2, byrow=T)
      }
      nu <- theta[1:2]
      beta <- theta[3:4]
      var <- theta[5:6]
      
      if( theta[1] < 0.0001 | theta[2] < 0.0001 | theta[3] < 0.0001 | theta[4] < 0.0001 | theta[5] < 0 | theta[6] < 0 |
          theta[7]< 0 |  theta[7] > 1 | theta[8] < 0 | theta[8] > 1 | theta[9] < 0 |  theta[9] > 1 | theta[10] < 0 | theta[10] > 1){
        return(Inf)
      }else{
        theo_list <- list()
        
        for(i in 1:2){
          theo_list_temp <- ifelse(h != 0, var[i]*(h/beta[i])^nu[i]*besselK(h/beta[i],nu[i])/(2^(nu[i]-1)*gamma(nu[i])),var[i]+nug[i])
          theo_list[[i]] <- theo_list_temp
        }
        
        for(j in 1:2){
          theo <- alpha[j,1]^2*theo_list[[1]] + alpha[j,2]^2*theo_list[[2]]
          if (weights == 1) 
            tloss <- sum((theo - emp_cov1[, 3+j])^2)
          if (weights == 2) 
            tloss <- sum(np*(emp_cov1[, 3+j] - theo)^2)
          if (weights == 3) 
            tloss <- sum(((emp_cov1[, 3+j] - theo)/(1.001 - theo))^2)
          
          loss <- loss + tloss
        }
        
        theo <- alpha[1,1]*alpha[2,1]*theo_list[[1]] + alpha[1,2]*alpha[2,2]*theo_list[[2]]
        if (weights == 1) 
          tloss <- sum((theo - emp_cov1[,6])^2)
        if (weights == 2) 
          tloss <- sum(np*(emp_cov1[,6] - theo)^2)
        if (weights == 3) 
          tloss <- sum(((emp_cov1[,6] - theo)/(max(emp_cov1[, 3+j]) - theo))^2)
        
        loss <- loss + tloss
        
        return(loss)
      }
    }else if (step == 2){
      
      if(nug_eff == T){
        nug <- est_param[6:7]
      }else{
        nug <- c(0, 0)
      }
      nu <- est_param[1:2]
      beta <- est_param[3:4]
      var <- est_param[5:6]
      
      loss <- 0
      
      alpha <-  matrix(est_param[7:10], ncol=2, byrow=T)
      
      w <- matrix(theta[1:4], ncol=2, byrow=T)
      
      if(!rand.vel == T){
        theo_list <- list()
        
        for(i in 1:2){
          
          if(meters == T){
            h <- sqrt((emp_cov1[,1] - emp_cov1[,3]*w[i,1])^2 + (emp_cov1[,2] - emp_cov1[,3]*w[i,2])^2)/1000
          }else{
            h <- sqrt((emp_cov1[,1] - emp_cov1[,3]*w[i,1])^2 + (emp_cov1[,2] - emp_cov1[,3]*w[i,2])^2)
          }
          theo_list_temp <- ifelse(h != 0, var[i]*(h/beta[i])^nu[i]*besselK(h/beta[i], nu[i])/(2^(nu[i]-1)*gamma(nu[i])), var[i] + nug[i])
          theo_list[[i]] <- theo_list_temp
        }
      }else{
        sigma <- list()
        sigma[[1]] <- matrix(c(theta[5:6], theta[6:7]), ncol=2, byrow=T)
        sigma[[2]] <- matrix(c(theta[8:9], theta[9:10]), ncol=2, byrow=T)
        
        if(theta[5] < 0 | theta[7] < 0 | min(eigen(sigma[[1]])$values) <= 0 | theta[8] < 0 | theta[10] < 0 | min(eigen(sigma[[2]])$values) <= 0 ){
          return(Inf)
        }else{
          theo_list <- list()
          
          for(i in 1:2){
            
            if(meters == T){
              h <- sqrt(diag((emp_cov1[,1:2] - matrix(c(emp_cov1[, 3]*w[i,1], emp_cov1[, 3]*w[i,2]), ncol=2))%*%solve(diag(2) + sigma[[i]])%*%t(emp_cov1[,1:2] - matrix(c(emp_cov1[, 3]*w[i,1], emp_cov1[, 3]*w[i,2]), ncol=2))))/1000
            }else{
              h <- sqrt(diag((emp_cov1[,1:2] - matrix(c(emp_cov1[, 3]*w[i,1], emp_cov1[, 3]*w[i,2]), ncol=2))%*%solve(diag(2) + sigma[[i]])%*%t(emp_cov1[,1:2] - matrix(c(emp_cov1[, 3]*w[i,1], emp_cov1[, 3]*w[i,2]), ncol=2))))
            }
            theo_list_temp <- ifelse(h != 0, var[i]*(h/beta[i])^nu[i]*besselK(h/beta[i], nu[i])/(2^(nu[i]-1)*gamma(nu[i])), var[i] + nug[i])
            theo_list[[i]] <- theo_list_temp/sqrt(det(diag(2) + sigma[[i]]))
          }
        }
      }
      
      for(j in 1:2){
        theo <- alpha[j,1]^2*theo_list[[1]] + alpha[j,2]^2*theo_list[[2]]
        if (weights == 1) 
          tloss <- sum((theo - emp_cov1[, 3+j])^2)
        if (weights == 2) 
          tloss <- sum(np* (emp_cov1[, 3+j] - theo)^2)
        if (weights == 3) 
          tloss <- sum(((emp_cov1[, 3+j] - theo)/(1.001 - theo))^2)
        
        loss <- loss+tloss
      }
      
      theo <- alpha[1,1]*alpha[2,1]*theo_list[[1]] + alpha[1,2]*alpha[2,2]*theo_list[[2]]
      if (weights == 1) 
        tloss <- sum((theo - emp_cov1[, 6])^2)
      if (weights == 2) 
        tloss <- sum(np* (emp_cov1[, 6] - theo)^2)
      if (weights == 3) 
        tloss <- sum(((emp_cov1[, 6] - theo)/(max(emp_cov1[, 3+j]) - theo))^2)
      
      loss <- loss+tloss
      
      return(loss)
      
    }
  }else if(mod == 'gneitingmatern'){
    
    if(nug_eff == T){
      nug <- est_param[6:7]
      rho <- est_param[8]
    }else{
      nug <- c(0, 0)
      rho <- est_param[6]
    }
    nu <- est_param[1:2]
    beta <- est_param[3]
    var <- est_param[4:5]
    
    alpha <- theta[1]
    b <- theta[2]
    
    nu1 <- nu[1]
    nu2 <- nu[2]
    nu3 <- (nu1 + nu2)/2
    
    loss<- 0
    
    new_h <- h/(alpha + 1)^(b/2)
    
    if( theta[1] < 0.0001 | theta[2] < 0.0001 | theta[2] > 1){
      return(Inf)
    }else{
      for(i in 1:2){
        theo <- ifelse(h != 0 & emp_cov1[, 3] != 0, var[i]*(new_h/beta)^nu[i] * besselK(new_h/beta, nu[i])/(2^(nu[i] - 1)*gamma(nu[i]))/(alpha + 1), (var[i] + nug[i])/(alpha + 1))
        if (weights == 1) 
          tloss <- sum((theo - emp_cov1[,i+3])^2)
        if (weights == 2) 
          tloss <- sum(np* (emp_cov1[,i+3] - theo)^2)
        if (weights == 3) 
          tloss <- sum(((emp_cov1[, i + 3] - theo)/(1.001 - theo))^2)
        
        loss <- loss + tloss
      }
      
      for(i in 3:3){
        
        theo <- ifelse(h !=0 & emp_cov1[, 3] != 0, (new_h/beta)^nu3 * besselK(new_h/beta, nu3)/(2^(nu3 - 1)*gamma(nu3))*sqrt(var[1] * var[2])*rho/((alpha + 1)), sqrt(var[1] * var[2])*rho/((alpha + 1)))
        if (weights == 1) 
          tloss <- sum((theo-emp_cov1[, i+3])^2)
        if (weights == 2) 
          tloss <- sum(np* (emp_cov1[, i+3] - theo)^2)
        if (weights == 3) 
          tloss <- sum(((emp_cov1[, i+3] - theo)/(1.001 - theo))^2)
        
        loss <- loss+tloss
      }
      return(loss)
    }
  }
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
