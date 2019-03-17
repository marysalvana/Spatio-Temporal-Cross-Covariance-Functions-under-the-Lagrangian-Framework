#---------STATIONARY---------#

fitting <- function(coordinates, obs1, obs2 = NULL){
  # FORMAT OF INPUT VALUES
  # coordinates: two-column matrix with longitude in column 1 and latitude in column 2
  # obs1: T-column matrix of observations where a row represents observations of one location and the columns are hourly/monthly/yearly observations
  # obs2: another covariate
  
  data_matrix <- data_format_into_matrix(data1 = obs1, data2 = obs2, temporal_replicates = ncol(obs1), simulated = F)
  conso_cor <- empirical_st_cov(data1 = t(data_matrix[insample_loc_index, 1:(ncol(data_matrix)/2)]), data2 = t(data_matrix[insample_loc_index, (ncol(data_matrix)/2 + 1):ncol(data_matrix)]), locations = coordinates[insample_loc_index,], max_time_lag = 5, simulated = F)
  binned <- empirical_covariance_dataframe(data_cov = conso_cor)
  
  #Define initial parameters
  
  R <- matrix(c(-1.1339125*cos(116.9811232), -1.1339125*sin(116.9811232), 0.8617141*sin(116.9811232),
                -0.8617141*cos(116.9811232)), ncol=2, byrow=T)
  
  theta_init <- c(3,4, 794.8, 1, 1, max(binned[which(binned[,3] == 0), 6]) - 0.001)
  emp_vel <- which.max(binned[which(binned[,3] == 1), 4])
  #kappa_ind <- which.max(binned[which(binned[,3] == 0), 6])
  #kappa <- binned[which(binned[,3] == 0)[kappa_ind], 1:2] %*% R
  kappa <- matrix(c(704400, 205700), ncol = 2)
  #w_init <- binned[which(binned[,3] == 1)[emp_vel], 1:2] %*% R
  w_init <- c(-2495995, -1783349) %*% R
  decay <- c(0.4, 0.1)
  t_init <- c(0.99, 0.99)
  other_init <- c(0.5, 2000, 0.5, 0.5, 6)

  #w_init=c(-756163.3, 579430.6)
  #w_var <- c(3905504484.2, 1093290456.9, 1545066174.0)
  #-1038390.4     278837.9 4455873832.1 1223028720.4 2222390454.9
  
  w_init <- c(-393048.2, 475609.2)
  w_var <- c(2014856.1, 121333.7, 2496301.5)
  #-756724.6  579775.5 3901641.7 1097350.8 1547412.4
  #-1039768.7   279262.1  4535528.4  1181704.7  2149263.4
  
  lmc_init <- c(1, 2, 1200.6314, 1591.1927, 1.21802, 1.033560, 0.8621720, 0.5367004, 0.0000002056025, 0.9935446)
  lmc_w_init <- c(-897386.7980, -438772.3458, -373454.1753, 167183.9760, 6258.984, -032.18460, 1654.979, 1380.273, -0070.52379, 0112.9060)
  
  M1 <- fit_model(init = theta_init, wind_init = c(w_init, w_var), mod = 'matern', randvel = T, weight = 3, empcov_spatial = binned[which(binned[,3] == 0),], empcov_st = binned[which(binned[,3] == 1),], nug_eff = F, meters = T, num_iter = 0, aniso = T, rotation_matrix = R, asym = kappa)
  M2 <- fit_model(init = theta_init, wind_init = c(w_init, decay), mod = 'matern', randvel = F, weight = 3, empcov_spatial = binned[which(binned[,3] == 0),], empcov_st = binned[which(binned[,3] == 1),], nug_eff = F, meters = T, num_iter = 10, aniso = T, rotation_matrix = R, asym = kappa)
  M3 <- fit_model(init = lmc_init, wind_init = lmc_w_init, mod = 'lmc', randvel = T, weight = 3, empcov_spatial = binned[which(binned[,3] == 0),], empcov_st = binned[which(binned[,3] == 1),], nug_eff = F, meters = T, num_iter = 0, aniso = T, rotation_matrix = R, asym = kappa)
  M4 <- fit_model(init = c(theta_init[-length(theta_init)], other_init, theta_init[length(theta_init)]), temp_init = c(t_init, w_init), mod = 'gneitingmatern_lagrangian', randvel = F, weight = 3, empcov_spatial = binned[which(binned[,3] == 0),], empcov_st = binned[which(binned[,3] == 1),], nug_eff = F, meters = T, num_iter = 10, aniso = T, rotation_matrix = R)
  M5 <- fit_model(init = theta_init, temp_init = t_init, mod = 'gneitingmatern', randvel = F, weight = 3, empcov_spatial = binned[which(binned[,3] == 0),], empcov_st = binned[which(binned[,3] == 1),], nug_eff = F, meters = T, num_iter = 0, aniso = T, rotation_matrix = R)
  
  hlag <- sqrt(binned[which(binned[,3] == 0), 1]^2 + binned[which(binned[,3] == 0), 2]^2)
  #display plots
  par(mfrow = c(1,3))
  plot(hlag/1000, binned[which(binned[,3] == 0), 4], pch = 3, ylab='', col = 1, xlab = 'Spatial Lag (km)', main='', col.main = "#4EC1DE")
  plot(hlag/1000, binned[which(binned[,3] == 0), 5], pch = 3, ylab='', col = 1, xlab = 'Spatial Lag (km)', main='', col.main = "#4EC1DE")
  plot(hlag/1000, binned[which(binned[,3] == 0), 6], pch = 3, ylab='', col = 1, xlab = 'Spatial Lag (km)', main='', col.main = "#4EC1DE")
}

fit_model <- function(init = NULL, wind_init = NULL, temp_init = NULL, mod, randvel, weight, empcov_spatial = NULL, empcov_st, nug_eff = F, meters = T, num_iter, 
                      est_param.temp = NULL, aniso = F, rotation_matrix = NULL, asym = NULL){
  
  # num_iter : number of loops to run optim
  
  if(mod == 'matern' | mod == 'gneitingmatern'){
    fit1.mod <- optim(par = init[-length(init)], wls, emp_cov1 = empcov_spatial, nug_eff = F, meters = T, weights = weight, step = 1, ani = aniso, rotmat = rotation_matrix, asymmetry = asym, model = mod, control=list(maxit = 10000, parscale = init[-length(init)], trace = 5))
    fit2.mod <- optim(par = init[length(init)], wls, emp_cov1 = empcov_spatial, nug_eff = F, meters = T, weights = weight, step = 2, est_param = fit1.mod$par, ani = aniso, rotmat = rotation_matrix, asymmetry = asym, model = mod, method='SANN', control = list(maxit = 3000, parscale = init[length(init)], trace = 5))
    
    if(mod == 'matern'){
      if(randvel == T){
        
        fit3.mod <- optim(par = wind_init, wls, emp_cov1 = empcov_st, nug_eff = F, meters = T, weights = weight, step = 3, est_param = c(fit1.mod$par, fit2.mod$par), ani = aniso, rotmat = rotation_matrix, asymmetry = asym, rand.vel = randvel, model = mod, control = list(maxit = 10000, parscale = wind_init, trace = 5))
        
        fit3.mod <- optim(par = wind_init, wls, emp_cov1 = empcov_st[whichpart(empcov_st[, 4], 10),], nug_eff = F, meters = T, weights = weight, step = 3, est_param = c(fit1.mod$par, fit2.mod$par), ani = aniso, rotmat = rotation_matrix, asymmetry = asym, rand.vel = randvel, model = mod, control = list(maxit = 10000, parscale = wind_init, trace = 5))
        
        new_wind_init <- fit3.mod$par
        fit3.mod <- optim(par = new_wind_init, wls, emp_cov1 = empcov_st[whichpart(empcov_st[, 4], 3000),], nug_eff = F, meters = T, weights = weight, step = 3, est_param = c(fit1.mod$par, fit2.mod$par), ani = aniso, rotmat = rotation_matrix, asymmetry = asym, rand.vel = randvel, model = mod, control = list(maxit = 10000, parscale = wind_init, trace = 5))
        
        lst <- list(parameters = c(fit1.mod$par, fit2.mod$par, fit3.mod$par), fn_value = fit1.mod$value + fit2.mod$value + fit3.mod$value)
        
        return(lst)
      }else if(randvel == F){
        
        fit3.mod <- optim(par = wind_init, wls, emp_cov1 = empcov_st, nug_eff = F, meters = T, weights = weight, step = 3, est_param = c(fit1.mod$par, fit2.mod$par), ani = aniso, rotmat = rotation_matrix, asymmetry = asym, rand.vel = randvel, model = mod, control = list(maxit = 10000, parscale = wind_init, trace = 5))
        
        if(num_iter > 0){
          for(iter3 in 1:num_iter){
            new_wind_init <- fit3.mod$par
            fit3.mod <- optim(par = new_wind_init, wls, emp_cov1 = empcov_st, nug_eff = F, meters = T, weights = weight, step = 3, est_param = c(fit1.mod$par, fit2.mod$par), ani = aniso, rotmat = rotation_matrix, asymmetry = asym, rand.vel = randvel, model = mod, control = list(maxit = 10000, parscale = new_wind_init, trace = 5))
          }
        }
        
        lst <- list(parameters = c(fit1.mod$par, fit2.mod$par, fit3.mod$par), fn_value = fit1.mod$value + fit2.mod$value + fit3.mod$value)
        
        return(lst)
        
      }
    }else if(mod == 'gneitingmatern'){
      fit3.mod <- optim(par = temp_init, wls, emp_cov1 = empcov_st, nug_eff = F, meters = T, weights = weight, step = 3, est_param = c(fit1.mod$par, fit2.mod$par), ani = aniso, rotmat = rotation_matrix, asymmetry = asym, model = mod, control = list(maxit = 10000, parscale = temp_init, trace = 5))
      
      if(num_iter > 0){
        for(iter3 in 1:num_iter){
          new_temp_init <- fit3.mod$par
          fit3.mod <- optim(par = new_temp_init, wls, emp_cov1 = empcov_st, nug_eff = F, meters = T, weights = weight, step = 3, est_param = c(fit1.mod$par, fit2.mod$par), ani = aniso, rotmat = rotation_matrix, asymmetry = asym, rand.vel = randvel, model = mod, control = list(maxit = 10000, parscale = new_wind_init, trace = 5))
        }
      }
      
      lst <- list(parameters = c(fit1.mod$par, fit2.mod$par, fit3.mod$par), fn_value = fit1.mod$value + fit2.mod$value + fit3.mod$value)
      
      return(lst)
    }
  }else if(mod == 'lmc'){
    if(randvel == T){
      fit1.mod <- optim(par = init, wls, emp_cov1 = empcov_spatial, nug_eff = F, meters = T, weights = weight, step = 1, ani = aniso, rotmat = rotation_matrix, asymmetry = asym, model = mod, control=list(maxit = 10000, parscale = init, trace = 5))
      
      for(aa in 1:2){
        new_init <- fit1.mod$par
        fit1.mod <- optim(par = new_init, wls, emp_cov1 = empcov_spatial, nug_eff = F, meters = T, weights = weight, step = 1, ani = aniso, rotmat = rotation_matrix, asymmetry = asym, model = mod, control=list(maxit = 10000, parscale = init, trace = 5))
        
      }
           
      fit2.mod <- optim(par = wind_init, wls, emp_cov1 = empcov_st, nug_eff = F, meters = T, weights = weight, step = 2, est_param = fit1.mod$par, ani = aniso, rotmat = rotation_matrix, asymmetry = asym, rand.vel = randvel, model = mod, control = list(maxit = 10000,parscale = wind_init, trace = 5))
      
      lst <- list(parameters = c(fit1.mod$par, fit2.mod$par), fn_value = fit1.mod$value + fit2.mod$value)
      
      return(lst)
    }else if(randvel == F){
      fit1.mod <- optim(par = init, wls, emp_cov1 = empcov_spatial, nug_eff = F, meters = T, weights = weight, step = 1, ani = aniso, rotmat = rotation_matrix, asymmetry = asym, model = mod, control=list(maxit = 10000, parscale = init, trace = 5))
      fit2.mod <- optim(par = wind_init, wls, emp_cov1 = empcov_st, nug_eff = F, meters = T, weights = weight, step = 2, est_param = fit1.mod$par, ani = aniso, rotmat = rotation_matrix, asymmetry = asym, rand.vel = randvel, model = mod, control = list(maxit = 10000,parscale = wind_init, trace = 5))
      
      lst <- list(parameters = c(fit1.mod$par, fit2.mod$par), fn_value = fit1.mod$value + fit2.mod$value)
      
      return(lst)
    }
  }else if(mod == 'gneitingmatern_lagrangian'){
    fit1.mod <- optim(par = init[-length(init)], wls, emp_cov1 = empcov_spatial, nug_eff = F, meters = T, weights = weight, step = 1, ani = aniso, rotmat = rotation_matrix, asymmetry = asym, model = mod, control=list(maxit = 10000, parscale = init[-length(init)], trace = 5))
    for(aa in 1:10){
      new_init <- fit1.mod$par
      fit1.mod <- optim(par = new_init, wls, emp_cov1 = empcov_spatial, nug_eff = F, meters = T, weights = weight, step = 1, ani = aniso, rotmat = rotation_matrix, asymmetry = asym, model = mod, control=list(maxit = 10000, parscale = init[-length(init)], trace = 5))
    }
    fit2.mod <- optim(par = init[length(init)], wls, emp_cov1 = empcov_spatial, nug_eff = F, meters = T, weights = weight, step = 2, est_param = fit1.mod$par, ani = aniso, rotmat = rotation_matrix, asymmetry = asym, model = mod, method = 'SANN', control=list(maxit = 3000, parscale = init[length(init)], trace = 5))
    fit3.mod <- optim(par = temp_init, wls, emp_cov1 = empcov_st, nug_eff = F, meters = T, weights = weight, step = 3, est_param = c(fit1.mod$par, fit2.mod$par), ani = aniso, rotmat = rotation_matrix, asymmetry = asym, rand.vel = randvel, model = mod, control = list(maxit = 10000, parscale = temp_init, trace = 5))
    
    lst <- list(parameters = c(fit1.mod$par, fit2.mod$par, fit3.mod$par), fn_value = fit1.mod$value + fit2.mod$value + fit3.mod$value)
    
    return(lst)
  }
}

kriging <- function(params, model, max_time_lag = 1, q = 2){
  
  obs <- matrix(, nrow = nrow(loc)*q*(max_time_lag + 1) - length(to_remove), ncol = (ncol(obs1) - 1))
  valid_pts <- matrix(, nrow = length(to_remove), ncol = (ncol(obs1) - 1))
  
  for(dd in 1:(ncol(obs1) - 1)){
    
    orig_var1 <- orig_var2 <- matrix(, ncol = (max_time_lag + 1), nrow = nrow(loc))
    
    for(aa in 1:time){
      orig_var1[,aa] <- obs1[, aa + dd - 1]
      orig_var2[,aa] <- obs2[, aa + dd - 1]
    }
    obs[, dd] <- c(c(orig_var1),c(orig_var2))[-to_remove]
    valid_pts[, dd] <- matrix(c(c(orig_var1), c(orig_var2))[to_remove], ncol=1, byrow=F)
  }
  
  if(model == 'M1'){
    sim.cov <- simulate_model(mod = 'matern', randvel = T, theta = params[1:6], wind = params[7:8], wind_var = params[9:11], maxtimelag = max_time_lag, p = q, locations = loc, meters = T, nugeff = nugget, kaps = c(704400, 205700))
  }
  
  if(model == 'M2'){
    th <- params[9:10]
    sim.cov <- simulate_model(mod = 'matern', randvel = F, theta = params[1:6], wind = params[7:8], maxtimelag = max_time_lag, p = q, locations = loc, meters = T, nugeff = nugget, kaps = c(704400, 205700))
    
    sim.cov[1:100,101:200] <- exp(-th[1])*sim.cov[1:100,101:200]
    sim.cov[101:200,1:100] <- exp(-th[1])*sim.cov[101:200,1:100] 
    sim.cov[201:300,301:400] <- exp(-th[2])*sim.cov[201:300,301:400]
    sim.cov[301:400,201:300] <- exp(-th[2])*sim.cov[301:400,201:300]
    sim.cov[1:100,301:400] <- exp(-(th[1] + th[2])/2)*sim.cov[1:100,301:400]
    sim.cov[301:400,1:100] <- exp(-(th[1] + th[2])/2)*sim.cov[301:400,1:100]
    sim.cov[101:200,201:300] <- exp(-(th[1] + th[2])/2)*sim.cov[101:200,201:300]
    sim.cov[201:300,101:200] <- exp(-(th[1] + th[2])/2)*sim.cov[201:300,101:200]
  }
  
  if(model == 'M3'){
    sim.cov <- simulate_model(mod = 'lmc', randvel = T, theta = params[1:10], wind = params[11:14], wind_var = list(matrix(c(params[15:16], params[16:17]), ncol=2), matrix(c(params[18:19], params[19:20]), ncol=2)), maxtimelag = max_time_lag, p = q, locations = loc, meters = T, nugeff = nugget)
  }
  
  if(model == 'M4'){
    sim.cov_old <- simulate_model(mod = 'M5', theta = c(params[1:5], params[11:13]), maxtimelag = 1, p = q, locations = loc, meters = T, nugeff = nugget)
    sim.cov <- simulate_model(mod = model, theta = c(params[8:10], params[7]), wind = params[14:15], maxtimelag = max_time_lag, p = q, locations = loc, meters = T, nugeff = nugget)
    full_cov1 <- (1-params[6])*sim.cov_old + params[6]*sim.cov
  }
  
  if(model == 'M5'){
    sim.cov <- simulate_model(mod = model, theta = params, maxtimelag = 1, p = q, locations = loc, meters = T, nugeff = nugget)
  }
  
  full_cov1 <- sim.cov
  given_cov1 <- full_cov1[-to_remove, -to_remove]
  unknown_cov1 <- full_cov1[to_remove, -to_remove]
  var_pred1 <- diag(full_cov1[to_remove, to_remove] - unknown_cov1 %*% solve(given_cov1) %*% t(full_cov1[to_remove, -to_remove]))
  
  test1 <- unknown_cov1 %*%chol2inv(chol(given_cov1))%*% obs
  
  diff1 <- test1 - valid_pts
  
  rmse <- sqrt(mean(diff1^2))
  
  mae <- mean(abs(diff1))
  
  z <- diff1/var_pred1 ## center and scale
  crps<- var_pred1 * (z*(2*pnorm(z, 0, 1) - 1) + 2*dnorm(z, 0, 1) - 1/sqrt(pi))
  crps_mat <- mean(crps)
  
  if(model == 'M1' | model == 'M2'){
    aic <- q*length(to_remove)*log(xx) + 2*(length(params) + 4 + 1)
  }else{
    aic <- q*length(to_remove)*log(mean(diff1^2)) + 2*(length(params) + 2 + 1)
  }

  mean(rmse)
  mean(crps_mat)
  mean(aic)
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

empirical_covariance_dataframe <- function(data_cov){
  
  empirical_var1 <- data_cov[[1]]
  
  binned.1 <- empirical_var1 %>% group_by(xlag, ylag) %>% summarize(avg1=mean(var1_cor))
  binned.2 <- empirical_var1 %>% group_by(xlag, ylag) %>% summarize(avg1=mean(var2_cor))
  binned.3 <- empirical_var1 %>% group_by(xlag, ylag) %>% summarize(avg1=mean(cross_cor))
  
  binned_orig <- cbind(binned.1$xlag, binned.1$ylag, rep(0,nrow(binned.1)), binned.1$avg1,
                       binned.2$avg1, binned.3$avg1)
  
  if(length(data_cov) > 1){
    for (i in 2:length(data_cov)){
      empirical_var1 <- data_cov[[i]]
      
      binned.1 <- empirical_var1 %>% group_by(xlag, ylag) %>% summarize(avg1=mean(var1_cor))
      binned.2 <- empirical_var1 %>% group_by(xlag, ylag) %>% summarize(avg1=mean(var2_cor))
      binned.3 <- empirical_var1 %>% group_by(xlag, ylag) %>% summarize(avg1=mean(cross_cor))
      
      binned_orig <- rbind(binned_orig,cbind(binned.1$xlag,binned.1$ylag,rep(i-1,nrow(binned.1)),binned.1$avg1,
                                             binned.2$avg1,binned.3$avg1))
    }
  }
  colnames(binned_orig) <- c('xlag', 'ylag', 'tlag', 'var1_cor', 'var2_cor', 'cross_cor')
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

simulate_model <- function(mod, theta, wind = NULL, wind_var = NULL, maxtimelag, p = 2, locations, meters = T, nugeff, randvel = NULL, kaps = NULL){
  if(!is.null(randvel)){
    if(mod == 'matern' & randvel == T){
      cov.mod <- matern_random_cov(theta, wind, wind_var, max_time_lag = maxtimelag, q = p, new_locations = locations, nug_eff = nugeff, kap = kaps)
    }else if(mod == 'matern' & randvel == F){
      cov.mod <- matern_cov(theta, wind, max_time_lag = maxtimelag, q = p, new_locations = locations, nug_eff = nugeff, kap = kaps)
    }else if(mod == 'lmc' & randvel == T){
      cov.mod <- lmc_random_cov(theta, wind, wind_var, max_time_lag = maxtimelag, q = p, new_locations = locations, nug_eff = nugeff )
    }else if(mod == 'lmc' & randvel == F){
      cov.mod <- lmc_cov(theta, wind, max_time_lag = maxtimelag, q = p, new_locations = locations, nug_eff = nugeff )
    }
  }else{
    if(mod == 'M5'){
      cov.mod <- matern_allard(theta, max_time_lag = maxtimelag, q = p, new_locations = locations, nug_eff = nugeff )
    }
  }
  return(cov.mod)
}

wls <- function(theta, emp_cov1, weights, nug_eff, step, est_param = NULL, meters, ani, rotmat, asymmetry, model, rand.vel) {
  
  # nug_eff: T means we also need to estimate nugget effect
  
  if(!ani == T){
    R <- diag(2)
  }else{
    R <- rotmat
  }
  
  hh <- emp_cov1[,1:2] %*% R
  
  if(!is.null(asymmetry)){
    hh2 <- cbind(hh[, 1] - asymmetry[1,1], hh[, 2] - asymmetry[1,2])
  }else{
    hh2 <- hh
  }
  
  if(meters == T){
    h <- sqrt(hh[, 1]^2 + hh[, 2]^2)/1000
    h2 <- sqrt(hh2[, 1]^2 + hh2[, 2]^2)/1000
  }else{
    h <- sqrt(hh[, 1]^2 + hh[, 2]^2)
    h2 <- sqrt(hh2[, 1]^2 + hh2[, 2]^2)
  }
  
  np <- ifelse(h != 0, 1/h, 1)
  
  if( model == 'matern' | model == 'gneitingmatern'){
    if(step == 1){
      
      loss<- 0
      
      if(nug_eff == T){
        nug <- theta[6:7]
      }else{
        nug <- c(0, 0)
      }
      nu <- theta[1:2]
      beta <- theta[3]
      var <- theta[4:5]
      
      sim.cov <- simulate_model(mod = 'M5', theta = c(theta, theta_init[length(theta_init)], t_init), maxtimelag = 0, p = q, locations = loc, meters = T, nugeff = nug_eff)
      
      #sim.cov <- matern_allard(theta = c(theta, theta_init[length(theta_init)], temp_init), max_time_lag = 0, q=2, new_locations = locations, meters = T, nug_eff = F)
      
      full_cov1 <- sim.cov
      given_cov1 <- full_cov1[-to_remove_old, -to_remove_old]
      unknown_cov1 <- full_cov1[to_remove_old, -to_remove_old]
      var_pred1 <- tryCatch(diag(full_cov1[to_remove_old, to_remove_old] - unknown_cov1 %*% solve(given_cov1) %*% t(full_cov1[to_remove_old, -to_remove_old])), error = function(e) NULL)
      
      
      if( theta[1] < 0.0001 | theta[2] < 0.0001 | theta[3] < 0.0001 | theta[4] < 0.0001 | theta[5] < 0.0001 | length(var_pred1) == 0 ){
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
      
      loss<- 0
      
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
      
      sim.cov <- simulate_model(mod = 'M5', theta = c(est_param, theta, t_init), maxtimelag = 0, p = q, locations = loc, meters = T, nugeff = nug_eff)
      full_cov1 <- sim.cov
      given_cov1 <- full_cov1[-to_remove_old, -to_remove_old]
      unknown_cov1 <- full_cov1[to_remove_old, -to_remove_old]
      var_pred1 <- tryCatch(diag(full_cov1[to_remove_old, to_remove_old] - unknown_cov1 %*% solve(given_cov1) %*% t(full_cov1[to_remove_old, -to_remove_old])), error = function(e) NULL)
      
      if( abs(rho) > 1 | length(var_pred1) == 0){
        return(Inf)
      }else{
        
        theo <- ifelse(h2 != 0, rho*sqrt(var[1]*var[2])*(h2/beta)^nu3*besselK(h2/beta, nu3)/(2^(nu3-1)*gamma(nu3)), rho*sqrt(var[1]*var[2]))
        if (weights == 1) 
          tloss <- sum((theo - emp_cov1[,6])^2)
        if (weights == 2) 
          tloss <- sum(np* (emp_cov1[,6] - theo)^2)
        if (weights == 3) 
          tloss <- sum(((emp_cov1[,6] - theo)/(1.001 - theo))^2)
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
      nu3 <- (nu[1] + nu[2])/2
      
      if(model == 'matern'){
        
        loss<- 0
        
        if(!rand.vel == T){
          w <- theta[1:2]
          c <- theta[3:4]
          
          c3 <- (c[1] + c[2])/2
          
          if(meters == T){
            h <- sqrt((hh[,1] - emp_cov1[,3]*w[1])^2 + (hh[,2] - emp_cov1[,3]*w[2])^2)/1000
            h2 <- sqrt((hh2[,1] - emp_cov1[,3]*w[1])^2 + (hh2[,2] - emp_cov1[,3]*w[2])^2)/1000
          }else{
            h <- sqrt((hh[,1] - emp_cov1[,3]*w[1])^2 + (hh[,2] - emp_cov1[,3]*w[2])^2)
            h2 <- sqrt((hh2[,1] - emp_cov1[,3]*w[1])^2 + (hh2[,2] - emp_cov1[,3]*w[2])^2)
          }
          
          np <- ifelse(h != 0, 1/h, 1)
          if(c[1] <= 0 | c[2] <= 0 ){
            return(Inf)
          }else{
            for(i in 1:2){
              
              theo <- ifelse(h!=0, var[i]*(h/beta)^nu[i]*besselK(h/beta,nu[i])/(2^(nu[i]-1)*gamma(nu[i])), var[i]+nug[i])
              if (weights == 1) 
                tloss <- sum((exp(-c[i])*theo - emp_cov1[,i+3])^2)
              if (weights == 2) 
                tloss <- sum(np* (emp_cov1[,i+3] - exp(-c[i])*theo)^2)
              if (weights == 3) 
                tloss <- sum(((emp_cov1[,i+3] - exp(-c[i])*theo)/(1.001 - exp(-c[i])*theo))^2)
              
              loss <- loss + tloss
            }
            
            for(i in 3:3){
              
              theo <- ifelse(h2 != 0, rho*sqrt(var[1]*var[2])*(h2/beta)^nu3*besselK(h2/beta, nu3)/(2^(nu3-1)*gamma(nu3)), rho*sqrt(var[1]*var[2]))
              if (weights == 1) 
                tloss <- sum((theo - exp(-c3)*emp_cov1[,i+3])^2)
              if (weights == 2) 
                tloss <- sum(np * (emp_cov1[,i+3] - exp(-c3)*theo)^2)
              if (weights == 3) 
                tloss<-sum(((emp_cov1[,i+3] - theo)/(1.001 - exp(-c3)*theo))^2)
              loss <- loss + tloss
            }
            return(loss)
          }
        }else{
          
          beta <- est_param[3]/1000
          
          if(meters == T){
            w <- theta[1:2]/1000000
            Sigma <- matrix(c(theta[3:4], theta[4:5]), ncol=2)/1000000
            h <- hh/1000000
            h2 <- hh2/1000000
          }else{
            w <- theta[1:2]
            Sigma <- matrix(c(theta[3:4], theta[4:5]), ncol=2)
            h <- hh
            h2 <- hh2
          }
          
          if(min(eigen(Sigma)$val) < 0 | theta[3] < 0 | theta[5] < 0){
            return(Inf)
          }else{
            for(i in 1:2){
              
              Int.func <- function(c,hvec){   
                y.fun  <- function(y) y^(nu[i])*exp(-y)*dmvn(X = hvec[1:2], mu = hvec[3]*w, sigma = (hvec[3]^2*Sigma + beta^2*2*y*diag(2)))
                sapply(c, y.fun)
              }
              lai <- function(xxxx) integrate(Int.func, lower = 0, upper = Inf, hvec = xxxx)$val
              
              theo <- apply(cbind(h, emp_cov1[, 3]), 1, lai)
              if(weight == 3){
                tloss <- sum(((emp_cov1[,i+3] - var[i]*4*pi*beta^2/gamma(nu[i])*theo)/(1.001 - var[i]*4*pi*beta^2/gamma(nu[i])*theo))^2)
              }else{
                tloss <- sum(((emp_cov1[,i+3] - var[i]*4*pi*beta^2/gamma(nu[i])*theo))^2)
              }
              loss <- loss + tloss
            }
            
            for(i in 3:3){
              
              Int.func <- function(c, hvec){   
                y.fun  <- function(y) y^(nu3)*exp(-y)*dmvn(X = hvec[1:2], mu = hvec[3]*w, sigma = (hvec[3]^2*Sigma + beta^2*2*y*diag(2)))
                sapply(c, y.fun)
              }
              lai <- function(xxxx) integrate(Int.func, lower = 0, upper=Inf, hvec = xxxx)$val
              theo <- apply(cbind(h2, emp_cov1[, 3]), 1, lai)
              if(weight == 3){
                tloss <- sum(((emp_cov1[,i+3] - sqrt(var[1] * var[2])*rho*4*pi*beta^2/gamma(nu3)*theo)/(1.001 - sqrt(var[1] * var[2])*rho*4*pi*beta^2/gamma(nu3)*theo))^2)
              }else{
                tloss <- sum(((emp_cov1[,i+3] - sqrt(var[1] * var[2])*rho*4*pi*beta^2/gamma(nu3)*theo))^2)
              }
              loss <- loss + tloss
            }
            return(loss)
          }
        }
        
      }else if(model == 'gneitingmatern'){
        
        alpha <- theta[1]
        b <- theta[2]
        
        loss<- 0
        
        new_h <- h/(alpha + 1)^(b/2)
        
        sim.cov <- simulate_model(mod = 'M5', theta = c(est_param, theta), maxtimelag = max_time_lag, p = q, locations = loc, meters = T, nugeff = nug_eff)
        full_cov1 <- sim.cov
        given_cov1 <- full_cov1[-to_remove, -to_remove]
        unknown_cov1 <- full_cov1[to_remove, -to_remove]
        var_pred1 <- tryCatch(diag(full_cov1[to_remove, to_remove] - unknown_cov1 %*% solve(given_cov1) %*% t(full_cov1[to_remove, -to_remove])), error = function(e) NULL)
        
        
        if( theta[1] < 0.000001 | theta[2] < 0.000001 | theta[2] > 1 | length(var_pred1) == 0){
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
            
            theo <- ifelse(h != 0 & emp_cov1[, 3] != 0, (new_h/beta)^nu3 * besselK(new_h/beta, nu3)/(2^(nu3 - 1)*gamma(nu3))*sqrt(var[1] * var[2])*rho/((alpha + 1)), sqrt(var[1] * var[2])*rho/((alpha + 1)))
            if (weights == 1) 
              tloss <- sum((theo - emp_cov1[, i+3])^2)
            if (weights == 2) 
              tloss <- sum(np* (emp_cov1[, i+3] - theo)^2)
            if (weights == 3) 
              tloss <- sum(((emp_cov1[, i+3] - theo)/(1.001 - theo))^2)
            
            loss <- loss + tloss
          }
          return(loss)
        }
      }
    }
  } else if(model == 'lmc'){
    loss <- 0
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
      
      sim.cov <- simulate_model(mod = 'lmc', randvel = T, theta = theta[1:10], wind = theta[1:4], wind_var = list(matrix(c(theta[5:6], theta[6:7]), ncol=2), matrix(c(theta[8:9], theta[9:10]), ncol=2)), maxtimelag = 0, p = q, locations = loc, meters = T, nugeff = nug_eff)
      
      full_cov1 <- sim.cov
      given_cov1 <- full_cov1[-to_remove_old, -to_remove_old]
      unknown_cov1 <- full_cov1[to_remove_old, -to_remove_old]
      var_pred1 <- tryCatch(diag(full_cov1[to_remove_old, to_remove_old] - unknown_cov1 %*% solve(given_cov1) %*% t(full_cov1[to_remove_old, -to_remove_old])), error = function(e) NULL)
      
      
      if( theta[1] < 0.0001 | theta[1] > 4 | theta[2] < 0.0001 | theta[2] < 0.0001 | theta[3] < 0.0001 | theta[4] < 0.0001 | theta[5] < 0.99 | theta[6] < 0.99 |
          theta[7]< 0 |  theta[7] > 1 | theta[8] < 0 | theta[8] > 1 | theta[9] < 0 |  theta[9] > 1 | theta[10] < 0 | theta[10] > 1 | length(var_pred1) == 0){
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
          tloss <- sum(((emp_cov1[,6] - theo)/(1.001 - theo))^2)
        
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
      
      w <- matrix(theta[1:4], ncol=2, byrow=T)/1000
      
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
        sigma[[1]] <- matrix(c(theta[5:6], theta[6:7]), ncol=2, byrow=T)/1000
        sigma[[2]] <- matrix(c(theta[8:9], theta[9:10]), ncol=2, byrow=T)/1000
        
        sim.cov <- simulate_model(mod = 'lmc', randvel = T, theta = est_param[1:10], wind = theta[1:4], wind_var = list(matrix(c(theta[5:6], theta[6:7]), ncol=2), matrix(c(theta[8:9], theta[9:10]), ncol=2)), maxtimelag = max_time_lag, p = q, locations = loc, meters = T, nugeff = nug_eff)
        
        full_cov1 <- sim.cov
        given_cov1 <- full_cov1[-to_remove, -to_remove]
        unknown_cov1 <- full_cov1[to_remove, -to_remove]
        var_pred1 <- tryCatch(diag(full_cov1[to_remove, to_remove] - unknown_cov1 %*% solve(given_cov1) %*% t(full_cov1[to_remove, -to_remove])), error = function(e) NULL)
        
        if(theta[5] < 0 | theta[7] < 0 | min(eigen(sigma[[1]])$values) <= 0 | theta[8] < 0 | theta[10] < 0 | min(eigen(sigma[[2]])$values) <= 0 |
           length(var_pred1) == 0){
          return(Inf)
        }else{
          theo_list <- list()
          
          for(i in 1:2){
            
            if(meters == T){
              hh1 <- hh/1000
              h <- sqrt(diag((hh1 - matrix(c(emp_cov1[, 3]*w[i,1], emp_cov1[, 3]*w[i,2]), ncol=2))%*%solve(diag(2) + sigma[[i]])%*%t(hh1 - matrix(c(emp_cov1[, 3]*w[i,1], emp_cov1[, 3]*w[i,2]), ncol=2))))
            }else{
              h <- sqrt(diag((hh1 - matrix(c(emp_cov1[, 3]*w[i,1], emp_cov1[, 3]*w[i,2]), ncol=2))%*%solve(diag(2) + sigma[[i]])%*%t(hh1 - matrix(c(emp_cov1[, 3]*w[i,1], emp_cov1[, 3]*w[i,2]), ncol=2))))
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
        tloss <- sum(((emp_cov1[, 6] - theo)/(1.001 - theo))^2)
      
      loss <- loss+tloss
      
      return(loss)
      
    }
  }else if(model == 'gneitingmatern_lagrangian'){
    if(step == 1){
      
      loss <- 0
      
      if(nug_eff == T){
        nug <- theta[6:7]
        lambda <- theta[8]
        scale <- theta[9]
        mu <- theta[10:11]
        nu.L <- theta[12]
      }else{
        nug <- c(0, 0)
        lambda <- theta[6]
        scale <- theta[7]
        mu <- theta[8:9]
        nu.L <- theta[10]
      }
      nu <- theta[1:2]
      beta <- theta[3]
      var <- theta[4:5]
      
      mu1 <- mu[1]
      mu2 <- mu[2]
      mu3 <- (mu[1] + mu[2])/2
      
      nu1 <- nu[1]
      nu2 <- nu[2]
      nu3 <- (nu1 + nu2)/2
      
      
      if( theta[1] < 0.0001 | theta[2] < 0.0001 | theta[3] < 0.0001 | theta[4] < 0.0001 | theta[5] < 0.0001 |theta[6]< 0.001 | theta[6] > 1 | theta[7] < 100 |
          theta[8] < 0.0001 | theta[9] < 0.0001 | theta[10] < 5 ){
        return(Inf)
      }else{
        for(i in 1:2){
          
          theo.temp <- ifelse(h != 0, var[i]*(h/beta)^nu[i] * besselK(h/beta,nu[i])/(2^(nu[i]-1)*gamma(nu[i])), var[i] + nug[i])
          lagrangian <- pmax((1 - h/scale), 0)^(nu.L + mu[i])
          
          theo <- (1 - lambda)*theo.temp + lambda*lagrangian
          if (weights == 1) 
            tloss <- sum((theo - emp_cov1[,i+3])^2)
          if (weights == 2) 
            tloss <- sum(np* (emp_cov1[,i+3] - theo)^2)
          if (weights == 3) 
            tloss <- sum(((emp_cov1[,i+3] - theo)/(1.001-theo))^2)
          
          loss <- loss+tloss
        }
        return(loss)
      }
    }else if (step == 2){
      
      if(nug_eff == T){
        nug <- est_param[6:7]
      }else{
        nug <- c(0, 0)
        lambda <- est_param[6]
        scale <- est_param[7]
        mu <- est_param[8:9]
        nu.L <- est_param[10]
      }
      nu <- est_param[1:2]
      beta <- est_param[3]
      var <- est_param[4:5]
      
      loss <- 0
      
      rho <- theta
      
      mu1 <- mu[1]
      mu2 <- mu[2]
      mu3 <- (mu[1] + mu[2])/2
      
      nu1 <- nu[1]
      nu2 <- nu[2]
      nu3 <- (nu1 + nu2)/2
      
      if( abs(rho) > 1){
        return(Inf)
      }else{
        for(i in 3:3){
          
          theo.temp <- ifelse(h !=0 ,(h/beta)^nu3 * besselK(h/beta, nu3)/(2^(nu3 - 1)*gamma(nu3))*sqrt(var[1] * var[2])*rho, sqrt(var[1] * var[2])*rho)
          
          beta2.3 <- (gamma(1 + mu3)/gamma(1 + nu.L + mu3))*sqrt((gamma(1 + nu.L + mu1)*gamma(1 + nu.L + mu2))/(gamma(1 + mu1)*gamma(1 + mu2)))
          lagrangian <- beta2.3*pmax((1 - h/scale), 0)^(nu.L + mu3)
          theo <- (1 - lambda)*theo.temp + lambda*lagrangian
          
          if (weights == 1) 
            tloss <- sum((theo - emp_cov1[,i+3])^2)
          if (weights == 2) 
            tloss <- sum(np* (emp_cov1[,i+3] - theo)^2)
          if (weights == 3) 
            tloss <- sum(((emp_cov1[,i+3] - theo)/(1.001 - theo))^2)
          loss <- loss + tloss
        }
        return(loss)
      }
    }else{
      
      if(nug_eff == T){
        nug <- est_param[6:7]
      }else{
        nug <- c(0, 0)
        lambda <- est_param[6]
        scale <- est_param[7]
        mu <- est_param[8:9]
        nu.L <- est_param[10]
        rho <- est_param[11]
      }
      nu <- est_param[1:2]
      beta <- est_param[3]
      var <- est_param[4:5]
      
      loss <- 0
    
      alpha <- theta[1]
      b <- theta[2]
      w <- theta[3:4]
      loss <- 0
      
      new_h <- sqrt(hh[,1]^2+hh[,2]^2)/1000
      h <- new_h/(alpha+1)^(b/2)
      t <- emp_cov1[,3]
      h2 <- sqrt((hh[,1]-w[1])^2+(hh[,2]-w[2])^2)/1000
      
      if(meters == T){
        h <- sqrt(hh[,1]^2 + hh[,2]^2)/1000
        h2 <- sqrt((hh[,1] - emp_cov1[,3]*w[1])^2 + (hh2[,2] - emp_cov1[,3]*w[2])^2)/1000
        h <- h/(alpha+1)^(b/2)
      }else{
        h <- sqrt(hh[,1]^2 + hh[,2]^2)
        h2 <- sqrt((hh2[,1] - emp_cov1[,3]*w[1])^2 + (hh2[,2] - emp_cov1[,3]*w[2])^2)
        h <- h/(alpha+1)^(b/2)
      }
      
      mu1 <- mu[1]
      mu2 <- mu[2]
      mu3 <- (mu[1] + mu[2])/2
      
      nu1 <- nu[1]
      nu2 <- nu[2]
      nu3 <- (nu1 + nu2)/2
      
      if( theta[1] < 0.0001 | theta[2] < 0.0001 | theta[2] > 1){
        return(Inf)
      }else{
        for(i in 1:2){
          theo.temp <- ifelse(h != 0, var[i]*(h/beta)^nu[i] * besselK(h/beta, nu[i])/(2^(nu[i] - 1)*gamma(nu[i]))/(alpha + 1), (var[i] + nug[i])/(alpha + 1))
          lagrangian <- pmax((1 - h2/scale), 0)^(nu.L + mu[i])
          
          theo <- (1 - lambda)*theo.temp + lambda*lagrangian
          if (weights == 1) 
            tloss <- sum((theo - emp_cov1[,i+3])^2)
          if (weights == 2) 
            tloss <- sum(np* (emp_cov1[,i+3] - theo)^2)
          if (weights == 3) 
            tloss <- sum(((emp_cov1[,i+3] - theo)/(1.001 - theo))^2)
          
          loss <- loss+tloss
        }
        for(i in 3:3){
          theo.temp <- ifelse(h != 0, (h/beta)^nu3 * besselK(h/beta, nu3)/(2^(nu3 - 1)*gamma(nu3))*sqrt(var[1] * var[2])*rho/((alpha + 1)), sqrt(var[1] * var[2])*rho/((alpha + 1)))
          
          beta2.3 <- (gamma(1 + mu3)/gamma(1 + nu.L + mu3))*sqrt((gamma(1 + nu.L + mu1)*gamma(1 + nu.L + mu2))/(gamma(1 + mu1)*gamma(1 + mu2)))
          lagrangian <- beta2.3*pmax((1 - h2/scale), 0)^(nu.L+mu3)
          theo <- (1 - lambda)*theo.temp + lambda*lagrangian
          
          if (weights == 1) 
            tloss <- sum((theo - emp_cov1[,i+3])^2)
          if (weights == 2) 
            tloss <- sum(np* (emp_cov1[,i+3] - theo)^2)
          if (weights == 3) 
            tloss <- sum(((emp_cov1[,i+3] - theo)/(1.001 - theo))^2)
          loss <- loss + tloss
        }
        return(loss)
      }
    }
  }
}

data_plots_functions <- function(filename, var1, var2){
  
  colfunc <- colorRampPalette(c("blue", "yellow", "red"))
  
  pdf(file = filename, width = 20, height = 9)
  
  starting_pt <- 0
  
  zr <- round(range(c(var1[, 1:(starting_pt + 1:5)[5]], var2[, 1:(starting_pt + 1:5)[5]])),1)
  endpt <- ceil(max(abs(min(zr)), abs(max(zr))))
  col_ind <- seq(-endpt, endpt, length.out = 1000)
  
  split.screen(rbind(c(0.1,0.85,0,1), c(.95,0.99,0,0.98)))
  screen(1)
  layout(mat = matrix(1:10, nrow = 2, ncol = 5, byrow = T))
  
  for(bb in starting_pt + 1:5){
    xx <- findInterval(var1[, bb], col_ind)
    
    if(bb==1){
      par(pty="s") 
      par(mai=c(0.4,0.4,0.4,0.4))
      map("worldHires", xlim = c(range(grid_locations[,1])[1]-1,range(grid_locations[,1])[2]+2), ylim = c(range(grid_locations[,2])[1]-1, range(grid_locations[,2])[2]+2), lwd=0.5)
      plot(spdf, add=TRUE, pch=16, col=colfunc(1000)[xx], lwd=0.5)
      mtext(expression(Z[1]), side = 2, line = 0, adj = 0.5, cex = 1.5, font=3,col="#0086FF")
      mtext(paste("t = ", bb - starting_pt, sep = ""), side = 3, line = 1, adj = 0.5, cex = 1.5, font=2,col="#4EC1DE")
    }else{
      par(pty="s") 
      par(mai=c(0.4,0.4,0.4,0.4))
      map("worldHires", xlim = c(range(grid_locations[,1])[1]-1,range(grid_locations[,1])[2]+2), ylim = c(range(grid_locations[,2])[1]-1,range(grid_locations[,2])[2]+2), lwd=0.5)
      plot(spdf, add=TRUE, pch=16, col=colfunc(1000)[xx], lwd=0.5)
      mtext(paste("t = ", bb - starting_pt,sep=""), side = 3, line = 1, adj = 0.5, cex = 1.5, font=2,col="#4EC1DE")
    }
    
  }
  for(bb in starting_pt + 1:5){
    xx <- findInterval(var2[, bb], col_ind)
    
    if(bb==1){
      par(pty="s") 
      par(mai=c(0.4,0.4,0.4,0.4))
      map("worldHires", xlim = c(range(grid_locations[,1])[1]-1,range(grid_locations[,1])[2]+2), ylim = c(range(grid_locations[,2])[1]-1,range(grid_locations[,2])[2]+2), lwd=0.5)
      plot(spdf, add=TRUE, pch=16, col=colfunc(1000)[xx], lwd=0.5)
      mtext(expression(Z[2]), side = 2, line = 0, adj = 0.5, cex = 1.5, font=3,col="#0086FF")
    }else{
      par(pty="s") 
      par(mai=c(0.4,0.4,0.4,0.4))
      map("worldHires", xlim = c(range(grid_locations[,1])[1]-1,range(grid_locations[,1])[2]+2), ylim = c(range(grid_locations[,2])[1]-1,range(grid_locations[,2])[2]+2), lwd=0.5)
      plot(spdf, add=TRUE, pch=16, col=colfunc(1000)[xx], lwd=0.5)
      axis(1, at=seq(0,1,by=0.2), labels = seq(0,1,by=0.2))
    }
  }
  screen(2)
  screen(2)
  x1 <- c(0.3, 0.4, 0.4, 0.3)
  y1 <- c(0.35, 0.35, 0.6, 0.6)
  legend.gradient2(cbind(x1, y1), cols = colorRampPalette(colors) (50), title = " ", 
                   limits = seq(-endpt, endpt, length.out = 5), cex=1)
  close.screen( all=TRUE)
  dev.off()
  
}

whichpart <- function(x, n=30) {
  nx <- length(x)
  p <- nx-n
  xp <- sort(x, partial=p)[p]
  which(x > xp)
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
