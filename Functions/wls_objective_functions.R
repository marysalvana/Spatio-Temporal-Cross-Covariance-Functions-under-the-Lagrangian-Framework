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

wls3_for_mod5<-function(theta, emp_cov1, weights) {
  
  nu=tempo_theta[1:2]
  beta=tempo_theta[3]
  nug <- c(0,0)
  var <- c(1,1)
  
  rho <- theta[1]
  
  loss<- 0
  
  R <- matrix(c(-1.1339125 *cos(116.9811232),-1.1339125 *sin(116.9811232),0.8617141*sin(116.9811232),
                -0.8617141*cos(116.9811232)),ncol=2,byrow=T)
  
  nu1 <- nu[1]
  nu2 <- nu[2]
  nu3 <- (nu[1]+nu[2])/2
  
  hh=emp_cov1[,1:2]%*%R
  hh <- cbind(hh[,1],hh[,2])
  h <- sqrt((hh[,1])^2+(hh[,2])^2)/1000
  
  np <- ifelse(h!=0,1/h,1)
  
  if( abs(rho)>1 ){
    return(Inf)
  }else{
    
    theo <- ifelse(h!=0,rho*sqrt(var[1]*var[2])*(h/beta)^nu3*besselK(h/beta,nu3)/(2^(nu3-1)*gamma(nu3)),rho*sqrt(var[1]*var[2]))
    if (weights == 1) 
      tloss <- sum((theo-emp_cov1[,3+3])^2)
    if (weights == 2) 
      tloss <- sum(np* (emp_cov1[,3+3] - theo)^2)
    if (weights == 3) 
      tloss<-sum(((emp_cov1[,3+3]-theo)/(1.001-theo))^2)
    loss <- loss+tloss
    
    return(loss)
  }
}

m2_with_decay<-function(theta, emp_cov1,weights){
  
  nu=tempo_theta[1:2]
  beta=tempo_theta[3]
  nug <- c(0,0)
  var <- c(1,1)
  rho <- mod1_parms[4]
  kappa <- mod1_parms[5:6]
  loss<- 0
  
  R <- matrix(c(-1.1339125 *cos(116.9811232),-1.1339125 *sin(116.9811232),0.8617141*sin(116.9811232),
                -0.8617141*cos(116.9811232)),ncol=2,byrow=T)
  
  w <- theta[1:2]
  c <- theta[3:4]
  
  hh= emp_cov1[,1:2]%*%R
  h <- sqrt((hh[,1]-w[1])^2+(hh[,2]-w[2])^2)/1000
  
  hh2 <- cbind(hh[,1]-kappa[1]-w[1],hh[,2]-kappa[2]-w[2])
  h2 <- sqrt((hh2[,1])^2+(hh2[,2])^2)/1000
  
  nu1 <- nu[1]
  nu2 <- nu[2]
  nu3 <- (nu[1]+nu[2])/2
  
  c3 <- (c[1]+c[2])/2
  
  if(c[1] <= 0 | c[2] <= 0 ){
    return(Inf)
  }else{
    for(i in 1:2){

      theo <- ifelse(h!=0,var[i]*(h/beta)^nu[i]*besselK(h/beta,nu[i])/(2^(nu[i]-1)*gamma(nu[i])),var[i]+nug[i])
      if (weights == 1) 
        tloss <- sum((theo-emp_cov1[,i+3])^2)
      if (weights == 2) 
        tloss <- sum(np* (emp_cov1[,i+3] - theo)^2)
      if (weights == 3) 
        tloss<-sum(((emp_cov1[,i+3]-exp(-c[i])*theo)/(1.001-exp(-c[i])*theo))^2)
      
      loss <- loss+tloss
    }
    for(i in 3:3){
      
      theo <- ifelse(h2!=0,rho*sqrt(var[1]*var[2])*(h2/beta)^nu3*besselK(h2/beta,nu3)/(2^(nu3-1)*gamma(nu3)),rho*sqrt(var[1]*var[2]))
      if (weights == 1) 
        tloss <- sum((theo-emp_cov1[,i+3])^2)
      if (weights == 2) 
        tloss <- sum(np * (emp_cov1[,i+3] - theo)^2)
      if (weights == 3) 
        tloss<-sum(((emp_cov1[,i+3]-exp(-c3)*theo)/(1.001-exp(-c3)*theo))^2)
      loss <- loss+tloss
    }
  }
  return(loss)
}

m2_no_decay<-function(theta, emp_cov1,weights){
  
  nu=tempo_theta[1:2]
  beta=tempo_theta[3]
  nug <- c(0,0)
  var <- c(1,1)
  rho <- mod1_parms[4]
  loss<- 0
  
  w <- theta[1:2]
  
  hh= emp_cov1[,1:2]
  h <- h2 <-sqrt((hh[,1]-w[1])^2+(hh[,2]-w[2])^2)/1000
  
  nu1 <- nu[1]
  nu2 <- nu[2]
  nu3 <- (nu[1]+nu[2])/2
  
  for(i in 1:2){
    
    theo <- ifelse(h!=0,var[i]*(h/beta)^nu[i]*besselK(h/beta,nu[i])/(2^(nu[i]-1)*gamma(nu[i])),var[i]+nug[i])
    if (weights == 1) 
      tloss <- sum((theo-emp_cov1[,i+3])^2)
    if (weights == 2) 
      tloss <- sum(np* (emp_cov1[,i+3] - theo)^2)
    if (weights == 3) 
      tloss<-sum(((emp_cov1[,i+3]-theo)/(1.001-theo))^2)
    
    loss <- loss+tloss
  }
  for(i in 3:3){
    
    theo <- ifelse(h2!=0,rho*sqrt(var[1]*var[2])*(h2/beta)^nu3*besselK(h2/beta,nu3)/(2^(nu3-1)*gamma(nu3)),rho*sqrt(var[1]*var[2]))
    if (weights == 1) 
      tloss <- sum((theo-emp_cov1[,i+3])^2)
    if (weights == 2) 
      tloss <- sum(np * (emp_cov1[,i+3] - theo)^2)
    if (weights == 3) 
      tloss<-sum(((emp_cov1[,i+3]-theo)/(1.001-theo))^2)
    loss <- loss+tloss
  }
  return(loss)
}

m4<-function(theta, emp_cov1, weights) {
  
  nu=theta[1:2]
  beta=theta[3]
  nug <- c(0,0)
  var <- c(1,1)
  
  lambda <- theta[4]
  scale <- theta[5]
  mu=theta[6:7]
  nu.L=theta[8]
  
  alpha <- theta[9]
  b <- theta[10]
  w <- theta[11:12]
  #rho <- (0.52-lambda)/(1-lambda)
  rho <- theta[13]
  loss<- 0
  
  R <- matrix(c(-1.1339125 *cos(116.9811232),-1.1339125 *sin(116.9811232),0.8617141*sin(116.9811232),
                -0.8617141*cos(116.9811232)),ncol=2,byrow=T)
  
  hh= emp_cov1[,1:2]%*%R
  
  new_h <- sqrt(hh[,1]^2+hh[,2]^2)/1000
  h <- new_h/(alpha+1)^(b/2)
  t <- emp_cov1[,3]
  h2 <- sqrt((hh[,1]-w[1])^2+(hh[,2]-w[2])^2)/1000
  
  mu1=mu[1]
  mu2=mu[2]
  mu3=(mu[1]+mu[2])/2
  
  nu1=nu[1]
  nu2=nu[2]
  nu3=(nu[1]+nu[2])/2
  
  if( theta[1] < 0.0001 | theta[2] < 0.0001 | theta[3] < 0.0001 |theta[4]<0.001 | theta[4]>1 | theta[5]<100 |
      theta[6] < 0.0001 | theta[7] < 0.0001 | theta[8] <5 | theta[9]<0.0001 | theta[10] < 0.0001 | theta[10]>1 |
      abs(rho)>1){
    return(Inf)
  }else{
    for(i in 1:2){
      h_star <- h
      theo.temp <- ifelse(h!=0,var[i]*(h_star/beta)^nu[i] * besselK(h_star/beta,nu[i])/(2^(nu[i]-1)*gamma(nu[i]))/(alpha+1),(var[i]+nug[i])/(alpha+1))
      lagrangian=pmax((1-h2/scale),0)^(nu.L+mu[i])
      
      theo <- (1-lambda)*theo.temp+lambda*lagrangian
      if (weights == 1) 
        tloss <- sum((theo-emp_cov1[,i+3])^2)
      if (weights == 2) 
        tloss <- sum(np* (emp_cov1[,i+3] - theo)^2)
      if (weights == 3) 
        tloss<-sum(((emp_cov1[,i+3]-theo)/(1.001-theo))^2)
      
      loss <- loss+tloss
    }
    for(i in 3:3){
      
      h_star <- h
      
      theo.temp <- ifelse(h!=0,(h_star/beta)^nu3 * besselK(h_star/beta,nu3)/(2^(nu3-1)*gamma(nu3))*sqrt(var[1] * var[2])*rho/((alpha+1)),sqrt(var[1] * var[2])*rho/((alpha+1)))
      
      beta2.3=(gamma(1+mu3)/gamma(1+nu.L+mu3))*sqrt((gamma(1+nu.L+mu1)*gamma(1+nu.L+mu2))/(gamma(1+mu1)*gamma(1+mu2)))
      lagrangian=beta2.3*pmax((1-h2/scale),0)^(nu.L+mu3)
      theo <- (1-lambda)*theo.temp+lambda*lagrangian
      
      if (weights == 1) 
        tloss <- sum((theo-emp_cov1[,i+3])^2)
      if (weights == 2) 
        tloss <- sum(np* (emp_cov1[,i+3] - theo)^2)
      if (weights == 3) 
        tloss<-sum(((emp_cov1[,i+3]-theo)/(1.001-theo))^2)
      loss <- loss+tloss
    }
    return(loss)
  }
}

m4_step1<-function(theta, emp_cov1, weights) {
  
  nu=theta[1:2]
  beta=theta[3]
  nug <- c(0,0)
  var <- c(1,1)
  
  lambda <- theta[4]
  scale <- theta[5]
  mu=theta[6:7]
  nu.L=theta[8]
  
  loss<- 0
  
  R <- matrix(c(-1.1339125 *cos(116.9811232),-1.1339125 *sin(116.9811232),0.8617141*sin(116.9811232),
                -0.8617141*cos(116.9811232)),ncol=2,byrow=T)
  
  hh= emp_cov1[,1:2]%*%R
  
  h <- sqrt(hh[,1]^2+hh[,2]^2)/1000
  
  mu1=mu[1]
  mu2=mu[2]
  mu3=(mu[1]+mu[2])/2
  
  nu1=nu[1]
  nu2=nu[2]
  nu3=(nu[1]+nu[2])/2
  
  if( theta[1] < 0.0001 | theta[2] < 0.0001 | theta[3] < 0.0001 |theta[4]<0.001 | theta[4]>1 | theta[5]<100 |
      theta[6] < 0.0001 | theta[7] < 0.0001 | theta[8] <5 ){
    return(Inf)
  }else{
    for(i in 1:2){
  
      theo.temp <- ifelse(h!=0,var[i]*(h/beta)^nu[i] * besselK(h/beta,nu[i])/(2^(nu[i]-1)*gamma(nu[i])),var[i]+nug[i])
      lagrangian=pmax((1-h/scale),0)^(nu.L+mu[i])
      
      theo <- (1-lambda)*theo.temp+lambda*lagrangian
      if (weights == 1) 
        tloss <- sum((theo-emp_cov1[,i+3])^2)
      if (weights == 2) 
        tloss <- sum(np* (emp_cov1[,i+3] - theo)^2)
      if (weights == 3) 
        tloss<-sum(((emp_cov1[,i+3]-theo)/(1.001-theo))^2)
      
      loss <- loss+tloss
    }
    return(loss)
  }
}

m4_step2<-function(theta, emp_cov1, weights) {
  
  nu=mod4_parms[1:2]
  beta=mod4_parms[3]
  nug <- c(0,0)
  var <- c(1,1)
  
  lambda <- mod4_parms[4]
  scale <- mod4_parms[5]
  mu=mod4_parms[6:7]
  nu.L=mod4_parms[8]
  
  #rho <- (0.52-lambda)/(1-lambda)
  rho <- theta
  loss<- 0
  
  R <- matrix(c(-1.1339125 *cos(116.9811232),-1.1339125 *sin(116.9811232),0.8617141*sin(116.9811232),
                -0.8617141*cos(116.9811232)),ncol=2,byrow=T)
  
  hh= emp_cov1[,1:2]%*%R
  
  h <- sqrt(hh[,1]^2+hh[,2]^2)/1000
  
  mu1=mu[1]
  mu2=mu[2]
  mu3=(mu[1]+mu[2])/2
  
  nu1=nu[1]
  nu2=nu[2]
  nu3=(nu[1]+nu[2])/2
  
  if( abs(rho)>1){
    return(Inf)
  }else{
    for(i in 3:3){
      
      theo.temp <- ifelse(h!=0,(h/beta)^nu3 * besselK(h/beta,nu3)/(2^(nu3-1)*gamma(nu3))*sqrt(var[1] * var[2])*rho,sqrt(var[1] * var[2])*rho)
      
      beta2.3=(gamma(1+mu3)/gamma(1+nu.L+mu3))*sqrt((gamma(1+nu.L+mu1)*gamma(1+nu.L+mu2))/(gamma(1+mu1)*gamma(1+mu2)))
      lagrangian=beta2.3*pmax((1-h/scale),0)^(nu.L+mu3)
      theo <- (1-lambda)*theo.temp+lambda*lagrangian
      
      if (weights == 1) 
        tloss <- sum((theo-emp_cov1[,i+3])^2)
      if (weights == 2) 
        tloss <- sum(np* (emp_cov1[,i+3] - theo)^2)
      if (weights == 3) 
        tloss<-sum(((emp_cov1[,i+3]-theo)/(1.001-theo))^2)
      loss <- loss+tloss
    }
    return(loss)
  }
}

m4_step3<-function(theta, emp_cov1, weights) {
  
  nu=mod4_parms[1:2]
  beta=mod4_parms[3]
  nug <- c(0,0)
  var <- c(1,1)
  
  lambda <- mod4_parms[4]
  scale <- mod4_parms[5]
  mu=mod4_parms[6:7]
  nu.L=mod4_parms[8]
  #mu=c(0.639,0.828)
  #nu.L=8
  
  alpha <- theta[1]
  b <- theta[2]
  w <- theta[3:4]
  rho <- mod4_parms[9]
  loss<- 0
  
  R <- matrix(c(-1.1339125 *cos(116.9811232),-1.1339125 *sin(116.9811232),0.8617141*sin(116.9811232),
                -0.8617141*cos(116.9811232)),ncol=2,byrow=T)
  
  hh= emp_cov1[,1:2]%*%R
  
  new_h <- sqrt(hh[,1]^2+hh[,2]^2)/1000
  h <- new_h/(alpha+1)^(b/2)
  t <- emp_cov1[,3]
  h2 <- sqrt((hh[,1]-w[1])^2+(hh[,2]-w[2])^2)/1000
  
  mu1=mu[1]
  mu2=mu[2]
  mu3=(mu[1]+mu[2])/2
  
  nu1=nu[1]
  nu2=nu[2]
  nu3=(nu[1]+nu[2])/2
  
  if(  theta[1]<0.0001 | theta[2] < 0.0001 | theta[2]>1){
    return(Inf)
  }else{
    for(i in 1:2){
      h_star <- h
      theo.temp <- ifelse(h!=0,var[i]*(h_star/beta)^nu[i] * besselK(h_star/beta,nu[i])/(2^(nu[i]-1)*gamma(nu[i]))/(alpha+1),(var[i]+nug[i])/(alpha+1))
      lagrangian=pmax((1-h2/scale),0)^(nu.L+mu[i])
      
      theo <- (1-lambda)*theo.temp+lambda*lagrangian
      if (weights == 1) 
        tloss <- sum((theo-emp_cov1[,i+3])^2)
      if (weights == 2) 
        tloss <- sum(np* (emp_cov1[,i+3] - theo)^2)
      if (weights == 3) 
        tloss<-sum(((emp_cov1[,i+3]-theo)/(1.001-theo))^2)
      
      loss <- loss+tloss
    }
    for(i in 3:3){
      
      h_star <- h
      
      theo.temp <- ifelse(h!=0,(h_star/beta)^nu3 * besselK(h_star/beta,nu3)/(2^(nu3-1)*gamma(nu3))*sqrt(var[1] * var[2])*rho/((alpha+1)),sqrt(var[1] * var[2])*rho/((alpha+1)))
      
      beta2.3=(gamma(1+mu3)/gamma(1+nu.L+mu3))*sqrt((gamma(1+nu.L+mu1)*gamma(1+nu.L+mu2))/(gamma(1+mu1)*gamma(1+mu2)))
      lagrangian=beta2.3*pmax((1-h2/scale),0)^(nu.L+mu3)
      theo <- (1-lambda)*theo.temp+lambda*lagrangian
      
      if (weights == 1) 
        tloss <- sum((theo-emp_cov1[,i+3])^2)
      if (weights == 2) 
        tloss <- sum(np* (emp_cov1[,i+3] - theo)^2)
      if (weights == 3) 
        tloss<-sum(((emp_cov1[,i+3]-theo)/(1.001-theo))^2)
      loss <- loss+tloss
    }
    return(loss)
  }
}

m4_step1_for_sim<-function(theta, emp_cov1, weights) {
  
  nu=theta[1:2]
  beta=theta[3]
  nug <- c(0,0)
  var <- c(1,1)
  
  lambda <- theta[4]
  scale <- theta[5]
  mu=theta[6:7]
  nu.L=theta[8]
  
  loss<- 0
  
  h <- emp_cov1[,1]
  
  mu1=mu[1]
  mu2=mu[2]
  mu3=(mu[1]+mu[2])/2
  
  nu1=nu[1]
  nu2=nu[2]
  nu3=(nu[1]+nu[2])/2
  
  if( theta[1] < 0.0001 | theta[2] < 0.0001 | theta[3] < 0.0001 |theta[4]<0.001 | theta[4]>1 | theta[5]<100 |
      theta[6] < 0.0001 | theta[7] < 0.0001 | theta[8] <5 ){
    return(Inf)
  }else{
    for(i in 1:2){
      
      theo.temp <- ifelse(h!=0,var[i]*(h/beta)^nu[i] * besselK(h/beta,nu[i])/(2^(nu[i]-1)*gamma(nu[i])),var[i]+nug[i])
      lagrangian=pmax((1-h/scale),0)^(nu.L+mu[i])
      
      theo <- (1-lambda)*theo.temp+lambda*lagrangian
      if (weights == 1) 
        tloss <- sum((theo-emp_cov1[,i+1])^2)
      if (weights == 2) 
        tloss <- sum(np* (emp_cov1[,i+1] - theo)^2)
      if (weights == 3) 
        tloss<-sum(((emp_cov1[,i+1]-theo)/(1.001-theo))^2)
      
      loss <- loss+tloss
    }
    return(loss)
  }
}

m4_step2_for_sim<-function(theta, emp_cov1, weights) {
  
  nu=mod4_parms[1:2]
  beta=mod4_parms[3]
  nug <- c(0,0)
  var <- c(1,1)
  
  lambda <- mod4_parms[4]
  scale <- mod4_parms[5]
  mu=mod4_parms[6:7]
  nu.L=mod4_parms[8]
  
  #rho <- (0.52-lambda)/(1-lambda)
  rho <- theta
  loss<- 0
  
  h <- emp_cov1[,1]
  
  mu1=mu[1]
  mu2=mu[2]
  mu3=(mu[1]+mu[2])/2
  
  nu1=nu[1]
  nu2=nu[2]
  nu3=(nu[1]+nu[2])/2
  
  if( abs(rho)>1){
    return(Inf)
  }else{
    for(i in 3:3){
      
      theo.temp <- ifelse(h!=0,(h/beta)^nu3 * besselK(h/beta,nu3)/(2^(nu3-1)*gamma(nu3))*sqrt(var[1] * var[2])*rho,sqrt(var[1] * var[2])*rho)
      
      beta2.3=(gamma(1+mu3)/gamma(1+nu.L+mu3))*sqrt((gamma(1+nu.L+mu1)*gamma(1+nu.L+mu2))/(gamma(1+mu1)*gamma(1+mu2)))
      lagrangian=beta2.3*pmax((1-h/scale),0)^(nu.L+mu3)
      theo <- (1-lambda)*theo.temp+lambda*lagrangian
      
      if (weights == 1) 
        tloss <- sum((theo-emp_cov1[,i+1])^2)
      if (weights == 2) 
        tloss <- sum(np* (emp_cov1[,i+1] - theo)^2)
      if (weights == 3) 
        tloss<-sum(((emp_cov1[,i+1]-theo)/(1.001-theo))^2)
      loss <- loss+tloss
    }
    return(loss)
  }
}

m4_step3_for_sim<-function(theta, emp_cov1, weights) {
  
  nu=mod4_parms[1:2]
  beta=mod4_parms[3]
  nug <- c(0,0)
  var <- c(1,1)
  
  lambda <- mod4_parms[4]
  scale <- mod4_parms[5]
  mu=mod4_parms[6:7]
  nu.L=mod4_parms[8]
  #mu=c(0.639,0.828)
  #nu.L=8
  
  alpha <- theta[1]
  b <- theta[2]
  w <- theta[3:4]
  rho <- mod4_params.temp2[9]
  loss<- 0
  
  hh= emp_cov1[,1:2]
  
  new_h <- sqrt(hh[,1]^2+hh[,2]^2)/1000
  h <- new_h/(alpha+1)^(b/2)
  t <- emp_cov1[,3]
  h2 <- sqrt((hh[,1]-w[1])^2+(hh[,2]-w[2])^2)/1000
  
  mu1=mu[1]
  mu2=mu[2]
  mu3=(mu[1]+mu[2])/2
  
  nu1=nu[1]
  nu2=nu[2]
  nu3=(nu[1]+nu[2])/2
  
  if(  theta[1]<0.0001 | theta[2] < 0.0001 | theta[2]>1){
    return(Inf)
  }else{
    for(i in 1:2){
      h_star <- h
      theo.temp <- ifelse(h!=0,var[i]*(h_star/beta)^nu[i] * besselK(h_star/beta,nu[i])/(2^(nu[i]-1)*gamma(nu[i]))/(alpha+1),(var[i]+nug[i])/(alpha+1))
      lagrangian=pmax((1-h2/scale),0)^(nu.L+mu[i])
      
      theo <- (1-lambda)*theo.temp+lambda*lagrangian
      if (weights == 1) 
        tloss <- sum((theo-emp_cov1[,i+3])^2)
      if (weights == 2) 
        tloss <- sum(np* (emp_cov1[,i+3] - theo)^2)
      if (weights == 3) 
        tloss<-sum(((emp_cov1[,i+3]-theo)/(1.001-theo))^2)
      
      loss <- loss+tloss
    }
    for(i in 3:3){
      
      h_star <- h
      
      theo.temp <- ifelse(h!=0,(h_star/beta)^nu3 * besselK(h_star/beta,nu3)/(2^(nu3-1)*gamma(nu3))*sqrt(var[1] * var[2])*rho/((alpha+1)),sqrt(var[1] * var[2])*rho/((alpha+1)))
      
      beta2.3=(gamma(1+mu3)/gamma(1+nu.L+mu3))*sqrt((gamma(1+nu.L+mu1)*gamma(1+nu.L+mu2))/(gamma(1+mu1)*gamma(1+mu2)))
      lagrangian=beta2.3*pmax((1-h2/scale),0)^(nu.L+mu3)
      theo <- (1-lambda)*theo.temp+lambda*lagrangian
      
      if (weights == 1) 
        tloss <- sum((theo-emp_cov1[,i+3])^2)
      if (weights == 2) 
        tloss <- sum(np* (emp_cov1[,i+3] - theo)^2)
      if (weights == 3) 
        tloss<-sum(((emp_cov1[,i+3]-theo)/(1.001-theo))^2)
      loss <- loss+tloss
    }
    return(loss)
  }
}


m5<-function(theta, emp_cov1, weights) {
  
  nu=mod1_parms[1:2]
  beta=mod1_parms[3]
  nug <- c(0,0)
  var <- c(1,1)
  rho <- 0.52
  
  nu1 <- nu[1]
  nu2 <- nu[2]
  nu3 <- (nu[1]+nu[2])/2
  
  alpha <- theta[1]
  b <- theta[2]
  loss<- 0
  
  R <- matrix(c(-1.1339125 *cos(116.9811232),-1.1339125 *sin(116.9811232),0.8617141*sin(116.9811232),
                -0.8617141*cos(116.9811232)),ncol=2,byrow=T)
  
  hh= emp_cov1[,1:2]%*%R
  
  new_h <- sqrt(hh[,1]^2+hh[,2]^2)/1000
  h <- new_h/(alpha+1)^(b/2)
  t <- emp_cov1[,3]
  
  np <- 1/h
  np[!is.finite(np)]=1
  
  if( theta[1]<0.0001 | theta[2]<0.0001 | theta[2]>1  ){
    return(Inf)
  }else{
    for(i in 1:2){
      h_star <- h
      theo <- ifelse(h!=0,var[i]*(h_star/beta)^nu[i] * besselK(h_star/beta,nu[i])/(2^(nu[i]-1)*gamma(nu[i]))/(alpha+1),(var[i]+nug[i])/(alpha+1))
      if (weights == 1) 
        tloss <- sum((theo-emp_cov1[,i+3])^2)
      if (weights == 2) 
        tloss <- sum(np* (emp_cov1[,i+3] - theo)^2)
      if (weights == 3) 
        tloss<-sum(((emp_cov1[,i+3]-theo)/(1.001-theo))^2)
      
      loss <- loss+tloss
    }
    for(i in 3:3){
      
      h_star <- h
      theo <- ifelse(h!=0,(h_star/beta)^nu3 * besselK(h_star/beta,nu3)/(2^(nu3-1)*gamma(nu3))*sqrt(var[1] * var[2])*rho/((alpha+1)),sqrt(var[1] * var[2])*rho/((alpha+1)))
      if (weights == 1) 
        tloss <- sum((theo-emp_cov1[,i+3])^2)
      if (weights == 2) 
        tloss <- sum(np* (emp_cov1[,i+3] - theo)^2)
      if (weights == 3) 
        tloss<-sum(((emp_cov1[,i+3]-theo)/(1.001-theo))^2)
      
      loss <- loss+tloss
    }
    return(loss)
  }
}

