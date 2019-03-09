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

