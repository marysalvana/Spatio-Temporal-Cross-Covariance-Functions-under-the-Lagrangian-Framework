WORKSTATION = 4

# 0.1 SET WORKING DIRECTORY WHERE DATA IS SAVED

if(WORKSTATION == 1){
  root <- '/Volumes/GoogleDrive/My Drive/Summer 2017_2018/Project 1/'
}else if (WORKSTATION == 2){
  root <- 'G://My Drive/Summer 2017_2018/Project 1/'
}else if (WORKSTATION == 3){
  root <- '/Users/salvanmo/Dropbox/05SalvanaMary/Spatio-Temporal Cross-Covariance Functions under the Lagrangian Framework/'
}else{
  root <- '/Volumes/GoogleDrive/My Drive/Summer 2017_2018/Project 2/Codes/FINAL'
}

source(file=paste(root,"/covariance_plots_functions.R",sep=''))
source(file=paste(root,"/matern_cov.R",sep=''))
source(file=paste(root,"/toeplitz_mat.R",sep=''))
source(file=paste(root,"/load_packages.R",sep=''))
source(file=paste(root,"/legend.gradient2.R",sep=''))

grid_x = seq(from=0.5, to=1, length.out = 11)
grid_y = seq(from=0.5, to=1, length.out = 11)
sim_grid_locations = expand.grid(grid_x, grid_y)
sim_grid_locations <- matrix(cbind(sim_grid_locations[,1],sim_grid_locations[,2]),ncol=2)

n <- nrow(sim_grid_locations)
loc_mat_x <- matrix(sim_grid_locations[,1],ncol=length(grid_x),nrow=length(grid_x),byrow=T)
loc_mat_y <- matrix(sim_grid_locations[,2],ncol=length(grid_x),nrow=length(grid_x),byrow=T)

n <- nrow(sim_grid_locations)

loc <- coords <- sim_grid_locations
t=5

if (t==1){
  loc <- loc
} else {
  for (tt in 1:(t-1)){
    temploc <- matrix(,ncol=2,nrow=nrow(coords))
    for(rr in 1:nrow(coords)){
      temploc[rr,] <- c(coords[rr,1]+tt*w[1],coords[rr,2]+tt*w[2])
    }
    loc <- rbind(loc, temploc)
  }
}

locations <- loc
##-----------------------fig1a-----------------------##

nu=c(0.5,1.5)
beta=2
nug <- c(0,0)
var <- c(1,1)
rot <- rho <- c(0.8)

nu1 <- nu[1]
nu2 <- nu[2]
nu3 <- (nu[1]+nu[2])/2
w=c(0.05,0.05)

sim_cov <- list()
sim_cov2 <- list()
sim_cov3 <- list()

theta2 <- function (n,beta0,beta1,beta2,beta3,beta4) {
  theta3 <- matrix(0,nrow=nrow(locations),ncol=1)
  
  for(i in (1 : nrow(locations))){
    x1 <- locations[i,1]
    x2 <- locations[i,2]
    theta3[i,] <- beta0 + beta1*(x1 - .5) + beta2*(x2-.5) + 
      beta3*(x1 - .5)^2 + beta4*(x2 - .5)^2
  }
  return(theta3)
}

log.lam1<-theta2(n,-3,0,0,-6,-7)
log.lam2<-theta2(n,-5,0,0,6,-4)
logit.phi<-theta2(n,0,1,-2,0,1)

kernel.local <- array(0, dim = c(2, 2, nrow(locations)))
for(i in 1:nrow(locations)){
  lam1 <- exp(log.lam1[i,])
  lam2 <- exp(log.lam2[i,])
  phi <- (pi/2)*exp(logit.phi[i,])/(1+exp(logit.phi[i,]))
  Pmat <- matrix(c(cos(phi), -sin(phi), sin(phi), cos(phi)), nrow = 2, byrow = T)
  Dmat <- diag(c(lam1, lam2))
  Sigma <- Pmat %*% Dmat %*% t(Pmat)
  kernel.local[, ,i] <-  Sigma
}

Sigma.mat <- matrix(rep(NA, (n*t)^2), nrow = n*t)
Q.mat <- matrix(rep(NA, (n*t)^2), nrow = n*t)
Inv_ij <- matrix(rep(NA,4),2,2)

for (i in 1:nrow(locations)) {
  #Sigma.mat[i, i] <- 1
  #Q.mat[i, i] <- 0
  Kernel_i <- kernel.local[, , i]
  det_i <- Kernel_i[1,1] * Kernel_i[2,2] - Kernel_i[1,2] * Kernel_i[2,1]
  for (j in i:nrow(locations)) {
    Kernel_j <- kernel.local[, , j]
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
    Sigma.mat[j, i] <- Sigma.mat[i, j]
    Q.mat[j, i] <- Q.mat[i, j]
  }
}

Sigma.mat1 <- Sigma.mat

for(aa in 1:5){
  Sigma.mat <- Sigma.mat1[1:121,1:121+(aa-1)*121]
  h <- Q.mat[1:121,1:121+(aa-1)*121]
  
  kk=1             
  sim_cov[[aa]] <- ifelse(h!=0,Sigma.mat*(h/beta)^nu[kk]*besselK(h/beta,nu[kk])/(2^(nu[kk]-1)*gamma(nu[kk])),Sigma.mat)
  kk=2
  sim_cov2[[aa]] <- ifelse(h!=0,Sigma.mat*(h/beta)^nu[kk]*besselK(h/beta,nu[kk])/(2^(nu[kk]-1)*gamma(nu[kk])),Sigma.mat)
  
  rho=rot*(gamma(nu1+3/2)/gamma(nu1))^(1/2) * (gamma(nu2+3/2)/gamma(nu2))^(1/2)*gamma(nu3)/(gamma(nu3+3/2))
  sim_cov3[[aa]] <- ifelse(h!=0,(h/beta)^nu3 * besselK(h/beta,nu3)/(2^(nu3-1)*gamma(nu3))*sqrt(Sigma.mat* Sigma.mat)*rho,sqrt(Sigma.mat * Sigma.mat)*rho)
}

covariance_plots_functions('nonstationary_fig1a.pdf',sim_cov,sim_cov2,sim_cov3)

##############################################################

source(file=paste(root,"/matern_cov.R",sep=''))

matern_theta_strong=c(0.5,1.5,5,0.8)
thets <- matern_theta_strong
materncov <- matern_cov_regular_grid_v2(thets,c(.1001,.1001),time=5)

set.seed(12) #12345
A <- mvrnorm(1, rep(0,dim(materncov)[1]),materncov)

var1_sim <- var2_sim <- matrix(,ncol=5,nrow=n)
for(aa in 1:5){
  var1_sim[,aa] <- A[((aa-1)*n+1):(n*aa)]-mean(A[((aa-1)*n+1):(n*aa)])
  var2_sim[,aa] <- A[(((aa-1)*n+1):(n*aa))+n*5]-mean(A[(((aa-1)*n+1):(n*aa))+n*5])
}

simulation_plots_functions('fig2a.pdf',var1_sim,var2_sim)

##############################################################

log.lam1.1<-theta2(n,-3,1,1,-6,-7)
log.lam1.2<-theta2(n,-5,1,1,2,-12)
logit.phi.1<-theta2(n,0,1,-2,0,1)

log.lam2.1<-theta2(n,-3,-1,-1,-6,-7)
log.lam2.2<-theta2(n,-5,-1,-1,2,-12)
logit.phi.2<-theta2(n,0,-1,-2,0,1)


source(file=paste(root,"/covariance_plots_functions.R",sep=''))

lam1.1 <- list()
lam1.2 <- list()
phi.1 <- list()

lam2.1 <- list()
lam2.2 <- list()
phi.2 <- list()

for(aa in 1:5){
  lam1.1[[aa]] <- matrix(log.lam1.1[1:121+(aa-1)*121,],ncol=sqrt(n),nrow=sqrt(n))
  lam1.2[[aa]] <- matrix(log.lam1.2[1:121+(aa-1)*121,],ncol=sqrt(n),nrow=sqrt(n))
  phi.1[[aa]] <- matrix(logit.phi.1[1:121+(aa-1)*121,],ncol=sqrt(n),nrow=sqrt(n))
  lam2.1[[aa]] <- matrix(log.lam2.1[1:121+(aa-1)*121,],ncol=sqrt(n),nrow=sqrt(n))
  lam2.2[[aa]] <- matrix(log.lam2.2[1:121+(aa-1)*121,],ncol=sqrt(n),nrow=sqrt(n))
  phi.2[[aa]] <- matrix(logit.phi.2[1:121+(aa-1)*121,],ncol=sqrt(n),nrow=sqrt(n))
}
#range(exp(lam2.1[[1]]),exp(lam2.1[[2]]),exp(lam2.1[[3]]),exp(lam2.1[[4]]),exp(lam2.1[[5]]))
sigma_plots_functions('params_var1.pdf',lam1.1,lam1.2,phi.1,1,lam2.1,lam2.2,phi.2)
sigma_plots_functions('params_var2.pdf',lam2.1,lam2.2,phi.2,2,lam1.1,lam1.2,phi.1)

###############    ESTIMATED  ################################

# -1.50380504  -0.07558557   4.82244196  -5.66094136   2.99578278  43.99618927
#  19.24018098 -14.63719610  -0.51861398  15.23767899   2.68615214  -3.00711322

##############################################################

source(file=paste(root,"/matern_cov.R",sep=''))

matern_theta_strong=c(0.5,1.5,2,0.8)
thets <- matern_theta_strong
materncov <- matern_cov_regular_grid_v4(thets,c(.05001,.05001),time=5)
#isSymmetric(materncov)
image.plot(materncov)

set.seed(12345678) #12, 1234, 112
A <- mvrnorm(1, rep(0,dim(materncov)[1]),materncov)

var1_sim <- var2_sim <- matrix(,ncol=5,nrow=n)
for(aa in 1:5){
  var1_sim[,aa] <- A[((aa-1)*n+1):(n*aa)]-mean(A[((aa-1)*n+1):(n*aa)])
  var2_sim[,aa] <- A[(((aa-1)*n+1):(n*aa))+n*5]-mean(A[(((aa-1)*n+1):(n*aa))+n*5])
}

simulation_plots_functions('fig2a.pdf',var1_sim,var2_sim)

################################################################
log.lam1.1<-theta2(n,-3,0,0,-6,-7)
log.lam1.2<-theta2(n,-5,0,0,6,-4)
logit.phi.1<-theta2(n,0,1,-2,0,1)

log.lam2.1<-theta2(n,3,0,0,0,-1)
log.lam2.2<-theta2(n,-5,0,0,3,-1)
logit.phi.2<-theta2(n,0,1,-2,1,1)

range(pi/2 * exp(logit.phi.2)/(1 + exp(logit.phi.2)))

#log.lam2.1<-theta2(n,-1,0,0,-1,-1)
#log.lam2.2<-theta2(n,-5,0,0,3,-1)
#logit.phi.2<-theta2(n,0,1,-2,0,1)

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
#KERNEL_LIST[[2]] <-KERNEL_LIST[[1]]
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
    for (j in i:nrow(locations)) {
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
      Sigma.mat[j, i] <- Sigma.mat[i, j]
      Q.mat[j, i] <- Q.mat[i, j]
    }
  }
  FIN_Sigma.mat[[KK]] <- Sigma.mat
  dist0[[KK]] <- Q.mat
}

par(mfrow=c(2,5))
for(aa in 1:5){
  image.plot(FIN_Sigma.mat[[1]][1:121,1:121+(aa-1)*121])
}
for(aa in 1:5){
  image.plot(FIN_Sigma.mat[[2]][1:121,1:121+(aa-1)*121])
}
