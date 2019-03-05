library(mvnfast)
library(dplyr)
# 0.1 SET WORKING DIRECTORY WHERE DATA IS SAVED

root <- '/scratch/dragon/intel/salvanmo/'

source(file=paste(root,"Functions/empirical_st_cov.R",sep=''))
source(file=paste(root,"Functions/wls_objective_functions.R",sep=''))

load(paste(root,'DATA/UTM_simulated_locations_only.RData',sep=''))
load(paste(root,'DATA/simulated_locations_only.RData',sep=''))

R <- matrix(c(-1.1339125 *cos(116.9811232),-1.1339125 *sin(116.9811232),0.8617141*sin(116.9811232),
              -0.8617141*cos(116.9811232)),ncol=2,byrow=T)

aggs_max <- 1 #the number of consecutive timesteps you want to average. Here we take the data as is and do not take averages since they are already Gaussian and stationary

#we use only 25 years of data: 1980-2004.

ave_u <- ave_v <- matrix(,ncol=floor(25*31*8/aggs_max),nrow=100)
for(locat in 1:100){
  new_uu <- new_uu2 <- matrix(,ncol=floor(25*31*8/aggs_max),nrow=aggs_max)
  
  for(gg in 1:floor(25*31*8/aggs_max)){
    new_uu[,gg] <- u[pts[locat],((gg-1)*aggs_max+1):(gg*aggs_max)]
    new_uu2[,gg] <- v[pts[locat],((gg-1)*aggs_max+1):(gg*aggs_max)]
  }
  ave_u[locat,]<- colMeans(new_uu)
  ave_v[locat,]<- colMeans(new_uu2)
}

var1_cov <- empirical_st_cov_UTM(t(ave_u[1:80,]))  ##this is coded for u=0:5
var2_cov <- empirical_st_cov_UTM(t(ave_v[1:80,]))
cross <- empirical_st_cov_cross_UTM(t(ave_u[1:80,]),t(ave_v[1:80,]))

empirical_var1 <- var1_cov[[1]]
colnames(empirical_var1) <- c('x','y','var1','x_old','y_old')
empirical_var2 <- var2_cov[[1]]
colnames(empirical_var2) <- c('x','y','var2','x_old','y_old')
empirical_var3 <- cross[[1]]
colnames(empirical_var3) <- c('x','y','var3','x_old','y_old')

binned.1 <- empirical_var1 %>% group_by(x,y) %>% summarize(avg1=mean(var1))
binned.2 <- empirical_var2 %>% group_by(x,y) %>% summarize(avg1=mean(var2))
binned.3 <- empirical_var3 %>% group_by(x,y) %>% summarize(avg1=mean(var3))

binned_orig <- cbind(binned.1$x,binned.1$y,rep(0,nrow(binned.1)),binned.1$avg1,
                     binned.2$avg1,binned.3$avg1)

dist0 <- sqrt(binned_orig[,1]^2+binned_orig[,2]^2)/1000

for (i in 2:6){
  empirical_var1 <- var1_cov[[i]]
  colnames(empirical_var1) <- c('x','y','var1','x_old','y_old')
  empirical_var2 <- var2_cov[[i]]
  colnames(empirical_var2) <- c('x','y','var2','x_old','y_old')
  empirical_var3 <- cross[[i]]
  colnames(empirical_var3) <- c('x','y','var3','x_old','y_old')
  
  binned.1 <- empirical_var1 %>% group_by(x,y) %>% summarize(avg1=mean(var1))
  binned.2 <- empirical_var2 %>% group_by(x,y) %>% summarize(avg1=mean(var2))
  binned.3 <- empirical_var3 %>% group_by(x,y) %>% summarize(avg1=mean(var3))
  
  binned_orig <- rbind(binned_orig,cbind(binned.1$x,binned.1$y,rep(i-1,nrow(binned.1)),binned.1$avg1,
                                         binned.2$avg1,binned.3$avg1))
}

binned0 <- cbind(dist0,binned_orig[which(binned_orig[,3]==0),4],
                 binned_orig[which(binned_orig[,3]==0),5],
                 binned_orig[which(binned_orig[,3]==0),6],
                 binned_orig[which(binned_orig[,3]==0),1],
                 binned_orig[which(binned_orig[,3]==0),2])

binned_orig2_0 <- binned_orig[which(binned_orig[,3]==0),]

binned_orig2_orig <- binned_orig[which(binned_orig[,3]<=1),]

binned_orig2_0 <- binned_orig[which(binned_orig[,3]==0),]

binned_orig2 <- binned_orig[which(binned_orig[,3]==1),]

mod1_parms <- c(3.218,3.736,794.8,0.59,704400,205700)
#theta <- c(-5969,-75.26,1,0.5,1)
theta <- c(-4147.2203286,-20.9478018,2.2028818,0.7302455,1.3907865)
start_time <- Sys.time()

fit.nsst<-optim(theta, m1, emp_cov1=binned_orig2, control=list(maxit=3000,parscale=theta,trace=5))

end_time <- Sys.time()

end_time <- end_time - start_time

mod1_params_motion <- fit.nsst$par
mod1_params_value <- fit.nsst$value

save(mod1_params_motion,end_time,mod1_params_value, file="/scratch/dragon/intel/salvanmo/Results/m1_v2.Rdata")

