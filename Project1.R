
WORKSTATION = 1

# 0.1 SET WORKING DIRECTORY WHERE DATA IS SAVED

if(WORKSTATION == 1){
  root <- '/Volumes/GoogleDrive/My Drive/Phd_Dissertation/ideal-happiness/'
}

source(file=paste(root,"Functions/load_packages.R",sep=''))
source(file=paste(root,"Functions/empirical_spacetime_covariance.R",sep=''))
source(file=paste(root,"Functions/data_format.R",sep=''))

setwd(paste(root,'Figures',sep=''))

load(paste(root,'Data/UTM_simulated_locations_only.RData',sep=''))
load(paste(root,'Data/simulated_locations_only.RData',sep=''))

############################################################################
############################################################################
###################     SIMULATED DATA ANALYSIS         ####################
############################################################################
############################################################################

outsample_loc_index <- 81:100 #location indices for validation

to_remove <- c()
for(aa in 1:length(outsample_loc_index)){
  to_remove <- c(to_remove,seq(outsample_loc_index[aa],time*100,by=100),seq(outsample_loc_index[aa],time*100,by=100)+time*100)
}

#------------------------------------------

kappa <- matrix(c(0,0,0,0),ncol=2,byrow=T)
thets <- c(3.218,3.736,(3.218+3.736)/2,rep(794.8,3),0,0,1,1,0.59)

materncov <- matern_cov_mod1_UTM(thets,c(-497426.7319,-39634.6099),time=2)

mod2_params <- matrix(,ncol=6,nrow=100)
for(samp in 1:100){
  A <- mvrnorm(n=1000,mu=rep(0,dim(materncov)[1]),Sigma=materncov)
  A1 <- A[,-to_remove]
  
  empcor <- cor(A1)
  
  empcov <- list()
  empcov[[1]] <- empcor[1:nrow(new_grid_locations),1:nrow(new_grid_locations)]
  empcov[[2]] <- empcor[(2*nrow(new_grid_locations)+1):(3*nrow(new_grid_locations)),(2*nrow(new_grid_locations)+1):(3*nrow(new_grid_locations))]
  empcov[[3]] <- empcor[1:nrow(new_grid_locations),(2*nrow(new_grid_locations)+1):(3*nrow(new_grid_locations))]
  
  index <- index1 <- index2 <- which(dist0 <= max(dist0),arr.ind = T)
  
  empirical <- data.frame(dist0[index],empcov[[1]][index],
                          empcov[[2]][index],empcov[[3]][index])
  colnames(empirical) <- c('dist','var1','var2','var1.2')
  
  binned.1 <- empirical %>% group_by(dist) %>% summarize(avg1=mean(var1),
                                                         avg2=mean(var2),
                                                         avg3=mean(var1.2))
  binned <- cbind(binned.1$dist,binned.1$avg1,binned.1$avg2,binned.1$avg3)
  binned0 <- binned
  
  empcov <- list()
  empcov[[1]] <- empcor[1:(time*nrow(new_grid_locations)),1:(time*nrow(new_grid_locations))]
  empcov[[2]] <- empcor[(time*nrow(new_grid_locations)+1):(2*time*nrow(new_grid_locations)),(time*nrow(new_grid_locations)+1):(2*time*nrow(new_grid_locations))]
  empcov[[3]] <- empcor[1:(time*nrow(new_grid_locations)),(time*nrow(new_grid_locations)+1):(2*time*nrow(new_grid_locations))]
  
  index1 <- which(dist0 <= max(dist0),arr.ind = T)
  lags <- lags1 <- new_grid_locations[index1[,2],]-new_grid_locations[index1[,1],]
  index <- cbind(index1[,1],index1[,2]+nrow(new_grid_locations))
  
  empirical <- data.frame(lags[,1],lags[,2],rep(seq(1,1),each=nrow(lags1)),empcov[[1]][index],
                          empcov[[2]][index],empcov[[3]][index])
  colnames(empirical) <- c('x','y','t','var1','var2','var1.2')
  
  binned.1 <- empirical %>% group_by(x,y,t) %>% summarize(avg1=mean(var1),
                                                          avg2=mean(var2),
                                                          avg3=mean(var1.2))
  binned <- cbind(binned.1$x,binned.1$y,binned.1$t,binned.1$avg1,binned.1$avg2,binned.1$avg3)
  
  binned_orig2 <- binned
  
  theta <- c(3.218,3.736,794.8)
  fit.nsst<-optim(theta, wls2_for_sim, emp_cov1=binned0, weights=3,control=list(maxit=3000,parscale=theta,trace=5))
  mod2_params.temp <- tempo_theta <- fit.nsst$par
  
  theta <- 0.5
  fit.nsst<-optim(theta, wls3_for_sim, emp_cov1=binned0, weights=3,method='SANN', control=list(maxit=3000,parscale=theta,trace=5))
  mod2_params.temp2 <-c(mod2_params.temp,fit.nsst$par)
  tempo_theta <- mod2_params.temp2
  
  theta <- c(-500400,-59740)
  fit.nsst<-optim(theta, m2_no_decay, emp_cov1=binned_orig2, weights=3, control=list(maxit=3000,parscale=theta,trace=5))
  mod2_params[samp,] <- c(mod2_params.temp2,fit.nsst$par/1000)
}


############################################################################
############################################################################
###################           REAL DATA ANALYSIS         ###################
############################################################################
############################################################################



# You can use raw netcdf data and pre-process it using preprocessing.R or you load this already preprocessed data

#we use only 25 years of data: 1980-2004 of January

insample_loc_index <- 1:80

nyears = 2004-1980+1
ndays = 31
nhours = 8

data_matrix <- data_format_into_matrix(data1 = u[pts,], data2 = v[pts,], temporal_replicates = nyears*ndays*nhours, simulated = F)

var1_cov <- empirical_st_cov(data1 = t(data_matrix[insample_loc_index,1:(ncol(data_matrix)/2)]), cross = F, locations = grid_locations_UTM[pts[insample_loc_index],], max_time_lag = 5)  
var2_cov <- empirical_st_cov(data1 = t(data_matrix[insample_loc_index,(ncol(data_matrix)/2+1):ncol(data_matrix)]), cross = F, locations = grid_locations_UTM[pts[insample_loc_index],], max_time_lag = 5)  
cross <- empirical_st_cov(data1 = t(data_matrix[insample_loc_index,1:(ncol(data_matrix)/2)]), data2 = t(data_matrix[insample_loc_index,(ncol(data_matrix)/2+1):ncol(data_matrix)]), cross = T, locations = grid_locations_UTM[pts[insample_loc_index],], max_time_lag = 5)  

binned <- empirical_covariance_dataframe(data1_cov = var1_cov, data2_cov = var2_cov, cross_cov = cross)

hlag <- sqrt(binned[which(binned[,3]==0),1]^2 + binned[which(binned[,3]==0),2]^2)

#Check plots

plot(hlag/1000, binned[which(binned[,3]==0),6], pch=3, ylab='', col=1,xlab='Spatial Lag (km)', main='', col.main= "#4EC1DE",ylim=c(0,1))
