
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

time <- 5
insample_loc_index <- 1:80
outsample_loc_index <- 81:100 #location indices for validation

to_remove <- c()
for(aa in 1:length(outsample_loc_index)){
  to_remove <- c(to_remove,seq(outsample_loc_index[aa],time*100,by=100),seq(outsample_loc_index[aa],time*100,by=100)+time*100)
}

#------------------------------------------

kappa <- matrix(c(0, 0, 0, 0), ncol=2, byrow=T)
thets <- c(3.218, 3.736, 794.8, 0.59, 1, 1, 0, 0) #set parameters to empirical estimates of real dataset

materncov <- matern_cov(thets, wind = c(-497426.7319,-39634.6099), max_time_lag = time - 1, p = 2, locations = grid_locations_UTM[pts,])

mod2_params <- matrix(, ncol=8, nrow=100)

for(iter in 1:100){
  
  A <- mvrnorm(n = 1000, mu = rep(0, dim(materncov)[1]), Sigma = materncov)
  A1 <- A[, -to_remove]
  
  conso_cor <- empirical_st_cov(data1 = A1, locations = grid_locations_UTM[pts[insample_loc_index],], max_time_lag = time - 1, simulated = T)  
  binned <- empirical_covariance_dataframe(data1_cov = conso_cor, simulated = T)
  
  if(iter == 1){
    hlag <- sqrt(binned[which(binned[,3]==0),1]^2 + binned[which(binned[,3]==0),2]^2)
    #display plots
    par(mfrow=c(1,3))
    plot(hlag/1000, binned[which(binned[,3]==0), 4], pch=3, ylab='', col=1,xlab='Spatial Lag (km)', main='', col.main= "#4EC1DE",ylim=c(0,1))
    plot(hlag/1000, binned[which(binned[,3]==0), 5], pch=3, ylab='', col=1,xlab='Spatial Lag (km)', main='', col.main= "#4EC1DE",ylim=c(0,1))
    plot(hlag/1000, binned[which(binned[,3]==0), 6], pch=3, ylab='', col=1,xlab='Spatial Lag (km)', main='', col.main= "#4EC1DE",ylim=c(0,1))
  }
  
  theta_init <- c(3.218, 3.736, 794.8, 1, 1, 0.5) #change this values to empirical
  w_init <- c(-500400, -59740)
  mod2 <- fit_model(init = theta_init, wind_init = w_init, mod = 2, weight = 3, empcov_spatial = binned[which(binned[,3]==0),], empcov_st = binned[which(binned[,3] > 0),], nug_eff = F, meters = T, num_iter = 10)
  mod2_params[iter, ] <- mod2$parameters
}


############################################################################
############################################################################
###################           REAL DATA ANALYSIS         ###################
############################################################################
############################################################################



# You can use raw netcdf data and pre-process it using preprocessing.R or you load this already preprocessed data

#we use only 25 years of data: 1980-2004 of January

nyears = 2004-1980+1
ndays = 31
nhours = 8

data_matrix <- data_format_into_matrix(data1 = u[pts,], data2 = v[pts,], temporal_replicates = nyears*ndays*nhours, simulated = F)

var1_cov <- empirical_st_cov(data1 = t(data_matrix[insample_loc_index,1:(ncol(data_matrix)/2)]), cross = F, locations = grid_locations_UTM[pts[insample_loc_index],], max_time_lag = 5, simulated = F)  
var2_cov <- empirical_st_cov(data1 = t(data_matrix[insample_loc_index,(ncol(data_matrix)/2+1):ncol(data_matrix)]), cross = F, locations = grid_locations_UTM[pts[insample_loc_index],], max_time_lag = 5, simulated = F)  
cross <- empirical_st_cov(data1 = t(data_matrix[insample_loc_index,1:(ncol(data_matrix)/2)]), data2 = t(data_matrix[insample_loc_index,(ncol(data_matrix)/2+1):ncol(data_matrix)]), cross = T, locations = grid_locations_UTM[pts[insample_loc_index],], max_time_lag = 5, simulated = F)  

binned <- empirical_covariance_dataframe(data1_cov = var1_cov, data2_cov = var2_cov, cross_cov = cross, simulated = F)

hlag <- sqrt(binned[which(binned[,3]==0),1]^2 + binned[which(binned[,3]==0),2]^2)

#Check plots

plot(hlag/1000, binned[which(binned[,3]==0),6], pch=3, ylab='', col=1,xlab='Spatial Lag (km)', main='', col.main= "#4EC1DE",ylim=c(0,1))
