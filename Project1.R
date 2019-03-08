
WORKSTATION = 1

# 0.1 SET WORKING DIRECTORY WHERE DATA IS SAVED

if(WORKSTATION == 1){
  root <- '/Volumes/GoogleDrive/My Drive/Phd_Dissertation/ideal-happiness/'
}

source(file=paste(root,"Functions/load_packages.R",sep=''))
source(file=paste(root,"Functions/empirical_spacetime_covariance.R",sep=''))
source(file=paste(root,"Functions/data_format.R",sep=''))
source(file=paste(root,"Functions/cov_func.R",sep=''))
source(file=paste(root,"Functions/toeplitz_mat.R",sep=''))
source(file=paste(root,"Functions/wls_objective_functions.R",sep=''))

setwd(paste(root,'Figures',sep=''))

load(paste(root,'Data/UTM_simulated_locations_only.RData',sep=''))
load(paste(root,'Data/simulated_locations_only.RData',sep=''))

############################################################################
############################################################################
###################     SIMULATED DATA ANALYSIS       ######################
############################################################################
############################################################################

time <- 2
insample_loc_index <- 1:80
outsample_loc_index <- 81:100 #location indices for validation

to_remove <- c()
for(aa in 1:length(outsample_loc_index)){
  to_remove <- c(to_remove,seq(outsample_loc_index[aa],time*100,by=100),seq(outsample_loc_index[aa],time*100,by=100)+time*100)
}

#------------------------------------------

kappa <- matrix(c(0, 0, 0, 0), ncol=2, byrow=T)
thets <- c(3.218, 3.736, 794.8, 0.59, 1, 1) #set parameters to empirical estimates of real dataset
thets_lmc <- c(3, 4, 281.6, 700, 1, 1, 0.838, 0.545, 0, 0.999)

m1 <- simulation_study(true_param_spatial = thets, true_param_velocity = c(-497426.7319, -39634.6099, 100, 0.00009, 100), sim_model = 'matern', rand_vel = T, num_sim = 1, max_u = 1, num_variables = 2, location = grid_locations_UTM[pts,])
m2 <- simulation_study(true_param_spatial = thets, true_param_velocity = c(-497426.7319, -39634.6099), sim_model = 'matern', rand_vel = F, num_sim = 1, max_u = 1, num_variables = 2, location = grid_locations_UTM[pts,])
m3 <- simulation_study(true_param_spatial = thets_lmc, true_param_velocity = c(-2600000, -661200, 3604000, 1947000, 100, 0, 100, 100, 0, 100), sim_model = 'lmc', rand_vel = T, num_sim = 1, max_u = 1, num_variables = 2, location = grid_locations_UTM[pts,])
m4 <- simulation_study(true_param_spatial = thets_lmc, true_param_velocity = c(-2600000, -661200, 3604000, 1947000), sim_model = 'lmc', rand_vel = F, num_sim = 1, max_u = 1, num_variables = 2, location = grid_locations_UTM[pts,])

############################################################################
############################################################################
###################           REAL DATA ANALYSIS         ###################
############################################################################
############################################################################

# You can use raw netcdf data and pre-process it using preprocessing.R or you load this already preprocessed data

#we use only 25 years of data: 1980-2004 of January

nyears = 2004 - 1980 + 1
ndays = 31
nhours = 8

data_matrix <- data_format_into_matrix(data1 = u[pts,], data2 = v[pts,], temporal_replicates = nyears*ndays*nhours, simulated = F)

var1_cov <- empirical_st_cov(data1 = t(data_matrix[insample_loc_index, 1:(ncol(data_matrix)/2)]), cross = F, locations = grid_locations_UTM[pts[insample_loc_index],], max_time_lag = 5, simulated = F)  
var2_cov <- empirical_st_cov(data1 = t(data_matrix[insample_loc_index, (ncol(data_matrix)/2 + 1):ncol(data_matrix)]), cross = F, locations = grid_locations_UTM[pts[insample_loc_index],], max_time_lag = 5, simulated = F)  
cross <- empirical_st_cov(data1 = t(data_matrix[insample_loc_index, 1:(ncol(data_matrix)/2)]), data2 = t(data_matrix[insample_loc_index, (ncol(data_matrix)/2 + 1):ncol(data_matrix)]), cross = T, locations = grid_locations_UTM[pts[insample_loc_index],], max_time_lag = 5, simulated = F)  

binned <- empirical_covariance_dataframe(data1_cov = var1_cov, data2_cov = var2_cov, cross_cov = cross, simulated = F)

hlag <- sqrt(binned[which(binned[,3] == 0), 1]^2 + binned[which(binned[,3] == 0), 2]^2)
#display plots
par(mfrow = c(1,3))
plot(hlag/1000, binned[which(binned[,3]==0), 4], pch=3, ylab='', col=1,xlab='Spatial Lag (km)', main='', col.main= "#4EC1DE", ylim=c(0,1))
plot(hlag/1000, binned[which(binned[,3]==0), 5], pch=3, ylab='', col=1,xlab='Spatial Lag (km)', main='', col.main= "#4EC1DE", ylim=c(0,1))
plot(hlag/1000, binned[which(binned[,3]==0), 6], pch=3, ylab='', col=1,xlab='Spatial Lag (km)', main='', col.main= "#4EC1DE", ylim=c(0,1))

rotation_matrix <- matrix(c(-1.1339125*cos(116.9811232), -1.1339125*sin(116.9811232),
                            0.8617141*sin(116.9811232), -0.8617141*cos(116.9811232)), ncol = 2, byrow = T)
