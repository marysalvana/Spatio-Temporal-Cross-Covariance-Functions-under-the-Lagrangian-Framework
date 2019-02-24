
WORKSTATION = 1

# 0.1 SET WORKING DIRECTORY WHERE DATA IS SAVED

if(WORKSTATION == 1){
  root <- '/Volumes/GoogleDrive/My Drive/Phd_Dissertation/ideal-happiness/'
}

source(file=paste(root,"Functions/load_packages.R",sep=''))
source(file=paste(root,"Functions/empirical_spacetime_covariance.R",sep=''))
source(file=paste(root,"Functions/data_format.R",sep=''))

setwd(paste(root,'Figures',sep=''))

############################################################################
############################################################################
###################     SIMULATED DATA ANALYSIS         ####################
############################################################################
############################################################################



############################################################################
############################################################################
###################           REAL DATA ANALYSIS         ###################
############################################################################
############################################################################

load(paste(root,'Data/UTM_simulated_locations_only.RData',sep=''))
load(paste(root,'Data/simulated_locations_only.RData',sep=''))

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
