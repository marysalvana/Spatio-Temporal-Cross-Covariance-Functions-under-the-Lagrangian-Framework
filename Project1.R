
WORKSTATION = 1

# 0.1 SET WORKING DIRECTORY WHERE DATA IS SAVED

if(WORKSTATION == 1){
  root <- '/Volumes/GoogleDrive/My Drive/Phd_Dissertation/ideal-happiness/'
}

source(file=paste(root,"Functions/empirical_spacetime_covariance.R",sep=''))

setwd(paste(root,'Figures',sep=''))

############################################################################
############################################################################
###################           REAL DATA ANALYSIS         ###################
############################################################################
############################################################################

load(paste(root,'Data/UTM_simulated_locations_only.RData',sep=''))
load(paste(root,'Data/simulated_locations_only.RData',sep=''))

# You can use raw netcdf data and pre-process it using preprocessing.R or you load this already preprocessed data

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

var1_cov <- empirical_st_cov(data1 = t(ave_u[1:80,]), cross = F, locations = grid_locations_UTM[pts[1:80],], max_time_lag = 5)  
var2_cov <- empirical_st_cov(data1 = t(ave_v[1:80,]), cross = F, locations = grid_locations_UTM[pts[1:80],], max_time_lag = 5)  
cross <- empirical_st_cov(data1 = t(ave_u[1:80,]), data2 = t(ave_v[1:80,]), cross = T, locations = grid_locations_UTM[pts[1:80],], max_time_lag = 5)  

binned <- empirical_covariance_dataframe(data1_cov = var1_cov, data2_cov = var2_cov, cross_cov = cross)

binned_orig2_0 <- binned_orig[which(binned_orig[,3]==0),]

binned_orig2_orig <- binned_orig[which(binned_orig[,3]<=1),]

binned_orig2_0 <- binned_orig[which(binned_orig[,3]==0),]

binned_orig2 <- binned_orig[which(binned_orig[,3]==1),]

plot(dist0,binned0[,4],pch=3,ylab='', col=1,xlab='Spatial Lag (km)',
     main='', col.main= "#4EC1DE",ylim=c(0,1))
