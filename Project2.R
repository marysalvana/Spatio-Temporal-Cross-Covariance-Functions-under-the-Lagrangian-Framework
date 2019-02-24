
WORKSTATION = 1

# 0.1 SET WORKING DIRECTORY WHERE DATA IS SAVED

if(WORKSTATION == 1){
  root <- '/Volumes/GoogleDrive/My Drive/Phd_Dissertation/ideal-happiness/'
}

source(file=paste(root,"Functions/load_packages.R",sep=''))
source(file=paste(root,"Functions/empirical_spacetime_covariance.R",sep=''))
source(file=paste(root,"Functions/data_format.R",sep=''))
source(file=paste(root,"Functions/NSconvo_fit_multi.R",sep=''))
source(file=paste(root,"Functions/matern_cov.R",sep=''))

setwd(paste(root,'Figures',sep=''))

############################################################################
############################################################################
###################     SIMULATED DATA ANALYSIS         ####################
############################################################################
############################################################################

n<-400
N<-sqrt(n)
t=1
grid_x = seq(from=0.5, to=1, length.out = N)
grid_y = seq(from=0.5, to=1, length.out = N)
sim_grid_locations = expand.grid(grid_x, grid_y)
sim_grid_locations <- matrix(cbind(sim_grid_locations[,1],sim_grid_locations[,2]),ncol=2)

##generate locations
matern_theta_strong=c(0.5,2,4,0.8)
thets <- matern_theta_strong
materncov <- matern_cov_regular_grid_v4(thets,c(.05001,.05001),time=t)

A <- mvrnorm(1000, rep(0,dim(materncov)[1]),materncov)

N.mc <- 4
fit.radius <- 0.29

locations <- coords <- as.matrix(sim_grid_locations)
N <- dim(coords)[1]
data <- matrix(A, nrow=2*N)

lon_min <- min(coords[,1])
lon_max <- max(coords[,1])
lat_min <- min(coords[,2])
lat_max <- max(coords[,2])

#=======================================
# mixture component knot locations
#=======================================
mc_x <- seq(from = lon_min + 0.5*(lon_max - lon_min)/floor(sqrt(N.mc)),
            to = lon_max - 0.5*(lon_max - lon_min)/floor(sqrt(N.mc)),
            length = floor(sqrt(N.mc)) )
mc_y <- seq(from = lat_min + 0.5*(lat_max - lat_min)/floor(sqrt(N.mc)),
            to = lat_max - 0.5*(lat_max - lat_min)/floor(sqrt(N.mc)),
            length = floor(sqrt(N.mc)) )
mc.locations <- expand.grid( mc_x, mc_y )
mc.locations <- matrix(c(mc.locations[,1], mc.locations[,2]), ncol=2, byrow=F)

K <- dim(mc.locations)[1]

#===========================================================================
# Check the mixture component locations
#===========================================================================
check.mc.locs <- mc_N( coords, mc.locations, fit.radius )
lambda.w <- ( 0.5*min(dist(mc.locations)) )^2

mean.model = data ~ 1
data <- A[1:400]
OLS.model1 <- lm(mean.model, x=TRUE)
data <- A[401:800]
OLS.model2 <- lm(mean.model, x=TRUE)
Xmat <- matrix( unname( OLS.model1$x ), nrow=2*N )

#===========================================================================
# Specify lower, upper, and initial parameter values for optim()
#===========================================================================

max.distance <- sqrt(sum((c(lon_min,lat_min) - c(lon_max,lat_max))^2))

resid.var <- (max(summary(OLS.model1)$sigma,summary(OLS.model2)$sigma))^2

coords_mc_dist <- StatMatch::mahalanobis.dist(data.x = coords, data.y = mc.locations, vc = diag(2))

k = 2
ind_local <- (coords_mc_dist[,k] <= fit.radius)

# Subset
temp.locations <- coords[ind_local,]
n.fit <- dim(temp.locations)[1]

if(!is.matrix(A)){
  temp.data <- A[rep(ind_local,2*t)]
}else{
  temp.data <- A[,rep(ind_local,2*t)]
}

Xtemp <- as.matrix(rep(1,length(which(ind_local==TRUE))),ncol=1)

#####################################################
# Local estimation

data_matrix = A[,rep(ind_local,2)]

var1_cov <- empirical_st_cov(data1 = data_matrix[,1:(ncol(data_matrix)/2)], cross = F, locations = temp.locations, max_time_lag = 0)  
var2_cov <- empirical_st_cov(data1 = data_matrix[,(ncol(data_matrix)/2+1):ncol(data_matrix)], cross = F, locations = temp.locations, max_time_lag = 0)  
cross <- empirical_st_cov(data1 = data_matrix[,1:(ncol(data_matrix)/2)], data2 = data_matrix[,(ncol(data_matrix)/2+1):ncol(data_matrix)], cross = T, locations = temp.locations, max_time_lag = 0)  

binned <- empirical_covariance_dataframe(data1_cov = var1_cov, data2_cov = var2_cov, cross_cov = cross)

hlag <- sqrt(binned[which(binned[,3]==0),1]^2 + binned[which(binned[,3]==0),2]^2)

#Check plots

plot(hlag, binned[which(binned[,3]==0),5], pch=3, ylab='', col=1,xlab='Spatial Lag (km)', main='', col.main= "#4EC1DE",ylim=c(0,1))



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
