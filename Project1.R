
WORKSTATION = 1

# 0.1 SET WORKING DIRECTORY WHERE DATA IS SAVED

if(WORKSTATION == 1){
  root <- '/Volumes/GoogleDrive/My Drive/Phd_Dissertation/ideal-happiness/'
}

source(file=paste(root,"Functions/load_packages.R",sep=''))
source(file=paste(root,"Functions/main_func.R",sep=''))
source(file=paste(root,"Functions/data_format.R",sep=''))
source(file=paste(root,"Functions/cov_func.R",sep=''))
source(file=paste(root,"Functions/toeplitz_mat.R",sep=''))
source(file=paste(root,"/legend.gradient2.R",sep=''))

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

dat1 <- dat2 <- matrix(, ncol = nyears*ndays*nhours, nrow = nrow(spdf))
for(aa in 1:(nyears*ndays*nhours)){
  
  data1 <- data.frame(grid_locations[pts,], u[pts, aa], v[pts, aa])
  colnames(data1) <- c('lon', 'lat', 'var1', 'var2')
  
  spdf <- SpatialPointsDataFrame(coords = data1[, c("lon", "lat")], data = data1,
                                 proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  X <- matrix(cbind(1, coordinates(spdf), coordinates(spdf)[,1]/coordinates(spdf)[,2]), nrow = nrow(spdf))
  coef1 <- solve(t(X) %*% X) %*% t(X) %*% spdf$var1
  coef2 <- solve(t(X) %*% X) %*% t(X) %*% spdf$var2
  res1 <- matrix(spdf$var1, ncol = 1) - X %*% coef1
  res2 <- matrix(spdf$var2, ncol = 1) - X %*% coef2
  
  dat1[,aa] <- (res1 - mean(res1))/sd(res1)
  dat2[,aa] <- (res2 - mean(res2))/sd(res2)
}

dat_analysis <- fitting(coordinates = grid_locations_UTM[pts,], obs1 = dat1, obs2 = dat2)

############################################################################
############################################################################
###################           APPENDIX         ###################
############################################################################
############################################################################

data_plots_functions('testdata.pdf', dat1, dat2)


