
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
source(file=paste(root,"/covariance_plots_functions.R",sep=''))


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
#nyears*ndays*nhours
dat1 <- dat2 <- matrix(, ncol = 5, nrow = nrow(grid_locations))
for(aa in 1:5){
  
  data1 <- data.frame(grid_locations_UTM/1000, u[, aa], v[, aa])
  colnames(data1) <- c('lon', 'lat', 'var1', 'var2')
  
  spdf <- SpatialPointsDataFrame(coords = data1[, c("lon", "lat")], data = data1,
                                 proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  #X <- matrix(cbind(1, coordinates(spdf), coordinates(spdf)[,1]/coordinates(spdf)[,2]), nrow = nrow(spdf))
  X <- matrix(cbind(1, coordinates(spdf)), nrow = nrow(spdf))
  coef1 <- solve(t(X) %*% X) %*% t(X) %*% spdf$var1
  coef2 <- solve(t(X) %*% X) %*% t(X) %*% spdf$var2
  res1 <- matrix(spdf$var1, ncol = 1) - X %*% coef1
  res2 <- matrix(spdf$var2, ncol = 1) - X %*% coef2
  
  dat1[,aa] <- (res1 - mean(res1))/sd(res1)
  dat2[,aa] <- (res2 - mean(res2))/sd(res2)
}

for(aa in 1:5){
  
  data1 <- data.frame(grid_locations, u[, aa], v[, aa])
  colnames(data1) <- c('lon', 'lat', 'var1', 'var2')
  
  spdf <- SpatialPointsDataFrame(coords = data1[, c("lon", "lat")], data = data1,
                                 proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  #X <- matrix(cbind(1, coordinates(spdf), coordinates(spdf)[,1]/coordinates(spdf)[,2]), nrow = nrow(spdf))
  X <- matrix(cbind(1, coordinates(spdf)), nrow = nrow(spdf))
  coef1 <- solve(t(X) %*% X) %*% t(X) %*% spdf$var1
  coef2 <- solve(t(X) %*% X) %*% t(X) %*% spdf$var2
  res1 <- matrix(spdf$var1, ncol = 1) - X %*% coef1
  res2 <- matrix(spdf$var2, ncol = 1) - X %*% coef2
  
  dat1[,aa] <- (res1 - mean(res1))/sd(res1)
  dat2[,aa] <- (res2 - mean(res2))/sd(res2)
}

dat1 <- dat2 <- matrix(, ncol = nyears*ndays*nhours, nrow = nrow(grid_locations[pts,]))
for(aa in 1:(nyears*ndays*nhours)){
  dat1[,aa] <- (u[pts,aa] - mean(u[pts,aa]))/sd(u[pts,aa])
  dat2[,aa] <- (v[pts,aa] - mean(v[pts,aa]))/sd(v[pts,aa])
}

dat1 <- dat2 <- matrix(, ncol = nyears*ndays*nhours, nrow = nrow(grid_locations))
for(aa in 1:(nyears*ndays*nhours)){
  dat1[,aa] <- (u[,aa] - mean(u[,aa]))/sd(u[,aa])
  dat2[,aa] <- (v[,aa] - mean(v[,aa]))/sd(v[,aa])
}

dat1 <- dat2 <- matrix(, ncol = nyears*ndays*nhours, nrow = nrow(grid_locations[pts,]))
for(aa in 1:(nyears*ndays*nhours)){
  dat1[,aa] <- u[pts,aa]
  dat2[,aa] <- v[pts,aa]
}

obs1 <- dat1
obs2 <- dat2

obs1 <- u[pts, 1:(nyears*ndays*nhours)]
obs2 <- v[pts, 1:(nyears*ndays*nhours)]

#dat1 <- dat2 <- matrix(, ncol = 5, nrow = nrow(spdf))
#for(aa in 1:5){
#  dat1[,aa] <- (u[pts,aa] - mean(u[pts,aa]))/sd(u[pts,aa])
#  dat2[,aa] <- (v[pts,aa] - mean(v[pts,aa]))/sd(v[pts,aa])
#}

dat_analysis <- fitting(coordinates = grid_locations_UTM[pts,], obs1 = dat1, obs2 = dat2)

############################################################################
############################################################################
###################           APPENDIX         ###################
############################################################################
############################################################################

data_plots_functions('testdata.pdf', dat1, dat2)


data1 <- data.frame(grid_locations_UTM[pts,] %*% R/1000000, u[pts, 1], v[pts, 1])
colnames(data1) <- c('lon', 'lat', 'var1', 'var2')
spdf <- SpatialPointsDataFrame(coords = data1[, c("lon", "lat")], data = data1,
                               proj4string = CRS("+proj=longlat +datum=WGS84"))

data1 <- cbind(grid_locations[pts,], 1, u[pts, 1], v[pts, 1])
for(bb in 2:6200){
  data1 <- rbind(data1, cbind(grid_locations[pts,], bb, u[pts, bb], v[pts, bb]))
}

#X <- matrix(cbind(1, coordinates(spdf), coordinates(spdf)[,1]/coordinates(spdf)[,2]), nrow = nrow(spdf))
X <- data1[,1:3]
coef1 <- solve(t(X) %*% X) %*% t(X) %*% data1[,4]
coef2 <- solve(t(X) %*% X) %*% t(X) %*% data1[,5]
res1 <- matrix(data1[,4], ncol = 1) - X %*% coef1
res2 <- matrix(data1[,5], ncol = 1) - X %*% coef2

dat1 <- dat2 <- matrix(, ncol = 5, nrow = nrow(grid_locations[pts,]))

for(aa in 1:5){
  dat1[,aa] <- (res1[(aa - 1)*100 + 1:100,] - mean(res1[(aa - 1)*100 + 1:100,]))/sd(res1[(aa - 1)*100 + 1:100,])
  dat2[,aa] <- (res2[(aa - 1)*100 + 1:100,] - mean(res2[(aa - 1)*100 + 1:100,]))/sd(res2[(aa - 1)*100 + 1:100,])
}

for(aa in 1:5){
  dat1[,aa] <- (u[pts,aa] - mean(u[pts,aa]))/sd(u[pts,aa])
  dat2[,aa] <- (v[pts,aa] - mean(v[pts,aa]))/sd(v[pts,aa])
}

colfunc <- colorRampPalette(c("blue","yellow","red"))

var1 = dat1
var2 = dat2

pdf(file=filename, width=20, height=9)

zr <- range(c(var1, var2))

split.screen( rbind(c(0.1,0.85,0,1), c(.95,0.99,0,0.98)))
screen(1)
layout(mat=matrix(1:10, nrow=2, ncol=5, byrow=T))

for(bb in 1:5){
  
  xx <- findInterval(var1[,bb], sort(c(var1, var2)))
  
  if(bb==1){
    par(pty="s") 
    par(mai=c(0.4,0.4,0.4,0.4))
    map("worldHires", xlim = c(range(grid_locations[,1])[1]-1,range(grid_locations[,1])[2]+2), ylim = c(range(grid_locations[,2])[1]-1, range(grid_locations[,2])[2]+2), lwd=0.5)
    plot(spdf, add=TRUE, pch=16, col=colfunc(length(c(var1, var2)))[xx], lwd=0.5)
    mtext(expression(Z[1]), side = 2, line = 0, adj = 0.5, cex = 1.5, font=3,col="#0086FF")
    mtext(paste("t=",bb,sep=""), side = 3, line = 1, adj = 0.5, cex = 1.5, font=2,col="#4EC1DE")
  }else{
    par(pty="s") 
    par(mai=c(0.4,0.4,0.4,0.4))
    map("worldHires", xlim = c(range(grid_locations[,1])[1]-1,range(grid_locations[,1])[2]+2), ylim = c(range(grid_locations[,2])[1]-1,range(grid_locations[,2])[2]+2), lwd=0.5)
    plot(spdf, add=TRUE, pch=16, col=colfunc(length(c(var1, var2)))[xx], lwd=0.5)
    mtext(paste("t=",bb,sep=""), side = 3, line = 1, adj = 0.5, cex = 1.5, font=2,col="#4EC1DE")
  }
  
}
for(bb in 1:5){
  xx <- findInterval(var2[,bb], sort(c(var1, var2)))
  
  if(bb==1){
    par(pty="s") 
    par(mai=c(0.4,0.4,0.4,0.4))
    map("worldHires", xlim = c(range(grid_locations[,1])[1]-1,range(grid_locations[,1])[2]+2), ylim = c(range(grid_locations[,2])[1]-1,range(grid_locations[,2])[2]+2), lwd=0.5)
    plot(spdf, add=TRUE, pch=16, col=colfunc(length(c(var1, var2)))[xx], lwd=0.5)
    mtext(expression(Z[2]), side = 2, line = 0, adj = 0.5, cex = 1.5, font=3,col="#0086FF")
  }else{
    par(pty="s") 
    par(mai=c(0.4,0.4,0.4,0.4))
    map("worldHires", xlim = c(range(grid_locations[,1])[1]-1,range(grid_locations[,1])[2]+2), ylim = c(range(grid_locations[,2])[1]-1,range(grid_locations[,2])[2]+2), lwd=0.5)
    plot(spdf, add=TRUE, pch=16, col=colfunc(length(c(var1, var2)))[xx], lwd=0.5)
    axis(1, at=seq(0,1,by=0.2), labels = seq(0,1,by=0.2))
  }
}
screen(2)
screen(2)
x1 <- c(0.3,0.4,0.4,0.3)
y1 <- c(0.35,0.35,0.6,0.6)
legend.gradient2(cbind(x1,y1),cols = colorRampPalette(colors) (colsteps), title = "", 
                 limits = seq(-3,3,length.out = 5), cex=1)
close.screen( all=TRUE)
dev.off()
