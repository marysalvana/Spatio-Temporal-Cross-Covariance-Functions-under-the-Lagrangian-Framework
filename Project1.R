
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

load(paste(root,'Data/UTM_simulated_locations_only.RData',sep='')) #contains u, v, pts
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
  to_remove <- c(to_remove,seq(outsample_loc_index[aa], time*100,by=100),seq(outsample_loc_index[aa],time*100,by=100)+time*100)
}

time <- 1

to_remove_old <- c()
for(aa in 1:length(outsample_loc_index)){
  to_remove_old <- c(to_remove_old,seq(outsample_loc_index[aa], time*100,by=100),seq(outsample_loc_index[aa],time*100,by=100)+time*100)
}
#------------------------------------------------------#

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

nyears <- 2004 - 1980 + 1
ndays <- 31
nhours <- 8

var1 <- dat1 <- u[pts,1:(nyears*ndays*nhours)]
var2 <- dat2 <- v[pts,1:(nyears*ndays*nhours)]
loc <- grid_locations_UTM[pts,]

dat_analysis <- fitting(coordinates = loc, obs1 = dat1, obs2 = dat2)

#--M1

  #$parameters
  # -1.019976e+06 -7.893983e+04  7.272666e+09 -1.496189e+09  1.049021e+09

  #$fn_value
  # 1440.368

#--M2

  #$parameters
  #  3.435771e+00  4.040689e+00  7.614931e+02  1.000070e+00  9.991663e-01  5.898672e-01
  # -2.417319e+05  7.370889e+04  7.940310e-01  4.693621e-01

  #$fn_value
  # 1463.647

####################################################

WORKSTATION = 4

# 0.1 SET WORKING DIRECTORY WHERE DATA IS SAVED

if(WORKSTATION == 1){
  root <- '/Volumes/GoogleDrive/My Drive/Summer 2017_2018/Project 1/'
}else if (WORKSTATION == 2){
  root <- 'G://My Drive/Summer 2017_2018/Project 1/'
}else if (WORKSTATION == 3){
  root <- '/Users/salvanmo/Dropbox/05SalvanaMary/Spatio-Temporal Cross-Covariance Functions under the Lagrangian Framework/'
}else{
  root <- '/Users/salvanmo/silver-octo-spork/'
}

source(file=paste(root,"/load_packages.R",sep=''))
source(file=paste(root,"/legend.gradient2.R",sep=''))
source(file=paste(root,"/empirical_st_cov.R",sep=''))
source(file=paste(root,"/wls_objective_functions.R",sep=''))
source(file=paste(root,"/covariance_plots_functions.R",sep=''))

load(paste(root,'DATA/UTM_simulated_locations_only.RData',sep=''))
load(paste(root,'DATA/simulated_locations_only.RData',sep=''))

## NOTE: Z1 is saved as u and Z1 is saved as v

setwd(paste(root,'Figures',sep=''))

#------------------------BINNING OF DATA------------------------#

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

plot(dist0,binned0[,4],pch=3,ylab='', col=1,xlab='Spatial Lag (km)',
     main='', col.main= "#4EC1DE",ylim=c(0,1))

##-----------------------fig5b-----------------------##

binned3 <- binned_orig

dt=3

colors=c("blue","yellow","red")
colsteps=100

binned <- list()
for(cc in 1:dt){
  binned[[cc]] <- binned3[binned3[,3]==cc-1,]
}

zr1 <- range(binned_orig[,4:6])

covariance_plots_functions_section5('fig5b.pdf',binned,F)

##-----------------------fig5c-----------------------##

R <- matrix(c(-1.1339125 *cos(116.9811232),-1.1339125 *sin(116.9811232),0.8617141*sin(116.9811232),
              -0.8617141*cos(116.9811232)),ncol=2,byrow=T)

covariance_plots_functions_section5('fig5c.pdf',binned,T)

##-----------------------fig6a-----------------------##

selected <- c(34,42,43)

binned <- list()

for(ss in 1:3){
  date_stamp = 1+selected[ss]
  
  to_include <- seq(date_stamp,6200,by= 248)
  to_include2 <- seq(date_stamp+1,6200,by= 248)
  
  var1_cov <- empirical_st_cov_cross_UTM_foradvection_v2(t(ave_u[1:80,to_include]),t(ave_u[1:80,to_include2]))
  var2_cov <- empirical_st_cov_cross_UTM_foradvection_v2(t(ave_v[1:80,to_include]),t(ave_v[1:80,to_include2]))
  cross <- empirical_st_cov_cross_UTM_foradvection_v2(t(ave_u[1:80,to_include]),t(ave_v[1:80,to_include2]))
  
  empirical_var1 <- var1_cov
  colnames(empirical_var1) <- c('x','y','var1','x_old','y_old')
  empirical_var2 <- var2_cov
  colnames(empirical_var2) <- c('x','y','var2','x_old','y_old')
  empirical_var3 <- cross
  colnames(empirical_var3) <- c('x','y','var3','x_old','y_old')
  
  binned.1 <- empirical_var1 %>% group_by(x,y) %>% summarize(avg1=mean(var1))
  binned.2 <- empirical_var2 %>% group_by(x,y) %>% summarize(avg1=mean(var2))
  binned.3 <- empirical_var3 %>% group_by(x,y) %>% summarize(avg1=mean(var3))
  
  binned_orig <- cbind(binned.1$x,binned.1$y,rep(1,nrow(binned.1)),binned.1$avg1,
                       binned.2$avg1,binned.3$avg1)
  binned[[ss]] <- binned_orig
}

colors=c("blue","yellow","red")

colsteps=100

new_kappa_coords <- matrix(c(0,-0.7,1.7,-0.7,-0.2,1.8),ncol=2,byrow=T)
transfo <- (binned[[1]][,1:2]%*%R)/1000000

zr1 <- range(min(binned[[1]][,4:6],binned[[2]][,4:6],binned[[3]][,4:6]),1)
dt=3

pdf(paste('fig6a.pdf',sep=''), width=8.7, height=8)

split.screen( rbind(c(0.1,0.93,0.1,1), c(.93,0.99,0.1,1)))
split.screen( figs = c( 3, 1 ), screen = 1 )
split.screen( figs = c( 1, dt ), screen = 3 )
split.screen( figs = c( 1, dt ), screen = 4 )
split.screen( figs = c( 1, dt ), screen = 5 )

for(bb in 1:dt){
  screen(bb+5)
  if(bb==1){
    par(pty="s") 
    par(mai=c(0.1,0.1,0.2,0.1))
    plot(transfo, col=colorRampPalette(colors) (colsteps) [ findInterval(binned[[bb]][,4], seq(-1,1, length.out=colsteps)) ]
         , pch = 20, xaxt="n", cex.axis=0.9,
         main=paste('sub-period ',bb,sep=' '),xlab='', ylab='',col.main= "#4EC1DE")
    mtext(expression(hat(C)[11]), side = 2, line = 3, adj = 0.5, font=3,col="#0086FF",cex=1.3)
    #mtext('subperiod 1', side = 2, line = 3.5, adj = 0.5, font=3,col="#5970AF",cex=0.8)
    arrows(0, 0, transfo[which.max(binned[[bb]][,4]),1], transfo[which.max(binned[[bb]][,4]),2],lwd=2)
  }else{
    par(pty="s") 
    par(mai=c(0.1,0.1,0.2,0.1))
    plot(transfo, col=colorRampPalette(colors) (colsteps) [ findInterval(binned[[bb]][,4], seq(zr1[1],zr1[2], length.out=colsteps)) ]
         , pch = 20, yaxt="n", xaxt="n", 
         main=paste('sub-period ',bb,sep=' '),xlab='', ylab='',col.main= "#4EC1DE")
    arrows(0, 0, transfo[which.max(binned[[bb]][,4]),1], transfo[which.max(binned[[bb]][,4]),2],lwd=2)
  }
  abline(h=0,v=0,lwd=1,lty=2,col=3)
}

for(bb in 1:dt){
  screen(bb+5+dt)
  if(bb==1){
    par(pty="s") 
    par(mai=c(0.1,0.1,0.2,0.1))
    plot(transfo, col=colorRampPalette(colors) (colsteps) [ findInterval(binned[[bb]][,5], seq(zr1[1],zr1[2], length.out=colsteps)) ]
         , pch = 20,xlab='',ylab='',col.main= "blue", xaxt="n",cex.axis=0.9)
    mtext(expression(hat(C)[22]), side = 2, line = 3, adj = 0.5, font=3,col="#0086FF",cex=1.3)
    arrows(0, 0, transfo[which.max(binned[[bb]][,5]),1], transfo[which.max(binned[[bb]][,5]),2],lwd=2)
    #mtext('subperiod 2', side = 2, line = 3.5, adj = 0.5, font=3,col="#5970AF",cex=0.8)
  }else{
    par(pty="s") 
    par(mai=c(0.1,0.1,0.2,0.1))
    plot(transfo, col=colorRampPalette(colors) (colsteps) [ findInterval(binned[[bb]][,5], seq(zr1[1],zr1[2], length.out=colsteps)) ]
         , pch = 20, xlab='',ylab='',col.main= "blue", yaxt="n",xaxt="n")
    arrows(0, 0,transfo[which.max(binned[[bb]][,5]),1], transfo[which.max(binned[[bb]][,5]),2],lwd=2)
  }
  abline(h=0,v=0,lwd=1,lty=2,col=3)
  
}

for(bb in 1:dt){
  screen(bb+5+dt*2)
  if(bb==1){
    par(pty="s") 
    par(mai=c(0.1,0.1,0.2,0.1))
    plot(transfo, col=colorRampPalette(colors) (colsteps) [ findInterval(binned[[bb]][,6], seq(zr1[1],zr1[2], length.out=colsteps)) ]
         , pch = 20, xlab='',ylab='',col.main= "blue",cex.axis=0.8)
    mtext(expression(hat(C)[12]), side = 2, line = 3, adj = 0.5, font=3,col="#0086FF",cex=1.3)
    arrows(new_kappa_coords[bb,1],new_kappa_coords[bb,2], transfo[which.max(binned[[bb]][,6]),1], transfo[which.max(binned[[bb]][,6]),2],lwd=2)
    #mtext('subperiod 3', side = 2, line = 3.5, adj = 0.5, font=3,col="#5970AF",cex=0.8)
  }else{
    par(pty="s") 
    par(mai=c(0.1,0.1,0.2,0.1))
    plot(transfo, col=colorRampPalette(colors) (colsteps) [ findInterval(binned[[bb]][,6], seq(zr1[1],zr1[2], length.out=colsteps)) ]
         , pch = 20, xlab='',ylab='',col.main= "blue", yaxt="n",cex.axis=0.8)
    arrows(new_kappa_coords[bb,1],new_kappa_coords[bb,2], transfo[which.max(binned[[bb]][,6]),1], transfo[which.max(binned[[bb]][,6]),2],lwd=2)
  }
  abline(h=0,v=0,lwd=1,lty=2,col=3)
}

screen(2)
x1 <- c(0,0.15,0.15,0)
y1 <- c(0.25,0.25,0.7,0.7)
legend.gradient2(cbind(x1,y1),cols = colorRampPalette(colors) (colsteps), title = "", 
                 limits = seq(0,1,by=0.25), cex=1)
mtext("E-W distance (x 1000 km)", side=1, line = 7, adj = 3.05, font=2)
mtext("N-S distance (x 1000 km)", side = 2, line = 41.3, adj = 0.5, font=2)

close.screen( all=TRUE)
dev.off()

##-----------------------fig6b-----------------------##

library(mixtools)

max_lag <- 124

data_matrix <- data_format_into_matrix(data1 = obs1, data2 = obs2, temporal_replicates = ncol(obs1), simulated = F)

vel_est1 <- vel_est2 <- vel_est3 <- matrix(, ncol = 2, nrow = 247)

for(date_stamp in 1:247){
  cat(date_stamp,'\n')
  to_include <- seq(date_stamp, 6200,by = 248)
  to_include2 <- seq(date_stamp + 1, 6200, by = 248)
  
  conso_cor <- empirical_st_cov(data1 = t(data_matrix[insample_loc_index, to_include]), data2 = t(data_matrix[insample_loc_index, ncol(data_matrix)/2 + to_include2]), locations = coordinates[insample_loc_index,], max_time_lag = 1, simulated = F)
  binned <- empirical_covariance_dataframe(data_cov = conso_cor)
  binned_orig <- binned[which(binned[,3] == 1), ]
  vel_est1[date_stamp, ] <- binned_orig[which.max(binned_orig[,4]), 1:2]
  vel_est2[date_stamp, ] <- binned_orig[which.max(binned_orig[,5]), 1:2]
  vel_est3[date_stamp, ] <- binned_orig[which.max(binned_orig[,6]), 1:2]
}

vel_est <- rbind(vel_est1, vel_est2, vel_est3)

advec_est_fin <- vel_est/1000000

tiff('new_fig6b.tiff', units = "in", width = 4, height = 4, res = 300)

par(mai = c(0.9, 0.9, 0.3, 0.3))
par(pty = "s") 
plot(advec_est_fin, ylab = expression(v[y]), xlab = expression(v[x]), xlim = c(-8,8), ylim = c(-8,8), pch = 20, lwd = 0.4, col = '#BDBDBD')
abline(h = 0, v = 0, lwd = 1, lty = 2, col = 3)

ellipse(colMeans(advec_est_fin), cov(advec_est_fin), alpha = .5, col = 1, lty = 2)
ellipse(colMeans(advec_est_fin), cov(advec_est_fin), alpha = .05, col = 1, lty = 2)

mu <- c(-1.0200, -0.0789)
sigs <- matrix(c(7272.6663, -1496.1894, -1496.1894, 1049.0212), ncol = 2, byrow = T)/1000

mu2 <- c(-2.393, -0.7131)
sigs2 <- matrix(c(0.0034, 0, 0, 0.0001), ncol = 2, byrow = T)

mu3 <- c(-0.3303, 0.1645)
sigs3 <- matrix(c(0.0013, 0, 0, 0.0001), ncol = 2, byrow = T)

ellipse(mu, sigs, alpha = .5, col = '#6967CE', lty = 1, lwd = 1)
ellipse(mu, sigs, alpha = .05, col = '#6967CE', lty = 1, lwd = 1)

ellipse(mu2, sigs2, alpha = .5, col = "#A7F432", lty = 1, lwd = 1)
ellipse(mu2, sigs2, alpha = .05, col = "#A7F432", lty = 1, lwd = 1)

ellipse(mu3, sigs3, alpha = .5, col = "#A7F432", lty = 1, lwd = 1)
ellipse(mu3, sigs3, alpha = .05, col = "#A7F432", lty = 1, lwd = 1)

points(cbind(mu[1], mu[2]), pch = 8, col = "#6967CE", lwd = 1)
points(cbind(mu2[1], mu2[2]), pch = 8,col = "#A7F432", lwd = 1)
points(cbind(mu3[1], mu3[2]), pch = 8, col = "#A7F432", lwd = 1)
points(cbind(-497.4, -39.63)/1000, pch = 8, col = "#FA5B3D", lwd = 1)
points(cbind(-839.9, 747.9)/1000, pch = 8, col = "#5DADEC", lwd = 1)

legend("topleft", legend=c(expression(hat(bold("\u03bc"))^(M1)),expression(hat(bold(v))^(M2)), expression(list(hat(tilde(bold("\u03bc")))[1]^(M3),hat(tilde(bold("\u03bc")))[2]^(M3))),expression(hat(bold(v))^(M4))),
       col=c("#6967CE","#FA5B3D","#A7F432","#5DADEC"), pch = rep(8,4), cex=0.6,box.lty=0, inset=.02)

dev.off()

##-----------------------fig7-----------------------##

nu <- nu_mod5 <- c(3.218,3.736)
nu3 <- (nu[1]+nu[2])/2
beta_mod5 <-beta <- 794.8
nug <- nug_mod3 <- nug_mod4 <- nug_mod5 <- c(0,0)
var <- var_mod3 <- var_mod4 <- var_mod5 <-c(1,1)
rho <- 0.59

nu_mod3=c(8.991,2.3374)
beta_mod3=c(417.0461,1097.5639)
alpha <-matrix(c(0.8380,0.5456,0.0000,0.9999),ncol=2,byrow=T)

h <- seq(0,max(binned0[,1]))

hh=binned0[,5:6]%*%R
new_h <- sqrt(hh[,1]^2+hh[,2]^2)/1000

h <- seq(0,max(new_h))

theo_list <- list()

for(pp in 1:2){
  theo_list_temp <- ifelse(h!=0,var_mod3[pp]*(h/beta_mod3[pp])^nu_mod3[pp]*besselK(h/beta_mod3[pp],nu_mod3[pp])/(2^(nu_mod3[pp]-1)*gamma(nu_mod3[pp])),var_mod3[pp]+nug_mod3[pp])
  theo_list[[pp]] <- theo_list_temp
}

nu_mod4 <- c(3.757,4.394,(3.757+4.394)/2)
beta_mod4 <- 723.9
mu=c(0.499,0.674)
nu.L=5.846
lambda <-0.001
scale <-   1948
rho_mod4 <- 0.55

#2.9901    5.9380  759.1722    0.0010 2905.6495    0.6396    0.8280    8.0012 0.4510    0.9994 -126.1124  -58.4599

pdf('fig7.pdf', width=8, height=3)

par(mfrow = c(1,3))
par(pty="s") 
plot(new_h,binned0[,2],pch=3,ylab='Correlation',xlab='Spatial Lag (km)',
     main='', col=1, ylim=c(0,1), col.main= "#4EC1DE")
mtext(expression(C[11]), side = 3, line = 0, adj = 0.5, cex = 0.65, font=2)
i=1

theo <- ifelse(h!=0,var[i]*(h/beta[i])^nu[i]*besselK(h/beta[i],nu[i])/(2^(nu[i]-1)*gamma(nu[i])),var[i]+nug[i])
lines(h,theo,col='#FA5B3D',lwd=2)

theo.temp <- ifelse(h!=0,var_mod4[i]*(h/beta_mod4)^nu_mod4[i]*besselK(h/beta_mod4,nu_mod4[i])/(2^(nu_mod4[i]-1)*gamma(nu_mod4[i])),var_mod4[i]+nug_mod4[i])
lagrangian=pmax((1-h/scale),0)^(nu.L+mu[i])
theo <- (1-lambda)*theo.temp+lambda*lagrangian
lines(h,theo,col='#5DADEC',lwd=2)

theo <- ifelse(h!=0,var_mod5[i]*(h/beta_mod5[i])^nu_mod5[i]*besselK(h/beta_mod5[i],nu_mod5[i])/(2^(nu_mod5[i]-1)*gamma(nu_mod5[i])),var_mod5[i]+nug_mod5[i])
lines(h,theo,col='#FFAA1D',lwd=2)

theo <- alpha[1,1]^2*theo_list[[1]]+alpha[1,2]^2*theo_list[[2]]
lines(h,theo,col='#A7F432',lwd=2)

legend("bottomleft", legend=c("M1 & M2", "M3","M4", "M5"),
       col=c("#FA5B3D","#A7F432","#5DADEC","#FFAA1D"), lty = rep(1,3), cex=0.8,box.lty=0, inset=.02)

####--------------------

plot(new_h,binned0[,3],pch=3,ylab='',xlab='Spatial Lag (km)',
     main='', col=1, ylim=c(0,1), col.main= "#4EC1DE")
mtext(expression(C[22]), side = 3, line = 0, adj = 0.5, cex = 0.65, font=2)

i=2
theo <- ifelse(h!=0,var[i]*(h/beta)^nu[i]*besselK(h/beta,nu[i])/(2^(nu[i]-1)*gamma(nu[i])),var[i]+nug[i])
lines(h,theo,col='#FA5B3D',lwd=2)

theo.temp <- ifelse(h!=0,var_mod4[i]*(h/beta_mod4)^nu_mod4[i]*besselK(h/beta_mod4,nu_mod4[i])/(2^(nu_mod4[i]-1)*gamma(nu_mod4[i])),var_mod4[i]+nug_mod4[i])
lagrangian=pmax((1-h/scale),0)^(nu.L+mu[i])
theo <- (1-lambda)*theo.temp+lambda*lagrangian
lines(h,theo,col='#5DADEC',lwd=2)

theo <- ifelse(h!=0,var_mod5[i]*(h/beta_mod5)^nu_mod5[i]*besselK(h/beta_mod5,nu_mod5[i])/(2^(nu_mod5[i]-1)*gamma(nu_mod5[i])),var_mod5[i]+nug_mod5[i])
lines(h,theo,col='#FFAA1D',lwd=2)

theo <- alpha[2,1]^2*theo_list[[1]]+alpha[2,2]^2*theo_list[[2]]
lines(h,theo,col='#A7F432',lwd=2)

legend("bottomleft", legend=c("M1 & M2", "M3","M4", "M5"),
       col=c("#FA5B3D","#A7F432","#5DADEC","#FFAA1D"), lty = rep(1,3), cex=0.8,box.lty=0, inset=.02)

####--------------------
plot(new_h,binned0[,4],pch=3,ylab='', col=1,xlab='Spatial Lag (km)',
     main='', col.main= "#4EC1DE",ylim=c(0,1))
#mtext(expression(list(nu[12]==2.46,a[12]==439)), side = 3, line = 0, adj = 0.5, cex = 0.65, font=2)
mtext(expression(C[12]), side = 3, line = 0, adj = 0.5, cex = 0.65, font=2)

#h <- seq(0,2800)
new_hh <- abs(seq(-734,length(h)-734-1,by=1))
theo <- ifelse(new_hh!=0,rho*sqrt(var[1]*var[2])*(new_hh/beta)^nu3*besselK(new_hh/beta,nu3)/(2^(nu3-1)*gamma(nu3)),rho*sqrt(var[1]*var[2]))
lines(h,theo,col='#FF3855',lwd=2)

theo.temp <- ifelse(h!=0,rho_mod4*sqrt(var_mod4[1]*var_mod4[2])*(h/beta_mod4)^nu_mod4[i]*besselK(h/beta_mod4,nu_mod4[i])/(2^(nu_mod4[i]-1)*gamma(nu_mod4[i])),rho_mod4*sqrt(var_mod4[1]*var_mod4[2]))
lagrangian=pmax((1-h/scale),0)^(nu.L+mu[i])
theo <- (1-lambda)*theo.temp+lambda*lagrangian
lines(h,theo,col='#5DADEC',lwd=2)

theo <- ifelse(h!=0,0.52*sqrt(var_mod5[1]*var_mod5[2])*(h/beta_mod5)^nu3*besselK(h/beta_mod5,nu3)/(2^(nu3-1)*gamma(nu3)),0.52*sqrt(var_mod5[1]*var_mod5[2]))
lines(h,theo,col='#FFAA1D',lwd=2)

theo <- alpha[1,1]*alpha[2,1]*theo_list[[1]]+alpha[1,2]*alpha[2,2]*theo_list[[2]]
lines(h,theo,col='#A7F432',lwd=2)

legend("topright", legend=c("M1 & M2", "M3","M4", "M5"),
       col=c("#FA5B3D","#A7F432","#5DADEC","#FFAA1D"), lty = rep(1,3), cex=0.8,box.lty=0, inset=.02)
dev.off()

############################################################################
############################################################################
########################         APPENDIX           ########################
############################################################################
############################################################################

#--*Figure 2
data_plots_functions('testdata.pdf', dat1, dat2)

filename = 'temporalstationarity.pdf'

pdf(file=filename, width=8, height=8)

par(mfrow=c(10,2))

lims <- ceil(max(abs(min(c(dat1[1:5,], dat2[1:5,]))), abs(max(c(dat1[1:5,], dat2[1:5,])))))
for(ll in 1:9){
  if(ll < 9){
    par(mai=c(0,0.5,0.2,0.1))
    plot(dat1[ll,] - mean(dat1[ll,]), type='l', ylim = c(-lims, lims), xaxt = 'n', xlab = '', ylab = '')
    if(ll == 1){
      mtext('1000 hPa', side = 3, line = 0, adj = 0.5, cex = 1, font = 2, col="#4EC1DE")
    }
    par(mai=c(0,0.1,0.2,0.4))
    plot(dat2[ll,] - mean(dat2[ll,]), type='l', ylim = c(-lims, lims), xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
    if(ll == 1){
      mtext('975 hPa', side = 3, line = 0, adj = 0.5, cex = 1, font = 2, col="#4EC1DE")
    }
    mtext(paste('Station', ll), side = 4, line = 0, adj = 0.5, cex = 0.5, font = 2, col="#0086FF")
  }else{
    par(mai=c(0,0.5,0.2,0.1))
    plot(dat1[ll,] - mean(dat1[ll,]), type='l', ylim = c(-lims, lims), xlab = '', ylab = '')
    par(mai=c(0,0.1,0.2,0.4))
    plot(dat2[ll,] - mean(dat2[ll,]), type='l', ylim = c(-lims, lims), yaxt = 'n', xlab = '', ylab = '')
    mtext(paste('Station', ll), side = 4, line = 0, adj = 0.5, cex = 0.5, font = 2, col="#0086FF")
  }
}

mtext('log PM Concentrations', side = 2, line = 29.5, adj = -2.1, cex = 1, font = 3, col="#0086FF")
mtext('Temporal Locations', side = 1, line = 3, adj = -0.5, cex = 1, font = 2, col="#4EC1DE")
dev.off()