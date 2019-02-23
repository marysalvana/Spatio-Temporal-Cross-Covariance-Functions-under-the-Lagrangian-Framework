

## LOAD ALL PACKAGES FROM [1]preprocessing.R

WORKSTATION = 1

# 0.1 SET WORKING DIRECTORY WHERE DATA IS SAVED

if(WORKSTATION == 1){
  root <- '/Volumes/GoogleDrive/My Drive/Summer 2017_2018/Project 1/'
}else if (WORKSTATION == 2){
  root <- 'G://My Drive/Summer 2017_2018/Project 1/'
}else{
  root <- '/Users/salvanmo/silver-octo-spork/'
}

load(paste(root,'DATA/UTM_simulated_locations_only.RData',sep=''))
load(paste(root,'DATA/simulated_locations_only.RData',sep=''))

WORKSTATION = 3

# 0.1 SET WORKING DIRECTORY WHERE DATA IS SAVED

if(WORKSTATION == 1){
  root <- '/Volumes/GoogleDrive/My Drive/Summer 2017_2018/Project 1/'
}else if (WORKSTATION == 2){
  root <- 'G://My Drive/Summer 2017_2018/Project 1/'
}else{
  root <- '/Volumes/GoogleDrive/My Drive/Phd_Dissertation/ideal-happiness/Functions'
}

setwd(paste(root,'/Figures',sep=''))

source(file=paste(root,"/empirical_st_cov.R",sep=''))

#You can use raw netcdf data and pre-process it using preprocessing.R or you load this already preprocessed data



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
