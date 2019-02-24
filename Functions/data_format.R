data_format_into_matrix <- function(data1, data2 = NULL, temporal_replicates, aggs_max = 1, simulated = T){
  
  # data1, data2: data matrix with temporal replicates as columns and each row represents the locations
  # aggs_max: the number of consecutive timesteps you want to average. Here we take the data as is and do not take averages since they are already Gaussian and stationary
  
  # output: data matrix with the first t columns for variable 1 and the t+1 column to 2t column for variable 2
  
  if(!simulated == T){
    ave_var1 <- ave_var2 <- matrix(,ncol=floor(temporal_replicates/aggs_max),nrow=nrow(data1))
    for(locat in 1:nrow(data1)){
      new_uu <- new_uu2 <- matrix(,ncol=floor(temporal_replicates/aggs_max),nrow=aggs_max)
      
      for(gg in 1:floor(temporal_replicates/aggs_max)){
        new_uu[,gg] <- data1[locat,((gg-1)*aggs_max+1):(gg*aggs_max)]
        new_uu2[,gg] <- data2[locat,((gg-1)*aggs_max+1):(gg*aggs_max)]
      }
      ave_var1[locat,]<- colMeans(new_uu)
      ave_var2[locat,]<- colMeans(new_uu2)
    }
    
    return(cbind(ave_var1, ave_var2))
  }
}