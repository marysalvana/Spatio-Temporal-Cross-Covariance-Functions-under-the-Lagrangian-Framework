toeplitz_mat <- function(S_list){
  k <- min(unlist(lapply(S_list, dim)))
  n <- length(S_list)
  #
  # Create the "strip".
  #
  strip <- array(NA, dim=c(k,k,2*n-1))
  for (i in 1:n) strip[,,i] <- S_list[[n+1-i]]
  if (n > 1) for (i in 2:n) strip[,,n+i-1] <- t(S_list[[i]])
  #
  # Assemble into "block-Toeplitz" form.
  #
  X <- array(NA, dim=c(k,k,n,n))
  # Blast the strip across X.
  #
  for (i in 1:n) X[,,,i] <- strip[,,(n+1-i):(2*n-i)]
  X <- matrix(aperm(X, c(1,3,2,4)), n*k)
}
