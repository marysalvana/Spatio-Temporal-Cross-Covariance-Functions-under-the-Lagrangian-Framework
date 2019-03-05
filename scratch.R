MLEs <- optim( par = MLEs.save[1,] , fn = f_loglik, control=list(maxit=50000,parscale=c(lam1.init, lam2.init, pi/4,lam1.init, lam2.init, pi/4, thets,1,1,1,1),trace=5))

new_thet <- MLEs$par
MLEs <- optim( par = new_thet, fn = f_loglik, control=list(maxit=50000,parscale=c(lam1.init, lam2.init, pi/4,lam1.init, lam2.init, pi/4, thets,1,1,1,1),trace=5))

materncov <- matern_cov_regular_grid_v4(thets,c(.05001,.05001),time=t)
A <- mvrnorm(1, rep(0,dim(materncov)[1]),materncov)

k=l=1
index <- (k-1)*NKR+ local.index + Nk*nkr*(l-1)

cov.chol <- chol(materncov[c(index,index+400),c(index,index+400)])
tmp1 <- backsolve(cov.chol,  A[c(index,index+400)], transpose = TRUE)
#tmp1 <- backsolve(cov.chol, c(data)-c(rep(b0[1],dim(data)[1]/2),rep(b0[2],dim(data)[1]/2)), transpose = TRUE)

ResCinvRes <- t(tmp1) %*% tmp1

new_cov_target <- materncov[c(index,index+400),c(index,index+400)]
image.plot(new_cov_target)

make.local.loglik <- local_loglik_v3(locations = temp.locations, 
                                  data = temp.data, p = p1[l+(k-1)*Nk,], time=t) 
test_theta <- rep(0.1,12)
MLEs.local <- optim(test_theta, make.local.loglik,control=list(maxit=10000, parscale=test_theta, trace=5))

test_theta <- c(MLEs.local$par,p[1:6])
MLEs.local <- optim(test_theta, make.local.loglik,control=list(maxit=10000, parscale=test_theta, trace=5))

for(lalai in 1:100){
  test_theta <- MLEs.local$par+0.01
  MLEs.local <- optim(test_theta, make.local.loglik,control=list(maxit=10000, parscale=test_theta, trace=5))
}

temp.data <- A[c(index,index + 400)]

make.local.loglik <- local_loglik_v4(locations = temp.locations, 
                                     data = temp.data, p = p1[l+(k-1)*Nk,], time=t) 
test_theta <- MLEs.local$par
MLEs.local <- optim(test_theta, make.local.loglik,control=list(maxit=10000, trace=5))

test_theta <- rep(0,18)
MLEs.local <- optim(test_theta, make.local.loglik,control=list(maxit=10000, trace=5))

test_theta <- MLEs.local$par
MLEs.local <- optim(test_theta, make.local.loglik,control=list(maxit=10000, trace=5))


test_theta <- test_theta_old
MLEs.local <- optim(test_theta, make.local.loglik,control=list(maxit=10000, trace=5))

# -1.59772094   1.04801227   8.02758410 -10.85316793  -1.30693969  45.16599423
# 4.83542356  -5.77314258   3.33519760   5.58629767   3.00216171  -7.75165332
# -3.51616760  -3.31513869  -8.09289960  -0.56906147   0.07647102   1.36922142

#  -1.5051897   1.0997939   8.0664497 -10.7594767  -1.6648476  56.1657259
#   4.9116717  -5.6707770   3.3128773   5.5537447   2.9105141  -7.7747117
#  -3.5055135  -3.3196077  -9.5537869  -0.5197488   0.2689982   1.4082707

make.local.loglik <- local_loglik(locations = temp.locations, 
                                  data = temp.data, p = p1[l+(k-1)*Nk,], time=t) 
test_theta <- rep(0,12)
MLEs.local <- optim(test_theta, make.local.loglik,control=list(maxit=87, parscale=test_theta, trace=5))


cov.chol1 <- chol(S)
cov.chol2 <- chol(new_cov_target)

cov.chol <- chol(S)
tmp1 <- backsolve(cov.chol,  A[c(index,index+400)], transpose = TRUE)

t(tmp1) %*% tmp1
