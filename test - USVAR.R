context("USVAR Functions Testing")

# test for findQ
# using both zero restrictions from first example in Binning's Paper
# 
test_that("findQ", {
  
  f <- c(0,0,1,1,
         1,1,1,1,
         0,1,1,1,
         1,1,1,1,
         0,0,0,1,
         1,1,1,1,
         1,1,1,1,
         1,1,1,1)
  m <- 4
  f <- matrix(f,nrow=m*2, ncol=m, byrow=TRUE)
  Qs.trial <- findQs(m, f)
  
  # True Qs, rank, and identification
  Qs.1 <- c(1,0,0,0,0,0,0,0,
            0,0,1,0,0,0,0,0,
            0,0,0,0,1,0,0,0,
            0,0,0,0,0,0,0,0)
  
  Qs.2 <- c(1,0,0,0,0,0,0,0,
            0,0,0,0,1,0,0,0,
            0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0)
  
  Qs.3 <- c(0,0,0,0,1,0,0,0,
            0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0)
  
  Qs.4 <- rep(0,m*m*2)
  
  r.index <- c(1,2,3,4)
  id <- 0
  
  mnames <- c("Qs.1","Qs.2","Qs.3","Qs.4")
  eval.mnames <- lapply(mnames,get)
  Qs <- lapply(eval.mnames, matrix, nrow=m,ncol=m*2, byrow=TRUE)
  
  Qs.true <- list(Q=Qs, index=r.index, flag=id)
  
  expect_equal(Qs.trial, Qs.true)
  
}
)


# test for findP
# using both zero restrictions and Qs from first example in Binning's Paper
#
test_that("findP", {
  m <- 4
  # generate a random covariance matrix, sigma
  rand.matrix <- matrix(rnorm(m^2), m, m)
  rm.eig <- eigen(rand.matrix+t(rand.matrix))
  # eigen decomposition
  rand.sigma <- rm.eig$vectors %*% diag(abs(rm.eig$values)) %*% t(rm.eig$vectors)
  # lags
  p <- 1
  # random coefficient matrix, B = transpose(B) in the context, 
  # where the intercepts estimates are the first row
  B <- matrix(rnorm(m*(m*p+1)), m*p+1, m)
  
  # take the first trial zero restrictions
  f <- c(0,0,1,1,
         1,1,1,1,
         0,1,1,1,
         1,1,1,1,
         0,0,0,1,
         1,1,1,1,
         1,1,1,1,
         1,1,1,1)
  f <- matrix(f,nrow=m*2, ncol=m, byrow=TRUE)
  
  Qs.1 <- c(1,0,0,0,0,0,0,0,
            0,0,1,0,0,0,0,0,
            0,0,0,0,1,0,0,0,
            0,0,0,0,0,0,0,0)
  
  Qs.2 <- c(1,0,0,0,0,0,0,0,
            0,0,0,0,1,0,0,0,
            0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0)
  
  Qs.3 <- c(0,0,0,0,1,0,0,0,
            0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0)
  
  Qs.4 <- rep(0,m*m*2)
  
  Qs.1 <- matrix(Qs.1, 4, 8, byrow=TRUE)
  Qs.2 <- matrix(Qs.2, 4, 8, byrow=TRUE)
  Qs.3 <- matrix(Qs.3, 4, 8, byrow=TRUE)
  Qs.4 <- matrix(Qs.4, 4, 8, byrow=TRUE)
  
  Q <- list(Qs.1,Qs.2,Qs.3,Qs.4)
  # # of lags
  # find P
  
  P <- findP(t(chol(rand.sigma)), B, Q, p, m, c(1,2,3,4))
  
  # form Fmat
  k <- m
  
  L0 <- t(chol(rand.sigma))
  
  beta_temp <- t(B[2:nrow(B),])
  
  beta <- matrix(rep(0,k*k), ncol = k)
  
  for(i in 1:p){
    
    beta <- beta + beta_temp[,(1:k) + (i-1)*k]
    
  }
  
  Linf <- solve(diag(rep(1,k)) - beta, L0)
  
  Fmat <- rbind(L0,Linf)
  
  FP <- Fmat%*%P
  
  expected <- which(f == 0)
  
  expect_equal(which(abs(FP) < 1e-10), expected)
  
})

# test for xlags
# use random data with 10 elements
test_that("xlags", {
  
  data <- data.frame(runif(10))
  
  trial <- xlags(data,1)
  
  true <- as.matrix(data[1:9,])
  
  expect_equal(trial,true)
})

