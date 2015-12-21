#' usvar
#'
#' Solve a underidentified Structural Autogressive models SVAR models
#'
#' @param p number of lags in the model
#' @param k number of dependent Variables
#' @param f zero restriction matrix
#' @param sr sign restriction matrix
#' @param B coefficient matrix
#' @param C1 initial impact response
#' @param t length of impulse response function
#' @param draws Number of draws to find
#' @param var_pos ...
#'
#' @return the function returns a list that stores the impulse response for every draw,
#'         and each draw is a nested list that contains the impulse response matrix
#'         for each response variable.
#' @export
usvar <- function(p, k, f, sr, B, C1, t, draws, var_pos){
  # algorithm to solve underidentified SVAR models
  # p = number of lags
  # k = number of dependent variables
  # f = zero restriction matrix, 2k-by-k
  # sr = sign restriction matrix to impose on impact
  # B = coefficient matrix
  # C1 = initial impact response
  # t = length of impulse response function
  # draws = number of draws

  shocks <- diag(1,k)

  Qs <- findQs(k,f)
  Q <- Qs$Q
  index <- Qs$index
  flag <- Qs$flag

  if (flag == 1)
    stop("Rank condition not satisfied, the model is overidentified")

  Btilde <- t(B[-1,])
  # build the companion matrix
  temp <- cbind(diag(1,k*(p-1)),matrix(0, k*(p-1),k))
  alpha <- rbind(Btilde, temp)

  counter <- 1

  # stores the impulse response function at all draws
  R <- list()

  while (counter < draws + 1){

    C <- generateDraw(C1,k) # generate the new short-run impact matrix

    P <- findP(C,B,Q,p,k,index) # find orthogonal rotation matrix

    W <- C%*%P # rotated short-run impact matrix

    # stores the shocks at each draw
    shock.list <- list()

    # ordering columns
    for (j in 1:k){

      if (W[var_pos[j],j] < 0)
        shock <- -shocks[,j]
      else
        shock <- shocks[,j]


      V <- matrix(0, k*p, t)

      V[1:k,1] <- W%*%shock

      chk <- W%*%shock

      sr_index <- which(!is.na(sr[,j]))
      # check the sign restrictions
      tmp <- sign(chk[sr_index]) - sr[sr_index,j]

      if (any(tmp!=0)){
        j <- 0
        break
      }

      for (i in 2:t){
        V[,i] <- alpha%*%V[,i-1]
      }

      # store the impulse response functions under each draw
      draw.name <- paste('draw',counter, sep='')
      shock.name <- paste('shock',j, sep='')

      shock.list[[shock.name]] <- V[1:k,]
      shock.list[['P.matrix']] <- P
      shock.list[['W.matrix']] <- W
      shock.list[['C.matrix']] <- C

      R[[draw.name]] <- shock.list

    }
    if (j == k)
      counter <- counter + 1
      print(counter)
  }
  return(R)
}

#' xlags
#'
#' Create lagged values for a data matrix where every column contains data for a variable
#'
#' @param data A matrix containing data for lagging, where each column contains time series data for a variable
#' @param p Number of lags
#'
#' @return the function returns a data matrix after lag operation
#' @export
xlags <- function(data, p){
  # create lagged values
  # input:
  # data = A matrix containing original predictor time series
  #        each column contains data for one predictor
  # p = Number of laggs
  # output:
  # xlags = lagged predictor time series data
  n <- nrow(data)
  m <- ncol(data)
  X <- matrix(0, n-p, m*p)

  for (j in 1:p){
    X[,(j-1)*m + 1:m] <- as.matrix(data[(p-j+1):(n-j),])
  }

  return(X)
}


#' findQs
#'
#' find the Q matrices
#'
#' @param f zero restriction matrix
#' @param k number of predictors
#'
#' @return the function returns a list which contains three elements.
#'         Q is a list containting the linear restrictions for each column
#'         index is ordered index
#'         flag is the status of the system, where "1" means Overidentified,
#'         "0" means exactly indentified and "-1" means underidentified
#
#' @export
findQs <- function(k,f){
  # find the Q matrices
  # input:
  # f = zero restriction matrix
  # k = number of predictors
  #
  # output:
  # res is a list which contains
  # Q = A list containting the linear restrictions for each column
  # index = ordered index
  # flag =
  #
  E = diag(rep(1,k));

  Q_init_mat <- vector("list",k)
  Q_init_rank <- rep(NA,k)

  for(i in 1:k){

    Q_init_mat[[i]] <- diag(as.numeric(ifelse(f %*% E[,i] == 0,1,0)))
    Q_init_rank[i] <- sum(svd(Q_init_mat[[i]])$d)

  }

  for(i in 1:k){

    temp <- Q_init_mat[[i]]
    temp2 <- matrix(rep(0,k*k*2), nrow = k)
    keep.index <- which(apply(temp,2,sum) > 0)
    if(length(keep.index) > 0){
      temp <- temp[keep.index,]
      temp2[1:length(keep.index),] <- temp
    }
    Q_init_mat[[i]] <- temp2

  }

  ordering <- order(Q_init_rank, decreasing = TRUE)
  Q_init_rank <- Q_init_rank[ordering]

  index <- rep(NA,k)

  for(i in 1:k){

    index[ordering[i]] <- i

  }

  Q <- vector("list",k)

  for(i in 1:k){

    Q[[i]] <- Q_init_mat[[ordering[i]]]

  }

  if(any(Q_init_rank - (k - (1:k)) > 0)){  #  over-identified
    flag = 1
  }else if(all(Q_init_rank  - (k - (1:k)) == 0)){ # exactly identified
    flag = 0
  }else if(any(Q_init_rank  - (k - (1:k)) < 0)) { #  under-identified
    flag = -1
  }

  res <- list(Q = Q, index = index, flag = flag)

  return(res)

}

#' findP
#'
#' find the orthogonal rotation matrix P
#'
#' @param C initial short run impact matrix, usually from a cholesky
#' @param B matrix of coefficients (including intercept estimates)
#' @param Q a list containting the linear restrictions for each column
#' @param p number of lags
#' @param k number of dependent variabels
#' @param ordering original column ordering in the matrix of restrictions
#' @return the function returns an orthogonal rotation matrix
#
#' @export
findP <- function(C,B,Q,p,k,ordering){
  # inputs:
  # C = initial short run impact matrix, usually from a cholesky
  # decomposition of the forecast error variance
  # B = Matrix of coefficients (including intercept estimates)
  # Q = A cell containing the linear restrictions for each columnn
  # p = number of lags
  # k = number of dependent variables
  # ordering = original column ordering in the matrix of restrictions
  #
  # outputs:
  # P = orthogonal rotation matrix
  L0 <- C

  beta_temp <- t(B[2:nrow(B),])

  beta <- matrix(rep(0,k*k), ncol = k)

  for(i in 1:p){

    beta <- beta + beta_temp[,(1:k) + (i-1)*k]

  }
  # solve , reverse
  Linf <- Matrix::solve(diag(rep(1,k)) - beta, C)
  ## cbind -> rbind
  Fmat <- rbind(L0,Linf)

  P <- matrix(rep(0,k*k),ncol = k)

  for(i in 1:k){

    if(i == 1){
      ## F -> Fmat
      Qtilde <- Q[[i]] %*% Fmat;
    }else {
      ## cbind -> rbind
      Qtilde = rbind(Q[[i]] %*% Fmat, t(P))
    }

    qr.mat <- qr(t(Qtilde))
    QQ <- qr.Q(qr.mat);
    P_temp <- QQ[,ncol(QQ)];
    P[,i] <- P_temp

  }

  # () -> []
  P <- P[,ordering]

  return(P)

}

#' generateDraw
#'
#' generate a orthogonal matrix by taking the QR
#' decomposition of a random matrix, and post multiply
#' the initial matrix by this random orthogonal draw
#'
#' @param C initial short run impact matrix
#' @param k number of dependent variables
#'
#' @return the function returns a randomised initial short-run impact matrix
#
#' @export
generateDraw <- function(C,k){
  # Input:
  # C = initial short run impact matrix
  # k = number of columns/number of predictors
  #
  # Output:
  # IRF after rotation
  newmatrix <- matrix(rnorm(k*k), ncol = k)

  qr.mat <- qr(newmatrix)
  Q <- qr.Q(qr.mat)
  R <- qr.R(qr.mat)

  for(i in 1:k){

    if(R[i,i] < 0) Q[,i] <- -1*Q[,i]

  }

  Z <- C %*% Q

  return(Z)

}

#' data.transform
#'
#' transform the list containing all the impulse response functions
#' into a dataframe with flat structure
#'
#' @param nos number of shocks
#' @param nov number of dependent variables
#' @param data.list a list containing all the IRF draws
#' @param fwd forward steps in impulse response function
#' @param variable.name a character vector containing name of dependent variables
#' @param shock.name a character vector containing shock names
#'
#' @return the function returns a dataframe containing the impulse response function with user-defined forward steps
#
#' @export
data.transform <- function(nos, nov, data.list, fwd, variable.name, shock.name){
  # transform the list containing all the IRFs into a dataframe
  # input:
  # nos = number of shocks
  # nov = number of variables
  # data.list = a list containing the data of IRFs
  # fwd = forward lags
  # variable.name = name of response variables
  # shock.name = name of shocks
  # output:
  # irf
  irf <- list()
  for (j in 1:nos){
    for (i in 1:nov){
      sn <- paste("shock", j, sep='')
      temp.shock <- lapply(data.list, function(x) x[[sn]][i,1:fwd])
      temp.shock <- reshape2::melt(data.frame(do.call(cbind, temp.shock)))
      step <- data.frame(step = rep(1:fwd,draws), shock_name = shock.name[j],
                         var_name = variable.name[i])
      irf[[(j-1)*nos+i]] <- cbind(step,temp.shock)
    }
  }
  irf <- do.call(rbind, irf)
  return(irf)
}

#' find.base.draw
#'
#' find the baseline draw whose corresponding impulse responses
#' has the minimimum sum of distance squared to the median value
#'
#' @param irf.all A list containing all the IRF draws
#'
#' @return the function returns the baseline draw index
#
#' @export
find.base.draw <- function(irf.all){
  # find the baseline draw by minimizing the distance to the median
  # input:
  # irf.all = data containing all the IRF draws
  # output:
  # draw.name = the baseline draw number

  # for each draw, each imf, compute the distance to median
  irf.all.median <- plyr::ddply(irf.all, c("var_name","shock_name","step"), summarise, value = median(value))
  irf.tmp <- merge(irf.all, irf.all.median, by=c("step", "var_name","shock_name"), all.x=TRUE)
  irf.tmp$meddis <- (irf.tmp$value.x-irf.tmp$value.y)^2

  # for each draw, compute the sum of squares
  dist.to.med <- plyr::ddply(irf.tmp,c("variable"), summarise, sos = sum(meddis))
  # find the min
  draw.ind <- which(dist.to.med[,2] == min(dist.to.med[,2]))
  # find the corresponding draw
  draw.name <- paste("draw",draw.ind,sep='')
  return(draw.name)
}


#' hist.decomp
#'
#' find the historical decomposition for model residuals.
#'
#' @param Z the matrix that transforms structral shocks to shocks
#' @param u model residuals
#' @param B estimated coefficient matrix
#' @param X.mat data for the predictors
#' @param nov number of variable
#' @param nos number of shocks
#' @param var.name variable names
#' @param shock.name shock names
#'
#' @return the function returns an array containing the historical decompositions 
#'         for the model residuals.
#'         
#'         The dimensions of the array are response variables, shocks, time respectively
#
#' @export
hist.decomp <- function(Z, u, B, X.mat, nov, nos, t.last, var.name, shock.name){
  # Z - transform structral shocks to shocks
  # u - residuals
  # B - coefficient matrix
  # X.mat - predictor matrix
  # nov - number of response variables
  # nos - number of shocks
  # t.last - steps ahead for decomposition
  # Date - date index for the data
  
  e <- Matrix::solve(Z, t(u))
  
  # initialize array to store shocks, each matrix from time t=1 to t=t.last
  # has response variables as its rows and shocks as its columns
  contrib <- array(data = NA, dim = c(nov, nos, t.last))
  const  <- matrix(NA, nrow = nov, ncol = t.last)
  
  A <- t(B)[,-1]
  A.power <- lapply(0:t.last, function(n,A) A %^% n, A = A)
  
  # recover historical decompositions
  for (t in 1:t.last){
    
    tmp <- 0
    
    for (k in 1:t){
      
      tmp <- (A.power[[t-k+1]] %*% Z) * matrix(rep(e[,k],2), nrow = nov, byrow = TRUE) + tmp
      
    }
    
    contrib[,,t] <- tmp
    
    tmp <- 0
    
    for(k in 1:t){
      
      tmp <- A.power[[t-k+1]] %*% B[1,] + tmp
      
    }
    
    const[,t] <- A.power[[t+1]] %*% X.mat[1,2:ncol(X.mat)] + tmp
    
  }
  
  rownames(contrib) <- var.name
  colnames(contrib) <- shock.name
  
  return(contrib)
}



