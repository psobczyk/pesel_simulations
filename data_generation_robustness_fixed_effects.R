#' 
#' Piotr Sobczyk, PWr
#' Simulation schemes
#' 

compare_methods <- function(simulation_function, numb.repetitions = 10, 
                            n = 50, SNRs = c(1), vars = c(50,100), k = 2, 
                            scale=TRUE, maxPC = NULL, id='A', 
                            pathToMatlab="/usr/local/MATLAB/R2013a/bin/matlab"){
  require(pesel)
  require(FactoMineR)
  if(is.null(maxPC)) maxPC = 2*k
  
  methods.names <- c("Minka BIC", "Laplace evidence", "Properly dervied BIC",
                     "GCV", "Rajan - penalized version", "Passemier", "CSV")
  
  results <- NULL
  
  tabela <- matrix(0, nrow=length(methods.names), ncol=numb.repetitions)

  for(i in 1:length(SNRs)){
    for(j in 1:length(vars)){
      tabela <- matrix(0, nrow=length(methods.names), ncol=numb.repetitions)
      SNR <- SNRs[i]
      var <- vars[j]
      
      SIGNAL = simulation_function(n = n, SNR = SNR, numb.vars = var, 
                                   k = k, scale = scale)$signal
      
      for(rep in 1:numb.repetitions){
        sim.data <- simulation_function(SIGNAL=SIGNAL, n = n, SNR = SNR, numb.vars = var, 
                                        k = k, scale = scale)
        
        BICs <- pesel(sim.data$X, 1, maxPC, method = "heterogenous", asymptotics = "n")$vals
        tabela[1,rep] <- which.max(BICs)
        BICs <- sapply(1:maxPC, function(j) pca.Laplace(sim.data$X, j))
        tabela[2,rep] <- which.max(BICs)
        BICs <- pesel(sim.data$X, 1, maxPC, method = "heterogenous", asymptotics = "p")$vals
        tabela[3,rep] <- which.max(BICs)
        estimatedNpc <- estim_ncp(sim.data$X, ncp.min = 1, ncp.max = maxPC, method =  "GCV")$ncp
        tabela[4,rep] <- estimatedNpc
        BICs <- pesel(sim.data$X, 1, maxPC, method = "homogenous", asymptotics = "p")$vals
        tabela[5,rep] <- which.max(BICs)
        
        xFileName = paste0("temp_data/x_", id, ".csv")
        write.table(sim.data$X, sep=",", dec=".", file = xFileName, col.names=FALSE, row.names = FALSE)
        command = paste0(pathToMatlab, " -nodisplay -r \"id='", id, 
                         "'; kmax=", maxPC, "; run('estimating_pcs_passemier2.m'); exit\"")
        system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)
        resFileName = paste0("temp_data/res_", id, ".csv")
        result = read.csv(resFileName)
      
        result = read.csv(resFileName)
        if(is.na(result$sigma)){
          tabela[6,rep] <- 1
        } else{
          tabela[6,rep] <- min(result$k,maxPC)
        }
        
        

        n.folds <- 10
        folds <- NULL
        p <- var
        fold.size <- n*p/n.folds
        free.indices <- 1:(n*p)
        selected <- NULL
        temp <- sample(free.indices, size = fold.size, replace = FALSE)
        folds[temp] <- 1
        selected <- c(selected, temp)
        for(fold.i in 2:n.folds){
          temp <- sample(free.indices[-selected], size = fold.size, replace = FALSE)
          folds[temp] <- fold.i
          selected <- c(selected, temp)
        }
      
        ### creating lambda grid
        lam0 <- lambda0(sim.data$X)
        lamseq=exp(seq(from=log(lam0),to=log(1),length=20))
      
        ### cross-validation
        error <- matrix(0, nrow=n.folds, ncol=20)
        for(fold.i in 1:n.folds){
          xs <- sim.data$X
          xs[which(folds==fold.i)] <- NA
          for(lambda.j in seq_along(lamseq)){
            lambda=lamseq[lambda.j]
            fits=softImpute(xs, rank.max=5, type="svd", lambda = lambda)
            if(fits$d[1]!=0){
              error[fold.i,lambda.j] <- sum((sim.data$X-fits$u%*%diag(fits$d)%*%t(fits$v))[which(folds==fold.i)]^2)
            } else {
              error[fold.i,lambda.j] <- sum((sim.data$X[which(folds==fold.i)])^2)
            }
          }
        }
      
        #estimating sigma
        fits=softImpute(sim.data$X, rank.max = 5, lambda = lamseq[which.min(colMeans(error))],type="svd")
        deg.free <- sum(fits$d>0)
        if(deg.free!=0){
          sigma.hat <- sum((sim.data$X-fits$u%*%diag(fits$d)%*%t(fits$v))^2)/(n*(p-deg.free))
        } else{
          sigma.hat <- sum((sim.data$X)^2)/(n*(p-deg.free))
        }
      
        command = paste0(pathToMatlab, " -nodisplay -r \"id='", id,
                         "'; maxPC=", maxPC, "; sigma=", sigma.hat, "; run('csv_function'); exit\"")
        system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)
      
        result = read.csv(resFileName)
        if(is.na(result$sigma)){
          tabela[7,rep] <- 1
        } else{
          tabela[7,rep] <- min(result$k,maxPC)
        }

        results[[paste0(SNR, "_", var)]] <- tabela
      }
    }
  }
  return(results)
}

#' Laplace evidence for PCA, as given by Minka
#' 
#' Computes the value of Laplace approximation for given data set and 
#' number of factors.
#' 
#' @param X a matrix with only continuous variables
#' @param k number of principal components fitted
#' @keywords internal
#' @references Automatic choice of dimensionality for PCA, Thomas P. Minka
#' @return L value of Laplace evidence
pca.Laplace <- function(X, k, alfa=1){
  X <- t(X)
  d <- dim(X)[1]
  N <- dim(X)[2]
  m <- d*k - k*(k+1)/2
  
  lambda <- abs(eigen(cov(t(X)), only.values = TRUE)$values)
  v <- sum(lambda[(k+1):d])/(d-k) 
  
  t1 <- -N/2*sum(log(lambda[0:k]))
  t2 <- -N*(d-k)/2*log(v)
  t3 <- -k/2*log(N)
  if(k>0){
    Az <- sapply(1:k, function(i) sum( log(1/lambda[(i+1):d] - 1/lambda[i] ) + log(lambda[i] - lambda[(i+1):d]) + log(N) ))
    if( any(is.nan(Az)) )
      warning(paste("Number of observations ", N, " is to little compared to number of variables ", d, 
                    " to perform a meaningful estimation"))
  } else {
    Az <- 0
  }
  
  t4 <- sum(Az)*(-1)/2
  t5 <- log(2*pi)*(m+k)/2
  t6 <- -k*log(2) + sum( lgamma( (d-1:k+1)/2 ) - (d-1:k+1)/2*log(pi) )
  t1+t2+t3+t4+t5+t6
}

data.simulation.model <- function(n = 100, SNR = 1, numb.vars = 30, k= 2, scale = TRUE){
  sigma <- 1/SNR
  
  X <- NULL
  factors <- NULL
  
  #factors are drawn from normal distribution
  factors <- replicate(k, rnorm(n, 0, 1))
  #coefficients are drawn from centralized (wrt mean) exponential distribution
  coeff <- replicate(numb.vars, rnorm(k, 0, 1))
  SIGNAL <- factors %*% coeff
  SIGNAL <- scale(SIGNAL)
  X <- SIGNAL + replicate(numb.vars, sigma*rnorm(n))
  if(scale){
    X = scale(X)
  }
  return(list(X = X,
              signal = SIGNAL,
              factors = factors))
}

data.simulation.student.noise <- function(SIGNAL=NULL, n = 100, SNR = 1, numb.vars = 30, k = 2, scale = TRUE, df = 2){
  sigma <- 1/SNR
  
  X <- NULL
  
  if(is.null(SIGNAL)){
    factors <- NULL
    
    #factors are drawn from normal distribution
    factors <- replicate(k, rnorm(n, 0, 1))
    #coefficients are drawn from uniform distribution
    coeff <- replicate(numb.vars, rnorm(k, 0, 1))
    SIGNAL <- factors %*% coeff
    
  }
  
  X <- SIGNAL + replicate(numb.vars, sigma*rt(n, df = df)*sqrt((df-2)/df))
  if(scale){
    X = scale(X)
  }
  return(list(X = X,
              signal = SIGNAL))
}


data.simulation.additional.variables <- function(SIGNAL=NULL, n = 100, SNR = 1, numb.vars = 30, k = 2, scale = TRUE, ratio = 0.5){
  sigma <- 1/SNR
  
  X <- NULL
  
  if(is.null(SIGNAL)){
    factors <- NULL
    
    #factors are drawn from normal distribution
    factors <- replicate(k, rnorm(n, 0, 1))
    #coefficients are drawn from uniform distribution
    coeff <- replicate(numb.vars, rnorm(k, 0, 1))
    SIGNAL <- factors %*% coeff
    SIGNAL <- scale(SIGNAL)
  }
  X <- SIGNAL + replicate(numb.vars, sigma*rnorm(n, 0, 1))
  X <- cbind(X, matrix(rnorm(floor(numb.vars*ratio)*n), nrow=n))
  if(scale){
    X = scale(X)
  }
  return(list(X = X,
              signal = SIGNAL))
}

data.simulation.lognormal.noise <- function(SIGNAL=NULL, n = 100, SNR = 1, numb.vars = 30, k = 2, scale = TRUE, mu = 0, sd = 1){
  sigma <- 1/SNR
  
  X <- NULL
  
  if(is.null(SIGNAL)){
    factors <- NULL
    
    #factors are drawn from normal distribution
    factors <- replicate(k, rnorm(n, 0, 1))
    #coefficients are drawn from lognormal distribution
    coeff <- replicate(numb.vars, rnorm(k, 0, 1))
    SIGNAL <- factors %*% coeff
  }
  
  lnsd = sqrt((exp(sd^2)-1)*exp(2*mu + sd^2))
  X <- SIGNAL + replicate(numb.vars, sigma/lnsd*rlnorm(n, meanlog = mu, sdlog = sd))
  if(scale){
    X = scale(X)
  }
  return(list(X = X,
              signal = SIGNAL))
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}