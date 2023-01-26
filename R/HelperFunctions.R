rmnvAnder <-  function(mean,sigma){
  drop(mean + t(chol(sigma))%*%rnorm(length(mean)))
}

# this leaves off the 2pi and evaluates it at 0
dmvn0Ander <- function(mean,sigma){
  Ch <- chol(sigma)
  return(-sum(log(diag(Ch))) -sum(backsolve(Ch,mean, transpose = TRUE)^2)/2)
}

# this leaves off the 2pi and evaluates it at 0
# accepts the cholesky decomp of the variance matrix
dmvn0ChAnder <- function(mean,Ch){
  return(-sum(log(diag(Ch))) -sum(backsolve(Ch,mean, transpose = TRUE)^2)/2)
}

