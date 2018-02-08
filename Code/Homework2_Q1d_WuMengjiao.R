x <- c(1.77, -0.23, 2.76, 3.80, 3.47, 56.75, -1.34, 4.24, -2.44,
       3.29, 3.71, -2.40, 4.53, -0.07, -1.05, -13.87, -2.53, -1.75)
n <- length(x)

loglik <- function(theta) {-n*log(pi)-sum(log(1+(x-theta)^2))}


# Calculate theta

start_time <- Sys.time()
mle <- function(theta0, toler=1e-08) {     
  theta <- theta0
  dloglik <- sum(2*(x-theta)/(1+(x-theta)^2))
  i <- 1
  
  while (abs(dloglik) > toler) {
      theta <- theta + dloglik/(n/2)
      dloglik <- sum(2*(x-theta)/(1+(x-theta)^2))
      i <- i+1
      if (i >= 1000) return(theta)
  }
  dloglik <- sum(2*(x-theta)/(1+(x-theta)^2))
  
  if (abs(dloglik) > toler)
    repeat {
      ddloglik <- 2*sum(((x-theta)^2-1)/(1+(x-theta)^2)^2)
      theta <- theta-dloglik/ddloglik
      dloglik <- sum(2*(x-theta)/(1+(x-theta)^2)) 
    }
  theta
}



# Output

IV <- c(-11, -1, 0, 1.5, 4, 4.7, 7, 8, 38)
Result <- as.matrix(sapply(IV, mle))
Result

end_time <- Sys.time()
time_taken <- end_time-start_time
time_taken

