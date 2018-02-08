x <- c(1.77, -0.23, 2.76, 3.80, 3.47, 56.75, -1.34, 4.24, -2.44,
       3.29, 3.71, -2.40, 4.53, -0.07, -1.05, -13.87, -2.53, -1.75)
n <- length(x)

loglik <- function(theta) {-n*log(pi)-sum (log(1+(x-theta)^2))}



# Calculate

start_time <- Sys.time()
mle <- function(theta0, alpha, toler=1e-08) {     
  theta <- theta0
  dloglik <- sum(2*(x-theta)/(1+(x-theta)^2))
  i <- 1 
  
  while (abs(dloglik) > toler) {
    ddloglik <- 2*sum(((x-theta)^2-1)/(1+(x-theta)^2)^2)
    theta <- theta+dloglik*alpha
    dloglik <- sum(2*(x-theta)/(1+(x-theta)^2))
    i <- i+1
    if (i >= 1000) return("Maximum iterations reached")
  }
  theta
}



# Output

Result <- matrix(NA, 9, 3)
IV <- as.matrix(c(-11, -1, 0, 1.5, 4, 4.7, 7, 8, 38))
alpha_set <- t(as.matrix(c(1, 0.64, 0.25)))
for (j in 1:3) {
  for (m in 1:9){
    Result[m,j] <- mle(IV[m,1], alpha_set[1,j])
  }
}
Result

end_time <- Sys.time()
time_taken <- end_time-start_time
time_taken

