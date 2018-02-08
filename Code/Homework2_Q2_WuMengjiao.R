x <- c(3.91, 4.85, 2.28, 4.06, 3.70, 4.04, 5.46, 3.53, 2.28, 1.96,
       2.53, 3.88, 2.22, 3.47, 4.82, 2.46, 2.99, 2.54, 0.52)

loglik <- function(theta) sum(log(1-cos(x-theta)))-log(2*pi)



# Part (a)
S <- seq(-pi, pi, len=1000)
loglik_S <- sapply(S, loglik)
plot(S, loglik_S, type ="l", xlab=expression(theta), ylab=expression(logL(theta)))



# Part (b)
mme <- asin(mean(x)-pi)



# Part (c) and (d)

mle <- function(theta0, toler=1e-08) {     
  theta <- theta0
  dloglik <- -sum(sin(x-theta)/(1-cos(x-theta)))
  
  while (abs(dloglik) > toler) {
    ddloglik <- -sum(1/(1-cos(x-theta)))
    theta <- theta-dloglik/ddloglik
    dloglik <- -sum(sin(x-theta)/(1-cos(x-theta)))
  }
  theta
}

Result <- as.matrix(c(mle(mme), mle(-2.7), mle(2.7)))
Result


# Parts (e)
M <- seq(-pi, pi, len=200)
N <- function(theta) round(mle(theta), 4)
theta_value <- sapply(M, N)
soln <- split(M, theta_value)
min_value <- sapply(soln, min)
max_value <- sapply(soln, max)
Result_partition <- cbind(round(min_value, 5), round(max_value, 5))

