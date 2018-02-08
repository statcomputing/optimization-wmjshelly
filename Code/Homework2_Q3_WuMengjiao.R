D <- c(0, 8, 28, 41, 63, 69, 97, 117, 135, 154)
B <- c(2, 47, 192, 256, 768, 896, 120, 896, 1184, 1024)
N0 <- B[1]
n <- length(B)
library(MASS)
library(bbmle)


# Part (a)

gauss_newton <- function(y, f, df, x, eps=1e-08, maxiter=1000, ...) {
  
  fx <- f(x, ...)
  dfx <- df(x, ...)
  i <- 0
  repeat {
    i <- i+1
    x_new <- x+ginv(t(dfx) %*% dfx) %*% t(dfx) %*% (y-fx)
    if(mean(abs(x_new-x)) < eps || i >= maxiter) {
      if(i >= maxiter) warning("Maximum number of iterations reached")
      break
    }
    x <- x_new
    fx <- f(x_new, ...)
    dfx <- df(x_new, ...)
  }
  return(list(x=as.numeric(x_new), fx=f(x_new, ...)))
}

f <- function(u) {
  K <- u[1]
  r <- u[2]
  return(K*N0/(N0+(K-N0)*exp(-r*D)))
}
df <- function(u) {
  fu <- f(u)
  K <- u[1]
  r <- u[2]
  o1 <- fu*(1-fu*exp(-r*D))/K
  o2 <- fu**2*(K-N0)*D*exp(-r*D)/K/N0
  return(matrix(c(o1, o2), nrow=n, ncol=2, byrow=FALSE))
}
o_gn <- gauss_newton(B, f, df, c(800, 0.1))
print(o_gn)



# Part (b)

g <- function(t) {o_gn$x[1] * N0 / (N0 + (o_gn$x[1] - N0) * exp(-o_gn$x[2] * t))}
k <- seq(750, 900, by=5)
r <- seq(0.05, 0.2, by=0.025)
Q <- 0 * outer(k, r)
for(i in seq_along(k)) {
  for(j in seq_along(r)) {
    Q[i, j] <- -sum((B - f(c(k[i], r[j])))**2)
    }
}
plot(x=0, y=0, xlim=range(k), ylim=range(r), xlab="K", ylab="r")
contour(k, r, Q, lty="solid", add=TRUE, nlevels=15)



# Part (c)

fit <- mle2(log(N)~dnorm(mean=log((K*N0)/(N0+(K-N0)*exp(-r*t))), sd=sqrt(var)), 
            start=list(K=500, r=0.5, var=0.5), data=list(t=D, N=B))
coef(summary(fit)) 


