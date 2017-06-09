######################################################
#### Reed Frost Model Inference ######################

rf_fit <- function(n0, n1, n2, niters, inits, priors) {
  q <- numeric(niters)
  q[1] <- inits$q
  n21 <- numeric(niters)
  n21[1] <- inits$n21
  for(i in 1:niters) {
    q[i + 1] <- rbeta(n = 1,
                   2 * (n0 + n1) + n21[i] + priors$alpha,
                   n1 + 2 * n2 + priors$delta)
    n21[i + 1] <- rbinom(n = 1, n2, 2 * q[i] / (2 * q[i] + 1))
  }
  return(data.frame(q = q, n21 = n21))
}

rf_fit1 <- function(n0, n1, n2, niters, inits, priors) {
  q <- inits$q
  n21 <- inits$n21
  out <- rbind(c(), c(q = q, n21 = n21))
  for(i in 1:niters) {
    q <- rbeta(n = 1,
               2 * (n0 + n1) + n21 + priors$alpha,
               n1 + 2 * n2 + priors$delta)
    n21 <- rbinom(n = 1, n2, 
                  2 * q / (2 * q + 1))
    out <- rbind(out, c(q, n21))
  }
  return(as.data.frame(out))
}

priors = list(alpha = 0.01, delta = 0.01)
inits = list(q = 0.2, n21 = 104)
df <- rf_fit(n0 = 14, n1 = 25, n2 = 104 + 257, 10000, inits, priors)
df1 <- rf_fit1(n0 = 14, n1 = 25, n2 = 104 + 257, 10000, inits, priors)

par(mfrow = c(2, 2))
plot(df$q, type = "l", xlab = "iter", ylab = "q")
abline(h = inits$q, col = "red")
plot(df$n21, type = "l", xlab = "iter", ylab = "n21")
abline(h = inits$n21, col = "red")

plot(df1$q, type = "l", xlab = "iter", ylab = "q")
abline(h = inits$q, col = "red")
plot(df1$n21, type = "l", xlab = "iter", ylab = "n21")
abline(h = inits$n21, col = "red")

pdf("rf_results.pdf", width = 8, height = 5)
par(mfrow = c(2, 2), mar=c(4,4,1,1)+0.1)
n_sim <- length(df$q)
n_mid <- floor(n_sim / 2)

plot(ecdf(df$q[1:n_mid]), main = "", xlab = "q")
lines(ecdf(df$q[(n_mid + 1):n_sim]), col = "blue", lty = "dashed")

plot(ecdf(df$n21[1:n_mid]), main = "", xlab = "n21", pch = ".")
lines(ecdf(df$n21[(n_mid + 1):n_sim]), col = "blue", pch = ".")

plot(density(df$q), xlab = "q", ylab = "density", main = "")
abline(v = inits$q, col = "red")

hist(df$n21, xlab = "n21", main = "")
abline(v = inits$n21, col = "red")
dev.off()