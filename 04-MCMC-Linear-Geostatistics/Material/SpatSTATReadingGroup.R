rm(list = ls())
################################################################################
#Metropolis within Gibbs for Linear geostatistical model
###############################################

n <- 50
xy <- cbind(runif(n, min = 0, max = 10), runif(n, min = 0, max = 6))
require(geoR)
U <- as.matrix(dist(xy))
sigma2 <- 2.5
phi <- 1.5
tau <- 0.3
V <- sigma2*exp(-U/phi) + tau*diag(n)
x <- rnorm(n, 2, 0.5)
X <- cbind(1, x)
p <- ncol(X)
beta0 <- 2
beta1 <- 0.5
require(MASS)
S <- mvrnorm(mu= rep(0, n),  Sigma = V)
y <- beta0 + beta1*x + S


loglikelihood<-function(theta){
  beta <- theta[1:p]
  sigma2 <- exp(theta[p+1])
  tau <- exp(theta[p+2])
  phi <- exp(theta[p+3])
  V <- exp(-U/phi) + (tau/ sigma2) * diag(n)
  logl <- -0.5*(n*log(2*pi) + n*log(sigma2) + as.numeric(determinant(V)$modulus) + 
    (1/sigma2)*(t(y-X%*%beta)%*%solve(V)%*%(y-X%*%beta)))
  return(logl)
}

require(pscl)
sigma.logprior = function(theta){
  p <- ncol(X)
  sigma2 = theta[p+1]
  sigma.prior = densigamma(exp(sigma2), 2, 1)
  return(log(sigma.prior))
}

tau.logprior = function(theta){
  p <- ncol(X)
  tau2=theta[p+2]
  tau.prior= densigamma(exp(tau2), 2, 1)
  return(log(tau.prior))
}

phi.logprior = function(theta){
  p <- ncol(X)
  phi=theta[p+3]
  phi.prior= densigamma(exp(phi), 2, 1)
  return(log(phi.prior))
}

# phi.logprior = function(theta){
#   p <- ncol(X)
#   phi=theta[p+3]
#   phi.prior= dunif(exp(phi), min=1, max=2)
#   return(log(phi.prior))
# }

beta.conditional.function = function(theta){
  beta=theta[1:p]
  sigma2=exp(theta[p+1])
  tau=exp(theta[p+2])
  phi <- exp(theta[p+3])
  V=sigma2 *exp(-U/phi) + tau * diag(n)
  beta.cond.Sigma=solve(t(X)%*% (solve(V))%*%X + diag(p))
  beta.cond.mu= beta.cond.Sigma%*%(t(X)%*%(solve(V))%*%y)
  return(mvrnorm(mu = beta.cond.mu, Sigma = beta.cond.Sigma))
}

theta.proposalfunction = function(theta){
  p <- ncol(X)
  return(mvrnorm(mu = theta[(p+1):(p+3)], Sigma = c(1, 1, 1) * diag(3)))
}

run_metropolis_MCMC = function(startvalue, iterations, burn.in, thinning){
  p <- ncol(X)
  chain = matrix(NA, nrow=iterations+1, ncol=length(startvalue))
  #startvalue[(p+1):(p+3)] <- exp(startvalue[(p+1):(p+3)]) 
  chain[1,] = startvalue
  sigma.acc <- 0
  tau.acc <- 0
  phi.acc <- 0
  for (i in 1:iterations){
    proposal=vector()
    proposal[1:p] = beta.conditional.function(chain[i,])
    proposal[(p+1):(p+3)] = theta.proposalfunction(chain[i,])
    print(proposal)
    chain[i+1,][1:p]=proposal[1:p]
    sigma.probab = loglikelihood(proposal)+ sigma.logprior(proposal) + proposal[(p+1)] - 
                         loglikelihood(chain[i,]) - sigma.logprior(chain[i,])- chain[i,][(p+1)] 
    print(sigma.probab)
    #Gibbs sampler for the beta 
    if (log(runif(1)) < sigma.probab){
      chain[i+1,][(p+1)] = proposal[(p+1)]
      sigma.acc = sigma.acc +1
    }
    else{
      chain[i+1,][(p+1)] = chain[i,][(p+1)]
    }
    tau.probab = loglikelihood(proposal)+ tau.logprior(proposal) + proposal[(p+2)] - 
                       loglikelihood(chain[i,]) - tau.logprior(chain[i,]) - chain[i,][(p+2)] 
    if (log(runif(1)) < tau.probab){
      chain[i+1,][(p+2)] = proposal[(p+2)]
      tau.acc <- tau.acc +1
    }
    else{
      chain[i+1,][(p+2)] = chain[i,][(p+2)]
    }

    phi.probab = loglikelihood(proposal)+ phi.logprior(proposal) + proposal[(p+3)] - 
                       loglikelihood(chain[i,]) - phi.logprior(chain[i,]) - chain[i,][(p+3)]
    if (log(runif(1)) < phi.probab){
      chain[i+1,][(p+3)] = proposal[(p+3)]
      phi.acc <- phi.acc +1
    }
    else{
      chain[i+1,][(p+3)] = chain[i,][(p+3)]
    }
  }
  # burn-in and thining the sample
  chain = apply(chain, 2, function(chain) chain[seq(from=burn.in,to=iterations,by=thinning)])
  require(coda)
  return(list(chain=mcmc(chain), sigma.acc= sigma.acc/iterations, 
              tau.acc=tau.acc/iterations, phi.acc=phi.acc/iterations)) 
}


startvalue = c(0.1, 0.1, log(c(2.6, 0.4, 1.6)))
itermax=10000
burn.in = 2000
thinning = 8
run = run_metropolis_MCMC(startvalue, itermax, burn.in, thinning)
chain <- run$chain
summary(chain)
plot(chain)
plot(exp(chain[,3]))
pairs(data.frame(chain))

#############
chain[, (3:5)] <- exp(chain[, (3:5)])

plot(acf(chain[, 5])$lag, acf(chain[, 5])$acf, type = "l", xlab = "lag", 
     ylab = "autocorrelation", ylim = c(-0.1, 1), 
     main = "Autocorrelogram of the simulated samples" )
abline(h = 0, lty = "dashed", col = 2)




