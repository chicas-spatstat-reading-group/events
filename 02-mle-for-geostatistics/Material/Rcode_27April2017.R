rm(list=ls())
library(geoR)

data("elevation")

D <- cbind(1,elevation$coords[,1],elevation$coords[,2],
           elevation$coords[,1]^2,elevation$coords[,2]^2,
           elevation$coords[,1]*elevation$coords[,2])
p <- ncol(D)
y <- elevation$data
n <- length(y)  
coords <- elevation$coords
  
U <- dist(coords)
  
log.lik <- function(par) {
  beta <- par[1:p]
  sigma2 <- exp(par[p+1])
  phi <- exp(par[p+2])
  nu2 <- exp(par[p+3])
  
  R <- varcov.spatial(dists.lowertri=U,kappa=kappa,
                      cov.pars=c(1,phi),nugget=nu2)$varcov
  R.inv <- solve(R)
  ldet.R <- determinant(R)$modulus
  mu <- as.numeric(D%*%beta)
  diff.y <- y-mu
  out <- -0.5*(n*log(sigma2)+ldet.R+
         t(diff.y)%*%R.inv%*%(diff.y)/sigma2)
  as.numeric(out)
}

data.aux <- data.frame(y,D[,-1])
lm.fit <- lm(y~.,data=data.aux)
beta.start <- coef(lm.fit)   
sigma2.start <- mean(lm.fit$residuals^2)


plot(variog(coords=coords,data=lm.fit$residuals))

tau2.start <- 10
nu2.start <- sigma2.start

curve(tau2.start+sigma2.start*(1-exp(-x/1)),add=TRUE)
phi.start <- 1

par.start <- c(beta.start,log(c(sigma2.start,phi.start,tau2.start/sigma2.start)))

matern.grad.phi <- function(U,phi,kappa) {
  n <- attr(U,"Size")
  grad.phi.mat <- matrix(NA,nrow=n,ncol=n)
  ind <- lower.tri(grad.phi.mat)
  grad.phi <- der.phi(as.numeric(U),phi,kappa)
  grad.phi.mat[ind] <-  grad.phi
  grad.phi.mat <- t(grad.phi.mat)
  grad.phi.mat[ind] <-  grad.phi
  diag(grad.phi.mat) <- rep(der.phi(0,phi,kappa),n)
  grad.phi.mat
}

der.phi <- function(u,phi,kappa) {
  u <- u+10e-16
  if(kappa==0.5) {
    out <- (u*exp(-u/phi))/phi^2
  } else {
    out <- ((besselK(u/phi,kappa+1)+besselK(u/phi,kappa-1))*
              phi^(-kappa-2)*u^(kappa+1))/(2^kappa*gamma(kappa))-
      (kappa*2^(1-kappa)*besselK(u/phi,kappa)*phi^(-kappa-1)*
         u^kappa)/gamma(kappa)
  }
  out
}

grad.log.lik <- function(par) {
  beta <- par[1:p]
  sigma2 <- exp(par[p+1])
  phi <- exp(par[p+2])
  nu2 <- exp(par[p+3])
  R <- varcov.spatial(dists.lowertri=U,cov.model="matern",
                      cov.pars=c(1,phi),
                      nugget=nu2,kappa=kappa)$varcov
  
  # Derivative of the paer
  R.inv <- solve(R)
  R1.phi <- matern.grad.phi(U,phi,kappa)
  m1.phi <- R.inv%*%R1.phi
  t1.phi <- -0.5*sum(diag(m1.phi))
  m2.phi <- m1.phi%*%R.inv; rm(m1.phi)
      

  t1.nu2 <- -0.5*sum(diag(R.inv))
  m2.nu2 <- R.inv%*%R.inv

  mu <- D%*%beta
  diff.y <- y-mu
  q.f <- t(diff.y)%*%R.inv%*%diff.y
  grad.beta <-  t(D)%*%R.inv%*%(diff.y)/sigma2
  grad.log.sigma2 <- (-n/(2*sigma2)+0.5*q.f/(sigma2^2))*sigma2
  grad.log.phi <- (t1.phi+0.5*as.numeric(t(diff.y)%*%m2.phi%*%(diff.y))/sigma2)*phi
  grad.log.nu2 <-  (t1.nu2+0.5*as.numeric(t(diff.y)%*%m2.nu2%*%(diff.y))/sigma2)*nu2
        
  out <- c(grad.beta,grad.log.sigma2,grad.log.phi,grad.log.nu2)
  return(out)
}
    
library(maxLik)

kappa <- 1.5

# How to chek if the derivates are ok

check <- compareDerivatives(
  log.lik,
  grad.log.lik,
  t0=par.start
)
    
der2.phi <- function(u,phi,kappa) {
  u <- u+10e-16
  if(kappa==0.5) {
    out <- (u*(u-2*phi)*exp(-u/phi))/phi^4
  } else {
    bk <- besselK(u/phi,kappa)
    bk.p1 <- besselK(u/phi,kappa+1)
    bk.p2 <- besselK(u/phi,kappa+2)
    bk.m1 <- besselK(u/phi,kappa-1)
    bk.m2 <- besselK(u/phi,kappa-2)
    out <- (2^(-kappa-1)*phi^(-kappa-4)*u^kappa*(bk.p2*u^2+2*bk*u^2+
           bk.m2*u^2-4*kappa*bk.p1*phi*u-4*
           bk.p1*phi*u-4*kappa*bk.m1*phi*u-4*bk.m1*phi*u+
           4*kappa^2*bk*phi^2+4*kappa*bk*phi^2))/(gamma(kappa))
  }
  out
}

matern.hessian.phi <- function(U,phi,kappa) {
  n <- attr(U,"Size")
  hess.phi.mat <- matrix(NA,nrow=n,ncol=n)
  ind <- lower.tri(hess.phi.mat)
  hess.phi <- der2.phi(as.numeric(U),phi,kappa)
  hess.phi.mat[ind] <-  hess.phi
  hess.phi.mat <- t(hess.phi.mat)
  hess.phi.mat[ind] <-  hess.phi
  diag(hess.phi.mat) <- rep(der2.phi(0,phi,kappa),n)
  hess.phi.mat
}

hess.log.lik <- function(par) {
  beta <- par[1:p]
  sigma2 <- exp(par[p+1])
  phi <- exp(par[p+2])
  mu <- D%*%beta
  nu2 <- exp(par[p+3])
  
  R <- varcov.spatial(dists.lowertri=U,cov.model="matern",
                      cov.pars=c(1,phi),
                      nugget=nu2,kappa=kappa)$varcov
      
  R.inv <- solve(R)
  R1.phi <- matern.grad.phi(U,phi,kappa)
  m1.phi <- R.inv%*%R1.phi
  t1.phi <- -0.5*sum(diag(m1.phi))
  m2.phi <- m1.phi%*%R.inv; rm(m1.phi)
      
  t1.nu2 <- -0.5*sum(diag(R.inv))
  m2.nu2 <- R.inv%*%R.inv
  t2.nu2 <- 0.5*sum(diag(m2.nu2))
  n2.nu2 <- 2*R.inv%*%m2.nu2
  t2.nu2.phi <- 0.5*sum(diag(R.inv%*%R1.phi%*%R.inv))
  n2.nu2.phi <- R.inv%*%(R.inv%*%R1.phi+
                         R1.phi%*%R.inv)%*%R.inv

      
  R2.phi <- matern.hessian.phi(U,phi,kappa)
  t2.phi <- -0.5*sum(diag(R.inv%*%R2.phi-R.inv%*%R1.phi%*%R.inv%*%R1.phi))
  n2.phi <- R.inv%*%(2*R1.phi%*%R.inv%*%R1.phi-R2.phi)%*%R.inv
    
  diff.y <- y-mu
  q.f <- t(diff.y)%*%R.inv%*%diff.y

  grad.log.sigma2 <- (-n/(2*sigma2)+0.5*q.f/(sigma2^2))*sigma2
  grad.log.phi <- (t1.phi+0.5*as.numeric(t(diff.y)%*%m2.phi%*%(diff.y))/sigma2)*phi
  grad.log.nu2 <-  (t1.nu2+0.5*as.numeric(t(diff.y)%*%m2.nu2%*%(diff.y))/sigma2)*nu2

      
  H <- matrix(0,nrow=length(par),ncol=length(par))
  H[1:p,1:p] <- -t(D)%*%R.inv%*%D/sigma2
  H[1:p,p+1] <- H[p+1,1:p] <- -t(D)%*%R.inv%*%(diff.y)/sigma2
  H[1:p,p+2] <- H[p+2,1:p] <- -phi*as.numeric(t(D)%*%m2.phi%*%(diff.y))/sigma2
  H[p+1,p+1] <- (n/(2*sigma2^2)-q.f/(sigma2^3))*sigma2^2+
                grad.log.sigma2
      
  H[p+1,p+2] <- H[p+2,p+1] <- (grad.log.phi/phi-t1.phi)*(-phi)
      
  H[p+2,p+2] <- (t2.phi-0.5*t(diff.y)%*%n2.phi%*%(diff.y)/sigma2)*phi^2+
                grad.log.phi
      
  H[1:p,p+3] <- H[p+3,1:p] <- -nu2*as.numeric(t(D)%*%m2.nu2%*%(diff.y))/sigma2
  H[p+2,p+3] <- H[p+3,p+2] <- (t2.nu2.phi-0.5*t(diff.y)%*%n2.nu2.phi%*%(diff.y)/sigma2)*phi*nu2
  H[p+1,p+3] <- H[p+3,p+1] <- (grad.log.nu2/nu2-t1.nu2)*(-nu2)
  H[p+3,p+3] <- (t2.nu2-0.5*t(diff.y)%*%n2.nu2%*%(diff.y)/sigma2)*nu2^2+
                grad.log.nu2
      
  return(H)
}
    

check <- compareDerivatives(
  grad.log.lik,
  hess.log.lik,
  t0=par.start+runif(p+3,-1,1)
)

estim <- nlminb(par.start,
                function(x) -log.lik(x),
                function(x) -grad.log.lik(x),
                function(x) -hess.log.lik(x),
                control=list(trace=1))
estim$gradient <- grad.log.lik(estim$par)
estim$hessian <- hess.log.lik(estim$par)
