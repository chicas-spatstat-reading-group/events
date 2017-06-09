######################################################
#### Reed Frost Simulation ###########################
sim_RF <- function(N, I_t, R_t, S_t, q) {
  S <- S_t
  I <- I_t
  R <- R_t
  while (I_t > 0) {
    # simulate number of new infectives
    I_new = rbinom(n = 1, size = S_t, prob = 1-q^I_t)
    
    # update number of infectives in each state
    S_t = S_t - I_new # I_new susceptibles are infected
    R_t = R_t + I_t   # previous infectives recover
    I_t = I_new       # I_new infectives for the next gen
    
    # save results
    S <- c(S, S_t)
    I <- c(I, I_t)
    R <- c(R, R_t)
  }
  return(as.matrix(data.frame(S = S, I = I, R = R), drop = F, ncol = 3))
}

rf <- sim_RF(N = 2000, 
             I_t = 2, R_t = 0, S_t = 1998, 
             q = 0.999)

rf

plot(1:length(rf[,"S"]), rf[,"S"], type = "l", 
     ylim = c(0, 1000),
     ylab = "number of individuals",
     xlab = "generation")
lines(1:length(rf[,"I"]), rf[,"I"], col = 2)
lines(1:length(rf[,"R"]), rf[,"R"], col = 3)


set.seed(1)
res <- lapply(1:400, function(i) {
  sim_RF(N = 3, I_t = 1, R_t = 0, S_t = 2,
         q = 0.2)
})

n0 <- n1 <- n2 <- n21 <- 0
count_n <- for(i in 1:400) {
  if (nrow(res[[i]]) == 3) {
    if (isTRUE(all.equal(res[[i]][2, ], c(S = 1, I = 1, R = 1)))) {
      n1 <- n1 + 1
    } else if (isTRUE(all.equal(res[[i]][2, ], c(S = 0, I = 2, R = 1)))) {
      n2 <- n2 + 1
    }
  } else if (nrow(res[[i]]) == 4) {
    n21 <- n21 + 1
  } else if (nrow(res[[i]]) == 2) {
    n0 <- n0 + 1
  }
} 

cat(n0, n1, n2, n21)

