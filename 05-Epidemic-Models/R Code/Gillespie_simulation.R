######################################################
#### Gillespie Simulation ############################

##
# a is the vector of rates of events
## simulate time to next event 
# (time to next event, tau, is epxonentially distributed 
# with rate equa to the sum of rates for all possible events)
# f(t)=(sum(a)exp(-tau*sum(a)))
## Simulate event type
# P(Event = v) = a[v]/sum(a)

S = 999
I = 1
R = 0
N = 1000

sim_glsp <- function(S, I, R, N, beta, gamma) {
  transitions <- rbind(c(-1, 1, 0), c(0, -1, 1))
  tau <- time <- 0
  out <- rbind(c(), c(event=1, time=time, S=S, I=I, R=R))
  while (I > 0) {
    rates <- c(beta * I * S,# / N, # rate infection
               gamma * I)          # rate removal
    tau <- rexp(1, sum(rates)) # Simulate next time and event type
    event_type <- sample(1:2, 1, prob = rates) #/sum(rates))
    S <- S + transitions[event_type, 1] # Update states
    I <- I + transitions[event_type, 2]
    R <- R + transitions[event_type, 3]
    time <- time + tau
    out <- rbind(out, c(event_type, time, S, I, R))
  }
  return(as.data.frame(out))
} 

glsp <- sim_glsp(S, I, R, beta = 0.75/N, gamma = 0.3)
repeat {
  u <- sample(1:1000, 1)
  print(u)
  set.seed(u)
  glsp <- sim_glsp(S=29,I=1,R=0, N=30, 0.2, 0.1)
  if (dim(glsp)[1] >= 20) break
}

set.seed(186)
glsp <- sim_glsp(S=29,I=1,R=0, N = 30, 0.2/30, 0.1)
glsp$time_0 <- glsp$time - glsp$time[which(glsp$event == 2)[1]]
paste(round(subset(glsp, event == 1)$time_0, 1), collapse = ", ")
paste(round(subset(glsp, event == 2)$time_0, 1), collapse = ", ")

pdf("Gillespie_sim_data.pdf", width = 8, height = 2.4)
par(mar=c(3,3,1,1)+0.1, mgp=c(2,1,0))
plot(glsp$time, glsp$S, type = "l", 
     ylim = c(0, 30),
     ylab = "number of individuals",
     xlab = "time")
lines(glsp$time, glsp$I, col = 2)
lines(glsp$time, glsp$R, col = 3)
dev.off()

plot(glsp$time_0, glsp$R, type = "s",
     ylab = "Number of removed individuals", xlab = "Time")

pdf("int_infected.pdf", width = 4, height = 3)
par(mar=c(3,3,1,1)+0.1, mgp=c(2,1,0))
plot(glsp$time_0, glsp$I, type = "s",
     ylab = "No. infected", xlab = "Time")
dev.off()

plot(glsp$time_0, glsp$S, type = "s",
     ylab = "Number of susceptible individuals", xlab = "Time")

plot(glsp$time_0, glsp$I * glsp$S, type = "s",
     ylab = "S*I", xlab = "Time")
