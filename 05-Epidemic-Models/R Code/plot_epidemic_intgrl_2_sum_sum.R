I_times <- c(0, 3.5, 2.5, 4, -1)
R_times <- c(7, 8, 6, 11, 3)
E_times <- c(I_times, R_times)[order(c(I_times, R_times))]
E_labels <- c(paste0("Y", 1:5), paste0("O", 1:5))[order(c(I_times, R_times))]
beta <- 0.2
ip_1 <- c(beta, beta)
t_1 <- c(E_times[1:2])
ip_2 <- c(beta, 2*beta, 3*beta, 2*beta, 2*beta)
t_2 <- c(E_times[1:5])
ip_3 <- c(beta, 2*beta, 2*beta)
t_3 <- c(E_times[1:3])
ip_4 <- c(beta, 2*beta, 3*beta, 2*beta, 3*beta, 3*beta)
t_4 <- c(E_times[1:6])

pdf("sumsum_IP.pdf", width = 6, height = 7)
par(mfrow = c(5, 1), mar=c(4, 4, 2, 0.5))
plot(0, type="n", xlim=c(-1,11), ylim=c(1,6), xlab="Time", ylab="ID", main="Epidemic Progression: beta = 0.2")
for(x in 1:length(I_times)) lines(c(I_times[x], R_times[x]), c(x, x))
for(x in 1:(length(I_times)*2)) lines(c(E_times[x], E_times[x]), c(-1, 6), lty = 2, col = "grey")
text(x = E_times, y = 5.5, labels = E_labels, col = "blue")
par(mar=c(2, 4, 2, 0.5))
## infectious pressure on each individual
plot(0, type="n", xlim=c(-1,11), ylim=c(0,0.6), xlab="Time", ylab="IP", main="Infectious pressure on ID 1")
lines(ip_1 ~ t_1, type = "s")
for(x in 1:(length(I_times)*2)) lines(c(E_times[x], E_times[x]), c(-1, 6), lty = 2, col = "grey")
text(x = E_times, y = 0.1, labels = E_labels, col = "blue")

plot(0, type="n", xlim=c(-1,11), ylim=c(0,0.6), xlab="Time", ylab="IP", main="Infectious pressure on ID 2")
lines(ip_2 ~ t_2, type = "s")
for(x in 1:(length(I_times)*2)) lines(c(E_times[x], E_times[x]), c(-1, 6), lty = 2, col = "grey")
text(x = E_times, y = 0.1, labels = E_labels, col = "blue")

plot(0, type="n", xlim=c(-1,11), ylim=c(0,0.6), xlab="Time", ylab="IP", main="Infectious pressure on ID 3")
lines(ip_3 ~ t_3, type = "s")
for(x in 1:(length(I_times)*2)) lines(c(E_times[x], E_times[x]), c(-1, 6), lty = 2, col = "grey")
text(x = E_times, y = 0.1, labels = E_labels, col = "blue")

plot(0, type="n", xlim=c(-1,11), ylim=c(0,0.6), xlab="Time", ylab="IP", main="Infectious pressure on ID 4")
lines(ip_4 ~ t_4, type = "s")
for(x in 1:(length(I_times)*2)) lines(c(E_times[x], E_times[x]), c(-1, 6), lty = 2, col = "grey")
text(x = E_times, y = 0.1, labels = E_labels, col = "blue")
dev.off()