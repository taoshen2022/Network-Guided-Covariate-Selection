library(latex2exp)
load("./Visualization/Sfig2-left-dcsbm.Rdata")
load("./Visualization/Sfig2-right-dcsbm.Rdata")

pdf("dcsbm.pdf", width = 15, height = 6)  
par(mfrow = c(1,2), mar = c(4, 5, 3, 9),                                  # Specify par parameters
    xpd = TRUE)

plot(aseq[5:21], tau_matrix[5:21,1], type = "n", ylim = c(min(tau_matrix[5:21,]), max(tau_matrix[5:21,])), 
     xlab = TeX("Network Assortativity $a$"), ylab = TeX("Aggregated Signal Strength $\\tau$"), cex = 1.5)
lines(aseq[5:21], tau_matrix[5:21,5], col = 1, lwd = 2, lty = 1)
lines(aseq[5:21], tau_matrix[5:21,1], col = 2, lwd = 2)
lines(aseq[5:21], tau_matrix[5:21,2], col = 2, lwd = 2, lty = 2)
lines(aseq[5:21], tau_matrix[5:21,3], col = 2, lwd = 1.5, lty = 3)
lines(aseq[5:21], tau_matrix[5:21,4], col = 2, lwd = 1.5, lty = 4)
legend("bottomleft", c("Oracle", "NGCS(A,K)", "NGCS(L,K)","NGCS(A,2K)", "NGCS(L,2K)"), 
       lty = c(1,1,2,3,4), col = c(1,2,2,2,2), lwd = c(2,2,2,1.5,1.5), cex = 1)


TP334 <- TP33
for (i in 1:6){
  TP334[,i] <- 1-TP33[,i]/pseq22
}
plot(pseq22, TP334[,1], type = "n", ylim = c(min(TP334[,-7]), max(TP334[,-7])), 
     xlab = TeX("Recovered Positives $|\\hat{S}|$"), ylab = TeX("False Discovery Rate $|\\hat{S}\\cap S^c|/|\\hat{S}|$"),cex=1.5)
lines(pseq22, TP334[,6], col = 1, lwd = 2)
lines(pseq22, TP334[,1], col = 2, lwd = 2)
lines(pseq22, TP334[,2], col = 2, lwd = 2, lty = 2)
lines(pseq22, TP334[,3], col = 3, lwd = 2)
lines(pseq22, TP334[,4], col = 4, lwd = 2)
lines(pseq22, TP334[,5], col = 5, lwd = 2)
abline(v = 49, col = 2, lwd = 1, lty = 2, xpd = FALSE)
points(49, 0.01020408,  lwd = 2, pch = 15, col =1)
text(70, 0.02, "HC:(49, 0.01)",col=1)  
legend("topleft", c("Oracle", "NGCS(A,K)", "NGCS(L,K)", "Chi", "FSCA", "SKmeans"), 
       col = c(1,2,2,3,4,5), lwd = c(2,2,2,2,2,2), lty = c(1,1,2,1,1,1), cex = 1)


dev.off() 