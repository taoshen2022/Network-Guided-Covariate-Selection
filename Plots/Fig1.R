load("./Plots/fig1-left-dcsbm.Rdata")
load("./Plots/fig1-right-dcsbm.Rdata")
load("./Plots/fig1-left-dcsbmm.Rdata")
load("./Plots/fig1-right-dcsbmm.Rdata")
load("./Plots/fig1-left-rgpg.Rdata")
load("./Plots/fig1-right-rdpg.Rdata")

library(latex2exp)
pdf("plot_with_par.pdf", width = 13, height = 8)  
par(mfrow = c(2,3), mar = c(4, 5, 3, 9),                                  # Specify par parameters
    xpd = TRUE)
plot(aseq[5:21], tau_matrix[5:21,1], type = "n", ylim = c(min(tau_matrix[5:21,]), max(tau_matrix[5:21,])), 
     xlab = TeX("Network Assortativity $a$"), ylab = TeX("Aggregated Signal Strength $\\tau$"), cex = 1.5)
lines(aseq[5:21], tau_matrix[5:21,5], col = 1, lwd = 2, lty = 1)
lines(aseq[5:21], tau_matrix[5:21,1], col = 2, lwd = 2)
lines(aseq[5:21], tau_matrix[5:21,2], col = 2, lwd = 2, lty = 2)
lines(aseq[5:21], tau_matrix[5:21,3], col = 2, lwd = 1.5, lty = 3)
lines(aseq[5:21], tau_matrix[5:21,4], col = 2, lwd = 1.5, lty = 4)

#legend(x = max(aseq)-0.15, y = 64, "X only", cex = 0.7, bty = "n", text.col = 3)
#legend(x = max(aseq)-0.17, y = 50, "Y known", cex = 0.7, bty = "n", text.col = 4)
legend("bottomleft", c("Oracle", "NGCS(A,K)", "NGCS(L,K)","NGCS(A,2K)", "NGCS(L,2K)"), 
       lty = c(1,1,2,3,4), col = c(1,2,2,2,2), lwd = c(2,2,2,1.5,1.5), cex = 1)

plot(tseq, tau_matrix22[,1], type = "n", ylim = c(min(tau_matrix22), max(tau_matrix22)), 
     xlab = "Mixture Parameter h", ylab = TeX("Aggregated Signal Strength $\\tau$"), cex=1.5)
lines(tseq, tau_matrix22[,5], col = 1, lwd = 2, lty = 1)
lines(tseq, tau_matrix22[,1], col = 2, lwd = 2)
lines(tseq, tau_matrix22[,2], col = 2, lwd = 2, lty = 2)
lines(tseq, tau_matrix22[,3], col = 2, lwd = 1.5, lty = 3)
lines(tseq, tau_matrix22[,4], col = 2, lwd = 1.5, lty = 4)
legend("bottomleft", c("Oracle", "NGCS(A,K)", "NGCS(L,K)","NGCS(A,2K)", "NGCS(L,2K)"), 
       lty = c(1,1,2,3,4), col = c(1,2,2,2,2), lwd = c(2,2,2,1.5,1.5), cex = 1)

plot(mseq, tau_matrix33[,1], type = "n", ylim = c(min(tau_matrix33), max(tau_matrix33)), 
     xlab = TeX(paste("Signal Strength", "$\\mu$")), ylab = TeX("Aggregated Signal Strength $\\tau$"),cex = 1.5)
lines(mseq, tau_matrix33[,5], col = 1, lwd = 2, lty = 1)
lines(mseq, tau_matrix33[,1], col = 2, lwd = 2)
lines(mseq, tau_matrix33[,2], col = 2, lwd = 2, lty = 2)
lines(mseq, tau_matrix33[,3], col = 2, lwd = 1.5, lty = 3)
lines(mseq, tau_matrix33[,4], col = 2, lwd = 1.5, lty = 4)
#legend(x = max(aseq)-0.15, y = 64, "X only", cex = 0.7, bty = "n", text.col = 3)
#legend(x = max(aseq)-0.17, y = 50, "Y known", cex = 0.7, bty = "n", text.col = 4)
legend("topleft",  c("Oracle", "NGCS(A,K)", "NGCS(L,K)","NGCS(A,2K)", "NGCS(L,2K)"), 
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

TP34 <- TP3
for (i in 1:7){
  TP34[,i] <- 1-TP3[,i]/pseq
}
plot(pseq, TP34[,1], type = "n", ylim = c(min(TP34), max(TP34)), 
     xlab = TeX("Recovered Positives $|\\hat{S}|$"), ylab = TeX("False Discovery Rate $|\\hat{S}\\cap S^c|/|\\hat{S}|$"),cex=1.5)
lines(pseq, TP34[,7], col = 1, lwd = 2)
lines(pseq, TP34[,1], col = 2, lwd = 2)
lines(pseq, TP34[,2], col = 2, lwd = 2, lty = 2)
lines(pseq, TP34[,3], col = 3, lwd = 2)
lines(pseq, TP34[,4], col = 4, lwd = 2)
lines(pseq, TP34[,5], col = 5, lwd = 2)
#lines(pseq, TP3[,6], col = 6, lwd = 2)
abline(v = 53, col = 2, lwd = 1, lty = 2, xpd = FALSE)
points(53, 0.0754717,  lwd = 1, pch = 15, col = 1)
text(140, 0.08, "HC:(53, 0.07)", col = 1)  
legend("bottomright", c("Oracle", "NGCS(A,K)", "NGCS(L,K)", "Chi","FSCA","SKmeans"), 
       col = c(1,2,2,3,4,5), lwd = c(2,2,2,2,2,2), lty = c(1,1,2,1,1,1), cex = 1)


TP24 <- TP2
for (i in 1:5){
  TP24[,i] <- 1-TP2[,i]/pseq44
}
plot(pseq44, TP24[,1], type = "n", ylim = c(min(TP24), max(TP24)), 
     xlab = TeX("Recovered Positives $|\\hat{S}|$"), ylab =TeX("False Discovery Rate $|\\hat{S}\\cap S^c|/|\\hat{S}|$"),cex=1.5)
lines(pseq44, TP24[,5], col = 1, lwd = 2)
lines(pseq44, TP24[,1], col = 2, lwd = 2)
lines(pseq44, TP24[,2], col = 2, lwd = 2, lty=2)
lines(pseq44, TP24[,3], col = 6, lwd = 2)
lines(pseq44, TP24[,4], col = 7, lwd = 2)
abline(v = 50, col = 2, lwd = 1, lty = 2, xpd = FALSE)
points(50, 0.026, lwd = 2, pch = 15,col=1)
text(30, 0.05, "HC:(50, 0.02)",col=1)  
legend("bottomright", c("Oracle", "NGCS(A,K)", "NGCS(L,K)", "Marginal", "RF"), 
       lty = c(1,1,2,1,1), col = c(1,2,2,6,7), lwd = c(2,2,2,2,2), cex = 1)
dev.off() 
