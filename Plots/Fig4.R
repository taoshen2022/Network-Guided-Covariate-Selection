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

TP334 <- TP33
for (i in 1:6){
  TP334[,i] <- 1-TP33[,i]/pseq22
}
plot(pseq22, TP334[,1], type = "n", ylim = c(min(TP334[,-7]), 1.1), 
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

plot(pseqq, error[,1], type = "n", ylim = c(min(error), 0.9), 
     xlab = "ln(p)", ylab = "Clustering Error Rate", cex = 1.5)
lines(pseqq, error[,2], col = 1, lwd = 2)
#lines(pseqq, error[,3], col = 1, lwd = 2, lty = 2)
lines(pseqq, error[,4], col = 2, lwd = 2)
lines(pseqq, error[,1], col = 3, lwd = 2)
lines(pseqq, error[,5], col = 4, lwd = 2)
#lines(ppseq, error[,6], col = 6, lwd = 2)
lines(pseqq, error[,7], col = "pink", lwd = 2)
#legend(x = max(aseq)-0.15, y = 64, "X only", cex = 0.7, bty = "n", text.col = 3)
#legend(x = max(aseq)-0.17, y = 50, "Y known", cex = 0.7, bty = "n", text.col = 4)
legend("topright", c("NG-Clu", "IF-PCA", "Spec", "SKmeans", "SAS"), 
       lty = c(1,1,1,1,1), col = c(1,2,3,4,"pink"), lwd = c(2,2,2,2,2), cex = 1)

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

plot(museq, mse[,1], type = "n", ylim = c(min(mse[,-4], na.rm = TRUE), 1.8), 
     xlab = TeX("Signal Strength $\\mu$"), ylab = "Mean Squared Error", cex = 1.5)
lines(museq, mse[,2], col = 1, lwd = 2)
lines(museq, mse[,1], col = 2, lwd = 2)
#lines(mseq, mse[,3], col = 3, lwd = 2)
#lines(museq, mse[,4], col = 4, lwd = 2)
lines(museq, mse[,5], col = 5, lwd = 2)
lines(museq, mse[,6], col = 6, lwd = 2)
#lines(museq, mse[,7], col = 6, lwd = 2)
lines(museq, mse[,8], col = 7, lwd = 2)
lines(museq, mse[,9], col = 8, lwd = 2)
#lines(museq, mse[,10], col = 9, lwd = 2, lty = 2)
legend("topright", c("NG-Reg", "Lasso", "CAR", "PCR", "MCP","SCAD"), 
       lty = c(1,1,1,1,1,1), col = c(1,2,5,6,7,8), lwd = c(2,2,2,2,2,2), cex = 1)
dev.off() 
