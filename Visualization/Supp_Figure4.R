library(latex2exp)
load("./Visualization/Sfig4-left-rdpg.Rdata")
load("./Visualization/Sfig4-right-rdpg.Rdata")

pdf("rdpg.pdf", width = 15, height = 6)  
par(mfrow = c(1,2), mar = c(4, 5, 3, 9),                                  # Specify par parameters
    xpd = TRUE)

plot(mseq, tau_matrix33[,1], type = "n", ylim = c(min(tau_matrix33), max(tau_matrix33)), 
     xlab = TeX(paste("Signal Strength", "$\\mu$")), ylab = TeX("Aggregated Signal Strength $\\tau$"),cex = 1.5)
lines(mseq, tau_matrix33[,5], col = 1, lwd = 2, lty = 1)
lines(mseq, tau_matrix33[,1], col = 2, lwd = 2)
lines(mseq, tau_matrix33[,2], col = 2, lwd = 2, lty = 2)
lines(mseq, tau_matrix33[,3], col = 2, lwd = 1.5, lty = 3)
lines(mseq, tau_matrix33[,4], col = 2, lwd = 1.5, lty = 4)
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

dev.off() 