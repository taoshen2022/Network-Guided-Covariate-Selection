library(latex2exp)
load("./Visualization/Sfig3-left-dcsbmm.Rdata")
load("./Visualization/Sfig3-right-dcsbmm.Rdata")

pdf("dcmm.pdf", width = 15, height = 6)  
par(mfrow = c(1,2), mar = c(4, 5, 3, 9),                                  # Specify par parameters
    xpd = TRUE)

plot(tseq, tau_matrix22[,1], type = "n", ylim = c(min(tau_matrix22), max(tau_matrix22)), 
     xlab = "Mixture Parameter h", ylab = TeX("Aggregated Signal Strength $\\tau$"), cex=1.5)
lines(tseq, tau_matrix22[,5], col = 1, lwd = 2, lty = 1)
lines(tseq, tau_matrix22[,1], col = 2, lwd = 2)
lines(tseq, tau_matrix22[,2], col = 2, lwd = 2, lty = 2)
lines(tseq, tau_matrix22[,3], col = 2, lwd = 1.5, lty = 3)
lines(tseq, tau_matrix22[,4], col = 2, lwd = 1.5, lty = 4)
legend("bottomleft", c("Oracle", "NGCS(A,K)", "NGCS(L,K)","NGCS(A,2K)", "NGCS(L,2K)"), 
       lty = c(1,1,2,3,4), col = c(1,2,2,2,2), lwd = c(2,2,2,1.5,1.5), cex = 1)



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
abline(v = 53, col = 2, lwd = 1, lty = 2, xpd = FALSE)
points(53, 0.0754717,  lwd = 1, pch = 15, col = 1)
text(140, 0.08, "HC:(53, 0.07)", col = 1)  
legend("bottomright", c("Oracle", "NGCS(A,K)", "NGCS(L,K)", "Chi","FSCA","SKmeans"), 
       col = c(1,2,2,3,4,5), lwd = c(2,2,2,2,2,2), lty = c(1,1,2,1,1,1), cex = 1)


dev.off() 