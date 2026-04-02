library(latex2exp)

pdf("ReF9.pdf", width = 11, height = 8)
##param setup
par(
  mfrow = c(3, 3),
  mar   = c(3.2, 3.2, 2.0, 2.0),  
  oma   = c(0.6, 0.6, 0.4, 0.4),
  mgp   = c(1.5, 0.5, 0),
  tcl   = -0.25,
  xaxs  = "i", yaxs = "i",
  xpd   = NA
)

museq = seq(0.1, 0.5, by = 0.05);
load("./Visualization/FDR_DCSBM_(a).Rdata")
plot(museq, FDR[,1], main = "DCSBM, Gaussian", type = "n", ylim = c(0, 1), 
     xlab = TeX("Signal Strength $\\mu$"), ylab = TeX("False Discovery Rate $|\\hat{S}\\cap S^c|/|\\hat{S}|$"), cex.lab = 1.0)
lines(museq, FDR[,1], col = 1, lwd = 2)
lines(museq, FDR[,4], col = 3, lwd = 2, lty = 2)
lines(museq, FDR[,5], col = 4, lwd = 2, lty = 2)
lines(museq, FDR[,7], col = 5, lwd = 2)
lines(museq, FDR[,8], col = 6, lwd = 2)
lines(museq, FDR[,3], col = 2, lwd = 3)


load("./Visualization/FDR_DCMM_(a).Rdata")
plot(museq, FDR[,1], main = "DCMM, Gaussian", type = "n", ylim = c(0, 1), 
     xlab = TeX("Signal Strength $\\mu$"), ylab = TeX("False Discovery Rate $|\\hat{S}\\cap S^c|/|\\hat{S}|$"), cex.lab = 1.0)
lines(museq, FDR[,1], col = 1, lwd = 2)
lines(museq, FDR[,4], col = 3, lwd = 2, lty = 2)
lines(museq, FDR[,5], col = 4, lwd = 2, lty = 2)
lines(museq, FDR[,7], col = 5, lwd = 2)
lines(museq, FDR[,8], col = 6, lwd = 2)
lines(museq, FDR[,3], col = 2, lwd = 3)


load("./Visualization/FDR_RDPG_(a).Rdata")
museq = seq(0.05, 0.30, by = 0.05);
plot(museq, FDR[,1], main = "RDPG, Gaussian", type = "n", ylim = c(0, 1), 
     xlab = TeX("Signal Strength $\\mu$"), ylab = TeX("False Discovery Rate $|\\hat{S}\\cap S^c|/|\\hat{S}|$"), cex.lab = 1.0)
lines(museq, FDR[,1], col = 1, lwd = 2)
lines(museq, FDR[,4], col = 3, lwd = 2, lty = 2)
lines(museq, FDR[,5], col = 4, lwd = 2, lty = 2)
lines(museq, FDR[,7], col = 5, lwd = 2)
lines(museq, FDR[,8], col = 6, lwd = 2)
lines(museq, FDR[,3], col = 2, lwd = 3)
legend("topright", c("NGCS-HCT-A", "NGCS-HCT-L", "NGCS-50", "NGCS-HW", "Chi", "SKmeans"), 
       col = c(2,3,1,4,5,6), lwd = c(3,2,2,2,2,2), lty = c(1,2,1,2,1,1), cex = 1.0)


load("./Visualization/FDR_DCSBM_(b).Rdata")
museq = seq(0.1, 0.5, by = 0.05);
plot(museq, FDR[,1], main = "DCSBM, Chi-square", type = "n", ylim = c(0, 1), 
     xlab = TeX("Signal Strength $\\mu$"), ylab = TeX("False Discovery Rate $|\\hat{S}\\cap S^c|/|\\hat{S}|$"), cex.lab = 1.0)
lines(museq, FDR[,1], col = 1, lwd = 2)
lines(museq, FDR[,4], col = 3, lwd = 2, lty = 2)
lines(museq, FDR[,5], col = 4, lwd = 2, lty = 2)
lines(museq, FDR[,7], col = 5, lwd = 2)
lines(museq, FDR[,8], col = 6, lwd = 2)
lines(museq, FDR[,3], col = 2, lwd = 3)

       
load("./Visualization/FDR_DCMM_(b).Rdata")
plot(museq, FDR[,1], main = "DCMM, Chi-square", type = "n", ylim = c(0, 1), 
     xlab = TeX("Signal Strength $\\mu$"), ylab = TeX("False Discovery Rate $|\\hat{S}\\cap S^c|/|\\hat{S}|$"), cex.lab = 1.0)
lines(museq, FDR[,1], col = 1, lwd = 2)
lines(museq, FDR[,4], col = 3, lwd = 2, lty = 2)
lines(museq, FDR[,5], col = 4, lwd = 2, lty = 2)
lines(museq, FDR[,7], col = 5, lwd = 2)
lines(museq, FDR[,8], col = 6, lwd = 2)
lines(museq, FDR[,3], col = 2, lwd = 3)


load("./Visualization/FDR_RDPG_(b).Rdata")
museq = seq(0.05, 0.30, by = 0.05);
plot(museq, FDR[,1], main = "RDPG, Chi-square", type = "n", ylim = c(0, 1), 
     xlab = TeX("Signal Strength $\\mu$"), ylab = TeX("False Discovery Rate $|\\hat{S}\\cap S^c|/|\\hat{S}|$"), cex.lab = 1.0)
lines(museq, FDR[,1], col = 1, lwd = 2)
lines(museq, FDR[,4], col = 3, lwd = 2, lty = 2)
lines(museq, FDR[,5], col = 4, lwd = 2, lty = 2)
lines(museq, FDR[,7], col = 5, lwd = 2)
lines(museq, FDR[,8], col = 6, lwd = 2)
lines(museq, FDR[,3], col = 2, lwd = 3)


museq = seq(0.1, 0.5, by = 0.05);
load("./Visualization/FDR_DCSBM_(c).Rdata")
plot(museq, FDR[,1], main = "DCSBM, Sub-Gaussian", type = "n", ylim = c(0, 1), 
     xlab = TeX("Signal Strength $\\mu$"), ylab = TeX("False Discovery Rate $|\\hat{S}\\cap S^c|/|\\hat{S}|$"), cex.lab = 1.0)
lines(museq, FDR[,1], col = 1, lwd = 2)
lines(museq, FDR[,4], col = 3, lwd = 2, lty = 2)
lines(museq, FDR[,5], col = 4, lwd = 2, lty = 2)
lines(museq, FDR[,7], col = 5, lwd = 2)
lines(museq, FDR[,8], col = 6, lwd = 2)
lines(museq, FDR[,3], col = 2, lwd = 3)


load("./Visualization/FDR_DCMM_(c).Rdata")
plot(museq, FDR[,1], main = "DCMM, Sub-Gaussian", type = "n", ylim = c(0, 1), 
     xlab = TeX("Signal Strength $\\mu$"), ylab = TeX("False Discovery Rate $|\\hat{S}\\cap S^c|/|\\hat{S}|$"), cex.lab = 1.0)
lines(museq, FDR[,1], col = 1, lwd = 2)
lines(museq, FDR[,4], col = 3, lwd = 2, lty = 2)
lines(museq, FDR[,5], col = 4, lwd = 2, lty = 2)
lines(museq, FDR[,7], col = 5, lwd = 2)
lines(museq, FDR[,8], col = 6, lwd = 2)
lines(museq, FDR[,3], col = 2, lwd = 3)


load("./Visualization/FDR_RDPG_(c).Rdata")
museq = seq(0.05, 0.30, by = 0.05);
plot(museq, FDR[,1], main = "RDPG, Sub-Gaussian", type = "n", ylim = c(0, 1), 
     xlab = TeX("Signal Strength $\\mu$"), ylab = TeX("False Discovery Rate $|\\hat{S}\\cap S^c|/|\\hat{S}|$"), cex.lab = 1.0)
lines(museq, FDR[,1], col = 1, lwd = 2)
lines(museq, FDR[,4], col = 3, lwd = 2, lty = 2)
lines(museq, FDR[,5], col = 4, lwd = 2, lty = 2)
lines(museq, FDR[,7], col = 5, lwd = 2)
lines(museq, FDR[,8], col = 6, lwd = 2)
lines(museq, FDR[,3], col = 2, lwd = 3)

dev.off() 