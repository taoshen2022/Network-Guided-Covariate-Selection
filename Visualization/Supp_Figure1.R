load("./Visualization/Supp_HCtab.Rdata")
library(latex2exp)

re_aseq <- c(0.9, 0.7, 0.5, 0.3, 0.1)

pdf("ReF9-supp.pdf", width = 10, height = 6)  
par(mfrow = c(1,2),                                   
    xpd = TRUE)

plot(re_aseq, HCtab[,1]/50, type = "n", ylim = c(min(HCtab/50)-0.1, 1), 
     xlab = TeX("1-Between Community Intensity a"), ylab = TeX("Proportion of Violation"), cex.lab = 1.0)
lines(re_aseq, HCtab[,3]/50, col = 2, lwd = 2)

plot(re_aseq, FDR[,1], type = "n", ylim = c(min(FDR), 1), 
     xlab = TeX("1-Between Community Intensity a"), ylab = TeX("False Discovery Rate $|\\hat{S}\\cap S^c|/|\\hat{S}|$"), cex.lab = 1.0)
lines(re_aseq, FDR[,3], col = 2, lwd = 2)

dev.off() 
