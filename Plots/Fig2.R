load("/Users/shentao/Desktop/fig2-left.Rdata")
load("/Users/shentao/Desktop/fig2-right.Rdata")

museq = seq(1.0, 2.0, 0.1)

pdf("clu_result_pre.pdf", width = 15, height = 6)  
par(mfrow = c(1,2), mar = c(5, 4, 4, 8),                                  # Specify par parameters
    xpd = TRUE)
plot(pseqq, error[,1], type = "n", ylim = c(min(error), max(error)), 
     xlab = "ln(p)", ylab = "Clustering Error Rate", cex.lab = 1)
lines(pseqq, error[,2], col = 1, lwd = 2)
#lines(pseqq, error[,3], col = 1, lwd = 2, lty = 2)
lines(pseqq, error[,4], col = 2, lwd = 2)
lines(pseqq, error[,1], col = 3, lwd = 2)
lines(pseqq, error[,5], col = 4, lwd = 2)
#lines(ppseq, error[,6], col = 6, lwd = 2)
lines(pseqq, error[,7], col = "pink", lwd = 2)
#legend(x = max(aseq)-0.15, y = 64, "X only", cex = 0.7, bty = "n", text.col = 3)
#legend(x = max(aseq)-0.17, y = 50, "Y known", cex = 0.7, bty = "n", text.col = 4)
legend("topright", inset = c(-0.25, 0), c("NGCS-Clu", "IF-PCA", "Spec", "SKmeans", "SAS"), 
       lty = c(1,1,1,1,1), col = c(1,2,3,4,"pink"), lwd = c(2,2,2,2,2), cex = 0.8)


plot(museq, mse[,1], type = "n", ylim = c(min(mse[,-4], na.rm = TRUE), 1.8), 
     xlab = TeX("Signal Strength $\\mu$"), ylab = "Mean Squared Error", cex.lab = 1)
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
legend("topright", inset = c(-0.26, 0), c("NGCS-Reg", "Lasso", "CAR", "PCR", "MCP","SCAD"), 
       lty = c(1,1,1,1,1,1), col = c(1,2,5,6,7,8), lwd = c(2,2,2,2,2,2), cex = 0.8)
dev.off()

