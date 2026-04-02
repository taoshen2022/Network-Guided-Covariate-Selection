load("./Visualization/Sfig5-left.Rdata")
load("./Visualization/Sfig5-right.Rdata")

museq <- seq(1.0, 2.0, 0.1)
pseqq <- seq(6.0, 12.0, 0.5)
pdf("clu_result_pre.pdf", width = 14, height = 6)  
par(mfrow = c(1,2), mar = c(5, 4, 4, 8),                                  # Specify par parameters
    xpd = TRUE)
plot(pseqq, error[,1], type = "n", ylim = c(min(error), 0.85), 
     xlab = "ln(p)", ylab = "Clustering Error Rate", cex.lab = 1)
lines(pseqq, error[,2], col = 1, lwd = 2)
lines(pseqq, error[,4], col = 2, lwd = 2)
lines(pseqq, error[,1], col = 3, lwd = 2)
lines(pseqq, error[,5], col = 4, lwd = 2)
lines(pseqq, error[,7], col = "pink", lwd = 2)
legend("topright", c("NG-clu", "IF-PCA", "Spec", "SKmeans", "SAS"), 
       lty = c(1,1,1,1,1), col = c(1,2,3,4,"pink"), lwd = c(2,2,2,2,2), cex = 0.8)


plot(museq, mse[,1], type = "n", ylim = c(min(mse[,-4], na.rm = TRUE), 1.8), 
     xlab = TeX("Signal Strength $\\mu$"), ylab = "Mean Squared Error", cex.lab = 1)
lines(museq, mse[,2], col = 1, lwd = 2)
lines(museq, mse[,1], col = 2, lwd = 2)
lines(museq, mse[,5], col = 5, lwd = 2)
lines(museq, mse[,6], col = 6, lwd = 2)
lines(museq, mse[,8], col = 7, lwd = 2)
lines(museq, mse[,9], col = 8, lwd = 2)
legend("topright", c("NG-reg", "Lasso", "CAR", "PCR", "MCP","SCAD"), 
       lty = c(1,1,1,1,1,1), col = c(1,2,5,6,7,8), lwd = c(2,2,2,2,2,2), cex = 0.8)
dev.off()

