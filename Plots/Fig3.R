load("/Users/shentao/Desktop/plot\ for\ NGFS/fig3-left.Rdata")
load("/Users/shentao/Desktop/plot\ for\ NGFS/fig3-right.Rdata")

pdf("clu_result_pre.pdf", width = 14, height = 6)  
par(mfrow = c(1,2), mar = c(5, 4, 4, 8),                                  # Specify par parameters
    xpd = TRUE)
plot(kseq, TPsina[,1], type = "n", ylim = c(min(TPsina), max(TPsina)), 
     xlab = TeX("\\hat{K}"), ylab = "")
lines(kseq, TPsina[,1], col = 2, lwd = 2, lty = 2)
lines(kseq, TPsina[,2], col = 1, lwd = 2)
legend("topright", inset = c(-0.32, 0), c("Positives","True Positives"), 
       col = c(1,2), lwd = c(2,2), lty = c(1,2), cex = 0.8)



plot(pseq_sina, TPsina2[,1], type = "n", ylim = c(0, 11), 
     xlab = "Positives", ylab = "True Positives")
lines(pseq_sina, TPsina2[,1], col = 2, lwd = 2)
#lines(pseq, TPsina2[,2], col = 2, lwd = 2)
lines(pseq_sina, TPsina2[,3], col = 3, lwd = 2, lty = 2)
lines(pseq_sina, TPsina2[,4], col = 4, lwd = 2, lty = 2)
#lines(pseq, TPsina2[,5], col = 5, lwd = 2)
lines(pseq_sina, TPsina2[,6], col = 1, lwd = 2, lty = 1)
points(13, 8.8,  lwd = 2, pch = 15, col = 'black')
text(13.5,9.2,"HC:(13,8.8)",col=1)
#points(length(fslap_K$selected), TPsina[which(kseq==K),2],  lwd = 2, pch = 15, col = 'red')
legend("topright", inset = c(-0.26, 0), c("Oracle","NGCS", "KS", "SKmeans"), 
       col = c(1,2,3,4), lwd = c(2,2,2,2), lty = c(1,1,2,2), cex = 0.8)

dev.off()