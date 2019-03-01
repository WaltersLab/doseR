#' @title dafs
#' @description This function filters the expression using DAFS (see ref).
#' This is a core function invoked by the DAFS wrapper.
#' @param VEC1 Vector 1, not intended for user interaction.
#' @param PLOT Boolean, toggles plotting.
#' @details This function filters the expression, using
#' Data Adaptive Flag method.
#' @return Returns vx[which.min(vv)] to wrapper function.
#' @examples
#' data(hmel.se)
#' f_se <- dafsFilter(se)
#' @author AJ Vaestermark, JR Walters.
#' @references BMC Bioinformatics, 2014, 15:92

dafs <- function(VEC1, PLOT) {

VEC1 <- VEC1[VEC1!=0]

#take log2 of data
log2xx <- log2(VEC1)

#vector to store Kolmogorov Smirnov distance statistics
vv <- rep(0,0)
vx <- rep(0,0)

#select start point
start <- length(log2xx[log2xx==min(log2xx)])/length(log2xx)

#set sequence
s <- seq(round(start,2),0.5,by=0.005)

#loop through cuts of the data to determine targeted K-S statistic
for(q in s) {

#select data greater than a quantile and run Mclust on that data to
#determine theoretical distribution
d <- log2xx[which(log2xx>quantile(log2xx,q,na.rm=TRUE))]
vx <- c(vx,quantile(log2xx,q,na.rm=TRUE))
out <- mclust::Mclust(d,G=1 , verbose=FALSE )
ks <- suppressWarnings( ks.test(d,"pnorm",out$parameter$mean,
sqrt(out$parameter$variance$sigmasq)) )
vv <- c(vv,ks$statistic)

}

if(PLOT) {
plot(density(log2xx),  main="", xlab="Expression",cex.axis=2)

lines(vx,vv,col="red")

legend(x="topright", legend=c("Data", "KS statistic", "Cutoff"),
col = c("black", "red", "red"),
text.col = c("black", "red", "red"),
lty = c(1, 1, 2), pch = c(NA,NA),
merge = TRUE, bg = "gray90")
abline(v = vx[which.min(vv)], col = "red", lty = 2)
}

return( vx[which.min(vv)] )
} # dafsFilter
