#' @title countReadsMara3 Internal function to calculate ratiosCorrelation.
#' @description This function generates ratiosCorrelation underlying data.
#' @usage countReadsMara3(scaffolds, G_RESULT2)
#' @param scaffolds Choice of scaffold.
#' @param G_RESULT2 G_RESULT file.
#' @details This function returns data from internal function to calculate ratiosCorrelation.
#' @return Returns ratiosCorrelation data.
#' @examples
#' data(G.RESULT)
#' countReadsMara3(scaffolds='DPSCF300002', G_RESULT2=G_RESULT)
#' @author AJ Vaestermark, JR Walters.
#' @references The "Zwyx" package, 2018 (in press).

countReadsMara3 <- function(scaffolds='DPSCF300001', G_RESULT2) {

  requireNamespace("stringr")
  requireNamespace("Rsamtools")


  SCAFFOLDS <- as.numeric(unlist(strsplit(scaffolds, "DPSCF30"))[2])

  M1 <-  unname(unlist(G_RESULT2[SCAFFOLDS]@unlistData@elementMetadata@listData[1]))
  F1<-  unname(unlist(G_RESULT2[SCAFFOLDS]@unlistData@elementMetadata@listData[2]))
  M2<-  unname(unlist(G_RESULT2[SCAFFOLDS]@unlistData@elementMetadata@listData[3]))
  F2<-  unname(unlist(G_RESULT2[SCAFFOLDS]@unlistData@elementMetadata@listData[4]))
  M3<-  unname(unlist(G_RESULT2[SCAFFOLDS]@unlistData@elementMetadata@listData[5]))
  F3<-  unname(unlist(G_RESULT2[SCAFFOLDS]@unlistData@elementMetadata@listData[6]))


  #Len <- suppressMessages( getReadCountsFromBAM(x[1], WL=1, refSeqName=scaffolds) )



  M_MW <- (M1+M2+M3) # That's mean windows?
  F_MW <- (F1+F2+F3)

  M_SR <- sum(M1+M2+M3)
  F_SR <- sum(F1+F2+F3)

  PLOT_MW <- mean(M_MW/F_MW, na.rm=T)  ## LOOKS RIGHT
  PLOT_SR <- (M_SR/F_SR)  ## LOOKS RIGHT

  if ( is.finite(PLOT_MW) && is.finite(PLOT_SR) ) { return(c(PLOT_MW,   PLOT_SR )) }

}

#############

#' @title ratiosCorrelation Function to plot ratiosCorrelation.
#' @description This function plots ratiosCorrelation.
#' @usage ratiosCorrelation(GRESULT3)
#' @param GRESULT3 G_RESULT file.
#' @details This function plots ratiosCorrelation.
#' @return Returns ratiosCorrelation value plot for samples.
#' @examples
#' data(G.RESULT)
#' ratiosCorrelation(G_RESULT)
#' @author AJ Vaestermark, JR Walters.
#' @references The "Zwyx" package, 2018 (in press).

ratiosCorrelation <- function(GRESULT3) {

  yyx<-NULL
  for(i in 1:5397) {
    j<-i
    j<-str_pad(j, 4, pad = "0")
    yyx[i]<- paste0('DPSCF30',j)
  }

  x_values<-NULL
  y_values<-NULL
  Sample <- sample (1:500, 300)
  pb <- txtProgressBar(max=300)
  for(j in 1:300){
    setTxtProgressBar(pb, j)
    x_values<-c(x_values,countReadsMara3 (scaffolds=yyx[Sample[j]], GRESULT3)[1])
    y_values<-c(y_values,countReadsMara3 (scaffolds=yyx[Sample[j]], GRESULT3)[2])
  }

  close(pb)

  txtplot(x_values-1, y_values-1, width=50, height=20, xlim=c(-5,2) , ylim=c(-5,2))



}
