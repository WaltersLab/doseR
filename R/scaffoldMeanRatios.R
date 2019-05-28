#' @title countReadsMara2 Internal function to calculate scaffoldMeanRatios.
#' @description This function generates scaffoldMeanRatios underlying data.
#' @usage countReadsMara2(scaffolds, G_RESULT2)
#' @param scaffolds Choice of scaffold.
#' @param G_RESULT2 G_RESULT file.
#' @details This function returns data from internal function to calculate scaffoldMeanRatios.
#' @return Returns scaffoldMeanRatios data.
#' @author AJ Vaestermark, JR Walters.
#' @references The "Zwyx" package, 2018 (in press).

countReadsMara2 <- function(scaffolds='DPSCF300001', G_RESULT2 ) {

  requireNamespace("stringr")
  requireNamespace("Rsamtools")

  x <- c('extdata/Danaus_poolC_Plex_BOS_HI023_M_CTTGGA_L003realigned.bam', 'MK_Plex_Mex_1742_F_GATCAG_L002realigned.bam', 'MK_Plex_NJ_116_M_GAGTGG_L008.realigned.bam','MK_Plex_StMFL_122_F_ATCACG_L008realigned.bam', 'MK_Plex_StMFL_109_M_ACTTGA_L007realigned.bam', 'MK_Plex_TX_T11_F_GGCTAC_L0013realigned.bam')
  header <- scanBamHeader(files=x[1])
  scaffLengths <- header[[1]][[1]]

  Len <- unname(scaffLengths[as.numeric(unlist(str_split(scaffolds, "DPSCF30"))[2])])

  SCAFFOLDS <- as.numeric(unlist(strsplit(scaffolds, "DPSCF30"))[2])

  M1 <-  unname(unlist(G_RESULT2[SCAFFOLDS]@unlistData@elementMetadata@listData[1]))
  F1<-  unname(unlist(G_RESULT2[SCAFFOLDS]@unlistData@elementMetadata@listData[2]))
  M2<-  unname(unlist(G_RESULT2[SCAFFOLDS]@unlistData@elementMetadata@listData[3]))
  F2<-  unname(unlist(G_RESULT2[SCAFFOLDS]@unlistData@elementMetadata@listData[4]))
  M3<-  unname(unlist(G_RESULT2[SCAFFOLDS]@unlistData@elementMetadata@listData[5]))
  F3<-  unname(unlist(G_RESULT2[SCAFFOLDS]@unlistData@elementMetadata@listData[6]))

  #M <- mean(M1+M2+M3) # That's mean windows?
  #F <- mean(F1+F2+F3)

  M <- sum(M1+M2+M3)
  F <- sum(F1+F2+F3)

  PLOT <- log2(F/M)  ## LOOKS RIGHT

  #return(c(PLOT,   Len@seqnames@lengths ))

  if ( is.finite(PLOT) && is.finite(Len) ) { return(c(PLOT,   Len )) }

  #return(c(PLOT,   Len ))

}


#################

#' @title scaffoldMeanRatios Function to plot scaffoldMeanRatios.
#' @description This function plots scaffoldMeanRatios.
#' @usage scaffoldMeanRatios(GRESULT3)
#' @param GRESULT3 G_RESULT file.
#' @details This function plots scaffoldMeanRatios.
#' @return Returns scaffoldMeanRatios value plot for samples.
#' @author AJ Vaestermark, JR Walters.
#' @references The "Zwyx" package, 2018 (in press).

scaffoldMeanRatios <- function(GRESULT3) {

  requireNamespace("doParallel")
  requireNamespace("parallel")
  requireNamespace("foreach")
  doParallel::registerDoParallel(cores=10)

  Xval <- NULL
  Yval <- NULL

  yyx<-NULL
  for(i in 1:5397) {
    j<-i
    j<-str_pad(j, 4, pad = "0")
    yyx[i]<- paste0('DPSCF30',j)
  }

  pb <- txtProgressBar(max=10)


  #globalVariables('i')
  #globalVariables('%dopar%')

  for(j in 1:10) {
    setTxtProgressBar(pb, j)
    SAMPLE <-  sample(1:5397, 5) ;




    xHI023_M <- (  foreach(i=1:5) %dopar% countReadsMara2(scaffolds=yyx[SAMPLE[i]]  , GRESULT3) ) ;
    for(i in 1:5)  Xval <- c(Xval, unlist(xHI023_M)[i*2])
    for(i in 1:5)  Yval <- c(Yval, unlist(xHI023_M)[i*2-1])
  }

  close(pb)


  #txtplot(log10(Xval), Yval, width=300, height=60) # for summed reads
  txtplot(log10(Xval), Yval, width=50, height=20)  # for mean windows



}

