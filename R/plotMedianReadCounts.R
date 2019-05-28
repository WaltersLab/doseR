#' @title countMedian Function to calculate median.
#' @description This function generates "median".
#' @usage countMedian(scaffolds, Sample, G_RESULT2)
#' @param scaffolds Choice of scaffold.
#' @param Sample Choice of sample.
#' @param G_RESULT2 G_RESULT file.
#' @details This function returns median expression for selected scaffold and sample.
#' @return Returns median value.
#' @examples
#' data(G.RESULT)
#' countMedian(scaffolds='DPSCF300112', Sample=1, G_RESULT)
#' @author AJ Vaestermark, JR Walters.
#' @references The "Zwyx" package, 2018 (in press).

# The following reads from the object, not the files...

countMedian <- function(scaffolds='DPSCF300001', Sample=1, G_RESULT2) {

 # data(G.RESULT)

  x <- c('Danaus_poolC_Plex_BOS_HI023_M_CTTGGA_L003realigned.bam', 'MK_Plex_Mex_1742_F_GATCAG_L002realigned.bam', 'MK_Plex_NJ_116_M_GAGTGG_L008.realigned.bam',
         'MK_Plex_StMFL_122_F_ATCACG_L008realigned.bam', 'MK_Plex_StMFL_109_M_ACTTGA_L007realigned.bam', 'MK_Plex_TX_T11_F_GGCTAC_L0013realigned.bam')

  i<-Sample

  SCAFFOLDS <- as.numeric(unlist(strsplit(scaffolds, "DPSCF30"))[2])
  sample2 <- x[Sample]

  #G_RESULT[1]@unlistData@elementMetadata@listData$Danaus_poolC_Plex_BOS_HI023_M_CTTGGA_L003realigned.bam
  G_RESULT2[SCAFFOLDS]@unlistData@elementMetadata@listData[Sample]

  #M3 <- unname( unlist( ( suppressMessages(    getReadCountsFromBAM(x[i], WL=500, refSeqName=scaffolds)))@elementMetadata@listData[1]));
  M3 <-  unname(unlist(G_RESULT2[SCAFFOLDS]@unlistData@elementMetadata@listData[Sample]))

  return( median(M3) )

}

############




############




##################

# highest allowed
#DPSCF305397



##################

# The following reads from the object, not the files...

#' @title countMedian Function to plot median read counts.
#' @description This function generates "median" read count plot.
#' @usage plotMedianReadCounts(G_RESULT3)
#' @param G_RESULT3 G_RESULT file.
#' @details This function plots median expression for selected sample.
#' @return Returns median value plot for samples.
#' @examples
#' data(G.RESULT)
#' plotMedianReadCounts(G_RESULT)
#' @author AJ Vaestermark, JR Walters.
#' @references The "Zwyx" package, 2018 (in press).


plotMedianReadCounts <- function(G_RESULT3) {

  requireNamespace("txtplot")
  requireNamespace("stringr")

  yyx<-NULL
  for(i in 1:5397) {
    j<-i
    j<-str_pad(j, 4, pad = "0")
    yyx[i]<- paste0('DPSCF30',j)
  }

  pb <- txtProgressBar(max=6)

  SAMPLE <-  sample(1:5397, 30) ;
xHI023_M <- NULL
for(i in 1:30) xHI023_M <- c( xHI023_M , countMedian(yyx[SAMPLE[i]], Sample=1, G_RESULT3) )
setTxtProgressBar(pb, 1)
x1742_F <- NULL
for(i in 1:30) x1742_F <- c( x1742_F , countMedian(yyx[SAMPLE[i]], Sample=2, G_RESULT3) )
setTxtProgressBar(pb, 2)
x116_M <- NULL
for(i in 1:30) x116_M <- c( x116_M , countMedian(yyx[SAMPLE[i]], Sample=3, G_RESULT3) )
setTxtProgressBar(pb, 3)
x122_F <- NULL
for(i in 1:30) x122_F <- c( x122_F , countMedian(yyx[SAMPLE[i]], Sample=4, G_RESULT3) )
setTxtProgressBar(pb, 4)
x109_M <- NULL
for(i in 1:30) x109_M <- c( x109_M , countMedian(yyx[SAMPLE[i]], Sample=5, G_RESULT3) )
setTxtProgressBar(pb, 5)
xT11_F <- NULL
for(i in 1:30) xT11_F <- c( xT11_F , countMedian(yyx[SAMPLE[i]], Sample=6, G_RESULT3) )
setTxtProgressBar(pb, 6)

close(pb)

txtboxplot(x109_M, x116_M ,xT11_F,x122_F,x1742_F,xHI023_M, width=50)

}

###################################



##################

#https://github.com/WaltersLab/Dplex_Changepoint/blob/master/Results/Dplex_v3/Monarch__MedianReadCounts.pdf


