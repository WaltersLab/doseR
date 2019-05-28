#' @title countReadsMara Internal function to find chimeric scaffold.
#' @description This function finds chimeric scaffold from underlying data.
#' @usage countReadsMara(scaffolds, G_RESULT2)
#' @param scaffolds Choice of scaffold.
#' @param G_RESULT2 G_RESULT file.
#' @details This function returns data from internal function to find chimeric scaffold.
#' @return Returns countReadsMara data.
#' @examples
#' data(G.RESULT)
#' @author AJ Vaestermark, JR Walters.
#' @references The "Zwyx" package, 2018 (in press).

############################################################



countReadsMara <- function(scaffolds='DPSCF300001', G_RESULT2) {

  # The following reads from the object, not the files...

  SCAFFOLDS <- as.numeric(unlist(strsplit(scaffolds, "DPSCF30"))[2])

  M1 <-  unname(unlist(G_RESULT2[SCAFFOLDS]@unlistData@elementMetadata@listData[1]))
  F1<-  unname(unlist(G_RESULT2[SCAFFOLDS]@unlistData@elementMetadata@listData[2]))
  M2<-  unname(unlist(G_RESULT2[SCAFFOLDS]@unlistData@elementMetadata@listData[3]))
  F2<-  unname(unlist(G_RESULT2[SCAFFOLDS]@unlistData@elementMetadata@listData[4]))
  M3<-  unname(unlist(G_RESULT2[SCAFFOLDS]@unlistData@elementMetadata@listData[5]))
  F3<-  unname(unlist(G_RESULT2[SCAFFOLDS]@unlistData@elementMetadata@listData[6]))

  M <- (M1+M2+M3)
  F <- (F1+F2+F3)

  PLOT <- log2(M/F)
  PLOT[is.infinite(PLOT)] <- -1
  PLOT[is.nan(PLOT)]      <- -1

  invisible(PLOT)

  # The following reads from the object, not the files...

}

############################################################

#' @title countCpt Function to calculate countCpt.
#' @description This function calculates countCpt.
#' @usage countCpt(temp.cpt)
#' @param temp.cpt temp.cpt file.
#' @details This function calculates countCpt.
#' @return Returns countCpt value for samples.
#' @examples
#' data(G.RESULT)
#' @author AJ Vaestermark, JR Walters.
#' @references The "Zwyx" package, 2018 (in press).

countCpt <- function(temp.cpt=temp.cpt) {

  #scaling <- (300-32) / temp.cpt@cpts[length( temp.cpt@cpts)]
  scaling <- (100-4) / temp.cpt@cpts[length( temp.cpt@cpts)]

  temp.cpt@cpts <- append( temp.cpt@cpts, 0,after = 0)

  temp.cpt@param.est$mean <- round(as.numeric(temp.cpt@param.est$mean))

  temp.cpt@param.est$mean[ temp.cpt@param.est$mean > 0 ] <- "Z"
  temp.cpt@param.est$mean[ temp.cpt@param.est$mean == 0 ] <- "a"
  temp.cpt@param.est$mean[ temp.cpt@param.est$mean < 0 ] <- "W"

  cat('   +------------+')

  for(i in 1:(length(temp.cpt@cpts)-1)) {

    for(j in 1:(         ceiling((temp.cpt@cpts[i+1]-temp.cpt@cpts[i])*scaling)    ) ) {
      cat(temp.cpt@param.est$mean[i])
    }

  }

  message()

}

############################################################



#############

#' @title findChimericScaffold Wrapper function to find chimeric scaffold.
#' @description This function finds chimeric scaffold.
#' @usage findChimericScaffold(G__RESULT)
#' @param G__RESULT G_RESULT file.
#' @details This function finds chimeric scaffold.
#' @return Returns findChimericScaffold value for samples.
#' @examples
#' data(G.RESULT)
#' @author AJ Vaestermark, JR Walters.
#' @references The "Zwyx" package, 2018 (in press).

findChimericScaffold <- function(G__RESULT) {

  requireNamespace("doParallel")
    doParallel::registerDoParallel(cores=10)

  yyx<-NULL
  for(i in 1:5397) {
    j<-i
    j<-str_pad(j, 4, pad = "0")
    yyx[i]<- paste0('DPSCF30',j)
  }

temp1 <- foreach(i=1:10) %dopar% { countReadsMara(scaffolds=yyx[i], G__RESULT)   }
txtplot(unlist(temp1[1]), width=100, height=40)
temp.cpt <- cpt.mean( unlist(temp1[1]) , method = "PELT", Q = 5, penalty= "SIC", minseglen = 5)
countCpt(temp.cpt)
}

