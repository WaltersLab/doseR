## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval = F-----------------------------------------------------------
#  
#  setwd("~/Documents/doseR_Project/doseR_vignette")
#  path.to.doser <- "./doseR_1.4.0.tar.gz"
#  install.packages(pkgs = path.to.doser, repos = NULL, type="source", lib = getwd() )

## ---- eval = F-----------------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("doseR")

## ------------------------------------------------------------------------
library(mclust)
library(doseR)
library(edgeR)

## ------------------------------------------------------------------------
#load("hmel.data.doser.Rdata") # loads hmel.dat
data(hmel.data.doser)  # loads hmel.dat list
lapply(hmel.dat, head)

## ------------------------------------------------------------------------
reps <- c("Male", "Male", "Male", "Female", "Female", "Female")

## ------------------------------------------------------------------------
annotxn <- data.frame("Chromosome" = factor(hmel.dat$chromosome, levels = 1:21))
annotxn$ZA <- factor(ifelse(hmel.dat$chromosome == 21, "Z", "A"), levels = c("A", "Z"))

## ------------------------------------------------------------------------
hmel.cd <- new("countDat", data = hmel.dat$readcounts, replicates = reps, seglens = hmel.dat$trxLength, annotation = annotxn)

hmel.cd

## ------------------------------------------------------------------------
libsizes(hmel.cd) <- getLibsizes2(hmel.cd, estimationType = "total")
libsizes(hmel.cd)

## ------------------------------------------------------------------------
hmel.cd@RPKM <- make_RPKM(hmel.cd)
hmel.cd

## ---- fig.align="center", fig.width=6------------------------------------
plotMA.CD(cD = hmel.cd, samplesA ="Male", samplesB = "Female", cex = .2 , pch = 19, col = rgb(0,0,0, .2), xlab = "Log2(Average RPKM)", ylab = "Log2(Male:Female)")

## ------------------------------------------------------------------------
hmel.filt <- simpleFilter(cD = hmel.cd, mean_cutoff = 0.01, counts = F)


## ---- fig.align="center", fig.width=5------------------------------------
plotExpr(cD = hmel.filt, groupings = "ZA", clusterby_grouping = F, col=c("grey80","red","grey80","red"), notch=T, ylab = "Log2(RPKM)")

## ------------------------------------------------------------------------
hmel.male <- hmel.filt[, hmel.filt@replicates == "Male"]
male_ZvA <- generateStats(cD = hmel.male , groupings = "ZA", LOG2 = F)
male_ZvA$summary  # distributional summary statistics
male_ZvA$kruskal  # htest class output from kruskal.test()
lapply(male_ZvA$data, head) # a record of values used for statistics.

## ---- fig.align="center"-------------------------------------------------
plotExpr(cD = hmel.filt, groupings = "Chromosome", col=c(rep("grey80", 20), "red"), notch=T, ylab = "Log2(RPKM)", las = 2, treatment = "Male",  clusterby_grouping = T )

## ---- fig.align="center", fig.width=4------------------------------------
plotRatioBoxes(cD = hmel.filt, groupings = "ZA", treatment1 = "Male", treatment2 = "Female", outline = F, col = c("grey80", "red"), ylab = "Log(Male:Female)" )

## ---- fig.align="center", fig.width = 5----------------------------------
plotRatioDensity(cD = hmel.filt, groupings = "ZA", treatment1 = "Male", treatment2 = "Female", type = "l", xlab = "Log(Male:Female)", ylab = "Density")

## ---- fig.align="center", fig.width=10-----------------------------------
par(mfrow = c(1,2))
plotRatioBoxes(cD = hmel.filt, groupings = "Chromosome", treatment1 = "Male", treatment2 = "Female", outline = F, col=c(rep("grey80", 20), "red"), ylab = "Log(Male:Female)", xlab = "Chromosome" )

plotRatioDensity(cD = hmel.filt, groupings = "Chromosome", treatment1 = "Male", treatment2 = "Female", type = "l", xlab = "Log(Male:Female)", ylab = "Density", col=c(rep("grey80", 20), "red"), lty = 1)

## ------------------------------------------------------------------------
za.ratios.test <- test_diffs(cD = hmel.filt, groupings = "ZA", treatment1 = "Male", treatment2 = "Female", LOG2 = F )
za.ratios.test$summary  # summary statistics for each grouping
za.ratios.test$kruskal  # htest class output from kruskal.test()
lapply(za.ratios.test$data, head) # values used for summaries and tests

