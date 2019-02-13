
library(doseR)

data(hmel.data.doser)
reps <- c("Male", "Male", "Male", "Female", "Female", "Female")
annotxn <- data.frame("Chromosome" = factor(hmel.dat$chromosome,
levels = 1:21))
hm.tr<-hmel.dat$trxLength
hm<-new("countDat",data=hmel.dat$readcounts,seglens=hm.tr,
annotation=annotxn)
replicates(hm) <- reps
libsizes(hm) <- getLibsizes2(hm, estimationType = "total")
rpkm(hm) <- make_RPKM(hm)
f_hm <- quantFilter(hm, lo.bound = 0.4, hi.bound = 0.5)
dm <- cD.DM(f_hm)
temp <- glSeq(dm, "-1 + replicate")

RUnit::checkEqualsNumeric(round(temp@optinfo$val[4], 1), round(6.008, 1))
RUnit::checkEqualsNumeric(round(temp@optinfo$val[5], 1), round(5.840, 1))
