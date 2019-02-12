library('RUnit')

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

checkEqualsNumeric(temp@optinfo$val[4], 6.008)
checkEqualsNumeric(temp@optinfo$val[5], 5.840)
