
library(doseR)
#library(SummarizedExperiment)

data(hmel.se)
f_se <- quantFilter(se, lo.bound = 0.4, hi.bound = 0.5)
dm <- se.DM(f_se)
temp <- glSeq(dm, "-1 + replicate")

RUnit::checkEqualsNumeric(round(temp@optinfo$val[4], 1), round(6.008, 1))
RUnit::checkEqualsNumeric(round(temp@optinfo$val[5], 1), round(5.840, 1))
