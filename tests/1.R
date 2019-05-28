library(Zwyx)

data(G.RESULT)
temp <- countMedian(scaffolds='DPSCF300112', Sample=1, G_RESULT)

RUnit::checkEqualsNumeric(temp, 117)


