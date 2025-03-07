library(Cardinal)
library(CardinalIO)
set.seed(2025, kind="L'Ecuyer-CMRG")
setCardinalParallel(6)

parseImzML('', ibd=TRUE)
