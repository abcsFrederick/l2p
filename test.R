
library(l2p)
genes <- c( "TP53", "PTEN", "APC" )
x = l2p(as.vector(genes))
options(max.print=1000000)
options(width=10000)
print(x)
