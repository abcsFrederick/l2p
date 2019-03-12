#install.packages("https://github.com/CCBR/l2p/raw/master/l2p_0.1-3.tar.gz", repos=NULL)
library(l2p)

# load a list of unique gene names from a file - you must create this 
universe <- scan("u.zhsu", what="", sep="\n")

# load a list of unique gene names from a file - you must create this 
upgenes <- scan("/up.zhsorted", what="", sep="\n")

# call l2p with a universe (u) ...
x2 = l2pu(as.vector(upgenes),universe)

options(max.print=1000000)
options(width=10000)
print(x2)

