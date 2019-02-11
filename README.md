# l2p program
<pre>
l2p = list to pathway

This an R package for "gene set enrichment".
The interface is a function "l2p()".  The output is a data frame with the following fields ...
 
     1  pval
     2  fdr
     3  ratio                      if positive, genes are OVER REPRESENTED, if negative genes are UNDER REPRESENTED
     4  pathwayhitcount            number of genes hit in pathway
     5  numberofgenesinpathway     number of genes in the pathway
     6  inputnumberofgenes         total count of user genes (user input)
     7  genesinpathwaysuniverse    total number of unique genes in all pathways
     8  pathwayaccessionidentifier canonical accession ( if available, otherwise assigned by us )
     9  source                     KEGG,REACTOME,GO,PANT(=PANTHER),PID=(pathway interaction database)
    10  pathwayname                Name of pathway
    11  genesinpathway             genes from user that hit the pathway (separated by spaces)
    
Example usage is:
    
library(l2p)
genes <- c( "TP53", "PTEN", "APC" )
x = l2p(as.vector(genes))
options(max.print=1000000)
options(width=10000)
print(x)

Installation:
Download package (l2p_0.1-1.tar.gz) : wget https://github.com/CCBR/l2p/raw/master/l2p_0.1-1.tar.gz
then  run R CMD INSTALL , i.e:

R CMD INSTALL l2p_0.1-1.tar.gz

or, from inside R, run this command:

install.packages("https://github.com/CCBR/l2p/raw/master/l2p_0.1-1.tar.gz", repos=NULL) 

Test program
R --vanilla < test.R

You can get the conda package with this command : wget https://github.com/CCBR/l2p/raw/master/r-l2p-0.0_1-r351_0.tar.bz2
</pre>

