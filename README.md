# l2p
l2p = list to pathway

This an R package for "gene set enrichment".
The interface is a function "l2p()".  The output is a data frame with the following fields ...
<pre>
     1  pval
     2  fdr
     3  ratio                      if postive, genes are OVER REPRESENTED, if negative genes are UNDER REPRESENTED
     4  pathwayhitcount            number of genes hit in pathway
     5  numberofgenesin pathway    number of genes in the pathway
     6  inputnumberofgenes          total count of user genes (user input)
     7  genesinpathwaysuniverse    total number of unique genes in all pathways
     8  pathwayaccessionidentifier canonical accession ( if availible, otherwise assigned by us )
     9  source                     KEGG,REACTOME,GO,PANT(=PANTHER),PID=(pathway interaciton database)
    10  pathwayname                Name of pathway
    11  pathwaytype genes_space_separated   HUGO genes from user that hit the pathway
    
Example usage is:
    
library(l2p)
genes <- c( "TP53", "PTEN", "APC" )
x = l2p(as.vector(genes))
options(max.print=1000000)
options(width=10000)
print(x)

Installation:
R CMD INSTALL l2p_0.1-1.tar.gz

Test program
R --vanilla < test.R
</pre>

