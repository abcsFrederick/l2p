# l2p program
<pre>
l2p = list to pathway

This an R package for "gene set enrichment".  It is optimized for speed.

Available functions:
l2p(genelist)            - return data frame with proabilities that arg (list of genes) matches a pathway
l2pu(list,universe)  - return data frame with proabilities with list of genes and user specified universe
l2pwcats(list,categeories)     - return data frome with categories specified
l2puwcats(list,universe,categories) - same las l2pwcats but also with a universe
l2pver - return l2p version
l2pgetlongdesc(acc)  - get the full (possibly very long) description for pathway accession identifer string
l2pgetgenes4acc(acc) - get the list all the genes for a pathway, use the accession.
m2h(mousegeneist)    - return list of human genes for input list of mouse gene names

The "msig" functions are obsoleted, the msigdb pathways are now provided with the "categories" parameter in the "wcats" functions.
Example categories parmeter is "KEGG,PID,C4", this is a string. User can mix and match as desired.

Available categories are :
BIOCYC  - organism specific Pathway/ Genome Databases (PGDBs)  - https://biocyc.org/
GO  - initiative to unify representation of gene and gene product attributes -  http://geneontology.org
KEGG - databases dealing with genomes, biological pathways, - https://www.kegg.jp/
PANTH - databases for protein analysis through evolutionary relationships - http://www.pantherdb.org/
PID  - Pathway interaction database: legacy database from Carl Schaefer & buddies at NCI
REACTOME - curated database of biological pathways - https://reactome.org/
WikiPathways - community resource for biological pathways - https://www.wikipathways.org
C1 - MSigDB only, positional gene sets for each human chromosome and cytogenetic band.
C2 - MSigDB only, curated gene sets from online pathway databases, publications in PubMed, and experts.
C3 - MSigDB only, motif gene sets based on conserved cis-regulatory motifs from comparative analysis
C4 - MSigDB only, computational gene sets defined by mining large collections of cancer-oriented microarray data.
C6 - MSigDB only, oncogenic gene sets defined directly from microarray data from cancer gene perturbations.
C7 - MSigDB only, immunologic gene sets  from microarray data from immunologic studies.

MSigDB "ARCHIVED" pathways are not provided.  MSigDB category "C5" is not there. Use "GO" category (from NCBI biosystems),instead.

An example function call is : x=l2pwcats(as.vector(genelist),"GO,WikiPathways,C4,C5,C6")

The output is a data frame with the following fields ...
 
     1  pval
     2  fdr
     3  ratio                      if positive, genes are OVER REPRESENTED, if negative genes are UNDER REPRESENTED
     4  pathwayhitcount            number of genes hit in pathway
     5  numberofgenesinpathway     number of genes in the pathway
     6  inputnumberofgenes         total count of user genes (user input)
     7  genesinpathwaysuniverse    total number of unique genes in all pathways
     8  pathwayaccessionidentifier canonical accession ( if available, otherwise assigned by us )
     9  category                   KEGG,REACTOME,GO,PANTH,PID(=PANTHER),PID=(pathway interaction database)
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
Download package (l2p_0.1-4.tar.gz) : wget https://github.com/CCBR/l2p/raw/master/l2p_0.1-4.tar.gz
then  run R CMD INSTALL , i.e:

R CMD INSTALL l2p_0.1-4.tar.gz

or, from inside R, run this command:

install.packages("https://github.com/CCBR/l2p/raw/master/l2p_0.1-4.tar.gz", repos=NULL) 

Test program
R --vanilla < test.R

You can get the conda package with this command : wget https://github.com/CCBR/l2p/raw/master/r-l2p-0.0_4-r332_0.tar.bz2
</pre>

