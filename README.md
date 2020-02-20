# l2p program
<pre>
l2p = list to pathway

This an R package for "gene set enrichment".  It is optimized for speed.

Installation:
Download package (l2p_0.1-7.tar.gz) : wget https://github.com/CCBR/l2p/raw/master/l2p_0.1-7.tar.gz
then  run R CMD INSTALL , i.e:

R CMD INSTALL l2p_0.1-7.tar.gz

or, from inside R, run this command:

install.packages("https://github.com/CCBR/l2p/raw/master/l2p_0.1-7.tar.gz", repos=NULL) 

You can get the conda package with this command :
 wget "https://github.com/CCBR/l2p/blob/master/r-l2p-0.0_7-r351_0.tar.bz2?raw=true" -O r-l2p-0.0_7-r351_0.tar.bz2

Available functions:
l2p(genelist)            - return data frame with proabilities that arg (list of genes) matches a pathway
l2pgetlongdesc(acc)  - get the full (possibly very long) description for pathway accession identifer string
l2pgetgenes4acc(acc) - get the list all the genes for a pathway, use the accession.
m2h(mousegeneist)    - return list of human genes for input list of mouse gene names

Convenience Functions:
l2pu(list,universe)  - return data frame with proabilities with list of genes and user specified universe
l2pwcats(list,categeories)     - return data frome with categories specified
l2puwcats(list,universe,categories) - same las l2pwcats but also with a universe
l2pver - return l2p version

l2p is now supporting R style argument passing.
Variable parameters are 

    universe   = list of gene names
    categories = see categories below
    custompathways = A list of vectors in GMT (gene matrix transpose style).  Each custom pw vector : pwname, desc, gene1, gene2 ...
    customfile = a gmt file
    universefile = list of genes one per line

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
     4  pwhitcount                 number of genes hit in pathway
     5  pwnohitcount               pathway number of genes in the pathway
     6  inputcount                  total count of user genes (user input)
     7  pwuniverseminuslist        total number of unique genes in all pathways
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



Test program
R --vanilla < test.R



How to do make a custom pathway:
vec1 = c("lall_ad.2","all_ad","AARS","ABCA1","ABCC9","ACTA1","ACTA2","ACTB","ACTC1","ACTG1","ACTN2","ACTN4","ACVR2B","ACVRL1","ADAR","AFG3L2","AFP","AIP","AK1","AKAP9")
vec2 = c("ACMG_2_0.2","ACMG_2_0","BRCA1","BRCA2","TP53","STK11","MLH1","MSH2","MSH6","PMS2")
vec3 = c("berg_ad.2","berg_ad","AARS","ABCC9","ACTA2","ACTB","ACTC1","ACTG1","ACTN2","ACTN4","ACVR2B","ACVRL1","ADAR","AFG3L2","AIP","AK1","AKAP9","AKT2","AMPD1","ANG","ANK2","ANKH","APC","APOA2","APOA5","APOB","APP","ATL1","ATP1A2","ATP2A2","ATP2C1","ATXN1","ATXN10","ATXN2","ATXN3","ATXN7","AXIN2","BAG3","BCO1","BEST1")
mylist <- list(vec1, vec2,vec3)

genes <- c( "TP53", "PTEN", "APC" , "CENPF" , "DLAT", "TP53" , "NOTAGENE" ,"ABCA1","ABCC9","ACTA1", "ADH1A" ,"ATXN3", "BEST1")
x = l2p(as.vector(genes),custompathways=mylist)
print(length(x))
print(x);
</pre>

