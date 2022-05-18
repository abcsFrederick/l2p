# l2p

[![DOI](https://zenodo.org/badge/164483194.svg)](https://zenodo.org/badge/latestdoi/164483194)
[![GitHub releases](https://img.shields.io/github/release/CCBR/l2p)](https://github.com/CCBR/l2p/releases) [![GitHub issues](https://img.shields.io/github/issues/CCBR/l2p)](https://github.com/CCBR/l2p/issues) [![GitHub license](https://img.shields.io/github/license/CCBR/l2p)](https://github.com/CCBR/l2p/blob/master/LICENSE)

The beta web version is here: https://ccbr.github.io/l2p/ . ( You can copy/paste a whole paper in to see an instant pathway analysis. )

List-to-pathway, or `l2p`, is an R package for gene set enrichment analysis that is _optimized for speed!_ 

`l2p` can be used to determine whether a biological process or function is over-represented in a user-defined gene list. This can be a list of differential expressed genes, or a list of annotated differential bound regions using a tool like [uropa](https://www.nature.com/articles/s41598-017-02464-y) or [homer](http://homer.ucsd.edu/homer/ngs/annotation.html). 

l2psupp is the "l2p supplemental" package which contains routines for converting gene symbols.  l2psupp 


## Installation

The latest package of `l2p` can be downloaded directly from Github. Here we describe each method in more detail.

**Option 1:** Download and install latest R package, `l2p_0.0-13.tar.gz`, from command-line:
```bash
# Get l2p from Github
wget https://github.com/CCBR/l2p/raw/master/l2p_0.0-13.tar.gz
# Install as a site package 
R CMD INSTALL l2p_0.0-13.tar.gz
# install l2psupp ( "l2p supplemental")
wget https://github.com/CCBR/l2p/raw/master/l2psupp_0.0-13.tar.gz
# Install as a site package 
R CMD INSTALL l2psupp_0.0-13.tar.gz
```
 
**Option 2:** Install `l2p` within an R console or RStudio session:
```R
# Install from R console or 
install.packages("https://github.com/CCBR/l2p/raw/master/l2p_0.0-13.tar.gz", repos=NULL) 
install.packages("https://github.com/CCBR/l2p/raw/master/l2psupp_0.0-13.tar.gz", repos=NULL) 

```

**Option 3:** Download and install `l2p` using conda:
```bash
# Download Package
wget https://github.com/CCBR/l2p/blob/master/r-l2p-0.0_13-r35_0.tar.bz2?raw=true -O r-l2p-0.0_13-r35_0.tar.bz2
# Install in a conda enviroment
conda install r-l2p-0.0_13-r35_0.tar.bz2
wget https://github.com/CCBR/l2p/blob/master/r-l2psupp-0.0_13-r35_0.tar.bz2?raw=true -O r-l2psupp-0.0_13-r35_0.tar.bz2
conda install r-l2psupp-0.0_13-r35_0.tar.bz2

```

> _**Please Note:**_ It is assumed [R](https://cran.r-project.org/doc/manuals/R-admin.html) is installed on the target system. 

## Usage
```
# Available functions in l2p
l2p(genelist)               # return data frame with proabilities that arg (list of genes) matches a pathway
l2pgetlongdesc(acc)         # get the full (possibly very long) description for pathway accession identifer string
l2pgetgenes4acc(acc)        # get the list all the genes for a pathway, use the accession.

# Convenience Functions:
l2pu(list,universe)         # return data frame with probabilities with list of genes and user specified universe
l2pwcats(list,categeories)  # return data frome with categories specified
l2puwcats(list,universe,categories) # same as l2pwcats but also with a universe
l2pver                              # return l2p version


# Available functions in l2psupp 
m2h(mousegenelist)          # return list of human genes for input list of mouse gene names
a2a(genelist,fromspecies,tospecies) # return list of source species genes and return list of orthologs for destination species
updategenes(genelist , [trust=1] , [ legitonly=0 ] )   # update old gene names to current HGNC (HUGO) gene names.
egid2hugos(entrez_gene_list) # get HGNC names for entrez gene ids
 
```

The `l2p` function supopors R style arguments    :
Here is a description of each parameter: 
 - universe: list of gene names
 - categories: see categories below
 - custompathways: A list of vectors in gene matrix transpose (GMT) format. Each custom pw vector: pwname, desc, gene1, gene2...
 - customfile: a GMT file
 - universefile: list of genes one per line

**Available categories**:
 - **BIOCYC**: organism specific Pathway/ Genome Databases (PGDBs)  - https://biocyc.org/
 - **GO**: initiative to unify representation of gene and gene product attributes -  http://geneontology.org
 - **KEGG**: databases dealing with genomes, biological pathways, - https://www.kegg.jp/
 - **PANTH**: databases for protein analysis through evolutionary relationships - http://www.pantherdb.org/
 - **PID**: Pathway interaction database: legacy database from Carl Schaefer & buddies at NCI
 - **REACTOME**: curated database of biological pathways - https://reactome.org/
 - **WikiPathways**: community resource for biological pathways 
 - **C1**: MSigDB positional gene sets for each human chromosome and cytogenetic band.
 - **C2**: MSigDB curated gene sets from online pathway databases, publications in PubMed, and experts.
 - **C3**: MSigDB motif gene sets based on conserved cis-regulatory motifs from comparative analysis
 - **C4**: MSigDB computational gene sets defined by mining large collections of cancer-oriented microarray data.
 - **C6**: MSigDB oncogenic gene sets defined directly from microarray data from cancer gene perturbations.
 - **C7**: MSigDB immunologic gene sets  from microarray data from immunologic studies.
 - **C8**: MsigDB markers identified in single-cell sequencing studies of human tissuetps://www.wikipathways.org

> _**Please Note:**_ MSigDB "ARCHIVED" pathways are not provided.  MSigDB category "C5" is not there. Use "GO" category (from NCBI biosystems), instead.

#### Example function call
```R 

```

#### Output
```
The output is a data frame with the following fields ...
 

1 pathway_name                  name of pathway
2 pval                          fisher's exact p-value
3 fdr                           false discovery rate: benjamini hochberg
4 enrichment_score              same as old but multiplied by 100 : ((number_hits /(number_hits+number_misses)) - (number_user_genes/(number_user_genes+total_gens_minus_input))) * 100
5 percent_gene_hits_per_pathway (number_hits/(number_hits+number_misses))
6 number_hits                   number of genes hit in pathway
7 number_misses                 pathway number of genes in the pathway
8 number_user_genes             total count of user genes (user input)
9 total_genes_minus_input       total number of unique genes in all pathways
10 pathway_id                   canonical accession ( if available, otherwise assigned by us )
11 category                     KEGG,REACTOME,GO,PANT(=PANTHER),PID=(pathway interaction database)  *was "source"*
12 pathway_type                 functional_set,pathway,structural_complex,custom
13 genesinpathway               HUGO genes from user that hit the pathway

#### More Examples    
```R
# Example usage for l2p
    
library(l2p)
genes <- c( "TP53", "PTEN", "APC" )
x = l2p(as.vector(genes))
options(max.print=1000000)
options(width=10000)
print(x)

# How to make a custom pathway
vec1 = c("lall_ad.2","all_ad","AARS","ABCA1","ABCC9","ACTA1","ACTA2","ACTB","ACTC1","ACTG1","ACTN2","ACTN4","ACVR2B","ACVRL1","ADAR","AFG3L2","AFP","AIP","AK1","AKAP9")
vec2 = c("ACMG_2_0.2","ACMG_2_0","BRCA1","BRCA2","TP53","STK11","MLH1","MSH2","MSH6","PMS2")
vec3 = c("berg_ad.2","berg_ad","AARS","ABCC9","ACTA2","ACTB","ACTC1","ACTG1","ACTN2","ACTN4","ACVR2B","ACVRL1","ADAR","AFG3L2","AIP","AK1","AKAP9","AKT2","AMPD1","ANG","ANK2","ANKH","APC","APOA2","APOA5","APOB","APP","ATL1","ATP1A2","ATP2A2","ATP2C1","ATXN1","ATXN10","ATXN2","ATXN3","ATXN7","AXIN2","BAG3","BCO1","BEST1")
mylist <- list(vec1, vec2,vec3)
genes <- c( "TP53", "PTEN", "APC" , "CENPF" , "DLAT", "TP53" , "NOTAGENE" ,"ABCA1","ABCC9","ACTA1", "ADH1A" ,"ATXN3", "BEST1")
x = l2p(as.vector(genes),custompathways=mylist)
print(length(x))
print(x);

# How to set a user universe
library(l2p)
options(width=10000)
options(max.print=999999)

fv<-c("ADH1A","CATSPERG","HLA-DQA2","HINT2P1","MIR3150A","OR5BS1P","LINC02338","C4orf48","PARD3B","CX3CR1","RPL21P121","ARHGAP1","GAPDHP36","CNBD1","C8orf48","HTR3D","LINC00396","HIGD1AP5","C16orf90","RNU1-134P","CKAP2P1","AP5M1","FFAR3","LAD1","RNU6-524P","TJP3","JRKL","CRADD","RN7SL333P","CYP4F26P","CD1A","B3GNT5","TACC1P1","LINC02763","LOC100505664","TEX15","RPSAP18","CHP2","TRAV8-3","PFDN5","RPL7P8","SERPINA9","DNTTIP1","MELTF","HESX1","LINC02277","SFSWAP","SLC7A11","NAA16","FAM171B","GMNN","ZBTB2","WNT6","LINC02799","MRPL4","MTND1P37","HMGN2P40","NMD3P1","MIR195","LINC02785","DYM","TADA3","CEACAMP5","FAM198B-AS1","FZD8","TTC39C-AS1","RN7SL470P","IQANK1","IGKV1OR9-1","RPL10AP3","BPI","RPL5P25","CARD16","LINC02415","UBE2Q2P10","MIR6761","RNU6-903P","LINC01559","ARL17A","MIR518F","BRAP","LINC01165","XPC","RNU6-505P","LRRIQ4","MIR192","CCL27","LAPTM4BP2","INVS","TMEM161B-AS1","FAM197Y6","HSPD1","UGT1A9","TOR3A","TAF15","MIR6726","TMEM87A","HMGB1","MEI4","NAGPA-AS1","MAPK6P5","HTRA2","HSPB1P1","DYRK1A","IFFO2","TACO1","PPP6C","OR5D14","RNU6-313P","LINC01940","BBS2","RN7SL435P","LINC02422","OR3B1P","ZZEF1","EARS2","LINC02558","LINC00265-2P","KCNH1","GSTP1P1","MIR8076","RNU6-370P","RNA5SP279","RN7SL752P","CXorf49B","ANKRD36P1","IDH3A","RNU6-644P","NUCB2","CHCHD4","FAM138C","MIR198","CDC23","BRCA1","LINC02681","TFB2M","PPIP5K2","MAP2K1","MTATP6P14","COX6B1P2","HDAC5","RAB11FIP2","VSIG4","RN7SL690P","DNAJC13","GOT2P1","GTF2H1","BIRC2","LOC100132202","GAGE4","MTRNR2L10","LINC02319","C8orf49","CCNG2","LINC01524","RN7SKP49","CLDN22","FXYD6","LINC00384","ZNF14","PCGF3","CCDC6","TM4SF20","PRPS1L1","PRORSD1P","SEPHS1P1","KCNA10","MGAT5","LINC02015","BSDC1","POTEM","PHAX","RNU4-65P","MTND1P16","GPRIN2","GALE","CALY","QTRT2","RNU2-18P","TNFRSF10A-AS1","NECTIN3","RNU7-84P","PCK2","BBS5","CEACAMP4","UBE2R2","ABCB9","INTS13","ZNF69","PLEKHM2","LDHA","PHKBP1","SLC9B2","HNRNPA3P9","ARGFXP1","IER5L","CAPRIN1","RNA5SP19","NOP9","COX6CP16")

genes<-c("ADH1A","CATSPERG","HLA-DQA2","HINT2P1","MIR3150A","OR5BS1P","LINC02338","C4orf48","PARD3B","CX3CR1","RPL21P121","ARHGAP1","GAPDHP36","CNBD1","C8orf48","HTR3D","LINC00396","HIGD1AP5")

# x=l2p(genes,categories="KEGG")
x=l2p(genes,universe=fv,categories="KEGG")
print(length(x))
x=l2pgetuniverse(categories="PID")
print(length(x))
x=l2pgetuniverse()
print(length(x))
x=l2pgetuniverse(categories="KEGG")
print(length(x))
x=l2pgetuniverse(categories="C2")
print(length(x))
x=l2pgetuniverse(categories="C3")
print(length(x))
x=l2pgetuniverse(categories="C2,C3")
print(length(x))
```
## Testing  
``` bash
# Running l2p's QA test program
R --vanilla < test.R
```

## Citation
<sup>Finney, R. & Nelson, G. (2020, July 13). List-to-Pathway: an ultrafast R package for gene set enrichment analysis (Version v0.0.3). Zenodo. http://doi.org/10.5281/zenodo.3942233</sup>
  
<hr>
  
<p align="center">
	<a href="#l2p">Back to Top</a>
</p>
