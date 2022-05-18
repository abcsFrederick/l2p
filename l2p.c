
/*
 * todo  :

1) make sure custom works.
2) updategenes()
3) ensembl2gene()
4) pretty print  , pvals  if want raw data frame DOUBLE CHECK AGAINST EMAIL
5) Another column will be called “%Gene Hits per Pathway”, explained as  (A/(A+B))
6) bench marks permute vs. non . plot results.

vi l2p2.c ; make
gcc --Wall -Os -o l2p2 pwgenes.c small.c

example command line debug: cat deglist | valgrind -v --tool=memcheck --leak-check=full --show-leak-kinds=all --track-origins=yes  ./l2p -permute

NUMBER_HITS … becomes Significant genes IN Pathway (A)
NUMBER_MISSES … becomes Non-Significant genes IN Pathway (B)
NUMBER_USER_GENES … becomes Significant genes NOT IN Pathway (C)
TOTAL_GENES_MINUS_INPUT … becomes Non-Significant Genes NOT IN Pathway (D)
     
      new_name                       old_name                   description
   1  pval                           pval                       fisher's 2x2 exact test
   2  fdr                            fdr                        false discovery rate: default is benjamini hochberg, GPCC method if permute=1 
   3  enrichment_score               ratio                      same as old but multiplied by 100 : ((number_hits /(number_hits+number_misses)) - (number_user_genes/(number_user_genes+total_gens_minus_input))) * 100
** 4  percent_gene_hits_per_pathway  NEW FIELD                  (number_hits/(number_hits+number_misses))
   5  number_hits                    pwhitcount                 number of genes hit in pathway
   6  number_misses                  pwnohitcount               pathway number of genes in the pathway
   7  number_user_genes              inputcount                 total count of user genes (user input)
   8  total_genes_minus_input        pwuniverseminuslist        total number of unique genes in all pathways
   9  pathway_id                     pathwayaccessionidentifier canonical accession ( if available, otherwise assigned by us )
  10  category                       category                   KEGG,REACTOME,GO,PANT(=PANTHER),PID=(pathway interaction database)  *was "source"*
  11  pathway_name                   pathwayname                Name of pathway
  12  pathway_type                   pathwaytype                functional_set,pathway,structural_complex,custom
  13  genes                          genes_space_separated      HUGO genes from user that hit the pathway


The Fisher’s exact test on the permutation table is basically stand-alone, it doesn’t need a permutation test initially.  But it looks like I’m testing something different from what is usually tested; for my test the total number of genes in the pathway isn’t relevant.  Good, I want to do something different.  

The contingency table: So I can think clearly I’ll set this up the way I’m used to (from epidemiology, first column is no disease, second is disease; first row in no exposure, second is exposure).  For us “disease” means occurring in the pathway in question; “exposure” means being on the investigator’s list of genes. 

First without the correction:
A:  Count of genes in the universe but not on the list and not in the pathway
B:  Count of genes in the universe but not on the list, in the pathway
C:  Count of genes on the list but not in the pathway
D:  Count of genes on the list and in the pathway                                            
c d b a
a     4  pwhitcount           number of genes hit in pathway
b     5  pwnohitcount         pathway number of genes in the pathway
c     6  inputcount           total count of user genes (user input)
d     7  pwuniverseminuslist  unique genes in universe minus inputcount

Now with the correction:
B and D are unchanged.
For A and C:  Sum over the genes, the number of pathways (in the universe of pathways) each gene occurs in.
For example suppose there are 3 genes in C (must be a very short list of genes!) so in the uncorrected table the entry in C is 3.
The first gene occurs in 5 pathways, the second in 2, the third in 10.  The entry in C is 17.  
Same for A, but this will always be a sum over a large number of genes.

The sum of the entries in the table will now be larger, so the table will need to be corrected, but ignore this for now. 

original                               GPCC
a=hits                                  d
b=pathway-hits                          b
c=list-hits                             c
d=universe-pathway-list+hits            a

so ...
   UUUUUUUUUUUUUUUUUUUUUUUUUUUUU 
   U                           U 
   U universe   PPPPPPPPPPPP   U 
   U            P          P   U 
   U            P pathway  P   U 
   U            P          P   U 
   U       LLLLLLLLLLLLLL  P   U 
   U       L    P       L  P   U 
   U       L    P hits  L  P   U 
   U       L    P       L  P   U 
   U       L    PPPPPPPPLPPP   U 
   U       L            L      U 
   U       L  list      L      U 
   U       L            L      U 
   U       LLLLLLLLLLLLLL      U 
   U                           U 
   UUUUUUUUUUUUUUUUUUUUUUUUUUUUU 

   UUUUUUUUUUUUUUUUUUUUUUUUUUUUU    UUUUUUUUUUUUUUUUUUUUUUUUUUUUU    UUUUUUUUUUUUUUUUUUUUUUUUUUUUU 
   U                           U    U                           U    U                           U 
   U universe                  U    U            PPPPPPPPPPPP   U    U                           U 
   U                           U    U            P          P   U    U                           U 
   U                           U    U            P pathway  P   U    U                           U 
   U                           U    U            P          P   U    U                           U 
   U                           U    U            P          P   U    U       LLLLLLLLLLLLLL      U 
   U                           U    U            P          P   U    U       L            L      U 
   U                           U    U            P          P   U    U       L            L      U 
   U                           U    U            P          P   U    U       L            L      U 
   U                           U    U            PPPPPPPPPPPP   U    U       L            L      U 
   U                           U    U                           U    U       L            L      U 
   U                           U    U                           U    U       L  list      L      U 
   U                           U    U                           U    U       L            L      U 
   U                           U    U                           U    U       LLLLLLLLLLLLLL      U 
   U                           U    U                           U    U                           U 
   UUUUUUUUUUUUUUUUUUUUUUUUUUUUU    UUUUUUUUUUUUUUUUUUUUUUUUUUUUU    UUUUUUUUUUUUUUUUUUUUUUUUUUUUU 

A = u - p - h
B = u - l 
C = l - h 
D = h

dcba

see https://www.biostars.org/p/252937/ on how DAVID does it.

list=300 hits=3  pw=40 genome=30001
         usergenes     genome
ipw         3-1         40
notinpw     297       29960

a=hits
b=list-hits
c=list
d=u-list

problem #1

problem #2
   One sided test - count

               diseased healthy
exposed           DE       HE
notexposed        DN       HN
OR=(DE/HE)/(DN/HN)

a/b / c/d

*/
#define NELSON_C 1
#define NELSON_TEST 0

// RICH:  for me running C program:
#if NELSON_TEST
#define DEGLISTNAME "brca_up_in_meta500.txt"
#define OUTPUTFILE "test_GWN_output.txt"
#define TECHOUTPUTFILE "test_GWN_tech_output.txt"
//#define DEGLISTNAME "deglistupinreactive.txt"
// RICH: will need to be passed paramter:
#endif

#define SCALE_UNIVERSE 1
 // try 1000000 for testing 
#define NUM_PERMUTES 200000
#define NUM_TEST_PERMUTES 200000
#define SIG_P 0.05

#ifdef L2P_USING_R
#include <R.h>
#include <Rdefines.h>
#endif


#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>
#include <math.h>

#include <inttypes.h>

#include "pathworks.h"

#ifdef WEBASSEMBLY
#include "small.h"
#else

#if L2PUSETHREADS
#include <pthread.h>
#endif

#include "big.h"
#endif


extern unsigned short int pwgenes[];
extern struct smallgenetype genes[];
extern struct pwtype pws[];
extern unsigned int numpws;
extern unsigned int numgenes;
extern unsigned int numpwgenes;

struct smallgenetype *by_egids; // a copy of "genes" ,but sorted by egid (entrez gene id )

// George's globals
unsigned int ingenecnt;
static struct tree_with_count **aug_treepointers;
static struct tree_with_count **treepointers;
static struct used_path_type **sig_path_ptr = (struct used_path_type **)0;
static unsigned int deg_count = 0;
static unsigned int *gene_path_cts;
static unsigned int *deg_path_cts;
static unsigned int aug_deg_count = 0; // number of times the userinput appears in all the pathways 
static unsigned int deg_count_sum;
static unsigned int real_genect = 0;
static unsigned int aug_gene_ct = 0;
static unsigned int *auguniverse;
static double mean_deg_gpcount;
static double mean_u_gpcount;
unsigned int ingenes[MAX_INGENES];
// static unsigned int n_sig_paths = 0;   not used, for debugging
// static struct used_path_type **augreverse = (struct used_path_type **)0;
// static struct tree_with_count *head = (struct tree_with_count *)0; Got moved to stack

#if NELSON_TEST
FILE *path_gene_countsfile;
// FILE *gene_path_countsfile;
#endif



#if 0
 // debug routine
void dump_used_paths(struct used_path_type *u, int num, char *msg)
{
    struct used_path_type *uptr; // used path pointer
    int i,j;
    char *z;

fprintf(stderr,"in dump_used_paths() num=%d msg=[%s]\n",num,msg); fflush(stderr); 
    for (i=0 ; i<100 ; i++)   // for each line of custom file
    {
        uptr = (u+i);
        fprintf(stderr,"acc %d %s\n",i,uptr->acc);  fflush(stderr); 
        fprintf(stderr,"name %d %s\n",i,uptr->name); 
        fprintf(stderr,"numgenes %d %u\n",i,uptr->numgenes); 
        fprintf(stderr,"numfixed %d %u\n",i,uptr->numfixedgenes); 
        fprintf(stderr,"genehits %d %p\n",i,uptr->genehits); 
        fprintf(stderr,"egids %d %p\n",i,uptr->egids); 
        for (j=0;j<10;j++)
       // for (j=0;j<uptr->numfixedgenes;j++)
        {
            z = egid2hugo(*(uptr->egids+j));
            fprintf(stderr,"gene %d %d = %u = %s [%s]\n",i,j,*(uptr->egids+j),z,msg);
        }
        fflush(stderr);
    }
    return;
}
#endif

#if RADIX_SORT
  // from awardofsky's github  . This is faster than qsort.
void radix_ui(register unsigned int vector[], register const unsigned int size) 
{

    /* Support for variable sized integers without overflow warnings */
    const int MAX_UINT__ = ((((1 << ((sizeof(unsigned int) << 3) - 2)) - 1) << 1) + 1);
//    const int LAST_EXP__ = (sizeof(unsigned int) - 1) << 3;
    
    /* Define std preliminary, constrain and expression to check if all bytes are sorted */
#define PRELIMINARY__ 100
#define MISSING_BITS__ exp < (sizeof(unsigned int) << 3) && (max >> exp) > 0
    /* Check for biggest integer in [a, b[ array segment */
#define LOOP_MAX__(a, b)				\
    for(s = &vector[a], k = &vector[b]; s < k; ++s) {	\
	if(*s > max)  {					\
	    max = *s;					\
	}						\
    }
    
    /* b = helper array pointer ; s, k and i = array iterators */
    /* exp = bits sorted, max = maximun range in array         */
    /* point = array of pointers to the helper array           */
    register unsigned int *b, *s, *k;
    register unsigned int exp = 0;
    register unsigned int max = exp;
    unsigned int i, *point[0x100];
    int swap = 0;
    
    /* Set preliminary according to size */
    const unsigned int preliminary = (size > PRELIMINARY__) ? PRELIMINARY__ : (size >> 3);
    
    /* If we found a integer with more than 24 bits in preliminar, */
    /* will have to sort all bytes either way, so max = MAX_UINT__ */
    LOOP_MAX__(1, preliminary);
    if(max <= (MAX_UINT__ >> 7)) {	
	LOOP_MAX__(preliminary, size);
    }
    
    /* Helper array initialization */
    b = (unsigned int *)malloc(sizeof(unsigned int) * size);
    
    /* Core algorithm: for a specific byte, fill the buckets array, */
    /* rearrange the array and reset the initial array accordingly. */
#define SORT_BYTE__(vec, bb, shift)					\
    unsigned int bucket[0x100] = {0};					\
    register unsigned char *n = (unsigned char *)(vec) + (exp >> 3),*m; \
    for(m = (unsigned char *)(&vec[size & 0xFFFFFFFC]); n < m;) {	\
	++bucket[*n]; n += sizeof(int);					\
	++bucket[*n]; n += sizeof(int);					\
	++bucket[*n]; n += sizeof(int);					\
	++bucket[*n]; n += sizeof(int);					\
    }									\
    for(n = (unsigned char *)(&vec[size & 0xFFFFFFFC]) + (exp >> 3),	\
	    m = (unsigned char *)(&vec[size]); n < m;) {		\
	++bucket[*n]; n += sizeof(int);					\
    }									\
    s = bb;								\
    int next = 0;							\
    for(i = 0; i < 0x100; ++i) {					\
	if(bucket[i] == size) {						\
	    next = 1;							\
	    break;							\
	}								\
    }									\
    if(next) {								\
	exp += 8;							\
	continue;							\
    }									\
    for(i = 0; i < 0x100; s += bucket[i++]) {				\
	point[i] = s;							\
    }									\
    for(s = vec, k = &vec[size]; s < k; ++s) {				\
	*point[(*s shift) & 0xFF]++ = *s;				\
    }									\
    swap = 1 - swap;							\
    exp += 8;
    
    /* Sort each byte (if needed) */
    while(MISSING_BITS__) {
	if(exp) {
	    if(swap) {
		SORT_BYTE__(b, vector, >> exp);
	    } else {
		SORT_BYTE__(vector, b, >> exp);
	    }
	} else {
	    SORT_BYTE__(vector, b, );
	}
    }

    if(swap) {
	memcpy(vector, b, sizeof(unsigned int) * size);
    }
    
    /* Free helper array */
    free(b);
    
    /* Undefine function scoped macros for eventual later use */
#undef PRELIMINARY__
#undef MISSING_BITS__
#undef LOOP_MAX__
#undef SORT_BYTE__
    
}
#endif


 // seems to be a bug in bsearch2 , valgrind reports invalid access. so user c stdlib bsearch for now
int bsearch2(const unsigned int key, const unsigned int  *base, size_t nmemb) 
{
    unsigned int  current_element;
    int medium;
    int first = 0;
    int last = nmemb;

    while (first <= last) 
    {
        medium = first + (last - first) / 2;
        current_element = *(base + medium);
        if (key < current_element)
            last = medium - 1;
        else if (key > current_element)
            first = medium + 1;
        else
            return 1;
    }
    return 0;
}


#if 0
int binsearch_6( unsigned int  array[], size_t size, unsigned int key, size_t *index )
{
  if( !array || !size ) return 0;
  arr_t *p=array;
  switch( size )
  {
#define C(n) case ((size_t)1<<n)-1: if( p[1<<(n-1)]<=key ) p+=1<<(n-1); 
#if SIZE_MAX == UINT64_MAX 
                                              C(63) C(62) C(61)
    C(60) C(59) C(58) C(57) C(56) C(55) C(54) C(53) C(52) C(51)
    C(50) C(49) C(48) C(47) C(46) C(45) C(44) C(43) C(42) C(41)
    C(40) C(39) C(38) C(37) C(36) C(35) C(34) C(33) C(32)
#endif 
                                                          C(31)
    C(30) C(29) C(28) C(27) C(26) C(25) C(24) C(23) C(22) C(21)
    C(20) C(19) C(18) C(17) C(16) C(15) C(14) C(13) C(12) C(11)
    C(10) C( 9) C( 8) C( 7) C( 6) C( 5) C( 4) C( 3) C( 2) C( 1)
#undef C 
      break;
    default:
      while( size > 0 ){
        size_t w=size/2;
        if( p[w] < key ){ p+=w+1; size-=w+1; } else size=w;
      }
  }
  *index=p-array; return p[0]==key;
}
#endif



char *type2string(int type)
{
   static char functional_set[] = "functional_set";
   static char pathway  [] = "pathway";
   static char structural_complex  [] = "structural_complex";
   static char custom_string  [] = "custom";
   static char null_info  [] = "NA";

   if (type == type_functional_set) return &functional_set[0];
   else if (type == type_pathway) return &pathway[0];
   else if (type == type_structural_complex) return &structural_complex[0];
   else if (type == type_custom) return &custom_string[0];
   return &null_info[0];
}

int cmp_hugo(const void *a, const void *b)
{
    return strcmp(((struct smallgenetype *)a)->hugo, ((struct smallgenetype *)b)->hugo);
}

int cmp_by_egid(const void *a, const void *b) // compare entrez gene id
{
    if      ( ((struct smallgenetype *)a)->egid < ((struct smallgenetype *)b)->egid ) return -1;
    else if ( ((struct smallgenetype *)a)->egid > ((struct smallgenetype *)b)->egid ) return 1;
    return 0;
}

char *egid2hugo(int egid)
{
    struct smallgenetype glocal;
    struct smallgenetype *gptr;

    glocal.egid = egid;
    gptr = (struct smallgenetype *)bsearch(&glocal,by_egids,numgenes,sizeof(struct smallgenetype),cmp_by_egid);
    if (gptr) return gptr->hugo;
    else return (void *)0;
}

unsigned int hugo2egid(char *h)
{
    struct smallgenetype glocal;
    struct smallgenetype *gptr;

    glocal.hugo = h;
    gptr = (struct smallgenetype *)bsearch(&glocal,genes,numgenes,sizeof(struct smallgenetype),cmp_hugo);
    if (gptr) return gptr->egid;
    else return (unsigned int)UINT_MAX;
}

struct smallgenetype *hugo2geneptr(char *h)
{
    struct smallgenetype glocal;
    struct smallgenetype *gptr;

    glocal.hugo = h;
    gptr = (struct smallgenetype *)bsearch(&glocal,genes,numgenes,sizeof(struct smallgenetype),cmp_hugo);
    return gptr;
}

unsigned short int hugo2usedindex(char *h)
{
    struct smallgenetype glocal;
    struct smallgenetype *gptr;
    unsigned short int ret;
    int idx;

    glocal.hugo = h;
    gptr = (struct smallgenetype *)bsearch(&glocal,genes,numgenes,sizeof(struct smallgenetype),cmp_hugo);
    if (gptr)
    {
        idx = gptr - &genes[0];
        if (idx > 20000) {fprintf(stderr,"ERROR in hub2usedindex %s gptr=%p idx=%d\n",h,gptr,idx); fflush(NULL); exit(0); }
        ret = (unsigned short int)idx;
    }
    else
    {
        ret = (unsigned short int)USHRT_MAX;
    }
    return ret;
}

int cmp_ui(const void *a, const void *b)
{
    if      ( *(unsigned int *)a < *(unsigned int *)b ) return -1;
    else if ( *(unsigned int *)a > *(unsigned int *)b ) return 1;
    return 0;
}

int cmp_double(const void *a, const void *b)
{
    if      ( *(double *)a < *(double *)b ) return -1;
    else if ( *(double *)a > *(double *)b ) return 1;
    return 0;
}

#if 0
int cmp_usi(const void *a, const void *b)
{
    if      ( *(unsigned short int *)a < *(unsigned short int *)b ) return -1;
    else if ( *(unsigned short int *)a > *(unsigned short int *)b ) return 1;
    return 0;
}


 not used
int cmp_float(const void *a, const void *b)
{
    if      ( *(float *)a < *(float *)b ) return -1;
    else if ( *(float *)a > *(float *)b ) return 1;
    return 0;
}
#endif

#if 0
/*
ReservoirSample(S[1..n], R[1..k])
  R[1] := S[1]
  for i from 2 to k do
      j := randomInteger(1, i)  // inclusive range
      R[i] := R[j]
      R[j] := S[i]
  for i from k + 1 to n do
      j := randomInteger(1, i)  // inclusive range
      if (j <= k)
          R[j] := S[i]
*/

static inline void reservoir(unsigned int s[], int n, unsigned int r[],int k)
{
    int i,j;

    r[0] = s[0];
    for (i=1;i<k;i++)
    {
        j = rand() % i;
        r[i] = r[j];
        r[j] = s[i];
    }
    for (i=k+1;i<n;i++)
    {
        j = rand() % i;
        if (j <= k)
            r[j] = s[i];
    }
    return;
}
#endif

// -- for benjamini hochberg FDR ...
struct ordertype
{
    double val;
    int order;
};

static int cmp_ordertype_by_order(const void *a, const void *b)
{
    struct ordertype *aa;
    struct ordertype *bb;
    aa = (void *)a;
    bb = (void *)b;
    if      (aa->order >  bb->order) return 1;
    else if (aa->order <  bb->order) return -1;
    return 0;
}

int cmp_ordertype_by_val_REV(const void *a, const void *b)
{
    struct ordertype *aa;
    struct ordertype *bb;
    aa = (void *)a;
    bb = (void *)b;
    if      (aa->val <  bb->val) return 1;
    else if (aa->val >  bb->val) return -1;
    // if (aa>bb) return 1; else if (aa<bb) return -1;
    return 0;
}

static void benjaminihochberg(int n,double pvals[], double returnpvals[])
{
/*
 here's the code from R that I re-imagined
                   i <- lp:1L
                   o <- order(p, decreasing = TRUE)
                   ro <- order(o)
                   pmin(1, cummin( n / i * p[o] ))[ro]
*/
    int j,k;
    struct ordertype *i;
    struct ordertype *o;
    struct ordertype *po;
    struct ordertype *cummin;
//    struct ordertype *ro;
//    struct ordertype *intermed;

// fprintf(stderr,"rpf in benjaminihochberg\n"); fflush(stderr); 

    i = (struct ordertype *)malloc(sizeof(struct ordertype)*n);
    for (k=n,j=0;j<n;j++,k--) (i+j)->order=k;

#define RDEBUG 0
#if RDEBUG
FILE *fp;
fp = fopen("test.pvals","w");
#endif
    o = (struct ordertype *)malloc(sizeof(struct ordertype)*n);
    for (j=0 ; j<n ; j++)
    {
#if RDEBUG
fprintf(fp,"%20.18f\n",pvals[j]);
#endif
        (o+j)->val=pvals[j];
        (o+j)->order=j+1;
    }
#if RDEBUG
fclose(fp);
#endif
    qsort(o,n,sizeof(struct ordertype),cmp_ordertype_by_val_REV);

#if 0
    ro = (struct ordertype *)malloc((sizeof(struct ordertype))*n);
    for (j=0;j<n;j++)
    {
        (ro+j)->val = (double)(o+j)->order;
        (ro+j)->order = j+1;
    }
    qsort(ro,n,sizeof(struct ordertype),cmp_ordertype_by_val);
#endif
    po = (struct ordertype *)malloc(sizeof(struct ordertype)*n);
    memset(po,0,sizeof(struct ordertype)*n);
    for (j=0;j<n;j++)
    {
        (po+j)->val = (double)pvals[j];
        (po+j)->order = (o->order); // why the hell isn't this ro? what the what?
    }
    qsort(po,n,sizeof(struct ordertype),cmp_ordertype_by_val_REV); // == p[o]

    cummin = (struct ordertype *)malloc((sizeof(struct ordertype))*n); // holds n / i * po
    for (j=0;j<n;j++)
    {
        (cummin+j)->val = (double)n / (double)(i+j)->order * ((po+j)->val) ;
    }
                   // Rcode: pmin(1, cummin( n / i * p[o] ))[ro]           ******************
    for (j=1;j<n;j++)
    {
        if ((cummin+j)->val > (cummin+j-1)->val)
            (cummin+j)->val = (cummin+j-1)->val;
    }
    for (j=0;j<n;j++)
    {
        if ((cummin+j)->val > 1)
            (cummin+j)->val = 1;
        (cummin+j)->order = (o+j)->order ;
    }
    qsort(cummin,n,sizeof(struct ordertype),cmp_ordertype_by_order);
#if RDEBUG
FILE *fp2;
fp2 = fopen("test.fdrs","w");
#endif
    for (j=0;j<n;j++)
    {
        returnpvals[j] = (cummin+j)->val;
#if RDEBUG
fprintf(fp2,"%20.18f\n",returnpvals[j]);
#endif
    }
#if RDEBUG
fclose(fp2);
#endif
    if (i) free(i);
    if (o) free(o);
    if (po) free(po);
    if (cummin) free(cummin);

    return;
}



#if 1
inline void subsamp(unsigned int s[], int n, int k)
{
    int i,j,dncnt;
    unsigned int  ui;

    for (dncnt=n,i=0;i<k;i++,dncnt--)
    {
        j = (rand() % dncnt) + i;
        ui = s[i];
        s[i] = s[j];
        s[j] = ui;
    }
    return;
}

void shuffle(unsigned int s[], int n)
{
    int i,j;
    unsigned int ui;

// for i from 0 to n−2 do
//      j ← random integer such that i ≤ j < n
//      exchange a[i] and a[j]
// or
// for i from n−1 downto 1 do
//      j ← random integer such that 0 ≤ j ≤ i
//      exchange a[j] and a[i]

    for (i=n-1;i>=1;i--)
    {
#ifdef L2P_USING_R
        j = (int)(unif_rand() * (double)i); // unif_rand() appears to return between 0.0 and 1.0
#else
        j = rand() % i;
#endif
        ui = s[i];
        s[i] = s[j];
        s[j] = ui;
    }
    return;
}



#if 0
int permute_test(struct used_path_type used_paths[],unsigned int num_used_paths,unsigned int *real_universe,unsigned int real_universe_cnt, int incnt)
{
    unsigned int r[40002];
    double   p[NUM_PERMUTES];
    struct used_path_type *uptr; // used path pointer 
    int j,k;
    unsigned int ui;
    unsigned int usi_j;
    double d;
    int tmphitcnt;
    unsigned int ui_uk,ui_ej;
    int kickat;

    
    for (ui=0 ; ui<num_used_paths ; ui++)
    {
        uptr = &used_paths[ui];
        uptr->pval4 = uptr->fdr4 = (double)1.0;
        if (uptr->pval == 1.0)
        {
            uptr->pval4 = 1.0;
            continue;
               // don't bother 
        }
        memset(p,0,sizeof(p));
        memcpy(r,real_universe,real_universe_cnt*sizeof(unsigned int));
        radix_ui(r,ingenecnt);
        shuffle(r,real_universe_cnt);
//        subsamp(r, real_universe_cnt,ingenecnt); // inlined if optimization on 
#if 0
    for (k=0;k<ingenecnt;k++)
    {
    fprintf(stderr,"k=%d %u ",k,*(r+k));
            z = egid2hugo(*(r+k));
            if (z) fprintf(stderr,"%s",z);
            else   fprintf(stderr,"ERROR");
            fprintf(stderr,"\n");
    }
#endif
        for (j=0 ; j<NUM_PERMUTES ; j++)
        {
            usi_j = k = tmphitcnt = 0;
            while ((usi_j<uptr->numfixedgenes) && (k < ingenecnt))
            {
                ui_ej = *(uptr->egids+usi_j);
                ui_uk = *(r+k);
                if (ui_ej == ui_uk)
                {
#if 0
    char *z;
            z = egid2hugo(ui_ej);
            if (z) fprintf(stderr,"got %s %d",z,ui_ej);
            else   fprintf(stderr,"ERROR");
            fprintf(stderr,"\n");
#endif
                    tmphitcnt++;
                    k++;
                    usi_j++;
                    continue;
                }
                else if (ui_ej < ui_uk) usi_j++;
                else                    k++;
            }
            if (tmphitcnt == 0) d = 1.0;
            else 
            {
                d = exact22((int)tmphitcnt,uptr->numfixedgenes-tmphitcnt,ingenecnt,real_universe_cnt-ingenecnt);
//                fprintf(stderr,"ex22=%f %d %d %d %d\n",d,(int)tmphitcnt,uptr->numfixedgenes-tmphitcnt,ingenecnt,real_universe_cnt-ingenecnt);
            }
            p[j] = d;
 if ((d < 0.05) && (d < uptr->pval))
 {
 fprintf(stderr,"upv:%f p[%d]=%f tmphitcnt=%d %d %d %d %d %s\n",uptr->pval,j,d,tmphitcnt,(int)tmphitcnt,uptr->numfixedgenes-tmphitcnt,ingenecnt,real_universe_cnt-ingenecnt,uptr->name); 
        for (k=0 ; k<NUM_PERMUTES ; k++) fprintf(stderr,"%f ",p[k]);
        fprintf(stderr,"\n"); 
 }
        }
        qsort(p,NUM_PERMUTES,sizeof(double),cmp_double);
        if (p[0] == 1.0)
        {
            uptr->pval4 = 1.0;
            continue;
        }
        kickat = NUM_PERMUTES;
        for (j=0;j<NUM_PERMUTES;j++) 
        {
            if (p[j] > uptr->pval) 
            { 
                kickat = j;
                break; 
            }
        }
fprintf(stderr,"kick at %d d: %f pv= %f\n",kickat,d,uptr->pval); 
        d = (double)kickat/(double)NUM_PERMUTES;
        uptr->pval4 = d;
 if (test9999)
 {
 fprintf(stderr,"test9999\n"); 
         for (j=0 ; j<NUM_PERMUTES ; j++) fprintf(stderr,"%f ",p[j]);
         fprintf(stderr,"\n"); 
 }
    // fprintf(stderr,"pv4=%f\n",d);
    }
    return 0;
}
#endif

struct otype  // order type
{
    unsigned int val;
    int order;
};

void FisherYates(unsigned int *p, int n) 
{ //implementation of Fisher Yates shuffle
    int i, j;
    unsigned int tmp;

    for (i = n - 1; i > 0; i--) 
    {
        j = random() % (i + 1); //randomize
        tmp = p[j];
        p[j] = p[i];
        p[i] = tmp;
    }
    return;
}

#if 1
int permute2(struct used_path_type used_paths[],unsigned int num_used_paths, unsigned int incnt, unsigned int num_permutes)
{
    unsigned int *tmphits = (unsigned int *)0;
    unsigned int *hyper_verse = (unsigned int *)0;
    unsigned int num_hyper_genes;
    struct used_path_type *uptr; // used path pointer 
    unsigned int i,j,k;
    unsigned int usi_j;
    unsigned int ll,ui;
    int tmphitcnt;
    unsigned int ui_uk,ui_ej;
    size_t sizet;
    int ret = -1;
    unsigned int perm_count;
    unsigned int sofar ;
    double d;
    unsigned int *zhits = (unsigned int *)0; // [num_used_paths][num_permutes];
    // unsigned int r[40002];
//    char *z;


fprintf(stderr,"rpf permute2 in permute2 incnt=%d num_use_paths=%u, num_permutes=%u\n",incnt,num_used_paths,num_permutes);  fflush(stderr);
    for (i=0 ; i<num_used_paths ; i++)
    {
        uptr = (used_paths+i);
        if (!uptr->egids) continue;
        j = k = ll = 0;
        while ((j<uptr->numfixedgenes) && (k < incnt))
        {
            ui_ej = *(uptr->egids+j);
            ui = ingenes[k];
            if (ui_ej == ui)
            {
                *((uptr->genehits) + (ll++)) = ui; // remember, because need to print out later
                k++;
                j++;
                // aug hit count for aug contingency table
                continue;
            }
            else if (ui_ej < ui) j++;
            else                 k++;
        }
        uptr->hitcnt = ll;
    }

    for (num_hyper_genes=i=0 ; i<num_used_paths ; i++)
        num_hyper_genes = num_hyper_genes + used_paths[i].numfixedgenes;    // get num_hyper_genes

fprintf(stderr,"rpf permute2 num_hyper_genes=%d\n",num_hyper_genes);  fflush(stderr);
    hyper_verse = malloc(num_hyper_genes*(sizeof (unsigned int)));
    if (!hyper_verse) goto PERM2_END;

    zhits = malloc(num_used_paths * num_permutes * sizeof (unsigned int) );
    if (!zhits) goto PERM2_END;

    sizet = sizeof(unsigned int) * incnt;
    tmphits = malloc(sizet);
    if (!tmphits) goto PERM2_END;

    for (j=i=0 ; i<num_used_paths ; i++)
    {   // get hyper_verse
        for (k=0 ; k < used_paths[i].numfixedgenes ; k++)
             *(hyper_verse + j++) = *(used_paths[i].egids+k);
    }
fprintf(stderr,"rpf permute2 hyper_verse setup , shuffle next , num_hyper_genes=%u\n",num_hyper_genes);  fflush(stderr);

#if 0
FILE *fp;
fp=fopen("test7","w");
for (i=0;i<num_hyper_genes;i++)
fprintf(fp,"%u\n",*(hyper_verse+i));
fclose(fp);
#endif

    FisherYates(hyper_verse,num_hyper_genes);

#if 0
fp=fopen("test8","w");
for (i=0;i<num_hyper_genes;i++)
fprintf(fp,"%u\n",*(hyper_verse+i));
fclose(fp);
#endif
    for (perm_count=0 ; perm_count < num_permutes ; perm_count++)
    {
if ((perm_count%1000)==0) { fprintf(stderr,"rpf permute2 loop %d\n",perm_count);  fflush(stderr); }
        for (sofar=0 ; sofar < incnt ; )          // get same number of unique genes and user incnt (input list) 
        {
            k = (unsigned int)(rand() % num_hyper_genes); 
            k = *(hyper_verse + k);
#if 0
fprintf(stderr,"rand = %u of %u\n",k,num_hyper_genes); 
   z = egid2hugo(k);
   if (z) fprintf(stderr,"got %s %d hits=%d",z,k,tmphitcnt);
   else   fprintf(stderr,"ERROR");
   fprintf(stderr,"\n");
#endif

            for (j=0 ; j<sofar ; j++)
                 if ((*tmphits+j) == k) break;
            if (j < sofar) // already there
               continue;
            *(tmphits+sofar) = k;
            sofar++;
        }
        qsort(tmphits,incnt,sizeof(unsigned int),cmp_ui);
        for (i=0 ; i<num_used_paths ; i++)
        {
            uptr = &used_paths[i];
            uptr->pval4 = uptr->fdr4 = (double)1.0;
            usi_j = k = tmphitcnt = 0;
            while ((usi_j<uptr->numfixedgenes) && (k < incnt))
            {
                ui_ej = *(uptr->egids+usi_j);
                ui_uk = *(tmphits+k);
                if (ui_ej == ui_uk)
                {
                    tmphitcnt++;
#if 0
            z = egid2hugo(ui_uk);
            if (z) fprintf(stderr,"got %s %d hits=%d",z,ui_ej,tmphitcnt);
            else   fprintf(stderr,"ERROR");
            fprintf(stderr,"\n");
#endif
                    k++;
                    usi_j++;
                    continue;
                }
                else if (ui_ej < ui_uk) usi_j++;
                else                    k++;
            }
            *(zhits + (i * num_permutes) + perm_count) = tmphitcnt;
        }
     }
     for (i=0 ; i<num_used_paths ; i++)
     {
        uptr = &used_paths[i];
        for (tmphitcnt=j=0 ; j < num_permutes ; j++)
        {
            k = *(zhits + (i * num_permutes) + j );
            if ( uptr->hitcnt > k ) tmphitcnt++;
            // else if ( uptr->hitcnt == k ) tmphitcnt++;
        }
// 3 5 3 3 3 0 10 5 0 0 0 0
// 2 100
        d = 1.0 - ((double)tmphitcnt/(double)num_permutes);
// fprintf(stderr,"kick tmphitcnt %d d: %f upval: %f u-hit= %u of genesinpw= %d mn= %u mx= %u avg=%f %s\n",tmphitcnt,d,uptr->pval,uptr->hitcnt,uptr->numfixedgenes,mmin,mmax,(double)sum/(double)num_permutes,uptr->name); 
        uptr->pval4 = d;
    }

    ret = 0;
PERM2_END:
    if (zhits)      free(zhits);
    if (tmphits)    free(tmphits);
    if (hyper_verse)free(hyper_verse);
    return ret;
}

int permute3(struct used_path_type used_paths[],unsigned int num_used_paths, unsigned int incnt, unsigned int num_permutes)
{
    // unsigned int r[40002];
    unsigned int *tmphits = (unsigned int *)0;
    unsigned int *hyper_verse = (unsigned int *)0;
    unsigned int num_hyper_genes;
    struct used_path_type *uptr; // used path pointer 
    unsigned int i,j,k;
    unsigned int usi_j;
    unsigned int tmphitcnt;
    unsigned int ui_uk,ui_ej;
    size_t sizet;
    int ret = -1;
    unsigned int perm_count;
    double d;
    unsigned int *zgenes = (unsigned int *)0; // [num_used_paths][incnt];
    unsigned int morecnt;
    // char *z;

fprintf(stderr,"rpf permute3 in permute3 start incnt=%d num_use_paths=%u, num_permutes=%u\n",incnt,num_used_paths,num_permutes);  fflush(stderr);

    for (num_hyper_genes=i=0 ; i<num_used_paths ; i++)
        num_hyper_genes = num_hyper_genes + used_paths[i].numfixedgenes;    // get num_hyper_genes

fprintf(stderr,"rpf permute3 num_hyper_genes=%d\n",num_hyper_genes);  fflush(stderr);
    hyper_verse = malloc(num_hyper_genes*(sizeof (unsigned int)));
    if (!hyper_verse) goto PERM3_END;

    zgenes = malloc(num_permutes * incnt * sizeof (unsigned int) );
    if (!zgenes) goto PERM3_END;

    sizet = sizeof(unsigned int) * incnt;
    tmphits = malloc(sizet);
    if (!tmphits) goto PERM3_END;

    for (j=i=0 ; i<num_used_paths ; i++)
    {   // get hyper_verse
        for (k=0 ; k < used_paths[i].numfixedgenes ; k++)
             *(hyper_verse + j++) = *(used_paths[i].egids+k);
    }
fprintf(stderr,"rpf permute3 hyper_verse setup , shuffle next , num_hyper_genes=%u\n",num_hyper_genes);  fflush(stderr);

#if 0
FILE *fp;
fp=fopen("test7","w");
for (i=0;i<num_hyper_genes;i++)
fprintf(fp,"%u\n",*(hyper_verse+i));
fclose(fp);
#endif

    FisherYates(hyper_verse,num_hyper_genes);

#if 0
fp=fopen("test8","w");
for (i=0;i<num_hyper_genes;i++)
fprintf(fp,"%u\n",*(hyper_verse+i));
fclose(fp);
#endif
fprintf(stderr,"rpf permute3 before loop num_permutes=%u\n",num_permutes);  fflush(stderr);
    for (perm_count=0 ; perm_count < num_permutes ; perm_count++)
    {
        i = 0;
        while (i < incnt )
        {
            k = (unsigned int)(rand() % num_hyper_genes); 
            k = *(hyper_verse + k);
            for (j=0 ; j<i ; j++)
                 if ( (*(zgenes+(perm_count*incnt)+j)) == k) break;
            if (i==j)
            {
                *(zgenes+(perm_count*incnt)+i) = k;
                i++;
            }
        }
        qsort(zgenes+(perm_count*incnt),incnt,sizeof(unsigned int),cmp_ui);
    }
fprintf(stderr,"rpf permute3 after  loop num_permutes=%u\n",num_permutes);  fflush(stderr);

if ((perm_count%1000)==0) { fprintf(stderr,"rpf permute3 loop %d\n",perm_count);  fflush(stderr); }
    for (i=0 ; i<num_used_paths ; i++)
    {
// fprintf(stderr,"rpf permute3 in 2nd loop i=%u to %u\n",i,num_used_paths);  fflush(stderr);
         morecnt = 0;
         uptr = &used_paths[i];
         uptr->pval4 = uptr->fdr4 = (double)1.0;
         for (perm_count=0 ; perm_count < num_permutes ; perm_count++)
         {
            usi_j = k = tmphitcnt = 0;
            while ((usi_j<uptr->numfixedgenes) && (k < incnt))
            {
                ui_ej = *(uptr->egids+usi_j);
                ui_uk = *(zgenes+(perm_count*incnt)+k);
                if (ui_ej == ui_uk)
                {
                    tmphitcnt++;
#if 0
            z = egid2hugo(ui_uk);
            if (z) fprintf(stderr,"got %s %d hits=%d",z,ui_ej,tmphitcnt);
            else   fprintf(stderr,"ERROR");
            fprintf(stderr,"\n");
#endif
                    k++;
                    usi_j++;
                    continue;
                }
                else if (ui_ej < ui_uk) usi_j++;
                else                    k++;
            }
            if (tmphitcnt > uptr->hitcnt)
            {
                morecnt++;
            }
        }
        d = ((double)morecnt/(double)num_permutes);
// 3 5 3 3 3 0 10 5 0 0 0 0
// 2 100
        // d = 1.0 - ((double)morecnt/(double)num_permutes);
// fprintf(stderr,"kick tmphitcnt %d d: %f upval: %f u-hit= %u of genesinpw= %d mn= %u mx= %u avg=%f %s\n",tmphitcnt,d,uptr->pval,uptr->hitcnt,uptr->numfixedgenes,mmin,mmax,(double)sum/(double)num_permutes,uptr->name); 
        uptr->pval4 = d;
    }
    ret = 0;
PERM3_END:
    if (tmphits)    free(tmphits);
    if (hyper_verse)free(hyper_verse);
    return ret;
}
#endif

#endif



void test22(int a,int b,int c,int d,double pv)
{
    FILE *fp;
    char cmd[512];
    char s[512];
    char junks[512];
    char junk2[512];
    char p[512];
    double diff;
    
    fp = fopen("tmp.R","w");
    
    fprintf(fp,"data = matrix(c(%d,%d,%d,%d), nrow = 2)\n",a,b,c,d);
//    fprintf(fp,"fisher.test(data,alternative=\"two.sided\")\n");
    fprintf(fp,"fisher.test(data,alternative=\"greater\")\n");
//    fprintf(fp,"fisher.test(data,alternative=\"less\")\n");
//    fprintf(fp,"library(\"exact2x2\")\n");
//    fprintf(fp,"exact2x2(data,alternative=\"two.sided\")\n");
//    fprintf(fp,"exact2x2(data,alternative=\"greater\")\n");
//    fprintf(fp,"exact2x2(data,alternative=\"less\")\n");
    fclose(fp);

    strcpy(cmd,"R -q --vanilla < tmp.R");
    fp = popen(cmd,"r");
    while (fgets(s,sizeof(s),fp))
    {
        if (strstr(s,"p-value"))
        {
            sscanf(s,"%s %s %s\n",junks,junk2,p);
            diff = pv - atof(p);
            if (diff < 0) diff = 0 - diff;
// fprintf(stderr,"%d %d %d %d %f %s %f\n",a,b,c,d,pv,p,diff);
        }
    }
    fclose(fp);
    return;
}


int do_pvals_and_bh(unsigned int incnt, struct used_path_type usedpaths[], unsigned int num_used_paths,unsigned int real_universe_cnt,int oneside)
{
    double pv;
    double pv2;
    double enrichment_score;
    unsigned int i;
    double *pvals;
    double *fdrs;
    unsigned int localhitcnt ;
    struct used_path_type *uptr;	// used path pointer 
    int a2 , b2 , c2 , d2;


// fprintf(stderr,"debug in do_pvals_and_bh() incnt=%d num_used_paths=%d real_universe_cnt=%d oneside=%d\n",incnt,num_used_paths,real_universe_cnt,oneside); fflush(stderr); 

    pv = 1.0;
// fprintf(stderr,"in do_pvals_and_bh() incnt=%d num_used_paths=%d real_universe_cnt=%d 2\n",incnt,num_used_paths,real_universe_cnt); fflush(stderr); 
    pvals = (double *)malloc((size_t)(sizeof (double)*(num_used_paths)));
// fprintf(stderr,"in do_pvals_and_bh() 3 pvals=%p\n",pvals); fflush(stderr); 
    if (!pvals) { fprintf(stderr,"ERROR: no memory\n"); fflush(stderr); return -1; }
// fprintf(stderr,"in do_pvals_and_bh() before loop to %d (num_used_paths)\n",num_used_paths); fflush(stderr); 
    for (i=0 ; i<num_used_paths;i++)
    {
        pv = pv2 = 1.0; // initialize 
// fprintf(stderr,"in do_pvals_and_bh() path=%d of %d\n",i,num_used_paths); fflush(stderr); 
        uptr = (usedpaths+i);
        localhitcnt = uptr->hitcnt;
        if (localhitcnt)
        {
            if (localhitcnt > (unsigned int)uptr->numfixedgenes) 
            {
                 fprintf(stderr,"ERROR: more hits than genes (localhitcnt %d > %d for %s)\n", 
                         localhitcnt , uptr->numfixedgenes,uptr->acc);
                 fprintf(stderr,"details: index=%d , hitcnt=%u\n", i,localhitcnt); 
                 fflush(stderr);
                 return -2;
            }
        }
/*
     4  pwhitcount            number of genes hit in pathway
     5  pwnohitcount          pathway    number of genes in the pathway
     6  inputlist                  total count of user genes (user input)
     7  pwuniverseminuslist    tinlistpathwaysuniverse    total number of unique genes in all pathways
a=hits
b=thispw-hits
c=list
d=u-list
        a2 = localhitcnt;
        b2 = uptr->numfixedgenes - localhitcnt;
        c2 = incnt;
        d2 = real_universe_cnt - incnt;
*/
a2 = localhitcnt;
b2 = uptr->numfixedgenes - localhitcnt;
c2 = incnt-localhitcnt;
d2 = real_universe_cnt - ( uptr->numfixedgenes ) - incnt + localhitcnt;

        if (oneside)
        {
#if 0
            pv = exact22_oneside((uint32_t)a, (uint32_t)b, (uint32_t)c, (uint32_t)d,0); // last arg is debug flag
 // testing 
double pv2 = exact22(a,b,c,d);
fprintf(stderr,"after exact22_oneside pv1 = %20.17f , pv2 = %20.17f %d %d %d %d\n",pv,pv2,a,b,c,d); fflush(NULL);
// test22(a,b,c,d,pv);
#endif
#if 1
/*
int u,p,h,l;
u = real_universe_cnt;
p = uptr->numfixedgenes;
h = localhitcnt;
l = incnt;
int a2=hits
int b2=pathway - hits
int c2=list - hits
int d2=universe - pathway - list + hits
*/
a2 = localhitcnt;
b2 = uptr->numfixedgenes - localhitcnt;
c2 = incnt-localhitcnt;
d2 = real_universe_cnt - ( uptr->numfixedgenes ) - incnt + localhitcnt;
            pv2 = exact22_oneside((uint32_t)a2, (uint32_t)b2, (uint32_t)c2, (uint32_t)d2,0); // last arg is debug flag
            pv = pv2;
#if 0
 if (strcmp(uptr->acc,"ko00010") == 0)
 {
 fprintf(stderr,"ko00010 %u %u %u %u  %18.16f numfixedgenes=%u\n",a2,b2,c2,d2,pv2,uptr->numfixedgenes); 
 }
// fprintf(stderr,"dbg pv = %f pv2=%f %d %d %d %d %d %d %d %d %s\n",pv,pv2,u,p,h,l,a2,b2,c2,d2,uptr->acc); 
#endif
#endif
        }
        else
        {
            pv = exact22(a2,b2,c2,d2);
        }
        uptr->pval2 = *(pvals+i) = pv2;
        uptr->pval = *(pvals+i) = pv;
        uptr->a = a2;
        uptr->b = b2;
        uptr->c = c2;
        uptr->d = d2;
        // enrichment_score = ( ((double)localhitcnt / (double)uptr->numfixedgenes) - ((double)incnt/(double)real_universe_cnt) );
        enrichment_score = ( ((double)a2 / (double)(a2+b2)) - ((double)c2/(double)(c2+d2)) );
        enrichment_score *= 100.0;
        uptr->enrichment_score = enrichment_score;
// if (strcmp(uptr->acc,"ko00010") == 0) { fprintf(stderr,"adbg : do_pvals_and_bh ko00010 enrichment_score is %f from %d %d %d %d\n",enrichment_score,a2,b2,c2,d2); }
    }
#if 0
for (i=0 ; i<num_used_paths ; i++)
{
 fprintf(stderr,"rpf checking 1 *(pvals+i) %d  %f\n",i,*(pvals+i) ); 
}
#endif

    fdrs = (double *)malloc((size_t)(sizeof (double)*num_used_paths)); 
    if (!fdrs) { free(pvals); /* clean up */ fprintf(stderr,"ERROR: no memory in do_pvals_and_bh() 2\n"); return -3; }

#if 0
for (i=0 ; i<num_used_paths ; i++)
{
fprintf(stderr,"rpf checking 2 *(pvals+i) %d  %f\n",i,*(pvals+i) ); 
}
#endif
// fprintf(stderr,"rpf before benjaminihochberg\n"); fflush(stderr); 
    benjaminihochberg(num_used_paths,pvals,fdrs);
// fprintf(stderr,"rpf after benjaminihochberg\n"); fflush(stderr); 
    for (i=0 ; i<num_used_paths ; i++)
    {
         usedpaths[i].fdr = *(fdrs+i);
// fprintf(stderr,"%d %f\n", i,usedpaths[i].fdr);
    }
    free(pvals);
    free(fdrs);
// fprintf(stderr,"in do_pvals_and_bh() end\n"); fflush(stderr); 
    return 0;
}


int do_pvals_and_bh_gpcc(unsigned int incnt, struct used_path_type usedpaths[], unsigned int num_used_paths,unsigned int real_universe_cnt,int oneside)
{
    double pv;
    double pv2;
    double enrichment_score;
    unsigned int i;
    double *pvals;
    double *fdrs;
    unsigned int localhitcnt ;
    struct used_path_type *uptr;	// used path pointer 
    int a2 , b2 , c2 , d2;


    pv = 1.0;
    pvals = (double *)malloc((size_t)(sizeof (double)*(num_used_paths)));
    if (!pvals) { fprintf(stderr,"ERROR: no memory\n"); fflush(stderr); return -1; }
    for (i=0 ; i<num_used_paths;i++)
    {
        pv = pv2 = 1.0; // initialize 
        uptr = (usedpaths+i);
        localhitcnt = uptr->hitcnt;
        if (localhitcnt)
        {
            if (localhitcnt > (unsigned int)uptr->numfixedgenes) 
            {
                 fprintf(stderr,"ERROR: more hits than genes (localhitcnt %d > %d for %s)\n", 
                         localhitcnt , uptr->numfixedgenes,uptr->acc);
                 fprintf(stderr,"details: index=%d , hitcnt=%u\n", i,localhitcnt); 
                 fflush(stderr);
                 return -2;
            }
        }
/*
     4  pwhitcount            number of genes hit in pathway
     5  pwnohitcount          pathway    number of genes in the pathway
     6  inputlist                  total count of user genes (user input)
     7  pwuniverseminuslist    tinlistpathwaysuniverse    total number of unique genes in all pathways
a=hits
b=thispw-hits
c=list
d=u-list
*/
        a2 = localhitcnt;
        b2 = uptr->numfixedgenes - localhitcnt;
        c2 = incnt;
        d2 = real_universe_cnt - incnt;

        if (oneside)
        {
#if 0
            pv = exact22_oneside((uint32_t)a, (uint32_t)b, (uint32_t)c, (uint32_t)d,0); // last arg is debug flag
 // testing 
double pv2 = exact22(a,b,c,d);
fprintf(stderr,"after exact22_oneside pv1 = %20.17f , pv2 = %20.17f %d %d %d %d\n",pv,pv2,a,b,c,d); fflush(NULL);
// test22(a,b,c,d,pv);
#endif
#if 1
/*
int u,p,h,l;
u = real_universe_cnt;
p = uptr->numfixedgenes;
h = localhitcnt;
l = incnt;
int a2=hits
int b2=pathway - hits
int c2=list - hits
int d2=universe - pathway - list + hits
*/
a2 = localhitcnt;
b2 = uptr->numfixedgenes - localhitcnt;
c2 = incnt-localhitcnt;
d2 = real_universe_cnt - ( uptr->numfixedgenes ) - incnt + localhitcnt;
            pv2 = exact22_oneside((uint32_t)a2, (uint32_t)b2, (uint32_t)c2, (uint32_t)d2,0); // last arg is debug flag
// fprintf(stderr,"dbg pv = %f pv2=%f %d %d %d %d %d %d %d %d %s\n",pv,pv2,u,p,h,l,a2,b2,c2,d2,uptr->acc); 
            pv = pv2;
#endif
        }
        else
        {
            pv = exact22(a2,b2,c2,d2);
        }
        uptr->pval2 = *(pvals+i) = pv2;
        uptr->pval = *(pvals+i) = pv;
        uptr->a = a2;
        uptr->b = b2;
        uptr->c = c2;
        uptr->d = d2;
        // enrichment_score = ( ((double)localhitcnt / (double)uptr->numfixedgenes) - ((double)incnt/(double)real_universe_cnt) );
        enrichment_score = ( ((double)a2 / (double)(a2+b2)) - ((double)c2/(double)(c2+d2)) );
        enrichment_score *= 100.0;
        uptr->enrichment_score = enrichment_score;
    }
#if 0
for (i=0 ; i<num_used_paths ; i++)
{
 fprintf(stderr,"rpf checking 1 *(pvals+i) %d  %f\n",i,*(pvals+i) ); 
}
#endif

    fdrs = (double *)malloc((size_t)(sizeof (double)*num_used_paths)); 
    if (!fdrs) { free(pvals); /* clean up */ fprintf(stderr,"ERROR: no memory in do_pvals_and_bh() 2\n"); return -3; }

#if 0
for (i=0 ; i<num_used_paths ; i++)
{
fprintf(stderr,"rpf checking 2 *(pvals+i) %d  %f\n",i,*(pvals+i) ); 
}
#endif
// fprintf(stderr,"rpf before benjaminihochberg\n"); fflush(stderr); 
    benjaminihochberg(num_used_paths,pvals,fdrs);
// fprintf(stderr,"rpf after benjaminihochberg\n"); fflush(stderr); 
    for (i=0 ; i<num_used_paths ; i++)
    {
         usedpaths[i].fdr = *(fdrs+i);
// fprintf(stderr,"%d %f\n", i,usedpaths[i].fdr);
    }
    free(pvals);
    free(fdrs);
// fprintf(stderr,"in do_pvals_and_bh() end\n"); fflush(stderr); 
    return 0;
}


int do_just_bh(unsigned int incnt, struct used_path_type usedpaths[], unsigned int num_used_paths,unsigned int real_universe_cnt)
{
    double enrichment_score;
    unsigned int i;
    double *pvals;
    double *fdrs;
    unsigned int localhitcnt ;
    struct used_path_type *uptr; // used path pointer 
    int a,b,c,d;


// fprintf(stderr,"debug in do_just_bh() incnt=%d num_used_paths=%d real_universe_cnt=%d \n",incnt,num_used_paths,real_universe_cnt); fflush(stderr); 

    pvals = (double *)malloc((size_t)(sizeof (double)*(num_used_paths)));
    if (!pvals) { fprintf(stderr,"ERROR: no memory\n"); fflush(stderr); return -1; }
    for (i=0 ; i<num_used_paths;i++)
    {
        uptr = (usedpaths+i);
        localhitcnt = uptr->hitcnt;
        if (localhitcnt)
        {
            if (localhitcnt > (unsigned int)uptr->numfixedgenes) 
            {
                 fprintf(stderr,"ERROR: more hits than genes (localhitcnt %d > %d for %s)\n", 
                         localhitcnt , uptr->numfixedgenes,uptr->acc);
                 fprintf(stderr,"details: index=%d , hitcnt=%u\n", i,localhitcnt); 
                 fflush(stderr);
                 return -2;
            }
        }
/*
     4  pwhitcount            number of genes hit in pathway
     5  pwnohitcount          pathway    number of genes in the pathway
     6  inputlist                  total count of user genes (user input)
     7  pwuniverseminuslist    tinlistpathwaysuniverse    total number of unique genes in all pathways
a=hits
b=thispw-hits
c=list
d=u-list
*/
        a = localhitcnt;
        b = uptr->numfixedgenes - localhitcnt;
        c = incnt;
        d = real_universe_cnt - incnt;

  /* new way : */
       *(pvals+i) = uptr->pval;
        uptr->a = a;
        uptr->b = b;
        uptr->c = c;
        uptr->d = d;
        // enrichment_score = ( ((double)localhitcnt / (double)uptr->numfixedgenes) - ((double)incnt/(double)real_universe_cnt) );
        enrichment_score = ( ((double)a / (double)(a+b)) - ((double)c/(double)(c+d)) );
        enrichment_score *= 100.0;
        uptr->enrichment_score = enrichment_score;
// if (strcmp(uptr->acc,"ko00010") == 0) { fprintf(stderr,"adbg : do_just_bh ko00010 enrichment_score is %f from %d %d %d %d\n",enrichment_score,a,b,c,d); }
    }
    fdrs = (double *)malloc((size_t)(sizeof (double)*num_used_paths)); 
    if (!fdrs) { free(pvals); /* clean up */ fprintf(stderr,"ERROR: no memory in do_just_bh() 2\n"); return -3; }
#if 0
 // debug
 fprintf(stderr,"before benjaminihochberg %d %p %p\n",num_used_paths,pvals,fdrs); fflush(NULL); 
 for (i=0;i<num_used_paths;i++)
 {
 fprintf(stderr," pv %d %f\n",i,*(pvals+i)); 
 }
#endif
    benjaminihochberg(num_used_paths,pvals,fdrs);
    for (i=0 ; i<num_used_paths ; i++)
    {
         usedpaths[i].fdr = *(fdrs+i);
    }
    free(pvals);
    free(fdrs);
    return 0;
}




#if 0
#if L2PUSETHREADS
#define NUM_THREADS 4
int permute_range( unsigned int num_used_paths,struct used_path_type usedpaths[],unsigned int real_universe_cnt,unsigned int *real_universe, int incnt, int lo, int hi)
{
    int i;

    for (i=lo ; i<=hi ; i++)
    {
// fprintf(stderr, "%d from %d to %d of %d %s \n",i,lo,hi,num_used_paths,usedpaths[i].name); fflush(stderr); 
        permute(i,&usedpaths[i],real_universe_cnt,real_universe,incnt);
    }
    return 0;
}

struct thread_data_t {
    int tid;
    int myid;
    int lo;
    int hi;
    unsigned int ingenecnt;
    struct used_path_type *usedpaths;
    unsigned int num_used_paths;
    unsigned int real_universe_cnt;
    unsigned int *real_universe;
    unsigned int *sampling;
    float *pspace;
};

#if 0
/* thread function */
void *thr_func(void *arg) 
{
    struct thread_data_t *data = (struct thread_data_t *)arg;
    unsigned int *r;
   // double   p[NUM_PERMUTES];
    int i,j,k;
    double d;
    int tmphitcnt;
    unsigned int ui_uk,ui_ej;
    // float *pspace = data->pspace;
    unsigned int ingenecnt = data->ingenecnt;
    struct used_path_type *uptr = data->usedpaths;

    r = data->sampling;
    for (i=data->lo;i<data->hi;i++)
    {
        j = k = tmphitcnt = 0;
        while ((j<uptr->numfixedgenes) && (k < ingenecnt))
        {
            ui_ej = *(uptr->egids+j);
            ui_uk = *(r+k);
            if (ui_ej == ui_uk)
            {
                tmphitcnt++;
                k++;
                j++;
                continue;
            }
            else if (ui_ej < ui_uk) j++;
            else                    k++;
        }

        d = exact22((int)tmphitcnt,uptr->numfixedgenes,ingenecnt,       data->real_universe_cnt);
	//fix here
        // *(pspace + (i*NUM_PERMUTES)+data->used_path_index; = d;
    }
#if 0
    qsort(p,NUM_PERMUTES,sizeof(double),cmp_double);
    d = 1.0;
    for (j=0;j<NUM_PERMUTES;j++) 
    {
        if (p[j] > uptr->pval) { d = p[j]; break; }
    }
    d = (double)j/(double)NUM_PERMUTES;
    uptr->fdr2 = d;
#endif
  pthread_exit(NULL);
}
#endif
#endif


#if 0
int oldthreaded_permutes(unsigned int ingenecnt, struct used_path_type usedpaths[], unsigned int num_used_paths,unsigned int real_universe_cnt,unsigned int *real_universe)
{
    pthread_t thr[NUM_THREADS];
    unsigned int r[40002];
    struct thread_data_t thr_data[NUM_THREADS];
    int i,j,k;
    int rc;
    unsigned int *pspace;

    pspace = malloc(sizeof(float)*NUM_PERMUTES*num_used_paths); // no dire need to memset
    memcpy(r,real_universe,real_universe_cnt*sizeof(unsigned int));
fprintf(stderr,"in oldthreaded_permutes()\n"); fflush(stderr); 
    for (i=0;i<NUM_PERMUTES;i++)
    {
        subsamp(r, real_universe_cnt,ingenecnt); // inlined if optimization on 
        quickSort(r,ingenecnt); // inlined if optimization on 
        k = NUM_PERMUTES/NUM_THREADS;
        for (j = 0; j < NUM_THREADS; ++j) 
        {     /* create threads */
            thr_data[j].lo = j*k;
            thr_data[j].hi = (j*k)+k-1;
            if (j == (NUM_THREADS-1))
                thr_data[j].hi = num_used_paths-1; // note "-1", 
fprintf(stderr,"%d of %d [ %d .. %d ] of %d size= %d\n",j,NUM_THREADS,thr_data[j].lo,thr_data[j].hi,num_used_paths,thr_data[j].hi-thr_data[j].lo); 
            thr_data[j].myid = j;
            thr_data[j].ingenecnt = ingenecnt;
            thr_data[j].sampling = r;
            // FIX thr_data[j].usedpathsindex = usedpaths;
            thr_data[j].usedpaths = usedpaths;
            thr_data[j].num_used_paths = num_used_paths;
fprintf(stderr,"creating thread %d\n",j); fflush(stderr); 
            if ((rc = pthread_create(&thr[j], NULL, thr_func, &thr_data[j]))) 
            {
              fprintf(stderr, "error: pthread_create, rc: %d\n", rc);
              return EXIT_FAILURE;
            }
        }
      /* block until all threads complete */
        for (j = 0; j < NUM_THREADS; ++j) 
        {
fprintf(stderr,"blocking %d \n",j);
            pthread_join(thr[j], NULL);
        }
fprintf(stderr,"after blocking  \n");
    }
  return EXIT_SUCCESS;
}
#endif

#endif



struct tree_type
{
    unsigned int val; // entrez gene id
    struct tree_type *left;
    struct tree_type *right;
};

void free_tree(struct tree_type *node)
{
    if (node==(void *)0) return;
    if (node->left != (void *)0)  free_tree(node->left);
    if (node->right != (void *)0) free_tree(node->right);
    free(node); 
    return;
}

void put_tree_to_array(struct tree_type *node, unsigned int puthere[], int *index_ptr)
{
    if (node != (void *)0)
    {
        put_tree_to_array(node->left,puthere,index_ptr);
        puthere[*index_ptr] = node->val;
        *index_ptr = *index_ptr + 1;
        put_tree_to_array(node->right,puthere,index_ptr);
    }
    return;
}
void free_tree_with_count(struct tree_with_count *node)
{
    if (node)
    {
        if (node->all_gene_paths) free(node->all_gene_paths);
        if (node->left) free_tree_with_count(node->left); // recursive calls
        if (node->right)free_tree_with_count(node->right);
        free(node);
    }
    return;
}
int put_to_tree(unsigned int this_gene, struct tree_with_count *head) // tree below head keeps growing; we will keep using it
{
    // has to:
    // 1. create a node (after the first, head) for each new gene (fill in gene number to node)
    //      --in left-right order
    // 2. count the number of unique genes (one for each new node, n.b. first node created above)
    // 3. add how many instances there are of each gene to ->count
     struct tree_with_count *z, *n;
     z = head;
     //printf("this_gene = %d\n", this_gene);
    while (1) // follow the tree down to the bottom switching right or left by gene number. New node if needed. Then outer loop takes to next gene
     {
          if (z->val == this_gene)
          {
              z->count += 1;
              //if (z->count == 1) ++real_genect; // shouldn't happen, we've seen this node before // this wastes time...
              //head->all_gene_paths = malloc(sizeof(struct used_path_type));// not yet
              //head->all_gene_paths = *(augreverse+ui);
              return(1); //gene has been placed, return
          }
          else if (this_gene < z->val)
          {
              if (z->left == (void *)0)
              {
                  n = z->left = malloc(sizeof(struct tree_with_count));
                  n->val=this_gene;
                  //n->all_gene_paths = &this_path; ? have to run genes and malloc to # repeats?
                  ++real_genect; // OK, it's a new node
                  n->count = 1;
                  n->deg = 0;
                  n->left=(void *)0;
                  n->right=(void *)0;
                  n->all_gene_paths=(void *)0;
                  n->pathindex = 0;
                  return(1);//gene has been placed, return
              }
              z = z->left; // and back to top of while
          }
          else  /// i.e.  (z->val > this_gene)
          {
              if (z->right == (void *)0)
              {
                  n = z->right = malloc(sizeof(struct tree_with_count));
                  n->val=this_gene;
                  ++real_genect;
                  n->count = 1;
                  n->deg = 0;
                  n->left=(void *)0;
                  n->right=(void *)0;
                  n->all_gene_paths=(void *)0;
                  n->pathindex = 0;
                  return(1);//gene has been placed, return
              }
              z = z->right; // and back to top of while
          }
     }
}



// RICH: Beginning of my code
void malloc_pathpointers(struct tree_with_count *node) // counts aligned with universe (real_universe)
{
    // int i;
// rpf     static int degindex = 0;
    if (node)
    {
        malloc_pathpointers(node->left);
        node->all_gene_paths = malloc(sizeof(struct used_path_type *)*(node->count));
               // node->count is number pathways this gene appears in.
        if (!node->all_gene_paths)
        {
            fprintf(stderr,"ERROR: not enough memory in malloc_pathpointers()\n"); fflush(stderr);
        }
        malloc_pathpointers(node->right);
    }
    return;
}
int path_pointers_to_tree(unsigned int this_gene, struct used_path_type *this_path, struct tree_with_count *head)
// following put_to_tree, but tree is already there
// I think right, have to find the unsorted genes here, this is guide for filling the pointer
// could have been done while the tree was created, but don't have space for the pointers malloc'd yet
{
    // all nodes exist,
    // has to:
     struct tree_with_count *z;
     // struct tree_with_count *n;
     z = head;
     while (1) // follow the tree down to the bottom switching right or left by gene number. Tree created already, fails if new node needed.
     {
          if (z->val == this_gene)
          {
              z->all_gene_paths[z->pathindex] = this_path; // all_gene_paths is malloc pointers to pathways[count] (aka: used_pathways )
              ++(z->pathindex);
              if (z->pathindex > z->count) // note this catches trouble after it's started (if z->pathindex > z->count, overwrote z->all_gene_paths above).  Could put above but wastes time.
              {
                  printf("bad pathway count for gene %ds, exiting\n", this_gene);
                  exit(0);
              }
              return(1);
          }
          else if (this_gene < z->val)
          {
              if (z->left == (void *)0)
              {
                  printf("gene %d not found while sorting pathpointers, exiting\n", this_gene);
                  exit(0);
              }
              z = z->left; // and back to top of while
          }
          else  // i.e.  (z->val > this_gene)
          {
            if (z->right == (void *)0)
            {
                  printf("gene %d not found while sorting pathpointers, exiting\n", this_gene);
            }
            z = z->right; // and back to top of while
          }
     }
}

static inline void put_count_tree_to_array(struct tree_with_count *node, unsigned int universe[], unsigned int *index_ptr, unsigned int *aug_ptr) // counts aligned with universe (real_universe)
{
    // no trap for going off the tree
    unsigned int ui;
    struct used_path_type *thispath;
    if (node != (void *)0) //  trap against going off the tree, but are we catching genes with no match?
    {
        //printf("left\t");
        put_count_tree_to_array(node->left,universe,index_ptr,aug_ptr);
        universe[*index_ptr] = node->val;
        // printf("tst tree %d\n",node->val);
        gene_path_cts[*index_ptr] = node->count;
        (*index_ptr) = (*index_ptr) + 1;
#if NELSON_TEST
        // debug
        // fprintf(gene_path_countsfile, "%d\t%d\n", node->val, node->count);
#endif
        for (ui = 0; ui < node->count; ui++)
        {
            if (*(aug_ptr) < aug_gene_ct) 
            {
// if not extra malloc, this messes up.  rpf fix
// fprintf(stderr,"*aug_ptr = %d  *index_ptr = %d\n",*aug_ptr,*index_ptr);  fflush(stderr);
                auguniverse[*aug_ptr] = node->val; // to randomly select genes for permutation test, weighted by pathway count
                aug_treepointers[*aug_ptr] = treepointers[*index_ptr] = node;   // multiple pointers to same node,
            }
            else
            {
                fprintf(stderr,"ERROR: *aug_ptr = %d , should be less than %d ui=%u node->count=%d \n",*aug_ptr,aug_gene_ct,ui,node->count); 
            }
            //printf("%d %d\n", *aug_ptr, node->val);
            //augreverse[*aug_ptr] = node; // replaced by
            // gp counts to the pathway
            thispath = node->all_gene_paths[ui];
            thispath->pathcountsum += node->count; // this sums the gene gp cts for all hits in the pathway
            *aug_ptr = *aug_ptr + 1;
        }
        //printf("right\t");
        put_count_tree_to_array(node->right,universe,index_ptr,aug_ptr);
    }
    //printf("up\t");
    return;
}


static inline void put_deg_tree_to_array(struct tree_with_count *node, unsigned int *index_ptr) // counts aligned with universe (real_universe)
{
    // no trap for going off the tree
    if (node != (void *)0) //  trap against going off the tree, but are we catching genes with no match?
    {
        put_deg_tree_to_array(node->left,index_ptr);
        if (node->deg)
        {
            deg_path_cts[*index_ptr] = node->count;
            aug_deg_count += node->count;
#if 0
struct tree_with_count
{
    unsigned int val; // entrez gene id
    unsigned int count; // times
    unsigned int deg; // 1 on deglist, 0 not on
    struct tree_with_count *left;
    struct tree_with_count *right;
    struct used_path_type **all_gene_paths; // all gene paths is an array of pointers
    unsigned int pathindex; // which array member gets the pointer to pathway?
};
#endif
            ingenes[*index_ptr] = node->val;
            ++*index_ptr;
        }
        put_deg_tree_to_array(node->right,index_ptr);
    }
    return;
}

 // gets called for every user in gene in their input list
int degs_to_tree(unsigned int this_gene, struct tree_with_count *head) // tree below head keeps growing; we will keep using it
{
    unsigned int i;
    struct tree_with_count *z, *node;   // fields: val:gene_id , count:times, deg:flag for in deglist, all_gene_paths: pointers to all pathways
#if 0
struct tree_with_count
{
    unsigned int val; // entrez gene id
    unsigned int count; // times
    unsigned int deg; // 1 on deglist, 0 not on
    struct tree_with_count *left;
    struct tree_with_count *right;
    struct used_path_type **all_gene_paths; // all gene paths is an array of pointers
    unsigned int pathindex; // which array member gets the pointer to pathway?
};
#endif

    struct used_path_type *thispath;
     z = head;
     while (1) // follow the tree down to the bottom switching right or left by gene number. If hit void, gene not in universe. Then outer loop takes to next gene
     {
          if (z->val == this_gene)
          {
              z->deg = 1;
              node = z;
              ++deg_count; // rpf: don't really need this, we know the deg_count.
              for (i=0; i<node->count; i++)
              {
                  thispath = node->all_gene_paths[i];
                  thispath->pathhits_gpsum += node->count;  // this sums the gene gp cts for all hits in the pathway
// rpf note: so this is the number of pathways that all deg genes in the pathway hit
              }
              return(1);
          }
          else if (this_gene < z->val)
          {
              if (z->left == (void *)0)
              {   // fprintf(stderr, "deg %d not in universe\n", this_gene);
                  return(0); // deg not in universe
              }
              z = z->left; // and back to top of while
          }
          else       /// i.e. (z->val > this_gene)
          {
              if (z->right == (void *)0)
              {    // fprintf(stderr, "deg %d not in universe\n", this_gene);
                  return(0); // deg not in universe
              }
              z = z->right; // and back to top of while
          } // but if coded wrong could go on forever, no trap; but ok, keeps going down till a null
     }
}


void chooseAugRandList(unsigned int n, unsigned int m, unsigned int *randlist) // n is size of augmented gene universe, m number to choose
{
    unsigned int i,j;
    unsigned int randint;
    unsigned int randgene;
    int already_got;
    // double scale = (double) (n-1)/(RAND_MAX);

// fprintf(stderr,"in chooseAugRandList() auggenecnt n=%d,  ingenecnt m=%d\n",n,m); fflush(stderr); 
    for (i=0;i<m;)
    {
#ifdef L2P_USING_R
        double randrslt;
        randrslt = unif_rand();
// wrong?        randrslt = (int)(unif_rand() * (double)i)+1;   // unif_rand() appears to return between 0.0 and 1.0
        // randint = (int)(randrslt*m);
        randint = (int)(randrslt*(double)n);
// fprintf(stderr,"in chooseAugRandList()3.1, m=%d randint=%d unif_rand randrslt=%f\n",m,randint,randrslt); fflush(stderr); 
#else
        randint = rand() % n;
#endif
        // char *zz;
        randgene = auguniverse[randint];
        // zz = egid2hugo(randint);
// fprintf(stderr,"in chooseAugRandList() 3.2 loop %d to %d i=%d randint=%d scale=%f zz=%s\n",i,m,i,randint,scale,zz); fflush(stderr); 
        for (already_got = j=0;j<i;j++)
        {
            if (randgene == randlist[j]) { already_got = 1; break; }
        }
        // zz = egid2hugo(randgene);
// fprintf(stderr,"in chooseAugRandList() 3.3 end loop %d to %d, alreadygot=%d, randgene=%u at j=%d (%s)\n",i,m,already_got,randgene,j,zz); fflush(stderr); 
        if (already_got) continue;
        randlist[i++] = randint;
            //fprintf(stdout, "i, j, chosen index, gene %f %d %d %d %d \n",randrslt, i, j, randint, randgene);
// fprintf(stderr,"in chooseAugRandList() 3.4 end loop %d to %d\n",i,m); fflush(stderr); 
    }
// fprintf(stderr,"in chooseAugRandList()end \n"); fflush(stderr); 
}


void chooseRandList(unsigned int n, unsigned int m,  unsigned int *randlist) // n is size of augmented gene universe, m number to choose
{
    unsigned int i;
    double randrslt;
    double scale = (double) (n-1)/(RAND_MAX);
    // unsigned int j;
    // unsigned int randint;
    // unsigned int randgene;
    //double scale = (double) n/RAND_MAX;

    for (i=0 ; i<m ; i++)
    {
#ifdef L2P_USING_R
        randrslt = (int)(unif_rand());
#else
        randrslt = rand();
#endif
        randlist[i] = randrslt*scale+1;
    }
}/**/


#if 0
 // *** rpf does not seem to be called 
int gpcount_correction(struct used_path_type usedpaths[], unsigned int num_used_paths) // ,int oneside)
{
    unsigned int i;
    double mean_gpct;
    struct used_path_type *uptr; // used path pointer
    for (i=0 ; i<num_used_paths;i++)
    {
        uptr = (usedpaths+i);
        mean_gpct =uptr->pathhits_gpsum/uptr->hitcnt;
        //if (uptr->hitcnt == 0) uptr->gpcc_p = 99;
        //else
            uptr->gpcc_p = pow(mean_gpct/mean_u_gpcount, log(uptr->hitcnt))*uptr->pval;
    }
    return 0;
}
#endif

int tablecalc(struct used_path_type usedpaths[], unsigned int num_used_paths, unsigned int real_universe_cnt) // ,int oneside)
{
    double  OR, pv, aug_pv;
    double aug_scale;
    unsigned int i;
    double *pvals;
    double *pvals2;
    double *fdrs;
    double *fdrs2;
    unsigned int localhitcnt;
    struct used_path_type *uptr; // used path pointer
    int a,b,c,d, A,B,C,D;
    int A_scaled, B_scaled, C_scaled, D_scaled;
    // double enrichment_score;
    
// fprintf(stderr,"gpccdbg: in tablecalc() ingenecnt=%d num_used_paths=%d real_universe_cnt=%d \n",ingenecnt,num_used_paths,real_universe_cnt); fflush(stderr);
    
    pvals = (double *)malloc((size_t)(sizeof (double)*(num_used_paths)));  // all in pathway struct, don't need
    if (!pvals) { fprintf(stderr,"ERROR: no memory\n"); fflush(stderr); return -1; }
    memset(pvals,0,sizeof (double)*(num_used_paths));
    pvals2 = (double *)malloc((size_t)(sizeof (double)*(num_used_paths)));  // all in pathway struct, don't need
    if (!pvals2) { fprintf(stderr,"ERROR: no memory\n"); fflush(stderr); return -1; }
    memset(pvals2,0,sizeof (double)*(num_used_paths));
    for (i=0 ; i<num_used_paths;i++)
    {
        uptr = (usedpaths+i);
        localhitcnt = uptr->hitcnt;
        //if (!localhitcnt) continue; //  Let's only look at overexpression
        if (localhitcnt > (unsigned int)uptr->numfixedgenes)
        {
             fprintf(stderr,"ERROR: more hits than genes (localhitcnt %d > %d for %s)\n",
                     localhitcnt , uptr->numfixedgenes,uptr->acc);
             fprintf(stderr,"details: index=%d , hitcnt=%u\n", i,localhitcnt);
             fflush(stderr);
             return -2;
        }
/*
     4  pwhitcount            number of genes hit in pathway
     5  pwnohitcount          pathway    number of genes in the pathway
     6  inputlist                  total count of user genes (user input)
     7  pwuniverseminuslist    tinlistpathwaysuniverse    total number of unique genes in all pathways
    a=hits
    b=thispw-hits
    c=list
    d=u-list
*/
        /*a = localhitcnt;
        b = uptr->numfixedgenes - localhitcnt;
        c = ingenecnt;
        d = real_universe_cnt - ingenecnt;*/
        
        d = localhitcnt;
        c = ingenecnt - localhitcnt;
        b = uptr->numfixedgenes - localhitcnt;
        a = real_universe_cnt - uptr->numfixedgenes - ingenecnt + localhitcnt;

        /* new way : */
        uptr->a = a;
        uptr->b = b;
        uptr->c = c;
        uptr->d = d;
        // enrichment_score = ( ((double)localhitcnt / (double)uptr->numfixedgenes) - ((double)ingenecnt/(double)real_universe_cnt) );
#if 0
  This is not right
        enrichment_score = ( ((double)a/(double)b) - ((double)c/(double)d) );
        uptr->enrichment_score = enrichment_score;
if (strcmp(uptr->acc,"ko00010") == 0) { fprintf(stderr,"adbg : tablecalc ko00010 enrichment_score is %f from %d %d %d %d\n",enrichment_score,a,b,c,d); }
#else
   // rpf 
        uptr->enrichment_score = ( ((double)a / (double)(a+b)) - ((double)c/(double)(c+d)) ) * 100.0;
#endif
        if ((c*b) != 0)
        {
            OR = ((double)a*(double)d)/((double)c*(double)b);
            uptr->OR = OR;
#if NELSON_TEST
            if (OR > 1000)
            {
                fprintf(stderr,"Note:big OR %f, path %s\n", OR, uptr->acc);
            }
#endif
        }
        else uptr->OR = OR = 99999;
        pv = exact22(a,b,c,d);
        //pv = exact22_oneside((uint32_t)a, (uint32_t)b, (uint32_t)c, (uint32_t)d,0); // last arg is debug flag
        /*if (!strncmp("basal lamina", uptr->name, 30))
        {
            printf("%s", uptr->name);
        }*/
        uptr->pval = *(pvals+i) = pv;
#if 0
        if ((pv < SIG_P) & (OR > 1)) ++n_sig_paths;
if (strcmp(uptr->acc,CHECKACC) == 0) { fprintf(stderr,"Got M40028 in GPCC, pathhits_gpsum=%f , pathcountsum=%u ,hits=%d, aug_gene_ct=%u aug_scale=%f\n",uptr->pathhits_gpsum,uptr->pathcountsum,uptr->hitcnt,aug_gene_ct,aug_scale); fflush(stderr); } 
#endif
//  VERY IMPORTANT CODE,  HERE'S WHERE IT HAPPENS ...
        D = uptr->pathhits_gpsum;     // wrong? number of pathways any and all genes in this pathway hit
        C = aug_deg_count - D;        // or d?   aug_deg_count is number of times the user input genes hit in any and all pathways.  
        B = uptr->pathcountsum - D;   // pathcountsum is the sum of the gene gp cts for all hits in the pathway
        A = aug_gene_ct - B - C - D;  //  aug_gene_ct is sum of all numfixedgenes in usedpaths[]
        aug_scale = (float)real_universe_cnt/((float)(aug_gene_ct + d)); // real_universe_cnt is count of unique genes in universe
        uptr->A_scaled = A_scaled = ceil(A*aug_scale);
        uptr->B_scaled = B_scaled = ceil(B*aug_scale);
        uptr->C_scaled = C_scaled = ceil(C*aug_scale);
        uptr->D_scaled = D_scaled = d; // not scaled!
        aug_pv = exact22(A_scaled,B_scaled,C_scaled,D_scaled);
        uptr->gpcc_p = aug_pv;        // heres the pvalue for GPCC ***
#define CHECKACC "R-HSA-5662702"
 if (strcmp(uptr->acc,CHECKACC) == 0) 
 { 
 fprintf(stderr,"Got %s in tablecalc, pathhits_gpsum=%f , pathcountsum=%u ,hit=%d aug_gene_ct=%u aug_scale=%f aug_deg_count=%u\n",
          CHECKACC,uptr->pathhits_gpsum,uptr->pathcountsum,uptr->hitcnt,aug_gene_ct,aug_scale,aug_deg_count); 
 fprintf(stderr,"Got %s in tablecalc aug_pv=%20.18f scaled= a.%d b.%d c.%d d.%d C=%d\n",CHECKACC,aug_pv,A_scaled,B_scaled,C_scaled,D_scaled,C); 
 fflush(stderr); 
 } 
        *(pvals2+i) = aug_pv;
        if (C_scaled*B_scaled != 0)
        {
            uptr->gpcc_OR = ((double)A_scaled*(double)D_scaled)/((double)C_scaled*(double)B_scaled);
        }
        else uptr->gpcc_OR = 99999;

    }
#if 0
        sig_path_ptr = malloc((sizeof(struct used_path_type *)*n_sig_paths));
#endif
    // HERE WE WOULD WANT TO SORT THE SIGNIFICANT P VALUES, FOR SPEED LATER
    // struct used_path_type *tstptr;
    //  I HAVEN'T LOOKED AT THE FDR CODE
    fdrs = (double *)malloc((size_t)(sizeof (double)*num_used_paths));
    if (!fdrs) { free(pvals); /* clean up */ fprintf(stderr,"ERROR: no memory in tablecalc() 2\n"); return -3; }
    fdrs2 = (double *)malloc((size_t)(sizeof (double)*num_used_paths));
    if (!fdrs2) { free(pvals); /* clean up */ fprintf(stderr,"ERROR: no memory in tablecalc() 2\n"); return -3; }
#if 0
     // debug
    fprintf(stderr,"before benjaminihochberg %d %p %p\n",num_used_paths,pvals,fdrs); fflush(NULL);
    for (i=0;i<num_used_paths;i++)
    {
    fprintf(stderr," pv %d %f\n",i,*(pvals+i));
    }
#endif
#if 0
 // rpf fix  : need to get gpcc pvals, not FE pvals into pvals 
        //j *(pvals+i) = uptr->pval;
    for (i=0;i<num_used_paths;i++)
    {
          *(pvals+i) = usedpaths[i].gpcc_p ;
    }
#endif
// fprintf(stderr,"rpf here before benjaminihochberg(%d,%p,%p), num_used_paths=%d\n",num_used_paths,pvals,fdrs,num_used_paths); fflush(NULL);
    benjaminihochberg(num_used_paths,pvals,fdrs);
    benjaminihochberg(num_used_paths,pvals2,fdrs2);
// fprintf(stderr,"rpf here after benjaminihochberg(%d,%p,%p), num_used_paths=%d\n",num_used_paths,pvals,fdrs,num_used_paths); fflush(NULL);
    for (i=0 ; i<num_used_paths ; i++)
    {
         usedpaths[i].fdr = *(fdrs+i);
         usedpaths[i].gpcc_fdr = *(fdrs2+i);
// fprintf(stderr,"gpccdbg gpcc_fdr = %f\n", usedpaths[i].gpcc_fdr);
    }
    free(pvals);
    free(fdrs);
    free(pvals2);
    free(fdrs2);
    return 0;
}

unsigned int GPCC(struct used_path_type usedpaths[], unsigned int num_used_paths, unsigned int real_universe_cnt, unsigned int *real_universe)
{
/*
     Create a vector of gene names, with each name appearing as often as the number of pathways the gene appears in.
     From this list randomly select as many genes as on your DEG list.
     Replace repeats, need a list of unique genes
     ii.     Create the contingency table for each pathway of interest (list of top pathways)
     iii.    Obtain and record its p value
     iv.     For each pathway of interest keep count of how many of the 200 permutations produced more significant p values than the real data
     
     Repeat a) etc. 200 times
     The counts/200 are approximate p values for the pathways.
*/
    unsigned int i,I,j,k, ll;
    int utilctr;
    unsigned int this_gene;
    unsigned int ui = 0;
    unsigned int ui_ej = 0;
    unsigned int *randgenesindex = (unsigned int *)0;
    unsigned int *pool = (unsigned int *)0;
    struct tree_with_count *head = (struct tree_with_count *)0;
    struct used_path_type **augreverse = (struct used_path_type **)0; // parallel to augmented universe ("auguniverse") , a "reverse lookup" mechanism
    struct used_path_type  *this_path = (struct used_path_type *)0;
    unsigned int *index_ptr = (unsigned int *)0;
    unsigned int *aug_ptr = (unsigned int *)0;
    struct used_path_type *uptr = (struct used_path_type *)0;    // used path pointer
    double *pspace = (double *)0;   // array of pvalues[num_used_paths*NUM_TEST_PERMUTES]
#if NELSON_TEST
    FILE *deg_path_count_file = (FILE *)0;
#endif
    FILE *pathgenes_fp = (FILE *)0;
#if 0
    unsigned int d_up = 0;
    usnigned int d_down = 0;
    unsigned int d_equal = 0;
    unsigned int a, rand_d, rand_c, rand_b;
    int var_d;
    double randOR;
    unsigned int sig_rand_p_ct = 0;
    double pv;
#endif


// fprintf(stderr,"rpf in GPCC(), real_universe_cnt=%u\n",real_universe_cnt); fflush(NULL);
#if 0
for(i=0;i<real_universe_cnt;i++)
 {
    char *zz;
fprintf(stderr,"gene %u = %u ",i,*(real_universe+i));
    zz = egid2hugo(*(real_universe+i));
    if (zz) fprintf(stderr,"%s ",zz);
    else    fprintf(stderr,"NA ");
fprintf(stderr,"\n");
 }
#endif


    //int *pool = (unsigned int *)0;
    // Set up augmented universe, count_tree, and pointers to pathways
#if NELSON_TEST
//    gene_path_countsfile = fopen("gene_pathway_counts.txt", "w");
    deg_path_count_file = fopen("deg_pathway_counts.txt", "w");
#endif
    //path_gene_countsfile = fopen("pathway_gene counts.txt", "w");
    if (!usedpaths)
    {
        fprintf(stderr,"ERROR: null used_path_type passed to GPCC()\n"); fflush(stderr);
        return (unsigned int)0 ;
    }
    // get count of all genes with repititions in pathways
#if 0
        pathgenes_fp = fopen("paths_genes_names.txt", "w");
#endif

// fprintf(stderr,"rpf in GPCC(), real_universe_cnt=%u 2\n",real_universe_cnt); fflush(NULL);
    for (i=aug_gene_ct=0 ; i<num_used_paths ; i++)
    {         // get number of genes in all pathways, count duplicated gene names
        uptr = (usedpaths+i);
        aug_gene_ct += uptr->numfixedgenes;
        //fprintf(path_gene_countsfile, "%s \t %d\n", uptr->name, uptr->numfixedgenes);
        /*for (j=0;j<uptr->numfixedgenes;j++)
        {
            fprintf(pathgenes_fp, "%s\t %s\t %d\n", uptr->name, uptr->acc, uptr->egids[j]);
        }*/
    }
// fprintf(stderr,"rpf debug aug_gene_ct = %d\n",aug_gene_ct); fflush(stderr); 
    if (pathgenes_fp) fflush(pathgenes_fp);
#if 1
    // fix fix fix  rpf
    auguniverse = malloc((aug_gene_ct+10) * sizeof(unsigned int)); // seems to be short by one uint , so just add it . what the heck? this is a bug. fix this.
    mean_u_gpcount = (double) aug_gene_ct/(double) real_universe_cnt;
#else
    auguniverse = malloc(aug_gene_ct * sizeof(unsigned int));
    mean_u_gpcount = (float) aug_gene_ct/(float) real_universe_cnt;
#endif
    if (!auguniverse)
    {
        fprintf(stderr,"ERROR: not enough memory. in GPCC()\n"); fflush(stderr);
        return (unsigned int)0 ;
    }
// fprintf(stderr,"rpf debug aug_gene_ct = %d mean_u_gpcount=%f\n",aug_gene_ct,mean_u_gpcount); fflush(stderr); 
    // for fast permute:  I need pointers from each gene (repeated) to all the pathways containing the gene
    // need to be sorted in gene number order
    // NOTE MANY VARIANTS ON PUTTING TO TREE AND GETTING FROM TREE FOLLOW HERE
// fprintf(stderr,"rpf in GPCC(), real_universe_cnt=%u 3\n",real_universe_cnt); fflush(NULL);
    augreverse = malloc(sizeof(struct used_path_type *)*(aug_gene_ct)); // temp, sorted contents go to ordrd_pathpointer
    for (i=k=0;i<num_used_paths;i++)
    {
        uptr = (usedpaths+i);
        for (j=0 ; j < uptr->numfixedgenes ; j++)
        {
            *(auguniverse+k) = *(uptr->egids + j);
            *(augreverse+k) = uptr; // since unlabelled, essential that this array of pointers has the same order as auguniverse
            k++;
            if (k > aug_gene_ct) // rpf: was just > -- 
            {
                fprintf(stderr,"ERROR: in GPCC() i=%d k=%d >= aug_gene_ct=%d . Overflow number of genes.\n",i,k,aug_gene_ct);  fflush(stderr);
                exit(0);
            }
        }
    }
// fprintf(stderr,"rpf in GPCC(), real_universe_cnt=%u 4\n",real_universe_cnt); fflush(NULL);
    head = malloc(sizeof(struct tree_with_count));
    head->val = *(auguniverse); // entrez gene id
    real_genect = 1;
    head->count = 1; // each node starts with count 1
    head->deg = 0;
    head->left = head->right = (void *)0;
    // run first in get used universe, or separate function
    // this gets counts
    // malloc to counts
    head->all_gene_paths = (void *)0;
    head->pathindex = 0;
    for (ui=1 ; ui<aug_gene_ct ; ui++)
    {
        this_gene = *(auguniverse+ui);
        this_path = *(augreverse+ui);    // this_path appears to not be used
        put_to_tree(this_gene, head); // for each gene, count the number of time it is in auguniverse
    }
    
// fprintf(stderr,"rpf in GPCC(), real_universe_cnt=%u 5\n",real_universe_cnt); fflush(NULL);
    index_ptr = malloc(sizeof(unsigned int));
    memset(index_ptr,0,sizeof(unsigned int));
    aug_ptr = malloc(sizeof(unsigned int)*2); // wtf?
    memset(aug_ptr,0,sizeof(unsigned int)*2);
    aug_treepointers = malloc(sizeof(struct tree_with_count **)*(aug_gene_ct+500)); // valgrind complains  fix bug
    memset(aug_treepointers ,0,sizeof(struct tree_with_count **)*(aug_gene_ct+500));
    treepointers = malloc(sizeof(struct tree_with_count **)*(real_universe_cnt+500));
    memset(treepointers ,0,sizeof(struct tree_with_count **)*(real_universe_cnt+500));
    // rpf treepointers = malloc(sizeof(struct tree_with_count **)*real_universe_cnt);
    gene_path_cts = malloc(sizeof(unsigned int)*(real_genect+1)); // valgrind complains so add 1
    memset(gene_path_cts ,0,sizeof(unsigned int)*(real_genect+1));

// fprintf(stderr,"rpf in GPCC(), real_universe_cnt=%u 6\n",real_universe_cnt); fflush(NULL);
    malloc_pathpointers(head);
     // fprintf(stderr,"test after malloc_pathpointers\n"); fflush(NULL);
    utilctr = 0;
    for (ui=0 ; ui<aug_gene_ct ; ui++)
    {
         this_gene = *(auguniverse+ui);
         this_path = *(augreverse+ui);
         path_pointers_to_tree(this_gene, this_path, head); // recall path_pointers are needed for permutation test
         if(this_gene == 7040) ++utilctr;   // rpf - what is this ?  APPEARS TO BE GN DEBUG
    }
    *index_ptr = 0; // fcn passes down pointer
    *aug_ptr = 0; // fcn passes down pointer
    put_count_tree_to_array(head,real_universe,index_ptr, aug_ptr); // n.b. auguniverse is now a global variable -- sorted by this recursive function
    deg_count = 0;
    for(i=0 ; i<ingenecnt ; i++)
    {
        this_gene = ingenes[i];
        degs_to_tree(this_gene, head);
    }
// fprintf(stderr,"rpf in GPCC(), real_universe_cnt=%u 7\n",real_universe_cnt); fflush(NULL);
    ingenecnt = deg_count; // should just use deg_count everywhere
    deg_path_cts = malloc(sizeof(unsigned int)*deg_count);
    memset(deg_path_cts,0,sizeof(unsigned int)*deg_count);
    *index_ptr = 0;
// fprintf(stderr,"in GPCC 2.0 head=%p index_ptr=%p\n",head,index_ptr); fflush(NULL);
    aug_deg_count = 0; // redundant
    put_deg_tree_to_array(head,index_ptr); // not getting deg count here, should I?;  note now could integrate with put_count_tree_to_array
    // no longer used I think
    for (i=0; i<deg_count; i++) deg_count_sum += deg_path_cts[i];
/*
fprintf(stderr,"before mean_deg_gpcount = deg_count_sum/deg_count;\n"); fflush(NULL);
fprintf(stderr,"ingenecnt = %d ",ingenecnt); fflush(NULL);
fprintf(stderr,"deg_count_sum = %d ",deg_count_sum); fflush(NULL);
fprintf(stderr,"deg_count = %d\n",deg_count); fflush(NULL);
*/

// fprintf(stderr,"rpf in GPCC(), real_universe_cnt=%u 8\n",real_universe_cnt); fflush(NULL);
    if (deg_count != 0) mean_deg_gpcount = deg_count_sum/deg_count;
    else                mean_deg_gpcount = 0;
// fprintf(stderr,"after mean_deg_gpcount = deg_count_sum/deg_count = %f\n",mean_deg_gpcount); fflush(NULL);
    for (i=j=0 ; i<num_used_paths ; i++)
    {
        uptr = (usedpaths+i);
        if (!uptr->egids) continue;
        uptr->aughitcnt = j = k = ll = 0;
        while ((j<uptr->numfixedgenes) && (k < deg_count))
        {
            ui_ej = *(uptr->egids+j);
            ui = ingenes[k];
            if (ui_ej == ui)
            {
                *((uptr->genehits) + (ll++)) = ui; // remember, because need to print out later
                uptr->aughitcnt += deg_path_cts[k]; // WRONG?  IF SO, WHY?
                k++;
                j++;
                // aug hit count for aug contingency table
                continue;
            }
            else if (ui_ej < ui) j++;
            else                 k++;
        }
        uptr->hitcnt = ll;
    }

// fprintf(stderr,"in GPCC 2.1 deg_count=%d \n",deg_count); fflush(NULL);
#if NELSON_TEST
    for (i=0;i<deg_count;i++) 
    {
         if (deg_path_count_file) fprintf(deg_path_count_file, "%d\n  ", deg_path_cts[i]);
    }
    fclose(deg_path_count_file); 
    deg_path_count_file = (FILE *)0;
#endif
    // now I will select randomly from auguniverse
    // this produces a set of gene numbers
    // For each gene number, in I have an array of pointers to the pathways--
    // --in the pathway count list
    // rows will correspond to pathways in ptr_from_unique. Need a consistent list of pathways (1-n); the row refers to this list
// fprintf(stderr,"in GPCC 2.2 deg_count=%d \n",deg_count); fflush(NULL);
    randgenesindex = malloc(sizeof(unsigned int *)*deg_count);
    //thispath = malloc(sizeof(struct used_path_type *));
//     int oneside = 1;
    // tablecalc(usedpaths, num_used_paths, real_universe, real_universe_cnt, oneside);
// fprintf(stderr,"in GPCC 2.3 deg_count=%d \n",deg_count); fflush(NULL);
    tablecalc(usedpaths, num_used_paths,  real_universe_cnt);
// fprintf(stderr,"in GPCC 2.4 deg_count=%d \n",deg_count); fflush(NULL);
    // FET has been done in tablecalc? Now do correction
    // gpcount_correction(usedpaths, num_used_paths);
    //*setseed = 12345; // for testing, but get an error; never malloc'd integer storage
#if NELSON_C
// fprintf(stderr,"in GPCC 2.5 deg_count=%d , NUM_TEST_PERMUTES=%d\n",deg_count,NUM_TEST_PERMUTES); fflush(NULL);

#if 1
fprintf(stderr,"doing TEST_PERMUTES\n"); 
    for (ui = 0; ui < NUM_TEST_PERMUTES; ui++)
    {
    struct used_path_type *thispath;
        // test: at ui = 0, put in the actual genes
        // Choose ingenecnt unique genes
        if ((ui % 10000) == 0) { fprintf(stderr,"permute = %u\n", ui); fflush(NULL); }
        for (I=0 ; I<num_used_paths ; I++) // test, remove
        {
            uptr = (usedpaths+I);
            if(uptr->randhits > 0)
                fprintf(stderr,"randhits not initialized\n");
        }
//fprintf(stderr,"rpfdbg here 2\n"); fflush(stderr); 
        // following returns the index of the chosen genes; this goes to pointer to all of paths for that gene
        // NEED R SWITCH HERE, FOR SEED
        if (SCALE_UNIVERSE) chooseAugRandList(aug_gene_ct, ingenecnt, randgenesindex);
        else chooseRandList(real_universe_cnt, ingenecnt, randgenesindex);
#if 0
        //for(j=0;j<ingenecnt;j++) printf("geneindex, gene, %d %d\n", randgenesindex[j], auguniverse[randgenesindex[j]]);
        /*if(ui == 1)
        for(j=0;j<ingenecnt;j++)
        {
            printf("%d\t%d\n", ui, aug_treepointers[randgenesindex[j]]->val);
        }
        printf("\n\n"); */
#endif
//fprintf(stderr,"rpfdbg here 3\n"); fflush(stderr); 
        for(j=0;j<ingenecnt;j++)
        {
    struct tree_with_count *thisnode = (struct tree_with_count *)0;
            if (SCALE_UNIVERSE) thisnode = aug_treepointers[randgenesindex[j]];
            else thisnode = treepointers[randgenesindex[j]];

            //thisnode = treepointers[randgenesindex[j]]; //no needs to point to tree;
            /*if (j < ingenecnt-1 && thisnode->val == treepointers[randgenesindex[j+1]]->val)
            {
                printf("bad\n");
            }*/
// fprintf(stderr,"j=%d of %d thisnode=%p\n",j,ingenecnt,thisnode); fflush(NULL);
// fprintf(stderr,"count=%d\n",thisnode->count); fflush(NULL);
            for(k=0;k<thisnode->count;k++)
            {
                thispath = thisnode->all_gene_paths[k];
                // (thispath+k)->randhits++;
                ++thispath->randhits;

                // printf("%s %d %d %d %d %d %d\n",thispath->acc, thispath->randhits, thisnode->count, thisnode->val, j,k, randgenesindex[j]);
                // if(!strcmp(thispath->acc, "GO:0034612") && ui == 14) printf("%d\n", thisnode->val);
            }
        }
//fprintf(stderr,"rpfdbg here 4"); fflush(stderr); 
        // now compare randhits with hits
        // could save time looking only at paths with randhits, hits
        for (I=0 ; I<num_used_paths ; I++)
        {         // get number of genes in all pathways, count duplicated gene names
            uptr = (usedpaths+I);
            /*if (!strcmp(uptr->acc, "P00001"))
            {
                printf("comparing counts   \n");
                printf("%d %d %s\n", uptr->hitcnt, uptr->randhits, uptr->acc);
            }*/
            //printf("hits and randhits %d %d\n", uptr->d, uptr->randhits);
            if (uptr->d < uptr->randhits) ++uptr->countover; // more random hits than real hits in this simulation instance; if never, p = 0
            //if(!strcmp(uptr->name, "basal lamina"))
                //printf("%11.2e %d %d %d\n", uptr->pval, uptr->hitcnt, uptr->randhits, uptr->countover);
            else if (uptr->d == uptr->randhits) ++uptr->countequal;
            else ++uptr->countunder;
            // quasi fdr; do for fewer permutes
            //printf("n, nsig %d %d\n", ui, sig_rand_p_ct);
            uptr->randhits = 0;
        }
//fprintf(stderr,"rpfdbg here 5"); fflush(stderr); 
    } // END PERMUTE LOOP
#endif
// fprintf(stderr,"in GPCC 2.6 deg_count=%d \n",deg_count); fflush(NULL);
#endif

// fprintf(stderr,"in GPCC 4.0.1\n"); fflush(NULL);
    for (I=0 ; I<num_used_paths ; I++)
    {
        uptr = (usedpaths+I);
        // trying without 0.5 x countequal
        uptr->p_permute_over  = (uptr->countover  + uptr->countequal + 0.5)/(float) NUM_TEST_PERMUTES; // low p value if few random hits > real hits
        uptr->p_permute_under = (uptr->countunder + uptr->countequal + 0.5)/(float) NUM_TEST_PERMUTES; // casts don't seem to be needed; shouldn't be for large ints
    }
// fprintf(stderr,"in GPCC 999\n"); fflush(NULL);

#if 0
  // rpf commented out because it does not do anything 
    unsigned int K;
    for (K=0; K<num_used_paths; K++) // need calcs combining path values calcd above, besides permutation calcs, no, ok to output scalar mean_deg_gpcount
    {
        uptr = (usedpaths+K);
    }
#endif
#ifdef L2P_USING_R
    //PutRNGstate(); // restore random state
#endif

    if (randgenesindex)   free(randgenesindex);
    if (head)             free_tree_with_count(head);
// fprintf(stderr,"freed %d nodes in head, num_used_paths=%u\n",hackcnt,num_used_paths);
    if (pool)             free(pool);
    if (pspace)           free(pspace);
    if (index_ptr)        free(index_ptr);
    if (aug_ptr)          free(aug_ptr);
    if (aug_treepointers) free(aug_treepointers);
    if (treepointers)     free(treepointers);
    if (gene_path_cts)    free(gene_path_cts);
    if (auguniverse)      free(auguniverse);
    if (augreverse)       free(augreverse);
    if (sig_path_ptr)     free(sig_path_ptr);
    if (deg_path_cts)     free (deg_path_cts);
// fprintf(stderr,"rpf dbug end GPCC\n"); fflush(NULL);

    return(0);
}


   // port of original gn code
static unsigned int gpcc2(struct used_path_type usedpaths[], unsigned int num_used_paths, unsigned int real_universe_cnt,  unsigned int user_incnt)
{
    struct used_path_type *uptr = (struct used_path_type *)0;
    struct used_path_type *uptr2 = (struct used_path_type *)0;
    unsigned int i,j,k,ui_ej,ui,jj;
    unsigned int A,B,C,D,a,b,c,d;
    double aug_scale = 0.0;
    double aug_pv = 1.0;
    int A_scaled, B_scaled, C_scaled, D_scaled;
    double *pvals2 = (double *)0;
    unsigned int local_aug_gene_ct = 0;
//    unsigned int newhitcnt;


fprintf(stderr,"in gpcc2(), user_incnt=%u\n",user_incnt);  fflush(stderr); 
    pvals2 = (double *)malloc((size_t)(sizeof (double)*num_used_paths) );
    if (!pvals2) { fprintf(stderr,"ERROR: no memory in gpcc2()\n"); fflush(stderr); return 1; }
    // memset(pvals2,0,sizeof(double)*num_used_paths); don't need to init 

    local_aug_gene_ct = 0;
    for ( i=0 ; i<num_used_paths ; i++ )
    {
// fprintf(stderr,"in gpcc2() rpf here 1 i=%d of %d\n",i,num_used_paths);  fflush(stderr); 
        uptr = (usedpaths+i);
        uptr->pathhits_gpsum = 0; // the number of pathways that all "HIT(!) genes" 
        uptr->pathcountsum = 0;
        local_aug_gene_ct = local_aug_gene_ct + uptr->numfixedgenes;
// rpf xxx
        if (uptr->hitcnt == 0) continue; // don't need to calculate if no hits
// fprintf(stderr,"in gpcc2() loop 1 %u to %u hit=%u\n",i,num_used_paths,uptr->hitcnt);  fflush(stderr); 
        for ( j=0 ; j<num_used_paths ; j++ )
        {         // get number of genes in all pathways, count duplicated gene names
            uptr2 = (usedpaths+j);
            jj = k = 0;
            while ( (jj<uptr->hitcnt) && (k < uptr2->numfixedgenes) )
            {
                ui_ej = *(uptr->genehits+jj);
                ui    = *(uptr2->egids+k);
                if (ui_ej == ui)
                {
                    uptr->pathhits_gpsum++;
                       // no *(uptr->genehits + newhitcnt++) = ui; already did this , for now.  *MAY* want to move calculation here, though.
                    jj++;
                    k++;
                }
                else if (ui_ej < ui) jj++;
                else                 k++;
            }
            jj = k = 0;
            while ( (jj<uptr->numfixedgenes) && (k < uptr2->numfixedgenes) )
            {
                ui_ej = *(uptr->egids+jj);
                ui    = *(uptr2->egids+k);
                if (ui_ej == ui)
                {
                    uptr->pathcountsum++;
                    jj++;
                    k++;
                }
                else if (ui_ej < ui) jj++;
                else                 k++;
            }
        }
    }
// fprintf(stderr,"in gpcc2() rpf here 2\n");  fflush(stderr); 

fprintf(stderr,"in gpcc2() 2\n");  fflush(stderr); 
    for (i=0 ; i<num_used_paths ; i++)
    {
        uptr = (usedpaths+i);

        d = uptr->hitcnt;
        c = ingenecnt - uptr->hitcnt;
        b = uptr->numfixedgenes - uptr->hitcnt;
        a = real_universe_cnt - uptr->numfixedgenes - ingenecnt + uptr->hitcnt;

        /* new way : */
        uptr->a = a;
        uptr->b = b;
        uptr->c = c;
        uptr->d = d;
        D = uptr->pathhits_gpsum; // wrong? number of pathways any and all genes in this pathway hit
        C = aug_deg_count - D; // or d?   aug_deg_count is number of times the user input genes hit in any and all pathways.  
        B = uptr->pathcountsum - D; // pathcountsum is the sum of the gene gp cts for all hits in the pathway
        A = local_aug_gene_ct - B - C - D;  //  local_aug_gene_ct is sum of all numfixedgenes in usedpaths[]
        aug_scale = (float)real_universe_cnt/((float)(local_aug_gene_ct + d)); // real_universe_cnt is count of unique genes in universe
        uptr->A_scaled = A_scaled = ceil(A*aug_scale);
        uptr->B_scaled = B_scaled = ceil(B*aug_scale);
        uptr->C_scaled = C_scaled = ceil(C*aug_scale);
        uptr->D_scaled = D_scaled = d; // not scaled!
        aug_pv = exact22(A_scaled,B_scaled,C_scaled,D_scaled);
        uptr->gpcc_p = aug_pv;     // here's the pvalue for GPCC ***
        uptr->pval = aug_pv;     // here's the pvalue for GPCC ***
 if (strcmp(uptr->acc,CHECKACC) == 0) 
 { 
 fprintf(stderr,"Got %s in gpcc2, pathhits_gpsum=%f , pathcountsum=%u ,hit=%d local_aug_gene_ct=%u aug_scale=%f\n",CHECKACC,uptr->pathhits_gpsum,uptr->pathcountsum,uptr->hitcnt,local_aug_gene_ct,aug_scale); 
 fprintf(stderr,"Got %s n gpcc2 aug_pv=%f scaled= %d %d %d %d\n",CHECKACC,aug_pv,A_scaled,B_scaled,C_scaled,D_scaled); 
 fflush(stderr); 
 } 
        *(pvals2+i) = aug_pv;
        if (C_scaled*B_scaled != 0)
        {
            uptr->gpcc_OR = ((double)A_scaled*(double)D_scaled)/((double)C_scaled*(double)B_scaled);
        }
        else uptr->gpcc_OR = 99999;
    }
    if (*pvals2) free (pvals2);
fprintf(stderr,"in gpcc2() rpf here 4\n");  fflush(stderr); 
fprintf(stderr,"in gpcc2() [end]\n");  fflush(stderr); 
    return 0;
}



struct gene_cnt_type
{
    unsigned int egid;
    unsigned int count;
    unsigned int left;
    unsigned int right;
};
struct gene_cnt_type *t;
int add_index = 0;

static void add_or_count (struct gene_cnt_type *t,unsigned int egid_arg, int idx)
{
    if (add_index == 0)
    {
        t[0].egid = egid_arg;
        t[0].count = 1;
        t[0].left = t[0].right = 0;
        add_index = 1;
        return;
    }
// fprintf(stderr,"in add_or_count() idx=%u cmp %u %u\n",idx,egid_arg,t[idx].egid);
    if (egid_arg == t[idx].egid) { t[idx].count++; return; }
    if (egid_arg<t[idx].egid)
    {
        if (t[idx].left == 0) 
        {
            t[add_index].egid = egid_arg;
            t[add_index].count = 1;
            t[add_index].left = t[add_index].right = 0;
            t[idx].left = add_index;
            add_index++;
            return;
        }
        add_or_count(t,egid_arg,t[idx].left);
    }
    else
    {
        if (t[idx].right == 0) 
        {
            t[add_index].egid = egid_arg;
            t[add_index].count = 1;
            t[add_index].left = t[add_index].right = 0;
            t[idx].right = add_index;
            add_index++;
            return;
        }
        add_or_count(t,egid_arg,t[idx].right);
    }
    return;
}


static int just_count(struct gene_cnt_type *t,unsigned int egid_arg, int idx)
{
//fprintf(stderr,"in just_count() idx=%u cmp %u %u\n",idx,egid_arg,t[idx].egid);
    if (egid_arg==t[idx].egid) 
    { 
//jfprintf(stderr,"in just_count() got %u, returning %u\n",egid_arg,t[idx].count);
          return t[idx].count;
    } 
    if (egid_arg<t[idx].egid)
    {
        return just_count(t,egid_arg,t[idx].left);
    }
    else
    {
        return just_count(t,egid_arg,t[idx].right);
    }
fprintf(stderr,"ERROR: in just_count() %u\n",egid_arg); 
    return -1;
}


   // optimized version
static unsigned int gpcc3(struct used_path_type usedpaths[], unsigned int num_used_paths, unsigned int real_universe_cnt,  unsigned int user_incnt, unsigned int aug_deg_count)
{
    struct used_path_type *uptr = (struct used_path_type *)0;
    unsigned int i,j,k;
    unsigned int A,B,C,D,a,b,c,d;
    double aug_scale = 0.0;
    double aug_pv = 1.0;
    unsigned int A_scaled, B_scaled, C_scaled, D_scaled;
    double *pvals2 = (double *)0;
    unsigned int local_aug_gene_ct = 0;
//    unsigned int newhitcnt;
    struct gene_cnt_type *gcnts;

// xxx
 fprintf(stderr,"gpcc3 augo_deg_count = %u\n",aug_deg_count);

    gcnts = (struct gene_cnt_type *)malloc(real_universe_cnt*sizeof(struct gene_cnt_type));

fprintf(stderr,"in gpcc3(), user_incnt=%u\n",user_incnt);  fflush(stderr); 
    local_aug_gene_ct = 0;
    for ( i=0 ; i<num_used_paths ; i++ )
    {
        uptr = (usedpaths+i);
        local_aug_gene_ct = local_aug_gene_ct + uptr->numfixedgenes;
        k = 0;
        while ( k < uptr->numfixedgenes)
        {
           add_or_count(gcnts, *(uptr->egids+k), 0 );
           k++;
        }
    }
    pvals2 = (double *)malloc((size_t)(sizeof (double)*num_used_paths) );
    if (!pvals2) { fprintf(stderr,"ERROR: no memory in gpcc3()\n"); fflush(stderr); free(gcnts); return 1; }
    // memset(pvals2,0,sizeof(double)*num_used_paths); don't need to init 

fprintf(stderr,"in gpcc3() 4\n");  fflush(stderr); 
    for (i=0 ; i<num_used_paths ; i++)
    {
        uptr = (usedpaths+i);
        for (a=j=0 ; j<uptr->numfixedgenes ; j++)
        {
            a = a + just_count(gcnts, *(uptr->egids+j), 0 );
        }
        uptr->pathcountsum = a;
        for (b=j=0 ; j<uptr->hitcnt ; j++)
        {
            b = b + just_count(gcnts, *((uptr->genehits)+j), 0 );
        }
        uptr->pathhits_gpsum = (double)b;

        d = uptr->hitcnt;
        c = ingenecnt - uptr->hitcnt;
        b = uptr->numfixedgenes - uptr->hitcnt;
        a = real_universe_cnt - uptr->numfixedgenes - ingenecnt + uptr->hitcnt;

        /* new way : */
        uptr->a = a;
        uptr->b = b;
        uptr->c = c;
        uptr->d = d;
        D = uptr->pathhits_gpsum; // wrong? number of pathways any and all genes in this pathway hit
        C = aug_deg_count - D; // or d?   aug_deg_count is number of times the user input genes hit in any and all pathways.  
        B = uptr->pathcountsum - D; // pathcountsum is the sum of the gene gp cts for all hits in the pathway
        A = local_aug_gene_ct - B - C - D;  //  local_aug_gene_ct is sum of all numfixedgenes in usedpaths[]
        aug_scale = (double)real_universe_cnt/((double)(local_aug_gene_ct + d)); // real_universe_cnt is count of unique genes in universe
        uptr->A_scaled = A_scaled = ceil(A*aug_scale);
        uptr->B_scaled = B_scaled = ceil(B*aug_scale);
        uptr->C_scaled = C_scaled = ceil(C*aug_scale);
        uptr->D_scaled = D_scaled = d; // not scaled!
        aug_pv = exact22(A_scaled,B_scaled,C_scaled,D_scaled);
if (strcmp(uptr->acc,CHECKACC) == 0) 
{ 
fprintf(stderr,"Got %s in gpcc3 pathhits_gpsum=%f , pathcountsum=%u ,hit=%d local_aug_gene_ct=%u aug_scale=%f\n",CHECKACC,uptr->pathhits_gpsum,uptr->pathcountsum,uptr->hitcnt,local_aug_gene_ct,aug_scale); 
fprintf(stderr,"Got %s in gpcc3 aug_pv=%20.18f scaled= a.%d b.%d c.%d d.%d ( C=%d D=%d ) aug_deg_count=%d\n",CHECKACC,aug_pv,A_scaled,B_scaled,C_scaled,D_scaled,C,D,aug_deg_count); 
fflush(stderr); 
} 
#if 0
if (d>0)
fprintf(stderr,"xxx %s : p=%f from scaled: %u %u %u %u scale=%f\n",uptr->acc,aug_pv,uptr->A_scaled , uptr->B_scaled, uptr->C_scaled, uptr->D_scaled ,aug_scale);
#endif
        uptr->gpcc_p = aug_pv;     // heres the pvalue for GPCC ***
// xxx rpf
        uptr->pval = aug_pv;     // heres the pvalue for GPCC ***
// xxx rpf
        *(pvals2+i) = aug_pv;
        if (C_scaled*B_scaled != 0)
        {
            uptr->gpcc_OR = ((double)A_scaled*(double)D_scaled)/((double)C_scaled*(double)B_scaled);
        }
        else uptr->gpcc_OR = 99999;
    }
    if (*pvals2) free (pvals2);
fprintf(stderr,"in gpcc3() [end]\n");  fflush(stderr); 
    return 0;
}

 // fast version
unsigned int gpcc4(struct used_path_type usedpaths[], unsigned int num_used_paths, unsigned int real_universe_cnt,  unsigned int user_incnt)
{
    struct used_path_type *uptr = (struct used_path_type *)0;
    unsigned int i,j,k;
    unsigned int A,B,C,D,a,b,c,d;
    double aug_scale = 0.0;
    double aug_pv = 1.0;
    int A_scaled, B_scaled, C_scaled, D_scaled;
    double *pvals2;
    unsigned int local_aug_gene_ct = 0;
//    unsigned int newhitcnt;
    struct gene_cnt_type *gcnts;

    gcnts = (struct gene_cnt_type *)malloc(real_universe_cnt*sizeof(struct gene_cnt_type));

fprintf(stderr,"in gpcc4(), user_incnt=%u\n",user_incnt);  fflush(stderr); 
    local_aug_gene_ct = 0;
    for ( i=0 ; i<num_used_paths ; i++ )
    {
        uptr = (usedpaths+i);
        local_aug_gene_ct = local_aug_gene_ct + uptr->numfixedgenes;
        k = 0;
        while ( k < uptr->numfixedgenes)
        {
           add_or_count(gcnts, *(uptr->egids+k), 0 );
           k++;
        }
    }
fprintf(stderr,"in gpcc4() 3\n");  fflush(stderr); 

    pvals2 = (double *)malloc((size_t)(sizeof (double)*num_used_paths) );
    if (!pvals2) { fprintf(stderr,"ERROR: no memory in gpcc4()\n"); fflush(stderr); free(gcnts); return 1; }
    // memset(pvals2,0,sizeof(double)*num_used_paths); don't need to init 

fprintf(stderr,"in gpcc4() 4\n");  fflush(stderr); 
    for (i=0 ; i<num_used_paths ; i++)
    {
        uptr = (usedpaths+i);
        for (a=j=0 ; j<uptr->numfixedgenes ; j++)
        {
            a = a + just_count(gcnts, *(uptr->egids+j), 0 );
        }
        uptr->pathhits_gpsum = a;
        for (b=j=0 ; j<uptr->hitcnt ; j++)
        {
            b = b + just_count(gcnts, *((uptr->genehits)+j), 0 );
        }
        uptr->pathcountsum = (double)b;
        d = uptr->hitcnt;
        c = ingenecnt - uptr->hitcnt;
        b = uptr->numfixedgenes - uptr->hitcnt;
        a = real_universe_cnt - uptr->numfixedgenes - ingenecnt + uptr->hitcnt;

        /* new way : */
        uptr->a = a;
        uptr->b = b;
        uptr->c = c;
        uptr->d = d;
        D = uptr->pathhits_gpsum; // wrong? number of pathways any and all genes in this pathway hit
        C = aug_deg_count - D; // or d?   aug_deg_count is number of times the user input genes hit in any and all pathways.  
        B = uptr->pathcountsum - D; // pathcountsum is the sum of the gene gp cts for all hits in the pathway
        A = local_aug_gene_ct - B - C - D;  //  local_aug_gene_ct is sum of all numfixedgenes in usedpaths[]
        aug_scale = (float)real_universe_cnt/((float)(local_aug_gene_ct + d)); // real_universe_cnt is count of unique genes in universe
        uptr->A_scaled = A_scaled = ceil(A*aug_scale);
        uptr->B_scaled = B_scaled = ceil(B*aug_scale);
        uptr->C_scaled = C_scaled = ceil(C*aug_scale);
        uptr->D_scaled = D_scaled = d; // not scaled!
        aug_pv = exact22(A_scaled,B_scaled,C_scaled,D_scaled);
 if (strcmp(uptr->acc,CHECKACC) == 0) 
 { 
 fprintf(stderr,"Got M40028 in gppc2, pathhits_gpsum=%f , pathcountsum=%u ,hit=%d local_aug_gene_ct=%u aug_scale=%f, hitcnt=%u aug_deg_count=%u\n",
 uptr->pathhits_gpsum,uptr->pathcountsum,uptr->hitcnt,local_aug_gene_ct,aug_scale,uptr->hitcnt,aug_deg_count); 
 fflush(stderr); 
 } 
        uptr->gpcc_p = aug_pv;     // heres the pvalue for GPCC ***
        *(pvals2+i) = aug_pv;
        if (C_scaled*B_scaled != 0)
        {
            uptr->gpcc_OR = ((double)A_scaled*(double)D_scaled)/((double)C_scaled*(double)B_scaled);
        }
        else uptr->gpcc_OR = 99999;
    }
fprintf(stderr,"in gpcc4() [end]\n");  fflush(stderr); 
    return 0;
}


#if 0
old
int fast_permutes(unsigned int ingenecnt, struct used_path_type usedpaths[], unsigned int num_used_paths,unsigned int real_universe_cnt,unsigned int *real_universe)
{
    // pthread_t thr[NUM_THREADS];
    // struct thread_data_t thr_data[NUM_THREADS];
    static unsigned int r[80002];
    double d;
    double oddsratio;
    unsigned int i,j,k;
    double *pspace;
    struct used_path_type *uptr;
    unsigned int usedpathindex;
    double *column;
    double *bestrow;
    unsigned int tmphitcnt;
    unsigned int ui_ej;
    unsigned int ui_uk;

    memcpy(r,real_universe,real_universe_cnt*sizeof(unsigned int));

    pspace = malloc(sizeof(double)*NUM_PERMUTES*num_used_paths); // no dire need to memset
    memset(pspace,0,sizeof(double)*NUM_PERMUTES*num_used_paths);

// fprintf(stderr,"debug in fast_permutes()\n"); fflush(stderr); 
    for (i=0;i<NUM_PERMUTES;i++)
    {
        subsamp(r, real_universe_cnt,ingenecnt); // inlined if optimization on 
        quickSort(r,ingenecnt);                  // inlined if optimization on 
// for (j=0 ; j<ingenecnt ; j++) fprintf(stderr,"ss: %d %d egid= %u\n",i,j,*(r+j));
        for (usedpathindex=0;usedpathindex<num_used_paths;usedpathindex++)
        {
            uptr = (usedpaths+usedpathindex);
            j = k = tmphitcnt = 0;
            while ((j<uptr->numfixedgenes) && (k < ingenecnt))
            {              // get the number of hits 
                ui_ej = *(uptr->egids+j); // entrez gene ids in the jth used pathway
                ui_uk = *(r+k);           // entrez gene from random sampling 
                if (ui_ej == ui_uk)
                {
                    tmphitcnt++;
                    k++;
                    j++;
                    continue;
                }
                else if (ui_ej < ui_uk) j++;
                else                    k++;
            }
	        /* d = exact22((int)tmphitcnt,uptr->numfixedgenes-tmphitcnt,ingenecnt,real_universe_cnt-ingenecnt);
                   *(pspace + (usedpathindex*NUM_PERMUTES) + i) = (double)d; */
	    d = ((double)tmphitcnt / (double)ingenecnt) / 
		 ((double)(uptr->numfixedgenes-tmphitcnt) / (double)(real_universe_cnt-ingenecnt));
// fprintf(stderr,"OR=%20.18f from (%u/%u)/(%u/%u)\n",oddsratio,tmphitcnt,ingenecnt,uptr->numfixedgenes-tmphitcnt,real_universe_cnt-ingenecnt);
            *(pspace + (usedpathindex*NUM_PERMUTES) + i) = (double)d;

/*
               diseased healthy
exposed           DE       HE            A   C
notexposed        DN       HN            B   D
OR=(DE/HE)/(DN/HN) = A/C / B/D

list=300 hits=3  pw=40 genome=30001
         usergenes     genome
ipw         3-1         40               tmphitcnt                        ingenecnt
notinpw     297       29960              uptr->numfixedgenes-tmphitcnt    real_universe_cnt-ingenecnt 
*/


#if 0
fprintf(stderr,"exact22 %f %d %d %d %d upi=%d i=%d %s\n",d,
		tmphitcnt,uptr->numfixedgenes-tmphitcnt,ingenecnt,real_universe_cnt-ingenecnt,
		usedpathindex,i,
		uptr->name);
#endif
        }

    }

#if 0
 for (i=0 ; i<NUM_PERMUTES ; i++) 
 { 
         fprintf(stderr,"perm:%d ",i);
         for (usedpathindex=0 ; usedpathindex<num_used_paths ; usedpathindex++)
         {
             d = *(pspace + (usedpathindex*NUM_PERMUTES) + i);
             fprintf(stderr,"%9.7f ",d);
         }
         fprintf(stderr,"\n");
 }
#endif

    column = malloc(sizeof(double)*num_used_paths); // no dire need to memset
    bestrow = malloc(sizeof(double)*NUM_PERMUTES); // no dire need to memset
    for (i=0;i<NUM_PERMUTES;i++)
    {
        for (usedpathindex=0;usedpathindex<num_used_paths;usedpathindex++)
        {
            d = *(pspace + (usedpathindex*NUM_PERMUTES) + i);
            *(column+usedpathindex) = d;
        }
        qsort(column,num_used_paths,sizeof(double),cmp_double);
        *(bestrow+i) = *(column+0);
    }
    qsort(bestrow,NUM_PERMUTES,sizeof(double),cmp_double);
#if 0
for (j=0 ; j<NUM_PERMUTES ; j++) { fprintf(stderr,"%d %f\n",j,*(bestrow+j)); }
#endif
    d = 1.0;
    for (usedpathindex=0 ; usedpathindex<num_used_paths ; usedpathindex++)
    {
        uptr = (usedpaths+usedpathindex);
        oddsratio = ((double)uptr->hitcnt / (double)ingenecnt) /
		 ((double)(uptr->numfixedgenes-(uptr->hitcnt)) / (double)(real_universe_cnt-ingenecnt));
        for (j=0 ; j<NUM_PERMUTES ; j++) 
        {
            if ( *(bestrow+j) < oddsratio )      // if (*(bestrow+j) > uptr->pval) 
            {
                d = *(bestrow+j); 
                break;
            }
        }
        d = (double)j/(double)NUM_PERMUTES;
        uptr->fdr2 = d;
    }

    free(pspace);
    free(column);
    free(bestrow);

    return 0;

#if 0
        k = NUM_PERMUTES/NUM_THREADS;
        for (j = 0; j < NUM_THREADS; ++j) 
        {     /* create threads */
            thr_data[j].lo = j*k;
            thr_data[j].hi = (j*k)+k-1;
            if (j == (NUM_THREADS-1))
                thr_data[j].hi = num_used_paths-1; // note "-1", 
fprintf(stderr,"%d of %d [ %d .. %d ] of %d size= %d\n",j,NUM_THREADS,thr_data[j].lo,thr_data[j].hi,num_used_paths,thr_data[j].hi-thr_data[j].lo); 
            thr_data[j].myid = j;
            thr_data[j].ingenecnt = ingenecnt;
            thr_data[j].sampling = r;
            thr_data[j].usedpathsindex = usedpaths;
            thr_data[j].usedpaths = usedpaths;
            thr_data[j].num_used_paths = num_used_paths;
fprintf(stderr,"creating thread %d\n",j); fflush(stderr); 
            if ((rc = pthread_create(&thr[j], NULL, thr_func, &thr_data[j]))) 
            {
              fprintf(stderr, "error: pthread_create, rc: %d\n", rc);
              return EXIT_FAILURE;
            }
        }
      /* block until all threads complete */
        for (j = 0; j < NUM_THREADS; ++j) 
        {
fprintf(stderr,"blocking %d \n",j);
            pthread_join(thr[j], NULL);
        }
fprintf(stderr,"after blocking  \n");
    return EXIT_SUCCESS;
#endif
}
#endif
 
int bh4( struct used_path_type usedpaths[], unsigned int num_used_paths)
{
    unsigned int i;
    double *pvals;
    double *fdrs;


// fprintf(stderr,"debug in bh4() ingenecnt=%d num_used_paths=%d real_universe_cnt=%d \n",ingenecnt,num_used_paths,real_universe_cnt); fflush(stderr); 

    pvals = (double *)malloc((size_t)(sizeof (double)*(num_used_paths)));
    if (!pvals) { fprintf(stderr,"ERROR: no memory\n"); fflush(stderr); return -1; }
    for (i=0 ; i<num_used_paths ; i++)
         *(pvals+i) = usedpaths[i].pval4; // nota bene
    fdrs = (double *)malloc((size_t)(sizeof (double)*num_used_paths)); 
    if (!fdrs) { free(pvals); /* clean up */ fprintf(stderr,"ERROR: no memory in bh4() 2\n"); return -3; }
    benjaminihochberg(num_used_paths,pvals,fdrs);
    for (i=0 ; i<num_used_paths ; i++)
    {
         usedpaths[i].fdr = *(fdrs+i);
    }
    free(pvals);
    free(fdrs);
    return 0;
}

          /* Be sure to free egids if you call setup_by_egids() !!! */
int setup_by_egids(void)
{
    size_t sz = sizeof(struct smallgenetype) * numgenes;

    by_egids = (struct smallgenetype *)malloc(sz);
    if (!by_egids) return -1;
    memcpy(by_egids,genes,sz);
    qsort(by_egids,numgenes,sizeof(struct smallgenetype),cmp_by_egid);
    return 0;
}

unsigned int  tmpgenes[MAX_INGENES];
 
int l2pfunc(struct used_path_type *usedpaths,unsigned int num_used_paths,unsigned int real_universe_cnt,
             unsigned int *real_universe, int calc_option, int *user_incnt_ptr, int oneside, unsigned int num_permutes)
{
    char s[512];
    unsigned int prev;
    unsigned int ui;
    unsigned int user_incnt = 0;
    unsigned int i,j;
    int ret = 0;
    int overflowflag = 0;
    struct used_path_type *uptr; // used path pointer 
    unsigned int ui_ej;
    unsigned int k,ll;
    double *pvals;
    double *fdrs;
// struct timespec time1, time2;
    // char *z;

fprintf(stderr,"rpf perm in l2pfunc() start, calc_option=%d num_permutes=%u\n",calc_option,num_permutes); fflush(stderr);  

    *user_incnt_ptr = j = 0;
    // ADDING INPUT FILE FOR C RUN, INPUT GENE LIST
    // struct timespec time1, time2;
    // char *z;
    // printf("input DEG list\n");
#if 0
 // # if NELSON_TEST
    FILE *deglist = (FILE *)0;
    if ((deglist = fopen (DEGLISTNAME, "r")) == NULL){
        printf ("ERROR: Can not open deglist file named \"%s\" exiting",DEGLISTNAME);
        exit (1);
    }
    // READING FROM DEGLIST
    while ( fgets(s, sizeof(s), deglist) ) // gets() function is deprecated
#else
    while ( fgets(s, sizeof(s), stdin) ) // gets() function is deprecated
#endif
    {
        for (i=0 ; s[i] ; i++) 
        { 
             if ((s[i] == '\n')||(s[i] == '\r')) s[i] = (char)0; 
             else if ((s[i] < ' ') || (s[i] >= 127)) s[i] = ' ';
        }
        ui = hugo2egid(s);
        if (ui == (unsigned int)UINT_MAX)
        {
             fprintf(stderr,"Note: invalid gene \"%s\" in user input\n",s);  
             continue;
        }
#if 0
        if (bsearch2(ui,real_universe,real_universe_cnt) == 0)
        {
             fprintf(stderr,"Note: gene not in universe : \"%s\" \n",s);  
             continue;
        }
        if (j < MAX_INGENES)
        {
            tmpgenes[j++] = ui;
        }
        else
        {
            if (overflowflag == 1) { fprintf(stderr,"NOTE: too many genes, input gene number limited to %d, rest are ignored, incnt=%d\n",MAX_INGENES,j); fflush(NULL); }
            overflowflag = 1;
        }
#else
        unsigned int *uiptr = (void *)0;
        uiptr = (unsigned int *)bsearch(&ui,real_universe,real_universe_cnt,sizeof(unsigned int),cmp_ui);
        if (!uiptr)
        {
             fprintf(stderr,"Note: gene not in universe : \"%s\" \n",s);  
             continue;
        }
        if (j < MAX_INGENES)
        {
            tmpgenes[j++] = ui;
        }
        else
        {
            if (overflowflag == 1) { fprintf(stderr,"NOTE: too many genes, input gene number limited to %d, rest are ignored, incnt=%d\n",MAX_INGENES,j); fflush(NULL); }
            overflowflag = 1;
        }
#endif
    }
#if 0
    if (deglist) { fclose(deglist); deglist = (FILE *)0; }
#endif
#if RADIX_SORT
    radix_ui(tmpgenes,j);
#else
    qsort(tmpgenes,j,sizeof(unsigned int),cmp_ui); // sort user input egids (entrez gene ids)
#endif
    prev = 0;
    for (i=user_incnt=0;i<j;i++)
    {                    // de-duplicate 
        ui = tmpgenes[i];
        if (ui != prev)
        {
            ingenes[user_incnt++] = ui;
        }
        prev = ui;
    }
    ingenecnt = user_incnt;

// clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
#if 0
clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
#define BILLION  1000000000.0;
double time_spent = (time2.tv_sec - time1.tv_sec) + (time2.tv_nsec - time1.tv_nsec) / BILLION;
fprintf(stderr,"Time elapsed is %f seconds\n", time_spent);  fflush(NULL);

fprintf(stderr,"here , before do_pvals %d %d %p %d \n",real_universe_cnt , user_incnt, usedpaths,num_used_paths); 
fflush(NULL);
#endif

    *user_incnt_ptr = user_incnt;

// fprintf(stderr,"rpf calc_option=%d\n",calc_option); 
    if (calc_option == CALC_OPTION_FE) // fe
    {
// fprintf(stderr,"rpf num_used_paths=%d\n",num_used_paths); 
        for (i=0 ; i<num_used_paths ; i++)
        {
            uptr = (usedpaths+i);
            if (!uptr->egids) continue;
            j = k = ll = 0;
            while ((j<uptr->numfixedgenes) && (k < user_incnt))
            {
                ui_ej = *(uptr->egids+j);
                ui = ingenes[k];
                if (ui_ej == ui)
                {
                    *((uptr->genehits) + (ll++)) = ui; // remember, because need to print out later
                    k++;
                    j++;
                    // aug hit count for aug contingency table
                    continue;
                }
                else if (ui_ej < ui) j++;
                else                 k++;
            }
            uptr->hitcnt = ll;
        }
/// where is lp2func ?
        do_pvals_and_bh(user_incnt, usedpaths,num_used_paths, real_universe_cnt,oneside);
    }
    else if (calc_option == CALC_OPTION_GPCC) // gpcc
    {
fprintf(stderr,"before nelson GPCC\n"); 
        GPCC(usedpaths,num_used_paths,real_universe_cnt, real_universe);
    // RICH: this code does all comparisons out of GPCC
        do_just_bh(user_incnt,usedpaths,num_used_paths,real_universe_cnt);
// rpf fix
    }
    else if (calc_option == CALC_OPTION_GPCC2) // gpcc2
    {
        for (i=0 ; i<num_used_paths ; i++)
        {
            uptr = (usedpaths+i);
            if (!uptr->egids) continue;
            j = k = ll = 0;
            while ((j<uptr->numfixedgenes) && (k < user_incnt))
            {
                ui_ej = *(uptr->egids+j);
                ui = ingenes[k];
                if (ui_ej == ui)
                {
                    *((uptr->genehits) + (ll++)) = ui; // remember, because need to print out later
                    k++;
                    j++;
                    // aug hit count for aug contingency table
                    continue;
                }
                else if (ui_ej < ui) j++;
                else                 k++;
            }
            uptr->hitcnt = ll;
        }
        gpcc2(usedpaths,num_used_paths,real_universe_cnt, user_incnt);
    // RICH: this code does all comparisons out of GPCC
        do_just_bh(user_incnt,usedpaths,num_used_paths,real_universe_cnt);
// rpf fix
    }
    else if (calc_option == CALC_OPTION_GPCC3) // gpcc3
    {
        unsigned int aug_deg_count = 0;
        for (i=0 ; i<num_used_paths ; i++)
        {
            uptr = (usedpaths+i);
            if (!uptr->egids) continue;
            j = k = ll = 0;
            while ((j<uptr->numfixedgenes) && (k < user_incnt))
            {
                ui_ej = *(uptr->egids+j);
                ui = ingenes[k];
                if (ui_ej == ui)
                {
                    *((uptr->genehits) + (ll++)) = ui; // remember, because need to print out later
                    k++;
                    j++;
                    aug_deg_count++;
                    // aug hit count for aug contingency table
                    continue;
                }
                else if (ui_ej < ui) j++;
                else                 k++;
            }
            uptr->hitcnt = ll;
        }
        gpcc3(usedpaths,num_used_paths,real_universe_cnt, user_incnt,aug_deg_count);
        fdrs = (double *)malloc((size_t)(sizeof (double)*num_used_paths)); 
//        benjaminihochberg(num_used_paths,pvals,fdrs);
//        do_just_bh(user_incnt,usedpaths,num_used_paths,real_universe_cnt);
        pvals = (double *)malloc((size_t)(sizeof (double)*(num_used_paths)));
        fdrs = (double *)malloc((size_t)(sizeof (double)*(num_used_paths)));
        //do_pvals_and_bh(user_incnt, usedpaths,num_used_paths, real_universe_cnt,oneside);
        for (i=0 ; i<num_used_paths ; i++) *(pvals+i) = usedpaths[i].gpcc_p;
        benjaminihochberg(num_used_paths,pvals,fdrs);
        for (i=0 ; i<num_used_paths ; i++) usedpaths[i].fdr = *(fdrs+i);
        free(pvals);
        free(fdrs);
fprintf(stderr,"rpf after benjaminihochberg\n"); fflush(stderr); 
// fprintf(stderr,"%d %f\n", i,usedpaths[i].fdr);
// rpf fix
    }
    else if (calc_option == CALC_OPTION_GPCC4) // gpcc4
    {
        for (i=0 ; i<num_used_paths ; i++)
        {
            uptr = (usedpaths+i);
            if (!uptr->egids) continue;
            j = k = ll = 0;
            while ((j<uptr->numfixedgenes) && (k < user_incnt))
            {
                ui_ej = *(uptr->egids+j);
                ui = ingenes[k];
                if (ui_ej == ui)
                {
                    *((uptr->genehits) + (ll++)) = ui; // remember, because need to print out later
                    k++;
                    j++;
                    // aug hit count for aug contingency table
                    continue;
                }
                else if (ui_ej < ui) j++;
                else                 k++;
            }
            uptr->hitcnt = ll;
        }
// unsigned int gpcc2(struct used_path_type usedpaths[], unsigned int num_used_paths, unsigned int real_universe_cnt, unsigned int *real_universe, const unsigned int user_ingenes[], unsigned int user_incnt)
        gpcc4(usedpaths,num_used_paths,real_universe_cnt, user_incnt);
//        do_just_bh(user_incnt,usedpaths,num_used_paths,real_universe_cnt);
        pvals = (double *)malloc((size_t)(sizeof (double)*(num_used_paths)));
        fdrs = (double *)malloc((size_t)(sizeof (double)*(num_used_paths)));
        //do_pvals_and_bh(user_incnt, usedpaths,num_used_paths, real_universe_cnt,oneside);
        for (i=0 ; i<num_used_paths ; i++) *(pvals+i) = usedpaths[i].gpcc_p;
        benjaminihochberg(num_used_paths,pvals,fdrs);
        for (i=0 ; i<num_used_paths ; i++) usedpaths[i].fdr = *(fdrs+i);
        free(pvals);
        free(fdrs);
fprintf(stderr,"rpf after benjaminihochberg\n"); fflush(stderr); 
// fprintf(stderr,"%d %f\n", i,usedpaths[i].fdr);
// rpf fix
    }
    else if (calc_option == CALC_OPTION_PERMUTE) // permute
    {
fprintf(stderr,"in l2pfunc() permute2 rpf num_used_paths=%d, calc_option=2 (CALC_OPTION_PERMUTE) num_permutes=%d\n",num_used_paths,num_permutes); 
// xxx
        if (num_permutes < 100) 
        {
fprintf(stderr,"Warning: overridd num_permutes of %d.  Setting to default %d (PERMUTE_DEFAULT_LOW).\n",num_permutes,PERMUTE_DEFAULT_LOW); fflush(NULL);
            num_permutes = PERMUTE_DEFAULT_LOW;
        }
fprintf(stderr,"rpf permute2 in lp2func(), after permute2, num_permutes=%d\n",num_permutes); fflush(NULL);
        permute2(usedpaths,num_used_paths,user_incnt,num_permutes);
fprintf(stderr,"rpf permute2 in lp2func(), after permute2, before bh4()\n"); fflush(NULL);
        bh4(usedpaths,num_used_paths);
fprintf(stderr,"rpf permute2 in lp2func(), after bh4\n"); fflush(NULL);
    }
    else if (calc_option == CALC_OPTION_PERMUTE3) // permute
    {
fprintf(stderr,"in l2pfunc() permute3 rpf num_used_paths=%d, calc_option=2\n",num_used_paths); 
        for (i=0 ; i<num_used_paths ; i++)
        {
            uptr = (usedpaths+i);
            if (!uptr->egids) continue;
            j = k = ll = 0;
            while ((j<uptr->numfixedgenes) && (k < user_incnt))
            {
                ui_ej = *(uptr->egids+j);
                ui = ingenes[k];
                if (ui_ej == ui)
                {
                    *((uptr->genehits) + (ll++)) = ui; // remember, because need to print out later
                    k++;
                    j++;
                    // aug hit count for aug contingency table
                    continue;
                }
                else if (ui_ej < ui) j++;
                else                 k++;
            }
            uptr->hitcnt = ll;
        }
// xxx
fprintf(stderr,"rpf permute3 in lp2func(), before permute3\n"); fflush(NULL);
        if (num_permutes < 100) 
        {
fprintf(stderr,"Warning: overridd num_permutes of %d.  Setting to default %d (PERMUTE_DEFAULT_LOW).\n",num_permutes,PERMUTE_DEFAULT_LOW); fflush(NULL);
            num_permutes = PERMUTE_DEFAULT_LOW;
        }
        permute3(usedpaths,num_used_paths,user_incnt,num_permutes);
fprintf(stderr,"rpf permute3 in lp2func(), after permute3, before bh4()\n"); fflush(NULL);
        bh4(usedpaths,num_used_paths);
fprintf(stderr,"rpf permute3 in lp2func(), after bh4\n"); fflush(NULL);
    }
    return ret;
}



int bitCount(int n)
{
    int cnt = 0;
    while (n)
    {
        cnt += n % 2;
        n >>= 1;
    }
    return cnt;
}

void usage(void)
{
fprintf(stderr,"l2p : \"list to pathways\" program.\n");
fprintf(stderr,"Usage: cat listofHUGOgenes_one_per_line.txt | l2p [optional args]\n");
fprintf(stderr,"possible optional args are ...\n");
fprintf(stderr," -help\n");
fprintf(stderr," -precise\n");
fprintf(stderr," -justheader\n");
fprintf(stderr," -noheader\n");
fprintf(stderr," -universe=Universefile_one_gene_per_line.txt\n");
fprintf(stderr," -categories=listofcatgories  (comma separated)\n");
fprintf(stderr,"    example: -categories=KEGG,REACTOME,BIOCYC,PANTH  (i.e. only use genes in those 4 categories\n");
fprintf(stderr,"    another -categories example: \"-categories=H,C6\"  (i.e. only use msigdb's Hallmark and C6 (cancer) category pathways\n");
fprintf(stderr,"    available categories are :\n");
fprintf(stderr,"    BIOCYC  - organism specific Pathway/ Genome Databases (PGDBs)  - https://biocyc.org/\n");
fprintf(stderr,"    GO  - initiative to unify representation of gene and gene product attributes -  http://geneontology.org\n");
fprintf(stderr,"    KEGG - databases dealing with genomes, biological pathways, - https://www.kegg.jp/\n");
fprintf(stderr,"    PANTH - databases for protein analysis through evolutionary relationships - http://www.pantherdb.org/\n");
fprintf(stderr,"    PID  - Pathway interaction database: legacy database from Carl Schaefer & buddies at NCI\n");
fprintf(stderr,"    REACTOME - curated database of biological pathways - https://reactome.org/\n");
fprintf(stderr,"    WikiPathways - community resource for biological pathways - https://www.wikipathways.org\n");
fprintf(stderr,"    C1 - MSigDB only, positional gene sets for each human chromosome and cytogenetic band.\n");
fprintf(stderr,"    C2 - MSigDB only, curated gene sets from online pathway databases, publications in PubMed, and experts.\n");
fprintf(stderr,"    C3 - MSigDB only, motif gene sets based on conserved cis-regulatory motifs from comparative analysis\n");
fprintf(stderr,"    C4 - MSigDB only, computational gene sets defined by mining large collections of cancer-oriented microarray data.\n");
fprintf(stderr,"    C5 - MSigDB only, gene sets consist of genes annotated by the same GO terms.\n");
fprintf(stderr,"    C6 - MSigDB only, oncogenic gene sets defined directly from microarray data from cancer gene perturbations.\n");
fprintf(stderr,"    C7 - MSigDB only, immunological signatures: represents cell states and perturbations within the immune system.\n");
fprintf(stderr,"    C8 - MSigDB only, markers identified in single-cell sequencing studies of human tissue.\n");
fprintf(stderr,"    H - MSigDB only, hallmark gene sets: signatures from MSigDB gene sets to represent biological processes.\n");
fprintf(stderr,"Example run : printf \"TP53\\nPTEN\\nAPC\\nKRAS\\nNRAS\\n\" | ./l2p -categories=PID | sort -k2,2n -k1,1n -k3,3nr | head\n");
            fflush(NULL);
}

static int parsecats(char *z, unsigned int *catspat)
{
    char ts[16][16];
    int bit;
    int j,k;
    int toks;

    for (j=0 ; *(z+j) ; j++)
    {
        if (*(z+j) == ',') *(z+j) = ' ';
    }
    memset(ts,0,sizeof(ts));
    toks = sscanf(z,"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s ",
          ts[0], ts[1], ts[2], ts[3], ts[4], ts[5], ts[6], ts[7], ts[8], ts[9],
          ts[10], ts[11], ts[12], ts[13], ts[14], ts[15]);
    for (k=0;k<toks;k++)
    {
        bit=string_to_category_code(ts[k]);
        if (bit)
            *catspat |= bit;
        else
        {
            fprintf(stderr,"ERROR: invalid category = %s\n",ts[k]);
            usage();
            exit(0);
        }
    }
    j = bitCount(*catspat) ;
    if (j == 0)
    {
#ifdef L2P_USING_R
        return 0;
#else
        fprintf(stderr,"ERROR: no categories specified\n");
        usage();
        return -1;
#endif
    }
    return 0;
}


void print_header(void)
{
 /* 1  */ printf("pval\t"); 
 /* 2  */ printf("fdr\t");
 /* 3  */ printf("enrichment_score\t");        // if positive, genes are OVER REPRESENTED, if negative genes are UNDER REPRESENTED
 /* 4  */ printf("pwhitcount\t");              // number of genes hit in pathway
 /* 5  */ printf("pwnohitcount\t");            // pathway number of genes in the pathway
 /* 6  */ printf("inputcount\t");              // total count of user genes (user input)
 /* 7  */ printf("pwuniverseminuslist\t");     // total number of unique genes in all pathways
 /* 8  */ printf("pathwayaccessionidentifier\t");// canonical accession ( if available, otherwise assigned by us )
 /* 9  */ printf("category\t");                // KEGG,REACTOME,GO,PANT(=PANTHER),PID=(pathway interaction database)  *was "source"*
 /* 10 */ printf("pathwayname\t");             // Name of pathway
            // last field, no tab
 /* 11 */ printf("genes");               // genes_space_separated   HUGO genes from user that hit the pathway
          printf("\n");
}

int l2p_init_C(int argc,char *argv[],unsigned int *catspatptr,int *precise_flag, int *calc_option,int *no_header_flag, 
       char universe_file[], char custom_file[], int *get_universe_flag, int *oneside, unsigned int *seed, unsigned int *num_permutes)
{
    char *z;
    int i;
    int oneerr = 0;

    *get_universe_flag = 0;
    *precise_flag = 0;
    *calc_option = CALC_OPTION_FE; // deafult is fisher's exact test
    *no_header_flag = 0;
    *oneside = 0; // default is two-sided fe
    *seed = 0;
    universe_file[0] = (char)0;
    custom_file[0] = (char)0;
    for (i=1 ; i<argc ; i++)
    {
        if ((strcmp(argv[i],"-help") == 0) || (strcmp(argv[i],"--help") == 0) || ((strcmp(argv[i],"help") == 0) ) )
        {
             usage();
             exit(1);
        }
        else if (strncmp(argv[i],"-categories=",12) == 0)
        {
            z = argv[i],
            z += 12;
            parsecats(z,catspatptr);
        }
        else if (strcmp(argv[i],"-precise") == 0)
              *precise_flag = 1; // print more digits out so user doesn't complain about "real pvals"
        else if (strcmp(argv[i],"-gpcc") == 0)
              *calc_option = CALC_OPTION_GPCC; 
        else if (strcmp(argv[i],"-permute") == 0)
              *calc_option = CALC_OPTION_PERMUTE;
        else if (strcmp(argv[i],"-permute3") == 0)
              *calc_option = CALC_OPTION_PERMUTE3;
        else if (strcmp(argv[i],"-gpcc2") == 0)
              *calc_option = CALC_OPTION_GPCC2;
        else if (strcmp(argv[i],"-gpcc3") == 0)
              *calc_option = CALC_OPTION_GPCC3;
        else if (strcmp(argv[i],"-noheader") == 0)
              *no_header_flag = 1;
        else if (strcmp(argv[i],"-justheader") == 0)
        {
            print_header();
            exit(0);
        }
        else if (strncmp(argv[i],"-oneside=",9) == 0)
        {
            oneerr = 0;
            if (argv[i] + 9)
            {
                *oneside = atoi(argv[i] + 9);
                if ( (*oneside < 0) || (*oneside > 2) )
                    oneerr = 1;
            }
            else 
            {
                oneerr = 1;
            }
            if (oneerr)
            {
                oneside = 0;
                fprintf(stderr,"ERROR: invalid one sided flag , using two-sided\n"); fflush(stderr);
            }
        }
        else if (strncmp(argv[i],"-universe=",10) == 0)
        {
            strcpy(universe_file,argv[i] + 10);
            // user_universe_flag = 1;
fprintf(stderr,"note: using \"%s\" as universe file\n",universe_file);
        }
        else if (strncmp(argv[i],"-seed=",6) == 0)
        {
            *seed = (unsigned int)atoi(argv[i]+6);
fprintf(stderr,"seed is %u\n",*seed); 
        }
        else if (strcmp(argv[i],"-getuniverse") == 0)
        {
            *get_universe_flag = 1;
        }
        else if (strncmp(argv[i],"-numpermutes=",13) == 0)
        {
            *num_permutes = atoi(argv[i]+13);
//             customflag++;
fprintf(stderr,"note: num_permutes = %u\n",*num_permutes);
        }
        else if (strncmp(argv[i],"-customfile=",12) == 0)
        {
            strcpy(custom_file,argv[i] + 12);
//             customflag++;
fprintf(stderr,"note: using \"%s\" as custom pathway gene listfile\n",custom_file);
        }
        else 
        {
fprintf(stderr,"Note: Invalid argument \"%s\" : ignored.\n",argv[i]);
        }
    }
    if (*catspatptr == 0) category_set_all(catspatptr);
    return 0;
}


void print_universe(unsigned int real_universe_cnt, unsigned int *ru)
{              // not a debug routine, used in command line
    unsigned int i;
    char *z;
// fprintf(stderr,"in print_universe real_universe_cnt = %d \n",real_universe_cnt);  fflush(stderr); 
    for (i=0;i<real_universe_cnt;i++)
    {
        z = egid2hugo(*(ru+i));
        if (z) printf("%s\n",z);
    }
    return;
}


unsigned int *get_used_universe(struct used_path_type *u, unsigned int num_used, unsigned int *real_universe_cnt_ptr)
{
    unsigned int i,j,k;
    size_t sz = 0;
    unsigned int prev = 0;
    unsigned int ui = 0;
    unsigned int maxelems = 0;
    unsigned int *tempverse  = (unsigned int *)0;
    unsigned int *real_universe  = (unsigned int *)0;
    struct used_path_type *uptr = (void *)0;    // used path pointer


    if (!u)
    {
        fprintf(stderr,"ERROR: null used_path_type passed to get_used_universe()\n"); fflush(stderr);
        return (unsigned int *)0 ;
    }

    maxelems = numpwgenes * 2;
    sz = (maxelems * sizeof(unsigned int)); //  maximum possible number of genes
    tempverse = malloc(sz);
    if (!tempverse)
    {
        fprintf(stderr,"ERROR: not enough memory in get_used_universe()\n"); fflush(stderr);
        return (unsigned int *)0 ;
    }
    memset(tempverse,0,sz);
    for (i=k=0;i<num_used;i++)
    {
        uptr = (u+i);
        for (j=0 ; j < uptr->numfixedgenes ; j++)
        {
            *(tempverse+k) = *(uptr->egids + j);
            k++;
            if (k >= maxelems)
            {
                fprintf(stderr,"in get_used_universe ERROR i=%d k=%d\n",i,k);  fflush(stderr);
                exit(0);
            }
        }
    }

#if RADIX_SORT
    radix_ui(tempverse,k);
#else
    qsort(tempverse,k,sizeof(unsigned int),cmp_ui);
#endif
    real_universe = (unsigned int *)malloc(sz);
    prev = 0;
    for (i=j=0;i<k;i++)
    {
       ui = *(tempverse+i);
       if (ui != prev)
       {
           *(real_universe+j) = ui;
           j++;
       }
       prev = ui;
    }
    free(tempverse);
    *real_universe_cnt_ptr = j;
    return (unsigned int *)real_universe;
}

unsigned int *load_user_universe_file(char universe_file[], unsigned int *universe_cnt_ret)
{
    char s[130000];   // 127570 is number of bytes in universe file, so should handle anything
    unsigned int ui;
    FILE *universe_fp;
    unsigned int *universe_genes = (unsigned int *)0;
    unsigned int universe_cnt;
    int i,j;

    *universe_cnt_ret = j = 0;
    if (universe_file[0])
    {                 // we have a universe file
        universe_fp = fopen(universe_file,"r");
        if (!universe_fp) 
        {
fprintf(stderr,"ERROR: Can't open universe_file named \"%s\"\n",universe_file); 
            return (void *)0;
        }
        universe_cnt = 0;
        while ( fgets(s, sizeof(s), universe_fp) )
            universe_cnt++;
        universe_genes = malloc(sizeof(unsigned int *)*universe_cnt);
        rewind( universe_fp);
        j = 0;
        while ( fgets(s, sizeof(s), universe_fp) )
        {
            for (i=0 ; s[i] ; i++) { if ((s[i] == '\n')||(s[i] == '\r')) s[i] = (char)0; }
// if (ui==29974) { fprintf(stderr,"got A1CF 29974\n"); exit(0); }
            ui = hugo2egid(s);
            if (ui == (unsigned int)UINT_MAX)
            {
                fprintf(stderr,"Note: invalid gene \"%s\" in universe file.\n",s);  
                continue;
            }
            *(universe_genes+j) = ui;
            j++;
        }
        fclose (universe_fp);
    }
    qsort(universe_genes,j,sizeof(unsigned int),cmp_ui);
    *universe_cnt_ret = j; // repair this , some genes may have been "invalid". initial count may have been wrong.
    return universe_genes;
}

 
int count_lines_in_file(char *fn)
{
    FILE *fp;
    int lines = 0;
    int ch = 0;

    fp = fopen(fn,"r");
    if (!fp) return 0;
    while(!feof(fp))
    {
        ch = fgetc(fp);
        if (ch == '\n') lines++;
    }
    fclose(fp);
    return lines;
}


struct used_path_type *setup_used_paths(unsigned int *num_used_paths, unsigned int catspat, char universe_file[], unsigned int in_universe_cnt,unsigned int *in_universe, char custom_file[], unsigned int *real_universe_cnt_ptr,unsigned int **real_universe,unsigned int lencust,struct custom_type *mycustompw)
{
    char s[100000];
    char tmps[5120];
    char custname[5120];
    unsigned int i,j,k,ll,used_index;
    int tokcnt;
    struct custom_type *kust = (struct custom_type *)0;
    char *z = (void *)0;
    int this_cust_pw_cnt = 0;
    unsigned int ui = 0;
    unsigned int *uiptr = (void *)0;
    unsigned short int usi = (unsigned short int)0;
    struct used_path_type *u = (struct used_path_type *)0;
    struct used_path_type *uptr; // used path pointer
    unsigned int custom_cnt = 0;
    unsigned int universe_cnt = 0;
    unsigned int *user_universe_genes = (unsigned int *)0;
    unsigned int *this_egids = (unsigned int *)0;
    unsigned int *tmpugenes = (unsigned int *)0;
    unsigned int prev = 0;
    unsigned int newcnt = 0;
    int keepgoing;
    FILE *fp;
// char **custom_data = (void *)0;
// struct timespec time1, time2;
// fprintf(stderr,"rpf debug mycustompw=%p lencust=%d\n",mycustompw,lencust); fflush(stderr);  
// fprintf(stderr,"in setup_used_paths() start catspat=%d = 0x%x\n",catspat,catspat); fflush(stderr); 
        // clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    *real_universe = (unsigned int *)0;

//fprintf(stderr,"in setup_used_paths\n"); fflush(stderr);

              // just need to know how many lines in gmt file
    if (custom_file[0]) 
    {
        custom_cnt = count_lines_in_file(custom_file);
// fprintf(stderr,"rpf debug custom_cnt=%d from %s\n",custom_cnt,custom_file); fflush(stderr);  
    }
    else
    {
        if (mycustompw) 
        {
// fprintf(stderr,"rpf debug mycustompw=%p\n",mycustompw); fflush(stderr);  
        }
    }

    if (universe_file[0])
    {
        user_universe_genes = load_user_universe_file(universe_file,&universe_cnt); // may return null 
    }
    else if (in_universe_cnt != 0)
    {
        user_universe_genes = in_universe;
	universe_cnt = in_universe_cnt;
    }
    u = (struct used_path_type *)malloc(sizeof(struct used_path_type) * (numpws+custom_cnt+lencust));
    if (!u) { fprintf(stderr,"ERROR: not enough memory.\n"); return (void *)0; }
    memset(u,0,sizeof(struct used_path_type) * (numpws+custom_cnt+lencust));

    used_index = 0;
    for (i=0 ; i<numpws ; i++)
    {
// fprintf(stderr,"in GOTEST i=%d pws[i].category=0x%x pat=0x%x\n",i,pws[i].category,catspat); 
        if ((pws[i].category & catspat) == 0) continue;
// fprintf(stderr,"in GOTEST KEEP i=%d pws[i].category=0x%x pat=0x%x\n",i,pws[i].category,catspat); 
        this_egids = (unsigned int *)malloc(sizeof(unsigned int)*pws[i].numgenes);
        if (this_egids == (void *)0)
        {
                 exit(0);
        }
        if (universe_cnt)   // get rid of unused genes - use only universe genes
        {
            for (newcnt = j = 0 ; j < pws[i].numgenes ; j++)
            {
                usi = pwgenes[pws[i].pwgenesindex+j];
                ui = genes[usi].egid;
#if 0
                if (bsearch2(ui,user_universe_genes,universe_cnt) == 1)
                {
                    *(this_egids+newcnt) = ui; // put this in the right place
                    newcnt++;
                }
                else
                {
                }
#else
                uiptr = (unsigned int *)bsearch(&ui,user_universe_genes,universe_cnt,sizeof(unsigned int),cmp_ui);
                if (uiptr)
                {
                    *(this_egids+newcnt) = ui; // put this in the right place
                    newcnt++;
                }
                else
                {
                }
#endif
            }
            if (newcnt < 3)  // not enough genes
            {
                free(this_egids);
                this_egids = (unsigned int *)0;
                continue;
            }
        }
        else
        {
            for (newcnt = j = 0 ; j<pws[i].numgenes ; j++)
            {
                ll = pwgenes[pws[i].pwgenesindex + j];
                if (ll == USHRT_MAX) // masked out as not in universe
                {
                    fprintf(stderr,"ERROR: %u \n",ll); exit(0);
                    continue;
                }
                *(this_egids+newcnt) = genes[ll].egid;
                newcnt++;
            }
        }
        uptr = (u+used_index);

#if 0
// do we need to sort?
#if RADIX_SORT
#if 0
    // test to make sure it is sorted 
unsigned int prev;
//fprintf(stderr,"testing if sorted\n"); 
        for (j=0 ; j<newcnt ; j++)
        {
           ui = *(this_egids+j) ;
// fprintf(stderr,"%u\n",ui);
           if (j)
           {
               if (ui <= prev) 
               {
                  fprintf(stderr,"%s not sorted %d <= %d \n",pws[i].acc ,ui,prev); 
               }
           }
           prev = ui;
        }
#endif
        if (newcnt < 48)
            radix_ui(this_egids,newcnt);
        else
            qsort(this_egids,newcnt,sizeof(unsigned int),cmp_ui);
#else
        qsort(this_egids,newcnt,sizeof(unsigned int),cmp_ui);
#endif
#endif
        uptr->egids = this_egids;
        uptr->numfixedgenes = newcnt;
        uptr->genehits = (unsigned int *)malloc(sizeof(unsigned int)*newcnt);
        uptr->category = pws[i].category;
        uptr->acc  = strdup(pws[i].acc);       // can't just assign (?)  . because need to free if custom
        uptr->name = strdup(pws[i].name);      // can't just assign
        uptr->numgenes = pws[i].numgenes;
        uptr->pwgenesindex = pws[i].pwgenesindex; // not going to care about this probably. info is now in egids
        uptr->hitcnt = 0;
        used_index++;
    }
// fprintf(stderr,"rpf used_index = %d\n",used_index); fflush(stderr); 

#define MAXUGENES 40000
    if (custom_cnt)
    {
        if (tmpugenes) { free(tmpugenes); tmpugenes = (unsigned int *)0; }
        tmpugenes=malloc(sizeof(unsigned int)*MAXUGENES);
        if (!tmpugenes) { fprintf(stderr,"can't malloc %u bytes \n",(unsigned int)(sizeof(unsigned int)*MAXUGENES));  fflush(stderr); }
        fp = fopen(custom_file,"r");
        if (!fp) { fprintf(stderr,"NOTE: can not open \"%s\" - ignoring\n",custom_file); fflush(stderr); }
        else
        {
            while ( fgets(s, sizeof(s)-2, fp) ) // for each line of custom file
            {
                for (k=0 ; s[k] ; k++) { if ((s[k] == '\n')||(s[k] == '\r')) s[k] = (char)0; } // strip newline
                uptr = (u+used_index); // used_index points to next available slot 
                z = &s[0];
                tmps[0] = custname[0] = (char)0;
                j = tokcnt = k = newcnt = 0;
                keepgoing = 1;
                while (keepgoing == 1)
                {
                    if ( (*(z+j) == '\t') || (*(z+j) == (char)0) )
                    {   // first token is name, 2nd is optional and rest are gene names separated by tabs
                       if (tokcnt == 0) strcpy(custname,tmps); 
                       else if (tokcnt == 1) { /* ignore */ } 
                       else
                       {
                           if (tmps[0])
                           {
                               ui = hugo2egid(tmps);
                               if (ui == (unsigned int)UINT_MAX)
                               {
                                   fprintf(stderr,"Note: invalid gene \"%s\" in custom pathway file.\n",tmps);  
                               }
                               else
                               {
                                   if (universe_cnt) // must be in user universe
                                   {
                                       uiptr = (unsigned int *)bsearch(&ui,user_universe_genes,universe_cnt,sizeof(unsigned int),cmp_ui);
                                       if (uiptr)
                                       {
                                           *(tmpugenes+newcnt) = ui; // put this in the right place
                                           newcnt++;
                                           if (newcnt == MAXUGENES) { fprintf(stderr,"ERROR: too many genes on custom pathway gmt file.\n"); fflush(stderr); }
                                       }
                                   }
                                   else 
                                   {
                                       *(tmpugenes+newcnt) = ui; // put this in the right place
                                       newcnt++;
                                   }
                               }
                           }
                       }
                       tokcnt++;
                       k = 0;
                       tmps[0] = (char)0;
                       if (*(z+j) == (char)0) keepgoing = 0;
                    }
                    else
                    {
                        tmps[k++] = *(z+j); 
                        tmps[k] = (char)0;
                    }
                    j++;
                }
                if (newcnt == 0) continue;
                qsort(tmpugenes,newcnt,sizeof(unsigned int),cmp_ui);
                this_cust_pw_cnt = 0;
                prev = 0;
                this_egids = (unsigned int *)malloc(sizeof(unsigned int)*newcnt);
                for (k=0;k<newcnt;k++) // de-duplicate, put only UNIQUE genes in
                {
                    ui = *(tmpugenes+k); 
                    if (ui != prev)
                    {
                        *(this_egids+this_cust_pw_cnt) = ui;
                        this_cust_pw_cnt++;
                    }
                    prev = ui;
                }
                uptr->acc = strdup(custname);  // remember to free this
                uptr->name = strdup(custname); // remember to free this
                uptr->egids = this_egids;      // remember to free this
                uptr->numgenes = uptr->numfixedgenes = this_cust_pw_cnt;
                uptr->pwgenesindex = USHRT_MAX;
                uptr->genehits = (unsigned int *)malloc(sizeof(unsigned int)*this_cust_pw_cnt); // remember to free this
                uptr->hitcnt = 0;
                uptr->category = CAT_CUSTOM | (string2type("custom") << 28); // "type" and "category" are CUSTOM 
                used_index++;
            }
        }
	if (fp ) { fclose(fp); fp = (FILE *)0; }
        if (tmpugenes) { free(tmpugenes); tmpugenes = (unsigned int *)0; }
    } // if custom cnt

    for (i=0;i<lencust;i++) // add R user custom lists passed as a list of vectors (with "gmt" lines)
    {
        uptr = (u+used_index); // used_index points to next available slot 
        kust = (mycustompw + i);
        uptr->acc = strdup(kust->name);
        uptr->name = strdup(kust->name);
        this_egids = (unsigned int *)malloc(sizeof(unsigned int) * kust->numgenes);
        for (j=0 ; j<kust->numgenes ; j++)   // for each gene in a custom pathway
        {
            *(this_egids+j) = *(kust->genes+j);
// fprintf(stderr,"rpf debug i=%d cust=%d egid=%d\n",i,j,*(this_egids+j)); fflush(stderr);  
        }
        uptr->egids = this_egids; // remember to free this
        uptr->numgenes = uptr->numfixedgenes = kust->numgenes;
        uptr->pwgenesindex = USHRT_MAX;
        uptr->genehits = (unsigned int *)malloc(sizeof(unsigned int)*kust->numgenes);
        uptr->hitcnt = 0;
        uptr->category = CAT_CUSTOM | (string2type("custom") << 28); // "type" and "category" are CUSTOM 
        used_index++;
    }
    *num_used_paths = used_index;
    if (user_universe_genes) free(user_universe_genes);

// fprintf(stderr,"uvcheck before *real_universe_cnt_ptr=%d\n",*real_universe_cnt_ptr); 
    *real_universe = get_used_universe(u, used_index, real_universe_cnt_ptr);
// fprintf(stderr,"uvcheck after *real_universe_cnt_ptr=%d\n",*real_universe_cnt_ptr); 
#if 0
 FILE *fptmp;
 fptmp=fopen("tmpu.chk","w");
 for (i=0;i<*real_universe_cnt_ptr;i++)
 {
 char *zz;
 zz = egid2hugo(*((*real_universe)+i));
 fprintf(fptmp,"%s\n",zz);
 }
 fclose(fptmp);
#endif

#if 0
clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
#define BILLION  1000000000.0;
double time_spent = (time2.tv_sec - time1.tv_sec) + (time2.tv_nsec - time1.tv_nsec) / BILLION;
fprintf(stderr,"Time elapsed is %f seconds\n", time_spent); 
#endif
// fprintf(stderr,"end setup_used_paths\n");
    return u;
}


void print_permute_header(void)
{
// fprintf(stderr,"in print_permute_header()\n"); fflush(NULL); 
    printf("a\t"); // uptr->a 
    printf("c\t"); // uptr->c 
    printf("b\t"); // uptr->b 
    printf("d\t"); // uptr->d 
    printf("countunder\t"); // uptr->countunder 
    printf("countequal\t"); // uptr->countequal 
    printf("countover\t"); // uptr->countover 
    printf("numfixedgenes\t"); // uptr->numfixedgenes 
    printf("OR\t"); // uptr->OR 
    printf("pval\t"); // uptr->pval 
    printf("p_permute_over\t"); // uptr->p_permute_over 
    printf("p_permute_under\t"); // uptr->p_permute_under 
    printf("p_pval\t");
    printf("acc\t"); // uptr->acc 
    printf("tmps\t"); // tmps 
    printf("name"); // uptr->name); 
    printf("\n");
    fflush(stdout);
// fprintf(stderr,"end print_permute_header()\n"); fflush(NULL); 
    return;
} 


#ifdef WEBASSEMBLY
#else

int main(int argc,char *argv[])
{
    char universe_file[PATH_MAX];
    char custom_file[PATH_MAX];
    char tmps[PATH_MAX];
    int user_incnt = 0;
    int no_header_flag = 0;
    int precise_flag = 0;
    int calc_option = 0; // 0=no, 1=nelson 2=rpf 
    int oneside = 0;
    unsigned int seed = 0;
    unsigned int catspat = 0;
    unsigned int num_used_paths = 0;
    struct used_path_type *u = (struct used_path_type *)0;
    struct used_path_type *uptr = (struct used_path_type *)0;
    int status = 0;
    unsigned int *real_universe = (unsigned int *)0;
    unsigned int real_universe_cnt = 0;
    int get_universe_flag = 0;
    unsigned int i = 0;
    unsigned int j = 0;
    char *z = (char *)0;
    double fdr_for_output;
    double d;
    unsigned int num_permutes = PERMUTE_DEFAULT_LOW;
#if NELSON_C
//    double corr_tst;
    FILE *pathcalls = (FILE *)0;
    FILE *pathcalls_tech = (FILE *)0;
#endif
#if NELSON_TEST
    double mean_gpct;
#endif


    srand(time(NULL) + getpid());
    srandom( time( NULL ) );
    (void)setup_by_egids();
    l2p_init_C(argc,argv,&catspat,&precise_flag,&calc_option,&no_header_flag,universe_file,custom_file,&get_universe_flag,&oneside,&seed,&num_permutes);
fprintf(stderr,"rpf calc_option=%d num_permutes=%u\n",calc_option,num_permutes); fflush(NULL);

    // prototype: struct used_path_type *setup_used_paths(unsigned int *num_used_paths, unsigned int catspat, char universe_file[], unsigned int in_universe_cnt,unsigned int *in_universe, char custom_file[], unsigned int *real_universe_cnt_ptr,unsigned int **real_universe,unsigned int lencust,struct custom_type *mycustompw);
    u = setup_used_paths(&num_used_paths, catspat,universe_file, 0,(void *)0,custom_file,&real_universe_cnt,&real_universe,0,(struct custom_type *)0);
    if (get_universe_flag == 1)
    {
        print_universe(real_universe_cnt,real_universe);
        return 0;
    }

fprintf(stderr,"in main(), before  l2pfunc num_permutes=%u\n",num_permutes); fflush(NULL);
    status = l2pfunc(u,num_used_paths,real_universe_cnt,real_universe,calc_option,&user_incnt,oneside,num_permutes);
// fprintf(stderr,"in main(), after  l2pfunc\n"); fflush(NULL);

#if NELSON_TEST
    pathcalls = fopen(OUTPUTFILE, "w");
    pathcalls_tech = fopen(TECHOUTPUTFILE, "w");

    fprintf(pathcalls, "GPCC_pval\tGPCC_fdr\tFET_pval\tFET_fdr\tGPCC_OR\tFET_OR\n");
    fprintf(pathcalls_tech, "a\tb\tc\td\ta\tb\tc\td\tcountunder\tcountequal\tcountover\tpwgenect\tpw_gpct\tOR\tgpcc_OR\tFET_p\tpermute_p\tmean_gpct\tmean_deg_gpct\tmean_u_gpcts\tgp_corr_p\tacc\tname\n");
#endif

    if (calc_option == CALC_OPTION_GPCC) // gpcc
    {
        if (no_header_flag == 0) print_permute_header();
        for (i=0 ; i<num_used_paths ; i++)
        {
            uptr = (u+i);
            if ( (uptr->a == 0) && (uptr->b == 0) ) d   = (double)0; else d = (double)((double)uptr->hitcnt / (double)(uptr->numfixedgenes));
            categories_pattern_to_strings(uptr->category,tmps);
            fdr_for_output = uptr->fdr;
            printf("%s\t%20.18f\t%20.18f\t%20.18f\ti%d\t%20.18f\t%s\t%s\t%d\t",uptr->name,uptr->enrichment_score,
                uptr->gpcc_p, uptr->gpcc_fdr, 
                uptr->hitcnt, d,
                uptr->acc,
                tmps, // uptr->category ,
                uptr->hitcnt); 
#if NELSON_TEST
// fprintf(stderr,"rpf here 3 %d in CALC_OPTION_GPCC\n",i); fflush(NULL); 
            printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%20.18f\t%11.2e\t%11.2e\t%11.2e\t%11.2e\t%s\t%s\t%s\t",
                          uptr->a, uptr->c, uptr->b, uptr->d,
                          uptr->countunder, uptr->countequal, uptr->countover, uptr->numfixedgenes, uptr->OR,
                          uptr->pval,  uptr->p_permute_over, uptr->p_permute_under, uptr->permute_pval,
                          uptr->acc, tmps, uptr->name);
            mean_gpct = (float) uptr->pathhits_gpsum/(float) uptr->hitcnt;
            categories_pattern_to_strings(uptr->category,tmps);
            fprintf(pathcalls, "%11.2e\t%11.2e\t%11.2e\t%11.2e\t%11.2f\t%11.2f\n",
                    uptr->gpcc_p, uptr->gpcc_fdr, uptr->pval, uptr->fdr, uptr->gpcc_OR, uptr->OR); // i, randhits diagnostic

            fprintf(pathcalls_tech, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%11.2f\t%11.2f\t%11.2e\t%11.2e\t%11.2f\t%11.2f\t%11.2f\t%11.2e\t%s\t%s\n",
                    uptr->a, uptr->b, uptr->c, uptr->d,
                    uptr->A_scaled, uptr->B_scaled, uptr->C_scaled, uptr->D_scaled,
                    uptr->countunder, uptr->countequal, uptr->countover, uptr->numfixedgenes, uptr->pathcountsum,
                    uptr->OR, uptr->gpcc_OR, uptr->pval, uptr->p_permute_over, mean_gpct,  mean_deg_gpcount, mean_u_gpcount,
                    uptr->gpcc_p, uptr->acc,  uptr->name);
#endif
// fprintf(stderr,"hitcnt = %d\n",uptr->hitcnt); 
            for (j=0 ; j<uptr->hitcnt ; j++)
            {
                z = egid2hugo(*((uptr->genehits)+j));
                printf("%s ",z);
            }
            printf("\n");
        }
    }
    else if ( (calc_option == CALC_OPTION_GPCC2) || (calc_option == CALC_OPTION_GPCC3) )
    {
        if (no_header_flag == 0) print_header();
        for (i=0 ; i<num_used_paths ; i++)
        {
           uptr = (u+i);
           uptr->enrichment_score = ( ((double)uptr->hitcnt / (double)uptr->numfixedgenes) - ((double)user_incnt/(double)real_universe_cnt) );
           uptr->enrichment_score *= 100.0;
           categories_pattern_to_strings(uptr->category,tmps);
           fdr_for_output = uptr->fdr;
           if (precise_flag)
           {
               printf("%20.18f\t%20.18f\t%11.9f\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t",
                   uptr->gpcc_p,     fdr_for_output, uptr->enrichment_score,
                   uptr->a, uptr->b, uptr->c, uptr->d,
                   uptr->acc, tmps,
                   uptr->name); // type2string(uptr->category >> 28))
           }
           else
           {
               printf("%11.9f\t%11.9f\t%11.9f\t%11.9f\t%11.7f\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t",
                   uptr->gpcc_p,     
                   fdr_for_output, uptr->pval4,uptr->fdr4,uptr->enrichment_score,
                   uptr->a, uptr->b, uptr->c, uptr->d,
                   uptr->acc, tmps,
                   uptr->name); // type2string(uptr->category >> 28))
#if 0
               printf("%11.9f\t%11.9f\t%11.7f\t%d\t%d\t%d\t%d\t%s\t%s\t%s\n",
                   uptr->pval,     fdr_for_output, uptr->enrichment_score,
                   uptr->a, uptr->b, uptr->c, uptr->d,
                   uptr->acc, tmps,
                   uptr->name, );
           printf("%11.9f\t%11.9f\t%11.7f\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n",
               uptr->pval,     fdr_for_output, uptr->enrichment_score,
               uptr->a, uptr->b, uptr->c, uptr->d,
               uptr->acc, tmps,
               uptr->name, type2string(uptr->category >> 28));
#endif
            }
            for (j=0 ; j<uptr->hitcnt ; j++)
            {
                z = egid2hugo(*((uptr->genehits)+j));
                printf("%s ",z);
            }
            printf("\n");
        }
    }
    else if ( calc_option == CALC_OPTION_PERMUTE)
    {
        if (no_header_flag == 0) print_header();
        for (i=0 ; i<num_used_paths ; i++)
        {
           uptr = (u+i);
           categories_pattern_to_strings(uptr->category,tmps);
           fdr_for_output = uptr->fdr;
           if (precise_flag)
           {
               printf("%20.18f\t%20.18f\t%11.9f\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t",
                   uptr->pval4,     fdr_for_output, uptr->enrichment_score,
                   uptr->a, uptr->b, uptr->c, uptr->d,
                   uptr->acc, tmps,
                   uptr->name, type2string(uptr->category >> 28));
           }
           else
           {
               printf("%11.9f\t%11.9f\t%11.9f\t%11.9f\t%11.7f\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t",
                   uptr->pval4,     
                   fdr_for_output, uptr->pval4,uptr->fdr4,uptr->enrichment_score,
                   uptr->a, uptr->b, uptr->c, uptr->d,
                   uptr->acc, tmps,
                   uptr->name, type2string(uptr->category >> 28));
#if 0
               printf("%11.9f\t%11.9f\t%11.7f\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n",
                   uptr->pval,     fdr_for_output, uptr->enrichment_score,
                   uptr->a, uptr->b, uptr->c, uptr->d,
                   uptr->acc, tmps,
                   uptr->name, type2string(uptr->category >> 28));
           printf("%11.9f\t%11.9f\t%11.7f\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n",
               uptr->pval,     fdr_for_output, uptr->enrichment_score,
               uptr->a, uptr->b, uptr->c, uptr->d,
               uptr->acc, tmps,
               uptr->name, type2string(uptr->category >> 28));
#endif
            }
            for (j=0 ; j<uptr->hitcnt ; j++)
            {
                z = egid2hugo(*((uptr->genehits)+j));
                printf("%s ",z);
            }
            printf("\n");
        }
    }
    else if ( calc_option == CALC_OPTION_PERMUTE3)
    {
        if (no_header_flag == 0) print_header();
        for (i=0 ; i<num_used_paths ; i++)
        {
           uptr = (u+i);
           uptr->a = uptr->hitcnt;
           uptr->b = uptr->numfixedgenes - uptr->hitcnt;
           uptr->c = user_incnt;
           uptr->d = real_universe_cnt - user_incnt;
           uptr->enrichment_score = ( ((double)uptr->hitcnt / (double)uptr->numfixedgenes) - ((double)user_incnt/(double)real_universe_cnt) );
           uptr->enrichment_score *= 100.0;
#if 0
   // might not be right , check this
a=hits
b=thispw-hits
c=list
d=u-list
        a2 = localhitcnt;
        b2 = uptr->numfixedgenes - localhitcnt;
        c2 = incnt;
        d2 = real_universe_cnt - incnt;
#endif

           categories_pattern_to_strings(uptr->category,tmps);
           fdr_for_output = uptr->fdr;
           if (precise_flag)
           {
// xxx aaa
               printf("%20.18f\t%20.18f\t%11.9f\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t",
                   uptr->pval4,     fdr_for_output, uptr->enrichment_score,
                   uptr->a, uptr->b, uptr->c, uptr->d,
                   uptr->acc, tmps,
                   uptr->name); 
           }
           else
           {
               printf("%11.9f\t%11.9f\t%11.9f\t%11.9f\t%11.7f\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t",
                   uptr->pval4,     
                   fdr_for_output, uptr->pval4,uptr->fdr4,uptr->enrichment_score,
                   uptr->a, uptr->b, uptr->c, uptr->d,
                   uptr->acc, tmps,
                   uptr->name); //  type2string(uptr->category >> 28))
#if 0
               printf("%11.9f\t%11.9f\t%11.7f\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n",
                   uptr->pval,     fdr_for_output, uptr->enrichment_score,
                   uptr->a, uptr->b, uptr->c, uptr->d,
                   uptr->acc, tmps,
                   uptr->name, type2string(uptr->category >> 28));
           printf("%11.9f\t%11.9f\t%11.7f\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n",
               uptr->pval,     fdr_for_output, uptr->enrichment_score,
               uptr->a, uptr->b, uptr->c, uptr->d,
               uptr->acc, tmps,
               uptr->name, type2string(uptr->category >> 28));
#endif
            }
            for (j=0 ; j<uptr->hitcnt ; j++)
            {
                z = egid2hugo(*((uptr->genehits)+j));
                printf("%s ",z);
            }
            printf("\n");
        }
    }
    else
    {
        if (no_header_flag == 0) print_header();
        for (i=0 ; i<num_used_paths ; i++)
        {
           uptr = (u+i);
           categories_pattern_to_strings(uptr->category,tmps);
           fdr_for_output = uptr->fdr;
           if (precise_flag)
           {
               printf("%20.18f\t%20.18f\t%11.9f\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t",
                   uptr->pval,     fdr_for_output, uptr->enrichment_score,
                   uptr->a, uptr->b, uptr->c, uptr->d,
                   uptr->acc, tmps,
                   uptr->name, type2string(uptr->category >> 28));
           }
           else
           {
               printf("%11.9f\t%11.9f\t%11.7f\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t",
                   uptr->pval,     
                   fdr_for_output, uptr->enrichment_score,
                   uptr->a, uptr->b, uptr->c, uptr->d,
                   uptr->acc, tmps,
                   uptr->name, type2string(uptr->category >> 28));
#if 0
               printf("%11.9f\t%11.9f\t%11.7f\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n",
                   uptr->pval,     fdr_for_output, uptr->enrichment_score,
                   uptr->a, uptr->b, uptr->c, uptr->d,
                   uptr->acc, tmps,
                   uptr->name, type2string(uptr->category >> 28));
           printf("%11.9f\t%11.9f\t%11.7f\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n",
               uptr->pval,     fdr_for_output, uptr->enrichment_score,
               uptr->a, uptr->b, uptr->c, uptr->d,
               uptr->acc, tmps,
               uptr->name, type2string(uptr->category >> 28));
#endif
            }
            for (j=0 ; j<uptr->hitcnt ; j++)
            {
                z = egid2hugo(*((uptr->genehits)+j));
                printf("%s ",z);
            }
            printf("\n");
        }
    }
// fprintf(stderr,"rpf near end no_header_flag=%d\n",no_header_flag); fflush(NULL); 
    
    fflush(NULL); 

    if (u)
    {
        for (i=0 ; i<num_used_paths ; i++)
        {
           uptr = (u+i);
           if (uptr->genehits) { free(uptr->genehits); uptr->genehits = (unsigned int *)0; }
           if (uptr->egids)    { free(uptr->egids);    uptr->egids = (unsigned int *)0; }
           if (uptr->acc)      { free(uptr->acc);      uptr->acc = (char *)0; }
           if (uptr->name)     { free(uptr->name);     uptr->name = (char *)0; }
        }
        free(u);
    }
    if (by_egids) free(by_egids);
    if (real_universe) free(real_universe);


#if NELSON_C
    if (pathcalls) fclose(pathcalls);
    if (pathcalls_tech) fclose(pathcalls_tech);
#endif

// fprintf(stderr,"rpf near end of main()\n"); fflush(NULL);

#if 0
       // next three lines are for valgrind
    fclose(stdin);
    fclose(stdout);
    fclose(stderr);
#endif

    return status;
}
#endif


/*

Hi Rich,

Create a vector of gene names, with each name appearing as often as the number of pathways the gene appears in.  
From this list randomly select as many genes as on your DEG list.
Replace repeats, need a list of unique genes
                                                             ii.     Create the contingency table for each pathway of interest (list of top pathways)

                                                            iii.     Obtain and record its p value

                                                            iv.     For each pathway of interest keep count of how many of the 200 permutations produced more significant p values than the real data

Repeat a) etc. 200 times
The counts/200 are approximate p values for the pathways.

REFERENCE NOTES on L2P:


1.Fisher’s Exact Test is as below
 
Odds ratio = (A*D)/(B*C)

A and D are high where there’s true significance in the pathway (a good pathway).

2. Adjusted p-value (GPCC): gene pathway count correction): another column for L2P


Calculate odds ratio for each pathway….    ->    Get all Odds ratios > 1

Select from a list of “augmented set of random genes” – represent same proportions as in the total pathways  

Recreate Fisher based on 500 “random genes” -  calculate the odds ratio 1000x or more

Is the Odds Ratio the same or more extreme in the same direction.

                Count these instances where (Odds Ratio)rand  > =  (Odds Ratio)real and divide by the number of permutations (ie. random tests).

The p value is the permutation p value

3. Multiple Comparison Correction (BH pvalue): another column for L2P

Rich’s own implementation of BH in C

4. Zeroes in the L2P output


?formatC
?apply
<0.001  -> 1x10-4
0.05 <- 5x10-2

colNumber new_name (l2p output) current_MAAPster_name old_name (l2p) description

1 Pathway_Name Description pathwayname Name of pathway 
2 Category Source category KEGG,REACTOME,GO,PANT(=PANTHER),PID=(pathway interaction database) *was soruce* 
3 P_Value P_Value pval fisher's 2x2 exact test 
4 FDR FDR fdr false discovery rate: default is benjamini hochberg, 
**5 Permutation_p-value Permute_pval See detailed methods 
6 Enrichment_Score Ratio ratio (A/B)/(C/D) * 100 
**7 Percent_Gene_Hits_Per_Pathway NEW FIELD (A/(A+B)) * 100 
8 Significant_Genes_IN_Pathway Number_Hits pwhitcount (A) 
9 Non-Significant_genes_IN_Pathway Number_Misses pwnohitcount (B) 
10 Significant_genes_NOT_IN_Pathway Number_User_Genes inputcount (C) 
11 Non-Significant_Genes_NOT_IN_Pathway Total_Genes_Minus_Input pwuniverseminuslist (D) 
12 Pathway_ID Pathway_ID pathwayaccessionidentifier canonical accession ( if available, otherwise assigned by us ) 
13 Pathway_Type pathwaytype functional_set, pathway, structural_complex , custom 
14 Gene_List Gene_List genes_space_separated HUGO genes from user list that hit the pathway

Table referenced from attached email:
Significant Genes Non-Significant Genes
IN Pathway A B
NOT IN Pathway C D
1 Pathway_name        Name of pathway
2 Enrichment_score    ((number_hits /(number_hits+number_misses)) - (number_user_genes/(number_user_genes+total_gens_minus_input))) * 100" <- tie-breaker #2
3 GPCC_pval           Gpcc (gene pathway count correction) pval <- sort on this first (#1)
4 GPCC_FDR            False discovery rate: benjamini hochberg (of GPCC p-value) <- this is what we recommend users use
5 FET_pval            Fisher's exact p-value (legacy)
6 FET_FDR             False discovery rate (FET)
7 Number_hits         Number of genes hit in pathway
8 Number_genes_in_pathway        Number of genes in the pathway (size of pathway)*
9 Percent_hits        Number_hits/total genes in pathway* <- tie-breaker #3
10 Pathway_id         Canonical accession ( if available, otherwise assigned by us )
11 Category           KEGG,REACTOME,GO,PANT(=PANTHER),PID=(pathway interaction database)  *was "source"*
12 Pathway_type       Functional_set,pathway,structural_complex,custom
13 Genes_in_pathway   HUGO genes from user that hit the pathway

VERSION 8 L2P:

Columns that will be output with L2P:

1 Pathway_name          Name of pathway
2 Enrichment_score      ((number_hits /(number_hits+number_misses)) - (number_user_genes/(number_user_genes+total_gens_minus_input))) * 100" <- tie-breaker #2
3 GPCC_pval             Gpcc (gene pathway count correction) pval <- sort on this first (#1)
4 GPCC_FDR              False discovery rate: benjamini hochberg (of GPCC p-value) <- this is what we recommend users use
5 FET_pval              Fisher's exact p-value (legacy)
6 FET_FDR               False discovery rate (FET)
7 Number_hits           Number of genes hit in pathway
8 Number_genes_in_pathway        Number of genes in the pathway (size of pathway)*
9 Percent_hits          Number_hits/total genes in pathway* <- tie-breaker #3
10 Pathway_id           Canonical accession ( if available, otherwise assigned by us )
11 Category             KEGG,REACTOME,GO,PANT(=PANTHER),PID=(pathway interaction database)  *was "source"*
12 Pathway_type         Functional_set,pathway,structural_complex,custom
13 Genes_in_pathway     HUGO genes from user that hit the pathway

Optional Columns: This is for Rich to put into a flag where optional columns can be added to the L2P output (Optionalcols = 1, by default Optionalcols = 0)

A – universe – (DEG count + pathway count) 
B – pathway count – number of hits
C – total DEG count - number of hits
D – number of hits
Odds_ratio              odds ratio
Sum_of_pathways         Sum of unique pathways where pathway genes are also found
*/

