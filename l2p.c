
/*
vim l2p.c ; gcc -Ofast -Wall -o l2p pathcommon.c l2p.c -lm  ; strip l2p
vim l2p.c ; gcc -O3 -march=native -flto -Wall -o l2p l2p.c -lm
vim l2p.c ; gcc -D_FORTIFY_SOURCE=2 -fstack-protector --param ssp-buffer-size=4 -fPIE -pie -Wl,-z,relro,-z,now -o l2p l2p.c -lm
vim l2p.c ; gcc -Wall -pg l2p.c -o l2p
#vim l2p.c ; gcc -D_FORTIFY_SOURCE=2 -fstack-protector --param ssp-buffer-size=4 -fPIE -pie -Wl,-z,relro,-z,now (ld -z relro and ld -z now) -o l2p l2p.c -lm

todo:
1) fix output
2) ensembl to

output:
     1  pval
     2  fdr
     3  ratio                      if postive, genes are OVER REPRESENTED, if negative genes are UNDER REPRESENTED
     4  pathwayhitcount            number of genes hit in pathway
     5  numberofgenesin pathway    number of genes in the pathway
     6  inputnumberofgenes         total count of user genes (user input)
     7  genesinpathwaysuniverse    total number of unique genes in all pathways
     8  pathwayaccessionidentifier canonical accession ( if availible, otherwise assigned by us )
     9  category                   KEGG,REACTOME,GO,PANT(=PANTHER),PID=(pathway interaction database)  *was "soruce"*
    10  pathwayname                Name of pathway
    11  pathwaytype genes_space_separated   HUGO genes from user that hit the pathway

    run with gdb

% gdb myprogram
gdb> run params ... < input.txt

    printf "TP53\nPTEN\nAPC\nKRAS\nNRAS\n"  > jokein
    printf "TP53\nPTEN\nAPC\nKRAS\nNRAS\n"  > jokein
    
    gdb ./l2p
    run  -categories=PID,KEGG < jokein
*/

#ifdef L2P_USING_R
#include <R.h>
#include <Rdefines.h>
#endif

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>
#include <errno.h>
#include <unistd.h>

#include "pathworks.h"

static int catspat = 0;   // categories pattern   , if 0, use all
static int no_header_flag = 0;
static int precise_flag = 0; // print out more precise 
static int user_universe_flag = 0;

unsigned char *ucz = (void *)0; // space
static struct binpathouttype binpath[MAXBSID];
static int numbinpaths = 0;
struct used_path_type
{
    double pval;
    double fdr;
    double ad; // ratio
    int index; // index to binpath
};

static struct used_path_type usedpaths[MAXBSID];
static int numusedpaths;

static struct bingenetype bingene[MAXGENE];   // int geneid; char hugo[]; char ensembl[]; int pathcount; int pathplace;
static int numbingenes;
static long int spacesize;

static char fn_pathbin[PATH_MAX];
static char fn_genesbin[PATH_MAX];
static char fn_spacebin[PATH_MAX];
static char universe_file[PATH_MAX];
static char exedir[PATH_MAX]; // executable directory 
static int numg = 0; // number of genes in universe  


// -- for benjamini hochberg FDR ...
struct ordertype
{
    double val;
    int order;
};

static int cmp_ordertype_by_val_REV(const void *a, const void *b)
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

#if 0
not used
static int cmp_ordertype_by_val(const void *a, const void *b)
{
    struct ordertype *aa;
    struct ordertype *bb;
    aa = (void *)a;
    bb = (void *)b;
    if      (aa->val >  bb->val) return 1;
    else if (aa->val <  bb->val) return -1;
    // if (aa>bb) return 1; else if (aa<bb) return -1;
    return 0;
}
static int cmp_ordertype_by_order_REV(const void *a, const void *b)
{
    struct ordertype *aa;
    struct ordertype *bb;
    aa = (void *)a;
    bb = (void *)b;
    if      (aa->order <  bb->order) return 1;
    else if (aa->order >  bb->order) return -1;
    return 0;
}
#endif

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


#define RDEBUG 0

static void benjaminihochberg(int n,double pvals[], double returnpvals[])
{
#if 0
// here's the code from R that I re-imagined
                   i <- lp:1L
                   o <- order(p, decreasing = TRUE)
                   ro <- order(o)
                   pmin(1, cummin( n / i * p[o] ))[ro]
#endif
    int j,k;
    struct ordertype *i;
    struct ordertype *o;
    struct ordertype *po;
    struct ordertype *cummin;
//    struct ordertype *ro;
//    struct ordertype *intermed;

    i = (struct ordertype *)malloc((sizeof(struct ordertype))*n);
    for (k=n,j=0;j<n;j++,k--) (i+j)->order=k;

#if RDEBUG
FILE *fp;
fp = fopen("test.pvals","w");
#endif
    o = (struct ordertype *)malloc((sizeof(struct ordertype))*n);
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
    po = (struct ordertype *)malloc((sizeof(struct ordertype))*n);
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


static double lngamm(double z)
// Reference: "Lanczos, C. 'A precision approximation
// of the gamma double ', J. SIAM Numer. Anal., B, 1, 86-96, 1964."
// Translation of  Alan Miller's FORTRAN-implementation
// See http://lib.stat.cmu.edu/apstat/245
{
  double x = 0;
  x += 0.1659470187408462e-06/(z+7);
  x += 0.9934937113930748e-05/(z+6);
  x -= 0.1385710331296526    /(z+5);
  x += 12.50734324009056     /(z+4);
  x -= 176.6150291498386     /(z+3);
  x += 771.3234287757674     /(z+2);
  x -= 1259.139216722289     /(z+1);
  x += 676.5203681218835     /(z);
  x += 0.9999999999995183;
//   return(Math.log(x)-5.58106146679532777-z+(z-0.5)*Math.log(z+6.5));
   return(log(x)-5.58106146679532777-z+(z-0.5)*log(z+6.5));
}

static double lnfact(int n)
{
  if(n<=1) return(0);
  return(lngamm(n+1));
}

static double lnbico(int n,int k)
{
  double ret;
  ret =lnfact(n)-lnfact(k)-lnfact(n-k);

  return(ret);
}

static double hyper_323(int n11,int n1_,int n_1,int n)
{
    double d;
    d = lnbico(n1_,n11) + lnbico(n-n1_,n_1-n11) - lnbico(n,n_1);
    d = exp(d);
    return d;
}

static int sn11,sn1_,sn_1,sn;
static double sprob;

static double hyper0(int n11i,int n1_i,int n_1i,int ni)
{

// printf("in hyper0 %d %d %d %d \n",n11i,n1_i,n_1i,ni);
  if(!(n1_i|n_1i|ni))
  {
//printf("in hyper0 NOT %d %d %d %d \n",n11i,n1_i,n_1i,ni);
    if(!(n11i % 10 == 0))
    {
      if(n11i==sn11+1)
      {
        sprob *= ((double)(sn1_-sn11)/(double)(n11i))*((double)(sn_1-sn11)/(double)(n11i+sn-sn1_-sn_1));
        sn11 = n11i;
        return sprob;
      }
      if(n11i==sn11-1)
      {
        sprob *= ((double)(sn11)/(double)(sn1_-n11i))*((double)(sn11+sn-sn1_-sn_1)/(double)(sn_1-n11i));
        sn11 = n11i;
        return sprob;
      }
    }
    sn11 = n11i;
  }
  else
  {
//printf("in hyper0 else %d %d %d %d \n",n11i,n1_i,n_1i,ni);
    sn11 = n11i;
    sn1_=n1_i;
    sn_1=n_1i;
    sn=ni;
  }
// printf("in hyper0 before hyper_323 %d %d %d %d\n",sn11,sn1_,sn_1,sn);
  sprob = hyper_323(sn11,sn1_,sn_1,sn);
// printf("hyper returns sprob = %10.15f after hyper_323\n",sprob);
  return sprob;
}

static double  hyper(int n11)
{
  return(hyper0(n11,0,0,0));
}

static double sleft,sright,sless,slarg;

static double exact(int n11,int n1_,int n_1,int n)
{
  int i,j;
  double p;
  double prob;
  int max=n1_;

// printf("in exact %d %d %d %d \n",n11,n1_, n_1,n);
  if(n_1<max) max=n_1;
  int min = n1_+n_1-n;
  if(min<0) min=0;
  if(min==max)
  {
    sless = 1;
    sright= 1;
    sleft = 1;
    slarg = 1;
    return 1;
  }

  prob=hyper0(n11,n1_,n_1,n);
// printf("in exact prob=%20.17f \n",prob);
  sleft=0;
  p=hyper(min);
  for(i=min+1; p<0.99999999*prob; i++)
  {
    sleft += p;
    p=hyper(i);
  }
  i--;
  if(p<1.00000001*prob) sleft += p;
  else i--;
  sright=0;
  p=hyper(max);
  for(j=max-1; p<0.99999999*prob; j--)
  {
    sright += p;
    p=hyper(j);
  }
  j++;
  if(p<1.00000001*prob) sright += p;
  else j++;
  if(abs(i-n11)<abs(j-n11))
  {
    sless = sleft;
    slarg = 1 - sleft + prob;
  }
  else
  {
    sless = 1 - sright + prob;
    slarg = sright;
  }
  return prob;
}


static double left,right,twotail;
//static int n11old=-1;
//static int n12old=-1;
//static int n21old=-1;
//static int n22old=-1;

static double exact22(int n11_,int n12_,int n21_,int n22_)
{
#if 0
  double prob = 0.0;
  double n11_ = parseInt("0"+n11,10);
  double n12_ = parseInt("0"+n12,10);
  double n21_ = parseInt("0"+n21,10);
  double n22_ = parseInt("0"+n22,10);
   if((n11old==n11_i) && (n12old==n12_) && (n21old==n21_) && (n22old==n22_)) return;
  n11old=n11_;
  n12old=n12_;
  n21old=n21_;
  n22old=n22_;
#endif
  if(n11_<0) n11_ *= -1;
  if(n12_<0) n12_ *= -1;
  if(n21_<0) n21_ *= -1;
  if(n22_<0) n22_ *= -1;

  int n1_ = n11_+n12_;
  int n_1 = n11_+n21_;
  int n   = n11_ +n12_ +n21_ +n22_;

  // prob = exact(n11_,n1_,n_1,n); don't need return value
  ( void ) exact(n11_,n1_,n_1,n);
// printf("prob after exact is  %30.25f\n",prob);

  left    = sless;
  right   = slarg;
  twotail = sleft+sright;
// printf("left=%20.15f right=%20.15f \n",sleft,sright);
  if(twotail>1) twotail=1;

  return twotail;
// printf("%d %d %d %d %12.8f prob=%20.15f twotail=%20.15f\n",(int)n11_, (int)n12_, (int)n21_, (int)n22_, twotail,prob,twotail);
/*
  document.form1.output.value +=
  newline+
  " TABLE = [ " +
  n11_+" , "+
  n12_+" , "+
  n21_+" , "+
  n22_+" ]" + newline +
  "Left   : p-value = "+ left + newline +
  "Right  : p-value = "+ right + newline +
  "2-Tail : p-value = "+ twotail +
  newline +   "------------------------------------------";
*/
}


#if 0
test code
int fish_exact_test()
{
    double d;
// TEST fisher's exact code
    int n11_;double n12_;double n21_;double n22_;

    n11_ = 3;
    n12_ = 40;
    n21_ = 297;
    n22_ = 29960; 

    d = exact22(n11_,n12_,n21_,n22_);
// 3 hits in p53 pathway which has  40 genes ,total 29960 genes , user input 300 genes  
    d = exact22(3,40,297,29960);
    return 0;
}
#endif


static struct hugo_type *hugos;  // parallel to bingenes, defined in .h file
static int numhugos;

static int cmp_hugo(const void *a, const void *b)
{
    return strcmp(((struct hugo_type *)a)->hugo, ((struct hugo_type *)b)->hugo);;
}

static int cmp_bingene(const void *a, const void *b)
{
    if      ( ((struct bingenetype *)a)->geneid < ((struct bingenetype *)b)->geneid ) return -1;
    else if ( ((struct bingenetype *)a)->geneid > ((struct bingenetype *)b)->geneid ) return 1;
    return 0;
}


// so we can search by hugo
static int setup_hugo_parallel_to_genebin(int n)
{
    int i;
    size_t sz;

    sz = ((size_t)(n)*sizeof(struct hugo_type));
    hugos = (struct hugo_type *)malloc(sz);
    if (!hugos) { fprintf(stderr,"ERROR: can't malloc in setup_hugo_parallel_to_genebin...()\n"); return -5; }
    for (i=0;i<n;i++)
    {
        (hugos+i)->hugo = strdup(bingene[i].hugo);
        (hugos+i)->generec_ptr = &bingene[i]; // from struct bingenetype bingene[MAXGENE];
        // (hugos+i)->generec_ptr = (struct bingenetype *)&bingene[i]; // from struct bingenetype bingene[MAXGENE];
    }
    numhugos = n;

    qsort(hugos,n,sizeof(struct hugo_type),cmp_hugo); // rearrange by hugo name for bsearch , keep link to original genebin rec 
#if 0
    for (i=0;i<n;i++)
    {
         fprintf(stderr,"cmp %s %s %s\n",(hugos+i)->hugo,(hugos+i)->generec_ptr->hugo,bingene[i].hugo);
    }
exit(0);
#endif

    return 0;
}

static void free_hugos(void)
{
    int i;
    if (hugos)
    {
        for (i=0;i<numhugos;i++)
        {
            if ((hugos+i)->hugo) free((hugos+i)->hugo);
        }
        free(hugos);
    }
    hugos = (void *)0;
    numhugos = 0;
    return;
}

static void free_binpaths(void)
{
    int i;
    for (i=0;i<numbinpaths;i++)
    {
        if (binpath[i].hits)
        {
            free(binpath[i].hits); 
            binpath[i].hits = (void *)0;
        }
    }
    numbinpaths = 0;
    return;
}



static char *type2string(int type)
{
   static char functional_set[] = "functional_set";
   static char pathway  [] = "pathway";  
   static char structural_complex  [] = "structural_complex";  
   static char null_info  [] = "NA";  

   if (type == 1) return &functional_set[0]; 
   else if (type == 2) return &pathway[0];
   else if (type == 3) return &structural_complex[0];
   return &null_info[0];
}


static void do_pvals_and_bh(int ng,int ingenecnt)
{
    struct binpathouttype *binpathptr;
    int skip = 1;
    double d;
    double ad;
    int i;
    double *pvals;
    double *fdrs;
    int hitcnt ;

    pvals = (double *)malloc((size_t)(sizeof (double)*numbinpaths)); 
    if (!pvals) { fprintf(stderr,"ERROR: can't malloc in do_pvals_and_bh() 1\n"); return; }

    for (numusedpaths=i=0 ; i<numbinpaths ; i++)
    {
        hitcnt = (unsigned int)0;
        binpathptr = &binpath[i];
        if (binpathptr->category & catspat) skip = 0; 
        else                                skip = 1; 
        if (skip == 1) continue;

        if (binpathptr->hits)
        {
            hitcnt = *(unsigned int *)(binpathptr->hits);
            if (hitcnt > (unsigned long int)binpathptr->numgenes) 
            {
                 fprintf(stderr,"ERROR: more hits than genes (hitcnt %d > %d for %s)\n", hitcnt , (binpath+i)->numgenes,(char *)(ucz+((binpathptr)->accession ) )) ;
                 fprintf(stderr,"details: i=%d , current hits=%u\n", i,hitcnt); 
                 fflush(stderr);
                 exit(0);
            }
        }
// test ...
//d = exact22(n11_,n12_,n21_,n22_); note:3 hits in p53 pathway which has  40 genes ,total 29960 genes , user input 300 genes 
// d = exact22(3,40,297,29960);
        d = exact22((int)hitcnt,(binpathptr)->numgenes,ingenecnt,ng);

        ad = ( ((double)hitcnt/ (double)(binpathptr)->numgenes) - ((double)ingenecnt/(double)ng) );
        *(pvals+numusedpaths) = usedpaths[numusedpaths].pval = d;
        usedpaths[numusedpaths].ad = ad;
        usedpaths[numusedpaths].index = i;
        numusedpaths++;
    }

    fdrs = (double *)malloc((size_t)(sizeof (double)*numbinpaths)); 
    if (!fdrs) { free(pvals); /* clean up */ fprintf(stderr,"ERROR: can't malloc in do_pvals_and_bh() 2\n"); exit(0); }
    benjaminihochberg(numusedpaths,pvals,fdrs);
    for (i=0 ; i<numusedpaths ; i++)
         usedpaths[i].fdr = *(fdrs+i);
    free(pvals); // path parallel
    free(fdrs);
    return;
}


#ifdef L2P_USING_R

static int R_deal_with_universe_list(SEXP lst, int numbingenes, int numbinpaths, 
                 int *numg_ptr, unsigned char *spaceptr) 
{
    char tmphugo[1024];
    struct hugo_type *hugoptr;
    struct hugo_type h;
    int found = 0;
    int geneid;
    int *iptr;
    int o2g;
    int i,j;
    int n = 0;
    int new_fix_gene_count;
    struct bingenetype Xgenerec; 
    struct bingenetype *Xgenerec_ptr; 

    n = length(lst);

//     fprintf(stderr,"rpf in R_deal_with_universe_list\n"); fflush(NULL); 
    for (i=0;i<numbingenes;i++) hugos[i].status = 0;  // assume "guilty"
    for (i=0;i<n;i++) 
    {
        memset(&h,0,sizeof(h));
        h.hugo = strdup(CHAR(STRING_ELT(lst, i)));
        hugoptr = bsearch(&h,hugos,numbingenes,sizeof(struct hugo_type),cmp_hugo);
        if (hugoptr) hugoptr->status = 1; // "innocent"
        free(h.hugo);
    }
    for (i=0 ; i<numbinpaths ; i++)         // for each pathway
    {
       new_fix_gene_count = binpath[i].numgenes;
       o2g = binpath[i].offset2geneids;
       if (o2g<=0) continue;
       iptr = (int *)(spaceptr + o2g);
       for (j=0 ; j<binpath[i].numgenes ; j++)  // for each gene for this pathway
       {
           geneid = *(iptr+j);
           found = 0;
           // memset(&Xgenerec,0,sizeof(struct bingenetype)); not needed
           Xgenerec.geneid = geneid;
           Xgenerec_ptr = bsearch(&Xgenerec,&bingene[0],numbingenes,sizeof(struct bingenetype),cmp_bingene);
           if (Xgenerec_ptr)
           {
               strcpy(tmphugo, Xgenerec_ptr->hugo);
               h.hugo = &tmphugo[0];
               hugoptr = bsearch(&h,hugos,numbingenes,sizeof(struct hugo_type),cmp_hugo);
               if (hugoptr)
               {
                   found = 1;
                   if (hugoptr->status == 0) // not in universe
                   {
// fprintf(stderr,"%s removed \n",h.hugo); fflush(stderr);  
                       new_fix_gene_count--;
                   }
               }
           }
           if (found == 0) 
           {
               fprintf(stderr,"Note: can not find geneid %d in deal with universe\n",geneid); 
               exit(0);
           }
       }
       if (binpath[i].numgenes != new_fix_gene_count)
            binpath[i].numgenes = new_fix_gene_count;
    }
//  fprintf(stderr,"in dealwithunivese R 3\n"); fflush(stderr); 
    j = 0;
    for (i=0 ; i<numbingenes ; i++) 
    {
       if (hugos[i].status == 1) j++;
    }
// fprintf(stderr,"Setting new universe to new %d ( from old universe  %d )n=%d\n",j,*numg_ptr,n); fflush(stderr);
    *numg_ptr = j;
    return 0;
}
 
#else
 // raw C version ..
static int deal_with_universe_file(char universe_fn[],int numbingenes, int numbinpaths, int *numg_ptr, unsigned char *spaceptr) 
{
    char s[20000];
    char tmphugo[100];
    int found = 0;
    int geneid;
    int *iptr;
    int o2g;
    int i,j;
    int new_fix_gene_count;
    int errorcode;
    FILE *fp;
    struct hugo_type *hugoptr;
    struct hugo_type h;
    struct bingenetype Xgenerec; 
    struct bingenetype *Xgenerec_ptr; 


    fp = fopen (universe_fn,"r");
    errorcode = errno;
    if (!fp) { fprintf(stderr,"Can not open user specified  \"universe\" file - %s , errno=%d\n",universe_fn,errorcode); return -1; }

    for (i=0;i<numbingenes;i++) hugos[i].status = 0;  // assume "guilty"
    s[0] = (char)0;
    while (fgets(s,20000,fp))
    {
        for (i=0;s[i];i++) { if (s[i] == '\n') s[i] = (char)0; if (s[i] == '\r') s[i] = (char)0; }
        h.hugo = strdup(s);
        hugoptr = bsearch(&h,hugos,numbingenes,sizeof(struct hugo_type),cmp_hugo);
        if (hugoptr) hugoptr->status = 1; // "innocent"
        free(h.hugo);
    }
    fclose(fp);

    for (i=0 ; i<numbinpaths ; i++) // for each pathway
    {
       new_fix_gene_count = binpath[i].numgenes;
       o2g = binpath[i].offset2geneids;
       if (o2g<=0) continue;
       iptr = (int *)(spaceptr + o2g);
       for (j=0 ; j<binpath[i].numgenes ; j++)  // for each gene for this pathway
       {
           geneid = *(iptr+j);
           found = 0;
           // memset(&Xgenerec,0,sizeof(struct bingenetype)); not needed
           Xgenerec.geneid = geneid;
           Xgenerec_ptr = bsearch(&Xgenerec,&bingene[0],numbingenes,sizeof(struct bingenetype),cmp_bingene);
           if (Xgenerec_ptr)
           {
               strcpy(tmphugo, Xgenerec_ptr->hugo);
               h.hugo = &tmphugo[0];
               hugoptr = bsearch(&h,hugos,numbingenes,sizeof(struct hugo_type),cmp_hugo);
               if (hugoptr)
               {
                   found = 1;
                   if (hugoptr->status == 0) // not in universe
                   {
//                fprintf(stderr,"%s removed \n",tmphugo); fflush(stderr);  
                       new_fix_gene_count--;
                   }
               }
           }
           if (found == 0) 
           {
               fprintf(stderr,"Error: can not find %d\n",geneid); 
               exit(0);
           }
       }
       if (binpath[i].numgenes != new_fix_gene_count)
            binpath[i].numgenes  = new_fix_gene_count;
    }
    j = 0;
    for (i=0;i<numbingenes;i++)
    {
       if (hugos[i].status == 1) j++;
    }
fprintf(stderr,"fixing new universe to new %d from old %d\n",j,*numg_ptr); fflush(stderr);
    *numg_ptr = j;
    return 0;
}
#endif

void usage(void)
{
fprintf(stderr,"l2p : \"list to pathways\" program.\n");
fprintf(stderr,"Usage: cat listofHUGOgenes_one_per_line.txt | l2p [optional args]\n");
fprintf(stderr,"possible optional args are ...\n");
fprintf(stderr," -help\n");
fprintf(stderr," -precise\n");
fprintf(stderr," -noheader\n");
fprintf(stderr," -universe=Universefile_one_gene_per_line.txt\n");
fprintf(stderr," -categories=listofcatgories  (comma separated)\n");
fprintf(stderr,"    example: -categories=KEGG,REACTOME,BIOCYC,PANTH  (i.e. only use genes in those 4 categories\n");
fprintf(stderr,"    another -categories example: \"-categories=H,C6\"  (i.e. only use msig's Hallmark and C6 (cancer) category pathways\n");
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
fprintf(stderr,"    C5 - MSigDB only, gene sets  consist of genes annotated by the same GO terms.\n");
fprintf(stderr,"    C6 - MSigDB only, oncogenic gene sets defined directly from microarray data from cancer gene perturbations.\n");
fprintf(stderr,"    C7 - MSigDB only, immunologic gene sets  from microarray data from immunologic studies.\n");
fprintf(stderr,"    H - MSigDB only, hallmark gene sets: signatures from MSigDB gene sets to represent biological processes.\n");
fprintf(stderr,"Example run : printf \"TP53\\nPTEN\\nAPC\\nKRAS\\nNRAS\\n\" | ./l2p -categories=PID | sort -k2,2n -k1,1n -k3,3nr | head\n");
            fflush(NULL);
}

void int2bin(int n, char s[]) 
{
    int i;
    for (i=0;i<40;i++) s[0] = (char)0;
    // determine the number of bits needed ("sizeof" returns bytes)
    int nbits = sizeof(n) * 8;
    // forcing evaluation as an unsigned value prevents complications
    // with negative numbers at the left-most bit
    unsigned int u = *(unsigned int*)&n;
    unsigned int mask = 1 << (nbits-1); // fill in values right-to-left
    for (i = 0; i < nbits; i++, mask >>= 1)
    {
            s[i] = ((u & mask) != 0) + '0';
            s[i+1] =  (char)0;
    }
    return;
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

static int parsecats(char *z)
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
            catspat |= bit;
        else
        {
            fprintf(stderr,"ERROR: invalid category = %s\n",ts[k]);
            usage();
            exit(0);
        }
    }
    j = bitCount(catspat) ;
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

int l2p_init_C(int Rflag,int argc,char *argv[])
{
    char *z;
    int i;

    universe_file[0] = 0;
    for (i=1 ; i<argc ; i++)
    {
        if ((strcmp(argv[i],"-help") == 0) || (strcmp(argv[i],"--help") == 0) || ((strcmp(argv[i],"help") == 0) ) )
        {
             usage();
             exit(1);
        }
        if (strncmp(argv[i],"-categories=",12) == 0)
        {
            z = argv[i], 
            z += 12;
            parsecats(z);
        }
        if (strcmp(argv[i],"-precise") == 0)
              precise_flag = 1; // print more digits out some user doesn't complain about "real pvals"
        if (strcmp(argv[i],"-noheader") == 0)
              no_header_flag = 1;
        if (strncmp(argv[i],"-universe=",10) == 0)
        {
            strcpy(universe_file,argv[i] + 10);
            user_universe_flag = 1;
fprintf(stderr,"note: using \"%s\" as universe file\n",universe_file);
        }
/*
        else if (strcmp(argv[i],"enc") == 0)
        {
            fprintf(stderr,"info: encode mode\n");
            // encode_flag = 1;
        }
*/
    } 
    if (catspat == 0) category_set_all(&catspat);
    return 0;
}

int l2p_init_R()
{
    char s[100];

//   fprintf(stderr,"in l2p_init_R\n"); fflush(NULL);
    int2bin(catspat,s);
//fprintf(stderr,"rpf in l2p_init_R , in catspat=%d=0x%x %s\n",catspat,catspat,s);  fflush(NULL);
    universe_file[0] = 0;
    if (catspat == 0) category_set_all(&catspat);

    return 0;
}

char* mystrcat( char* dest, char* src )
{
     int kick = 0;
     while (*dest) dest++;
     while (1)
     {
         if (*src == (char)0) kick = 1;
         *dest++ = *src++;
         if (kick) break;
     }
     return --dest;
}
 
int load_bin_dat(void)
{
    FILE *fp;
    int system_errorcode;
    int i,j;
    long int li;
    size_t sz;

    sprintf(fn_pathbin,"%s%s",exedir, "pathworks.bin");
    sprintf(fn_genesbin,"%s%s",exedir,"pathworksgenes.bin");
    sprintf(fn_spacebin,"%s%s",exedir,"pathworksspace.bin");

    fp = fopen (fn_pathbin,"r");
    system_errorcode = errno;
    if (!fp) 
    { 
        fprintf(stderr,"Can not open \"path binary data\" - %s - errno=%d\n",fn_pathbin,system_errorcode);  fflush(NULL);
        goto LOAD_ERROR;
    }
    fseek(fp, 0L, SEEK_END);
    li = ftell(fp);
    numbinpaths = li / sizeof(struct binpathouttype);
    rewind(fp);
    for (j=i=0;i<numbinpaths;i++)
    {
        sz = fread(&binpath[i],sizeof(struct binpathouttype),1,fp);
        if (sz != 1) {fprintf(stderr,"ERROR reading %dth record of %s sz=%zu recsize=%zu\n",j,fn_pathbin,sz,sizeof(struct binpathouttype)); fflush(NULL); exit(0); }
        if (catspat & binpath[i].category)
           j++;
        else 
        {
/*
char ss1[100];
char ss2[100];
int2bin(catspat,ss1);
int2bin(binpath[i].category,ss2);
catspat,ss1,i,binpath[i].category,binpath[i].category,ss2); 
*/
        }
    }
    numusedpaths= j;
    fclose(fp);
    fp = (FILE *)0;
// fprintf(stderr,"rpf in l2p_core numusedpaths=%d , numbinpaths=%d \n",numusedpaths,numbinpaths); fflush(stderr);

    fp = fopen (fn_genesbin,"r");
    system_errorcode = errno;
    if (!fp) { fprintf(stderr,"Can not open \"genes data\". filename: %s, errno=%d\n",fn_genesbin,system_errorcode); exit(0); }
    fseek(fp, 0L, SEEK_END);
    li = ftell(fp);
    numbingenes = (int)(li / sizeof(struct bingenetype));
// fprintf(stderr,"rpf in l2p_core numbingenes=%d \n",numbingenes); fflush(stderr);
    rewind(fp);
    numg = 0;
    for (i=0 ; ((i<numbingenes)&&(i<MAXGENE)) ; i++)
    {
        sz = fread(&bingene[i],sizeof(struct bingenetype),1,fp);
        system_errorcode = errno;
        if (sz != 1) {fprintf(stderr,"ERROR reading genes at %dth record of \"%s\". errno=%d\n",i,fn_genesbin,system_errorcode); fflush(NULL); exit(0); }
        if (catspat == 0) numg++;
        else if (catspat & bingene[i].categories) numg++;
    }
    fclose(fp);
    fp = (FILE *)0;


// -- for "spill over" , path and genes have some fields which are pointers (actually offsets!) to variable sized data in space
// these are accessed by integer offset (i.e. pointers) to null terminated strings or known sized array of ints
    fp = fopen (fn_spacebin,"r");
    system_errorcode = errno;
    if (!fp) { fprintf(stderr,"Can not open \"path space data\" - %s - errno=%d\n",fn_spacebin,system_errorcode); exit(0); }
    fseek(fp, 0L, SEEK_END);
    spacesize = ftell(fp);
    rewind(fp);
    ucz = malloc((size_t)spacesize); if (!ucz) { fprintf(stderr,"ERROR: can't malloc in space for \"spill\" data\n"); return 0; }
    sz = fread(ucz,spacesize,1,fp);
    if (sz != 1) {fprintf(stderr,"ERROR reading \"spill space\" file %s\n",fn_spacebin); fflush(NULL); exit(0); }
    fclose(fp);
    fp = (FILE *)0;

    return 0;
LOAD_ERROR:
    if (ucz) { free(ucz); ucz = (unsigned char *)0; }
    return -1;
}



#ifdef L2P_USING_R
SEXP l2p_core( int Rflag, int numingenes,char *genelist[], SEXP ulist, int msigflagarg)
#else
static int l2p_core( int Rflag, int numingenes,char *genelist[], int msigflagarg)
#endif
{
    char tmps_cat[40]; // temp string for "category" 
    struct hugo_type *hugoptr;
    struct hugo_type h;
    struct bingenetype *generec_ptr; 
    struct binpathouttype *binpathptr;
    unsigned int *usintptr = (void *)0;
    int status;
    FILE *fp;
    int j;
    int hitcnt;
    int i;
    int idx;
    int system_errorcode;
    size_t sz;
    long int li;
    int ingenecnt;

#ifdef L2P_USING_R
   char *p;
   char *p2;
   int k;
   SEXP pval; // 1 
   SEXP fdr; // 2
   SEXP ratio;                      // 3 if postive, genes are OVER REPRESENTED, if negative genes are UNDER REPRESENTED
   SEXP pathwayhitcount;            // 4 number of genes hit in pathway
   SEXP numberofgenesinpathway;     // 5 number of genes in the pathway
   SEXP inputnumberofgenes;         // 6 total count of user genes (user input)
   SEXP genesinpathwaysuniverse;    // 7 total number of unique genes in all pathways
   SEXP pathwayaccessionidentifier; // 8 canonical accession ( if availible, otherwise assigned by us )
   SEXP category;                     // 9 KEGG,REACTOME,GO,PANT(=PANTHER),PID=(pathway interaciton database)
   SEXP pathwayname;                // 10 Name of pathway
   SEXP genesinpathway;             // 11 genes_space_separated   HUGO genes from user that hit the pathway
   int maxflds = 11;
   SEXP Rret;
   SEXP cls; // R class 
   SEXP nam; // R name 
   SEXP rownam; // row names
   int protect_cnt = 0;

//    fprintf(stderr,"rpf inl2p_core numingens = %d \n",numingenes); 
#else
    char s[20000];
// fprintf(stderr,"rpf core pats=0x%x\n",catspat); fflush(NULL);
#endif

    // numg = 28992
#if 0
    if (encode_flag)
    {
// not used
        // numg = 26568;
        sprintf(fn_pathbin,"%s%s",exedir,"encpath.bin");
        sprintf(fn_genesbin,"%s%s",exedir,"encgenes.bin");
        sprintf(fn_spacebin,"%s%s",exedir,"encspace.bin");
    }
    else if (panther_only_flag)
    {
// not used
        // numg = 2238;
        sprintf(fn_pathbin, "%s%s",exedir,"pantherworks.bin");
        sprintf(fn_genesbin,"%s%s",exedir,"pantherworksgenes.bin");
        sprintf(fn_spacebin,"%s%s",exedir,"pantherworksspace.bin");
    }
    else if (msig_flag)
    {
// numg = 25785;
        sprintf(fn_pathbin, "%s%s",exedir,"msigpath.bin");
        sprintf(fn_genesbin,"%s%s",exedir,"msiggenes.bin");
        sprintf(fn_spacebin,"%s%s",exedir,"msigspace.bin");
    }
    else
#endif

    if (load_bin_dat() < 0)
    {
        return (SEXP)0;
    }

             // -- need hugo name access, so set up a parallel array
    status = setup_hugo_parallel_to_genebin( numbingenes);
    if (status != 0)
        return 0;

    if (user_universe_flag)
    {
#ifdef L2P_USING_R
        if (R_deal_with_universe_list(ulist,numbingenes,numbinpaths,&numg,ucz) != 0)
        {
// fprintf(stderr,"rpf after R_deal_with_universe_list(), return badness\n");
            return 0;
        }
#else
        if (deal_with_universe_file(universe_file,numbingenes,numbinpaths,&numg,ucz) != 0)
            return 0;
#endif
    }

    ingenecnt = 0;

#ifdef L2P_USING_R
    PROTECT(Rret = Rf_allocVector(VECSXP, 11)); // a list with 11 elements
    protect_cnt++;
    for (i=0 ; i<maxflds ; i++) // maxflds = 11 for now
    {
        PROTECT(pval=Rf_allocVector(REALSXP, numusedpaths ));
        PROTECT(fdr=Rf_allocVector(REALSXP, numusedpaths));
        PROTECT(ratio=Rf_allocVector(REALSXP, numusedpaths));
        PROTECT(pathwayhitcount=Rf_allocVector(INTSXP, numusedpaths));
        PROTECT(numberofgenesinpathway=Rf_allocVector(INTSXP, numusedpaths));
        PROTECT(inputnumberofgenes=Rf_allocVector(INTSXP, numusedpaths));
        PROTECT(genesinpathwaysuniverse=Rf_allocVector(INTSXP, numusedpaths));
        PROTECT(pathwayaccessionidentifier=Rf_allocVector(STRSXP, numusedpaths));
        PROTECT(category=Rf_allocVector(STRSXP, numusedpaths));   // is natively an int, but convert to string
        PROTECT(pathwayname=Rf_allocVector(STRSXP, numusedpaths));
        PROTECT(genesinpathway=Rf_allocVector(STRSXP, numusedpaths));
        protect_cnt+=11;
    }

// fprintf(stderr,"rpf after numingenes=%d\n",numingenes);
    for (k=0;k<numingenes;k++)
#else
    while ( fgets(s, sizeof(s), stdin) ) // gets() function is deprecated
#endif
    {
#ifdef L2P_USING_R
        h.hugo = strdup(genelist[k]);
#else
        for (i=0;s[i];i++) { if (s[i] == '\n') s[i] = (char)0; if (s[i] == '\r') s[i] = (char)0; }
        h.hugo = strdup(s);
#endif
        hugoptr = bsearch(&h,  hugos,numbingenes,sizeof(struct hugo_type),cmp_hugo);
        if (hugoptr)
        {
            generec_ptr = hugoptr->generec_ptr;
            if (!generec_ptr) { fprintf(stderr,"ERROR: null generecptr\n"); fflush(stderr); exit(0); }
            if ((generec_ptr)->pathcount)
            {
                 ingenecnt++;
                 if (generec_ptr->pathplace == -1) { fprintf(stderr,"ERROR: should not get pathplace of -1 with a pathcount\n");  exit(0); }
                 if (generec_ptr->pathplace == 0) { fprintf(stderr,"ERROR: should not get pathplace of 0 with a pathcount\n");  exit(0); }
                 for (i=0; i < generec_ptr->pathcount ; i++) // go through the paths for this gene, add hit to paths
                 {
                     idx = *(int *)(ucz + generec_ptr->pathplace + (i*4)); // get INDEX of record in binpath[]
                     binpathptr = &binpath[idx];
                     if ((binpathptr->hits) == (void *)0)
                     {
                         binpathptr->hits = (void *)malloc((size_t)(binpathptr->numgenes+2) *4);
                         usintptr = (unsigned int *)(binpathptr->hits);
                         *(usintptr+1) = (unsigned int)(generec_ptr - &bingene[0] ); // idx
                         *(usintptr) = (unsigned int) 1;
                     }
                     else
                     {
                         usintptr = (unsigned int *)(binpathptr->hits); // point to already allocated array
                         unsigned int curcount = *usintptr;
                         curcount = curcount + 1;
                         *usintptr = (unsigned int )(curcount);
                         *(usintptr+curcount) = (unsigned int)(generec_ptr - &bingene[0] ); // idx
                     }
                }
            }
        }
        if (h.hugo) { free(h.hugo); h.hugo = (char *)0; }
    }
// fprintf(stderr,"rpf before do_pvals_and_bh numg=%d\n",numg);
    do_pvals_and_bh(numg,ingenecnt);

#ifdef L2P_USING_R
#else
    if (no_header_flag == 0)
        printf("pval\tfdr\tratio\tpathwayhitcount\tnumgenesinpw\tpathway\tinputnumofgenes\tgenesinpathwaysuniverse\tpathwayaccessionid\tcategory\tpathwayname\tpathwaytype\tgenes\n");
#endif
    for (i=0 ; i<numusedpaths; i++)
    {
            binpathptr = &binpath[usedpaths[i].index];
            tmps_cat[0] = (char)0;
            category_code_to_string( (binpathptr)->category , tmps_cat);
            if (binpathptr->hits) hitcnt = *(unsigned int *)(binpathptr->hits);
            else hitcnt = 0;

#ifdef L2P_USING_R
            REAL(pval)[i] = usedpaths[i].pval;
            REAL(fdr)[i] = usedpaths[i].fdr;
            REAL(ratio)[i] = usedpaths[i].ad;
            INTEGER(pathwayhitcount)[i] = hitcnt;
            INTEGER(numberofgenesinpathway)[i] = (binpathptr)->numgenes;
            INTEGER(inputnumberofgenes)[i] = ingenecnt;
            INTEGER(genesinpathwaysuniverse)[i] = numg;
            SET_STRING_ELT(pathwayaccessionidentifier, i, mkChar((char *)(ucz+((binpathptr)->accession ) ) ) ); 
            SET_STRING_ELT(category, i, mkChar(tmps_cat));
            SET_STRING_ELT(pathwayname, i, mkChar( (char *)(ucz+((binpathptr)->name ) ) ) );
            SET_STRING_ELT(genesinpathway, i, mkChar( type2string((binpathptr)->type)) );

            p2 =  malloc((hitcnt+2) * 34);  // plus some extra space
            memset(p2,0,(hitcnt+2) * 34);
            if (binpathptr->hits)
            {
                p = p2;
                usintptr = (unsigned int *)(binpathptr->hits);
                for (j=1;j<(hitcnt + 1);j++)
                {
                     generec_ptr =  &bingene[*(usintptr+j)];
                     if (j>1) p = mystrcat(p," ");
                     p = mystrcat(p,generec_ptr->hugo);
                }
            }
            SET_STRING_ELT(genesinpathway, i, mkChar( p2 ) );
            if (p2) { free(p2); p2 = (char *)0; }
#else
            if (precise_flag == 1)
            {
                printf("%20.18f\t%20.18f\t%11.9f\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t",
                    usedpaths[i].pval, usedpaths[i].fdr, usedpaths[i].ad,
                             hitcnt, (binpathptr)->numgenes, ingenecnt, numg,
                             (char *)(ucz+(binpathptr)->accession), tmps_cat,
                             (ucz+((binpathptr)->name)), type2string((binpathptr)->type));
            }
            else
            {
                printf("%11.9f\t%11.9f\t%11.9f\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t",
                    usedpaths[i].pval, usedpaths[i].fdr, usedpaths[i].ad,
                             hitcnt, (binpathptr)->numgenes, ingenecnt, numg,
                             (char *)(ucz+(binpathptr)->accession), tmps_cat,
                             (ucz+((binpathptr)->name)), type2string((binpathptr)->type));
            }
            if (binpathptr->hits)
            {
                usintptr = (unsigned int *)(binpathptr->hits);
                for (j=1;j<(hitcnt + 1);j++)
                {
                     generec_ptr =  &bingene[*(usintptr+j)];
                     printf("%s ",generec_ptr->hugo);
                }
            }
            printf("\n");
#endif
    }

#ifdef L2P_USING_R
   SET_VECTOR_ELT( Rret,0, pval);
   SET_VECTOR_ELT( Rret,1, fdr);
   SET_VECTOR_ELT( Rret,2, ratio);
   SET_VECTOR_ELT( Rret,3, pathwayhitcount);
   SET_VECTOR_ELT( Rret,4, numberofgenesinpathway);
   SET_VECTOR_ELT( Rret,5, inputnumberofgenes);
   SET_VECTOR_ELT( Rret,6, genesinpathwaysuniverse);
   SET_VECTOR_ELT( Rret,7, pathwayaccessionidentifier);
   SET_VECTOR_ELT( Rret,8, category);
   SET_VECTOR_ELT( Rret,9, pathwayname);
   SET_VECTOR_ELT( Rret,10, genesinpathway);

   PROTECT(cls = allocVector(STRSXP, 1)); // class attribute
   protect_cnt++;

   SET_STRING_ELT(cls, 0, mkChar("data.frame"));
   classgets(Rret, cls);

   PROTECT(nam = allocVector(STRSXP, 11)); // names attribute (column names)
   protect_cnt++;

   SET_STRING_ELT( nam ,                      0, mkChar("pval"));
   SET_STRING_ELT( nam,                        1, mkChar("fdr"));
   SET_STRING_ELT( nam,                      2, mkChar("ratio"));
   SET_STRING_ELT( nam,            3, mkChar("pathwayhitcount"));
   SET_STRING_ELT( nam,     4, mkChar("numberofgenesinpathway"));
   SET_STRING_ELT( nam,         5, mkChar("inputnumberofgenes"));
   SET_STRING_ELT( nam,    6, mkChar("genesinpathwaysuniverse"));
   SET_STRING_ELT( nam, 7, mkChar("pathwayaccessionidentifier"));
   SET_STRING_ELT( nam,                   8, mkChar("category"));
   SET_STRING_ELT( nam,                9, mkChar("pathwayname"));
   SET_STRING_ELT( nam,            10, mkChar("genesinpathway"));
   namesgets(Rret, nam);

   PROTECT(rownam = allocVector(STRSXP, numusedpaths )); // row.names attribute
   protect_cnt++;
   for (i=0;i<numusedpaths;i++)
   {
       binpathptr = &binpath[usedpaths[i].index];
       SET_STRING_ELT(rownam, i, mkChar( (char *)(ucz+((binpathptr)->accession ) ) ) );
   }
   setAttrib(Rret, R_RowNamesSymbol, rownam);

// CORE_END:
   UNPROTECT(protect_cnt);
   for (i=0;i<numingenes;i++)
     if (genelist[i]) {free(genelist[i]); }
   free(genelist);

   free_hugos();
   free_binpaths();

   if (ucz) { free(ucz); ucz = (void *)0; }
// printf("in l2p_core() return Rret = %p\n",Rret); fflush(stdout);
   return Rret;
#else
   for (i=0;i<numingenes;i++)
     if (genelist[i]) { free(genelist[i]); }
   free(genelist);

   free_hugos();
   free_binpaths();
   if (ucz) { free(ucz); ucz = (void *)0; }
   return 0;
#endif

CORE_ERROR:
   for (i=0;i<numingenes;i++)
     if (genelist[i]) {free(genelist[i]); }
   free(genelist);

   free_hugos();
   free_binpaths();
   if (ucz) { free(ucz); ucz = (void *)0; }
   return 0;
}


#ifdef L2P_USING_R

SEXP l2p(SEXP lst, SEXP fpath)
{
   char path[PATH_MAX];
   int lastslash;
   int i;
   int len;
   char **z;

// fprintf(stderr,"in  l2p() 0 \n"); fflush(stderr);
   user_universe_flag = 0;
   catspat=0;
   len = length(lst);
   strncpy(path,CHAR(STRING_ELT(fpath, 0)),PATH_MAX-2);
   z = (char **)malloc(sizeof(char *)*len);
   if (!z)
   {
        return (SEXP) -1;
   }
   for (i = 0; i < len; i++)
   {
       *(z+i) = strdup(CHAR(STRING_ELT(lst, i)));
   }
  lastslash = -1;
  for (i=0;path[i];i++)
  {
      if (path[i] == '/') lastslash = i;
  }
  if (lastslash > 0) path[lastslash] = (char)0;
  strcpy(exedir,path);
  strcat(exedir,"/");

  l2p_init_R();
  return l2p_core(1,len, z, (void *)0,0);
}

SEXP l2pgetlongdesc(SEXP accarg, SEXP fpath)
{
    char path[PATH_MAX];
    char acc[PATH_MAX];
    int i;	  
    struct binpathouttype *binpathptr;
    int lastslash;
    char *z;
 
    strncpy(acc,CHAR(STRING_ELT(accarg, 0)),PATH_MAX-2);
    strncpy(path,CHAR(STRING_ELT(fpath, 0)),PATH_MAX-2);
    lastslash = -1;
    for (i=0;path[i];i++)
    {
        if (path[i] == '/') lastslash = i;
    }
    if (lastslash > 0) path[lastslash] = (char)0;
    strcpy(exedir,path);
    strcat(exedir,"/");
    if (load_bin_dat() < 0)
    {
        return (SEXP)0;
    }
    for (i=0;i<numbinpaths   ;i++)
    {
        binpathptr = &binpath[i];
        z = (char *)(ucz+((binpathptr)->desc ) );
// fprintf(stderr,"l2pgetlongdesc %s %s\n",acc,z); fflush(NULL);
	if (z)
	{
            if (strcmp(acc,z) == 0)
               break;
	}
	else z = (char *)0;
    }
    SEXP ret = PROTECT(allocVector(STRSXP, 1));
    if (z)
        SET_STRING_ELT(ret, 0, mkChar(z));
    else
        SET_STRING_ELT(ret, 0, mkChar(""));
    UNPROTECT(1);
    return ret;
}

SEXP l2pgetgenes4acc(SEXP accarg, SEXP fpath)
{
    char path[PATH_MAX];
    char acc[PATH_MAX];
    char tmphugo[1024];
    struct hugo_type *hugoptr;
    struct hugo_type h;
    int *iptr;
    int i,j;	  
    int o2g;
    struct binpathouttype *binpathptr;
    char *z;
    char *z2;
    int geneid;
    struct bingenetype Xgenerec; 
    struct bingenetype *Xgenerec_ptr; 
    int lastslash;
 
    strncpy(acc,CHAR(STRING_ELT(accarg, 0)),PATH_MAX-2);
    strncpy(path,CHAR(STRING_ELT(fpath, 0)),PATH_MAX-2);
    lastslash = -1;
    for (i=0;path[i];i++)
    {
        if (path[i] == '/') lastslash = i;
    }
    if (lastslash > 0) path[lastslash] = (char)0;
    strcpy(exedir,path);
    strcat(exedir,"/");
    if (load_bin_dat() < 0)
    {
        return (SEXP)0;
    }
    z2 = (char *)0;
    for (i=0;i<numbinpaths   ;i++)
    {
        binpathptr = &binpath[i];
        z = (char *)(ucz+((binpathptr)->accession ) );
        if (strcmp(acc,z) == 0)
	{
	    o2g = binpathptr->offset2geneids;
            if (o2g<=0) { z2 = (char *)0; break; }
            SEXP ret = PROTECT(allocVector(STRSXP, binpathptr->numgenes));
            iptr = (int *)(ucz + o2g);
            for (j=0 ; j<binpath[i].numgenes ; j++)  // for each gene for this pathway
            {
		tmphugo[0] = (char)0;
                geneid = *(iptr+j);
                Xgenerec.geneid = geneid;
                Xgenerec_ptr = bsearch(&Xgenerec,&bingene[0],numbingenes,sizeof(struct bingenetype),cmp_bingene);
                if (Xgenerec_ptr)
                    strcpy(tmphugo, Xgenerec_ptr->hugo);
                SET_STRING_ELT(ret, j, mkChar(tmphugo));
            }
            UNPROTECT(1);
	    return ret;
        }
    }
    return (SEXP)0;
}

SEXP l2pwcats (SEXP lst, SEXP catsarg, SEXP fpath)
{
   char path[PATH_MAX];
   char cats[PATH_MAX];
   int lastslash;
   int i;
   int len;
   char **z;

   user_universe_flag = 0;
   catspat=0;
   user_universe_flag = 0;
   len = length(lst);
   strncpy(path,CHAR(STRING_ELT(fpath, 0)),PATH_MAX-2);
   strncpy(cats,CHAR(STRING_ELT(catsarg, 0)),PATH_MAX-2);
   parsecats(cats); // set catpats
   z = (char **)malloc(sizeof(char *)*len);
   if (!z)
        return (SEXP) -1;
   for (i = 0; i < len; i++)
       *(z+i) = strdup(CHAR(STRING_ELT(lst, i)));
   lastslash = -1;
   for (i=0;path[i];i++)
   {
       if (path[i] == '/') lastslash = i;
   }
   if (lastslash > 0) path[lastslash] = (char)0;

   strcpy(exedir,path);
   strcat(exedir,"/");

   l2p_init_R();
   return l2p_core(1,len, z, (void *)0,0);
}


SEXP l2pmsig(SEXP lst, SEXP fpath)
{
   char path[PATH_MAX];
   int lastslash;
   int i;
   int len;
   char **z;


   user_universe_flag = 0;
   len = length(lst);
   strncpy(path,CHAR(STRING_ELT(fpath, 0)),PATH_MAX-2);
   z = (char **)malloc(sizeof(char *)*len);
   if (!z)
   {
        return (SEXP) -1;
   }
   for (i = 0; i < len; i++)
   {
       *(z+i) = strdup(CHAR(STRING_ELT(lst, i)));
   }
  lastslash = -1;
  for (i=0;path[i];i++)
  {
      if (path[i] == '/') lastslash = i;
  }
  if (lastslash > 0) path[lastslash] = (char)0;

  strcpy(exedir,path);
  strcat(exedir,"/");

  l2p_init_R();
  return l2p_core(1,len, z, (void *)0,1);
}

SEXP l2pmsigwcats(SEXP lst, SEXP catsarg, SEXP fpath)
{
   char cats[PATH_MAX];
   char path[PATH_MAX];
   int lastslash;
   int i;
   int len;
   char **z;


   user_universe_flag = 0;
   catspat=0;
   len = length(lst);
   strncpy(path,CHAR(STRING_ELT(fpath, 0)),PATH_MAX-2);
   strncpy(cats,CHAR(STRING_ELT(catsarg, 0)),PATH_MAX-2);
   parsecats(cats); // set catpats
   z = (char **)malloc(sizeof(char *)*len);
   if (!z)
   {
        return (SEXP) -1;
   }
   for (i = 0; i < len; i++)
   {
       *(z+i) = strdup(CHAR(STRING_ELT(lst, i)));
   }
  lastslash = -1;
  for (i=0;path[i];i++)
  {
      if (path[i] == '/') lastslash = i;
  }
  if (lastslash > 0) path[lastslash] = (char)0;

// hardwire for tesing strcpy(exedir,"/Users/finneyr/R/libs/l2p/extdata/");
  strcpy(exedir,path);
  strcat(exedir,"/");

  l2p_init_R();
  return l2p_core(1,len, z, (void *)0,1);
}

SEXP l2pu(SEXP lst, SEXP ulst, SEXP fpath )
{
   char path[PATH_MAX];
   int lastslash;
   int i;
   int len;
   char **z;

   user_universe_flag = 1;
   catspat=0;
   if (!lst)
       return (SEXP)-1;
   if (!ulst)
       return (SEXP)-2;
   len = length(lst);
   strncpy(path, CHAR(STRING_ELT(fpath, 0)),PATH_MAX-2);
   z = (char **)malloc(sizeof(char *)*len);
   if (!z)
       return (SEXP)-3;
   for (i = 0; i < len; i++)
   {
       *(z+i) = strdup(CHAR(STRING_ELT(lst, i)));
   }
   lastslash = -1;
   for (i=0;path[i];i++)
   {
      if (path[i] == '/') lastslash = i;
   }
   if (lastslash > 0) path[lastslash] = (char)0;

   strcpy(exedir,path);
   strcat(exedir,"/");

    char s[100];
    memset(s,0,sizeof(s)); 
    int2bin(catspat,s);
// fprintf(stderr,"rpf BEFORE  l2p_init_R , in catspat=%d=0x%x %s\n",catspat,catspat,s); 
    l2p_init_R();
    memset(s,0,sizeof(s)); 
    int2bin(catspat,s);
// fprintf(stderr,"rpf after  l2p_init_R , in catspat=%d=0x%x %s\n",catspat,catspat,s); 
   return l2p_core(1,len, z,ulst,0);
}

SEXP l2puwcats(SEXP lst, SEXP ulst, SEXP catsarg, SEXP fpath )
{
   char path[PATH_MAX];
   char cats[PATH_MAX];
   int lastslash;
   int i;
   int len;
   char **z;

   user_universe_flag = 1;
   catspat=0;
//   fprintf(stderr,"in l2puwcats 1\n"); fflush(NULL);
   if (!lst)
       return (SEXP)-1;
   if (!ulst)
       return (SEXP)-2;
   len = length(lst);
   strncpy(path, CHAR(STRING_ELT(fpath, 0)),PATH_MAX-2);
   strncpy(cats,CHAR(STRING_ELT(catsarg, 0)),PATH_MAX-2);
   parsecats(cats); // set catpats
   z = (char **)malloc(sizeof(char *)*len);
   if (!z)
       return (SEXP)-3;
   for (i = 0; i < len; i++)
   {
       *(z+i) = strdup(CHAR(STRING_ELT(lst, i)));
   }
   lastslash = -1;
   for (i=0;path[i];i++)
   {
      if (path[i] == '/') lastslash = i;
   }
   if (lastslash > 0) path[lastslash] = (char)0;

   strcpy(exedir,path);
   strcat(exedir,"/");

   l2p_init_R();
   return l2p_core(1,len, z,ulst,0);
}


#if 0
SEXP l2pumsig(SEXP lst, SEXP ulst, SEXP fpath )
{
   char path[PATH_MAX];
   int lastslash;
   int i;
   int len;
   char **z;

   if (!lst)
   {
       return (SEXP)-1;
   }
   if (!ulst)
   {
       return (SEXP)-2;
   }
   user_universe_flag = 1;
   len = length(lst);
   strncpy(path, CHAR(STRING_ELT(fpath, 0)),PATH_MAX-2);
   z = (char **)malloc(sizeof(char *)*len);
   if (!z)
   {
       return (SEXP)-3;
   }
   for (i = 0; i < len; i++)
   {
       *(z+i) = strdup(CHAR(STRING_ELT(lst, i)));
   }
   lastslash = -1;
   for (i=0;path[i];i++)
   {
      if (path[i] == '/') lastslash = i;
   }
   if (lastslash > 0) path[lastslash] = (char)0;

   strcpy(exedir,path);
   strcat(exedir,"/");

   l2p_init_R();
   return l2p_core(1,len, z,ulst,1);
}
#endif

#if 0
SEXP l2pmsigwcatsu(SEXP lst, SEXP ulst, SEXP catsarg, SEXP fpath )
{
   char cats[PATH_MAX];
   char path[PATH_MAX];
   int lastslash;
   int i;
   int len;
   char **z;

   if (!lst)
       return (SEXP)-1;
   if (!ulst)
       return (SEXP)-2;
   user_universe_flag = 1;
   len = length(lst);
   strncpy(path, CHAR(STRING_ELT(fpath, 0)),PATH_MAX-2);
   strncpy(cats,CHAR(STRING_ELT(catsarg, 0)),PATH_MAX-2);
   parsecats(cats); // set catpats
   z = (char **)malloc(sizeof(char *)*len);
   if (!z)
       return (SEXP)-3;
   for (i = 0; i < len; i++)
       *(z+i) = strdup(CHAR(STRING_ELT(lst, i)));
   lastslash = -1;
   for (i=0;path[i];i++)
   {
      if (path[i] == '/') lastslash = i;
   }
   if (lastslash > 0) path[lastslash] = (char)0;

   strcpy(exedir,path);
   strcat(exedir,"/");

   l2p_init_R();
   return l2p_core(1,len, z,ulst,1);
}
#endif

#else
#if __MACH__
// using mac os 
#include <libproc.h>
#endif
static void get_this_executable_path(char *puthere,int size)
{
    int i;
    int lastslash = 0;

    *puthere = (char)0;
 // if apple ... else linux ...
#if __MACH__
// /proc_pidpath/
    pid_t pid;
    pid = getpid();
    int ret;
    ret = proc_pidpath (pid, puthere, PATH_MAX);
    if ( ret <= 0 ) {
        fprintf(stderr, "PID %d: proc_pidpath ();\n", pid);
        fprintf(stderr, "    %s\n", strerror(errno));
    } else { // success
        // printf("proc %d: %s\n", pid, puthere);
    }
// does not work    strcpy(rlinfo,"/proc/curproc/file");
//  alternate method:  (_NSGetExecutablePath(path, &size) == 0)
#else
    char rlinfo[PATH_MAX];
    strcpy(rlinfo,"/proc/self/exe");
    if (readlink(rlinfo, puthere,PATH_MAX) == -1)
    {
         fprintf(stderr,"ERROR can not access full path of (this!) executable\n");
                     // what to do ? not sure.
         return;
    }
#endif
    for (i=0;puthere[i];i++)
        if (puthere[i] == '/') lastslash = i;
    puthere[lastslash+1] = (char)0;
    return;
// readlink("/proc/curproc/file", buf, bufsize) (FreeBSD)
// readlink("/proc/self/path/a.out", buf, bufsize)
// On Windows: use GetModuleFileName(NULL, buf, bufsize)
}

static int l2p_run_for_C(int argc,char *argv[])
{
    l2p_init_C(0,argc,argv);
    get_this_executable_path(exedir,PATH_MAX);
// fprintf(stderr,"rpf in l2p_run_for_C before get_this_executable_path pats=0x%x\n",catspat); fflush(NULL);
    l2p_core(0, 0,(void *)0,0);
    return 0;
}

int main(int argc,char *argv[])
{
    return l2p_run_for_C(argc,argv);
}

#endif



void category_set_all(int *pat)
{
    *pat = ( 
   CAT_NCBI_BIOCYC| CAT_NCBI_GO    | CAT_NCBI_KEGG  | CAT_NCBI_PANTH | CAT_NCBI_PID | CAT_NCBI_REACTOME     | CAT_NCBI_WikiPathways | 
   CAT_MSIG_C1   | CAT_MSIG_C2   | CAT_MSIG_C3   | CAT_MSIG_C4   | CAT_MSIG_C5   | CAT_MSIG_C6   | CAT_MSIG_C7   | CAT_MSIG_H  ) ;
}
void category_code_to_string(int cat,char puthere[])
{
         if (cat & CAT_NCBI_BIOCYC) strcpy(puthere,"BIOCYC"); 
    else if (cat & CAT_NCBI_GO    ) strcpy(puthere,"GO"); 
    else if (cat & CAT_NCBI_KEGG  ) strcpy(puthere,"KEGG"); 
    else if (cat & CAT_NCBI_PANTH ) strcpy(puthere,"PANTH"); 
    else if (cat & CAT_NCBI_PID ) strcpy(puthere,"PID");
    else if (cat & CAT_NCBI_REACTOME     ) strcpy(puthere,"REACTOME");
    else if (cat & CAT_NCBI_WikiPathways ) strcpy(puthere,"WikiPathways");
    else if (cat & CAT_MSIG_C1   ) strcpy(puthere,"C1");
    else if (cat & CAT_MSIG_C2   ) strcpy(puthere,"C2");
    else if (cat & CAT_MSIG_C3   ) strcpy(puthere,"C3");
    else if (cat & CAT_MSIG_C4   ) strcpy(puthere,"C4");
    else if (cat & CAT_MSIG_C5   ) strcpy(puthere,"C5");
    else if (cat & CAT_MSIG_C6   ) strcpy(puthere,"C6");
    else if (cat & CAT_MSIG_C7   ) strcpy(puthere,"C7");
    else if (cat & CAT_MSIG_H    ) strcpy(puthere,"H");
    else strcpy(puthere,"UNKNOWN");
//    else if (cat & CAT_MSIG_ARCHIVED     )   strcpy(puthere,"ARCHIVED"); 
}


int string_to_category_code(char cats[])
{
if  (strcmp(cats,"BIOCYC") == 0)    return CAT_NCBI_BIOCYC; 
else if  (strcmp(cats,"GO") == 0)   return CAT_NCBI_GO; 
else if  (strcmp(cats,"KEGG") == 0) return CAT_NCBI_KEGG; 
else if  (strcmp(cats,"PANTH") == 0) return CAT_NCBI_PANTH; 

else if  (strcmp(cats,"PID") == 0) return CAT_NCBI_PID; 
else if  (strcmp(cats,"Pathway Interaction Database") == 0) return CAT_NCBI_PID;  // dupe (see previous line)

else if  (strcmp(cats,"REACTOME") == 0) return CAT_NCBI_REACTOME; 
else if  (strcmp(cats,"WikiPathways") == 0) return CAT_NCBI_WikiPathways; 
else if  (strcmp(cats,"C1") == 0) return CAT_MSIG_C1; 
else if  (strcmp(cats,"C2") == 0) return CAT_MSIG_C2; 
else if  (strcmp(cats,"C3") == 0) return CAT_MSIG_C3; 
else if  (strcmp(cats,"C4") == 0) return CAT_MSIG_C4; 
else if  (strcmp(cats,"C5") == 0) return CAT_MSIG_C5; 
else if  (strcmp(cats,"C6") == 0) return CAT_MSIG_C6; 
else if  (strcmp(cats,"C7") == 0) return CAT_MSIG_C7; 
else if  (strcmp(cats,"H") == 0) return CAT_MSIG_H; 
else return 0;
// else if  (strcmp(cats,"ARCHIVED") == 0)   return CAT_MSIG_ARCHIVED; 
}

void categories_pattern_to_strings(int cat,char puthere[])
{
    puthere[0] = (char)0;
    if (cat & CAT_NCBI_BIOCYC) strcat(puthere,"BIOCYC "); 
    if (cat & CAT_NCBI_GO    ) strcat(puthere,"GO "); 
    if (cat & CAT_NCBI_KEGG  ) strcat(puthere,"KEGG "); 
    if (cat & CAT_NCBI_PANTH ) strcat(puthere,"PANTH "); 
    if (cat & CAT_NCBI_PID ) strcat(puthere,"PID ");
    if (cat & CAT_NCBI_REACTOME     ) strcat(puthere,"REACTOME ");
    if (cat & CAT_NCBI_WikiPathways ) strcat(puthere,"WikiPathways ");
//     if (cat & CAT_MSIG_ARCHIVED     )   strcat(puthere,"ARCHIVED "); 
    if (cat & CAT_MSIG_C1    ) strcat(puthere,"C1 ");
    if (cat & CAT_MSIG_C2   ) strcat(puthere,"C2 ");
    if (cat & CAT_MSIG_C3   ) strcat(puthere,"C3 ");
    if (cat & CAT_MSIG_C4   ) strcat(puthere,"C4 ");
    if (cat & CAT_MSIG_C5   ) strcat(puthere,"C5 ");
    if (cat & CAT_MSIG_C6   ) strcat(puthere,"C6 ");
    if (cat & CAT_MSIG_C7   ) strcat(puthere,"C7 ");
    if (cat & CAT_MSIG_H    ) strcat(puthere,"H ");
}

