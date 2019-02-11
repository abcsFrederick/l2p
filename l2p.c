
/*
vim l2p.c ; gcc -O3 -Wall -o l2p l2p.c -lm 
vim l2p.c ; gcc -O3 -march=native -flto -Wall -o l2p l2p.c -lm 
vim l2p.c ; gcc -D_FORTIFY_SOURCE=2 -fstack-protector --param ssp-buffer-size=4 -fPIE -pie -Wl,-z,relro,-z,now -o l2p l2p.c -lm 
#vim l2p.c ; gcc -D_FORTIFY_SOURCE=2 -fstack-protector --param ssp-buffer-size=4 -fPIE -pie -Wl,-z,relro,-z,now (ld -z relro and ld -z now) -o l2p l2p.c -lm 

output:
     1	pval
     2	fdr
     3	ratio                      if postive, genes are OVER REPRESENTED, if negative genes are UNDER REPRESENTED
     4	pathwayhitcount            number of genes hit in pathway
     5	numberofgenesin pathway    number of genes in the pathway
     6	inputnumberofgenes          total count of user genes (user input)
     7	genesinpathwaysuniverse    total number of unique genes in all pathways
     8	pathwayaccessionidentifier canonical accession ( if availible, otherwise assigned by us )
     9	source                     KEGG,REACTOME,GO,PANT(=PANTHER),PID=(pathway interaciton database)
    10	pathwayname                Name of pathway
    11	pathwaytype genes_space_separated   HUGO genes from user that hit the pathway
*/

#if 1
#define L2P_USING_R 1
#endif 

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

#include "pathworks.h"

static int header_flag = 0;
static int precise_flag = 0;
static int user_universe_flag = 0;
static struct binpathouttype binpath[MAXBSID];
static int numbinpaths = 0;

struct binpathparallel
{
    double pval;
    double fdr;
    double ad;
    double hitcnt;
};
static struct binpathparallel pathparallel[MAXBSID];        // uses "numbinpaths" variable as index

static struct bingenetype bingene[MAXGENE];
static int numbingenes;
static int spacesize;

// -- for benjamini hochberg FDR ...
struct ordertype
{
    double val;
    int order;
};

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

static void benjaminihochberg(double pvals[],int n, double returnpvals[])
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

    o = (struct ordertype *)malloc((sizeof(struct ordertype))*n);
    for (j=0;j<n;j++) 
    {
        (o+j)->val=pvals[j];
        (o+j)->order=j+1;
    }
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
    for (j=0;j<n;j++)
        returnpvals[j] = (cummin+j)->val;
    if (i) free(i);
    if (o) free(o);
//    if (ro) free(ro);
    if (po) free(po);
    if (cummin) free(cummin);
//    if (intermed) free(intermed);

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
  double prob = 0.0;
#if 0
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

  prob = exact(n11_,n1_,n_1,n);
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

static int cmp_hugo(const void *a, const void *b)
{
    int ret;
    struct hugo_type *aa; 
    struct hugo_type *bb; 
    aa = (struct hugo_type *)a;
    bb = (struct hugo_type *)b;

    ret = strcmp(aa->hugo, bb->hugo);
    return ret;
}

static void setup_hugo_parallel_to_genebin(size_t total_array_size_in_bytes,int n)
{
    int i;

    hugos = (struct hugo_type *)malloc((size_t)total_array_size_in_bytes);
    if (!hugos) { fprintf(stderr,"ERROR: can't malloc in setup_hugo_parallel_to_genebin...()\n"); exit(0); }
    for (i=0;i<n;i++)
    {
        (hugos+i)->hugo = strdup(bingene[i].hugo);
        (hugos+i)->generec_ptr = &bingene[i]; // from struct bingenetype bingene[MAXGENE];
        // (hugos+i)->generec_ptr = (struct bingenetype *)&bingene[i]; // from struct bingenetype bingene[MAXGENE];
    }

    qsort(hugos,n,sizeof(struct hugo_type),cmp_hugo); // rearrange by hugo name for bsearch , keep link to original genebin rec 
#if 0
    for (i=0;i<n;i++)
    {
         fprintf(stderr,"cmp %s %s %s\n",(hugos+i)->hugo,(hugos+i)->generec_ptr->hugo,bingene[i].hugo);
    }
exit(0);
#endif

    return;
}



#if 0
static void dump_hits( int index)
{
    struct hit_type *hitptr;
    struct hit_type *trav;
    int i = 0;

    hitptr =(struct hit_type *)binpath[index].hits;
    trav = hitptr;
    while (trav->genename)
    {
        if (i != 0) printf(" ");
        printf("%s",trav->genename);
        if (trav->n == (void *)0) break;
        trav = trav->n;
        i++;
    }
    return;
}

static int cnt_hits(int index_into_array_of_struct)
{
    int ret = 0;
    struct hit_type *hitptr;
    struct hit_type *trav;

    hitptr =(struct hit_type *)binpath[index_into_array_of_struct].hits;
    if (!hitptr) return 0 ;
    trav = hitptr;
    while (trav->genename)
    {
        ret++;
        if (trav->n == (void *)0) break;
        trav = trav->n;
    }
    return ret;
}


static int cmp_binpath_by_bsid(const void *a, const void *b)
{
    struct binpathouttype *aa;
    struct binpathouttype *bb;
    aa = (void *)a;
    bb = (void *)b;
    if      (aa->bsid <  bb->bsid) return -1;
    else if (aa->bsid >  bb->bsid) return 1;
    return 0;
}
#endif

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


static void get_this_executable_path(char *puthere,int size)
{
    int i;
    int lastslash = 0;

    *puthere = (char)0;
    (void)readlink("/proc/self/exe", puthere,size);
    for (i=0;puthere[i];i++)
        if (puthere[i] == '/') lastslash = i;
    puthere[lastslash+1] = (char)0;
    return;
// readlink("/proc/curproc/file", buf, bufsize) (FreeBSD)
// readlink("/proc/self/path/a.out", buf, bufsize) 
// On Windows: use GetModuleFileName(NULL, buf, bufsize)
}

static int encode = 0;
static int msig = 0;
static int panther_only = 0;


static void do_exact_and_benjamini_hochberg(int numg,int incnt)
{
    struct binpathouttype *binpathptr;
    double d;
    double ad;
    int i;
    double *z;
    double *fdrs;
    int hitcnt ;

    z = (double *)malloc((size_t)(sizeof (double)*numbinpaths)); if (!z) { fprintf(stderr,"ERROR: can't malloc in do_exact_and_benjamini_hochberg() 1\n"); exit(0); }
    fdrs = (double *)malloc((size_t)(sizeof (double)*numbinpaths)); if (!fdrs) { free(z); /* clean up */ fprintf(stderr,"ERROR: can't malloc in do_exact_and_benjamini_hochberg() 2\n"); exit(0); }
    for (i=0 ; i<numbinpaths ; i++)
    {
        hitcnt = (unsigned int)0;
        binpathptr = &binpath[i];
        if (binpathptr->hits)
        {
            hitcnt = *(unsigned int *)(binpathptr->hits);
            if (hitcnt > (unsigned long int)binpathptr->numgenes) 
            {
                 fprintf(stderr,"ERROR: more hits than genes hitcnt %d > %d for %s \n", hitcnt , (binpath+i)->numgenes,(binpath+i)->accession); 
                 fprintf(stderr,"ERROR: i=%d , path=%s,%s, current hits=%u\n", i,binpathptr->accession,binpathptr->source,hitcnt); 
                 fflush(stderr);
                 exit(0); 
            }
        }
        d = exact22((int)hitcnt,(binpathptr)->numgenes,incnt,numg);
        ad = ( ((double)hitcnt/ (double)(binpathptr)->numgenes) - ((double)incnt/(double)numg) );
        pathparallel[i].pval = d;
        pathparallel[i].ad = ad;
        *(z+i) = d;
    }
    benjaminihochberg(z,numbinpaths,fdrs);
    for (i=0 ; i<numbinpaths ; i++)
    {
         pathparallel[i].fdr = *(fdrs+i);
    }
    free(z);
    free(fdrs);
    return;
}


static void deal_with_universe_file(char universe_fn[],int n, int numbinpaths, int *numg_ptr, unsigned char *spaceptr) 
{
    char s[20000];
    char tmphugo[100];
    int found = 0;
    int geneid;
    int *iptr;
    int o2g;
    int i,j,k;
    int incnt;
    int new_fix_gene_count;
    int errorcode;
    FILE *fp;
    struct hugo_type *hugoptr;
    struct hugo_type h;

    fp = fopen (universe_fn,"r");
    errorcode = errno;
    if (!fp) { fprintf(stderr,"Can not open user specified  \"universe\" file - %s , errno=%d\n",universe_fn,errorcode); exit(0); }

    for (i=0;i<n;i++) hugos[i].status = 0;  // assume "guilty"
    incnt = 0;
    while (fgets(s,20000,fp))
    {
        for (i=0;s[i];i++) { if (s[i] == '\n') s[i] = (char)0; if (s[i] == '\r') s[i] = (char)0; }
        h.hugo = strdup(s);
        hugoptr = bsearch(&h,hugos,n,sizeof(struct hugo_type),cmp_hugo);
        if (hugoptr) hugoptr->status = 1; // "innocent"
        free(h.hugo);
        incnt++;
    }
    fclose(fp);

    for (i=0 ; i<numbinpaths ; i++) // for each path
    {
       new_fix_gene_count = binpath[i].numgenes;
       o2g = binpath[i].offset2geneids;
       if (o2g<=0) continue;
       iptr = (int *)(spaceptr + o2g);
       for (j=0 ; j<binpath[i].numgenes ; j++) // for each gene for this pathway
       {
           geneid = *(iptr+j);
           found = 0;
           for (k=0;k<numbingenes;k++)
           {
               if (geneid == bingene[k].geneid) 
               {
                  found = 1;
                   strcpy(tmphugo, bingene[k].hugo);
                   h.hugo = &tmphugo[0];
                   hugoptr = bsearch(&h,hugos,n,sizeof(struct hugo_type),cmp_hugo);
                   if (hugoptr)
                   {
                       if (hugoptr->status == 0) // not in universe
                       {
                           new_fix_gene_count--;
                       }
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
       {
            binpath[i].numgenes  = new_fix_gene_count;
       }
    }
    j = 0;
    for (i=0;i<n;i++) 
    {
       if (hugos[i].status == 1) j++;
    }
fprintf(stderr,"fixing new universe to new %d from old %d\n",j,*numg_ptr);
    *numg_ptr = j;
    return;
}

static char fn_pathbin[PATH_MAX];
static char fn_genesbin[PATH_MAX];
static char fn_spacebin[PATH_MAX];
static char universe_file[PATH_MAX];
static char exedir[PATH_MAX]; // executable directory 
static int numg = 0;

int l2p_init(int Rflag,int argc,char *argv[])
{
    int i;

// old: ingenes = 35058358, binpathrecsize=72 , genrecsize=80 numbsids=53772 numgenes=60094, hitgenes=19670 outpathreccnt=18418

    numg = 19657;  // universe , get withe "cat pathworks.txt | cut -f10 | tr " " "\n" | sort | uniq"

// no!  numg = numbingenes;   no? why?

    universe_file[0] = 0;

    encode = msig = 0;
    for (i=1;i<argc;i++)
    {
        if (strcmp(argv[i],"-help") == 0)
        {
fprintf(stderr,"l2p : \"list to pathways\" program.\n");
fprintf(stderr,"Example Usage: cat listofHUGOgenes_one_per_line.txt | l2p [optional args]\n");
fprintf(stderr,"possible optional args are ...\n");
fprintf(stderr," -help\n");
fprintf(stderr," -precise\n");
fprintf(stderr," -header\n");
fprintf(stderr," -universe=Universefile.txt\n");
fprintf(stderr," enc                        (use encode, no dash)\n");
fprintf(stderr," msig                       (use msig, no dash)\n");
fprintf(stderr," panther                    (use panther, no dash)\n");
            fflush(stderr);
            return 0;
        }
        if (strcmp(argv[i],"-precise") == 0)
        {
              precise_flag = 1;
        }
        if (strcmp(argv[i],"-header") == 0)
        {
              header_flag = 1;
        }
        if (strncmp(argv[i],"-universe=",10) == 0)
        {
            strcpy(universe_file,argv[i] + 10);
            user_universe_flag = 1;
fprintf(stderr,"note: using \"%s\" as universe file\n",universe_file);
        }
        else if (strcmp(argv[i],"enc") == 0)
        {
            fprintf(stderr,"info: encode mode\n");
            encode = 1;
        }
        else if (strcmp(argv[i],"msig") == 0)
        {
            fprintf(stderr,"info: msig mode\n");
            msig = 1;
        }
        else if (strcmp(argv[i],"panther") == 0)
        {
            fprintf(stderr,"info: panther only mode\n");
            panther_only = 1;
        }
    }
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
 
#ifdef L2P_USING_R

SEXP l2p_core( int Rflag, int numingenes,char *genelist[])
#else
int l2p_core( int Rflag, int numingenes,char *genelist[])
#endif
{
    char s[20000];
    struct hugo_type *hugoptr;
    struct hugo_type h;
    struct bingenetype *generec_ptr; 
    struct binpathouttype *binpathptr;
    unsigned int *usintptr = (void *)0;
    unsigned char *ucz = (void *)0; // space
    FILE *fp;
    int j;
    int hitcnt;
    int i;
    int index;
    int system_errorcode;
    int sz;
    int incnt;
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
   SEXP source;                     // 9 KEGG,REACTOME,GO,PANT(=PANTHER),PID=(pathway interaciton database)
   SEXP pathwayname;                // 10 Name of pathway
   SEXP genesinpathway;             // 11 genes_space_separated   HUGO genes from user that hit the pathway
   int maxflds = 11;
   SEXP Rret;
   SEXP cls; // R class 
   SEXP nam; // R name 
   SEXP rownam; // row names
   int protect_cnt = 0;
#endif

// printf("in l2p_core Rflag=%d numgingenes=%d, genelist=%p\n",Rflag,numingenes,genelist); fflush(stdout); 
    if (encode)
    {
        numg = 26568;  // brave new universe
        sprintf(fn_pathbin,"%s%s",exedir,"encpath.bin");
        sprintf(fn_genesbin,"%s%s",exedir,"encgenes.bin");
        sprintf(fn_spacebin,"%s%s",exedir,"encspace.bin");
    }
    else if (panther_only)
    {
        numg = 2238;
        sprintf(fn_pathbin, "%s%s",exedir,"pantherworks.bin");
        sprintf(fn_genesbin,"%s%s",exedir,"pantherworksgenes.bin");
        sprintf(fn_spacebin,"%s%s",exedir,"pantherworksspace.bin");
    }
    else if (msig)
    {
        numg = 25785;  // brave new universe
        sprintf(fn_pathbin, "%s%s",exedir,"msigpath.bin");
        sprintf(fn_genesbin,"%s%s",exedir,"msiggenes.bin");
        sprintf(fn_spacebin,"%s%s",exedir,"msigspace.bin");
    }
    else
    {
        numg = 19688;  // universe
        sprintf(fn_pathbin,"%s%s",exedir, "pathworks.bin");
        sprintf(fn_genesbin,"%s%s",exedir,"pathworksgenes.bin");
        sprintf(fn_spacebin,"%s%s",exedir,"pathworksspace.bin");
    }


    fp = fopen (fn_pathbin,"r");
    system_errorcode = errno;
    if (!fp) 
    { 
        fprintf(stderr,"Can not open \"path binary data\" - %s - errno=%d\n",fn_pathbin,system_errorcode); 
        exit(0); 
    }
    fseek(fp, 0L, SEEK_END);
    sz = ftell(fp);
    numbinpaths = sz / sizeof(struct binpathouttype); 
    rewind(fp);
    (void)fread(&binpath[0],(size_t)sz,1,fp);
    fclose(fp);
    fp = (FILE *)0;

     
    // setup_hits_parallel_to_pathways(sizeof(struct binpathouttype) * numbinpaths);

    fp = fopen (fn_genesbin,"r");
    system_errorcode = errno;
    if (!fp) { fprintf(stderr,"Can not open \"path genes data\" - %s , errno=%d\n",fn_genesbin,system_errorcode); exit(0); }
    fseek(fp, 0L, SEEK_END);
    sz = ftell(fp);
    numbingenes = (int)(sz / sizeof(struct bingenetype)); 
    rewind(fp);
    (void)fread(&bingene[0],sizeof(struct bingenetype),numbingenes,fp);
    fclose(fp);
    fp = (FILE *)0;

            // -- need hugo name access, so set up a parallel array 
    setup_hugo_parallel_to_genebin( ((size_t)(numbingenes)*sizeof(struct hugo_type)) ,numbingenes);

// -- for "spill over" , path and genes have pointers to variable sized data in space 
// these are accessed by integer offset (i.e. pointers) to null terminated strings or known sized array of ints
    fp = fopen (fn_spacebin,"r");
    system_errorcode = errno;
    if (!fp) { fprintf(stderr,"Can not open \"path space data\" - %s - errno=%d\n",fn_spacebin,system_errorcode); exit(0); }
    fseek(fp, 0L, SEEK_END);
    spacesize = ftell(fp);
    rewind(fp);
    ucz = malloc((size_t)spacesize); if (!ucz) { fprintf(stderr,"ERROR: can't malloc in l2p_core()\n"); exit(0); }
    (void)fread(ucz,spacesize,1,fp);
    fclose(fp);
    fp = (FILE *)0;

    if (user_universe_flag)
    {
fprintf(stderr,"before deal_with_universe_file\n"); fflush(stderr); 
        deal_with_universe_file(universe_file,numbingenes,numbinpaths,&numg,ucz); 
fprintf(stderr,"after deal_with_universe_file\n"); fflush(stderr); 
    }

    incnt = 0;
    s[0] = (char)0;

#ifdef L2P_USING_R
    PROTECT(Rret = Rf_allocVector(VECSXP, 11)); // a list with 11 elements
    protect_cnt++;
    for (i=0 ; i<maxflds ; i++)
    {
        PROTECT(pval=Rf_allocVector(REALSXP, numbinpaths )); 
        PROTECT(fdr=Rf_allocVector(REALSXP, numbinpaths)); 
        PROTECT(ratio=Rf_allocVector(REALSXP, numbinpaths));
        PROTECT(pathwayhitcount=Rf_allocVector(INTSXP, numbinpaths));
        PROTECT(numberofgenesinpathway=Rf_allocVector(INTSXP, numbinpaths));
        PROTECT(inputnumberofgenes=Rf_allocVector(INTSXP, numbinpaths));
        PROTECT(genesinpathwaysuniverse=Rf_allocVector(INTSXP, numbinpaths));
        PROTECT(pathwayaccessionidentifier=Rf_allocVector(STRSXP, numbinpaths));
        PROTECT(source=Rf_allocVector(STRSXP, numbinpaths));
        PROTECT(pathwayname=Rf_allocVector(STRSXP, numbinpaths));
        PROTECT(genesinpathway=Rf_allocVector(STRSXP, numbinpaths));
        protect_cnt+=11;
    }

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
                 incnt++;
                 if (generec_ptr->pathplace == -1) { fprintf(stderr,"ERROR: should not get pathplace of -1 and pathcount\n");  exit(0); }
                 if (generec_ptr->pathplace == 0) { fprintf(stderr,"ERROR: should not get pathplace of 0 and pathcount\n");  exit(0); }
                 for (i=0; i < generec_ptr->pathcount ; i++) // go through the paths for this gene, add hit to to path
                 {
                     index = *(int *)(ucz + generec_ptr->pathplace + (i*4)); // get INDEX of record in binpath[]
                     binpathptr = &binpath[index];
                     if ((binpathptr->hits) == (void *)0)
#if 1
                     {
                         binpathptr->hits = (void *)malloc((size_t)(binpathptr->numgenes+2) *4);
                         usintptr = (unsigned int *)(binpathptr->hits);
                         *(usintptr+1) = (unsigned int)(generec_ptr - &bingene[0] ); // index
                         *(usintptr) = (unsigned int) 1;
                     }
                     else
                     {
                         usintptr = (unsigned int *)(binpathptr->hits); // point to already allocated array
                         unsigned int curcount = *usintptr;
			 curcount = curcount + 1;
                         *usintptr = (unsigned int )(curcount); 
                         *(usintptr+curcount) = (unsigned int)(generec_ptr - &bingene[0] ); // index
                     }
#else
                     index = *(int *)(ucz + generec_ptr->pathplace + (i*4)); // get INDEX of record in binpath[]
                     binpathptr = &binpath[index];
                     if ((binpathptr->hits) == (void *)0)
                     {
                         binpathptr->hits = (void *)malloc((size_t)(binpathptr->numgenes+2) *8);
                         uslongintptr = (unsigned long int *)(binpathptr->hits);
                         *(uslongintptr+1) = (unsigned long int)generec_ptr;
                         *(uslongintptr) = (unsigned long int) 1;
                     }
                     else
                     {
                         uslongintptr = (unsigned long int *)(binpathptr->hits); // point to already allocated array
                         unsigned long int curcount = *uslongintptr;
			 curcount = curcount + 1;
                         *uslongintptr = (unsigned long int )(curcount); 
                         *(uslongintptr+curcount) = (unsigned long int)generec_ptr; // cast pointer as unsigned long int
                     }
#endif
                }
            }
        }
        free(h.hugo);
    }

    do_exact_and_benjamini_hochberg(numg,incnt);

    //d = exact22(n11_,n12_,n21_,n22_); note:3 hits in p53 pathway which has  40 genes ,total 29960 genes , user input 300 genes  
    // d = exact22(3,40,297,29960);
    if  (header_flag)
        printf("pval\tfdr\tratio\tpathwayhitcount\tnumberofgenesin\tpathway\tinputnumberofgens\tgenesinpathwaysuniverse\tpathwayaccessionidentifier\tsource\tpathwayname\tpathwaytype\tgenes\n");
    for (i=0 ; i<numbinpaths ; i++)
    {
            binpathptr = &binpath[i];
            if (binpathptr->hits) hitcnt = *(unsigned int *)(binpathptr->hits);
            else hitcnt = 0;
#ifdef L2P_USING_R
            REAL(pval)[i] = pathparallel[i].pval;
            REAL(fdr)[i] = pathparallel[i].fdr;
            REAL(ratio)[i] = pathparallel[i].ad;
            INTEGER(pathwayhitcount)[i] = hitcnt;
            INTEGER(numberofgenesinpathway)[i] = (binpathptr)->numgenes;
            INTEGER(inputnumberofgenes)[i] = incnt;
            INTEGER(genesinpathwaysuniverse)[i] = numg;
            SET_STRING_ELT(pathwayaccessionidentifier, i, mkChar((binpathptr)->accession));
            SET_STRING_ELT(source, i, mkChar((binpathptr)->source));
            SET_STRING_ELT(pathwayname, i, mkChar( (char *)(ucz+((binpathptr)->name ) ) ) );
            SET_STRING_ELT(pathwayname, i, mkChar( (char *)(ucz+((binpathptr)->name ) ) ) );
            SET_STRING_ELT(genesinpathway, i, mkChar( type2string((binpathptr)->type)) );

            p2 =  malloc((hitcnt+2) * 34);  // plus some extras space 
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
// printf("%d genes = %s\n",i,p2);
            }
            SET_STRING_ELT(genesinpathway, i, mkChar( p2 ) );
            if (p2) free(p2); p2 = (char *)0;
#else
            if (precise_flag == 1)
            {
                printf("%20.18f\t%20.18f\t%11.9f\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t",
                    pathparallel[i].pval, pathparallel[i].fdr, pathparallel[i].ad,
                             hitcnt, (binpathptr)->numgenes, incnt, numg,
                             (binpathptr)->accession, (binpathptr)->source,
                             (ucz+((binpathptr)->name)), type2string((binpathptr)->type));
            }
            else
            {
                printf("%11.9f\t%11.9f\t%11.9f\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t",
                    pathparallel[i].pval, pathparallel[i].fdr, pathparallel[i].ad,
                             hitcnt, (binpathptr)->numgenes, incnt, numg,
                             (binpathptr)->accession, (binpathptr)->source,
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
   SET_VECTOR_ELT(Rret, 0, pval);
   SET_VECTOR_ELT( Rret,1, fdr);

   SET_VECTOR_ELT( Rret,2, ratio);
   SET_VECTOR_ELT( Rret,3, pathwayhitcount);
   SET_VECTOR_ELT( Rret,4, numberofgenesinpathway);
   SET_VECTOR_ELT( Rret,5, inputnumberofgenes);
   SET_VECTOR_ELT( Rret,6, genesinpathwaysuniverse);
   SET_VECTOR_ELT( Rret,7, pathwayaccessionidentifier);
   SET_VECTOR_ELT( Rret,8, source);
   SET_VECTOR_ELT( Rret,9,pathwayname); 
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
   SET_STRING_ELT( nam,                     8, mkChar("source"));
   SET_STRING_ELT( nam,                9, mkChar("pathwayname"));
   SET_STRING_ELT( nam,            10, mkChar("genesinpathway"));
   namesgets(Rret, nam);

   PROTECT(rownam = allocVector(STRSXP, numbinpaths )); // row.names attribute
   protect_cnt++;
   for (i=0;i<numbinpaths;i++)
   {
       binpathptr = &binpath[i];
       SET_STRING_ELT(rownam, i, mkChar((binpathptr)->accession));
   }
   setAttrib(Rret, R_RowNamesSymbol, rownam);

   UNPROTECT(protect_cnt);
   for (i=0;i<numingenes;i++)
     free(genelist[i]);
   free(genelist);
printf("in l2p_core() return Rret = %p\n",Rret); fflush(stdout);
   return Rret;
#else
    return 0;
#endif
}


#ifdef L2P_USING_R

SEXP l2p(SEXP lst, SEXP fpath)
{
   char path[PATH_MAX];
   int lastslash;
   int i,j;
   int len;
   SEXP tmp;
   const char *f;
   int n;
   char *t;
   char **z;

   len = length(lst);
   strcpy(path,CHAR(STRING_ELT(fpath, 0)));
   z = (char **)malloc(sizeof(char *)*len);
   for (i = 0; i < len; i++) 
   {
       *(z+i) = strdup(CHAR(STRING_ELT(lst, i)));
   }
#if 0
   n = length(fpath);
printf("in l2p 2.1 n=%d\n",n); fflush(stdout); 
   t = malloc(n + 1);
   if (t != NULL) {
printf("in l2p 3\n"); fflush(stdout); 
      for (i = 0; i < n; i++) {
printf("in l2p 4\n"); fflush(stdout); 
// printf("Returned: {%s}\n", CHAR(STRING_ELT(val, 0)));
         t[i] = *CHAR(STRING_ELT(fpath, i));
      }
      t[n] = '\0';
   }
#endif

  lastslash = -1;
  for (i=0;path[i];i++)
  {
      if (path[i] == '/') lastslash = i;
  }
  if (lastslash > 0) path[lastslash] = (char)0;

// hardwire for tesing strcpy(exedir,"/Users/finneyr/R/libs/l2p/extdata/");
  strcpy(exedir,path);
  strcat(exedir,"/");

  l2p_init(1,0,(void *)0);
  return l2p_core(1,len, z);
}


#else
int l2p_run_for_C(int argc,char *argv[])
{
    l2p_init(0,argc,argv);
    get_this_executable_path(exedir,PATH_MAX);
    l2p_core(0, 0,(void *)0);
    return 0;
}

int main(int argc,char *argv[])
{
    return l2p_run_for_C(argc,argv);
}

#endif

