
/*
R interface
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

#include "pathworks.h"

#if 0
#include "small.h"
#endif

#include "big.h"

#define RDEBUG 0

extern unsigned short int pwgenes[];
extern struct smallgenetype genes[];
extern struct pwtype pws[];
extern int numpws;
extern int numgenes;
extern int numpwgenes;
extern struct smallgenetype *by_egids;  // a copy of "genes" ,but sorted by egid (entrez gene id )
extern unsigned int ingenes[MAX_INGENES];
extern unsigned int ingenecnt;

void category_code_to_string(unsigned int cat,char puthere[])
{
         if (cat & CAT_NCBI_BIOCYC)        strcpy(puthere,"BIOCYC"); 
    else if (cat & CAT_NCBI_GO    )        strcpy(puthere,"GO"); 
    else if (cat & CAT_NCBI_KEGG  )        strcpy(puthere,"KEGG"); 
    else if (cat & CAT_NCBI_PANTH )        strcpy(puthere,"PANTH"); 
    else if (cat & CAT_NCBI_PID )          strcpy(puthere,"PID");
    else if (cat & CAT_NCBI_REACTOME     ) strcpy(puthere,"REACTOME");
    else if (cat & CAT_NCBI_WikiPathways ) strcpy(puthere,"WikiPathways");
    else if (cat & CAT_MSIG_C1   ) strcpy(puthere,"C1");
    else if (cat & CAT_MSIG_C2   ) strcpy(puthere,"C2");
    else if (cat & CAT_MSIG_C3   ) strcpy(puthere,"C3");
    else if (cat & CAT_MSIG_C4   ) strcpy(puthere,"C4");
    else if (cat & CAT_MSIG_C5   ) strcpy(puthere,"C5");
    else if (cat & CAT_MSIG_C6   ) strcpy(puthere,"C6");
    else if (cat & CAT_MSIG_C7   ) strcpy(puthere,"C7");
    else if (cat & CAT_MSIG_C8   ) strcpy(puthere,"C8");
    else if (cat & CAT_MSIG_H    ) strcpy(puthere,"H");
    else if (cat & CAT_CUSTOM    ) strcpy(puthere,"CUSTOM");
    else strcpy(puthere,"UNKNOWN");
//    else if (cat & CAT_MSIG_ARCHIVED     )   strcpy(puthere,"ARCHIVED"); 
}

struct ordertype             // used for benjamini hochberg FDR 
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


void benjaminihochberg(int n,double pvals[], double returnpvals[])
{
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


#if 0
double lngamm(double z)
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

double left,right,twotail;
#endif

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

static int bitCount_unsigned(unsigned int n) 
{
    int cnt = 0;
    while (n) 
    {
        cnt += n % 2;
        n >>= 1;
    }
    return cnt;
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
            *catspat = bit | *catspat;
        else
        {
            fprintf(stderr,"ERROR: invalid category = %s\n",ts[k]);
            return 0;
        }
    }
    j = bitCount_unsigned(*catspat) ;
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


int l2pfunc_R(unsigned int *user_in_genes, unsigned int user_incnt, struct used_path_type *usedpaths,unsigned int num_used_paths,unsigned int real_universe_cnt, unsigned int *real_universe,int permute_flag,int gpcc_flag)
{
    struct used_path_type *uptr; // used path pointer 
    unsigned int ui_uk;
    unsigned int ui_ej;
    unsigned int i;
    unsigned int j,k,ll;
    int ret = 0;

   // this shuts up "unused" compiler warning
#define UNUSED(x) (void)(x)
UNUSED(permute_flag);

    
// fprintf(stderr,"l2pfunc_R : %p user_incnt=%u %p %u %u %p %u %u\n", user_in_genes, user_incnt, usedpaths,num_used_paths,real_universe_cnt, real_universe,permute_flag,gpcc_flag);

    for (i=0 ; i<num_used_paths ; i++)
    {
        uptr = (usedpaths+i);
        if (!uptr->egids) continue;
        j = k = ll = 0;
        while ((j<uptr->numfixedgenes) && (k < user_incnt))
        {
            ui_ej = *(uptr->egids+j);
            ui_uk = *(user_in_genes+k);
            if (ui_ej == ui_uk)
            {
                *((uptr->genehits) + (ll++)) = ui_uk; // remember, because need to print out later
                k++;
                j++;
                continue;
            }
            else if (ui_ej < ui_uk) j++;
            else                    k++;
        }
        uptr->hitcnt = ll;
    }
    if (gpcc_flag == 1)
    {
//fprintf(stderr,"gpccdbg in GPCC flag, user_incnt=%d\n",user_incnt); fflush(stderr); 
        for (i=0;i<user_incnt;i++)
            ingenes[i] = *(user_in_genes+i);
        ingenecnt = user_incnt;
// fprintf(stderr,"gpccdbg before GPCC %p %u %u %p\n",usedpaths,num_used_paths,real_universe_cnt,real_universe); fflush(stderr); 
        GPCC(usedpaths,num_used_paths,real_universe_cnt,real_universe);
//fprintf(stderr,"gpccdbg after GPCC\n"); fflush(stderr); 
    }
    return ret;
}


SEXP l2p(SEXP lst, SEXP categories, SEXP universe, SEXP custompws, SEXP customfn, SEXP universefn, SEXP permute_arg, SEXP oneside_arg, SEXP gpcc_arg, 
    SEXP legacy_arg, SEXP extra_arg)
{
    char tmps2[PATH_MAX];  // temp string for "category" 
    char tmps[512];
    char universe_file[PATH_MAX];
    char custom_file[PATH_MAX];
    unsigned int i,k,k2;
    unsigned int j;
    unsigned int len_of_user_pws = 0;  // length of list of lists , from Rf_length() which return R_len_t ( which in int? )
    unsigned int len_of_vector = 0;
    struct used_path_type *uptr;
    unsigned int *user_in_univ_ptr = (unsigned int *)0;
    unsigned int *in_universe_original = (unsigned int *)0;
    unsigned int in_universe_cnt = 0;
    unsigned int *real_universe = (unsigned int *)0;
    unsigned int real_universe_cnt = 0;
    int extra_flag = 0;
    int legacy_flag = 0;
    int permute_flag = 0;
    int gpcc_flag = 0;
    int user_universe_flag = 0;
    unsigned int user_incnt = 0;
    int oneside = 0;
    unsigned int num_used_paths = 0;
    unsigned int len = 0;
    unsigned int *user_in_genes = (unsigned int *)0;
    unsigned int *user_in_genes_original = (unsigned int *)0;
    char *p2 = (char *)0;
    char *zz = (char *)0;
    unsigned int catspat = 0;   // categories pattern   , if 0, use all
    unsigned int ui = 0;
    struct used_path_type *u = (struct used_path_type *)0;
    int protect_cnt = 0;
    unsigned int maxflds = 13;
    struct custom_type *mycustompw = (struct custom_type *)0;
    struct custom_type *mycustompwptr = (struct custom_type *)0;
    double fdr_for_output;
    SEXP list = (SEXP)0;
    SEXP pval = (SEXP)0;     // 1 
    SEXP fdr = (SEXP)0;      // 2
    SEXP GPCC_pval = (SEXP)0; 
    SEXP GPCC_FDR = (SEXP)0; 
    SEXP number_genes_in_pathway = (SEXP)0; 
    SEXP percent_hits = (SEXP)0; 
    SEXP enrichment_score = (SEXP)0;  // 3 if positive, genes are OVER REPRESENTED, if negative genes are UNDER REPRESENTED
    SEXP percent_gene_hits_per_pathway = (SEXP)0;  // 4
    SEXP number_hits = (SEXP)0;       // 5 number of genes hit in pathway
    SEXP number_misses = (SEXP)0;     // 6 number of genes in the pathway
    SEXP number_user_genes = (SEXP)0;       // 7 total count of user genes (user input)
    SEXP total_genes_minus_input = (SEXP)0;    // 8 total number of unique genes in all pathways
    SEXP pathway_id = (SEXP)0; // 9 canonical accession ( if available, otherwise assigned by us )
    SEXP pathway_type = (SEXP)0;
    SEXP category = (SEXP)0;                   // 10 KEGG,REACTOME,GO,PANT(=PANTHER),PID=(pathway interaciton database)
    SEXP pathway_name = (SEXP)0;                // 11 Name of pathway
    SEXP genesinpathway = (SEXP)0;             // 12 genes_space_separated   HUGO genes from user that hit the pathway
    SEXP permute_pval = (SEXP)0;     
    SEXP universe_count = (SEXP)0; 
    SEXP user_input_count = (SEXP)0; 
    SEXP odds_ratio = (SEXP)0; 
    SEXP Rret = (SEXP)0;
    SEXP cls = (SEXP)0;       // R class 
    SEXP nam = (SEXP)0;       // R name 
    SEXP rownam = (SEXP)0;    // row names


    (void)setup_by_egids();

    if (Rf_isNull(permute_arg)) { permute_flag = 0; } else { permute_flag = asInteger(permute_arg); }
    if (Rf_isNull(extra_arg))   { extra_flag = 0; }   else { extra_flag = asInteger(extra_arg); }
    if (Rf_isNull(legacy_arg))  { legacy_flag = 0; }  else { legacy_flag = asInteger(legacy_arg); }
    if (Rf_isNull(gpcc_arg))    { gpcc_flag = 0; }    else { gpcc_flag = asInteger(gpcc_arg); }

    len = length(lst);
    user_in_genes = (unsigned int *)malloc(sizeof(unsigned int)*len); // remember to free this
    if (!user_in_genes)
        return (SEXP) -1; // why not 0 ?
    user_in_genes_original = (unsigned int *)malloc(sizeof(unsigned int)*len); // remember to free this
    if (!user_in_genes_original)
    {
        free(user_in_genes);
        user_in_genes = (unsigned int *)0;
        return (SEXP) -1; // why not 0 ?
    }
    for (k = i = 0; i < len; i++)
    {
        strcpy(tmps,CHAR(STRING_ELT(lst, i)));
        ui = hugo2egid(tmps);
        if (ui != UINT_MAX) { *(user_in_genes_original+k) = ui; k++; }
    }
#if RADIX
    radix_ui(user_in_genes_original,k);
#else
    qsort(user_in_genes_original,k,sizeof(unsigned int),cmp_ui);
#endif

// fprintf(stderr,"adbg legacy_flag is %d",legacy_flag); 
// fprintf(stderr,"gpccdbg , in l2p, gpcc_arg=%p\n",gpcc_arg); 

    if (Rf_isNull(oneside_arg)) 
    {
        oneside = 0;
    }
    else 
    { 
        oneside = asInteger(oneside_arg);
        if ((oneside < 0) || (oneside > 2)) oneside = 0;      // just give them twosided fisher's exact test
    }

    universe_file[0] = custom_file[0] = tmps[0] = tmps2[0] = (char)0;
    if (Rf_isNull(custompws))      //  if (custompws == (SEXP)0)
    {
//        custom_flag = 0; not used
        len_of_user_pws = 0;
    }
    else
    {
/* struct custom_type { char *name; char *optional; unsigned int numgenes; unsigned int *genes; }; */
//        custom_flag = 1; not used
        len_of_user_pws = (int)length(custompws);
        mycustompw = (struct custom_type *)malloc(sizeof(struct custom_type )*len_of_user_pws); // remember to free this
        memset(mycustompw,0,sizeof(struct custom_type )*len_of_user_pws);
        for (i=0;i<len_of_user_pws;i++)
        {
            mycustompwptr = mycustompw + i;
            list = VECTOR_ELT(custompws, i);
            len_of_vector = length(list);
            mycustompwptr->genes = (unsigned int *)malloc(sizeof(unsigned int)*len_of_vector);
            memset( mycustompwptr->genes,0,(sizeof(unsigned int)*len_of_vector));
            mycustompwptr->numgenes = 0;   // will set this properly later , initialize to zero for now
            for (j=0 ; j<len_of_vector ; j++)
            {
                memset(tmps2,0,sizeof(tmps2));
                strncpy(tmps2,CHAR(STRING_ELT(list, j)),PATH_MAX-2);
                if (j == 0) mycustompwptr->name = strdup(tmps2);    // name. check: must free this!
                else if (j == 1) mycustompwptr->optional = strdup(tmps2); //  gmt says second field meaning is unspecified. could this be a URL someday?
                else 
                {
                    ui = hugo2egid(tmps2);
// fprintf(stderr,"custgene is i=%d j=%d %u %s\n",i,j,ui,tmps2);  fflush(stderr); 
                    if (ui != UINT_MAX) 
                    { 
                        *(mycustompwptr->genes + j - 2) = ui;  // nb:-2 . remember fields 1 and 2 are not the genes, fields3 and greater are
                        mycustompwptr->numgenes++;
                    }
                }
            }
#if RADIX
            radix_ui(mycustompwptr->genes,mycustompwptr->numgenes);
#else
            qsort(mycustompwptr->genes,mycustompwptr->numgenes,sizeof(unsigned int),cmp_ui);
#endif
        }
    }
    if (Rf_isNull(universe)) //   if (universe == (SEXP)0)
        user_universe_flag = 0;
    else
        user_universe_flag = 1;
    if (Rf_isNull(categories ))
    {
        category_set_all(&catspat);
// fprintf(stderr,"no parsecats after category_set_all: %u (null categories) \n",catspat);
    }
    else
    {
        if (isVector(categories))  // it's really a dam string
        {
             tmps[0] = (char)0;
             len = length(categories);
             for (k2 = 0; k2 < len; k2++)
             {
                 if (k2) strcat(tmps,",");
                 strncpy(tmps2,CHAR(STRING_ELT(categories, k2)),PATH_MAX-2);
                 strcat(tmps,tmps2);
             }
// fprintf(stderr,"before parsecats : %s\n",tmps);
            (void)parsecats(tmps,&catspat);     // set catpats
        }
        else 
        {
fprintf(stderr,"NOTE: Not sure what's up. Can't parse categroies \n");  fflush(NULL);
        }
    }
categories_pattern_to_strings(catspat,tmps);
// fprintf(stderr,"GOTEST rpfdbg categories here 9, catspats=%x = \"%s\"\n",catspat,tmps);   fflush(NULL); 

    if (Rf_isNull(customfn)) {} else strncpy(custom_file,CHAR(STRING_ELT(customfn, 0)),PATH_MAX-2);
    if (Rf_isNull(universefn)) {} else strncpy(universe_file,CHAR(STRING_ELT(universefn, 0)),PATH_MAX-2);

    for (j=i=0;i<k;i++)
    {                // de duplicate list 
        if (i > 0) 
        {
            if ( *(user_in_genes_original+i) == *(user_in_genes_original+i-1) )
               continue;
            *(user_in_genes+j) = *(user_in_genes_original+i);
            j++;
        }
        else 
        {
            *(user_in_genes+j) = *(user_in_genes_original+i);
            j++;
        }
    }
    user_incnt = (unsigned int )j;

    len = length(universe);
    if (user_universe_flag == 1)
    {
        user_in_univ_ptr = (unsigned int *)malloc(sizeof(unsigned int)*len);        // remember to free this
        if (!user_in_univ_ptr) return (SEXP) -1;
        in_universe_original = (unsigned int *)malloc(sizeof(unsigned int)*len);    // remember to free this
        if (!in_universe_original) 
	{
            free(user_in_univ_ptr);
            return (SEXP) -1;
	}
        for (k = i = 0; i < len; i++)
        {
            strcpy(tmps,CHAR(STRING_ELT(universe, i)));
            ui = hugo2egid(tmps);
            if (ui != UINT_MAX) { *(in_universe_original+k) = ui; k++; }
        }
#if RADIX
        radix_ui( in_universe_original,k);
#else
        qsort(in_universe_original,k,sizeof(unsigned int),cmp_ui);
#endif
        for (j=i=0;i<k;i++)
        {  // de duplicate universe 
            if (i > 0) 
            {
                if ( *(in_universe_original+i) == *(in_universe_original+i-1) )
                   continue;
                *(user_in_univ_ptr+j) = *(in_universe_original+i);
                j++;
            }
            else 
            {
                *(user_in_univ_ptr+j) = *(in_universe_original+i);
                j++;
            }
        }
        in_universe_cnt = j;
    }
    if (in_universe_original) { free(in_universe_original); in_universe_original = (void *)0; }

// xxx fix
// fprintf(stderr,"GOTEST before setup_used_paths\n");   fflush(NULL); 
    u = setup_used_paths(&num_used_paths, catspat,universe_file,in_universe_cnt,user_in_univ_ptr,custom_file,&real_universe_cnt,&real_universe,len_of_user_pws,mycustompw);
// fprintf(stderr,"GOTEST after setup_used_paths\n");   fflush(NULL); 
//fprintf(stderr,"in gpccdbg  after setup_used_path cats=%x after setup_used_paths() \n",catspat);  fflush(stderr);
// NO, freed in setup_used_paths    if (user_in_univ_ptr) { free(user_in_univ_ptr); user_in_univ_ptr = (void *)0; }
    
#if 0
 FILE *fptmp=fopen("tmp.uig.preint","w");
 for (ui=0;ui<user_incnt;ui++)
 {
 char *zz;
 zz = egid2hugo(*(user_in_genes+ui));
 fprintf(fptmp,"%s\n",zz);
 }
 fclose(fptmp);
#endif
    
// xxx
    unsigned int *tmp_in_genes;
    unsigned int uj,tmpegid;
    unsigned int *uiptr;
    tmp_in_genes = (unsigned int *)malloc(sizeof(unsigned int)*user_incnt);    // remember to free this
    for (ui=uj=0;ui<user_incnt;ui++)
    {
         tmpegid = *(user_in_genes+ui);
         uiptr = (unsigned int *)bsearch(&tmpegid,real_universe,real_universe_cnt,sizeof(unsigned int),cmp_ui);
         if (uiptr)
         {
             *(tmp_in_genes+uj) = *uiptr;
             uj++;
         }
    }
    user_incnt = uj;
    memcpy(user_in_genes,tmp_in_genes,user_incnt*sizeof(unsigned int));
    free(tmp_in_genes);

#if 0
 FILE *fptmp2=fopen("tmp.uig.afterint","w");
 for (ui=0;ui<user_incnt;ui++)
 {
 char *zzz;
 zzz = egid2hugo(*(user_in_genes+ui));
 fprintf(fptmp2,"%s\n",zzz);
 }
 fclose(fptmp2);
#endif

// fprintf(stderr,"gpccdbg in l2p, before l2pfunc_R gpcc_flag=%d user_incnt=%u\n",gpcc_flag,user_incnt); fflush(stderr); 
    GetRNGstate();
    (void)l2pfunc_R(user_in_genes,user_incnt,u,num_used_paths,real_universe_cnt,real_universe,permute_flag,gpcc_flag);
    PutRNGstate();
//fprintf(stderr,"gpccdbg in l2p, after l2pfunc_R\n"); fflush(stderr); 
    if (gpcc_flag)
    {
// fprintf(stderr,"gpccdbg: in l2p - in gpcc_flag output\n"); fflush(stderr); 
#if 0
    double gpcc_p; double gpcc_fdr; double enrichment_score;                  // ratio
Columns that will be output with L2P:
     8	"genes_in_pathway"
     9	"percent_hits"
    10	"pathway_id"
    11	"category"
    12	"hits"
    13	"genesinpathway"

8 Number_genes_in_pathway        Number of genes in the pathway (size of pathway)*
9 Percent_hits          Number_hits/total genes in pathway* <- tie-breaker #3
10 Pathway_id           Canonical accession ( if available, otherwise assigned by us )
11 Category             KEGG,REACTOME,GO,PANT(=PANTHER),PID=(pathway interaction database)  note: was ,source,
12 Pathway_type         Functional_set,pathway,structural_complex,custom
13 Genes_in_pathway     HUGO genes from user that hit the pathway

Optional Columns: This is for Rich to put into a flag where optional columns can be added to the L2P output (Optionalcols = 1, by default Optionalcols = 0)

* A – universe – (DEG count + pathway count) 
B – pathway count – number of hits XXX already
C – total DEG count - number of hits
D – number of hits
Odds_ratio              odds ratio
Sum_of_pathways         Sum of unique pathways where pathway genes are also found

universe
inputlistcount
or

#endif
        if (extra_flag == 0) maxflds = 13; else maxflds = 16; 
        PROTECT(Rret = Rf_allocVector(VECSXP, maxflds)); // a list
        protect_cnt++;
// fprintf(stderr,"gpccdbg: in l2p - in gpcc_flag output 1\n"); fflush(stderr); 
        PROTECT(pathway_name=Rf_allocVector(STRSXP, num_used_paths));              // 1 
        PROTECT(enrichment_score=Rf_allocVector(REALSXP, num_used_paths));         // 2 
        PROTECT(GPCC_pval=Rf_allocVector(REALSXP, num_used_paths ));               // 3 
        PROTECT(GPCC_FDR=Rf_allocVector(REALSXP, num_used_paths ));                // 4 
        PROTECT(pval=Rf_allocVector(REALSXP, num_used_paths ));                    // 5 
        PROTECT(fdr=Rf_allocVector(REALSXP, num_used_paths));                      // 6 
        PROTECT(number_hits=Rf_allocVector(INTSXP, num_used_paths));               // 7 
        PROTECT(number_genes_in_pathway =Rf_allocVector(INTSXP, num_used_paths)); // 8 
        PROTECT(percent_hits=Rf_allocVector(REALSXP, num_used_paths));              // 9 
        PROTECT(pathway_id=Rf_allocVector(STRSXP, num_used_paths));                // 10
        PROTECT(category=Rf_allocVector(STRSXP, num_used_paths));    // 11 is "l2p internal" an integer, but convert to string
        PROTECT(pathway_type=Rf_allocVector(STRSXP, num_used_paths)); // 12
        PROTECT(genesinpathway=Rf_allocVector(STRSXP, num_used_paths)); // 13
        if (extra_flag == 1)
        {
             PROTECT(universe_count=Rf_allocVector(INTSXP, num_used_paths));   // 14
             PROTECT(user_input_count=Rf_allocVector(INTSXP, num_used_paths));   // 15
             PROTECT(odds_ratio=Rf_allocVector(REALSXP, num_used_paths));   // 16
        }
        protect_cnt += maxflds;
// fprintf(stderr,"GOTEST here 1\n");  fflush(stderr); 
        for (ui=0 ; ui<num_used_paths ; ui++)
        {
            uptr = (u+ui);
            SET_STRING_ELT(pathway_name, ui, mkChar(uptr->name) ); // 1
            REAL(enrichment_score)[ui] = uptr->enrichment_score;   // 2 
// if (strcmp(uptr->acc,"ko00010") == 0) { fprintf(stderr,"adbg : l2p (rinterface.c)  ko00010 enrichment_score is %f \n",uptr->enrichment_score); }
            REAL(GPCC_pval)[ui] = uptr->gpcc_p;                    // 3 
            REAL(GPCC_FDR)[ui] = uptr->gpcc_fdr;                   // 4 
            REAL(pval)[ui] =  uptr->pval;                          // 5 
            REAL(fdr)[ui] =  uptr->fdr;                            // 6 
            INTEGER(number_hits)[ui] = uptr->hitcnt;               // 7 
            INTEGER(number_genes_in_pathway)[ui] = uptr->numfixedgenes;                    // 8 
            if ( (uptr->a == 0) && (uptr->b == 0) ) REAL(percent_hits)[ui] = (double)0;
            else 
            {
                  REAL(percent_hits)[ui] = (double)((double)uptr->hitcnt / (double)(uptr->numfixedgenes));
            }
            SET_STRING_ELT(pathway_id, ui, mkChar(uptr->acc) );
    
            category_code_to_string( uptr->category , tmps);
// fprintf(stderr,"GOTEST tmps=%s\n",tmps);  fflush(stderr); 
            SET_STRING_ELT(category, ui, mkChar(tmps));
    
            INTEGER(number_hits)[ui] = uptr->hitcnt; 
            p2 =  malloc((uptr->hitcnt +2) * 34);  // plus some extra space
            memset(p2,0,(uptr->hitcnt +2) * 34);
            for (j=0 ; j<uptr->hitcnt ; j++)
            {
               zz = egid2hugo(*((uptr->genehits)+j));
               if (j == 0) sprintf(tmps,"%s",zz); 
               else sprintf(tmps," %s",zz); 
               strcat(p2,tmps);
            }
            SET_STRING_ELT( genesinpathway, ui, mkChar( p2 ) );
            if (p2) { free(p2); p2 = (char *)0; }
            if (extra_flag == 1)
            {
                 INTEGER(universe_count)[ui] = uptr->d  ;   // 14
                 INTEGER(user_input_count)[ui] = uptr->c  ;   // 15
                 REAL(odds_ratio)[ui] = uptr->gpcc_OR  ;   // 16
// fprintf(stderr,"rpfdbg OR = %f\n",uptr->gpcc_OR);
            }
        }
// fprintf(stderr,"gpccdbg: in l2p - in gpcc_flag output 3\n"); fflush(stderr); 
#if 0
1 Pathway_name          Name of pathway
2 Enrichment_score      ((number_hits /(number_hits+number_misses)) - (number_user_genes/(number_user_genes+total_gens_minus_input))) * 100 <- tie-breaker #2
3 GPCC_pval             Gpcc (gene pathway count correction) pval <- sort on this first (#1)
4 GPCC_FDR              False discovery rate: benjamini hochberg (of GPCC p-value) <- this is what we recommend users use
5 FET_pval              Fishers exact p-value (legacy)
6 FET_FDR               False discovery rate (FET)
7 Number_hits           Number of genes hit in pathway
8 Number_genes_in_pathway        Number of genes in the pathway (size of pathway)*
9 Percent_hits          Number_hits/total genes in pathway* <- tie-breaker #3
10 Pathway_id           Canonical accession ( if available, otherwise assigned by us )
11 Category             KEGG,REACTOME,GO,PANT(=PANTHER),PID=(pathway interaction database)  -- was ,source,
12 Pathway_type         Functional_set,pathway,structural_complex,custom
13 Genes_in_pathway     HUGO genes from user that hit the pathway
Optional Columns: This is for Rich to put into a flag where optional columns can be added to the L2P output (Optionalcols = 1, by default Optionalcols = 0)
A – universe – (DEG count + pathway count) 
B – pathway count – number of hits
C – total DEG count - number of hits
D – number of hits
Odds_ratio              odds ratio
Sum_of_pathways         Sum of unique pathways where pathway genes are also found
#endif
        SET_VECTOR_ELT( Rret,0, pathway_name);
        SET_VECTOR_ELT( Rret,1, enrichment_score);
        SET_VECTOR_ELT( Rret,2, GPCC_pval);
        SET_VECTOR_ELT( Rret,3, GPCC_FDR);
        SET_VECTOR_ELT( Rret,4, pval);
        SET_VECTOR_ELT( Rret,5, fdr);
        SET_VECTOR_ELT( Rret,6, number_hits);
        SET_VECTOR_ELT( Rret,7, number_genes_in_pathway);
        SET_VECTOR_ELT( Rret,8, percent_hits);
        SET_VECTOR_ELT( Rret,9, pathway_id);
        SET_VECTOR_ELT( Rret,10,category);
        SET_VECTOR_ELT( Rret,11,number_hits);
        SET_VECTOR_ELT( Rret,12,genesinpathway);
        if (extra_flag == 1)
        {
            SET_VECTOR_ELT( Rret,13,universe_count);
            SET_VECTOR_ELT( Rret,14,user_input_count);
            SET_VECTOR_ELT( Rret,15,odds_ratio);
        }
    
        PROTECT(cls = allocVector(STRSXP, 1)); // class attribute
        protect_cnt++;
    
        SET_STRING_ELT(cls, 0, mkChar("data.frame"));
        classgets(Rret, cls);
    
        PROTECT(nam = allocVector(STRSXP, maxflds));     // names attribute (column names)
        protect_cnt++;

// xxx
        SET_STRING_ELT( nam, 0,  mkChar("pathway_name"));
        SET_STRING_ELT( nam, 1,  mkChar("enrichment_score"));
        SET_STRING_ELT( nam, 2,  mkChar("GPCC_pval"));
        SET_STRING_ELT( nam, 3,  mkChar("GPCC_FDR"));
        SET_STRING_ELT( nam, 4,  mkChar("pval"));
        SET_STRING_ELT( nam, 5,  mkChar("fdr"));
        SET_STRING_ELT( nam, 6,  mkChar("number_hits"));
        SET_STRING_ELT( nam, 7,  mkChar("genes_in_pathway"));
        SET_STRING_ELT( nam, 8,  mkChar("percent_hits"));
        SET_STRING_ELT( nam, 9,  mkChar("pathway_id"));
        SET_STRING_ELT( nam, 10, mkChar("category"));
        SET_STRING_ELT( nam, 11, mkChar("hits"));
        SET_STRING_ELT( nam, 12, mkChar("genesinpathway"));
        if (extra_flag == 1)
        {
            SET_STRING_ELT( nam, 13, mkChar("universe_count"));
            SET_STRING_ELT( nam, 14, mkChar("user_input_count"));
            SET_STRING_ELT( nam, 15, mkChar("odds_ratio"));
        }
    
        namesgets(Rret, nam);
    
        PROTECT(rownam = allocVector(STRSXP, num_used_paths )); // row.names attribute
        protect_cnt++;
        for (i=0 ; i<num_used_paths ; i++)
        {
            uptr = (u+i);
            SET_STRING_ELT(rownam, i, mkChar( (char *)uptr->acc) );
        }
        setAttrib(Rret, R_RowNamesSymbol, rownam);
    }
    else if (permute_flag)
    {
        do_just_bh(user_incnt, u,num_used_paths,real_universe_cnt);
        PROTECT(Rret = Rf_allocVector(VECSXP, 15)); // a list with 15 elements
        protect_cnt++;
        maxflds = 15;
// fprintf(stderr,"pmf 1\n"); fflush(stderr);
// xxx
        PROTECT(pathway_name=Rf_allocVector(STRSXP, num_used_paths));
        PROTECT(pval=Rf_allocVector(REALSXP, num_used_paths ));
        PROTECT(fdr=Rf_allocVector(REALSXP, num_used_paths));
        PROTECT(odds_ratio=Rf_allocVector(REALSXP, num_used_paths ));
        PROTECT(permute_pval=Rf_allocVector(REALSXP, num_used_paths ));
        PROTECT(enrichment_score=Rf_allocVector(REALSXP, num_used_paths));
        PROTECT(percent_gene_hits_per_pathway =Rf_allocVector(REALSXP, num_used_paths));
        PROTECT(number_hits=Rf_allocVector(INTSXP, num_used_paths));
        PROTECT(number_misses=Rf_allocVector(INTSXP, num_used_paths));
        PROTECT(number_user_genes=Rf_allocVector(INTSXP, num_used_paths));
        PROTECT(total_genes_minus_input=Rf_allocVector(INTSXP, num_used_paths));
        PROTECT(pathway_id=Rf_allocVector(STRSXP, num_used_paths));
        PROTECT(category=Rf_allocVector(STRSXP, num_used_paths));// is "l2p internal" an integer, but convert to string
        PROTECT(pathway_type=Rf_allocVector(STRSXP, num_used_paths));
        PROTECT(genesinpathway=Rf_allocVector(STRSXP, num_used_paths));
        protect_cnt += 15;
        for (ui=0 ; ui<num_used_paths ; ui++)
        {
            uptr = (u+ui);
    
            REAL(pval)[ui] = uptr->pval;                           // 0 
            fdr_for_output = uptr->fdr;
            REAL(fdr)[ui] =   fdr_for_output;                      // 1 
            REAL(pval)[ui] = uptr->OR;                             // 2 
            REAL(pval)[ui] = uptr->permute_pval;                          // 3 
            REAL(enrichment_score)[ui] = uptr->enrichment_score;   // 4 
            if ( (uptr->a == 0) && (uptr->b == 0) )
                REAL(percent_gene_hits_per_pathway)[ui] = (double)0;
            else
                REAL(percent_gene_hits_per_pathway)[ui] = (double)uptr->a/(double)((uptr->a)+(uptr->b)); // 5
            INTEGER(number_hits)[ui] =  uptr->a;                   // 6 
            INTEGER(number_misses)[ui] =  uptr->b;                 // 7 
            INTEGER(number_user_genes)[ui] = uptr->c;              // 8 
            INTEGER(total_genes_minus_input)[ui] = uptr->d;        // 9 
    
            SET_STRING_ELT(pathway_id, ui, mkChar(uptr->acc) );    // 10 
    
            category_code_to_string( uptr->category , tmps);   // 11 
            SET_STRING_ELT(category, ui, mkChar(tmps));        // 12 
    
            SET_STRING_ELT(pathway_name, ui, mkChar(uptr->name) ); // 13
            p2 =  malloc((uptr->hitcnt +2) * 34);  // plus some extra space
            memset(p2,0,(uptr->hitcnt +2) * 34);
            for (j=0 ; j<uptr->hitcnt ; j++)
            {
               zz = egid2hugo(*((uptr->genehits)+j));
               if (j == 0)
                   sprintf(tmps,"%s",zz); 
               else
                   sprintf(tmps," %s",zz); 
               strcat(p2,tmps);
            }
            SET_STRING_ELT( genesinpathway, ui, mkChar( p2 ) );
            if (p2) { free(p2); p2 = (char *)0; }
        }
//fprintf(stderr,"pmf 2\n"); fflush(stderr);
    
        SET_VECTOR_ELT( Rret,0,pathway_name);
        SET_VECTOR_ELT( Rret,1, category);
        SET_VECTOR_ELT( Rret,2, pval);
        SET_VECTOR_ELT( Rret,3, fdr);
        SET_VECTOR_ELT( Rret,4, odds_ratio);
        SET_VECTOR_ELT( Rret,5, permute_pval);
        SET_VECTOR_ELT( Rret,6, enrichment_score);
        SET_VECTOR_ELT( Rret,7, percent_gene_hits_per_pathway);
        SET_VECTOR_ELT( Rret,8, number_hits);
        SET_VECTOR_ELT( Rret,9, number_misses);
        SET_VECTOR_ELT( Rret,10, number_user_genes);
        SET_VECTOR_ELT( Rret,11, total_genes_minus_input);
        SET_VECTOR_ELT( Rret,12, pathway_id);
        SET_VECTOR_ELT( Rret,13, pathway_type);
        SET_VECTOR_ELT( Rret,14, genesinpathway);
    
        PROTECT(cls = allocVector(STRSXP, 1)); // class attribute
        protect_cnt++;
    
        SET_STRING_ELT(cls, 0, mkChar("data.frame"));
        classgets(Rret, cls);
    
        PROTECT(nam = allocVector(STRSXP, maxflds));     // names attribute (column names)
        protect_cnt++;
/*
     1	"pathway_name"
     2	"category"
     3	"pval"
     4	"fdr"
     5	"enrichment_score"
     6	"percent_gene_hits_per_pathway"
     7	"number_hits"
     8	"number_misses"
     9	"number_user_genes"
    10	"total_genes_minus_input"
    11	"pathway_id"
    12	"pathway_type"
    13	"genesinpathway"
*/
    
        SET_STRING_ELT( nam, 0,  mkChar("pathway_name"));
        SET_STRING_ELT( nam, 1,  mkChar("category"));
// "pathway_name"	"category"	"pval"	"fdr"	"enrichment_score"	"percent_gene_hits_per_pathway"	"number_hits"	"number_misses"	"number_user_genes"	"total_genes_minus_input"	"pathway_id"	"pathway_type"	"genesinpathway"
// "ko00010"	"Glycolysis / Gluconeogenesis"	"KEGG"	0.642755317070867	1	-0.0270562770562772	0.0238095238095238	1	41	89	3607	"ko00010"	""	"PGAM2"
        SET_STRING_ELT( nam, 2,  mkChar("pval"));
        SET_STRING_ELT( nam, 3,  mkChar("fdr"));
        SET_STRING_ELT( nam, 4,  mkChar("OR"));
        SET_STRING_ELT( nam, 5,  mkChar("permute_pval"));
        SET_STRING_ELT( nam, 6,  mkChar("enrichment_score"));
        SET_STRING_ELT( nam, 7,  mkChar("percent_gene_hits_per_pathway"));
        SET_STRING_ELT( nam, 8,  mkChar("number_hits"));
        SET_STRING_ELT( nam, 9,  mkChar("number_misses"));
        SET_STRING_ELT( nam, 10, mkChar("number_user_genes"));
        SET_STRING_ELT( nam, 11, mkChar("total_genes_minus_input"));
        SET_STRING_ELT( nam, 12, mkChar("pathway_id"));
        SET_STRING_ELT( nam, 13, mkChar("pathway_type"));
        SET_STRING_ELT( nam, 14, mkChar("genesinpathway"));
    
        namesgets(Rret, nam);
    
        PROTECT(rownam = allocVector(STRSXP, num_used_paths )); // row.names attribute
        protect_cnt++;
        for (i=0 ; i<num_used_paths ; i++)
        {
            uptr = (u+i);
            SET_STRING_ELT(rownam, i, mkChar( (char *)uptr->acc) );
        }
        setAttrib(Rret, R_RowNamesSymbol, rownam);
    }
    else if (legacy_flag == 1)
    {
// fprintf(stderr,"adbg in legacy_flag\n");  fflush(NULL); 
        do_pvals_and_bh(user_incnt, u, num_used_paths,real_universe_cnt,oneside);

        PROTECT(Rret = Rf_allocVector(VECSXP, 13)); // a list with 13 elements
        protect_cnt++;
        PROTECT(pathway_name=Rf_allocVector(STRSXP, num_used_paths));
        PROTECT(pval=Rf_allocVector(REALSXP, num_used_paths ));
        PROTECT(fdr=Rf_allocVector(REALSXP, num_used_paths));
        PROTECT(enrichment_score=Rf_allocVector(REALSXP, num_used_paths));
        PROTECT(percent_gene_hits_per_pathway =Rf_allocVector(REALSXP, num_used_paths));
        PROTECT(number_hits=Rf_allocVector(INTSXP, num_used_paths));
        PROTECT(number_misses=Rf_allocVector(INTSXP, num_used_paths));
        PROTECT(number_user_genes=Rf_allocVector(INTSXP, num_used_paths));
        PROTECT(total_genes_minus_input=Rf_allocVector(INTSXP, num_used_paths));
        PROTECT(pathway_id=Rf_allocVector(STRSXP, num_used_paths));
        PROTECT(category=Rf_allocVector(STRSXP, num_used_paths));// is "l2p internal" an integer, but convert to string
        PROTECT(pathway_type=Rf_allocVector(STRSXP, num_used_paths));
        PROTECT(genesinpathway=Rf_allocVector(STRSXP, num_used_paths));
        protect_cnt += 13;
        for (ui=0 ; ui<num_used_paths ; ui++)
        {
            uptr = (u+ui);
    
            REAL(pval)[ui] = uptr->pval;                           // 0 
            fdr_for_output = uptr->fdr;
            REAL(fdr)[ui] =   fdr_for_output;                      // 1 
            REAL(enrichment_score)[ui] = uptr->enrichment_score;                 // 2 
            if ( (uptr->a == 0) && (uptr->b == 0) )
                REAL(percent_gene_hits_per_pathway)[ui] = (double)0;
            else
                REAL(percent_gene_hits_per_pathway)[ui] = (double)((double)uptr->a/(double)((uptr->a)+(uptr->b))*100.0); // 3
            INTEGER(number_hits)[ui] =  uptr->a;                   // 4 
            INTEGER(number_misses)[ui] =  uptr->b;                 // 5 
            INTEGER(number_user_genes)[ui] = uptr->c;              // 6 
            INTEGER(total_genes_minus_input)[ui] = uptr->d;        // 7 
    
            SET_STRING_ELT(pathway_id, ui, mkChar(uptr->acc) );    // 8 
    
            category_code_to_string( uptr->category , tmps);   // 9 
            SET_STRING_ELT(category, ui, mkChar(tmps));        // 10 
    
            SET_STRING_ELT(pathway_name, ui, mkChar(uptr->name) ); // 11
            p2 =  malloc((uptr->hitcnt +2) * 34);  // plus some extra space
            memset(p2,0,(uptr->hitcnt +2) * 34);
            for (j=0 ; j<uptr->hitcnt ; j++)
            {
               zz = egid2hugo(*((uptr->genehits)+j));
               if (j == 0)
                   sprintf(tmps,"%s",zz); 
               else
                   sprintf(tmps," %s",zz); 
               strcat(p2,tmps);
            }
            SET_STRING_ELT( genesinpathway, ui, mkChar( p2 ) );
            if (p2) { free(p2); p2 = (char *)0; }
        }
    
        SET_VECTOR_ELT( Rret,0,pathway_name);
        SET_VECTOR_ELT( Rret,1, category);
        SET_VECTOR_ELT( Rret,2, pval);
        SET_VECTOR_ELT( Rret,3, fdr);
        SET_VECTOR_ELT( Rret,4, enrichment_score);
        SET_VECTOR_ELT( Rret,5, percent_gene_hits_per_pathway);
        SET_VECTOR_ELT( Rret,6, number_hits);
        SET_VECTOR_ELT( Rret,7, number_misses);
        SET_VECTOR_ELT( Rret,8, number_user_genes);
        SET_VECTOR_ELT( Rret,9, total_genes_minus_input);
        SET_VECTOR_ELT( Rret,10, pathway_id);
        SET_VECTOR_ELT( Rret,11, pathway_type);
        SET_VECTOR_ELT( Rret,12, genesinpathway);
    
        PROTECT(cls = allocVector(STRSXP, 1)); // class attribute
        protect_cnt++;
    
        SET_STRING_ELT(cls, 0, mkChar("data.frame"));
        classgets(Rret, cls);
    
        PROTECT(nam = allocVector(STRSXP, maxflds));     // names attribute (column names)
        protect_cnt++;
    
        SET_STRING_ELT( nam, 0,  mkChar("pathway_name"));
        SET_STRING_ELT( nam, 1,  mkChar("category"));
        SET_STRING_ELT( nam, 2,  mkChar("pval"));
        SET_STRING_ELT( nam, 3,  mkChar("fdr"));
        SET_STRING_ELT( nam, 4,  mkChar("enrichment_score"));
        SET_STRING_ELT( nam, 5,  mkChar("percent_gene_hits_per_pathway"));
        SET_STRING_ELT( nam, 6,  mkChar("number_hits"));
        SET_STRING_ELT( nam, 7,  mkChar("number_misses"));
        SET_STRING_ELT( nam, 8,  mkChar("number_user_genes"));
        SET_STRING_ELT( nam, 9,  mkChar("total_genes_minus_input"));
        SET_STRING_ELT( nam, 10, mkChar("pathway_id"));
        SET_STRING_ELT( nam, 11, mkChar("pathway_type"));
        SET_STRING_ELT( nam, 12, mkChar("genesinpathway"));
    
        namesgets(Rret, nam);
    
        PROTECT(rownam = allocVector(STRSXP, num_used_paths )); // row.names attribute
        protect_cnt++;
        for (i=0 ; i<num_used_paths ; i++)
        {
            uptr = (u+i);
            SET_STRING_ELT(rownam, i, mkChar( (char *)uptr->acc) );
        }
        setAttrib(Rret, R_RowNamesSymbol, rownam);
    }
    if (user_in_genes) { free(user_in_genes); user_in_genes = (unsigned int *)0; }
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
    if (by_egids) { free(by_egids); by_egids = (void *)0; }
    if (real_universe) { free(real_universe); real_universe = (void *)0; }
    if (user_in_genes_original) free (user_in_genes_original);
    if (mycustompw) 
    {
        for (i=0;i<len_of_user_pws;i++) 
        { 
            mycustompwptr = mycustompw + i;
            if (mycustompwptr->name) free (mycustompwptr->name); 
            if (mycustompwptr->optional) free (mycustompwptr->optional); 
            if (mycustompwptr->genes) free (mycustompwptr->genes); 
            mycustompwptr->genes = (void *)0;
        }
        free(mycustompw); 
    }
    UNPROTECT(protect_cnt);
// fprintf(stderr,"gpccdbg in l2pgetuniverseR after frees \n");  fflush(stderr);
    return Rret;
}


SEXP l2pgetuniverseR(SEXP categories)
{
    char tmps_cat[PATH_MAX]; // temp string for "category" 
    char tmps[256];          // max gene name length is 22 = DTX2P1-UPK3BP1-PMS2P11
    char junks[32];
    struct used_path_type *uptr; // used path pointer 
    struct used_path_type *u = (struct used_path_type *)0;
    unsigned int num_used_paths = 0;
    unsigned int real_universe_cnt = 0;
    unsigned int *real_universe = (unsigned int *)0;
    unsigned int catspat = 0;
    unsigned int i;
    char *zz;
    int protect_cnt = 0;
    SEXP Rret;

// fprintf(stderr,"in l2pgetuniverseR 1, categories=%p\n",categories);  fflush(stderr);
    (void)setup_by_egids();
    tmps_cat[0] = junks[0] = tmps[1] = (char)0;
    if (Rf_isNull(categories))
    {
        category_set_all(&catspat);
    }
    else
    {
        strncpy( tmps_cat,CHAR(STRING_ELT(categories, 0)),PATH_MAX-2);
        (void)parsecats(tmps_cat,&catspat);     // set catpats
// fprintf(stderr,"in l2pgetuniverseR 2 cats=%x\n",catspat);  fflush(stderr);
    }

    u = setup_used_paths(&num_used_paths,catspat, junks,0,(unsigned int *)0,junks, &real_universe_cnt,&real_universe,0,(struct custom_type *)0);

// fprintf(stderr,"in l2pgetuniverseR 3\n");  fflush(stderr);
    PROTECT(Rret = allocVector(STRSXP, real_universe_cnt));
    protect_cnt++;
    for (i=0 ; i<real_universe_cnt ; i++)
    {
        zz = egid2hugo(*(real_universe+i));
        if (zz) sprintf(tmps,"%s",zz); 
        else strcpy(tmps,"NA"); 
        SET_STRING_ELT( Rret, i, mkChar(tmps) );
    }
// fprintf(stderr,"in l2pgetuniverseR 4\n");  fflush(stderr);
    if (u)  //free
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
// fprintf(stderr,"in l2pgetuniverseR 5\n");  fflush(stderr);
    if (real_universe) free(real_universe);
    if (by_egids) { free(by_egids); by_egids = (void *)0; }
    UNPROTECT(protect_cnt);
    return (SEXP)Rret;
}

SEXP l2pgetlongdesc(SEXP accarg, SEXP fpath)
{
    char s[10000];
    char path[PATH_MAX];
    char acc[PATH_MAX];
    char datadir[PATH_MAX];
    char fn[PATH_MAX*2]; 
    int lastslash;
    int i,j;
    char *z =(char *)0;
    long int offset;
    FILE *fp;
    SEXP Rret;
 
    (void)setup_by_egids();
    strncpy(acc,CHAR(STRING_ELT(accarg, 0)),PATH_MAX-2);
    strncpy(path,CHAR(STRING_ELT(fpath, 0)),PATH_MAX-2);
    lastslash = -1;
    for (i=0;path[i];i++)
    {
        if (path[i] == '/') lastslash = i;
    }
    if (lastslash > 0) path[lastslash] = (char)0;
    strcpy(datadir,path);
    strcat(datadir,"/");
    
    for (i=0 ; i<numpws ; i++)
    {
        if  (strcmp(acc,pws[i].acc) == 0)
        {
            sprintf(fn,"%s%s",datadir,"longdata.txt");
            fp = fopen(fn,"r");
            if (!fp) { return (SEXP)-1; }
            offset = (long int)pws[i].longdesc;
            fseek(fp,offset,SEEK_SET);
            s[0] = (char)0; 
            if (!fgets(s,(int)(sizeof(s)-1),fp)) {s[0] = (char)0; }
            fclose(fp);
            z = &s[0];
            for (j=0;s[j];j++) 
            {
                if ((s[j] == '\r') || (s[j] == '\n'))  // get rid of carriage return
                {
                  s[j] = (char)0; 
                  break;
                }
                else if ((s[j] < 32)  || (s[j] >= 127))
                  s[j] = (char)' ';  // hook me dudes up with no junk in the string
            }
            break;
         }
    }
    if (by_egids) { free(by_egids); by_egids = (void *)0; }
    Rret = PROTECT(allocVector(STRSXP, 1));
    if (z)
        SET_STRING_ELT(Rret, 0, mkChar(z));
    else
        SET_STRING_ELT(Rret, 0, mkChar(""));
    UNPROTECT(1);
    return Rret;
}


SEXP l2pgetgenes4acc(SEXP accarg)
{
    char acc[PATH_MAX];
    int i;
    unsigned int j;
    unsigned short int usi;
    char *z = (char *)0;
    SEXP Rret;

    acc[0] = (char)0;
    (void)setup_by_egids();
    strncpy(acc,CHAR(STRING_ELT(accarg, 0)),PATH_MAX-2);
    for (i=0 ; i<numpws ; i++)
    {
// fprintf(stderr,"in l2pgetgenes4acc() looping for acc=%s i=%d to %d\n",acc,i,numpws); fflush(stderr);
        if (strcmp(acc,pws[i].acc) == 0)
        {
            PROTECT(Rret = allocVector(STRSXP, pws[i].numgenes));
            for (j = 0 ; j<pws[i].numgenes ; j++)
            {
                usi = pwgenes[pws[i].pwgenesindex + j]; // okay 
                if (usi == USHRT_MAX) // already masked out as not in universe. (but no masking for this function?)
                {
fprintf(stderr,"ERROR: masked gene at i=%d j=%d \n",i,j); fflush(stderr);
                    continue;
                }
                z = egid2hugo(genes[usi].egid);
                if (z)
                    SET_STRING_ELT( Rret, j, mkChar(z) );
                else
                {
fprintf(stderr,"ERROR: cant find gene egid=%u gene at i=%d j=%d \n",usi,i,j); fflush(stderr);
                    SET_STRING_ELT( Rret, j, mkChar("") );
                }
            }
            UNPROTECT(1);
            if (by_egids) { free(by_egids); by_egids = (void *)0; }
            return (SEXP)Rret;
        }
    }
    if (by_egids) { free(by_egids); by_egids = (void *)0; }
// fprintf(stderr,"in l2pgetgenes4acc() after loop. returning null \n"); fflush(stderr);
    return (SEXP)R_NilValue;

}


