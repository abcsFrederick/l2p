
/*
     1  pval
     2  fdr
     3  ratio                      if postive, genes are OVER REPRESENTED, if negative genes are UNDER REPRESENTED
     4  pwhitcount                 number of genes hit in pathway
     5  pwnohitcount               pathway    number of genes in the pathway
     6  inputcount                  total count of user genes (user input)
     7  pwuniverseminuslist        total number of unique genes in all pathways
     8  pathwayaccessionidentifier canonical accession ( if available, otherwise assigned by us )
     9  category                   KEGG,REACTOME,GO,PANT(=PANTHER),PID=(pathway interaction database)  *was "soruce"*
    10  pathwayname                Name of pathway
    11  pathwaytype genes_space_separated   HUGO genes from user that hit the pathway

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
#ifdef WEBASSEMBLY
#include "small.h"
#else
#include "big.h"
#endif

#define RDEBUG 0

extern unsigned short int pwgenes[];
extern struct smallgenetype genes[];
extern struct pwtype pws[];
extern int numpws;
extern int numgenes;
extern int numpwgenes;
extern struct smallgenetype *by_egids; // a copy of "genes" ,but sorted by egid (entrez gene id )

void category_code_to_string(unsigned int cat,char puthere[])
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
    else if (cat & CAT_CUSTOM    ) strcpy(puthere,"CUSTOM");
    else strcpy(puthere,"UNKNOWN");
//    else if (cat & CAT_MSIG_ARCHIVED     )   strcpy(puthere,"ARCHIVED"); 
}

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

int bitCount_unsigned(unsigned int n) 
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




// int l2pfunc_R(unsigned int *user_in_genes, unsigned int user_incnt, struct used_path_type *usedpaths,unsigned int num_used_paths,unsigned int real_universe_cnt,unsigned int *real_universe,  int permute_flag)

int l2pfunc_R(unsigned int *user_in_genes, unsigned int user_incnt, struct used_path_type *usedpaths,unsigned int num_used_paths,unsigned int real_universe_cnt, int permute_flag)
{
    struct used_path_type *uptr; // used path pointer 
    unsigned int ui_uk;
    unsigned int ui_ej;
    unsigned int i;
    unsigned int j,k,ll;
    int ret = 0;

    
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
    if (permute_flag)
    {
#if 1
        nelson_method(user_incnt, usedpaths,num_used_paths,real_universe_cnt);
#else
        fast_permutes(user_incnt, usedpaths,num_used_paths,real_universe_cnt,real_universe);
#endif
    }
    return ret;
}

unsigned int hugo2egid(char *h);
int cmp_ui(const void *a, const void *b);


SEXP l2p(SEXP lst, SEXP categories, SEXP universe, SEXP custompws, SEXP customfn, SEXP universefn, SEXP permuteval, SEXP oneside_arg)
{
    char tmps_cat[PATH_MAX];  // temp string for "category" 
    char tmps[512];
    char universe_file[PATH_MAX];
    char custom_file[PATH_MAX];
    unsigned int i,k;
    unsigned int j;
    unsigned int len_of_list = 0;  // length of list of lists , from Rf_length() which return R_len_t ( which in int? )
    unsigned int len_of_vector = 0;
    struct used_path_type *uptr;
    unsigned int *in_universe = (unsigned int *)0;
    unsigned int *in_universe_original = (unsigned int *)0;
    unsigned int in_universe_cnt = 0;
    unsigned int *real_universe = (unsigned int *)0;
    unsigned int real_universe_cnt = 0;
    int permute_flag = 0;
    int user_universe_flag = 0;
    int user_incnt = 0;
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
    unsigned int maxflds = 11;
    struct custom_type *mycustompw = (struct custom_type *)0;
    struct custom_type *mycustompwptr = (struct custom_type *)0;
    double fdr_for_output;
    SEXP list = (SEXP)0;
    SEXP pval = (SEXP)0;     // 1 
    SEXP fdr = (SEXP)0;      // 2
    SEXP ratio = (SEXP)0;         // 3 if positive, genes are OVER REPRESENTED, if negative genes are UNDER REPRESENTED
    SEXP pwhitcount = (SEXP)0;       // 4 number of genes hit in pathway
    SEXP pwnohitcount = (SEXP)0;     // 5 number of genes in the pathway
    SEXP inputcount = (SEXP)0;       // 6 total count of user genes (user input)
    SEXP pwuniverseminuslist = (SEXP)0;    // 7 total number of unique genes in all pathways
    SEXP pathwayaccessionidentifier = (SEXP)0; // 8 canonical accession ( if available, otherwise assigned by us )
    SEXP category = (SEXP)0;                   // 9 KEGG,REACTOME,GO,PANT(=PANTHER),PID=(pathway interaciton database)
    SEXP pathwayname = (SEXP)0;                // 10 Name of pathway
    SEXP genesinpathway = (SEXP)0;             // 11 genes_space_separated   HUGO genes from user that hit the pathway
    SEXP Rret = (SEXP)0;
    SEXP cls = (SEXP)0;       // R class 
    SEXP nam = (SEXP)0;       // R name 
    SEXP rownam = (SEXP)0;    // row names
//    SEXP tmpsexp = (SEXP)0;
//     int status = 0; not used
//    int lastslash = 0; not used
//    int custom_flag = 0;
    //int precise_flag = 0;


    (void)setup_by_egids();
    if (Rf_isNull(permuteval)) 
    {
        permute_flag = 0;
    }
    else 
    { 
        permute_flag = asInteger(permuteval);
    }

    if (Rf_isNull(oneside_arg)) 
    {
        oneside = 0;
    }
    else 
    { 
        oneside = asInteger(oneside_arg);
        if ((oneside < 0) || (oneside > 2)) oneside = 0;      // just give them twosided fisher's exact test
    }

    universe_file[0] = custom_file[0] = tmps[0] = tmps_cat[0] = (char)0;
    if (Rf_isNull(custompws))      //  if (custompws == (SEXP)0)
    {
//        custom_flag = 0; not used
    }
    else
    {
/* struct custom_type { char *name; char *optional; unsigned int numgenes; unsigned int *genes; }; */
//        custom_flag = 1; not used
        len_of_list = (int)length(custompws);
        mycustompw = (struct custom_type *)malloc(sizeof(struct custom_type )*len_of_list); // remember to free this
        memset(mycustompw,0,sizeof(struct custom_type )*len_of_list);
        for (i=0;i<len_of_list;i++)
        {
            mycustompwptr = mycustompw + i;
            list = VECTOR_ELT(custompws, i);
            len_of_vector = length(list);
            mycustompwptr->genes = (unsigned int *)malloc(sizeof(unsigned int)*len_of_vector);
            mycustompwptr->numgenes = 0;   // will set this properly later
            for (j=0;j<len_of_vector;j++)
            {
                tmps_cat[0] = (char)0;
                strncpy(tmps_cat,CHAR(STRING_ELT(list, j)),PATH_MAX-2);
                if (j == 0) mycustompwptr->name = strdup(tmps_cat);    // name
                else if (j == 1) mycustompwptr->optional = strdup(tmps_cat); // optional,  might be a URL
                else 
                {
                    ui = hugo2egid(tmps_cat);
// fprintf(stderr,"xxx custgene is i=%d j=%d %u %s\n",i,j,ui,tmps_cat);  fflush(stderr); 
                    if (ui != UINT_MAX) { *(mycustompwptr->genes + j - 2) = ui; mycustompwptr->numgenes++; } // nb:-2
                }
            }
            qsort(mycustompwptr->genes,mycustompwptr->numgenes,sizeof(unsigned int),cmp_ui);
        }
    }
    if (Rf_isNull(universe)) //   if (universe == (SEXP)0)
        user_universe_flag = 0;
    else
        user_universe_flag = 1;
    if (Rf_isNull(categories ))
    {
        category_set_all(&catspat);
    }
    else
    {
        strncpy( tmps_cat,CHAR(STRING_ELT(categories, 0)),PATH_MAX-2);
        (void)parsecats(tmps_cat,&catspat);     // set catpats
    }
    if (Rf_isNull(customfn)) {} else strncpy(custom_file,CHAR(STRING_ELT(customfn, 0)),PATH_MAX-2);
    if (Rf_isNull(universefn)) {} else strncpy(universe_file,CHAR(STRING_ELT(universefn, 0)),PATH_MAX-2);

    len = length(lst);
    user_in_genes = (unsigned int *)malloc(sizeof(unsigned int)*len); // remember to free this
    user_in_genes_original = (unsigned int *)malloc(sizeof(unsigned int)*len); // remember to free this

    if (!user_in_genes)
        return (SEXP) -1; // why not 0 ?
    for (k = i = 0; i < len; i++)
    {
        strcpy(tmps,CHAR(STRING_ELT(lst, i)));
        ui = hugo2egid(tmps);
        if (ui != UINT_MAX) { *(user_in_genes_original+k) = ui; k++; }
    }
    qsort(user_in_genes_original,k,sizeof(unsigned int),cmp_ui);

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
    user_incnt = j;

    len = length(universe);
    if (user_universe_flag == 1)
    {
        in_universe = (unsigned int *)malloc(sizeof(unsigned int)*len);             // remember to free this
        if (!in_universe) return (SEXP) -1;
        in_universe_original = (unsigned int *)malloc(sizeof(unsigned int)*len);    // remember to free this
        if (!in_universe_original) 
	{
            free(in_universe);
            return (SEXP) -1;
	}
        for (k = i = 0; i < len; i++)
        {
            strcpy(tmps,CHAR(STRING_ELT(universe, i)));
            ui = hugo2egid(tmps);
            if (ui != UINT_MAX) { *(in_universe_original+k) = ui; k++; }
        }
        qsort(in_universe_original,k,sizeof(unsigned int),cmp_ui);
        for (j=i=0;i<k;i++)
        {  // de duplicate universe 
            if (i > 0) 
            {
                if ( *(in_universe_original+i) == *(in_universe_original+i-1) )
                   continue;
                *(in_universe+j) = *(in_universe_original+i);
                j++;
            }
            else 
            {
                *(in_universe+j) = *(in_universe_original+i);
                j++;
            }
        }
        in_universe_cnt = j;
    }
// xxx 
    if (in_universe_original) { free(in_universe_original); in_universe_original = (void *)0; }

    u = setup_used_paths(&num_used_paths, catspat,universe_file,in_universe_cnt,in_universe,custom_file,&real_universe_cnt,&real_universe,len_of_list,mycustompw);
// fprintf(stderr,"in l2pgetuniverseR 2.1 cats=%x xxx after setup_used_paths() \n",catspat);  fflush(stderr);
// NO, freed in setup_used_paths    if (in_universe) { free(in_universe); in_universe = (void *)0; }

    // (void)l2pfunc_R(user_in_genes,user_incnt,u,num_used_paths,real_universe_cnt,real_universe,permute_flag);
    (void)l2pfunc_R(user_in_genes,user_incnt,u,num_used_paths,real_universe_cnt,permute_flag);
    do_pvals_and_bh(user_incnt, u,num_used_paths,real_universe_cnt,oneside);

    PROTECT(Rret = Rf_allocVector(VECSXP, 11)); // a list with 11 elements
    protect_cnt++;
    for (i=0 ; i<maxflds ; i++) // maxflds = 11 for now
    {
        PROTECT(pval=Rf_allocVector(REALSXP, num_used_paths ));
        PROTECT(fdr=Rf_allocVector(REALSXP, num_used_paths));
        PROTECT(ratio=Rf_allocVector(REALSXP, num_used_paths));
        PROTECT(pwhitcount=Rf_allocVector(INTSXP, num_used_paths));
        PROTECT(pwnohitcount=Rf_allocVector(INTSXP, num_used_paths));
        PROTECT(inputcount=Rf_allocVector(INTSXP, num_used_paths));
        PROTECT(pwuniverseminuslist=Rf_allocVector(INTSXP, num_used_paths));
        PROTECT(pathwayaccessionidentifier=Rf_allocVector(STRSXP, num_used_paths));
        PROTECT(category=Rf_allocVector(STRSXP, num_used_paths));   // is natively an int, but convert to string
        PROTECT(pathwayname=Rf_allocVector(STRSXP, num_used_paths));
        PROTECT(genesinpathway=Rf_allocVector(STRSXP, num_used_paths));
        protect_cnt+=11;

    }
// fprintf(stderr,"aaa 1 ,real_universe_cnt = %d, user_incnt = %d \n",real_universe_cnt,user_incnt); fflush(stderr);
    for (ui=0 ; ui<num_used_paths ; ui++)
    {
        uptr = (u+ui);
        REAL(pval)[ui] = uptr->pval;                                      // 0 
        if (permute_flag) fdr_for_output = uptr->fdr2; // use new method
        else              fdr_for_output = uptr->fdr;
        REAL(fdr)[ui] =   fdr_for_output;                                 // 1 
        REAL(ratio)[ui] = uptr->ad;                                       // 2 
        INTEGER(pwhitcount)[ui] =  uptr->a;                               // 3 
        INTEGER(pwnohitcount)[ui] =  uptr->b;                             // 4 
        INTEGER(inputcount)[ui] = uptr->c;                                // 5 
        INTEGER(pwuniverseminuslist)[ui] = uptr->d;                       // 6 

        SET_STRING_ELT(pathwayaccessionidentifier, ui, mkChar(uptr->acc) );     // 7 

        category_code_to_string( uptr->category , tmps_cat);             // 8 
        SET_STRING_ELT(category, ui, mkChar(tmps_cat));                   // 9 

        SET_STRING_ELT(pathwayname, ui, mkChar(uptr->name) );             // 10 
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

   SET_VECTOR_ELT( Rret,0, pval);
   SET_VECTOR_ELT( Rret,1, fdr);
   SET_VECTOR_ELT( Rret,2, ratio);
   SET_VECTOR_ELT( Rret,3, pwhitcount);
   SET_VECTOR_ELT( Rret,4, pwnohitcount);
   SET_VECTOR_ELT( Rret,5, inputcount);
   SET_VECTOR_ELT( Rret,6, pwuniverseminuslist);
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
   SET_STRING_ELT( nam,                 3, mkChar("pwhitcount"));
   SET_STRING_ELT( nam,               4, mkChar("pwnohitcount"));
   SET_STRING_ELT( nam,                 5, mkChar("inputcount"));
   SET_STRING_ELT( nam,        6, mkChar("pwuniverseminuslist"));
   SET_STRING_ELT( nam, 7, mkChar("pathwayaccessionidentifier"));
   SET_STRING_ELT( nam,                   8, mkChar("category"));
   SET_STRING_ELT( nam,                9, mkChar("pathwayname"));
   SET_STRING_ELT( nam,            10, mkChar("genesinpathway"));
   namesgets(Rret, nam);

   PROTECT(rownam = allocVector(STRSXP, num_used_paths )); // row.names attribute
   protect_cnt++;
   for (i=0 ; i<num_used_paths ; i++)
   {
       uptr = (u+i);
       SET_STRING_ELT(rownam, i, mkChar( (char *)uptr->acc) );
   }
   setAttrib(Rret, R_RowNamesSymbol, rownam);

// aaa
// fprintf(stderr,"in l2pgetuniverseR before frees \n");  fflush(stderr);
    if (user_in_genes) free(user_in_genes);
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
// fprintf(stderr,"aaa returning %p \n",Rret); fflush(stderr); 
    if (user_in_genes_original) free (user_in_genes_original);
    if (mycustompw) 
    {
        for (i=0;i<len_of_list;i++) 
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
// fprintf(stderr,"in l2pgetuniverseR after frees \n");  fflush(stderr);
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
    char fn[PATH_MAX];
    int lastslash;
    int i,j;
    char *z =(char)0;
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
    Rret = PROTECT(allocVector(STRSXP, 1));
    if (z)
        SET_STRING_ELT(Rret, 0, mkChar(z));
    else
        SET_STRING_ELT(Rret, 0, mkChar(""));
    if (by_egids) { free(by_egids); by_egids = (void *)0; }
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


