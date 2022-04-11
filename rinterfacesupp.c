
/*
l2psupp R interface
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

SEXP updategenes(SEXP ingenelist, SEXP trust_flag_arg)
{
    int protect_cnt = 0;
    SEXP cls = (SEXP)0;       // R class 
    SEXP Rret = (SEXP)0;
    SEXP nam = (SEXP)0;
    SEXP rownam = (SEXP)0;
    SEXP newname = (SEXP)0;
    SEXP oldname = (SEXP)0;
    SEXP change_flag = (SEXP)0;
    SEXP status = (SEXP)0;
    SEXP is_legit_name = (SEXP)0;
    unsigned int len;
    unsigned int ui;
    char **genes;
    struct updated_genes_type *updated_genes;
    struct updated_genes_type *ug;
    unsigned int maxflds = 5;
    int trust_flag;
    char tmps[5012];


    if (Rf_isNull(ingenelist))
    {
        return (SEXP)R_NilValue;
    }
    if (Rf_isNull(trust_flag_arg))   { trust_flag = 1; }   else { trust_flag = asInteger(trust_flag_arg); }
//fprintf(stderr,"trust_flag = %d (%p)\n",trust_flag,trust_flag_arg); 

    len = length(ingenelist);
    genes = malloc(sizeof(char *)*len);
    if (!genes) return (SEXP)R_NilValue;
    for (ui=0;ui<len;ui++)
    {
        strcpy(tmps,CHAR(STRING_ELT(ingenelist, ui)));
        *(genes+ui) = strdup(tmps); // be sure to free these
    }
    updated_genes = updategenesR(genes,len);
    for (ui=0;ui<len;ui++)
    {
        if (*(genes+ui)) free(*(genes+ui));
    }
    free(genes);
    if (trust_flag)
    {
       Rret = PROTECT(allocVector(STRSXP, len));
       for (ui=0;ui<len;ui++)
       {
           SET_STRING_ELT(Rret, ui, mkChar(updated_genes[ui].newname));
       }
       for (ui=0;ui<len;ui++)
       {
           if (updated_genes[ui].newname) free(updated_genes[ui].newname); 
           if (updated_genes[ui].oldname) free(updated_genes[ui].oldname); 
       }
       free(updated_genes);
       UNPROTECT(1);
       return Rret;
    }

#if 0
struct updated_genes_type { char *newname; char *oldname; int change_flag; int status; int is_legit_name; };
#endif
   maxflds = 5; 
   PROTECT(Rret = Rf_allocVector(VECSXP, maxflds)); // a list
   protect_cnt++;
   PROTECT(newname=Rf_allocVector(STRSXP, len));              // 1 
   protect_cnt++;
   PROTECT(oldname=Rf_allocVector(STRSXP, len));              // 2 
   protect_cnt++;
   PROTECT(change_flag=Rf_allocVector(INTSXP, len));       // 3 
   protect_cnt++;
   PROTECT(status=Rf_allocVector(INTSXP, len));               // 4 
   protect_cnt++;
   PROTECT(is_legit_name=Rf_allocVector(INTSXP, len));       // 5 
   protect_cnt++;
   for (ui=0 ; ui<len ; ui++)
   {
       ug = (updated_genes+ui);
       SET_STRING_ELT(newname, ui, mkChar(ug->newname) ); // 1
       SET_STRING_ELT(oldname, ui, mkChar(ug->oldname) ); // 2
       INTEGER(change_flag)[ui] = ug->change_flag;   // 3 
       INTEGER(status)[ui] = ug->status;               // 4 
       INTEGER(is_legit_name)[ui] = ug->is_legit_name;    // 5 
    }

// xxx
    SET_VECTOR_ELT( Rret,0, newname);
    SET_VECTOR_ELT( Rret,1, oldname);
    SET_VECTOR_ELT( Rret,2, change_flag);
    SET_VECTOR_ELT( Rret,3, status);
    SET_VECTOR_ELT( Rret,4, is_legit_name);
    
    PROTECT(cls = allocVector(STRSXP, 1)); // class attribute
    protect_cnt++;
    
    SET_STRING_ELT(cls, 0, mkChar("data.frame"));
    classgets(Rret, cls);
    
    PROTECT(nam = allocVector(STRSXP, maxflds));     // names attribute (column names)
    protect_cnt++;

    SET_STRING_ELT( nam, 0,  mkChar("newname"));
    SET_STRING_ELT( nam, 1,  mkChar("oldname"));
    SET_STRING_ELT( nam, 2,  mkChar("change_flag"));
    SET_STRING_ELT( nam, 3,  mkChar("status"));
    SET_STRING_ELT( nam, 4,  mkChar("is_legit_name"));

    namesgets(Rret, nam);
    
    PROTECT(rownam = allocVector(STRSXP, len )); // row.names attribute
    protect_cnt++;
    for (ui=0 ; ui<len ; ui++)
    {
        sprintf(tmps,"%u",ui);
        SET_STRING_ELT(rownam, ui, mkChar( (char *)&tmps[0]) );
    }
    setAttrib(Rret, R_RowNamesSymbol, rownam);
// fprintf(stderr,"gpccdbg in l2pgetuniverseR after frees \n");  fflush(stderr);
    for (ui=0;ui<len;ui++)
    {
        if (updated_genes[ui].newname) free(updated_genes[ui].newname); 
        if (updated_genes[ui].oldname) free(updated_genes[ui].oldname); 
    }
    free(updated_genes);
    UNPROTECT(protect_cnt);
    return (SEXP)Rret;
}


SEXP egids2hugosR(SEXP ingenelist) //, SEXP trust_flag_arg)
{
    SEXP Rret = (SEXP)0;
    int j;
    unsigned int len;
    unsigned int newlen;
    unsigned int ui;
    unsigned int *input_egids;
    struct entrez_hugo_ensemble_type *es2hs; // "entrez's to hugos" 
#if 0
    SEXP cls = (SEXP)0;       // R class 
    SEXP nam = (SEXP)0;
    SEXP rownam = (SEXP)0;
    SEXP newname = (SEXP)0;
    int protect_cnt = 0;
    int trust_flag;
    unsigned int maxflds = 5;
    char tmps[5012];
#endif

// fprintf(stderr,"in egids2hugosR start\n"); fflush(stderr); 
    if (Rf_isNull(ingenelist))
    {
        return (SEXP)R_NilValue;
    }
    len = length(ingenelist);
// fprintf(stderr,"in egids2hugosR lenght ingenelist is  = %d ( at %p) \n",len,ingenelist); fflush(stderr); 
    input_egids = malloc(sizeof(unsigned int)*len);
    if (!input_egids) return (SEXP)R_NilValue;
// fprintf(stderr,"in egids2hugosR getting input len = %d\n",len); fflush(stderr); 
    newlen = 0;
    for (ui=0; ui<len ;ui++)
    {
// fprintf(stderr,"in egids2hugosR in loop %u to %u\n",ui,len); fflush(stderr); 
#if 1
        j =  INTEGER(ingenelist)[ui];
#else
    double d;
        d = REAL(ingenelist)[ui];
        j = (int)d;
#endif
// fprintf(stderr,"in egids2hugosR in loop %u to %u j=%d from %f\n",ui,len,j,d); fflush(stderr); 
        if (j > 0)
            input_egids[newlen++] = j;
    }
// fprintf(stderr,"in egids2hugosR before egids2hugos(), len now %d\n",newlen); fflush(stderr); 
    es2hs = egids2hugos(input_egids, newlen);
// fprintf(stderr,"in egids2hugosR after egids2hugos() ret=%p\n",es2hs); fflush(stderr); 
    free(input_egids);
    if (!es2hs)
    {
        return (SEXP)R_NilValue;
    }
    newlen = 0;
    for (ui=0;    ;ui++)
    {
        if (es2hs[ui].gene_id == 0) break;
        newlen++;
    }
// fprintf(stderr,"in egids2hugosR near end newlen now =%d\n",newlen); fflush(stderr); 
    Rret = PROTECT(allocVector(STRSXP, newlen));
    for (ui=0;ui<newlen;ui++)
    {
        SET_STRING_ELT(Rret, ui, mkChar(es2hs[ui].hugo));
    }
    free(es2hs); 
    UNPROTECT(1);
// fprintf(stderr,"in egids2hugosR at end , will return %p\n",Rret); fflush(stderr); 
    return (SEXP)Rret;
}
    

