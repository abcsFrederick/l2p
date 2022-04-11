
/*
vi updategenes.c ; gcc -minline-all-stringops -Wall -o updategenes updategenessupport.c updategenes.c
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <limits.h>
#include <math.h>
#include <ctype.h>
#include <errno.h>

#include "pathworks.h"


#ifdef L2P_USING_R
#include <R.h>
#include <Rdefines.h>
#endif

extern struct synonym_type synonyms[];
extern char *hugo[];
extern int numhugo;
extern int numehe;
extern struct entrez_hugo_ensemble_type ehe[];

int cmp_sort_by_syn_name(const void *a, const void *b)
{
   struct synonym_type *aa;
   struct synonym_type *bb;
   aa = (struct synonym_type *)a;
   bb = (struct synonym_type *)b;
   return(strcmp(aa->Synonym, bb->Synonym));
}

int myStrCmp(const void *s1, const void *s2) {
  const char *key = s1;
  const char * const *arg = s2;
  // printf("myStrCmp: s1(%p): %s, s2(%p): %s\n", s1, key, s2, *arg);
  return strcmp(key, *arg);
}


#ifdef L2P_USING_R
struct updated_genes_type *updategenesR(char *genes[], const int len)
{
    char nm[5120]; // gene name
    char oldname[5120]; // gene name
    char newname[5120]; // gene name
    char hugo_string[512];
    int i,status;
    int change_flag = 0;
    int num_syns = 0;
    int is_legit_name = 0;
    struct synonym_type *syptr;
    struct synonym_type syspace;
    char *z;
    struct updated_genes_type *ret = (struct updated_genes_type *)0;

    for (num_syns=0; synonyms[num_syns].Synonym; num_syns++) { /* get number of synonyms */ }

// struct updated_genes_type { char *newname; char *oldname; int change_flag; int status; int is_legit_name; };

    ret = (struct updated_genes_type *)malloc(sizeof(struct updated_genes_type) * len);
    if (!ret) return ret;

    for (i=0;i<len;i++)
    {
        strcpy(nm,genes[i]);
        strcpy(oldname,nm);
        strcpy(newname,nm);
        is_legit_name = 0;
        strcpy(hugo_string,nm);
        is_legit_name = 0;
#if 0
        for (i=0;i<numhugo;i++)
        { 
// fprintf(stderr,"%s %s\n",hugo_string,hugo[i]); 
        if (strcmp(hugo_string,hugo[i]) == 0) 
        { 
             is_legit_name = 1; 
             break; 
        } 
        }
#else
// fprintf(stderr,"search for %s\n",hugo_string); 
        z = (char *)bsearch(&hugo_string[0],&hugo[0], numhugo, sizeof(char *), myStrCmp);
        if (z) is_legit_name = 1;
        // else fprintf(stderr,"not found for %s\n",hugo_string); 
#endif
        memset(&syspace,0,sizeof(syspace));
        syspace.Synonym = &nm[0];
        syptr = bsearch(&syspace,synonyms,num_syns,sizeof(struct synonym_type),cmp_sort_by_syn_name);
#if 0
        if (syptr) 
        {
            strcpy(newname,syptr->Symbol);
            change_flag = 1;
        }
// modes: 1=all (new and old)  2=just changes(just differences) 3=trustus(updated genes)
        if (mode == 1)
        {
            printf("%s %s\n",newname,oldname);
        }
        else if (mode == 2)
        {
            if (change_flag == 1)
                printf("%s %s\n",newname,oldname);
        }
        else if (mode == 3)
        {
            printf("%s\n",newname);
        }
#endif
        status = change_flag = 0;
        if (syptr)
        {
// printf("dbg: %s %s %d\n",syptr->Symbol,syptr->Synonym,syptr->status);
            if (syptr->status == 1)
                 status = 1;
            else
            {
                strcpy(newname,syptr->Symbol);
                change_flag = 1;
            }
        }
        (ret+i)->newname = strdup(newname); 
        (ret+i)->oldname = strdup(oldname); 
        (ret+i)->change_flag = change_flag;
        (ret+i)->status = status;
        (ret+i)->is_legit_name = is_legit_name;
        // printf("%s\t%s\t%d\t%d\t%d\n",newname,oldname,change_flag,status,is_legit_name);
    }
    return ret;
}


struct entrez_hugo_ensemble_type *egids2hugos(unsigned int egids[], const int len) // free if return is not null
{
    int i,j,k;
    struct entrez_hugo_ensemble_type *ret = (struct entrez_hugo_ensemble_type *)0;

// xxx 
    ret = (struct entrez_hugo_ensemble_type *)malloc(sizeof(struct entrez_hugo_ensemble_type) * (len + 1)); // will "null" terminate this
    if (!ret) return ret;

    for (i=k=0; i<len;i++)
    {
// fprintf(stderr,"egids2hugos() input %d of %d is %d\n",i,len,egids[i]);
        for (j=0;j<numehe;j++)
        {
            if (ehe[j].gene_id == egids[i]) 
            {
                ret[k++] = ehe[j]; // note: don't need to free hugo or ensembl field
                break;
            }
        }
    }
    ret[k].gene_id = 0; // "null terminate"
/*
    for (i=0;ret[i].gene_id;i++)
    {
        fprintf(stderr,"egids2hugos() got %s for %d\n",ret[i].hugo,ret[i].gene_id);
    }
*/
    return ret;
}
#else

int geneupdate_command_line(const int mode)
{
    char nm[5120]; // gene name
    char oldname[5120]; // gene name
    char newname[5120]; // gene name
    char hugo_string[512];
    int i,status;
    int change_flag = 0;
    int num_syns = 0;
    int is_legit_name = 0;
    struct synonym_type *syptr;
    struct synonym_type syspace;

    for (num_syns=0; synonyms[num_syns].Synonym; num_syns++) { /* get number of synonyms */ }

    while ( fgets(nm, sizeof(nm), stdin) )
    {
        for (i=0;nm[i];i++) { if ((nm[i] == '\r') || (nm[i] == '\n')) nm[i] = (char)0; }
        strcpy(oldname,nm);
        strcpy(newname,nm);
        is_legit_name = 0;
        strcpy(hugo_string,nm);


        is_legit_name = 0;
#if 0
        for (i=0;i<numhugo;i++)
        { 
// fprintf(stderr,"%s %s\n",hugo_string,hugo[i]); 
        if (strcmp(hugo_string,hugo[i]) == 0) 
        { 
             is_legit_name = 1; 
             break; 
        } 
        }
#else
    char *z;
// fprintf(stderr,"search for %s\n",hugo_string); 
        z = (char *)bsearch(&hugo_string[0],&hugo[0], numhugo, sizeof(char *), myStrCmp);
        if (z) is_legit_name = 1;
        // else fprintf(stderr,"not found for %s\n",hugo_string); 
#endif

        memset(&syspace,0,sizeof(syspace));
        syspace.Synonym = &nm[0];
        syptr = bsearch(&syspace,synonyms,num_syns,sizeof(struct synonym_type),cmp_sort_by_syn_name);
#if 0
        if (syptr) 
        {
            strcpy(newname,syptr->Symbol);
            change_flag = 1;
        }
// modes: 1=all (new and old)  2=just changes(just differences) 3=trustus(updated genes)
        if (mode == 1)
        {
            printf("%s %s\n",newname,oldname);
        }
        else if (mode == 2)
        {
            if (change_flag == 1)
                printf("%s %s\n",newname,oldname);
        }
        else if (mode == 3)
        {
            printf("%s\n",newname);
        }
#endif
        status = 0;
        change_flag = 0;
        if (syptr)
        {
// printf("dbg: %s %s %d\n",syptr->Symbol,syptr->Synonym,syptr->status);
            if (syptr->status == 1)
                 status = 1;
            else
            {
                strcpy(newname,syptr->Symbol);
                change_flag = 1;
            }
        }
// char *hugoptr; hugoptr = bsearch(&hugospace,synonyms,num_syns,sizeof(struct synonym_type),cmp_sort_by_syn_name);
        printf("%s\t%s\t%d\t%d\t%d\n",newname,oldname,change_flag,status,is_legit_name);
    }
    return 0;
}

int main(int argc,char *argv[])
{
    int bad;
    int mode = 1; // default =1 , 1=all 2=just changes 3=trust 

    mode = 1;
    bad = 0;
    if (argc > 2) bad = 1;
    if (argc == 2)
    {
        if ((strcmp(argv[1],"help")==0) ||
            (strcmp(argv[1],"-help")==0) ||
            (strcmp(argv[1],"--help")==0) ||
            (strcmp(argv[1],"-h")==0) ||
            (strcmp(argv[1],"--h")==0)
           )
         {
             bad = 1;
         }
         else
         {
             mode = atoi(argv[1]);
             if ((mode < 1) || (mode > 3))
             {
                 bad = 1;
             }
         }
    }
    if (bad)
    { 
        fprintf(stderr,"usage: updategenes [mode]\n");
        fprintf(stderr,"modes: 1=all (new and old)  2=just changes(just differences) 3=trustus(updated genes)\n");
        return 0;
    }
    geneupdate_command_line(mode);

    return 0;
}

#endif
