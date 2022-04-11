
#include <string.h>
#include "pathworks.h"

void category_set_all(unsigned int *pat)
{
    *pat = (
   CAT_NCBI_BIOCYC| CAT_NCBI_GO    | CAT_NCBI_KEGG  | CAT_NCBI_PANTH | CAT_NCBI_PID | CAT_NCBI_REACTOME     | CAT_NCBI_WikiPathways |
   CAT_MSIG_C1   | CAT_MSIG_C2   | CAT_MSIG_C3   | CAT_MSIG_C4   | CAT_MSIG_C5   | CAT_MSIG_C6   | CAT_MSIG_C7 | CAT_MSIG_C8  | CAT_MSIG_H  ) ;
}


int string_to_category_code(char cats[])
{
     if  (strcmp(cats,"BIOCYC") == 0)  return CAT_NCBI_BIOCYC;
else if  (strcmp(cats,"GO") == 0)      return CAT_NCBI_GO;
else if  (strcmp(cats,"KEGG") == 0)    return CAT_NCBI_KEGG;
else if  (strcmp(cats,"PANTH") == 0)   return CAT_NCBI_PANTH;

else if  (strcmp(cats,"PID") == 0)     return CAT_NCBI_PID;
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
else if  (strcmp(cats,"C8") == 0) return CAT_MSIG_C8;
else if  (strcmp(cats,"H") == 0) return CAT_MSIG_H;
else if  (strcmp(cats,"CUSTOM") == 0) return CAT_CUSTOM;
else return 0;

// else if  (strcmp(cats,"ARCHIVED") == 0)   return CAT_MSIG_ARCHIVED; 
}

void categories_pattern_to_strings(unsigned int cat,char puthere[])
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
    if (cat & CAT_MSIG_C8   ) strcat(puthere,"C8 ");
    if (cat & CAT_MSIG_H    ) strcat(puthere,"H ");
    if (cat & CAT_CUSTOM    ) strcat(puthere,"CUSTOM ");
}

unsigned int string2type(char *s)
{
        if (strcmp("functional_set",s) == 0)     return type_functional_set;
   else if (strcmp("pathway",s) == 0)            return type_pathway;
   else if (strcmp("structural_complex",s) == 0) return type_structural_complex;
/* 4,5 are scope */
   else if (strcmp("custom",s) == 0)             return type_custom;
   else if (strcmp("unknown",s) == 0)            return type_unknown;
   return 9;
}

