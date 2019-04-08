
#ifndef   PATWORKS_H 
#define   PATWORKS_H 
#endif

#define MAXBSID 25000
#define MAXBSIDPOSSIBLE 851568
#define MAXGENE 70000

#define CAT_NCBI_BIOCYC    (1<<0)     
                                           // count = 294
#define CAT_NCBI_GO        (1<<1)     
                                           // count = 13515 
#define CAT_NCBI_KEGG      (1<<2)   
                                           // count = 485 
#define CAT_NCBI_PANTH     (1<<3)     
                                           // count = 129 
#define CAT_NCBI_PID       (1<<4)                 
                                           // count = 183 
#define CAT_NCBI_REACTOME      (1<<5)  
                                           // count = 1548 
#define CAT_NCBI_WikiPathways  (1<<6)  
                                           // count = 345 
#define CAT_MSIG_C1            (1<<7)  
                                           // count = 325 
#define CAT_MSIG_C2            (1<<8)  
                                           // count = 3777 
#define CAT_MSIG_C3            (1<<9)  
                                           // count = 836 
#define CAT_MSIG_C4            (1<<10)  
                                           // count = 858 
#define CAT_MSIG_C5            (1<<11)  
                                           // count = 5871 
#define CAT_MSIG_C6            (1<<12)  
                                           // count = 144 
#define CAT_MSIG_C7            (1<<13)  
                                           // count = 1888 
#define CAT_MSIG_H             (1<<14)  
                                           // count = 50 

#if 0
#define CAT_MSIG_ARCHIVED      (1<<15)  
                                           // count = 858 
#endif
 
#define NCBI_PAT = (CAT_NCBI_BIOCYC|CAT_NCBI_GO|CAT_NCBI_KEGG|CAT_NCBI_PANTH|CAT_NCBI_PID|CAT_NCBI_REACTOME|CAT_NCBI_WikiPathways);
#define MSIG_PAT = (CAT_MSIG_C1|CAT_MSIG_C2|CAT_MSIG_C3|CAT_MSIG_C4|CAT_MSIG_C5|CAT_MSIG_C6|CAT_MSIG_C7|CAT_MSIG_H);

struct binpathouttype     // the "binary" pathway information file  
{ 
    int bsid;         //32 bit int
    int category;     // one categoery (particular bit set) examples:  CAT_NCBI_GO        
    int accession;    // spill from char *
    int name;  // spill for char *
    char type; // 1 byte, from char *
    int scope; // 32 bit int  from char *
    int taxid; // 32 bit int , straightforard
    int desc;  // spill from char * 
    int numgenes;

// little tricky here: "hits" is not set at record creation (it is set to null), then, later (i.e when running l2p) ,  
// it is used in processing when the binpath[] data is read in. use this "hits" field for the count of 
// genes that hit this pathway
    int offset2geneids;   // pointer to "numgenes" geneids 

//  hits is a ptr to an array with 
//    firstelemen=[0]=numhits, then the rest of the arry is [1...n] ptrs to struct of generecs (bingentype? right?)
    void *hits;   // used in  l2p for user "hits" to this pathway
};

#define MAXGENENAME 26
// maximum length ARHGAP27P1-BPTFP1-KPNA2P3 = 25

struct bingenetype
{
    int geneid; // entrez gene id
    char hugo[MAXGENENAME];
    char ensembl[MAXGENENAME];
    int pathcount;   // count of paths, ids are in int array famous at "pathplace"
    int pathplace;   // index to path (to a struct binpathouttype record,see above)
    int categories;  // bit patterns 
};

struct genelisttype // used by harvest programs 
{
    int geneid;
    struct genelisttype *n;
};

struct bstype // biosystems id and info - input into this array  -- used by harvest programs 
{
    int bsid;
    int  category; // use CAT_ bitpattern defines (above) 
    char *accession;
    char *name;
    char *type;
    char *scope;
    int taxid;
    char *desc;
    int redundant;       // flag for checking to see if this pathway is duplicated by another pathway 
// next two fields get values from other file
    int numgenes;       // "count of" in next line of code line (i.e. number of genes)
    struct genelisttype *geneslinkedlist; // a linked list
};

struct hit_type // a linked list
{
    char *genename; 
    struct hit_type *n;
};


struct hugo_type
{
    char *hugo;
    struct bingenetype *generec_ptr; 
    int status;   // this can be used for various purposes, initial reason is to use for "universe" masking 
};

struct genetype // from ncbi
{
    int geneid;
    char *hugo;
    char *ensembl;
    int categories;
};



// pathway commons
#define chemical_affects                             1
#define in_complex_with                              2
#define catalysis_precedes                           4
#define controls_expression_of                       8
#define controls_state_change_of                     16
#define controls_production_of                       32
#define consumption_controlled_by                    64
#define controls_phosphorylation_of                  128
#define used_to_produce                              256
#define transport                                    512
#define reacts_with                                  1024
#define interacts_with                               2048
#define reference                                    4096
#define other                                        8192

#define MAXPC 1578297
        // latest 1578297  16637893 850016130 PathwayCommons9.All.hgnc.txt

struct pctype
{
    int ID_Interactor_A;
    int ID_Interactor_B;
    int interaction_type;
};

#define MAXBIOGRID 303568

// bits for "interaction_type" field ...                         // count name  
#define association                                           1  // 8931 psi-mi:"MI:0914(association)"
#define colocalization                                        2  // 44101 psi-mi:"MI:0403(colocalization)"
#define synthetic_genetic_interaction_defined_by_inequality   4  // 50045 psi-mi:"MI:0794(synthetic genetic interaction defined by inequality)"
#define suppressive_genetic_interaction_defined_by_inequality 8  // 197811 psi-mi:"MI:0796(suppressive genetic interaction defined by inequality)"
#define direct_interaction                                    16 // 206875 psi-mi:"MI:0407(direct interaction)"
#define physical_association                                  32 // 329721 psi-mi:"MI:0915(physical association)"
#define additive_genetic_interaction_defined_by_inequality    64 // 535593 psi-mi:"MI:0799(additive genetic interaction defined by inequality)"
struct biogridtype
{
    int ID_Interactor_A;
    int ID_Interactor_B;
    int interaction_type;
};

void category_set_all(int *pat);
void category_code_to_string(int cat,char puthere[]);
int string_to_category_code(char cats[]);
void categories_pattern_to_strings(int cat,char puthere[]);



