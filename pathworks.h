
#define MAXBSID 851568
#define MAXGENE 100000 

#define NEWWAY 1

struct binpathouttype     // the "binary" pathway information file 
{
    int bsid;         //32 bit int
    char source[10];  // fixed size - null terminated
    char accession[20];  // fixed size - null terminated 
    int name;  // spill for char *
    char type; // 1 byte, from char *
    int scope; // 32 bit int  from char *
    int taxid; // 32 bit int , straightforard
    int desc;  // spill from char * 
    int numgenes; 
    int offset2geneids; // pointer to "numgenes" geneids 
// little tricky here: "hits" is not really set at output record creation (set to null), then, later,  
// it is used in processing when the binpath[] data is read in. use this "hits" field for the count of 
// genes that hit this pathway
// ALSO: hits is a ptr to an array with 
//    firstelemen=[0]=numhits, then the rest of the arry is [1...n] ptrs to struct of generecs (bingentype? right?)
    void *hits;   // used in  l2p for user "hits" to this pathway
};

#define MAXGENENAME 34
struct bingenetype
{
    int geneid;
    char hugo[MAXGENENAME];
    char ensembl[MAXGENENAME];
    int pathcount;   // count of paths, ids are in int array famous at "pathplace"
    int pathplace;   // index to path (to a struct binpathouttype record,see above)
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
              // was 713209  24751583 592549368 PathwayCommons.8.Detailed.EXTENDED_BINARY_SIF.hgnc.txt

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

