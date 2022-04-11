
struct pwtype
{
    unsigned int category; // top 2 bits are type, rest is category
    char *acc;
    char *name;
    unsigned int numgenes;
    unsigned int pwgenesindex;
#if 0
    unsigned int longdesc;  // offset to overflow text file - reason, not always accessed so not in memory
#endif
};
