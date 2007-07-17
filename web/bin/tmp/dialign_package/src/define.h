#define INT16 short int
#ifdef __alpha
#define INT32 int
#else
#define INT32 long int
#endif
#define REAL32 float
#define REAL64 double

#ifndef TRUE
#define TRUE  1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef MAC
#define Malloc malloc
#define Calloc calloc
#define Realloc realloc
#define Free free
#endif

#define STATUS 0x10000002
