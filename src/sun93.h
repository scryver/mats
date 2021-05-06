
#define FLT_UWORD_IS_FINITE(x)     ((x) <  0x7f800000L)
#define FLT_UWORD_IS_NAN(x)        ((x) >  0x7f800000L)
#define FLT_UWORD_IS_INFINITE(x)   ((x) == 0x7f800000L)
#define FLT_UWORD_MAX              0x7fffffff
#define FLT_LARGEST_EXP            (FLT_UWORD_MAX >> 23)
#define FLT_SMALLEST_EXP           -22

#define FLT_UWORD_IS_ZERO(x)       ((x) == 0)
#define FLT_UWORD_IS_SUBNORMAL(x)  ((x) <  0x00800000L)

#define OVERFLOW_INT   50000
#define UNDERFLOW_INT -50000
