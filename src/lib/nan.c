#include <config.h>

#if defined(FORTRANUNDERSCORE)

#define nan nan_

#elif defined(FORTRANDOUBLEUNDERSCORE)

#define nan nan__

#elif defined(FORTRANCAPS)

#define nan NAN

#endif

/* this is the NaN macro from the SGI library
   with the test value modified from 0x7ff to 0x7fb
   to work (apparently) with stuff passed from fortran
*/



typedef union 
{
         struct 
         {
#ifdef _MIPSEL
            unsigned fraction_low:32;
            unsigned bits:20;
            unsigned exponent :11;
            unsigned sign     : 1;
#else
            unsigned sign     : 1;
            unsigned exponent :11;
            unsigned bits:20;
            unsigned fraction_low:32;
#endif
         } inf_parts;
         struct 
         {
#ifdef _MIPSEL
            unsigned fraction_low: 32;
            unsigned bits     :19;
            unsigned qnan_bit : 1;
            unsigned exponent :11;
            unsigned sign     : 1;
#else
            unsigned sign     : 1;
            unsigned exponent :11;
            unsigned qnan_bit : 1;
            unsigned bits     :19;
            unsigned fraction_low: 32;
#endif
         } nan_parts;
         double d;

} dnan; 


#ifdef rs6000
#define IsNANorINF(X)  (((dnan *)&(X))->nan_parts.exponent == 0x7f8 || ((dnan *)&(X))->nan_parts.exponent == 0x7fc)
#endif
#ifdef hpux 
#define IsNANorINF(X)  (((dnan *)&(X))->nan_parts.exponent == 0x7f8 || ((dnan *)&(X))->nan_parts.exponent == 0x7fa)
#endif
#ifdef sun4
#define IsNANorINF(X)  (((dnan *)&(X))->nan_parts.exponent == 0x7f8 || ((dnan *)&(X))->nan_parts.exponent == 0x7ff)
#endif

int nan(x)
float *x;
{
/*
  printf("%f %d %x\n",*x, IsNANorINF(*x),((dnan *)&(*x))->nan_parts.exponent );
*/
  return(IsNANorINF(*x));
}


