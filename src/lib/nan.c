#include <config.h>

/* #include <nan.h> */
/* this is the NaN macro from the SGI library
   with the test value modified from 0x7ff to 0x7fb
   to work (apparently) with stuff passed from fortran
*/

/* defined for the IBM jpe */
#if defined(FORTRANUNDERSCORE)
#  define nan nan_
#elif defined(FORTRANDOUBLEUNDERSCORE)
#  define nan nan__
#elif defined(FORTRANCAPS)
#  define nan NAN
#endif


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
#ifdef IRIX
#include <ieeefp.h>
#include <fp_class.h>
#define IsNANorINF(X)  (!finite((double) *x))
#endif
#ifdef CRAYT3E
#include <fp.h>
#undef NAN
#define IsNANorINF(X) (!isfinite( *x))
#endif
#ifndef IsNANorINF
#define IsNANorINF(X)  (((dnan *)&(X))->nan_parts.exponent == 0x7f8)
#endif


int nan(float *x)
{
  /*
  printf("%e %d %x\n",*x, IsNANorINF(*x),((dnan *)&(*x))->nan_parts.exponent );
#ifdef IRIX
  if(IsNANorINF(*x)){
        printf("fp class: %d\n",fp_class_f(*x));
  }
#endif
  */
  return(IsNANorINF(*x));
}







