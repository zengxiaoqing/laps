      SUBROUTINE CNVVXDB(INTEGRIN1,INTEGRIN2,REALOUT)

C***********************************************************************
C  Purpose: Properly interpret the bit string that represents a VAX 
C  double precision number to obtain the correct number on the target 
C  system.   
C
C  Method: The bit string that represents the VAX double precision number
C  is passed into the routine as two integers.  Note the bytes have not
C  been swapped.  The sign bit, the bits that represent the exponent or 
C  characteristic, and the bits that represent the fraction are picked off
C  using the standard Fortran function IBITS.  The following table show how
C  the bits and what they represent (sign, exponent etc) map to most unix 
C  systems.
C    
C                        first integer   
C          byte1       byte2       byte3     byte4
C   VAX    07 06-00    15 14-8     23-16     31-24 
C          e1 f7-f1     s e8-e2    f15-f8    f23-f16
C   UNIX   31 30-24    23 22-16    15-08     07-00
C
C                        second integer   
C          byte1      byte2        byte3     byte4
C   VAX    07-00      15-08        23-16     31-24
C          f31-24     f39-f32      f47-40    f55-48  
C   UNIX   31-24      23-16        15-08     07-00
C
C  Here f1 stands for the first bit that is used to represent the 
C  fraction, f2 for the second bit used to represent the fraction, etc., 
C  e1 represents the first exponent bit, e2 the second, and s represents 
C  the sign bit.  A VAX floating point number is represented by :
C
C  floating point number = -1**signbit * fraction * 2 **(EXPONENT-128)
C    
C  !!!!!!Note: to save space the VAX system assumes the most significant
C              fraction bit is set.  This fact is taken into account when
C              computing the fraction.  Also, as appears to be standard
C              practice, the exponent is biased by a fix number(128 in 
C              this case). 
C    
C  References:
C  1.  Final Interface Specification for the Satellite Data Handling 
C  System Communications Network (SDHS-COMNET), 01 February 1988
C  2.  AFGWC/DONS PV-Wave program auto_convert.pro 
C  3.  VAX Fortran Reference Manual from the SDHS programmer library  
C***********************************************************************

      INTEGER INTEGRIN1,INTEGRIN2 !input integers that contain the 
                                    !the bit strings
      REAL*8 REALOUT     !The output double precision number
      INTEGER SIGNBIT  !the sign bit
      INTEGER CHARACTS !the integer representation of the exponent
                         !before it is biased


      REAL*8 M1,M2,M3,M4 ! components of the fraction;M1 is most 
                         !significant

      REAL*8 M5,M6,M7    ! components of the fraction;M1 is most 
                         !significant

      REAL*8 FRACTION    !the fraction portion of the double precision
                         !number 

      REAL*8 EXPONENT    !the exponent portion of the double precision 
                         !number

C***********************************************************************
C  Pick off the sign bit and convert it to an integer (a zero or a one)
C***********************************************************************

      SIGNBIT = IBITS(INTEGRIN1,23,1)
      

C***********************************************************************
C  Pick off the exponent bits (7 bits from one byte, 1 bit from another)
C  convert the two bit strings to their integer representation and add
C  to get the unbiased exponent value  
C***********************************************************************

      CHARACTS = IBITS(INTEGRIN1,16,7)*2 + IBITS(INTEGRIN1,31,1)

C***********************************************************************
C  Bias the exponent
C***********************************************************************
     
      EXPONENT = DFLOAT (CHARACTS - 128)

C***********************************************************************
C Pick off the bits strings from the three difference bytes that contain
C fraction bits.  Convert the bit strings to the appropriate fraction.
C To account for the fact that the most significant fraction bit is
C assumed to be set, 2**7 is added to the most significant fractional 
C component M1. 
C***********************************************************************
     
      M1 = DFLOAT( (IBITS( INTEGRIN1,24,7 ) + 2**7 ) )/(2.0D0)**8
      M2 = DFLOAT(IBITS(INTEGRIN1,8,8))/(2.0D0)**24 
      M3 = DFLOAT(IBITS(INTEGRIN1,0,8))/(2.0D0)**16
      M4 = DFLOAT( IBITS(INTEGRIN2,16,8) )/(2.0D0)**32
      M5 = DFLOAT(IBITS(INTEGRIN2,24,8))/(2.0D0)**40
      M6 = DFLOAT(IBITS(INTEGRIN2,0,8))/(2.0D0)**48
      M7 = DFLOAT(IBITS(INTEGRIN2,8,8))/(2.0D0)**56

C***********************************************************************
C  Compute the fraction portion of the double precision number by adding
C  the fractional components
C***********************************************************************

      FRACTION = M1+M2+M3+M4+M5+M6+M7

C***********************************************************************
C  Compute the double precision number according to the defintion of a 
C  VAX double precision number.
C***********************************************************************
     
      REALOUT= DFLOAT(1-2*SIGNBIT)*FRACTION*2.0D0**EXPONENT

C***********************************************************************
C  Since one must assume that the most significant fraction bit is set.
C  Zero in the above computation will be computed as .5 i.e. all fraction
C  bits were zero but the implied bit resulted in the fraction being .5.
C  To ensure zero is properly represented the following code checks for 
C  when the fraction is .5 and set the floating point number to zero
C***********************************************************************

      IF(FRACTION .EQ. .5)THEN
        REALOUT=0.0
      ENDIF
 
      RETURN
      END

