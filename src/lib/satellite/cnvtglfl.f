      SUBROUTINE CNVTGLFL(INTEGRIN,REALOUT)

C***********************************************************************
C  Purpose: Properly interpret the bit string that represents a Gould 
C  floating point number to obtain the correct floating point number on 
C  the target system.   
C
C  Method: The bit string that represents the Gould floating point number
C  is passed into the routine as an integer. The sign bit, the bits that
C  represent the exponent or characteristic, and the bits that represent 
C  the fraction are picked off using the standard Fortran function IBITS.
C  The following table show the bits and what they represent (sign, 
C  exponent etc).
C    
C          byte1       byte2    byte3     byte4
C          31 30-24    23-16    15-08     07-00
C          s  e7-e1    f8-f1    f16-f9    f24-f17
C
C    
C
C  Here f1 stands for the first bit that is used to represent the 
C  fraction, f2 for the second bit used to represent the fraction, etc., 
C  e1 represents the first exponent bit, e2 the second, and s represents 
C  the sign bit.  A Gould floating point number is represented by :
C
C  floating point number = -1**signbit * fraction * 2 **((EXPONENT-64)*4)
C    
C  !!!!!!Note: for negative numbers the bit string is stored in 2's
C              complement form.  This form must be inverted before the
C              bit string can be properly interpreted
C    
C  Reference: Operations Ground Equipment Interface Specification DRL 
C             504-02-1 Part1, Document No. E007020;Space Systems/Loral
C   
C***********************************************************************

      INTEGER INTEGRIN !The input integer that contains the bit string
      REAL*8 REALOUT !The output floating point number
      INTEGER SIGNBIT !the sign bit
      INTEGER CHARACTS !the integer representation of the exponent 
                         !before it is biased

      REAL*8 M1,M2,M3 !components of the fraction;M1 is most significant
      REAL*8 FRACTION !the fraction portion of the floating point number
      REAL*8 EXPONENT !the exponent portion of the floating point number

C***********************************************************************
C  Pick off the sign bit and convert it to an integer (a zero or a one)
C***********************************************************************
      
      SIGNBIT = IBITS(INTEGRIN,31,1)

C***********************************************************************
C  Check to see if the number is negative, if so, invert the 2's 
C  complement notation.  The inversion is done by clearing the sign bit,
C  subtracting one and flipping the bits
C***********************************************************************

      IF (SIGNBIT .EQ. 1) THEN
        INTEGRIN = IBCLR(INTEGRIN,31)
        INTEGRIN = NOT((INTEGRIN-1))
      ENDIF

C***********************************************************************
C  Pick off the exponent bits (7 bits) convert the two bit strings to
C  their integer representation and add to get the unbiased exponent
C  value  
C***********************************************************************

      CHARACTS = IBITS(INTEGRIN,24,7)

C***********************************************************************
C  Bias the exponent
C***********************************************************************

      EXPONENT = DFLOAT ((CHARACTS - 64)*4 )

C***********************************************************************
C Pick off the bits strings from the three difference bytes that contain
C fraction bits.  Convert the bit strings to the appropriate fraction.
C***********************************************************************

      M1 = DFLOAT( (IBITS( INTEGRIN,16,8 ) ) )/(2.0D0)**8
      M2 = DFLOAT(IBITS(INTEGRIN,8,8))/(2.0D0)**16 
      M3 = DFLOAT(IBITS(INTEGRIN,0,8))/(2.0D0)**24

C***********************************************************************
C  Compute the fraction portion of the floating point number by adding
C  the three fractional components
C***********************************************************************
 
      FRACTION = M1 +M2 + M3

C***********************************************************************
C  Compute the floating point number according to the defintion of a 
C  Gould floating point number.
C***********************************************************************

      REALOUT= FLOAT(1-2*SIGNBIT)*FRACTION*2.0**EXPONENT


      RETURN
      END
