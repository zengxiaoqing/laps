      SUBROUTINE CNVTINT(INTIN)

C***********************************************************************
C  Purpose: Properly interpret the bit string that represents a negative 
C  integer in 2's complement form.   
C
C  Method: If the sign bit of the input integer is set, then the number
C  is negative and is assumed to be in 2's complement form.  The 
C  appropriate representation of the integer on the target system is 
C  obtained by inverting the 2's complement form using the standard 
C  fortran functions IBCLR (clear a bit) and NOT (a bitwise complement)
C    
C  Reference: Operations Ground Equipment Interface Specification DRL 
C             504-02-1 Part1, Document No. E007020;Space Systems/Loral
C   
C***********************************************************************

      INTEGER INTIN
      INTEGER SIGNBIT

C***********************************************************************
C  Pick off the sign bit and convert it to an integer (a zero or a one)
C***********************************************************************

      SIGNBIT = IBITS(INTIN,31,1)

C***********************************************************************
C  Check to see if the number is negative, if so, invert the 2's 
C  complement notation.  The inversion is done by clearing the sign bit,
C  subtracting one , flipping the bits, and reclearing the sign bit. 
C  The result is the magnitude of the negative number
C***********************************************************************

      IF(SIGNBIT .EQ. 1) THEN
        INTIN = IBCLR(INTIN,31)
        INTIN = NOT((INTIN-1))
        INTIN = IBCLR(INTIN,31)
      ENDIF

C***********************************************************************
C  Compute the integer based on the sign bit and the magnitude of the
C  input number. 
C***********************************************************************

      INTIN = (1-2*SIGNBIT)*INTIN

      RETURN
      END

