      FUNCTION BYTESWP4(VAXINT)
C***********************************************************************
C  Purpose:  Perform the functional equivalent of a byte swap of a vax 
C  four byte integer
C
C  Method:  The bit strings contained in the bytes of the input integer
C  VAXINT are parsed into the integer variables BYTE1, BYTE2, BYTE3, 
C  and BYTE4.  These bits strings are used to compute the integer that 
C  they represent.  This method is used rather than just swapping bytes
C  to make the program system independent.  No knowledge of how the
C  target system stores its integers is needed.  One only needs to know 
C  how the VAX system stores integers.  For a VAX system the least
C  significant bits are stored in byte 1 and the most significant in byte
C  4--The reverse is true for most unix system.      
C 
C  References:
C  1.  Final Interface Specification for the Satellite Data Handling 
C  System Communications Network (SDHS-COMNET), 01 February 1988
C  2.  AFGWC/DONS PV-Wave program auto_convert.pro 
C  3.  VAX Fortran Reference Manual from the SDHS programmer library  
C***********************************************************************

      IMPLICIT NONE
      integer byteswp4
      INTEGER VAXINT
      INTEGER BYTE1
      INTEGER BYTE2
      INTEGER BYTE3
      INTEGER BYTE4

C***********************************************************************
C  Pick off 8-bit bytes, using the standard Fortran routine IBITS, from
C  the input variable and store the information in integer variables 
C  The most significant bits are stored in BYTE4, the next most 
C  significant in BYTE3, the next in BYTE2, and the least significant in
C  in BYTE1.  The input VAX bit string has the most significant bits in
C  the fourth byte, the next most significant in the third byte, the 
C  the next in the second byte and the least significant in the first
C  byte.  Since most unix sytem expect the bytes to be stored in reverse
C  order, a call to IBITS to get bits 0 to 8 on a unix machine means the
C  routine goes to the 4th word to read the the bit string.  Becuase the 
C  the 4th byte on a vax machine has the most significant bits a call to 
C  IBITS with the argument 0, 8 gets the most significant bits of the 
C  byte.  Sound complicated? It is. The bit string diagram given below
C  may help.
C   
C        byte1   byte2   byte3   byte4
C   VAX  07-00   15-08   23-16   31-24 
C   UNIX 31-24   23-16   15-08   07-00
C     
C***********************************************************************

      byte1=0
      byte2=0
      byte3=0
      byte4=0

      byte4 = IBITS(VAXINT,0,8)  
      byte3 = IBITS(VAXINT,8,8)
      byte2 = IBITS(VAXINT,16,8)
      byte1 = IBITS(VAXINT,24,8)

C***********************************************************************
C  Recontruct the integer by multipling the integer representation of an
C  8-bit byte by the appropriate power of 2 and adding the four bytes.
C  For example:  BYTE2 represents the second least significant string of
C  8 bits and must be multiplied by 2**8 to obtain the correct integer
C  representation. BYTE3 must be multiplied by 2**16 and BYTE4 by 2**24
C  to obtain the proper integer representation.  Their sum represents 
C  the integer.  Note: if negative values are stored in 2's or 1's 
C  complement form then the integer computed will not be a correct
C  representation, but the bit string will be equal to the input VAX
C  string with the bytes swapped.
C***********************************************************************
      
      byteswp4 = BYTE1 + BYTE2*256 + BYTE3*65536 + BYTE4*16777216

      RETURN
      END

