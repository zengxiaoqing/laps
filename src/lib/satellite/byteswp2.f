      function BYTESWP2(VAXINT)

C***********************************************************************
C  Purpose:  Perform the functional equivalent of a byte swap of a VAX 
C  two byte integer
C
C  Method:  The bit strings contained in the bytes of the input integer
C  VAXINT are parsed into the integer variables BYTE1,  and BYTE2 
C  These bits strings are used to compute the integer that 
C  they represent.  This method is used rather than just swapping bytes
C  to make the program system independent.  No knowledge of how the
C  target system stores its integers is needed.  One only needs to know 
C  how the VAX system stores integers.  For a VAX system the least
C  significant bits are stored in byte 1 and the most significant in byte
C  2--The reverse is true for most unix system. 
C  
C  References:
C  1.  Final Interface Specification for the Satellite Data Handling 
C  System Communications Network (SDHS-COMNET), 01 February 1988
C  2.  AFGWC/DONS PV-Wave program auto_convert.pro 
C  3.  VAX Fortran Reference Manual from the SDHS programmer library  
C***********************************************************************
      IMPLICIT NONE
      integer byteswp2
      INTEGER VAXINT
      INTEGER BYTE1
      INTEGER BYTE2

C***********************************************************************
C  Pick off 8-bit bytes, using the standard Fortran routine IBITS, from
C  the input variable and store the information in an integer variables 
C  The most significant bits are stored in BYTE2, and the
C  least significant in BYTE1.  The input VAX bit string has the least
C  significant bits in the first byte and the most significant in
C  the second byte.  Since most unix sytem expect the bytes to be stored
C  in reverse order, a call to IBITS to get bits 0 to 8 on a unix machine
C  means the routine goes to the 2th word to read the the bit string.  
C  Becuase the the 2th byte on a vax machine has the most significant 
C  bits a call to IBITS with the argument 0, 8 gets the most significant
C  bits of the integer.  Sound complicated? It is. The bit string 
C  diagram given below may help.
C   
C        byte1   byte2   
C   VAX  07-00   15-08    
C   UNIX 15-08   07-00 

C***********************************************************************

cc      byte2 = IBITS(VAXINT,0,8)
cc      byte1 = IBITS(VAXINT,8,8)
      byte1=0
      byte2=0
      byte2 = IBITS(VAXINT,16,8)
      byte1 = IBITS(VAXINT,24,8)
     

***********************************************************************
C  Recontruct the integer by multipling the integer representation of an
C  8-bit byte by the appropriate power of 2 and adding the two bytes.
C  For example:  BYTE2 represents the most significant string of
C  8 bits and must be multiplied by 2**8 to obtain the correct integer
C  representation.   The  sum of BYTE1 and BYTE2 represents the 
C  integer.  Note: if negative values are stored in 2's or 1's 
C  complement form then the integer computed will not be a correct
C  representation, but the bit string will be equal to the input VAX
C  string with the bytes swapped.
C*********************************************************************** 
      if(byte2.lt.128) then
         BYTESWP2 = BYTE1 + BYTE2*256 
      else
         write(6,*) 'WARNING: BYTESWP2 incountered value < 0'
     +             ,' the value of this translation is doubtful'
         BYTESWP2 = -BYTE1 + (128-BYTE2)*256
      endif

      RETURN
      END
