      SUBROUTINE convert_to_real(INTEGRIN,REALOUT)

C***********************************************************************
C  Purpose:  properly interpret the bit string for a input
C  floating point number to obtain the correct floating point number on 
C  the target system.   
C
C  Method: The bit string that represents the VAX floating point number is
C  passed into the routine as an integer.  Note the bytes have not been 
C  swapped.  The data is assumed to be in Big Endian IEEE format.  The 
C  4 bytes of the input integer are swaped, using MVBITS, so the result 
C  is a Little Endian representation.  This bit string is then transfered
C  to a  variable of type real using the function transfer.  Once the bit 
C  string is transfered it is properly interpreted as a real number. 
C
C  References:
C  1.  Final Interface Specification for the Satellite Data Handling 
C  System Communications Network (SDHS-COMNET), 01 February 1988
C  2.  AFGWC/DONS PV-Wave program auto_convert.pro 
C  3.  VAX Fortran Reference Manual from the SDHS programmer library  
C***********************************************************************

      INTEGER   INTEGRIN   !The input integer that contains the bit string
      REAL      REALOUT    !The output floating point number
c     INTEGER*4 INTEGRTEMP  !A temporary integer that contains the swapped
                            !bytes of the input interger INTEGRIN 

C***********************************************************************
C  Call the intrinsic routine mvbits to swap bytes 4 to 1, 1 to 4, 2 to 3
C  and 3 to 2.  Call the intrinsic function to transfer the bit string
C  from the integer where the bytes have been swaped to a variable of 
C  type real.  
C***********************************************************************

c     call mvbits(integrin,24,8,integrtemp,0)
c     call mvbits(integrin,16,8,integrtemp,8)
c     call mvbits(integrin,8,8,integrtemp,16) 
c     call mvbits(integrin,0,8,integrtemp,24)
c     realout=transfer(integrtemp,realout)
      realout=transfer(integrin,realout)

      RETURN
      END

