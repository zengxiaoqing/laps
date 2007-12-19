      SUBROUTINE convert_to_double(INTEGRIN1,INTEGRIN2,REALOUT)

C***********************************************************************
C  Purpose: Properly interpret the bit string that represents a VAX 
C  double precision number to obtain the correct number on the target 
C  system.   
C
C  Method: The bit string that represents the VAX double precision 
C  number is passed into the routine as two integers.  Note the bytes 
C  have not been swapped.  The data is assumed to be in Big Endian 
C  IEEE format.  The 4 bytes of the input integers are swaped, using 
C  MVBITS, so the result is a Little Endian representation stored in an
C  integer array of dimension two.  This bit string is then transfered 
C  to a  variable of type double precision (REAL*8) using the routine 
C  transfer.  Once the bit string is transfered it is properly 
C  interpreted as a double precision number. 
C
C    
C  References:
C  1.  Final Interface Specification for the Satellite Data Handling 
C  System Communications Network (SDHS-COMNET), 01 February 1988
C  2.  AFGWC/DONS PV-Wave program auto_convert.pro 
C  3.  VAX Fortran Reference Manual from the SDHS programmer library  
C***********************************************************************

      INTEGER   INTEGRIN1,INTEGRIN2 !input integers that contain the 
                                    !the bit strings
      REAL*8 REALOUT     !The output double precision number
      INTEGER,dimension(2) ::  integrtemp !A working integer array 
                                          !that holds the double
                                          !precision bit string after
                                          !swaping of the words and 
                                          
C***********************************************************************
C  On Big Endian Systems the first word contains the least significant
C  bits.  On a Little Endian System the second (high address 32 bit word)
C  contains the least significant bits.  So in addtion to swaping bytes
C  within a word the word themselves must be swapped.  Both actions are
C  accomplished with the use of MVBITS
C***********************************************************************                                          !bytes

c     call mvbits(integrin1,24,8,integrtemp(2),0)
c     call mvbits(integrin1,16,8,integrtemp(2),8)
c     call mvbits(integrin1,8,8,integrtemp(2),16)
c     call mvbits(integrin1,0,8,integrtemp(2),24)
c     call mvbits(integrin2,24,8,integrtemp(1),0)
c     call mvbits(integrin2,16,8,integrtemp(1),8)
c     call mvbits(integrin2,8,8,integrtemp(1),16)
c     call mvbits(integrin2,0,8,integrtemp(1),24)

      integrtemp(2)=integrin1
      integrtemp(1)=integrin2
      realout=transfer(integrtemp,realout)
 
      RETURN
      END

