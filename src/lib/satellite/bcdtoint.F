#include <config.h>
        subroutine bcd_to_int (re8,int4)

c       routine takes real*8 and converts into BCD form

        implicit none

        integer re8(2),int4(2),i,j
        character*8 string1

        j=2
        do i=1,2
           
           write(string1,'(4z2.2)') ibits(re8(i),24,8),
     +                              ibits(re8(i),16,8),
     +                              ibits(re8(i),8,8),
     +                              ibits(re8(i),0,8)

#ifdef WORDS_BIGENDIAN
           read(string1,*) int4 (i)
#else
           print*,'re8 =',re8,ibits(re8(i),0,8),ibits(re8(i),8,8)
           print*,'strint1 =',string1

           read(string1,*) int4(j)
           j=1
#endif
        enddo
        return
        end
