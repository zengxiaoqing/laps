
      SUBROUTINE CLOCK(H,c5_string)
      REAL*8 H,RHOUR
      CHARACTER*5 c5_string
      INTEGER HOUR,MINUTE

!     Hour is the input time in radians

      RHOUR=DMOD(H*3.819718634D0+48.D0,24.D0)
      HOUR=IDINT(RHOUR)

      MINUTE=NINT((RHOUR-HOUR)*60.D0)

      if(minute .eq. 60)then
          minute = 0
          hour = hour + 1
      endif

      MIN1 = MINUTE/10
      MIN2 = MINUTE - MIN1 * 10

      write(c5_string,1)HOUR,MIN1,MIN2
1     format(i2,':',i1,i1)

      RETURN
      END
