
      subroutine afwa_julhr_i4time(I_A1JUL,I_A1MIN,i4time)

!     I_A1JUL is number of hours since Dec 31, 1967 at 00z
!     This is converted to i4time, number of sec since Jan 1, 1960 at 00z
      i4time_hr  = I_A1JUL * 3600 + (8*365 - 1 + 2) * 86400
      i4time_min = I_A1MIN*60
      i4time     = i4time_hr + i4time_min

      return
      end

