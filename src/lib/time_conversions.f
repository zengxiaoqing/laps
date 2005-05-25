
      subroutine afwa_julhr_i4time(I_A1JUL,I_A1MIN,i4time)

!     I_A1JUL is number of hours since Dec 31, 1967 at 00z
!     This is converted to i4time, number of sec since Jan 1, 1960 at 00z
      i4time_hr  = I_A1JUL * 3600 + (8*365 - 1 + 2) * 86400
      i4time_min = I_A1MIN*60
      i4time     = i4time_hr + i4time_min

      return
      end

      subroutine jd_to_i4time(jd,i4time,istatus)

      double precision jd, jd1960
      integer i4time

cdoc  Converts Julian Day (Number of days since Jan 1 4713BCE to i4time)

!     Author: Steve Albers 2005

      jd_1960 = 2436934.5

      days_since_1960 = jd - jd_1960

      i4time = int(days_since_1960 * 86400.)

      istatus = 1

      return
      end
