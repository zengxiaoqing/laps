c
c
      function wfo_fname13_to_fname9(wfo_fname13)
c
c..... Routine converts wfo filename 'yyyymmdd_hhmm' to 'yyjjjhhmm' filename.
c
c  JSmart   8-96  Stole this from JSmart who previously stole it from
c                 PStamus and modified to generate the wfo filename time.
c
      character*2 cyr
      character*3 cjjj
      character wfo_fname13*13, wfo_fname13_to_fname9*9
      integer imon_a(12), imon(12)
      data imon_a/0,31,59,90,120,151,181,212,243,273,304,334/
c
c..... First, read the LAPS time and get julian days.
c
      read(wfo_fname13(3:4),11) iyr
      read(wfo_fname13(5:6),11) imm
      read(wfo_fname13(7:8),11) idy
 11   format(i2)
c
c..... Check for leap year.
c
      do i=1,12
         imon(i) = imon_a(i)
      enddo !i
      if(iyr .eq. (iyr/4)*4) then
         do i=3,12
            imon(i) = imon(i) + 1
         enddo !i
      endif
c
c..... Convert month, day to julian day.
c
      jjj = imon(imm) + idy
c
c..... Now write out the time.
c
      write(cyr,12) iyr
 12   format(i2)
      write(cjjj,13)jjj
 13   format(i3)

      if(cyr(1:1) .eq. ' ') cyr(1:1)='0'
      if(cyr(2:2) .eq. ' ') cyr(2:2)='0' 
      if(cjjj(1:1).eq. ' ') cjjj(1:1)= '0'
      if(cjjj(2:2).eq. ' ') cjjj(2:2)= '0'
      if(cjjj(3:3).eq. ' ') cjjj(3:3)= '0'

      wfo_fname13_to_fname9 = cyr//cjjj//wfo_fname13(10:13)
c
      return
      end
c
c
c ----------------------------------------------------------------
c
      function fname9_to_wfo_fname13(a9_time_in)
c
c..... Routine to convert LAPS 'yyjjjhhmm' time to 'yyyymmdd_hhmm' time.
c
c  JSmart   8-96  Stole this from PStamus and modified to generate the
c                 wfo filename time.
c
      character*2 cyr,cmm,cdy,chh,cmin
      character*4 cyyyy
      character*9 a9_time
      character*9 a9_time_in
      character*13 fname9_to_wfo_fname13
      integer imon_a(12), imon(12)
      data imon_a/0,31,59,90,120,151,181,212,243,273,304,334/

      a9_time = a9_time_in
c
c..... First, read the LAPS time and get julian days.
c
      read(a9_time,11) iyr, jjj, ihh, imin
 11   format(i2,i3,2i2)
c
c..... Check for leap year.
c
      do i=1,12
         imon(i) = imon_a(i)
      enddo !i
      if(iyr .eq. (iyr/4)*4) then
         do i=3,12
            imon(i) = imon(i) + 1
         enddo !i
      endif
c
c..... Convert julian day to month, day.
c
      do i=12,1,-1
         kk = jjj - imon(i)
         if(kk .gt. 0) then
            imm = i
            idy = kk
            go to 200
         elseif(kk .eq. 0) then
            imm = i - 1
            idy = jjj - imon(imm)
            go to 200
         endif
      enddo !i
      imm = 1
      idy = jjj
c
 200  continue
c
c..... Now write out the time.
c
      write(cyr,12) iyr
      write(cmm,12) imm
      write(cdy,12) idy
      write(chh,12) ihh
      write(cmin,12)imin
 12   format(i2)

      if(cmm(1:1) .eq. ' ') cmm(1:1) = '0'
      if(cdy(1:1) .eq. ' ') cdy(1:1) = '0'
      if(chh(1:1) .eq. ' ') chh(1:1) = '0'
      if(chh(2:2) .eq. ' ') chh(2:2) = '0'
      if(cmin(1:1).eq. ' ') cmin(1:1)= '0'
      if(cmin(2:2).eq. ' ') cmin(2:2)= '0'
      if(iyr.lt.90)then
         cyyyy='20'//cyr
      else
         cyyyy='19'//cyr
      endif

      fname9_to_wfo_fname13 = cyyyy//cmm//cdy//'_'//chh//cmin
c
      return
      end
c
c ----------------------------------------------------------------
c
       function cvt_wfo_fname13_i4time(wfo_fname13)

       character*13 wfo_fname13
       character*9  wfo_fname13_to_fname9
       character*9  fname9
       integer*4    cvt_wfo_fname13_i4time
       integer*4    i4time
       INTEGER*4 INT_FILE(9), I, NYEAR, JDAY, NHOUR, MIN, MONTH, NDAY
       INTEGER*4 I4TIME_INT_LP

c
c first convert wfo filename *13 to filename 'yyjjjhhmm'
c
       fname9 = wfo_fname13_to_fname9(wfo_fname13)
c
c convert fname9 to i4time
c
       call i4time_fname_lp (fname9, i4time, istatus)
       if(istatus.ne.1)then
          write(6,*)'Error converting fname9 to i4time'
          write(6,*)'Conversion incomplete in cvt_wfo_fname13_i4time'
       endif

       cvt_wfo_fname13_i4time = i4time
C
C      Successful return
C
       ISTATUS = 1
       RETURN
       end
c
c ---------------------------------------------------------------
c
       function cvt_i4time_wfo_fname13(i4time)

       character*9  fname9
       character*13 cvt_i4time_wfo_fname13
       character*13 fname9_to_wfo_fname13
       integer*4    i4time
       integer*4    i4time_temp
c
c first convert i4time to ascii filename *9
c
       i4time_temp=i4time
       call make_fnam_lp (i4time_temp, fname9, istatus)
       if(istatus.ne.1)then
          write(6,*)'Error converting i4time to fname9'
          write(6,*)'Conversion incomplete in cvt_i4time_wfo_fname13'
       endif
c
c convert fname9 to wfo_fname13
c
       cvt_i4time_wfo_fname13 = fname9_to_wfo_fname13(fname9)

       return
       end
