
        function a10_to_a9(a10_time,istatus)

!       Convert a10_time (yyMMddhhmm) to a9_time (yydddhhmm)
!       Steve Albers 1998

        character*10 a10_time
        character*9 a10_to_a9, a8_to_a9, a9_time
        character*8 a8_time

        a8_time = a10_time(1:8)
        a9_time = a8_to_a9(a8_time) 
        a9_time = a9_time(1:7)//a10_time(9:10)

        a10_to_a9 = a9_time

        istatus = 1

        return
        end


        function a8_to_a9(a8_time)

!       Convert a8_time (yyMMddhh) to a9_time (yydddhhmm)
!       Steve Albers 1998

        character*9 a8_to_a9
        character*9 a9
        character*8 a8_time

        integer*4 imon_a(12)

        data imon_a/0,31,59,90,120,151,181,212,243,273,304,334/

        read(a8_time,1)iyr,imn,idy,ih
1       format(i2,i2,i2,i2)

        id = imon_a(imn)

        idays = id + idy

!       Decide whether to add a day for leap year.
        if(iyr .eq. (iyr / 4) * 4 )then
            if(imn .ge. 3)then
                idays = idays + 1
            endif
        endif

        write(a9,2)iyr,idays,ih,im
2       format(i2,i3,i2,i2)

        if(a9(1:1) .eq. ' ')a9(1:1) = '0'    ! y
        if(a9(3:3) .eq. ' ')a9(3:3) = '0'    ! d
        if(a9(4:4) .eq. ' ')a9(4:4) = '0'    ! d
        if(a9(6:6) .eq. ' ')a9(6:6) = '0'    ! h
        if(a9(8:8) .eq. ' ')a9(8:8) = '0'    ! m

        a8_to_a9 = a9

        return
        end
c
c
      function a9_to_a8(a9_time)
c
c..... Routine to convert LAPS 'yyjjjhhmm' time to 'yymmddhh' time.
c..... Corrected for Y2K.  P. Stamus, NOAA/FSL   Oct 1998
c
      character a9_time*9, a9_to_a8*8, a8*8
      integer imon_a(12), imon(12)
      data imon_a/0,31,59,90,120,151,181,212,243,273,304,334/
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
      write(a8,12) iyr, imm, idy, ihh
 12   format(4i2)
c
      if(a8(1:1) .eq. ' ') a8(1:1) = '0'
      if(a8(3:3) .eq. ' ') a8(3:3) = '0'
      if(a8(5:5) .eq. ' ') a8(5:5) = '0'
      if(a8(7:7) .eq. ' ') a8(7:7) = '0'
c
      a9_to_a8 = a8
c
      return
      end
