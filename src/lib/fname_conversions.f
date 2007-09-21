        function a10_to_a9(a10_time,istatus)

cdoc    Convert a10_time (yyMMddhhmm) to a9_time (yydddhhmm)
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

        function yr_a10_to_a9(a10_time)

cdoc    Convert a10_time (yyyyMMddhh) to a9_time (yydddhhmm)
!       John Smart 2000

        character*10 a10_time
        character*10 new_a10_time
        character*9  yr_a10_to_a9
        character*9  a10_to_a9, a9_time

        new_a10_time=a10_time(3:10)//'00'         !warning: always assuming 00 minutes
        a9_time = a10_to_a9(new_a10_time,istatus)

        yr_a10_to_a9 = a9_time

        return
        end

        function a9_to_yr_a10_time(a9_time,istatus)

cdoc    Convert a9_time (yydddhhmm) to a10_time (yyyyMMddhh)
!       John Smart 2000

        integer i4time_sys,ilen
        integer len

        character*10 a9_to_yr_a10_time
        character*9  a9_time, a9_time_sys
        character*8  a8_time,a9_to_a8
        character*5  yr_jday
        character*16 asctim_str
        character*2  anal_hr,anal_min

        call get_systime_all(i4time_sys,a9_time_sys,anal_hr,
     1             anal_min,asctim_str,yr_jday,istatus)

        if(istatus.ne.1)then
           print*,'Error returned from get_systime_all'
           return
        endif

        a8_time=a9_to_a8(a9_time)
        ilen=len(asctim_str)
        a9_to_yr_a10_time=asctim_str(ilen-8:ilen-5)//a8_time(3:8)

        istatus = 1

        return
        end


        function a8_to_a9(a8_time)

cdoc    Convert a8_time (yyMMddhh) to a9_time (yydddhhmm)
!       Steve Albers 1998

        character*9 a8_to_a9
        character*9 a9
        character*8 a8_time

        integer imon_a(12)

        data imon_a/0,31,59,90,120,151,181,212,243,273,304,334/

        read(a8_time,1,err=2)iyr,imn,idy,ih
1       format(i2,i2,i2,i2)
        go to 3

2       write(6,*)' Input Error in a8_to_a9, a8_time = ',a8_time
        stop     

3       id = imon_a(imn)

        idays = id + idy

!       Decide whether to add a day for leap year.
        if(iyr .eq. (iyr / 4) * 4 )then
            if(imn .ge. 3)then
                idays = idays + 1
            endif
        endif

        im = 0

        write(a9,4)iyr,idays,ih,im
4       format(i2,i3,i2,i2)

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
cdoc  Routine to convert LAPS 'yyjjjhhmm' time to 'yymmddhh' time.
c     Corrected for Y2K.  P. Stamus, NOAA/FSL   Oct 1998
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

        function rsa13_to_a9(rsa13_time)

cdoc    Convert rsa13_time (yyyyjjjhhmmss) to a9_time (yyjjjhhmm)
!       Steve Albers 1998

        character*13 rsa13_time
        character*9 rsa13_to_a9

        rsa13_to_a9 = rsa13_time(3:11)

        istatus = 1

        return
        end

        function a9_to_rsa13(a9_time)

cdoc    Convert a9_time (yyjjjhhmm) to rsa13_time (yyyyjjjhhmmss)
!       Steve Albers 1998

        character*9 a9_time
        character*13 a9_to_rsa13
        character*2 cc

        if(a9_time(1:1) .eq. '9')then
            cc = '19'
        else
            cc = '20'
        endif

        a9_to_rsa13 = cc//a9_time//'00'

        istatus = 1

        return
        end
 
c---------------------------------------------------------------------------- 
      function a7_to_a9_time(a7_time)

c
cdoc  Routine takes FA (CWB) model intial time, as described below, and converts
cdoc  it to 9 character ascii format for laps. Routines below this 
cdoc  (make_fa_valtime, a9_to_a7_time, make_fa_ext, and fname13_to_FA_filename)
cdoc  take the FA model filename extension and convert it to a 4 character valid
cdoc  time, or vice versa.
c
c===========================================================================
cFA Model File description :
c===========================================================================
c   nf9093000.12m : Sample data file 1.
c                   9093000(YMMDDhh) => 1999/09/30/00z
c                   12m => 12hr fcst major-run(first cut)
c   nf9093012.12m : Sample data file 2.
c                   9093012(YMMDDhh) => 1999/09/30/12z
c===========================================================================

      implicit none

      integer istatus
      integer i4time_sys
c
c note: fname_in is 7 characters (ymmddhh)
c
      character*(*) a7_time
      character*8   a8_time
      character*9   a8_to_a9
      character*9   a7_to_a9_time
      character*9   a9_time

      call get_systime(i4time_sys,a9_time,istatus)

      a8_time=a9_time(1:1)//a7_time
      a9_time=a8_to_a9(a8_time)

      a7_to_a9_time = a9_time

      return
      end
c
c ---------------------------------------------------------------------------
c
      function a9_to_a7_time(a9_time)

      implicit none

      character*9  a9_time
      character*8  a9_to_a8
      character*8  a8_time
      character*7  a9_to_a7_time

      a8_time=a9_to_a8(a9_time)
      a9_to_a7_time=a8_time(2:8)

      return
      end
c
c ---------------------------------------------------------------------------
c
      function fname13_to_FA_filename(fname13,cmodel)

      implicit none

      character*13 fname13
      character*16 fname13_to_FA_filename
      character*10 fname10,a9_to_yr_a10_time
      character*9  a9_string
      character*4  a4_string
      character*(*) cmodel
      character*3  c3_FA_ext,c3_ext
      integer      lenc
      integer      istatus

      a9_string=fname13(1:9)
      a4_string=fname13(10:13)

c     fname7=a9_to_a7_time(a9_string)
      fname10=a9_to_yr_a10_time(a9_string,istatus)
      if(istatus.ne.1)then
         print*,'Error in conversion: a9_to_yr_a10_time'
         return
      endif
      c3_ext=c3_FA_ext(a4_string)
      call s_len(cmodel,lenc)
      if(cmodel(1:lenc).eq.'CWB_20FA_LAMBERT_NF')then
         fname13_to_FA_filename='nf'//fname10//'.'//c3_ext
      elseif(cmodel(1:lenc).eq.'CWB_20FA_LAMBERT_RE')then
         fname13_to_FA_filename='re'//fname10//'.'//c3_ext
      endif

      return
      end

c ---------------------------------------------------------------------------
c
      function c4_FA_valtime(c_FA_ext)

      implicit none

      character*3 c_FA_ext
      character*4 c4_FA_valtime

      c4_FA_valtime='0000'
      c4_FA_valtime(3:4)=c_FA_ext(1:2)

      return
      end

c ---------------------------------------------------------------------------
c
      function c3_FA_ext(a4_string)

      implicit none

      integer number
      character*4 a4_string
      character*3 c3_FA_ext

      read(a4_string,'(i4)')number
      if(number.lt.10)then
         c3_FA_ext='0'//a4_string(4:4)//'m'
      else
         c3_FA_ext=a4_string(3:4)//'m'
      endif

      return
      end

c ---------------------------------------------------------------------------
c
      function wrftolaps_c6_maprojname(map_proj_name)

      implicit none

      character*32  map_proj_name
      character*6   c6_maproj
      character*6   wrftolaps_c6_maprojname
 
      integer       len
      

      call s_len(map_proj_name,len)
      if(map_proj_name(1:len).eq.'polar')then
         c6_maproj='plrstr'
      elseif(map_proj_name(1:len).eq.'lambert')then
         c6_maproj='lambrt'
      elseif(map_proj_name(1:len).eq.'mercator')then
         c6_maproj='merctr'
       elseif(map_proj_name(1:len).eq.'rotlat')then
          c6_maproj='rotlat'
      else
         print*,'Error: Unknown map projection: ',' map_proj =
     1 ',map_proj_name(1:len),'.  Check the WRF namelist wrfsi.nl'
         return
      endif

      wrftolaps_c6_maprojname=c6_maproj
      return
      end
