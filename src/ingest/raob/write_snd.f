           
      subroutine write_snd(lun_out                         ! I
     1                    ,maxsnd,maxlvl,nsnd              ! I
     1                    ,iwmostanum                      ! I
     1                    ,stalat,stalon,staelev           ! I
     1                    ,c5_staid,a9time_ob,c8_obstype   ! I
     1                    ,nlvl                            ! I
     1                    ,height_m                        ! I
     1                    ,pressure_mb                     ! I
     1                    ,temp_c                          ! I
     1                    ,dewpoint_c                      ! I
     1                    ,dir_deg                         ! I
     1                    ,spd_mps                         ! I
     1                    ,istatus)                        ! O

!     Steve Albers FSL    2001

!     Write routine for 'snd' file

!     For missing data values, 'r_missing_data' should be passed in 

!.............................................................................

      integer iwmostanum(maxsnd),nlvl(maxsnd)
      real stalat(maxsnd,maxlvl),stalon(maxsnd,maxlvl),staelev(maxsnd)       
      character c5_staid(maxsnd)*5,a9time_ob(maxsnd,maxlvl)*9
     1         ,c8_obstype(maxsnd)*8,c_line*200

      character*9 a9_time

      character*5 c5_sta

      real height_m(maxsnd,maxlvl)
      real pressure_mb(maxsnd,maxlvl)
      real temp_c(maxsnd,maxlvl)
      real dewpoint_c(maxsnd,maxlvl)
      real dir_deg(maxsnd,maxlvl)
      real spd_mps(maxsnd,maxlvl)

!............................................................................

!     Get RAOB Time Window
      call get_windob_time_window('RAOB',i4_wind_ob,istatus)
      if(istatus .ne. 1)goto 990

      call get_tempob_time_window('RAOB',i4_temp_ob,istatus)
      if(istatus .ne. 1)goto 990

      i4_raob_window = max(i4_wind_ob,i4_temp_ob)

!     Get systime
      call GETENV('LAPS_A9TIME',a9_time)
      call s_len(a9_time,ilen)

      if(ilen .eq. 9)then
!       write(6,*)' systime (from env) = ',a9_time
        call i4time_fname_lp(a9_time,i4time_sys,istatus)
      else
        call get_systime(i4time_sys,a9_time,istatus)
        if(istatus .ne. 1)go to 990
!       write(6,*)' systime = ',a9_time
      endif

      do isnd = 1,nsnd

!       Reject observation times outside time window        
        call i4time_fname_lp(a9time_ob(isnd,1),i4time_raob,status)
        if(abs(i4time_raob - i4time_sys) .gt. i4_raob_window)then
            write(6,*)a9time_ob(isnd,1),
     1                ' is outside time window with sounding #',isnd
            goto 900
        endif

        call s_len(c5_staid(isnd),len_sta)
        if(len_sta .gt. 0)then
            c5_sta = c5_staid(isnd)
        else
            if(iwmostanum(isnd) .gt. 0 .and. 
     1         iwmostanum(isnd) .le. 99999)then
                write(6,*)' Filling in blank staid with WMOID'
                write(c5_sta,101)iwmostanum(isnd)
 101            format(i5)
            endif
        endif

!       Write Sounding Header

        write(6,511,err=990)
     1             iwmostanum(isnd),nlvl(isnd)
     1            ,stalat(isnd,1),stalon(isnd,1),staelev(isnd)
     1            ,c5_sta,a9time_ob(isnd,1),c8_obstype(isnd)

        write(lun_out,511,err=990)
     1             iwmostanum(isnd),nlvl(isnd)
     1            ,stalat(isnd,1),stalon(isnd,1),staelev(isnd)
     1            ,c5_sta,a9time_ob(isnd,1),c8_obstype(isnd)

  511   format(i12,i12,f11.4,f15.4,f15.0,1x,a5,3x,a9,1x,a8)

        do lvl = 1,nlvl(isnd)

!         Write Sounding Level (the character array helps keep everything
!                               in one line when using free format)

          write(c_line,*)height_m(isnd,lvl)," "
     1              ,pressure_mb(isnd,lvl)," "
     1              ,temp_c(isnd,lvl)," "
     1              ,dewpoint_c(isnd,lvl)," "
     1              ,dir_deg(isnd,lvl)," "
     1              ,spd_mps(isnd,lvl)," "
     1              ,a9time_ob(isnd,lvl)," "
     1              ,stalat(isnd,lvl)," "
     1              ,stalon(isnd,lvl)
          call s_len2(c_line,len_line)
          write(lun_out,521)c_line(1:len_line)
  521     format(a)

          if(isnd .le. 100)then
              write(6,521)c_line(1:len_line)
          endif

        enddo ! lvl

 900    continue

      enddo ! isnd

      go to 999

 990  write(6,*)' ERROR in write_snd'
      istatus=0
      return

 999  istatus = 1
      return
      end

