
      subroutine get_satsnd_afwa(i4time_sys,i4_satsnd_window
     1                                     ,NX_L,NY_L
     1                                     ,filename,istatus)

!     Steve Albers FSL    May 1999

      character*(*) filename

!.............................................................................

      character*6 C6_A1ACID
      character*9 a9_timeObs,a9_recptTime,a9time_ob 
      character*1000 c_line
      character*6 c_read
      character*5 c5_staid
      character*8 c8_obstype

      integer MAX_LEVELS

      parameter (MAX_LEVELS = 50)
      real*4 RHEIGHT(MAX_LEVELS)
      real*4 RH(MAX_LEVELS)
      real*4 PRESSURE(MAX_LEVELS)
      real*4 TEMP(MAX_LEVELS)

      real*4 lat_a(NX_L,NY_L)
      real*4 lon_a(NX_L,NY_L)
      real*4 topo_a(NX_L,NY_L)

!............................................................................

      open(11,file=filename,status='old')

      call get_domain_perimeter(NX_L,NY_L,'nest7grid',lat_a,lon_a, 
     1            topo_a,1.0,rnorth,south,east,west,istatus)
      if(istatus .ne. 1)then
          write(6,*)' Error in get_laps_perimeter'
          return
      endif
  
      i = 0

      do while (.true.)

          read(11,*,err=890,end=999)c_line

          read(c_line,101,err=890,end=999) !    NAME             UNITS & FACTOR
     1         I_A1CYCC,                       
     1         I_A1GWC,
     1         I_A1JUL,                ! JULIAN HOUR          HR since 673650000
     1         I_A1LAT,                ! LATITUDE             DEG * 100
     1         I_A1LON,                ! LONGITUDE            DEG * -100 
     1         I_A1NLVL

!    1         I_A1TYPE,
!    1         I_A1MIN,                ! TIME-REPORT-MINUTES
!    1         I_A1KIND,
!    1         I_A1ALT,                ! FLIGHT-ALT (true)    Meters MSL 
!    1         I_A1PLA,
!    1         I_A1DVAL,
!    1         I_A1HOLM,
!    1         I_A1FLTP,               ! TEMPERATURE          KELVINS * 10
!    1         I_A1WD,                 ! WIND-DIRECTION       DEG
!    1         I_A1WFLS,               ! WIND-SPEED           M/S * 10
!    1         C6_A1ACID 

 101      format(6(i9,2x))

          nlvls = I_A1NLVL

          iblk = 1
          do lvl = 1,nlvls
              istart = 66 + (iblk-1)*nlvls*11 + (lvl-1)*11 + 1
              iend = istart+8
              c_read = c_line(istart:iend) 
              read(c_read,102,err=890)I_HEIGHT
 102          format(i9)
              rheight(lvl) = I_HEIGHT
          enddo ! lvl

          iblk = 2
          do lvl = 1,nlvls
              istart = 66 + (iblk-1)*nlvls*11 + (lvl-1)*11 + 1
              iend = istart+8
              c_read = c_line(istart:iend) 
              read(c_read,102,err=890)I_RH
              rh(lvl) = I_RH
          enddo ! lvl

          iblk = 3
          do lvl = 1,nlvls
              istart = 66 + (iblk-1)*nlvls*11 + (lvl-1)*11 + 1
              iend = istart+8
              c_read = c_line(istart:iend) 
              read(c_read,102,err=890)I_PRESSURE
              pressure(lvl) = I_PRESSURE
          enddo ! lvl

          iblk = 4
          do lvl = 1,nlvls
              istart = 66 + (iblk-1)*nlvls*11 + (lvl-1)*11 + 1
              iend = istart+8
              c_read = c_line(istart:iend) 
              read(c_read,102,err=890)I_TEMP
              temp(lvl) = float(I_TEMP)/10.
          enddo ! lvl

          i = i + 1

          write(6,*)
          write(6,*)' satsnd #',i

          stalat =  float(I_A1LAT)/100.
          stalon = -float(I_A1LON)/100.

          write(6,*)' location = ',stalat,stalon

          if(stalat .le. rnorth .and. stalat .ge. south .and.
     1       stalon .ge. west   .and. stalon .le. east        )then       
              continue
          else ! Outside lat/lon perimeter - reject
              write(6,*)' lat/lon - reject'       
!    1                 ,stalat,stalon
              goto 900
          endif

          call afwa_jul_i4time(I_A1JUL,I_A1MIN,i4time_ob)

          call make_fnam_lp(i4time_ob,a9_time_ob,istatus)
          if(istatus .ne. 1)goto900

          a9_recptTime = '         '

!         call cv_asc_i4time(a9_timeObs,i4time_ob)

          i4_resid = abs(i4time_ob - i4time_sys)
          if(i4_resid .gt. i4_satsnd_window)then ! outside time window
              write(6,*)' time - reject '
     1           ,a9_time_ob,i4_resid,i4_satsnd_window
              goto 900        
          endif

          write(6,1)a9_time_ob,a9_recptTime 
          write(11,1)a9_time_ob,a9_recptTime 
 1        format(' Time - prp/rcvd:'/1x,a9,2x,a9) 

          write(6,2)stalat,stalon
          write(11,2)stalat,stalon
 2        format(' Lat, lon '/f8.3,f10.3,f8.0)  

          staelev = 0.

          c8_obstype = 'SATSND'

          write(6,511,err=990)
     1             iwmostanum,nlvls,stalat
     1            ,stalon,staelev,c5_staid,a9time_ob,c8_obstype
          write(11,511,err=990)
     1             iwmostanum,nlvls,stalat
     1            ,stalon,staelev,c5_staid,a9time_ob,c8_obstype

  511     format(i12,i12,f11.4,f15.4,f15.0,1x,a5,3x,a9,1x,a8)

          dewpoint = r_missing_data
          dir = r_missing_data
          spd = r_missing_data

          do lvl = 1,nlvls
              if(abs(temp(lvl)) .gt. 400.)then
                  temp(lvl) = r_missing_data
              endif

              write(6,*)rheight(lvl),pressure(lvl)
     1              ,temp(lvl)
     1              ,ilvl

              write(11,*)rheight(lvl),pressure(lvl)
     1              ,temp(lvl)
     1              ,dewpoint
     1              ,dir,spd

          enddo ! lvl

          go to 900

 890      write(6,*)' Warning: read/write error'

 900  enddo ! read line of AFWA file

!............................................................................

 990  write(6,*)' ERROR in ingest_satsnd_af'
      istatus=0
      return
  
 999  write(6,*)' End of AFWA file detected'
      istatus = 1
      return
      end


      subroutine afwa_jul_i4time(I_A1JUL,I_A1MIN,i4time)

!     I_A1JUL is number of hours since Dec 31, 1967 at 00z
!     This is converted to i4time, number of sec since Jan 1, 1960 at 00z
      i4time_hr  = I_A1JUL * 3600 + (8*365 - 1 - 2) * 86400
      i4time_min = I_A1MIN*60
      i4time     = i4time_hr + i4time_min

      return
      end

