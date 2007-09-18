           
      subroutine get_satsnd_afwa(i4time_sys,i4_satsnd_window
     1                          ,NX_L,NY_L
     1                          ,lun_in,filename,lun_out,istatus)

!     Steve Albers FSL    May 1999

      character*(*) filename

!.............................................................................

!     character*6 C6_A1ACID
      character*9 a9_timeObs,a9_recptTime,a9time_ob 
      character*1000 c_line
      character*9 c_read
      character*5 c5_staid
      character*8 c8_obstype

      integer MAX_LEVELS

      parameter (MAX_LEVELS = 50)
      real RHEIGHT(MAX_LEVELS)
      real RH(MAX_LEVELS)
      real PRESSURE(MAX_LEVELS)
      real TEMP_K(MAX_LEVELS)

      real lat_a(NX_L,NY_L)
      real lon_a(NX_L,NY_L)
      real topo_a(NX_L,NY_L)

!............................................................................

      open(lun_in,file=filename,status='old')

      call get_r_missing_data(r_missing_data,istatus)
      if(istatus .ne. 1)then
          write(6,*)' Error in get_r_missing_data'
          return
      endif

      call get_domain_perimeter(NX_L,NY_L,'nest7grid',lat_a,lon_a, 
     1            topo_a,1.0,rnorth,south,east,west,istatus)
      if(istatus .ne. 1)then
          write(6,*)' Error in get_laps_perimeter'
          return
      endif
  
      i = 0

      do while (.true.)

          read(lun_in,51,err=890,end=999)c_line
 51       format(a)

          read(c_line,101,err=890)     !    NAME             UNITS & FACTOR
     1         I_A1CYCC,                       
     1         I_A1GWC,
     1         I_A1JUL,                ! JULIAN HOUR          HR since 673650000
     1         I_A1LAT,                ! LATITUDE             DEG * 100
     1         I_A1LON,                ! LONGITUDE            DEG * -100 
     1         I_A1NLVL

!    1         I_A1TYPE,
!    1         I_A1MIN,                ! TIME-REPORT-MINUTES
!    1         I_A1KIND,
!    1         I_A1PLA,
!    1         I_A1DVAL,
!    1         I_A1HOLM,
!    1         I_A1FLTP,               ! TEMPERATURE          KELVINS * 10
!    1         I_A1WD,                 ! WIND-DIRECTION       DEG
!    1         I_A1WFLS,               ! WIND-SPEED           M/S * 10
!    1         C6_A1ACID 

 101      format(6(i9,2x))

!         nlvls = I_A1NLVL
          nlvls = 16          ! Hardwired as per AFWA documentation

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
              pressure(lvl) = float(I_PRESSURE) / 10.
          enddo ! lvl

          iblk = 4
          do lvl = 1,nlvls
              istart = 66 + (iblk-1)*nlvls*11 + (lvl-1)*11 + 1
              iend = istart+8
              c_read = c_line(istart:iend) 
              read(c_read,102,err=890)I_TEMP
              temp_k(lvl) = float(I_TEMP) / 10.
          enddo ! lvl

!         Read hours & minutes
          iblk = 5
          lvl = 2
          istart = 66 + (iblk-1)*nlvls*11 + (lvl-1)*11 + 1
          iend = istart+8
          c_read = c_line(istart:iend) 
          read(c_read,112,err=890)I_HR,I_MIN
 112      format(5x,2i2)

!         Read satellite ID
          iblk = 5
          lvl = 3
          istart = 66 + (iblk-1)*nlvls*11 + (lvl-1)*11 + 1
          iend = istart+8
          c_read = c_line(istart:iend) 
          c5_staid = c_read(1:2)//c_read(4:6)

          i = i + 1

          write(6,*)
          write(6,*)' satsnd #',i

          stalat =  float(I_A1LAT)/100.
          stalon = +float(I_A1LON)/100.

          write(6,2)stalat,stalon
 2        format(' Lat, lon '/f8.3,f10.3)  

          if(stalat .le. rnorth .and. stalat .ge. south .and.
     1       stalon .ge. west   .and. stalon .le. east        )then       
              continue
          else ! Outside lat/lon perimeter - reject
              write(6,*)' lat/lon - reject'       
!    1                 ,stalat,stalon
              goto 900
          endif

          I_HR_JUL = I_A1JUL - (I_A1JUL/24)*24

          write(6,*)' I_HR, I_HR_JUL, I_MIN ', I_HR, I_HR_JUL, I_MIN       

          if(I_HR .ne. I_HR_JUL)then
              write(6,*)' WARNING: I_HR discrepancy - reject '
     1                 ,I_HR,I_HR_JUL
              goto900
          endif

          call afwa_julhr_i4time(I_A1JUL,I_MIN,i4time_ob)

          call make_fnam_lp(i4time_ob,a9time_ob,istatus)
          if(istatus .ne. 1)goto900

          a9_recptTime = '         '

          i4_resid = abs(i4time_ob - i4time_sys)
          if(i4_resid .gt. i4_satsnd_window)then ! outside time window
              write(6,*)' time - reject '
     1           ,a9time_ob,i4_resid,i4_satsnd_window
              goto 900        
          endif

          write(6,1)a9time_ob
 1        format(' Time:'/1x,a9) 

          staelev = 0.

          c8_obstype = 'SATSND'

          iwmostanum = 0

          write(6,511,err=990)
     1             iwmostanum,nlvls,stalat
     1            ,stalon,staelev,c5_staid,a9time_ob,c8_obstype
          write(lun_out,511,err=990)
     1             iwmostanum,nlvls,stalat
     1            ,stalon,staelev,c5_staid,a9time_ob,c8_obstype

  511     format(i12,i12,f11.4,f15.4,f15.0,1x,a5,3x,a9,1x,a8)

          dir = r_missing_data
          spd = r_missing_data

          do lvl = 1,nlvls
              if(abs(rheight(lvl)) .gt. 100000.)then
                  rheight(lvl) = r_missing_data
              endif

              if(abs(pressure(lvl)) .gt. 100000.)then
                  pressure(lvl) = r_missing_data
              endif

              if(temp_k(lvl) .ge. 150. .and. temp_k(lvl) .le. 400.)then       
                  temp_c = temp_k(lvl) - 273.15
              else
                  temp_c = r_missing_data
              endif

!             Convert rh(lvl) to dewpoint
              if(rh(lvl) .gt. 0. .and. rh(lvl) .le. 100.
     1                           .and. temp_c .ne. r_missing_data)then       
                  dewpoint_c = dwpt(temp_c,rh(lvl))
              else
                  dewpoint_c = r_missing_data
              endif

              if(i .le. 100)write(6,*)rheight(lvl),pressure(lvl)
     1              ,temp_c
     1              ,dewpoint_c
     1              ,lvl

              write(lun_out,*)rheight(lvl),pressure(lvl)
     1              ,temp_c
     1              ,dewpoint_c
     1              ,dir,spd

          enddo ! lvl

          go to 900

 890      write(6,*)' Warning (get_satsnd_afwa): read/write error'

 900  enddo ! read line of AFWA file

!............................................................................

 990  write(6,*)' ERROR in ingest_satsnd_af'
      istatus=0
      return
  
 999  write(6,*)' End of AFWA file detected'
      istatus = 1
      return
      end

