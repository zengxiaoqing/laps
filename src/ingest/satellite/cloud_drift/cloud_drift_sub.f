
      subroutine get_cloud_drift_data(i4time_sys,i4_window
     1                                    ,NX_L,NY_L
     1                                    ,filename,istatus)

!     Steve Albers     Jan-1998

!.............................................................................

      character*170 filename
      character*10 a10_time
      character*9 a9_timeObs,a10_to_a9 
      character*6 a6_time
      character*4 a4_time
      character*2 c2_sat_type

      real*4 lat_a(NX_L,NY_L)
      real*4 lon_a(NX_L,NY_L)
      real*4 topo_a(NX_L,NY_L)

      call get_domain_perimeter(NX_L,NY_L,'nest7grid',lat_a,lon_a, 
     1            topo_a,1.0,rnorth,south,east,west,istatus)
      if(istatus .ne. 1)then
          write(6,*)' Error in get_laps_perimeter'
          return
      endif

      open(1,file=filename,status='old',err=990)
      i = 0

!     Read the file...
      read(1,*)               ! Skip header line
 40   read(1,11,err=50,end=990)c2_sat_type,a6_time,a4_time
     1                        ,rlat,rlon,ipres,spd,idir
 11   format(a2,13x,a6,2x,a4,f8.2,f9.2,i6,f7.1,i5)
 50   continue

      i = i+1
      pres_pa = float(ipres) * 100. ! Pascals
      dir = float(idir)
      rlon = -rlon

      if(rlat .le. rnorth .and. rlat .ge. south .and.
     1   rlon .ge. west   .and. rlon .le. east        )then       
          write(6,*)
          write(6,*)' cloud_drift #',i

      else ! Outside lat/lon perimeter - reject
!         write(6,*)' lat/lon - reject'       
          goto 900
      endif

      if(c2_sat_type .ne. 'IR' .and. c2_sat_type .ne. 'VI')then       
          write(6,*)' Bad Sat Type ',c2_sat_type
          goto 900
      endif

      a10_time = a6_time//a4_time
      a9_timeObs = a10_to_a9(a10_time,istatus)
      if(istatus .ne. 1)then
          write(6,*)' Bad observation time - reject ',a10_time       
          goto 900
      endif

      call cv_asc_i4time(a9_timeObs,i4time_ob)
      i4_resid = abs(i4time_ob - i4time_sys)
      if(i4_resid .gt. i4_window)then ! outside time window
          write(6,*)' time - reject ',a9_timeObs
     1                               ,i4_resid,i4_window
          goto 900        
      endif

      write(6 ,21)rlat,rlon,pres_pa,dir,spd,a9_timeObs
      write(11,21)rlat,rlon,pres_pa,dir,spd,a9_timeObs
 21   format(f8.3,f10.3,f8.0,f6.0,f6.1,2x,a9)

 900  continue
      
      go to 40 ! Loop back to read next line

 990  close(1)

      return
      END
