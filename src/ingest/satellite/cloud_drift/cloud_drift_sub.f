cdis   
cdis    Open Source License/Disclaimer, Forecast Systems Laboratory
cdis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
cdis    
cdis    This software is distributed under the Open Source Definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    In particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - Redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - Redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - All modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - If significant modifications or enhancements are made to this
cdis    software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
cdis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
cdis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
cdis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
cdis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
cdis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
cdis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
cdis   
cdis
cdis
cdis   
cdis

      subroutine get_cloud_drift_data(i4time_sys,i4_window
     1                                    ,NX_L,NY_L
     1                                    ,filename,istatus)

!     Steve Albers     Jan-1998

!.............................................................................

      character*170 filename
      character*10 a10_time
      character*9 a9_timeObs,a10_to_a9 
      character*6 a6_time
      character*5 a5_time
      character*4 a4_time
      character*2 c2_sat_type

      real lat_a(NX_L,NY_L)
      real lon_a(NX_L,NY_L)
      real topo_a(NX_L,NY_L)

      logical l_new_fmt, l_parse
      
      l_new_fmt = .true. ! .false.
      if(l_parse(filename,'goes11'))l_new_fmt = .true.
      if(l_parse(filename,'goes12'))l_new_fmt = .true.

      call get_domain_perimeter(NX_L,NY_L,'nest7grid',lat_a,lon_a, 
     1            topo_a,1.0,rnorth,south,east,west,istatus)
      if(istatus .ne. 1)then
          write(6,*)' Error in get_laps_perimeter'
          return
      endif

      open(1,file=filename,status='old',err=990)
      i = 0

!     Read the file...
      if(l_new_fmt)then
          write(6,*)' Reading ASCII file in new format ',filename
          nheader = 2
      else
          write(6,*)' Reading ASCII file in old format ',filename
          nheader = 1
      endif

      do iheader = 1,nheader
          read(1,*)               ! Skip header lines
      enddo ! iheader

 40   if(.not. l_new_fmt)then
          read(1,11,err=50,end=990)c2_sat_type,a6_time,a4_time
     1                            ,rlat,rlon,ipres,spd,idir
 11       format(a2,15x,a6,2x,a4,f8.2,f9.2,i6,f7.1,i5)
      else ! new format
          read(1,12,err=50,end=990)c2_sat_type,a5_time,a4_time
     1                            ,rlat,rlon,ipres,spd,idir
 12       format(a2,17x,a5,4x,a4,2x,f10.0,f10.0,i10,f10.0,i10)
      endif

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

      if(c2_sat_type .ne. 'IR' .and. c2_sat_type .ne. 'VI'
     1                         .and. c2_sat_type .ne. 'WV')then       
          write(6,*)' Bad Sat Type ',c2_sat_type
          goto 900
      endif

      if(l_new_fmt)then
          a9_timeObs = a5_time//a4_time
      else
          a10_time = a6_time//a4_time
          a9_timeObs = a10_to_a9(a10_time,istatus)
          if(istatus .ne. 1)then
              write(6,*)' Bad observation time - reject ',a10_time       
              goto 900
          endif
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
