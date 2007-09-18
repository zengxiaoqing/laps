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

      subroutine get_cloud_drift_afwa(i4time_sys,i4_ob_window
     1                                ,NX_L,NY_L
     1                                ,filename,istatus)

      character*(*) filename

!.............................................................................

      character*9 a9_timeObs
!     character*7 c7_skycover
      real lat_a(NX_L,NY_L)
      real lon_a(NX_L,NY_L)
      real topo_a(NX_L,NY_L)
      real latitude,longitude

!............................................................................

      lun_in = 21
      lun_out = 11

      open(lun_in,file=filename,status='old')

      r_mspkt = .518  

      i = 0

      do while (.true.)

          read(lun_in,101,err=890,end=999) !    NAME          UNITS & FACTOR
     1         I_G9CYC,                ! DB-CYCLE       
     1         I_G9GWCG,               ! GWC-REGION
     1         I_G9JULG,               ! JULIAN HOUR          HR since 673650000
     1         I_G9LATG,               ! LATITUDE             DEG * 100
     1         I_G9LONG,               ! LONGITUDE            DEG * 100 
     1         I_G9OPTC,               ! OPTIONS-CODE-GOES
     1         I_G9OTG,                ! HHMM
     1         I_G9PRGO,               ! PRESSURE             PA * 10 (HPA*1000)
     1         I_G9RTGO,               ! REPORT-TYPE-ID
     1         I_G9SATG,               ! SATELLITE OPERATOR
     1         I_G9SNG,                ! SATELLITE NUMBER
     1         I_G9TMGO,               ! TEMPERATURE
     1         I_G9WDGO,               ! WIND-DIRECTION       DEG
     1         I_G9WSGO,               ! WIND-SPEED           M/S * 10
     1         I_G9WUIG                ! WIND-SPEED UNITS     1=KT, 2=M/S
 101      format(14(i9,2x),i9)

          I_G9MIN = 0

          i = i + 1

          write(6,*)
          write(6,*)' cdw #',i

          latitude  =  float(I_G9LATG)/100.
          longitude = +float(I_G9LONG)/100.
          pres_pa   =  I_G9PRGO / 10.

          write(6,2)latitude,longitude,pres_pa
 2        format(' Lat, lon, pres_pa'/f8.3,f10.3,f8.0)  

          if(pres_pa .gt. 1e10)then
              write(6,*)' pres_pa is suspect - reject',pres_pa
              goto 900
          endif

          I_HR_OB  = I_G9OTG / 100
          I_MIN_OB = I_G9OTG - I_HR_OB*100
          I_HR_JUL = I_G9JULG - (I_G9JULG/24)*24

          if(I_HR_OB .ne. I_HR_JUL)then
              write(6,*)' Error, I_HR_OB .ne. I_HR_JUL',I_HR_OB,I_HR_JUL
              goto900
          endif

          call afwa_julhr_i4time(I_G9JULG,I_MIN_OB,i4time_ob)

          call make_fnam_lp(i4time_ob,a9_timeObs,istatus)
          if(istatus .ne. 1)goto900

          i4_resid = abs(i4time_ob - i4time_sys)
          if(i4_resid .gt. i4_ob_window)then ! outside time window
              write(6,*)' time - reject '
     1           ,a9_timeObs,i4_resid,i4_ob_window
              goto 900        
          endif

          windDIR = I_G9WDGO

!         Calculate Wind Speed in M/S
          if(I_G9WUIG .eq. 1)then     ! KT
!             write(6,*)' Convert Wind speed from KT to M/S ',I_G9WUIG       
              windSpeed = (float(I_G9WSGO) / 10.) * r_mspkt

          elseif(I_G9WUIG .eq. 2)then ! M/S
              windSpeed = float(I_G9WSGO) / 10.

          else ! invalid units identifier
              write(6,*)' warning: wind speed units are undefined '
     1                 ,I_G9WUIG,I_G9WSGO     
              goto900

          endif

          if(abs(windSpeed) .gt. 250. .or. windDir .gt. 360.)then
              write(6,*)' wind is suspect - reject',windDir,windSpeed
              goto900

          else ! write out valid wind
              write(6      ,21)latitude,longitude,pres_pa
     1                        ,winddir,windspeed,a9_timeObs       
              write(lun_out,21)latitude,longitude,pres_pa
     1                        ,winddir,windspeed,a9_timeObs       
 21           format(f8.3,f10.3,f8.0,f6.0,f6.1,2x,a9)

          endif

          go to 900

 890      write(6,*)' Warning: read error'

 900  enddo ! read line of AFWA file

!............................................................................

 999  write(6,*)' End of AFWA file detected'

      close(lun_in)
      istatus = 1
      return
      end
