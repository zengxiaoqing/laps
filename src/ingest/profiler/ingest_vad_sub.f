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

      subroutine get_vad_data(i4time_sys,ilaps_cycle_time
     1                                    ,NX_L,NY_L
     1                                    ,filename,istatus)

      character*170 filename

!.............................................................................

      include 'netcdf.inc'
      integer maxLevels, recNum,nf_fid, nf_vid, nf_status
C
C  Open netcdf File for reading
C
      nf_status = NF_OPEN(filename,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'NF_OPEN',filename
      endif
C
C  Fill all dimension values
C
C
C Get size of maxLevels
C
      nf_status = NF_INQ_DIMID(nf_fid,'maxLevels',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim maxLevels'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,maxLevels)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim maxLevels'
      endif
C
C Get size of recNum
C
      nf_status = NF_INQ_DIMID(nf_fid,'recNum',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim recNum'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,recNum)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim recNum'
      endif
      call main_sub(nf_fid, maxLevels, recNum,
!.............................................................................
     1              i4time_sys,ilaps_cycle_time,NX_L,NY_L,istatus)
!.............................................................................

      return
      end
C
C
      subroutine main_sub(nf_fid, maxLevels, recNum,
!.............................................................................
     1              i4time_sys,ilaps_cycle_time,NX_L,NY_L,istatus)

      real lat_a(NX_L,NY_L)
      real lon_a(NX_L,NY_L)
      real topo_a(NX_L,NY_L)

      character*9 a9time_ob
      character*4 c4_staname

!.............................................................................

      include 'netcdf.inc'
      integer maxLevels, recNum,nf_fid
      integer numLevels(recNum), windDir( maxLevels, recNum)
      real obHeight( maxLevels, recNum), staElev(recNum),
     +     staLat(recNum), staLon(recNum), windSpeed( maxLevels,
     +     recNum)
      double precision timeObs(recNum)
      character*4 staName(recNum)
      character*350 rawMsg(recNum)
      character rmsError( maxLevels, recNum)

!............................................................................

      call read_netcdf(nf_fid, maxLevels, recNum, numLevels, windDir, 
     +     obHeight, staElev, staLat, staLon, windSpeed, timeObs, 
     +     rawMsg, rmsError, staName)
C
C The netcdf variables are filled - your code goes here
C
!............................................................................

      call get_latlon_perimeter(NX_L,NY_L,1.0
     1                           ,lat_a,lon_a,topo_a
     1                           ,rnorth,south,east,west,istatus)
      if(istatus .ne. 1)then
          write(6,*)' Error reading LAPS perimeter'
          return
      endif

      num_vad = recNum
      write(6,*)' # of vad = ',recNum

      do ista = 1,num_vad

          write(6,*)
          write(6,*)' vad #',ista

!         write(6,*)' location = '
!    1             ,stalat(ista),stalon(ista),staelev(ista)

          if(stalat(ista) .le. rnorth .and. stalat(ista) .ge. south 
     1                                .and.
     1       stalon(ista) .ge. west   .and. stalon(ista) .le. east      
     1                                                             )then
              continue
          else ! Outside lat/lon perimeter - reject
              write(6,*)' lat/lon - reject'       
!    1                 ,stalat(ista),stalon(ista)
              goto 900
          endif


          if(abs(timeObs(ista))      .lt. 3d9)then
              call c_time2fname(nint(timeObs(ista)),a9time_ob)
          else
              write(6,*)' Bad observation time - reject',timeObs(ista)
              goto 900
          endif

          call cv_asc_i4time(a9time_ob,i4time_ob)
          i4_resid = abs(i4time_ob - i4time_sys)
          if(i4_resid .gt. (ilaps_cycle_time / 2) )then ! outside time window
              write(6,*)' time - reject '
     1           ,a9time_ob,i4_resid,ilaps_cycle_time / 2
              goto 900        
          endif

          write(6,1)a9time_ob
 1        format(' Ob Time:',1x,a9)

          id_num = 0

          c4_staname = 'K'//staname(ista)(1:3)

          if(abs(staelev(ista)) .gt. 999999.)then
              staelev_out = 999999.
          else
              staelev_out = staelev(ista)
          endif

          write(6,401)id_num,numlevels(ista),
     1                stalat(ista),stalon(ista),staelev_out,
     1                c4_staname,a9time_ob,'VAD     '
          write(1,401)id_num,numlevels(ista),
     1                stalat(ista),stalon(ista),staelev_out,
     1                c4_staname,a9time_ob,'VAD     '
401       format(i12,i12,f11.3,f15.3,f15.0,5x,a4,5x,a9,1x,a8)

          do ilvl = 1,numlevels(ista)
!             Set the rms error for the level
              rms = 10.
              if(rmserror(ilvl,ista)(1:1) .eq. 'A')rms = 1.03
              if(rmserror(ilvl,ista)(1:1) .eq. 'B')rms = 2.06
              if(rmserror(ilvl,ista)(1:1) .eq. 'C')rms = 3.09
              if(rmserror(ilvl,ista)(1:1) .eq. 'D')rms = 4.12
              if(rmserror(ilvl,ista)(1:1) .eq. 'E')rms = 5.14
              if(rmserror(ilvl,ista)(1:1) .eq. 'F')rms = 6.17
              if(rmserror(ilvl,ista)(1:1) .eq. 'G')rms = 7.20

              write(1,301,err=303)obheight(ilvl,ista),
     1                            float(winddir(ilvl,ista)), 
     1                            windspeed(ilvl,ista),rms       
              write(6,301,err=303)obheight(ilvl,ista),
     1                            float(winddir(ilvl,ista)), 
     1                            windspeed(ilvl,ista),rms       
301           format(1x,f6.0,f6.0,2f6.1)
303           continue
          enddo ! ilvl

 900  enddo ! i

!............................................................................

      return
      end
C
C  Subroutine to read the file "VAD Wind data : selected by receipt time : time range from 883602000 to 883603800" 
C
      subroutine read_netcdf(nf_fid, maxLevels, recNum, numLevels, 
     +     windDir, obHeight, staElev, staLat, staLon, windSpeed, 
     +     timeObs, rawMsg, rmsError, staName)
C
      include 'netcdf.inc'
      integer maxLevels, recNum,nf_fid, nf_vid, nf_status
      integer numLevels(recNum), windDir( maxLevels, recNum)
      real obHeight( maxLevels, recNum), staElev(recNum),
     +     staLat(recNum), staLon(recNum), windSpeed( maxLevels,
     +     recNum)
      double precision timeObs(recNum)
      character*4 staName(recNum)
      character*350 rawMsg(recNum)
      character rmsError( maxLevels, recNum)


C   Variables of type REAL
C
C     Variable        NETCDF Long Name
C      obHeight     "Observation height above MSL"
C
        nf_status = NF_INQ_VARID(nf_fid,'obHeight',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var obHeight'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,obHeight)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var obHeight'
      endif
C
C     Variable        NETCDF Long Name
C      staElev      "Elevation of station above MSL"
C
        nf_status = NF_INQ_VARID(nf_fid,'staElev',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var staElev'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,staElev)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var staElev'
      endif
C
C     Variable        NETCDF Long Name
C      staLat       "Station latitude"
C
        nf_status = NF_INQ_VARID(nf_fid,'staLat',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var staLat'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,staLat)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var staLat'
      endif
C
C     Variable        NETCDF Long Name
C      staLon       "Station longitude"
C
        nf_status = NF_INQ_VARID(nf_fid,'staLon',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var staLon'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,staLon)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var staLon'
      endif
C
C     Variable        NETCDF Long Name
C      windSpeed    "Wind speed"
C
        nf_status = NF_INQ_VARID(nf_fid,'windSpeed',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windSpeed'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,windSpeed)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windSpeed'
      endif

C   Variables of type INT
C
C
C     Variable        NETCDF Long Name
C      numLevels    "Number of observations for station"
C
        nf_status = NF_INQ_VARID(nf_fid,'numLevels',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var numLevels'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,numLevels)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var numLevels'
      endif
C
C     Variable        NETCDF Long Name
C      windDir      "Wind direction"
C
        nf_status = NF_INQ_VARID(nf_fid,'windDir',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDir'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,windDir)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDir'
      endif

C   Variables of type DOUBLE
C
C
C     Variable        NETCDF Long Name
C      timeObs      "Time of observation"
C
        nf_status = NF_INQ_VARID(nf_fid,'timeObs',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var timeObs'
      endif
        nf_status = NF_GET_VAR_DOUBLE(nf_fid,nf_vid,timeObs)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var timeObs'
      endif


C   Variables of type CHAR
C
C
C     Variable        NETCDF Long Name
C      rawMsg       "Receipt format VAD Winds message"
C
        nf_status = NF_INQ_VARID(nf_fid,'rawMsg',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var rawMsg'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,rawMsg)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var rawMsg'
      endif
C
C     Variable        NETCDF Long Name
C      rmsError     "RMS vector wind error"
C
        nf_status = NF_INQ_VARID(nf_fid,'rmsError',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var rmsError'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,rmsError)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var rmsError'
      endif
C
C     Variable        NETCDF Long Name
C      staName      "Three letter WMO site identifier"
C
        nf_status = NF_INQ_VARID(nf_fid,'staName',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var staName'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,staName)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var staName'
      endif

      nf_status = nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
      endif

      return
      end


