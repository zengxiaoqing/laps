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
C
C  Subroutine to read the file 
C
      subroutine read_gps(nf_fid, recNum, 
     +     pressure, staElev, staLat, staLon, temperature, 
     +     relativeHumidity, timeObs, staNam)
C
      include 'netcdf.inc'
      integer recNum,nf_fid, nf_vid, nf_status

      real dryDelay(recNum), formalError(recNum), pressure(recNum),
     +     staElev(recNum), staLat(recNum), staLon(recNum),
     +     temperature(recNum), relativeHumidity(recNum), 
     +     totalDelay(recNum),
     +     waterVapor(recNum), wetDelay(recNum)
      double precision timeObs(recNum)
      character*5 staNam(recNum)
      character*80 staLongNam(recNum)


C   Variables of type REAL
C
C     Variable        NETCDF Long Name
C      dryDelay     "Dry component GPS signal delay"
C
        nf_status = NF_INQ_VARID(nf_fid,'dryDelay',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dryDelay'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,dryDelay)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dryDelay'
      endif
C
C     Variable        NETCDF Long Name
C      formalError  "Formal Error"
C
        nf_status = NF_INQ_VARID(nf_fid,'formalError',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var formalError'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,formalError)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var formalError'
      endif
C
C     Variable        NETCDF Long Name
C      pressure     "Pressure used for PWV calculation"
C
        nf_status = NF_INQ_VARID(nf_fid,'pressure',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var pressure'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,pressure)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var pressure'
      endif
C
C     Variable        NETCDF Long Name
C      staElev      "Station Elevation (above MSL)"
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
C      staLat       "Station Latitude"
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
C      staLon       "Station Longitude"
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
C      temperature  "Temperature used for calculation"
C
        nf_status = NF_INQ_VARID(nf_fid,'temperature',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var temperature'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,temperature)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var temperature'
      endif
C
C     Variable        NETCDF Long Name
C      RH  "PerCent"
C
        nf_status = NF_INQ_VARID(nf_fid,'relativeHumidity',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var relativeHumidity'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,relativeHumidity)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var relativeHumidity'
      endif
C
C     Variable        NETCDF Long Name
C      totalDelay   "Total GPS signal delay"
C
        nf_status = NF_INQ_VARID(nf_fid,'totalDelay',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var totalDelay'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,totalDelay)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var totalDelay'
      endif
C
C     Variable        NETCDF Long Name
C      waterVapor   "Water Vapor"
C
        nf_status = NF_INQ_VARID(nf_fid,'waterVapor',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var waterVapor'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,waterVapor)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var waterVapor'
      endif
C
C     Variable        NETCDF Long Name
C      wetDelay     "Wet component GPS signal delay"
C
        nf_status = NF_INQ_VARID(nf_fid,'wetDelay',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var wetDelay'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,wetDelay)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var wetDelay'
      endif

C   Variables of type INT
C

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
C      staLongNam   "Station Location"
C
        nf_status = NF_INQ_VARID(nf_fid,'staLongNam',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var staLongNam'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,staLongNam)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var staLongNam'
      endif
C
C     Variable        NETCDF Long Name
C      staNam       "Alphanumeric station name"
C
        nf_status = NF_INQ_VARID(nf_fid,'staNam',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var staNam'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,staNam)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var staNam'
      endif

      nf_status = nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
      endif

      return
      end
