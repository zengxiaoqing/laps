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
c
c
      subroutine read_metar(nf_fid , maxSkyCover, recNum, altimeter,
     &     autoStationType, dewpoint, dpFromTenths, elevation,
     &     latitude, longitude, maxTemp24Hour, minTemp24Hour,
     &     precip1Hour, precip24Hour, precip3Hour, precip6Hour,
     &     presWeather, pressChange3Hour, pressChangeChar,
     &     reportType, seaLevelPress, skyCover, skyLayerBase,
     &     snowCover, stationName, tempFromTenths, temperature,
     &     timeObs, visibility, windDir, windGust, windSpeed, wmoId,
     &     badflag, istatus)
c
c======================================================================
c
c     Routine to read the METAR NetCDF files at FSL.
c     Code created with 'xgennet.pl' by J. Edwards, NOAA/FSL.
c     
c     Original:  P. Stamus, NOAA/FSL  12 Mar 1998
c     Changes:
c        P. Stamus, NOAA/FSL  01 Sep 1998  Fix for differences in
c           /public vs WFO cdls (no snowCover in WFO at this time.)
c
c======================================================================
c
      include 'netcdf.inc'
c
      integer maxSkyCover, recNum, nf_fid, nf_vid, nf_status
      integer ifilval
c
      character*6 autoStationType(recNum)
      character*25 presWeather(recNum)
      character*6 reportType(recNum)
      character*8 skyCover( maxSkyCover, recNum)
      character*5 stationName(recNum)
      integer   pressChangeChar(recNum), wmoId(recNum)

      double precision timeObs(recNum), dfilval

      real altimeter(recNum), dewpoint(recNum), dpFromTenths(recNum)
      real elevation(recNum), latitude(recNum), longitude(recNum)
      real maxTemp24Hour(recNum), minTemp24Hour(recNum)
      real precip1Hour(recNum), precip24Hour(recNum)
      real precip3Hour(recNum), precip6Hour(recNum)
      real pressChange3Hour(recNum), seaLevelPress(recNum)
      real skyLayerBase( maxSkyCover, recNum), snowCover(recNum)
      real tempFromTenths(recNum), temperature(recNum)
      real visibility(recNum), windDir(recNum), windGust(recNum)
      real windSpeed(recNum)

      real filval
c
c
c..... Start.
c
      istatus = 0
C
C     Variable        NETCDF Long Name
C      autoStationType"automated station type" 
C
        nf_status = NF_INQ_VARID(nf_fid,'autoStationType',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var autoStationType'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,autoStationType)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ autoStationType '
      endif
C
C     Variable        NETCDF Long Name
C      presWeather  "present weather" 
C
        nf_status = NF_INQ_VARID(nf_fid,'presWeather',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var presWeather'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,presWeather)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ presWeather '
      endif
C
C     Variable        NETCDF Long Name
C      reportType   "report type" 
C
        nf_status = NF_INQ_VARID(nf_fid,'reportType',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var reportType'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,reportType)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ reportType '
      endif
C
C     Variable        NETCDF Long Name
C      skyCover     "sky cover" 
C
        nf_status = NF_INQ_VARID(nf_fid,'skyCover',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var skyCover'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,skyCover)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ skyCover '
      endif
C
C     Variable        NETCDF Long Name
C      stationName  "alphanumeric station identification" 
C
        nf_status = NF_INQ_VARID(nf_fid,'stationName',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationName'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,stationName)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ stationName '
      endif
C
C     Variable        NETCDF Long Name
C      pressChangeChar"character of pressure change" 
C
        nf_status = NF_INQ_VARID(nf_fid,'pressChangeChar',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var pressChangeChar'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,pressChangeChar)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ pressChangeChar '
      endif
        nf_status = NF_GET_ATT_INT(nf_fid,nf_vid,'_FillValue',ifilval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var pressChangeChar'
      endif
      do i=1,recNum
         if(pressChangeChar(i) .eq. ifilval) 
     &                         pressChangeChar(i) = int(badflag)
      enddo !i
C
C     Variable        NETCDF Long Name
C      wmoId        "numeric WMO identification" 
C
        nf_status = NF_INQ_VARID(nf_fid,'wmoId',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var wmoId'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,wmoId)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ wmoId '
      endif
        nf_status = NF_GET_ATT_INT(nf_fid,nf_vid,'_FillValue',ifilval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var wmoId'
      endif
      do i=1,recNum
         if(wmoId(i) .eq. ifilval) wmoId(i) = int(badflag)
      enddo !i
C
C     Variable        NETCDF Long Name
C      altimeter    "altimeter setting" 
C
        nf_status = NF_INQ_VARID(nf_fid,'altimeter',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var altimeter'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,altimeter)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ altimeter '
      endif
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var altimeter'
      endif
      call ck_array_real(altimeter, recNum, filval, badflag)
C
C     Variable        NETCDF Long Name
C      dewpoint     "dewpoint" 
C
        nf_status = NF_INQ_VARID(nf_fid,'dewpoint',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dewpoint'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,dewpoint)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ dewpoint '
      endif
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var dewpoint'
      endif
      call ck_array_real(dewpoint, recNum, filval, badflag)
C
C     Variable        NETCDF Long Name
C      dpFromTenths "dewpoint from tenths of a degree Celsius" 
C
        nf_status = NF_INQ_VARID(nf_fid,'dpFromTenths',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dpFromTenths'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,dpFromTenths)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ dpFromTenths '
      endif
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
        print *,'in var dpFromTenths'
      endif
      call ck_array_real(dpFromTenths, recNum, filval, badflag)
C
C     Variable        NETCDF Long Name
C      elevation    "elevation" 
C
        nf_status = NF_INQ_VARID(nf_fid,'elevation',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var elevation'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,elevation)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ elevation '
      endif
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *,'in var elevation'
      endif
      call ck_array_real(elevation, recNum, filval, badflag)
C
C     Variable        NETCDF Long Name
C      latitude     "latitude" 
C
        nf_status = NF_INQ_VARID(nf_fid,'latitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var latitude'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,latitude)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ latitude '
      endif
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var latitude' 
      endif
      call ck_array_real(latitude, recNum, filval, badflag)
C
C     Variable        NETCDF Long Name
C      longitude    "longitude" 
C
        nf_status = NF_INQ_VARID(nf_fid,'longitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var longitude'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,longitude)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ longitude '
      endif
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var longitude' 
      endif
      call ck_array_real(longitude, recNum, filval, badflag)
C
C     Variable        NETCDF Long Name
C      maxTemp24Hour"24 hour max temperature" 
C
        nf_status = NF_INQ_VARID(nf_fid,'maxTemp24Hour',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var maxTemp24Hour'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,maxTemp24Hour)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ maxTemp24Hour '
      endif
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *,'in var maxTemp24Hour'
      endif
      call ck_array_real(maxTemp24Hour, recNum, filval, badflag)
C
C     Variable        NETCDF Long Name
C      minTemp24Hour"24 hour min temperature" 
C
        nf_status = NF_INQ_VARID(nf_fid,'minTemp24Hour',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var minTemp24Hour'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,minTemp24Hour)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ minTemp24Hour '
      endif
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *,'in var minTemp24Hour'
      endif
      call ck_array_real(minTemp24Hour, recNum, filval, badflag)
C
C     Variable        NETCDF Long Name
C      precip1Hour  "1 hour precipitation" 
C
        nf_status = NF_INQ_VARID(nf_fid,'precip1Hour',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precip1Hour'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,precip1Hour)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ precip1Hour '
      endif
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var precip1Hour'
      endif
      call ck_array_real(precip1Hour, recNum, filval, badflag)
C
C     Variable        NETCDF Long Name
C      precip24Hour "24 hour precipitation" 
C
        nf_status = NF_INQ_VARID(nf_fid,'precip24Hour',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precip24Hour'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,precip24Hour)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ precip24Hour '
      endif
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var precip24Hour'
      endif
      call ck_array_real(precip24Hour, recNum, filval, badflag)
C
C     Variable        NETCDF Long Name
C      precip3Hour  "3 hour precipitation" 
C
        nf_status = NF_INQ_VARID(nf_fid,'precip3Hour',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precip3Hour'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,precip3Hour)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ precip3Hour '
      endif
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var precip3Hour'
      endif
      call ck_array_real(precip3Hour, recNum, filval, badflag)
C
C     Variable        NETCDF Long Name
C      precip6Hour  "6 hour precipitation" 
C
        nf_status = NF_INQ_VARID(nf_fid,'precip6Hour',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precip6Hour'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,precip6Hour)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ precip6Hour '
      endif
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var precip6Hour'
      endif
      call ck_array_real(precip6hour, recNum, filval, badflag)
C
C     Variable        NETCDF Long Name
C      pressChange3Hour"3 hour pressure change" 
C
        nf_status = NF_INQ_VARID(nf_fid,'pressChange3Hour',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var pressChange3Hour'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,pressChange3Hour)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ pressChange3Hour '
      endif
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var pressChange3Hour'
      endif
      call ck_array_real(pressChange3Hour, recNum, filval, badflag)
C
C     Variable        NETCDF Long Name
C      seaLevelPress"sea level pressure" 
C
        nf_status = NF_INQ_VARID(nf_fid,'seaLevelPress',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var seaLevelPress'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,seaLevelPress)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ seaLevelPress '
      endif
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var seaLevelPress'
      endif
      call ck_array_real(seaLevelPress, recNum, filval, badflag)
C
C     Variable        NETCDF Long Name
C      skyLayerBase "sky cover layer base" 
C
        nf_status = NF_INQ_VARID(nf_fid,'skyLayerBase',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var skyLayerBase'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,skyLayerBase)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ skyLayerBase '
      endif
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var skyLayerBase'
      endif
      call ck_array_real(skyLayerBase, recNum, filval, badflag)
C
C     Variable        NETCDF Long Name
C      snowCover    "snow cover" 
C
        nf_status = NF_INQ_VARID(nf_fid,'snowCover',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var snowCover'
        do i=1,recNum
           snowCover(i) = badflag
        enddo !i
        go to 500
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,snowCover)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ snowCover '
      endif
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var snowCover'
      endif
      call ck_array_real(snowCover, recNum, filval, badflag)
C
C     Variable        NETCDF Long Name
C      tempFromTenths"temperature from tenths of a degree Celsius" 
C
 500  nf_status = NF_INQ_VARID(nf_fid,'tempFromTenths',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var tempFromTenths'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,tempFromTenths)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ tempFromTenths '
      endif
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var tempFromTenths'
      endif
      call ck_array_real(tempFromTenths, recNum, filval, badflag)
C
C     Variable        NETCDF Long Name
C      temperature  "temperature" 
C
        nf_status = NF_INQ_VARID(nf_fid,'temperature',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var temperature'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,temperature)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ temperature '
      endif
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var temperature'
      endif
      call ck_array_real(temperature, recNum, filval, badflag)
C
C     Variable        NETCDF Long Name
C      visibility   "visibility" 
C
        nf_status = NF_INQ_VARID(nf_fid,'visibility',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var visibility'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,visibility)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ visibility '
      endif
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var visibility'
      endif
      call ck_array_real(visibility, recNum, filval, badflag)
C
C     Variable        NETCDF Long Name
C      windDir      "wind direction" 
C
        nf_status = NF_INQ_VARID(nf_fid,'windDir',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDir'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,windDir)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ windDir '
      endif
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var windDir'
      endif
      call ck_array_real(windDir, recNum, filval, badflag)
C
C     Variable        NETCDF Long Name
C      windGust     "wind gust" 
C
        nf_status = NF_INQ_VARID(nf_fid,'windGust',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windGust'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,windGust)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ windGust '
      endif
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var windGust'
      endif
      call ck_array_real(windGust, recNum, filval, badflag)
C
C     Variable        NETCDF Long Name
C      windSpeed    "wind speed" 
C
        nf_status = NF_INQ_VARID(nf_fid,'windSpeed',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windSpeed'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,windSpeed)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ windSpeed '
      endif
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var windSpeed'
      endif
      call ck_array_real(windSpeed, recNum, filval, badflag)
C
C     Variable        NETCDF Long Name
C      timeObs      "time of observation" 
C
        nf_status = NF_INQ_VARID(nf_fid,'timeObs',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var timeObs'
      endif
        nf_status = NF_GET_VAR_DOUBLE(nf_fid,nf_vid,timeObs)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ timeObs '
      endif
       nf_status = NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'_FillValue',dfilval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var timeObs'
      endif
      do i=1,recNum
         if(timeObs(i) .eq. dfilval) timeObs(i) = badflag
      enddo !i
c
c
      nf_status = nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
      endif
c
c..... That's it.
c
      istatus = 1
c
      return
      end
