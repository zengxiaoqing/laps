cdis    Forecast Systems Laboratory
cdis    NOAA/OAR/ERL/FSL
cdis    325 Broadway
cdis    Boulder, CO     80303
cdis 
cdis    Forecast Research Division
cdis    Local Analysis and Prediction Branch
cdis    LAPS 
cdis 
cdis    This software and its documentation are in the public domain and 
cdis    are furnished "as is."  The United States government, its 
cdis    instrumentalities, officers, employees, and agents make no 
cdis    warranty, express or implied, as to the usefulness of the software 
cdis    and documentation for any purpose.  They assume no responsibility 
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis    
cdis    Permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  All modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making 
cdis    the modifications.  If significant modifications or enhancements 
cdis    are made to this software, the FSL Software Policy Manager  
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
c
c
	subroutine write_surface_obs(btime,outfile,n_obs_g,
     &    n_obs_b,wmoid,stations,provider,wx,reptype,autostntype,
     &    store_1,store_2,store_3,store_4,store_5,store_6,store_7,
     &    store_2ea,store_3ea,store_4ea,store_5ea,store_6ea,
     &    store_cldamt,store_cldht,maxsta,jstatus)
c
c*****************************************************************************
c
c	Routine to write the LAPS surface data file.   The data is passed
c       to this routine via the 'store' array.
c
c	Changes:
c		P. Stamus  03-27-98  Original version (from old format
c                                      version of write_surface_obs).
c                          05-01-98  Added soil moisture to 'store_5' & '_5ea'
c                          09-04-98  Final adjustments for operational use.
c
c
c               J. Edwards 09-16-98  moved all format definitions to 
c                                    src/include/lso_formats.inc
c                                    changed 909 definition to allow for 
c                                    missing data
c
c*****************************************************************************
c
	real*4 store_1(maxsta,4), 
     &         store_2(maxsta,3), store_2ea(maxsta,3),
     &         store_3(maxsta,4), store_3ea(maxsta,2),
     &         store_4(maxsta,5), store_4ea(maxsta,2),
     &         store_5(maxsta,4), store_5ea(maxsta,4),
     &         store_6(maxsta,5), store_6ea(maxsta,2),
     &         store_7(maxsta,3),
     &         store_cldht(maxsta,5)
c
	integer*4 jstatus, wmoid(maxsta)
c
	character btime*24, outfile*(*), 
     &         stations(maxsta)*20, provider(maxsta)*11,
     &         wx(maxsta)*25,reptype(maxsta)*6, 
     &         autostntype(maxsta)*6,store_cldamt(maxsta,5)*4
c
c
c.....	Write the file.
c
	open(11,file=outfile,status='unknown')
c
c.....	Write the header.
c
	write(11,900) btime,		! time
     &               n_obs_g,		! # of obs in the laps grid
     &               n_obs_b		! # of obs in the box
c
c.....	Write the station data.
c
	do k=1,n_obs_b
c
	   write(11,901) stations(k),           !station id
     &                   wmoid(k),              !WMO id number
     &                   provider(k),           !data provider
     &                   (store_1(k,i),i=1,3),  !lat, lon, elev
     &                   nint(store_1(k,4))     !obs time
c
	  write(11,903)  reptype(k),            !station report type
     &                   autostntype(k),        !station type (manual/auto)
     &                   wx(k)                  !present weather
c
	  write(11,905) store_2(k,1), store_2ea(k,1),   !temp, temp expected accuracy
     &                  store_2(k,2), store_2ea(k,2),   !dew point, dew point exp. accuracy
     &                  store_2(k,3), store_2ea(k,3)    !Rel hum, rh expected accuracy
c
	  write(11,907) store_3(k,1), store_3(k,2),     !wind dir, wind speed
     &                  store_3(k,3), store_3(k,4),     !wind gust dir, wind gust speed
     &                  store_3ea(k,1), store_3ea(k,2)  !dir expected accuracy, spd exp accuracy
c
	  write(11,909) store_4(k,1),                   !altimeter
     &                  store_4(k,2),                   !station pressure
     &                  store_4(k,3),                   !MSL pressure
     &                  nint(store_4(k,4)),             !3-h press change character
     &                  store_4(k,5),                   !3-h pressure change
     &                  store_4ea(k,1), store_4ea(k,2)  !pressure exp accuracy, alt exp accuracy


c
	  write(11,911) store_5(k,1), store_5ea(k,1),   !visibility, vis exp accuracy
     &                  store_5(k,2), store_5ea(k,2),   !solar, solar exp accuracy
     &                  store_5(k,3), store_5ea(k,3),   !soil/water temp, soil/water temp exp accuracy
     &                  store_5(k,4), store_5ea(k,4)    !soil moisture, soil moist exp accuracy
c
	  write(11,913)  store_6(k,1),                  !1-h precipitation
     &                   store_6(k,2),                  !3-h precipitation
     &                   store_6(k,3),                  !6-h precipitation
     &                   store_6(k,4),                  !24-h precipitation
     &                   store_6(k,5),                  !snow depth
     &                   store_6ea(k,1), store_6ea(k,2) !precip and snow exp accuracy
c
	  kkk_s = int(store_7(k,1))
	  write(11,915) kkk_s,                     !num cld layers (store_7(k,1))
     &                  store_7(k,2),              !24-h max temperature
     &                  store_7(k,3)               !24-h min temperature
c
c.....	Write the cloud data if we have any.
c
	  if(kkk_s .gt. 0) then
	    do ii=1,kkk_s
  	      write(11,917) store_cldamt(k,ii), store_cldht(k,ii)   !layer cloud amount and height
	    enddo !ii
	  endif
c
	enddo !k
c
	endfile(11)
	close(11)	

c
c..... End of data writing.  Let's go home...
c
	return
	include 'lso_formats.inc'
	end

        subroutine get_sfc_badflag(badflag_out,istatus)

cdoc    Returns "badflag" used in surface code

        include 'laps_sfc.inc'

        badflag_out = badflag

        istatus = 1
        return
        end

        subroutine get_ibadflag(ibadflag,istatus)

cdoc    Returns "ibadflag" used in surface code

        ibadflag = -99
        istatus = 1

        return
        end

c
c
      subroutine read_local(nf_fid , recNum, altimeter, dataProvider,
     +     dewpoint, elevation, latitude, longitude, observationTime,
     +     presWeather, relHumidity, rhChangeTime, seaLevelPressure,
     +     stationId, stationPressChangeTime, stationPressure,
     +     stationType, tempChangeTime, temperature, visibility,
     +     windDir, windDirChangeTime, windDirMax, windGust,
     +     windGustChangeTime, windSpeed, windSpeedChangeTime,
     &     badflag, istatus)
c
c**********************************************************************
c
c     Routine to read the LDAD NetCDF mesonet observation files at FSL.
c     Code created with 'xgennet.pl' by J. Edwards, NOAA/FSL.
c     
c     Original:  P. Stamus, NOAA/FSL  28 Aug 1998
c	         09-30-98  P. Stamus
c                     Housekeeping changes.
c                10-18-99  P. Stamus
c                     Add check for missing values.
c
c**********************************************************************
c
      include 'netcdf.inc'
      integer recNum, nf_fid, nf_vid, nf_status !, ifilval

      character*11 dataProvider(recNum)
      character*25 presWeather(recNum)
      character*6 stationId(recNum)
      character*11 stationType(recNum)
      real altimeter(recNum), dewpoint(recNum), elevation(recNum),
     +     latitude(recNum), longitude(recNum), relHumidity(recNum),
     +     seaLevelPressure(recNum), stationPressure(recNum),
     +     temperature(recNum), visibility(recNum), windDir(recNum),
     +     windDirMax(recNum), windGust(recNum), windSpeed(recNum)
      real filval, misval

      double precision observationTime(recNum), rhChangeTime(recNum),
     +     stationPressChangeTime(recNum), tempChangeTime(recNum),
     +     windDirChangeTime(recNum), windGustChangeTime(recNum),
     +     windSpeedChangeTime(recNum), dfilval, dmisval
c
c..... Start here.
c
      istatus = 0
C
C     Variable        NETCDF Long Name
C      dataProvider "LDAD data provider" 
C
        nf_status = NF_INQ_VARID(nf_fid,'dataProvider',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dataProvider'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,dataProvider)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ dataProvider '
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
C      stationId    "alphanumeric station Id" 
C
        nf_status = NF_INQ_VARID(nf_fid,'stationId',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationId'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,stationId)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ stationId '
      endif
C
C     Variable        NETCDF Long Name
C      stationType  "LDAD station type" 
C
        nf_status = NF_INQ_VARID(nf_fid,'stationType',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationType'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,stationType)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ stationType '
      endif
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
       nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var altimeter'
      endif
      call ck_array_real(altimeter, recNum, misval, badflag)
C
C     Variable        NETCDF Long Name
C      dewpoint     "dew point temperature" 
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
       nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var dewpoint'
      endif
      call ck_array_real(dewpoint, recNum, misval, badflag)
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
         print *, ' in var elevation'
      endif
      call ck_array_real(elevation, recNum, filval, badflag)
       nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var elevation'
      endif
      call ck_array_real(elevation, recNum, misval, badflag)
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
       nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var latitude'
      endif
      call ck_array_real(latitude, recNum, misval, badflag)
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
       nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var longitude'
      endif
      call ck_array_real(longitude, recNum, misval, badflag)
C
C     Variable        NETCDF Long Name
C      relHumidity  "relative humidity" 
C
        nf_status = NF_INQ_VARID(nf_fid,'relHumidity',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var relHumidity'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,relHumidity)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ relHumidity '
      endif
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var relHumidity'
      endif
      call ck_array_real(relHumidity, recNum, filval, badflag)
       nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var relHumidity'
      endif
      call ck_array_real(relHumidity, recNum, misval, badflag)
C
C     Variable        NETCDF Long Name
C      seaLevelPressure"sea level pressure" 
C
        nf_status = NF_INQ_VARID(nf_fid,'seaLevelPressure',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var seaLevelPressure'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,seaLevelPressure)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ seaLevelPressure '
      endif
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var seaLevelPressure'
      endif
      call ck_array_real(seaLevelPressure, recNum, filval, badflag)
       nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var seaLevelPressure'
      endif
      call ck_array_real(seaLevelPressure, recNum, misval, badflag)
C
C     Variable        NETCDF Long Name
C      stationPressure"station pressure" 
C
        nf_status = NF_INQ_VARID(nf_fid,'stationPressure',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationPressure'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,stationPressure)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ stationPressure '
      endif
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var stationPressure'
      endif
      call ck_array_real(stationPressure, recNum, filval, badflag)
       nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var stationPressure'
      endif
      call ck_array_real(stationPressure, recNum, misval, badflag)
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
       nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var temperature'
      endif
      call ck_array_real(temperature, recNum, misval, badflag)
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
       nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var visibility'
      endif
      call ck_array_real(visibility, recNum, misval, badflag)
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
       nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var windDir'
      endif
      call ck_array_real(windDir, recNum, misval, badflag)
C
C     Variable        NETCDF Long Name
C      windDirMax   "wind direction at gust" 
C
        nf_status = NF_INQ_VARID(nf_fid,'windDirMax',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDirMax'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,windDirMax)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ windDirMax '
      endif
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var windDirMax'
      endif
      call ck_array_real(windDirMax, recNum, filval, badflag)
       nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var windDirMax'
      endif
      call ck_array_real(windDirMax, recNum, misval, badflag)
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
       nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var windGust'
      endif
      call ck_array_real(windGust, recNum, misval, badflag)
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
       nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var windSpeed'
      endif
      call ck_array_real(windSpeed, recNum, misval, badflag)
C
C     Variable        NETCDF Long Name
C      observationTime"time of observation" 
C
        nf_status = NF_INQ_VARID(nf_fid,'observationTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var observationTime'
      endif
        nf_status = NF_GET_VAR_DOUBLE(nf_fid,nf_vid,observationTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ observationTime '
      endif
       nf_status = NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'_FillValue',dfilval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var observationTime'
      endif
      do i=1,recNum
        if(observationTime(i) .eq. dfilval) observationTime(i) = badflag
      enddo !i
      nf_status=NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'missing_value',dmisval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var observationTime'
      endif
      do i=1,recNum
        if(observationTime(i) .eq. dmisval) observationTime(i) = badflag
      enddo !i
C
C     Variable        NETCDF Long Name
C      rhChangeTime "relative humidity time of last change" 
C
        nf_status = NF_INQ_VARID(nf_fid,'rhChangeTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var rhChangeTime'
      endif
        nf_status = NF_GET_VAR_DOUBLE(nf_fid,nf_vid,rhChangeTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ rhChangeTime '
      endif
       nf_status = NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'_FillValue',dfilval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var rhChangeTime'
      endif
      do i=1,recNum
         if(rhChangeTime(i) .eq. dfilval) rhChangeTime(i) = badflag
      enddo !i
      nf_status=NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'missing_value',dmisval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var rhChangeTime'
      endif
      do i=1,recNum
         if(rhChangeTime(i) .eq. dmisval) rhChangeTime(i) = badflag
      enddo !i
C
C     Variable        NETCDF Long Name
C      stationPressChangeTime"station press time of last change" 
C
        nf_status = NF_INQ_VARID(nf_fid,'stationPressChangeTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationPressChangeTime'
      endif
        nf_status = NF_GET_VAR_DOUBLE(nf_fid,nf_vid,
     &                                        stationPressChangeTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ stationPressChangeTime '
      endif
       nf_status = NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'_FillValue',dfilval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var stationPressChangeTime'
      endif
      do i=1,recNum
         if(stationPressChangeTime(i) .eq. dfilval) 
     &                             stationPressChangeTime(i) = badflag
      enddo !i
      nf_status=NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'missing_value',dmisval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var stationPressChangeTime'
      endif
      do i=1,recNum
         if(stationPressChangeTime(i) .eq. dmisval) 
     &                             stationPressChangeTime(i) = badflag
      enddo !i
C
C     Variable        NETCDF Long Name
C      tempChangeTime"temperature time of last change" 
C
        nf_status = NF_INQ_VARID(nf_fid,'tempChangeTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var tempChangeTime'
      endif
        nf_status = NF_GET_VAR_DOUBLE(nf_fid,nf_vid,tempChangeTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ tempChangeTime '
      endif
       nf_status = NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'_FillValue',dfilval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var tempChangeTime'
      endif
      do i=1,recNum
         if(tempChangeTime(i) .eq. dfilval) tempChangeTime(i) = badflag
      enddo !i
      nf_status=NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'missing_value',dmisval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var tempChangeTime'
      endif
      do i=1,recNum
         if(tempChangeTime(i) .eq. dmisval) tempChangeTime(i) = badflag
      enddo !i
C
C     Variable        NETCDF Long Name
C      windDirChangeTime"wind direction time of last change" 
C
        nf_status = NF_INQ_VARID(nf_fid,'windDirChangeTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDirChangeTime'
      endif
        nf_status = NF_GET_VAR_DOUBLE(nf_fid,nf_vid,windDirChangeTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ windDirChangeTime '
      endif
       nf_status = NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'_FillValue',dfilval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var windDirChangeTime'
      endif
      do i=1,recNum
         if(windDirChangeTime(i) .eq. dfilval) 
     &                                 windDirChangeTime(i) = badflag
      enddo !i
      nf_status=NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'missing_value',dmisval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var windDirChangeTime'
      endif
      do i=1,recNum
         if(windDirChangeTime(i) .eq. dmisval) 
     &                                 windDirChangeTime(i) = badflag
      enddo !i
C
C     Variable        NETCDF Long Name
C      windGustChangeTime"wind gust time of last change" 
C
        nf_status = NF_INQ_VARID(nf_fid,'windGustChangeTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windGustChangeTime'
      endif
        nf_status = NF_GET_VAR_DOUBLE(nf_fid,nf_vid,windGustChangeTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ windGustChangeTime '
      endif
       nf_status = NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'_FillValue',dfilval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var windGustChangeTime'
      endif
      do i=1,recNum
         if(windGustChangeTime(i) .eq. dfilval) 
     &                               windGustChangeTime(i) = badflag
      enddo !i
      nf_status=NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'missing_value',dmisval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var windGustChangeTime'
      endif
      do i=1,recNum
         if(windGustChangeTime(i) .eq. dmisval) 
     &                               windGustChangeTime(i) = badflag
      enddo !i
C
C     Variable        NETCDF Long Name
C      windSpeedChangeTime"wind speed time of last change" 
C
        nf_status = NF_INQ_VARID(nf_fid,'windSpeedChangeTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windSpeedChangeTime'
      endif
        nf_status = NF_GET_VAR_DOUBLE(nf_fid,nf_vid,windSpeedChangeTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ windSpeedChangeTime '
      endif
       nf_status = NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'_FillValue',dfilval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var windSpeedChangeTime'
      endif
      do i=1,recNum
         if(windSpeedChangeTime(i) .eq. dfilval) 
     &                                windSpeedChangeTime(i) = badflag
      enddo !i
      nf_status=NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'missing_value',dmisval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var windSpeedChangeTime'
      endif
      do i=1,recNum
         if(windSpeedChangeTime(i) .eq. dmisval) 
     &                                windSpeedChangeTime(i) = badflag
      enddo !i
c
      nf_status = nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
      endif
c
      istatus = 1
c
      return
      end

         subroutine get_filetime_range(
     1                i4time_ob_b,i4time_ob_a                   ! I
     1               ,i4_contains_early,i4_contains_late        ! I
     1               ,intvl                                     ! I
     1               ,i4time_file_b,i4time_file_a)              ! O

cdoc     Determine the range of needed filetimes, given observation time range
cdoc     and other info about the files.

         integer i4time_ob_b        ! Earliest ob we want
         integer i4time_ob_a        ! Latest ob we want
         integer i4_contains_early  ! Earliest contained ob relative to filetime
         integer i4_contains_late   ! Latest contained ob relative to filetime
         integer intvl              ! Regular time interval of files

!        Range of file times we want to read
         i4time_file_b = i4time_ob_b - i4_contains_late
         i4time_file_a = i4time_ob_a + i4_contains_early

!        Range of filenames at fixed intervals
         i4time_file_b = ( (i4time_file_b + intvl - 1) / intvl)*intvl     
         i4time_file_a = ( (i4time_file_a            ) / intvl)*intvl

         return
         end

c
c
      subroutine ck_array_real(var, recNum, filval, badflag)
c
      integer recNum
      real*4 var(recNum), filval, badflag
c
      do i=1,recNum
         if(var(i) .eq. filval) var(i) = badflag
      enddo !i
c
      return
      end
