      subroutine obs_driver_gps
      include 'netcdf.inc'
      integer recNum,nf_fid, nf_vid, nf_status
C
C  Open netcdf File for reading
C
      nf_status = NF_OPEN('fsl001061515.nc',NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'NF_OPEN fsl001061515.nc'
      endif
C
C  Fill all dimension values
C
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
      call get_gps_dum(nf_fid, recNum)

      end
C
C
      subroutine get_gps_dum(nf_fid, recNum)

      include 'netcdf.inc'
      integer recNum,nf_fid, nf_vid, nf_status

      real dryDelay(recNum), formalError(recNum), pressure(recNum),
     +     staElev(recNum), staLat(recNum), staLon(recNum),
     +     temperature(recNum), totalDelay(recNum),
     +     waterVapor(recNum), wetDelay(recNum)
      double precision timeObs(recNum)
      character*5 staNam(recNum)
      character*80 staLongNam(recNum)

      call read_gps(nf_fid, recNum, dryDelay, formalError, 
     +     pressure, staElev, staLat, staLon, temperature, 
     +     totalDelay, waterVapor, wetDelay, timeObs, staLongNam, 
     +     staNam)
C
C The netcdf variables are filled - your code goes here
C
      return
      end
C
C  Subroutine to read the file 
C
      subroutine read_gps(nf_fid, recNum, 
     +     pressure, staElev, staLat, staLon, temperature, 
     +     timeObs, staNam)
C
      include 'netcdf.inc'
      integer recNum,nf_fid, nf_vid, nf_status

      real dryDelay(recNum), formalError(recNum), pressure(recNum),
     +     staElev(recNum), staLat(recNum), staLon(recNum),
     +     temperature(recNum), totalDelay(recNum),
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
