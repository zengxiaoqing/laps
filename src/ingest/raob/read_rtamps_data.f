C
C  Subroutine to read the file 
C
      subroutine read_rtamps_netcdf_dum(nfid, manLevel, maxStaticIds,        
     +     nInventoryBins, rawLevel, recNum, stdLevel, termLevel, 
     +     tropLevel, editFlag, firstInBin, firstOverflow, 
     +     globalInventory, indxRefr, invTime, inventory, irMan, 
     +     isOverflow, lastInBin, lastRecord, nStaticIds, oiMan, 
     +     optIndxRefr, prevRecord, storedObs, absHumidity, 
     +     airDensity, baromPressure, bpMan, bpStd, bpTerm, bpTrop, 
     +     dewPt, direction, dpMan, dpStd, dpTrop, drMan, drStd, 
     +     elevation, geomHeight, geopHeight, ghMan, ghTerm, ghTrop, 
     +     gpStd, gpTerm, latitude, longitude, precipWater, 
     +     relHumidity, rhMan, rhStd, riseRate, shear, shearDir, 
     +     shearMagX, shearMagY, spMan, spStd, speed, temperature, 
     +     tpMan, tpStd, tpTrop, vaporPressure, velError, velSound, 
     +     observationTime, receivedTime, reportTime, dataProvider, 
     +     providerId, staticIds, stationName)
C
      include 'netcdf.inc'
      integer manLevel, maxStaticIds, nInventoryBins, rawLevel, 
     +     recNum, stdLevel, termLevel, tropLevel,nfid, nf_vid, 
     +     nf_status
      integer editFlag( rawLevel, recNum), firstInBin(nInventoryBins),
     +     firstOverflow, globalInventory, indxRefr( rawLevel,
     +     recNum), invTime(recNum), inventory(maxStaticIds), irMan(
     +     manLevel, recNum), isOverflow(recNum),
     +     lastInBin(nInventoryBins), lastRecord(maxStaticIds),
     +     nStaticIds, oiMan( manLevel, recNum), optIndxRefr(
     +     rawLevel, recNum), prevRecord(recNum), storedObs(recNum)
      real absHumidity( rawLevel, recNum), airDensity( rawLevel,
     +     recNum), baromPressure( rawLevel, recNum), bpMan(
     +     manLevel, recNum), bpStd( stdLevel, recNum), bpTerm(
     +     termLevel, recNum), bpTrop( tropLevel, recNum), dewPt(
     +     rawLevel, recNum), direction( rawLevel, recNum), dpMan(
     +     manLevel, recNum), dpStd( stdLevel, recNum), dpTrop(
     +     tropLevel, recNum), drMan( manLevel, recNum), drStd(
     +     stdLevel, recNum), elevation(recNum), geomHeight(
     +     rawLevel, recNum), geopHeight( rawLevel, recNum), ghMan(
     +     manLevel, recNum), ghTerm( termLevel, recNum), ghTrop(
     +     tropLevel, recNum), gpStd( stdLevel, recNum), gpTerm(
     +     termLevel, recNum), latitude(recNum), longitude(recNum),
     +     precipWater( rawLevel, recNum), relHumidity( rawLevel,
     +     recNum), rhMan( manLevel, recNum), rhStd( stdLevel,
     +     recNum), riseRate( rawLevel, recNum), shear( rawLevel,
     +     recNum), shearDir( rawLevel, recNum), shearMagX( rawLevel,
     +     recNum), shearMagY( rawLevel, recNum), spMan( manLevel,
     +     recNum), spStd( stdLevel, recNum), speed( rawLevel,
     +     recNum), temperature( rawLevel, recNum), tpMan( manLevel,
     +     recNum), tpStd( stdLevel, recNum), tpTrop( tropLevel,
     +     recNum), vaporPressure( rawLevel, recNum), velError(
     +     rawLevel, recNum), velSound( rawLevel, recNum)
      double precision observationTime(recNum), receivedTime(recNum),
     +     reportTime(recNum)
      character*12 providerId(recNum)
      character*24 staticIds(maxStaticIds)
      character*51 stationName(recNum)
      character*11 dataProvider(recNum)


C   Variables of type REAL
C
C     Variable        NETCDF Long Name
C      absHumidity  "Absolute Humidity"
C
        nf_status = NF_INQ_VARID(nfid,'absHumidity',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var absHumidity'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,absHumidity)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var absHumidity'
      endif
C
C     Variable        NETCDF Long Name
C      airDensity   "Density of Air"
C
        nf_status = NF_INQ_VARID(nfid,'airDensity',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var airDensity'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,airDensity)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var airDensity'
      endif
C
C     Variable        NETCDF Long Name
C      baromPressure"Pressure"
C
        nf_status = NF_INQ_VARID(nfid,'baromPressure',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var baromPressure'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,baromPressure)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var baromPressure'
      endif
C
C     Variable        NETCDF Long Name
C      bpMan        "Pressure - Mandatory Level"
C
        nf_status = NF_INQ_VARID(nfid,'bpMan',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var bpMan'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,bpMan)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var bpMan'
      endif
C
C     Variable        NETCDF Long Name
C      bpStd        "Pressure - Standard Level"
C
        nf_status = NF_INQ_VARID(nfid,'bpStd',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var bpStd'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,bpStd)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var bpStd'
      endif
C
C     Variable        NETCDF Long Name
C      bpTerm       "Pressure - Termination"
C
        nf_status = NF_INQ_VARID(nfid,'bpTerm',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var bpTerm'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,bpTerm)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var bpTerm'
      endif
C
C     Variable        NETCDF Long Name
C      bpTrop       "Pressure - Tropopause Level"
C
        nf_status = NF_INQ_VARID(nfid,'bpTrop',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var bpTrop'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,bpTrop)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var bpTrop'
      endif
C
C     Variable        NETCDF Long Name
C      dewPt        "Dew Point Temperature"
C
        nf_status = NF_INQ_VARID(nfid,'dewPt',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dewPt'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,dewPt)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dewPt'
      endif
C
C     Variable        NETCDF Long Name
C      direction    "Wind Direction"
C
        nf_status = NF_INQ_VARID(nfid,'direction',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var direction'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,direction)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var direction'
      endif
C
C     Variable        NETCDF Long Name
C      dpMan        "Dew Point Temperature - Mandatory Level"
C
        nf_status = NF_INQ_VARID(nfid,'dpMan',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dpMan'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,dpMan)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dpMan'
      endif
C
C     Variable        NETCDF Long Name
C      dpStd        "Dew Point Temperature - Standard Level"
C
        nf_status = NF_INQ_VARID(nfid,'dpStd',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dpStd'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,dpStd)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dpStd'
      endif
C
C     Variable        NETCDF Long Name
C      dpTrop       "Dew Point Temperature - Tropopause Level"
C
        nf_status = NF_INQ_VARID(nfid,'dpTrop',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dpTrop'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,dpTrop)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dpTrop'
      endif
C
C     Variable        NETCDF Long Name
C      drMan        "Wind Direction - Mandatory Level"
C
        nf_status = NF_INQ_VARID(nfid,'drMan',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var drMan'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,drMan)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var drMan'
      endif
C
C     Variable        NETCDF Long Name
C      drStd        "Wind Direction - Standard Level"
C
        nf_status = NF_INQ_VARID(nfid,'drStd',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var drStd'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,drStd)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var drStd'
      endif
C
C     Variable        NETCDF Long Name
C      elevation    "Station Elevation"
C
        nf_status = NF_INQ_VARID(nfid,'elevation',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var elevation'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,elevation)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var elevation'
      endif
C
C     Variable        NETCDF Long Name
C      geomHeight   "Geometric Height"
C
        nf_status = NF_INQ_VARID(nfid,'geomHeight',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var geomHeight'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,geomHeight)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var geomHeight'
      endif
C
C     Variable        NETCDF Long Name
C      geopHeight   "Geopotential Height"
C
        nf_status = NF_INQ_VARID(nfid,'geopHeight',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var geopHeight'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,geopHeight)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var geopHeight'
      endif
C
C     Variable        NETCDF Long Name
C      ghMan        "Geometric - Mandatory Level"
C
        nf_status = NF_INQ_VARID(nfid,'ghMan',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var ghMan'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,ghMan)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var ghMan'
      endif
C
C     Variable        NETCDF Long Name
C      ghTerm       "Geometric - Termination"
C
        nf_status = NF_INQ_VARID(nfid,'ghTerm',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var ghTerm'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,ghTerm)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var ghTerm'
      endif
C
C     Variable        NETCDF Long Name
C      ghTrop       "Geometric - Tropopause Level"
C
        nf_status = NF_INQ_VARID(nfid,'ghTrop',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var ghTrop'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,ghTrop)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var ghTrop'
      endif
C
C     Variable        NETCDF Long Name
C      gpStd        "Geopotential - Standard Level"
C
        nf_status = NF_INQ_VARID(nfid,'gpStd',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var gpStd'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,gpStd)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var gpStd'
      endif
C
C     Variable        NETCDF Long Name
C      gpTerm       "Geopotential - Termination"
C
        nf_status = NF_INQ_VARID(nfid,'gpTerm',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var gpTerm'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,gpTerm)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var gpTerm'
      endif
C
C     Variable        NETCDF Long Name
C      latitude     "Station Latitude"
C
        nf_status = NF_INQ_VARID(nfid,'latitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var latitude'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,latitude)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var latitude'
      endif
C
C     Variable        NETCDF Long Name
C      longitude    "Station Longitude"
C
        nf_status = NF_INQ_VARID(nfid,'longitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var longitude'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,longitude)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var longitude'
      endif
C
C     Variable        NETCDF Long Name
C      precipWater  "Precipitable Water"
C
        nf_status = NF_INQ_VARID(nfid,'precipWater',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipWater'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,precipWater)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipWater'
      endif
C
C     Variable        NETCDF Long Name
C      relHumidity  "Relative Humidity"
C
        nf_status = NF_INQ_VARID(nfid,'relHumidity',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var relHumidity'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,relHumidity)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var relHumidity'
      endif
C
C     Variable        NETCDF Long Name
C      rhMan        "Relative Humidity - Mandatory Level"
C
        nf_status = NF_INQ_VARID(nfid,'rhMan',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var rhMan'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,rhMan)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var rhMan'
      endif
C
C     Variable        NETCDF Long Name
C      rhStd        "Relative Humidity - Standard Level"
C
        nf_status = NF_INQ_VARID(nfid,'rhStd',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var rhStd'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,rhStd)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var rhStd'
      endif
C
C     Variable        NETCDF Long Name
C      riseRate     "Rise Rate"
C
        nf_status = NF_INQ_VARID(nfid,'riseRate',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var riseRate'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,riseRate)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var riseRate'
      endif
C
C     Variable        NETCDF Long Name
C      shear        "Shear"
C
        nf_status = NF_INQ_VARID(nfid,'shear',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var shear'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,shear)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var shear'
      endif
C
C     Variable        NETCDF Long Name
C      shearDir     "Shear Direction"
C
        nf_status = NF_INQ_VARID(nfid,'shearDir',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var shearDir'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,shearDir)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var shearDir'
      endif
C
C     Variable        NETCDF Long Name
C      shearMagX    "Shear Magnitude X-direction"
C
        nf_status = NF_INQ_VARID(nfid,'shearMagX',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var shearMagX'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,shearMagX)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var shearMagX'
      endif
C
C     Variable        NETCDF Long Name
C      shearMagY    "Shear Magnitude Y-direction"
C
        nf_status = NF_INQ_VARID(nfid,'shearMagY',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var shearMagY'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,shearMagY)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var shearMagY'
      endif
C
C     Variable        NETCDF Long Name
C      spMan        "Wind Speed - Mandatory Level"
C
        nf_status = NF_INQ_VARID(nfid,'spMan',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var spMan'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,spMan)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var spMan'
      endif
C
C     Variable        NETCDF Long Name
C      spStd        "Wind Speed - Standard Level"
C
        nf_status = NF_INQ_VARID(nfid,'spStd',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var spStd'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,spStd)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var spStd'
      endif
C
C     Variable        NETCDF Long Name
C      speed        "Wind Speed"
C
        nf_status = NF_INQ_VARID(nfid,'speed',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var speed'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,speed)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var speed'
      endif
C
C     Variable        NETCDF Long Name
C      temperature  "Temperature"
C
        nf_status = NF_INQ_VARID(nfid,'temperature',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var temperature'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,temperature)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var temperature'
      endif
C
C     Variable        NETCDF Long Name
C      tpMan        "Temperature - Mandatory Level"
C
        nf_status = NF_INQ_VARID(nfid,'tpMan',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var tpMan'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,tpMan)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var tpMan'
      endif
C
C     Variable        NETCDF Long Name
C      tpStd        "Temperature - Standard Level"
C
        nf_status = NF_INQ_VARID(nfid,'tpStd',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var tpStd'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,tpStd)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var tpStd'
      endif
C
C     Variable        NETCDF Long Name
C      tpTrop       "Temperature - Tropopause Level"
C
        nf_status = NF_INQ_VARID(nfid,'tpTrop',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var tpTrop'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,tpTrop)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var tpTrop'
      endif
C
C     Variable        NETCDF Long Name
C      vaporPressure"Pressure"
C
        nf_status = NF_INQ_VARID(nfid,'vaporPressure',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var vaporPressure'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,vaporPressure)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var vaporPressure'
      endif
C
C     Variable        NETCDF Long Name
C      velError     "Velocity Error"
C
        nf_status = NF_INQ_VARID(nfid,'velError',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var velError'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,velError)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var velError'
      endif
C
C     Variable        NETCDF Long Name
C      velSound     "Velocity of Sound"
C
        nf_status = NF_INQ_VARID(nfid,'velSound',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var velSound'
      endif
        nf_status = NF_GET_VAR_REAL(nfid,nf_vid,velSound)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var velSound'
      endif

C   Variables of type INT
C
C
C     Variable        NETCDF Long Name
C      editFlag     "Edit Flag"
C
        nf_status = NF_INQ_VARID(nfid,'editFlag',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var editFlag'
      endif
        nf_status = NF_GET_VAR_INT(nfid,nf_vid,editFlag)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var editFlag'
      endif
C
C     Variable        NETCDF Long Name
C      firstInBin   
C
        nf_status = NF_INQ_VARID(nfid,'firstInBin',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var firstInBin'
      endif
        nf_status = NF_GET_VAR_INT(nfid,nf_vid,firstInBin)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var firstInBin'
      endif
C
C     Variable        NETCDF Long Name
C      firstOverflow
C
        nf_status = NF_INQ_VARID(nfid,'firstOverflow',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var firstOverflow'
      endif
        nf_status = NF_GET_VAR_INT(nfid,nf_vid,firstOverflow)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var firstOverflow'
      endif
C
C     Variable        NETCDF Long Name
C      globalInventory
C
        nf_status = NF_INQ_VARID(nfid,'globalInventory',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var globalInventory'
      endif
        nf_status = NF_GET_VAR_INT(nfid,nf_vid,globalInventory)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var globalInventory'
      endif
C
C     Variable        NETCDF Long Name
C      indxRefr     "Microwave Index of Refraction"
C
        nf_status = NF_INQ_VARID(nfid,'indxRefr',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var indxRefr'
      endif
        nf_status = NF_GET_VAR_INT(nfid,nf_vid,indxRefr)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var indxRefr'
      endif
C
C     Variable        NETCDF Long Name
C      invTime      
C
        nf_status = NF_INQ_VARID(nfid,'invTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var invTime'
      endif
        nf_status = NF_GET_VAR_INT(nfid,nf_vid,invTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var invTime'
      endif
C
C     Variable        NETCDF Long Name
C      inventory    
C
        nf_status = NF_INQ_VARID(nfid,'inventory',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var inventory'
      endif
        nf_status = NF_GET_VAR_INT(nfid,nf_vid,inventory)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var inventory'
      endif
C
C     Variable        NETCDF Long Name
C      irMan        "Microwave Index of Refraction - Mandatory level"
C
        nf_status = NF_INQ_VARID(nfid,'irMan',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var irMan'
      endif
        nf_status = NF_GET_VAR_INT(nfid,nf_vid,irMan)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var irMan'
      endif
C
C     Variable        NETCDF Long Name
C      isOverflow   
C
        nf_status = NF_INQ_VARID(nfid,'isOverflow',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var isOverflow'
      endif
        nf_status = NF_GET_VAR_INT(nfid,nf_vid,isOverflow)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var isOverflow'
      endif
C
C     Variable        NETCDF Long Name
C      lastInBin    
C
        nf_status = NF_INQ_VARID(nfid,'lastInBin',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var lastInBin'
      endif
        nf_status = NF_GET_VAR_INT(nfid,nf_vid,lastInBin)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var lastInBin'
      endif
C
C     Variable        NETCDF Long Name
C      lastRecord   
C
        nf_status = NF_INQ_VARID(nfid,'lastRecord',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var lastRecord'
      endif
        nf_status = NF_GET_VAR_INT(nfid,nf_vid,lastRecord)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var lastRecord'
      endif
C
C     Variable        NETCDF Long Name
C      nStaticIds   
C
        nf_status = NF_INQ_VARID(nfid,'nStaticIds',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var nStaticIds'
      endif
        nf_status = NF_GET_VAR_INT(nfid,nf_vid,nStaticIds)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var nStaticIds'
      endif
C
C     Variable        NETCDF Long Name
C      oiMan        "Optical Index of Refraction - Mandatory Level"
C
        nf_status = NF_INQ_VARID(nfid,'oiMan',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var oiMan'
      endif
        nf_status = NF_GET_VAR_INT(nfid,nf_vid,oiMan)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var oiMan'
      endif
C
C     Variable        NETCDF Long Name
C      optIndxRefr  "Optical Index of Refraction"
C
        nf_status = NF_INQ_VARID(nfid,'optIndxRefr',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var optIndxRefr'
      endif
        nf_status = NF_GET_VAR_INT(nfid,nf_vid,optIndxRefr)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var optIndxRefr'
      endif
C
C     Variable        NETCDF Long Name
C      prevRecord   
C
        nf_status = NF_INQ_VARID(nfid,'prevRecord',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var prevRecord'
      endif
        nf_status = NF_GET_VAR_INT(nfid,nf_vid,prevRecord)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var prevRecord'
      endif
C
C     Variable        NETCDF Long Name
C      storedObs    "Stored \'Raw\' Profile Observations"
C
        nf_status = NF_INQ_VARID(nfid,'storedObs',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var storedObs'
      endif
        nf_status = NF_GET_VAR_INT(nfid,nf_vid,storedObs)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var storedObs'
      endif

C   Variables of type DOUBLE
C
C
C     Variable        NETCDF Long Name
C      observationTime"Observation Time"
C
        nf_status = NF_INQ_VARID(nfid,'observationTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var observationTime'
      endif
        nf_status = NF_GET_VAR_DOUBLE(nfid,nf_vid,observationTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var observationTime'
      endif
C
C     Variable        NETCDF Long Name
C      receivedTime "Received Time"
C
        nf_status = NF_INQ_VARID(nfid,'receivedTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var receivedTime'
      endif
        nf_status = NF_GET_VAR_DOUBLE(nfid,nf_vid,receivedTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var receivedTime'
      endif
C
C     Variable        NETCDF Long Name
C      reportTime   "Report Time"
C
        nf_status = NF_INQ_VARID(nfid,'reportTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var reportTime'
      endif
        nf_status = NF_GET_VAR_DOUBLE(nfid,nf_vid,reportTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var reportTime'
      endif


C   Variables of type CHAR
C
C
C     Variable        NETCDF Long Name
C      dataProvider "Local data provider"
C
        nf_status = NF_INQ_VARID(nfid,'dataProvider',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dataProvider'
      endif
        nf_status = NF_GET_VAR_TEXT(nfid,nf_vid,dataProvider)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dataProvider'
      endif
C
C     Variable        NETCDF Long Name
C      providerId   "Data Provider station Id"
C
        nf_status = NF_INQ_VARID(nfid,'providerId',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var providerId'
      endif
        nf_status = NF_GET_VAR_TEXT(nfid,nf_vid,providerId)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var providerId'
      endif
C
C     Variable        NETCDF Long Name
C      staticIds    
C
        nf_status = NF_INQ_VARID(nfid,'staticIds',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var staticIds'
      endif
        nf_status = NF_GET_VAR_TEXT(nfid,nf_vid,staticIds)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var staticIds'
      endif
C
C     Variable        NETCDF Long Name
C      stationName  "alphanumeric station name"
C
        nf_status = NF_INQ_VARID(nfid,'stationName',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationName'
      endif
        nf_status = NF_GET_VAR_TEXT(nfid,nf_vid,stationName)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationName'
      endif

      nf_status = nf_close(nfid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
      endif

      return
      end
