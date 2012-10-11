
      subroutine get_vol_netcdf_hdr(filename,
     +        gateR, gateR_HI, gateV, gateV_HI, radialR, radialR_HI,
     +        radialV, radialV_HI, scanR, scanR_HI, scanV,
     +        scanV_HI,nf_fid, nf_vid, nf_status)

!     Argument List
      character*(*) filename

      include 'netcdf.inc'
      integer gateR, gateR_HI, gateV, gateV_HI, radialR, radialR_HI,
     +     radialV, radialV_HI, scanR, scanR_HI, scanV,
     +     scanV_HI,nf_fid, nf_vid, nf_status
C
C  Open netcdf File for reading
C
      nf_status = NF_OPEN(filename,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'NF_OPEN ',filename
      endif
C
C  Fill all dimension values
C
C
C Get size of gateR
C
      nf_status = NF_INQ_DIMID(nf_fid,'gateR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim gateR'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,gateR)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim gateR'
        gateR = 0
      endif
C
C Get size of gateR_HI
C
      nf_status = NF_INQ_DIMID(nf_fid,'gateR_HI',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim gateR_HI'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,gateR_HI)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim gateR_HI'
        gateR_HI = 0
      endif
C
C Get size of gateV
C
      nf_status = NF_INQ_DIMID(nf_fid,'gateV',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim gateV'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,gateV)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim gateV'
        gateV = 0
      endif
C
C Get size of gateV_HI
C
      nf_status = NF_INQ_DIMID(nf_fid,'gateV_HI',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim gateV_HI'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,gateV_HI)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim gateV_HI'
        gateV_HI = 0
      endif
C
C Get size of radialR
C
      nf_status = NF_INQ_DIMID(nf_fid,'radialR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim radialR'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,radialR)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim radialR'
        radialR = 0
      endif
C
C Get size of radialR_HI
C
      nf_status = NF_INQ_DIMID(nf_fid,'radialR_HI',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim radialR_HI'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,radialR_HI)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim radialR_HI'
        radialR_HI = 0
      endif
C
C Get size of radialV
C
      nf_status = NF_INQ_DIMID(nf_fid,'radialV',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim radialV'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,radialV)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim radialV'
        radialV = 0
      endif
C
C Get size of radialV_HI
C
      nf_status = NF_INQ_DIMID(nf_fid,'radialV_HI',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim radialV_HI'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,radialV_HI)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim radialV_HI'
        radialV_HI = 0
      endif
C
C Get size of scanR
C
      nf_status = NF_INQ_DIMID(nf_fid,'scanR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim scanR'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,scanR)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim scanR'
        scanR = 0
      endif
C
C Get size of scanR_HI
C
      nf_status = NF_INQ_DIMID(nf_fid,'scanR_HI',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim scanR_HI'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,scanR_HI)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim scanR_HI'
        scanR_HI = 0
      endif
C
C Get size of scanV
C
      nf_status = NF_INQ_DIMID(nf_fid,'scanV',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim scanV'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,scanV)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim scanV'
        scanV = 0
      endif
C
C Get size of scanV_HI
C
      nf_status = NF_INQ_DIMID(nf_fid,'scanV_HI',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim scanV_HI'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,scanV_HI)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim scanV_HI'
        scanV_HI = 0
      endif
!     call get_vol_netcdf_data(nf_fid, gateR, gateR_HI, gateV, gateV_HI, radialR,
!    +     radialR_HI, radialV, radialV_HI, scanR, scanR_HI, scanV,
!    +     scanV_HI)

      return
      end
C
C
      subroutine get_vol_netcdf_data(nf_fid, gateR, gateR_HI, gateV, 
     +     gateV_HI,
     +     radialR, radialR_HI, radialV, radialV_HI, scanR, scanR_HI,
     +     scanV, scanV_HI,
     +     Reflectivity, Reflectivity_HI,
     +     RadialVelocity, RadialVelocity_HI, 
     +     elevationR, elevationR_HI,
     +     elevationV, elevationV_HI,
     +     azimuthR, azimuthR_HI,
     +     azimuthV, azimuthV_HI,
     +     distanceR, distanceR_HI,
     +     distanceV, distanceV_HI,
     +     nyquistVelocityV, nyquistVelocityV_HI)

      include 'netcdf.inc'
      integer gateR, gateR_HI, gateV, gateV_HI, radialR, radialR_HI,
     +     radialV, radialV_HI, scanR, scanR_HI, scanV,
     +     scanV_HI,nf_fid, nf_vid, nf_status
      integer RadialVelocity( gateV,  radialV, scanV),
     +     RadialVelocity_HI( gateV_HI,  radialV_HI, scanV_HI),
     +     Reflectivity( gateR,  radialR, scanR), Reflectivity_HI(
     +     gateR_HI,  radialR_HI, scanR_HI), SpectrumWidth( gateV, 
     +     radialV, scanV), SpectrumWidth_HI( gateV_HI,  radialV_HI,
     +     scanV_HI), numGatesR(scanR), numGatesR_HI(scanR_HI),
     +     numGatesV(scanV), numGatesV_HI(scanV_HI),
     +     numRadialsR(scanR), numRadialsR_HI(scanR_HI),
     +     numRadialsV(scanV), numRadialsV_HI(scanV_HI), timeR(
     +     radialR, scanR), timeR_HI( radialR_HI, scanR_HI), timeV(
     +     radialV, scanV), timeV_HI( radialV_HI, scanV_HI)
      real azimuthR( radialR, scanR), azimuthR_HI( radialR_HI,
     +     scanR_HI), azimuthV( radialV, scanV), azimuthV_HI(
     +     radialV_HI, scanV_HI), distanceR(gateR),
     +     distanceR_HI(gateR_HI), distanceV(gateV),
     +     distanceV_HI(gateV_HI), elevationR( radialR, scanR),
     +     elevationR_HI( radialR_HI, scanR_HI), elevationV( radialV,
     +     scanV), elevationV_HI( radialV_HI, scanV_HI),
     +     nyquistVelocityV(scanV), nyquistVelocityV_HI(scanV_HI)


      call read_netcdf_vol(nf_fid, gateR, gateR_HI, gateV, gateV_HI, 
     +     radialR, radialR_HI, radialV, radialV_HI, scanR, scanR_HI, 
     +     scanV, scanV_HI, RadialVelocity, RadialVelocity_HI, 
     +     Reflectivity, Reflectivity_HI, SpectrumWidth, 
     +     SpectrumWidth_HI, numGatesR, numGatesR_HI, numGatesV, 
     +     numGatesV_HI, numRadialsR, numRadialsR_HI, numRadialsV, 
     +     numRadialsV_HI, timeR, timeR_HI, timeV, timeV_HI, 
     +     azimuthR, azimuthR_HI, azimuthV, azimuthV_HI, distanceR, 
     +     distanceR_HI, distanceV, distanceV_HI, elevationR, 
     +     elevationR_HI, elevationV, elevationV_HI,
     +     nyquistVelocityV, nyquistVelocityV_HI)
C
C The netcdf variables are filled - your code goes here
C
      return
      end
C
C  Subroutine to read the file 
C
      subroutine read_netcdf_vol(nf_fid, gateR, gateR_HI, gateV, 
     +     gateV_HI, 
     +     radialR, radialR_HI, radialV, radialV_HI, scanR, scanR_HI, 
     +     scanV, scanV_HI, RadialVelocity, RadialVelocity_HI, 
     +     Reflectivity, Reflectivity_HI, SpectrumWidth, 
     +     SpectrumWidth_HI, numGatesR, numGatesR_HI, numGatesV, 
     +     numGatesV_HI, numRadialsR, numRadialsR_HI, numRadialsV, 
     +     numRadialsV_HI, timeR, timeR_HI, timeV, timeV_HI, 
     +     azimuthR, azimuthR_HI, azimuthV, azimuthV_HI, distanceR, 
     +     distanceR_HI, distanceV, distanceV_HI, elevationR, 
     +     elevationR_HI, elevationV, elevationV_HI,
     +     nyquistVelocityV, nyquistVelocityV_HI)
C
   
      use mem_namelist, ONLY: r_missing_data

      include 'netcdf.inc'
      integer gateR, gateR_HI, gateV, gateV_HI, radialR, radialR_HI, 
     +     radialV, radialV_HI, scanR, scanR_HI, scanV, 
     +     scanV_HI,nf_fid, nf_vid, nf_status
      integer RadialVelocity( gateV,  radialV, scanV),
     +     RadialVelocity_HI( gateV_HI,  radialV_HI, scanV_HI),
     +     Reflectivity( gateR,  radialR, scanR), Reflectivity_HI(
     +     gateR_HI,  radialR_HI, scanR_HI), SpectrumWidth( gateV, 
     +     radialV, scanV), SpectrumWidth_HI( gateV_HI,  radialV_HI,
     +     scanV_HI), numGatesR(scanR), numGatesR_HI(scanR_HI),
     +     numGatesV(scanV), numGatesV_HI(scanV_HI),
     +     numRadialsR(scanR), numRadialsR_HI(scanR_HI),
     +     numRadialsV(scanV), numRadialsV_HI(scanV_HI), timeR(
     +     radialR, scanR), timeR_HI( radialR_HI, scanR_HI), timeV(
     +     radialV, scanV), timeV_HI( radialV_HI, scanV_HI)
      real azimuthR( radialR, scanR), azimuthR_HI( radialR_HI,
     +     scanR_HI), azimuthV( radialV, scanV), azimuthV_HI(
     +     radialV_HI, scanV_HI), distanceR(gateR),
     +     distanceR_HI(gateR_HI), distanceV(gateV),
     +     distanceV_HI(gateV_HI), elevationR( radialR, scanR),
     +     elevationR_HI( radialR_HI, scanR_HI), elevationV( radialV,
     +     scanV), elevationV_HI( radialV_HI, scanV_HI),
     +     nyquistVelocityR(scanR), nyquistVelocityR_HI(scanR_HI),
     +     nyquistVelocityV(scanV), nyquistVelocityV_HI(scanV_HI)


      i_missing_data = 255

C   Variables of type REAL
C
C     Variable        NETCDF Long Name
C      azimuthR     "azimuth angle in degrees: 0 = true north, 90 = east"
C
      nf_status = NF_INQ_VARID(nf_fid,'azimuthR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var azimuthR'
        azimuthR = r_missing_data
      else
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,azimuthR)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var azimuthR'
          azimuthR = r_missing_data
        endif
      endif
C
C     Variable        NETCDF Long Name
C      azimuthR_HI  "azimuth angle in degrees: 0 = true north, 90 = east"
C
      nf_status = NF_INQ_VARID(nf_fid,'azimuthR_HI',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var azimuthR_HI'
        azimuthR_HI = r_missing_data
      else
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,azimuthR_HI)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var azimuthR_HI'
          azimuthR_HI = r_missing_data
        endif
      endif
C
C     Variable        NETCDF Long Name
C      azimuthV     "azimuth angle in degrees: 0 = true north, 90 = east"
C
      nf_status = NF_INQ_VARID(nf_fid,'azimuthV',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var azimuthV'
        azimuthV = r_missing_data
      else
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,azimuthV)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var azimuthV'
          azimuthV = r_missing_data
        endif
      endif
C
C     Variable        NETCDF Long Name
C      azimuthV_HI  "azimuth angle in degrees: 0 = true north, 90 = east"
C
      nf_status = NF_INQ_VARID(nf_fid,'azimuthV_HI',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var azimuthV_HI'
        azimuthV_HI = r_missing_data
      else
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,azimuthV_HI)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var azimuthV_HI'
          azimuthV_HI = r_missing_data
        endif
      endif
C
C     Variable        NETCDF Long Name
C      distanceR    "radial distance to start of gate"
C
      nf_status = NF_INQ_VARID(nf_fid,'distanceR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var distanceR'
        distanceR = r_missing_data
      else
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,distanceR)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var distanceR'
          distanceR = r_missing_data
        endif
      endif
C
C     Variable        NETCDF Long Name
C      distanceR_HI "radial distance to start of gate"
C
      nf_status = NF_INQ_VARID(nf_fid,'distanceR_HI',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var distanceR_HI'
        distanceR_HI = r_missing_data
      else
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,distanceR_HI)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var distanceR_HI'
          distanceR_HI = r_missing_data
        endif
      endif
C
C     Variable        NETCDF Long Name
C      distanceV    "radial distance to start of gate"
C
      nf_status = NF_INQ_VARID(nf_fid,'distanceV',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var distanceV'
        distanceV = r_missing_data
      else
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,distanceV)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var distanceV'
          distanceV = r_missing_data
        endif
      endif
C
C     Variable        NETCDF Long Name
C      distanceV_HI "radial distance to start of gate"
C
      nf_status = NF_INQ_VARID(nf_fid,'distanceV_HI',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var distanceV_HI'
        distanceV_HI = r_missing_data
      else
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,distanceV_HI)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var distanceV_HI'
          distanceV_HI = r_missing_data
        endif
      endif
C
C     Variable        NETCDF Long Name
C      elevationR   "elevation angle in degres: 0 = parallel to pedestal base, 90 = perpendicular"
C
      nf_status = NF_INQ_VARID(nf_fid,'elevationR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var elevationR'
        elevationR = r_missing_data
      else
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,elevationR)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var elevationR'
          elevationR = r_missing_data
        endif
      endif
C
C     Variable        NETCDF Long Name
C      elevationR_HI"elevation angle in degres: 0 = parallel to pedestal base, 90 = perpendicular"
C
      nf_status = NF_INQ_VARID(nf_fid,'elevationR_HI',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var elevationR_HI'
        elevationR_HI = r_missing_data
      else
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,elevationR_HI)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var elevationR_HI'
          elevationR_HI = r_missing_data
        endif
      endif
C
C     Variable        NETCDF Long Name
C      elevationV   "elevation angle in degres: 0 = parallel to pedestal base, 90 = perpendicular"
C
      nf_status = NF_INQ_VARID(nf_fid,'elevationV',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var elevationV'
        elevationV = r_missing_data
      else
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,elevationV)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var elevationV'
          elevationV = r_missing_data
        endif
      endif
C
C     Variable        NETCDF Long Name
C      elevationV_HI"elevation angle in degres: 0 = parallel to pedestal base, 90 = perpendicular"
C
      nf_status = NF_INQ_VARID(nf_fid,'elevationV_HI',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var elevationV_HI'
        elevationV_HI = r_missing_data
      else
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,elevationV_HI)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var elevationV_HI'
          elevationV_HI = r_missing_data
        endif
      endif
C
C     Variable        NETCDF Long Name
C      nyquistVelocityR"Nyquist Velocity"
C
      nf_status = NF_INQ_VARID(nf_fid,'nyquistVelocityR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var nyquistVelocityR'
        nyquistVelocityR = r_missing_data
      else
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,nyquistVelocityR)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var nyquistVelocityR'
          nyquistVelocityR = r_missing_data
        endif
      endif
C
C     Variable        NETCDF Long Name
C      nyquistVelocityR_HI"Nyquist Velocity"
C
      nf_status = NF_INQ_VARID(nf_fid,'nyquistVelocityR_HI',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var nyquistVelocityR_HI'
        nyquistVelocityR_HI = r_missing_data
      else
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,nyquistVelocityR_HI)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var nyquistVelocityR_HI'
          nyquistVelocityR_HI = r_missing_data
        endif
      endif
C
C     Variable        NETCDF Long Name
C      nyquistVelocityV"Nyquist Velocity"
C
      nf_status = NF_INQ_VARID(nf_fid,'nyquistVelocityV',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var nyquistVelocityV'
        nyquistVelocityV = r_missing_data
      else
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,nyquistVelocityV)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var nyquistVelocityV'
          nyquistVelocityV = r_missing_data
        endif
      endif
C
C     Variable        NETCDF Long Name
C      nyquistVelocityV_HI"Nyquist Velocity"
C
      nf_status = NF_INQ_VARID(nf_fid,'nyquistVelocityV_HI',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var nyquistVelocityV_HI'
        nyquistVelocityV_HI = r_missing_data
      else
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,nyquistVelocityV_HI)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var nyquistVelocityV_HI'
          nyquistVelocityV_HI = r_missing_data
        endif
      endif

C   Variables of type INT
C
C
C     Variable        NETCDF Long Name
C      RadialVelocity"Radial Velocity"
C
      nf_status = NF_INQ_VARID(nf_fid,'RadialVelocity',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var RadialVelocity'
        RadialVelocity = i_missing_data
      else
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,RadialVelocity)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var RadialVelocity'
          RadialVelocity = i_missing_data
        endif
      endif
C
C     Variable        NETCDF Long Name
C      RadialVelocity_HI"Radial Velocity_HI"
C
      nf_status = NF_INQ_VARID(nf_fid,'RadialVelocity_HI',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var RadialVelocity_HI'
        RadialVelocity_HI = i_missing_data
      else
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,RadialVelocity_HI)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var RadialVelocity_HI'
          RadialVelocity_HI = i_missing_data
        endif
      endif
C
C     Variable        NETCDF Long Name
C      Reflectivity "Reflectivity"
C
      nf_status = NF_INQ_VARID(nf_fid,'Reflectivity',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Reflectivity'
        Reflectivity = i_missing_data
      else
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,Reflectivity)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var Reflectivity'
          Reflectivity = i_missing_data
        endif
      endif
C
C     Variable        NETCDF Long Name
C      Reflectivity_HI"Reflectivity_HI"
C
      nf_status = NF_INQ_VARID(nf_fid,'Reflectivity_HI',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Reflectivity_HI'
        Reflectivity_HI = i_missing_data
      else
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,Reflectivity_HI)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var Reflectivity_HI'
          Reflectivity_HI = i_missing_data
        endif
      endif
C
C     Variable        NETCDF Long Name
C      SpectrumWidth"Radial Spectrum"
C
      nf_status = NF_INQ_VARID(nf_fid,'SpectrumWidth',nf_vid)
      if(nf_status.ne.NF_NOERR .or. .true.) then
        print *, NF_STRERROR(nf_status)
        print *,'in var SpectrumWidth'
        SpectrumWidth = i_missing_data
      else
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,SpectrumWidth)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var SpectrumWidth'
          SpectrumWidth = i_missing_data
        endif
      endif
C
C     Variable        NETCDF Long Name
C      SpectrumWidth_HI"Radial Spectrum_HI"
C
      nf_status = NF_INQ_VARID(nf_fid,'SpectrumWidth_HI',nf_vid)
      if(nf_status.ne.NF_NOERR .or. .true.) then
        print *, NF_STRERROR(nf_status)
        print *,'in var SpectrumWidth_HI'
        SpectrumWidth_HI = i_missing_data
      else
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,SpectrumWidth_HI)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var SpectrumWidth_HI'
          SpectrumWidth_HI = i_missing_data
        endif
      endif
C
C     Variable        NETCDF Long Name
C      numGatesR    "number of valid gates in this scan"
C
      nf_status = NF_INQ_VARID(nf_fid,'numGatesR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var numGatesR'
        numGatesR = i_missing_data
      else
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,numGatesR)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var numGatesR'
          numGatesR = i_missing_data
        endif
      endif
C
C     Variable        NETCDF Long Name
C      numGatesR_HI "number of valid gates in this scan"
C
      nf_status = NF_INQ_VARID(nf_fid,'numGatesR_HI',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var numGatesR_HI'
        numGatesR_HI = i_missing_data
      else
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,numGatesR_HI)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var numGatesR_HI'
          numGatesR_HI = i_missing_data
        endif
      endif
C
C     Variable        NETCDF Long Name
C      numGatesV    "number of valid gates in this scan"
C
      nf_status = NF_INQ_VARID(nf_fid,'numGatesV',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var numGatesV'
        numGatesV = i_missing_data
      else
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,numGatesV)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var numGatesV'
          numGatesV = i_missing_data
        endif
      endif
C
C     Variable        NETCDF Long Name
C      numGatesV_HI "number of valid gates in this scan"
C
      nf_status = NF_INQ_VARID(nf_fid,'numGatesV_HI',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var numGatesV_HI'
        numGatesV_HI = i_missing_data
      else
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,numGatesV_HI)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var numGatesV_HI'
          numGatesV_HI = i_missing_data
        endif
      endif
C
C     Variable        NETCDF Long Name
C      numRadialsR  "number of valid radials in this scan"
C
      nf_status = NF_INQ_VARID(nf_fid,'numRadialsR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var numRadialsR'
        numRadialsR = i_missing_data
      else
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,numRadialsR)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var numRadialsR'
          numRadialsR = i_missing_data
        endif
      endif
C
C     Variable        NETCDF Long Name
C      numRadialsR_HI"number of valid radials in this scan"
C
      nf_status = NF_INQ_VARID(nf_fid,'numRadialsR_HI',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var numRadialsR_HI'
        numRadialsR_HI = i_missing_data
      else
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,numRadialsR_HI)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var numRadialsR_HI'
          numRadialsR_HI = i_missing_data
        endif
      endif
C
C     Variable        NETCDF Long Name
C      numRadialsV  "number of valid radials in this scan"
C
      nf_status = NF_INQ_VARID(nf_fid,'numRadialsV',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var numRadialsV'
        numRadialsV = i_missing_data
      else
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,numRadialsV)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var numRadialsV'
          numRadialsV = i_missing_data
        endif
      endif
C
C     Variable        NETCDF Long Name
C      numRadialsV_HI"number of valid radials in this scan"
C
      nf_status = NF_INQ_VARID(nf_fid,'numRadialsV_HI',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var numRadialsV_HI'
        numRadialsV_HI = i_missing_data
      else
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,numRadialsV_HI)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var numRadialsV_HI'
          numRadialsV_HI = i_missing_data
        endif
      endif
C
C     Variable        NETCDF Long Name
C      timeR        "time since base date"
C
      nf_status = NF_INQ_VARID(nf_fid,'timeR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var timeR'
        timeR = i_missing_data
      else
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,timeR)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var timeR'
          timeR = i_missing_data
        endif
      endif
C
C     Variable        NETCDF Long Name
C      timeR_HI     "time since base date"
C
      nf_status = NF_INQ_VARID(nf_fid,'timeR_HI',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var timeR_HI'
        timeR_HI = i_missing_data
      else
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,timeR_HI)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var timeR_HI'
          timeR_HI = i_missing_data
        endif
      endif
C
C     Variable        NETCDF Long Name
C      timeV        "time since base date"
C
      nf_status = NF_INQ_VARID(nf_fid,'timeV',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var timeV'
        timeV = i_missing_data
      else
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,timeV)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var timeV'
          timeV = i_missing_data
        endif
      endif
C
C     Variable        NETCDF Long Name
C      timeV_HI     "time since base date"
C
      nf_status = NF_INQ_VARID(nf_fid,'timeV_HI',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var timeV_HI'
        timeV_HI = i_missing_data
      else
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,timeV_HI)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var timeV_HI'
          timeV_HI = i_missing_data
        endif
      endif

C   Variables of type DOUBLE
C


C   Variables of type CHAR
C

      nf_status = nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
      endif

      return
      end
