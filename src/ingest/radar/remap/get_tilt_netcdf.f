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

      subroutine get_tilt_netcdf_hdr(filename,nf_fid
     1                               ,radarName
     1                               ,siteLat                        
     1                               ,siteLon                        
     1                               ,siteAlt                        
     1                               ,elevationAngle
     1                               ,numRadials 
     1                               ,elevationNumber
     1                               ,VCP
     1                               ,radialAzim
     1                               ,resolutionV
     1                               ,gateSizeV,gateSizeZ
     1                               ,firstGateRangeV,firstGateRangeZ
     1                               ,V_bin_max, Z_bin_max, radial_max ! I
     1                               ,V_bin,     Z_bin,     radial     ! O
     1                               ,istatus)

!     Argument List
      character*(*) filename
      character*5  radarName
      integer V_bin_max, Z_bin_max, radial_max
      real radialAzim(radial_max)

!     Local
      real radialElev(radial_max)
      character*132 siteName
      double precision esEndTime, esStartTime, radialTime(radial_max)

!.............................................................................

      include 'netcdf.inc'
      integer V_bin, Z_bin, radial,nf_fid, nf_vid, nf_status

!     include 'remap_constants.dat' ! for debugging only
!     include 'remap.cmn' ! for debugging only

      write(6,*)' get_tilt_netcdf_hdr: reading ',filename
C
C  Open netcdf File for reading
C
      nf_status = NF_OPEN(filename,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'NF_OPEN ',filename
        istatus = 0
        return
      endif
C
C  Fill all dimension values
C
C
C Get size of V_bin
C
      nf_status = NF_INQ_DIMID(nf_fid,'V_bin',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim V_bin'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,V_bin)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim V_bin'
      endif
C
C Get size of Z_bin
C
      nf_status = NF_INQ_DIMID(nf_fid,'Z_bin',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim Z_bin'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,Z_bin)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim Z_bin'
      endif
C
C Get size of radial
C
      nf_status = NF_INQ_DIMID(nf_fid,'radial',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim radial'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,radial)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim radial'
      endif

!.....Test whether dimensions of NetCDF file are within bounds...............

      if(V_bin .gt. V_bin_max)then
          write(6,*)' V_bin > permitted dimensions ',V_bin,V_bin_max
          stop
      endif

      if(Z_bin .gt. Z_bin_max)then
          write(6,*)' Z_bin > permitted dimensions ',Z_bin,Z_bin_max     
          stop
      endif

      if(radial .gt. radial_max)then
          write(6,*)' radial > permitted dimensions ',radial,radial_max       
          stop
      endif

      istatus = 1
      return

      end

      subroutine get_tilt_netcdf_data(filename,nf_fid
     1                               ,radarName
     1                               ,siteLat                        
     1                               ,siteLon                        
     1                               ,siteAlt                        
     1                               ,elevationAngle
     1                               ,numRadials 
     1                               ,elevationNumber
     1                               ,VCP
     1                               ,nyquist
     1                               ,radialAzim
     1                               ,Z  
     1                               ,V
     1                               ,resolutionV
     1                               ,gateSizeV,gateSizeZ
     1                               ,firstGateRangeV,firstGateRangeZ
     1                               ,Z_scale, Z_offset
     1                               ,V_scale, V_offset
     1                               ,V_bin_in, Z_bin_in, radial_in    ! I
     1                               ,istatus)

!     Argument List
      character*(*) filename
      character*5  radarName
      integer V_bin_in, Z_bin_in, radial_in
      integer V(V_bin_in,radial_in), Z(Z_bin_in,radial_in)
      real radialAzim(radial_in)

!     Local
      real radialElev(radial_in)
      character*132 siteName
      double precision esEndTime, esStartTime, radialTime(radial_in)

!.............................................................................

      include 'netcdf.inc'
      integer V_bin, Z_bin, radial,nf_fid, nf_vid, nf_status

!     include 'remap_constants.dat' ! for debugging only
!     include 'remap.cmn' ! for debugging only

      write(6,*)' get_tilt_netcdf_data: reading ',filename
      write(6,*)'                       v/z/rad = '
     1         ,v_bin_in,z_bin_in,radial_in
C
C  Open netcdf File for reading
C
!     nf_status = NF_OPEN(filename,NF_NOWRITE,nf_fid)
!     if(nf_status.ne.NF_NOERR) then
!       print *, NF_STRERROR(nf_status)
!       print *,'NF_OPEN ',filename
!       istatus = 0
!       return
!     endif
C
C  Fill all dimension values
C
C
C Get size of V_bin
C
      nf_status = NF_INQ_DIMID(nf_fid,'V_bin',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim V_bin'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,V_bin)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim V_bin'
      endif
C
C Get size of Z_bin
C
      nf_status = NF_INQ_DIMID(nf_fid,'Z_bin',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim Z_bin'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,Z_bin)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim Z_bin'
      endif
C
C Get size of radial
C
      nf_status = NF_INQ_DIMID(nf_fid,'radial',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim radial'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,radial)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim radial'
      endif

!.....Test whether dimensions of NetCDF file are within bounds...............

      if(V_bin .ne. V_bin_in)then
          write(6,*)' V_bin != permitted dimensions ',V_bin,V_bin_in
          stop
      endif

      if(Z_bin .ne. Z_bin_in)then
          write(6,*)' Z_bin != permitted dimensions ',Z_bin,Z_bin_in     
          stop
      endif

      if(radial .ne. radial_in)then
          write(6,*)' radial != permitted dimensions ',radial,radial_in       
          stop
      endif

      call read_netcdf(nf_fid, V_bin_in, Z_bin_in, radial_in,        ! I
!............................................................................
     +     V, VCP, 
     +     Z, elevationNumber, numGatesV, numGatesZ, numRadials, 
     +     atmosAttenFactor, calibConst, elevationAngle, 
     +     firstGateRangeV, firstGateRangeZ, gateSizeV, gateSizeZ, 
     +     Z_scale, Z_offset, V_scale, V_offset,
     +     nyquist, powDiffThreshold, radialAzim, radialElev, 
     +     resolutionV, siteAlt, siteLat, siteLon, unambigRange, 
     +     esEndTime, esStartTime, radialTime, radarName, siteName)

      if(sitealt .eq. 0. .or. sitelat .eq. 0. .or. sitelon .eq. 0.)then       
          write(6,*)' Warning, no site info in get_tilt_netcdf_data'       
          istatus = 0
          return
      else
          write(6,*)' Site info:',siteAlt, siteLat, siteLon       
          istatus = 1
          return
      endif

      end
C
C
C
C  Subroutine to read the file "WSR-88D Wideband Data" 
C
      subroutine read_netcdf(nf_fid, V_bin, Z_bin, radial, V, VCP, 
     +     Z, elevationNumber, numGatesV, numGatesZ, numRadials, 
     +     atmosAttenFactor, calibConst, elevationAngle, 
     +     firstGateRangeV, firstGateRangeZ, gateSizeV, gateSizeZ, 
     +     Z_scale, Z_offset, V_scale, V_offset,
     +     nyquist, powDiffThreshold, radialAzim, radialElev, 
     +     resolutionV, siteAlt, siteLat, siteLon, unambigRange, 
     +     esEndTime, esStartTime, radialTime, radarName, siteName)
C
      include 'netcdf.inc'
!     include 'remap_constants.dat' ! for debugging only
!     include 'remap.cmn' ! for debugging only
      integer V_bin, Z_bin, radial,nf_fid, nf_vid, nf_status
      integer V( V_bin, radial), VCP, Z( Z_bin,
     +     radial), elevationNumber, numGatesV, numGatesZ, numRadials
      real atmosAttenFactor, calibConst, elevationAngle,
     +     firstGateRangeV, firstGateRangeZ, gateSizeV, gateSizeZ,
     +     nyquist, powDiffThreshold, radialAzim(radial),
     +     radialElev(radial), resolutionV, siteAlt, siteLat,
     +     siteLon, unambigRange
      double precision esEndTime, esStartTime, radialTime(radial)
      character*5 radarName
      character*132 siteName


C   Variables of type REAL
C
C     Variable        NETCDF Long Name
C      atmosAttenFactor"Atmospheric attenuation factor"
C
        nf_status = NF_INQ_VARID(nf_fid,'atmosAttenFactor',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var atmosAttenFactor'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,atmosAttenFactor)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var atmosAttenFactor'
      endif
C
C     Variable        NETCDF Long Name
C      calibConst   "System gain calibration constant"
C
        nf_status = NF_INQ_VARID(nf_fid,'calibConst',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var calibConst'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,calibConst)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var calibConst'
      endif
C
C     Variable        NETCDF Long Name
C      elevationAngle"Elevation angle"
C
        nf_status = NF_INQ_VARID(nf_fid,'elevationAngle',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var elevationAngle'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,elevationAngle)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var elevationAngle'
      endif
C
C     Variable        NETCDF Long Name
C      firstGateRangeV"Range to 1st Doppler gate"
C
        nf_status = NF_INQ_VARID(nf_fid,'firstGateRangeV',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var firstGateRangeV'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,firstGateRangeV)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var firstGateRangeV'
      endif
C
C     Variable        NETCDF Long Name
C      firstGateRangeZ"Range to 1st Reflectivity gate"
C
        nf_status = NF_INQ_VARID(nf_fid,'firstGateRangeZ',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var firstGateRangeZ'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,firstGateRangeZ)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var firstGateRangeZ'
      endif
C
C     Variable        NETCDF Long Name
C      gateSizeV    "Doppler gate spacing"
C
        nf_status = NF_INQ_VARID(nf_fid,'gateSizeV',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var gateSizeV'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,gateSizeV)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var gateSizeV'
      endif
C
C     Variable        NETCDF Long Name
C      gateSizeZ    "Reflectivity gate spacing"
C
        nf_status = NF_INQ_VARID(nf_fid,'gateSizeZ',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var gateSizeZ'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,gateSizeZ)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var gateSizeZ'
      endif
C
C     Variable        NETCDF Long Name
C      Z_scale  "Reflectivity scale value"
C
      nf_status = NF_INQ_VARID(nf_fid,'Z_scale',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Z_scale'
      else
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,Z_scale)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var Z_scale'
        endif
      endif
C
C     Variable        NETCDF Long Name
C      Z_offset  "Reflectivity offset value"
C
      nf_status = NF_INQ_VARID(nf_fid,'Z_offset',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Z_offset'
      else
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,Z_offset)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var Z_offset'
        endif
      endif
C
C     Variable        NETCDF Long Name
C      nyquist      "Nyquist velocity"
C
        nf_status = NF_INQ_VARID(nf_fid,'nyquist',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var nyquist'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,nyquist)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var nyquist'
      endif
C
C     Variable        NETCDF Long Name
C      powDiffThreshold"Range de-aliasing threshold"
C
        nf_status = NF_INQ_VARID(nf_fid,'powDiffThreshold',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var powDiffThreshold'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,powDiffThreshold)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var powDiffThreshold'
      endif
C
C     Variable        NETCDF Long Name
C      radialAzim   "Radial azimuth angle"
C
        nf_status = NF_INQ_VARID(nf_fid,'radialAzim',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var radialAzim'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,radialAzim)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var radialAzim'
      endif
C
C     Variable        NETCDF Long Name
C      radialElev   "Radial elevation angle"
C
        nf_status = NF_INQ_VARID(nf_fid,'radialElev',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var radialElev'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,radialElev)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var radialElev'
      endif
C
C     Variable        NETCDF Long Name
C      resolutionV  "Doppler velocity resolution"
C
        nf_status = NF_INQ_VARID(nf_fid,'resolutionV',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var resolutionV'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,resolutionV)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var resolutionV'
      endif
C
C     Variable        NETCDF Long Name
C      V_scale  "Velocity scale value"
C
      nf_status = NF_INQ_VARID(nf_fid,'V_scale',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var V_scale'
      else
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,V_scale)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var V_scale'
        endif
      endif
C
C     Variable        NETCDF Long Name
C      V_offset  "Velocity offset value"
C
      nf_status = NF_INQ_VARID(nf_fid,'V_offset',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var V_offset'
      else
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,V_offset)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var V_offset'
        endif
      endif
C
C     Variable        NETCDF Long Name
C      siteAlt      "Altitude of site above mean sea level"
C
        nf_status = NF_INQ_VARID(nf_fid,'siteAlt',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var siteAlt'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,siteAlt)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var siteAlt'
      endif
C
C     Variable        NETCDF Long Name
C      siteLat      "Latitude of site"
C
        nf_status = NF_INQ_VARID(nf_fid,'siteLat',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var siteLat'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,siteLat)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var siteLat'
      endif
C
C     Variable        NETCDF Long Name
C      siteLon      "Longitude of site"
C
        nf_status = NF_INQ_VARID(nf_fid,'siteLon',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var siteLon'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,siteLon)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var siteLon'
      endif
C
C     Variable        NETCDF Long Name
C      unambigRange "Unambiguous range"
C
        nf_status = NF_INQ_VARID(nf_fid,'unambigRange',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var unambigRange'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,unambigRange)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var unambigRange'
      endif

C   Variables of type INT
C
C
C     Variable        NETCDF Long Name
C      V            "Velocity"
C
        nf_status = NF_INQ_VARID(nf_fid,'V',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var V'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,V)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var V'
      endif
C
C     Variable        NETCDF Long Name
C      VCP          "Volume Coverage Pattern"
C
        nf_status = NF_INQ_VARID(nf_fid,'VCP',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var VCP'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,VCP)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var VCP'
      endif
C
C     Variable        NETCDF Long Name
C      W            "Spectrum Width"
C
!        nf_status = NF_INQ_VARID(nf_fid,'W',nf_vid)
!      if(nf_status.ne.NF_NOERR) then
!        print *, NF_STRERROR(nf_status)
!        print *,'in var W'
!      endif
!        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,W)
!      if(nf_status.ne.NF_NOERR) then
!        print *, NF_STRERROR(nf_status)
!        print *,'in var W'
!      endif
C
C     Variable        NETCDF Long Name
C      Z            "Reflectivity"
C
        nf_status = NF_INQ_VARID(nf_fid,'Z',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Z'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,Z)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Z'
      endif
C
C     Variable        NETCDF Long Name
C      elevationNumber"Elevation number"
C
        nf_status = NF_INQ_VARID(nf_fid,'elevationNumber',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var elevationNumber'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,elevationNumber)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var elevationNumber'
      endif
C
C     Variable        NETCDF Long Name
C      numGatesV    "Number of Doppler gates"
C
        nf_status = NF_INQ_VARID(nf_fid,'numGatesV',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var numGatesV'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,numGatesV)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var numGatesV'
      endif
C
C     Variable        NETCDF Long Name
C      numGatesZ    "Number of reflectivity gates"
C
        nf_status = NF_INQ_VARID(nf_fid,'numGatesZ',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var numGatesZ'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,numGatesZ)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var numGatesZ'
      endif
C
C     Variable        NETCDF Long Name
C      numRadials   "Number of radials"
C
        nf_status = NF_INQ_VARID(nf_fid,'numRadials',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var numRadials'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,numRadials)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var numRadials'
      endif

C   Variables of type DOUBLE
C
C
C     Variable        NETCDF Long Name
C      esEndTime    "End time of elevation scan"
C
        nf_status = NF_INQ_VARID(nf_fid,'esEndTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var esEndTime'
      endif
        nf_status = NF_GET_VAR_DOUBLE(nf_fid,nf_vid,esEndTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var esEndTime'
      endif
C
C     Variable        NETCDF Long Name
C      esStartTime  "Start time of elevation scan"
C
        nf_status = NF_INQ_VARID(nf_fid,'esStartTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var esStartTime'
      endif
        nf_status = NF_GET_VAR_DOUBLE(nf_fid,nf_vid,esStartTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var esStartTime'
      endif
C
C     Variable        NETCDF Long Name
C      radialTime   "Time of radial"
C
        nf_status = NF_INQ_VARID(nf_fid,'radialTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var radialTime'
      endif
        nf_status = NF_GET_VAR_DOUBLE(nf_fid,nf_vid,radialTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var radialTime'
      endif


C   Variables of type CHAR
C
C
C     Variable        NETCDF Long Name
C      radarName    "Official name of the radar"
C
        nf_status = NF_INQ_VARID(nf_fid,'radarName',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var radarName'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,radarName)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var radarName'
      endif
C
C     Variable        NETCDF Long Name
C      siteName     "Long name of the radar site"
C
        nf_status = NF_INQ_VARID(nf_fid,'siteName',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var siteName'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,siteName)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var siteName'
      endif

      nf_status = nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
      endif

      return
      end
