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
cdis cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis

      subroutine read_gps (path, filename, time_diff,
     1     gps_tpw, gps_error, gps_lat, 
     1     gps_lon,gps_num,
     1     gps_n, istatus)


      implicit none

      include 'netcdf.inc'

c     parameter list variables
      integer gps_n, gps_num
      real gps_tpw(gps_n)
      real gps_error(gps_n)
      real gps_lat(gps_n)
      real gps_lon(gps_n)
      integer time_diff
      character*256 path 
      character*9 filename
      character*9 filefound
      

c     internal
      integer istatus, ptg_index
      integer file_name_length

      integer recNum, nf_fid, nf_vid, nf_status
      character*120 extension
      integer extension_index
      character*120 desired_ext
      integer de_index


c     prep code
      call s_len(path, ptg_index)

      file_name_length = 14
      desired_ext = 'nc'


!GFORTRAN modifications begin
      extension_index = 0
!GFORTRAN modifications end
      call get_newest_file (filename, time_diff,file_name_length,
     1     path,ptg_index,filefound,desired_ext, de_index,
     1     extension, extension_index, istatus)

!GFORTRAN modifications begin
      if (istatus.ne.1) then    !failure, return with failure code
         write(6,*) 'No GPS data'
         return
      endif
!GFORTRAN modifications end
c     prep filename to fit


C
C  Open netcdf File for reading
C
      nf_status = NF_OPEN(path(1:ptg_index)//'/'//
     1     filefound//'0030o.'//extension(1:extension_index),
     1     NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         istatus = 0
         write(6,*) 'failure getting GPS data'
         return
      else
         istatus = 1
   
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
        istatus = 0
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,recNum)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim recNum'
        istatus = 0
      endif
      call read_gps_data (nf_fid , recNum, gps_tpw, gps_error, 
     1     gps_lat, gps_lon,
     1     gps_n,gps_num)

      return
      end

 
C
C
      subroutine read_gps_data (nf_fid, recNum, gps_tpw, gps_error, 
     1     gps_lat, gps_lon,
     1     gps_n,gps_num)
      include 'netcdf.inc'
c     parameter list variables
      integer gps_n, gps_num
      real gps_tpw(gps_n)
      real gps_error(gps_n)
      real gps_lat(gps_n)
      real gps_lon(gps_n)
c
      integer recNum, nf_fid, nf_vid, nf_status,i
      character*80 staLongNam(recNum)
      character*5 staNam(recNum)
      real formalError(recNum), staLat(recNum), staLon(recNum), 
     +   waterVapor(recNum)
      call read_gps_basics (nf_fid , recNum, formalError, staLat,
     +    staLon, staLongNam, staNam, waterVapor)
C
C The netcdf variables are filled - your code goes here
C

      gps_num = recNum
      do i = 1, recNum
         gps_tpw(i) = waterVapor(i)
         gps_error(i) = formalError(i)
         gps_lat(i) = staLat(i)
         gps_lon(i) = staLon(i)
      enddo

      return
      end
      subroutine read_gps_basics (nf_fid , recNum, formalError, staLat,
     +    staLon, staLongNam, staNam, waterVapor)
      include 'netcdf.inc'
      integer recNum, nf_fid, nf_vid, nf_status

      character*80 staLongNam(recNum)
      character*5 staNam(recNum)
      real formalError(recNum), staLat(recNum), staLon(recNum), 
     +   waterVapor(recNum)
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
        print *,'in NF_GET_VAR_ staLongNam '
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
        print *,'in NF_GET_VAR_ staNam '
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
        print *,'in NF_GET_VAR_ formalError '
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
        print *,'in NF_GET_VAR_ staLat '
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
        print *,'in NF_GET_VAR_ staLon '
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
        print *,'in NF_GET_VAR_ waterVapor '
      endif
      nf_status = nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
      endif

      return
      end
