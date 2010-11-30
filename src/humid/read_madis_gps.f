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

      subroutine read_madis_gps (path, filename, time_diff,
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
      character*13 filefound, cvt_i4time_wfo_fname13
      

c     internal
      integer istatus, ptg_index,i4time
      integer file_name_length

      integer recNum, nf_fid, nf_vid, nf_status
      character*120 extension
      integer extension_index
      character*120 desired_ext
      integer de_index


c     prep code
      call s_len(path, ptg_index)
c     create i4time locally from input filename (time) reference variable
      call i4time_fname_lp (filename, i4time, istatus)
c     adjust time to read prior hour (madis data of current hour will NOT
c     contain any data
      i4time = i4time - 3600 ! 3600 = 1 hour, subtraction, one hour earlier
c     create filefound (wfo mode name) from local i4time just generated
      filefound = cvt_i4time_wfo_fname13(i4time)
C
C  Open desired netcdf File for reading
C
      nf_status = NF_OPEN(path(1:ptg_index)//'/'//
     1     filefound,
     1     NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         istatus = 0
         write(6,*) 'failure getting MADIS GPS data'
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
      call read_madis_gps_data (nf_fid , recNum, gps_tpw, gps_error, 
     1     gps_lat, gps_lon,
     1     gps_n,gps_num)

      

      return
      end

 
C
C
      subroutine read_madis_gps_data (nf_fid,recNum,gps_tpw,gps_error, 
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
      integer recNum, nf_fid, nf_vid, nf_status,i,j
      character*80 staLongNam(recNum)
      character*5 staNam(recNum)
      real formalError(recNum), staLat(recNum), staLon(recNum), 
     +   waterVapor(recNum)
      double precision observationTime(recNum),latestTime

      call read_gps_madis_basics (nf_fid , recNum, formalError, staLat,
     +    staLon,  waterVapor, observationTime)
C
C The netcdf variables are filled - your code goes here
C
c     okay for the madis problem we are almost home.  we now have to 
c     make a subtle mod to the old reader to make it look like the old reader
c     in the old reader, the recNum was the total number of data, now
c     that is not the case.  since madis uses recNum differently.

c     but in the old code we did introduce gps_num and now we will actually 
c     use it as intended.

      latestTime = 0.

      do j = 1, recNum
         if(waterVapor(j) .lt. 10000.) then ! good data
            if (latestTime .le. observationTime(j)) then
               latestTime = observationTime(j)
            endif
         endif
      enddo

c     at this point latestTime is the desired time to trap

      i = 0
      do j = 1, recNum
         
         if (waterVapor(j) .lt. 10000. .and.
     +        latestTime .eq. observationTime(j)) then ! presume good data
            i = i +1
            gps_tpw(i) = waterVapor(j)
            gps_error(i) = formalError(j)
            gps_lat(i) = staLat(j)
            gps_lon(i) = staLon(j)
         endif
      enddo

      gps_num = i
      
      return
      end



      subroutine read_gps_madis_basics (nf_fid , recNum, formalError, 
     +     latitude, longitude, totalColumnPWV, observationTime)

      include 'netcdf.inc'
      integer recNum, nf_fid, nf_vid, nf_status



      real formalError(recNum), totalColumnPWV(recNum),
     +     latitude(recNum), longitude(recNum)
      double precision observationTime(recNum)
C


C
C     Variable        NETCDF Long Name
C      formalError  "Formal Error" 
C

c     gather observation time

        nf_status = NF_INQ_VARID(nf_fid,'observationTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var observationTime'
      endif
        nf_status = NF_GET_VAR_Double(nf_fid,nf_vid,observationTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ observationTime '
      endif


c     get formal error


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
        nf_status = NF_INQ_VARID(nf_fid,'latitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var staLat'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,latitude)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ staLat '
      endif
C
C     Variable        NETCDF Long Name
C      staLon       "Station Longitude" 
C
        nf_status = NF_INQ_VARID(nf_fid,'longitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var staLon'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,longitude)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ staLon '
      endif
C
C     Variable        NETCDF Long Name
C      waterVapor   "Water Vapor" 
C
        nf_status = NF_INQ_VARID(nf_fid,'totalColumnPWV',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var waterVapor'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,totalColumnPWV)
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
