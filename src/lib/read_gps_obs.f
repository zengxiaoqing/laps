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

      subroutine read_gps_obs (lun, path, i4beg, i4end,
     1     imax, jmax, latgrid, longrid, bad_sfc,
     1     gps_tpw, gps_wet, gps_error, gps_xy, gps_elv, 
     1     gps_tim, gps_indomain, gps_n, istatus)
c     This code is originated from Dan Birkenheuer.

c     On July 16, 2009, Yuanfu Xie modified it to read an additional 
c     variable, wet delays, station elevation and gps obs time.

c     On August 2, 2009, Yuanfu Xie modified it to read all GPS data
c     files between i4begin and i4end.

c     Reworked a bit by Steve Albers in 2011

      implicit none

      include 'netcdf.inc'

c     parameter list variables
      integer max_files
      parameter (max_files=3000)
      integer lun, gps_n, gps_num, gps_indomain, verbose
      integer i4times(max_files),nfiles,gps_i,ifile,i,istat_latlon
      integer imax, jmax
      real latgrid(imax,jmax)
      real longrid(imax,jmax)
      real gps_tpw(gps_n)
      real gps_wet(gps_n)
      real gps_error(gps_n)
      real gps_xy(2,gps_n)
      real gps_elv(gps_n)
      real gps_tim(gps_n)

      real gps_lat(gps_n)
      real gps_lon(gps_n)
      integer i4beg, i4end
      character*256 path,filespec,c_filenames(max_files) 
      
c     internal
      integer istatus, ptg_index
      integer file_name_length

      integer recNum, nf_fid, nf_vid, nf_status
      real x, y, bad_sfc

      verbose = 0

c     Set file specs for get_file_times:
      call s_len(path, ptg_index)
      filespec(1:ptg_index) = path(1:ptg_index)
      filespec(ptg_index+1:ptg_index+9) = '*0030o.nc'  ! Hardcode for now by Yuanfu

c     Get all filenames under the path:
      call get_file_times(filespec(1:ptg_index+9),max_files,c_filenames,
     1                    i4times,nfiles,istatus)

      gps_i = 0    ! Total GPS data read
c     Loop through all available files: Yuanfu
      do ifile=1,nfiles
      
C
C       Open netcdf File for reading
C
        ! Read data file between i4beg and i4end:
        if (i4times(ifile) .ge. i4beg .and. 
     1      i4times(ifile) .le. i4end) then
        nf_status = NF_OPEN(c_filenames(ifile),NF_NOWRITE,nf_fid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          istatus = 0
          write(6,*) 'failure getting GPS data'
          nf_status = NF_CLOSE(nf_fid)
          return
        else
          istatus = 1
          write(6,*)' Opened gps file: ',trim(c_filenames(ifile))
        endif
C
C       Fill all dimension values
C
C
C       Get size of recNum
C
        nf_status = NF_INQ_DIMID(nf_fid,'recNum',nf_vid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'dim recNum'
          istatus = 0
          goto 900  
        endif
        nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,recNum)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'dim recNum'
          istatus = 0
          goto 900  
        endif

        ! Check if recNum is larger than space allocated:
        if (gps_i+recNum .gt. gps_n) then
          print *,' Too many GPS obs, increase your gps_n and rerun!'
     1	,gps_i+recNum,' > ',gps_n
          stop
        endif

        call read_gps_data (nf_fid , recNum, gps_tpw(gps_i+1), 
     1       gps_wet(gps_i+1), gps_error(gps_i+1), gps_lat(gps_i+1), 
     1       gps_lon(gps_i+1), gps_elv(gps_i+1), gps_tim(gps_i+1), 
     1       gps_n,gps_num)

c       Accumulate all:
        gps_i = gps_i+gps_num
        write(6,*)' gps_num / gps_i = ',gps_num,gps_i

        ! Finish read this qualified file:
        endif

 900  enddo
      ! Check if there is any gps obs:
      if (gps_i .ge. 1) istatus = 1
      if (gps_i .le. 0) istatus = 0

      ! Write the GPS wetdelay information into LAPS HMG file:
      gps_indomain = 0
      do i=1,gps_i ! loop through all obs from the files
        CALL LATLON_TO_RLAPSGRID(gps_lat(i),gps_lon(i),latgrid,
     1             longrid, imax, jmax, x, y, istat_latlon)

        if(verbose .ge. 1 .or. i .le. 10)then
          write(6,101)i,gps_lat(i),gps_lon(i),x,y,gps_wet(i)
 101      format(i7,4f9.2,e16.6)
        endif

        if (x .ge. 1 .and. x .le. imax .and.
     1      y .ge. 1 .and. y .le. jmax .and. 
     1    gps_wet(i) .ne. bad_sfc) then
          gps_indomain = gps_indomain + 1
          gps_tpw(gps_indomain) = gps_tpw(i)
          gps_wet(gps_indomain) = gps_wet(i)
          gps_error(gps_indomain) = gps_error(i)
          gps_elv(gps_indomain) = gps_elv(i)
          gps_tim(gps_indomain) = gps_tim(i)
          gps_xy(1,gps_indomain) = x
          gps_xy(2,gps_indomain) = y

          if(lun .gt. 0)then! Write wet delays and tpw into hmg file:
            write(lun, *) x, y,0, gps_wet(i), 'GPSWET'
            write(lun, *) x, y,0, gps_tpw(i), 'GPSTPW'
          endif

          if(verbose .ge. 1)then
            write(6, *) x, y,0, gps_wet(i), 'GPSWET'
            write(6, *) x, y,0, gps_tpw(i), 'GPSTPW'
          endif

        endif
      enddo

      return

      end

C
C
      subroutine read_gps_data (nf_fid, recNum, gps_tpw, 
     1     gps_wet, gps_error, gps_lat, gps_lon, gps_elv,
     1     gps_tim, gps_n,gps_num)
      include 'netcdf.inc'

c     This code is originated from Dan Birkenheuer.
c     On July 16, 2009, Yuanfu Xie modified it to read an additional
c     variable, wet delays, station elevation and gps obs time.

c     parameter list variables
      integer gps_n, gps_num
      real gps_tpw(gps_n)
      real gps_wet(gps_n)
      real gps_error(gps_n)
      real gps_lat(gps_n)
      real gps_lon(gps_n)
      real gps_elv(gps_n)
      real gps_tim(gps_n)
c
      integer recNum, nf_fid, nf_vid, nf_status,i
      character*80 staLongNam(recNum)
      character*5 staNam(recNum)
      real formalError(recNum), staLat(recNum), staLon(recNum), 
     +     staElv(recNum), timObs(recNum), waterVapor(recNum), 
     +     wetDelay(recNum)
      call read_gps_basics (nf_fid , recNum, formalError, staLat,
     +    staLon, staElv, timObs, staLongNam, staNam, waterVapor, 
     +    wetDelay) 
      ! Yuanfu: wetDelay, staElv -- station elevation, and gps obs time
C
C The netcdf variables are filled - your code goes here
C

      gps_num = recNum
      do i = 1, recNum
         gps_tpw(i) = waterVapor(i)
         gps_wet(i) = wetDelay(i)
         gps_error(i) = formalError(i)
         gps_lat(i) = staLat(i)
         gps_lon(i) = staLon(i)
         gps_elv(i) = staElv(i)
         gps_tim(i) = timObs(i)
      enddo

      return
      end
      subroutine read_gps_basics (nf_fid , recNum, formalError, 
     +           staLat,staLon, staElv, timObs, staLongNam, 
     +           staNam, waterVapor,wetDelay)	
      ! Yuanfu: add wetDelay and staElv -- station elevation

c     This code is originated from Dan Birkenheuer.
c     On July 16, 2009, Yuanfu Xie modified it to read an additional
c     variable, wet delays and station elevation.

      include 'netcdf.inc'
      integer recNum, nf_fid, nf_vid, nf_status

      character*80 staLongNam(recNum)
      character*5 staNam(recNum)
      real formalError(recNum), staLat(recNum), staLon(recNum), 
     +     staElv(recNum), timObs(recNum), waterVapor(recNum), 
     +     wetDelay(recNum)	
      ! Yuanfu: add wetDelay and staElv
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
C      staElv       "Station Elevation" 
C
        nf_status = NF_INQ_VARID(nf_fid,'staElev',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var staElev'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,staElv)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ staElev '
      endif
C
C     Variable        NETCDF Long Name
C      timObs       "Time of observation" 
C
        nf_status = NF_INQ_VARID(nf_fid,'timeObs',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var timeObs'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,timObs)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ timeObs '
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
C
C     Variable        NETCDF Long Name
C      wetDelay   "Wet component GPS signal delay" 
C
        nf_status = NF_INQ_VARID(nf_fid,'wetDelay',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var wetDelay'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,wetDelay)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ wetDelay '
      endif


C     Close file:
      nf_status = nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
      endif

      return
      end


      subroutine get_gps_path(path_to_gps_out,istatus)

      use mem_namelist, ONLY: path_to_gps

      character*256 path_to_gps_out

      logical l_exist

      path_to_gps_out = path_to_gps 

      inquire(file=trim(path_to_gps_out),exist=l_exist)
      if(l_exist .eqv. .true.)then
          write(6,*)' Path to GPS from namelist exists'         
      else
          write(6,*)
     1              ' Path to GPS non-existent, trying /public'
          path_to_gps_out = '/public/data/gpsmet/netcdf/'
      endif

      istatus = 1

      return
      end
