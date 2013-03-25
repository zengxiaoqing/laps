!dis
!dis    Open Source License/Disclaimer, Forecast Systems Laboratory
!dis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
!dis
!dis    This software is distributed under the Open Source Definition,
!dis    which may be found at http://www.opensource.org/osd.html.
!dis
!dis    In particular, redistribution and use in source and binary forms,
!dis    with or without modification, are permitted provided that the
!dis    following conditions are met:
!dis
!dis    - Redistributions of source code must retain this notice, this
!dis    list of conditions and the following disclaimer.
!dis
!dis    - Redistributions in binary form must provide access to this
!dis    notice, this list of conditions and the following disclaimer, and
!dis    the underlying source code.
!dis
!dis    - All modifications to this software must be clearly documented,
!dis    and are solely the responsibility of the agent making the
!dis    modifications.
!dis
!dis    - If significant modifications or enhancements are made to this
!dis    software, the FSL Software Policy Manager
!dis    (softwaremgr@fsl.noaa.gov) should be notified.
!dis
!dis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
!dis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
!dis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
!dis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
!dis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
!dis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
!dis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
!dis
!dis
!  Module to handle I/O from a WRF v2 netCDF file
!

MODULE wrf_netcdf

! F90 module to deal with reading WRF output files (WRFv1 netcdf)

  USE grid_utils
  include 'netcdf.inc'
  PRIVATE
    
  PUBLIC open_wrfnc
  PUBLIC close_wrfnc
  PUBLIC get_wrfsi_map
  PUBLIC get_wrf2_map
  PUBLIC get_wrf2_timeinfo
  PUBLIC get_wrf_misc
  PUBLIC get_wrf2_misc
  PUBLIC get_wrf_scalar
  PUBLIC get_wrf_1d
  PUBLIC get_wrfnc_2d
  PUBLIC get_wrfnc_3d
  PUBLIC make_wrf_file_name
  PUBLIC make_wrf2_file_name
  PUBLIC wrfio_wait

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE open_wrfnc(fname, lun, status)

    ! Opens a WRF netCDF file, returning the CDF id of the file

    IMPLICIT NONE
    include 'netcdf.inc'  
    ! Arguments

    CHARACTER(LEN=255),INTENT(IN)          :: fname
    INTEGER, INTENT(OUT)                   :: lun
    INTEGER, INTENT(OUT)                   :: status

    status = 0
    lun = NCOPN(fname, NCNOWRIT,   status)

    IF (status .ne. 0) then
       print *, 'Error opening netCDF file: ',TRIM(fname)
    ENDIF

    RETURN
  END SUBROUTINE open_wrfnc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE close_wrfnc(cdfid)

    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: cdfid
    INTEGER :: status 
    include 'netcdf.inc'

    status = nf_close(cdfid)
    if (status .ne. nf_noerr) THEN
      print *, 'Problem closing the netCDF file.'
    endif
    return
 END SUBROUTINE close_wrfnc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_wrfsi_map(lun, projection, lat1, lon1, stdlon, &
                         truelat1, truelat2, dx_m, dy_m, &
                         nx, ny, status)
  
    ! Reads from a WRFSI static file to get mapping information for a WRF domain
    ! NOTE:  The parameters returned are the non-staggered grid parameters,
    ! as defined in the wrfsi.nl file.  The dimensions found in the wrfout
    ! netCDF files require have one less element in each direction.

    ! Assumes that the static file has been opened with the open_wrfnc routine
    ! and the CDF ID is the "lun". 

    IMPLICIT NONE
    include 'netcdf.inc' 
    ! Args
    INTEGER, INTENT(IN)                    :: lun  
    CHARACTER (LEN=32), INTENT(OUT)        :: projection
    REAL, INTENT(OUT)                      :: lat1
    REAL, INTENT(OUT)                      :: lon1 
    REAL, INTENT(OUT)                      :: stdlon
    REAL, INTENT(OUT)                      :: truelat1
    REAL, INTENT(OUT)                      :: truelat2
    REAL, INTENT(OUT)                      :: dx_m
    REAL, INTENT(OUT)                      :: dy_m
    INTEGER, INTENT(OUT)                   :: nx
    INTEGER, INTENT(OUT)                   :: ny
    INTEGER, INTENT(OUT)                   :: status
  
    !Locals
    INTEGER :: vid, rcode
    CHARACTER(LEN=132) :: dum
    CHARACTER(LEN=132) :: grid_type
    INTEGER, DIMENSION(2) :: startc, countc
    INTEGER, DIMENSION(4) :: start, count
    status = 0

    ! Get x-y dimensions

    vid = NCVID(lun, 'x', rcode)
    CALL NCDINQ(lun,vid,dum,nx,rcode)
    IF (rcode .NE. 0) THEN
       PRINT *, 'Error getting nx.'
       status = 1
    ENDIF

    vid = NCVID(lun, 'y', rcode)
    CALL NCDINQ(lun,vid,dum,ny,rcode)
    IF (rcode .NE. 0) THEN
       PRINT *, 'Error getting nx.'
       status = 1
    ENDIF

    !  Get projection
    startc = (/1,1/)
    countc = (/132,1/)
    vid = NCVID(lun, 'grid_type', rcode)
    CALL NCVGTC(lun,vid,startc,countc,grid_type, 132, rcode)

    SELECT CASE(grid_type)
      CASE ('mercator') 
        projection = 'MERCATOR                        '
      CASE ('secant lambert conformal')
        projection = 'LAMBERT CONFORMAL               '
      CASE ('tangential lambert conformal')
        projection = 'LAMBERT CONFORMAL               '
      CASE ('polar stereographic')
        projection = 'POLAR STEREOGRAPHIC             '
      CASE DEFAULT
        print *, 'Unrecognized grid type: ', grid_type
        status = 1
    END SELECT

    ! Get SW corner point
    vid = NCVID(lun, 'lat', rcode)
    CALL NCVGT1(lun, vid, 1, lat1, rcode)
    vid = NCVID(lun, 'lon', rcode)
    CALL NCVGT1(lun, vid, 1, lon1, rcode) 

    ! Get Truelat1/trulat2
    vid = NCVID(lun, 'Latin1', rcode)
    CALL NCVGT1(lun, vid, 1, truelat1, rcode)
    vid = NCVID(lun, 'Latin2', rcode)
    CALL NCVGT1(lun, vid, 1, truelat2, rcode)

    ! Get standard longitude
    vid = NCVID(lun, 'LoV', rcode)
    CALL NCVGT1(lun, vid, 1, stdlon, rcode)

    ! Get dx/dy
    vid = NCVID(lun, 'Dx', rcode)
    CALL NCVGT1(lun, vid, 1, dx_m, rcode)
    vid = NCVID(lun, 'Dy', rcode)
    CALL NCVGT1(lun, vid, 1, dy_m, rcode)

    ! Clean up
    if (lon1 .gt. 180.) lon1 = lon1 - 360.
    if (stdlon .gt. 180.) stdlon = stdlon - 360.
 
    RETURN

  END SUBROUTINE get_wrfsi_map 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_wrf2_timeinfo(lun,reftime,dt,itimestep,tau_hr,tau_min,tau_sec, &
                               status)

    IMPLICIT NONE
    include 'netcdf.inc'

    INTEGER, INTENT(IN)                  :: lun
    CHARACTER(LEN=19), INTENT(OUT)       :: reftime
    REAL, INTENT(OUT)                    :: dt
    REAL                                 :: xtime
    INTEGER, INTENT(OUT)                 :: itimestep
    INTEGER, INTENT(OUT)                 :: tau_hr, tau_min, tau_sec, status

    INTEGER   :: rcode, vid

    reftime = '????-??-??_??:??:??' 
    tau_hr = 0
    tau_min = 0
    tau_sec = 0
    status = 0

    rcode = NF_GET_ATT_TEXT(lun,0,"SIMULATION_START_DATE",reftime)
    rcode = NF_GET_ATT_REAL(lun,0,"DT",dt)
    rcode = NF_INQ_VARID(lun,"ITIMESTEP",vid)
    rcode = NF_GET_VAR_INT(lun,vid,itimestep)

    !!!!! Modified by Wei-Ting (20130307) :            !!!!!
    !!!!!     Use XTIME instead of ITIMESTEP FOR WRFV3 !!!!!
    if (rcode .lt. 0) then
      write(6,*)' Warning: ITIMESTEP not found, looking for XTIME in get_wrf2_timeinfo'
      rcode = NF_INQ_VARID(lun,"XTIME",vid)
      rcode = NF_GET_VAR_REAL(lun,vid,xtime)
      itimestep = INT(xtime*60/dt)
    endif
    !!!!! End of Modifying !!!!!
    
    tau_hr = INT(FLOAT(itimestep)*dt)/3600
    tau_sec = MOD(INT(FLOAT(itimestep)*dt),3600)
    tau_min = tau_sec / 60
    tau_sec = MOD(tau_sec,60) 

    RETURN
END SUBROUTINE get_wrf2_timeinfo
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_wrf2_map(lun, stag, projcode, lat1, lon1, stdlon, &
                          truelat1, truelat2, dx_m, dy_m, &
                          nx, ny, nz, status)

    IMPLICIT NONE
    include 'netcdf.inc'
    ! Args
    INTEGER, INTENT(IN)                    :: lun
    CHARACTER (LEN=1)                      :: stag
    INTEGER, INTENT(OUT)                   :: projcode  
    REAL, INTENT(OUT)                      :: lat1
    REAL, INTENT(OUT)                      :: lon1
    REAL, INTENT(OUT)                      :: stdlon
    REAL, INTENT(OUT)                      :: truelat1
    REAL, INTENT(OUT)                      :: truelat2
    REAL, INTENT(OUT)                      :: dx_m
    REAL, INTENT(OUT)                      :: dy_m
    INTEGER, INTENT(OUT)                   :: nx
    INTEGER, INTENT(OUT)                   :: ny
    INTEGER, INTENT(OUT)                   :: nz
    INTEGER, INTENT(OUT)                   :: status

  
    ! Local variables
    INTEGER   :: rcode,dimid,vid
    INTEGER   :: projcode_wrf
    INTEGER   :: wrf_version
    real, allocatable, dimension(:,:) :: dum_2d

    ! Initialization
    lat1 = -999.
    lon1 = -999.
    projcode = 0     
    projcode_wrf = 0
    stdlon = -999.
    truelat1 = -999.
    truelat2 = -999.
    dx_m  = -999.
    dy_m  = -999.
    nx = 0
    ny = 0
    nz =0
    status = 0

    rcode = NF_GET_ATT_INT(lun, 0, "MAP_PROJ", projcode_wrf)
    SELECT CASE (projcode_wrf)
      CASE(1)
        projcode = 3
      CASE(2)
        projcode = 5
      CASE(3)
        projcode = 1
      CASE DEFAULT
        projcode = 99
    END SELECT 
    rcode = NF_GET_ATT_REAL(lun, 0, "STAND_LON",stdlon)
    rcode = NF_GET_ATT_REAL(lun, 0, "TRUELAT1",truelat1)
    rcode = NF_GET_ATT_REAL(lun, 0, "TRUELAT2",truelat2)
    rcode = NF_GET_ATT_REAL(lun, 0, "DX", dx_m)
    rcode = NF_GET_ATT_REAL(lun, 0, "DY", dy_m)
    
    ! Dimenions
    rcode = NF_INQ_DIMID(lun, "west_east", dimid)
    rcode = NF_INQ_DIMLEN(lun, dimid, nx)
    rcode = NF_INQ_DIMID(lun, "south_north", dimid)
    rcode = NF_INQ_DIMLEN(lun, dimid, ny)

!   Allocate 
    ALLOCATE(dum_2d (nx,ny))

    SELECT CASE (stag)
      CASE('T')
        rcode = NF_INQ_VARID(lun,"LAT_LL_T",vid)
        if(rcode .eq. -49)then ! assume version 3 WRF
          write(6,*)' Warning: LAT_LL_T not found, looking for XLAT'
          rcode = NF_INQ_VARID(lun,"XLAT",vid)
          rcode = NF_GET_VAR_REAL(lun,vid,dum_2d)
          lat1 = dum_2d(1,1)
          wrf_version = 3
        else
          rcode = NF_GET_VAR_REAL(lun,vid,lat1)
          wrf_version = 2
        endif

        rcode = NF_INQ_VARID(lun,"LON_LL_T",vid)
        if(rcode .eq. -49)then ! assume version 3 WRF
          write(6,*)' Warning: LON_LL_T not found, looking for XLONG'
          rcode = NF_INQ_VARID(lun,"XLONG",vid)
          rcode = NF_GET_VAR_REAL(lun,vid,dum_2d)
          lon1 = dum_2d(1,1)
        else
          rcode = NF_GET_VAR_REAL(lun,vid,lon1)
        endif

      CASE('U')
        rcode = NF_INQ_VARID(lun,"LAT_LL_U",vid)
        if(rcode .eq. -49)then ! assume version 3 WRF
          write(6,*)' Warning: LAT_LL_U not found, looking for XLAT_U'
          rcode = NF_INQ_VARID(lun,"XLAT_U",vid)
          rcode = NF_GET_VAR_REAL(lun,vid,dum_2d)
          lat1 = dum_2d(1,1)
          wrf_version = 3
        else
          rcode = NF_GET_VAR_REAL(lun,vid,lat1)
          wrf_version = 2
        endif

        rcode = NF_INQ_VARID(lun,"LON_LL_U",vid)
        if(rcode .eq. -49)then ! assume version 3 WRF
          write(6,*)' Warning: LON_LL_U not found, looking for XLONG_U'
          rcode = NF_INQ_VARID(lun,"XLONG_U",vid)
          rcode = NF_GET_VAR_REAL(lun,vid,dum_2d)
          lon1 = dum_2d(1,1)
        else
          rcode = NF_GET_VAR_REAL(lun,vid,lon1)
        endif
        nx = nx + 1
 
      CASE('V')
        rcode = NF_INQ_VARID(lun,"LAT_LL_V",vid)
        if(rcode .eq. -49)then ! assume version 3 WRF
          write(6,*)' Warning: LAT_LL_V not found, looking for XLAT_V'
          rcode = NF_INQ_VARID(lun,"XLAT_V",vid)
          rcode = NF_GET_VAR_REAL(lun,vid,dum_2d)
          lat1 = dum_2d(1,1)
          wrf_version = 3
        else
          rcode = NF_GET_VAR_REAL(lun,vid,lat1)
          wrf_version = 2
        endif

        rcode = NF_INQ_VARID(lun,"LON_LL_V",vid)
        if(rcode .eq. -49)then ! assume version 3 WRF
          write(6,*)' Warning: LON_LL_V not found, looking for XLONG_V'
          rcode = NF_INQ_VARID(lun,"XLONG_V",vid)
          rcode = NF_GET_VAR_REAL(lun,vid,dum_2d)
          lon1 = dum_2d(1,1)
        else
          rcode = NF_GET_VAR_REAL(lun,vid,lon1)
        endif
        ny = ny + 1

    END SELECT

    DEALLOCATE(dum_2d)

!   if(wrf_version .eq. 2)then
!       write(6,*)' Version 2, looking for bottom_top for nz '
        rcode = NF_INQ_DIMID(lun, "bottom_top", dimid)
        rcode = NF_INQ_DIMLEN(lun, dimid, nz)
!   else
!       write(6,*)' Version 3, looking for bottom_top_stag for nz '
!       rcode = NF_INQ_DIMID(lun, "bottom_top_stag", dimid)
!       rcode = NF_INQ_DIMLEN(lun, dimid, nz)
!   endif

    write(6,*)' WRF nx,ny,nz = ',nx,ny,nz

    RETURN
  END SUBROUTINE get_wrf2_map

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_wrf_misc(cdfid, nzh, nzf, ptop, clwflag, iceflag, &
                          graupelflag)

    ! Subroutine to get a few key parameters from the WRF model for the
    ! model post processor.

    IMPLICIT NONE
    include 'netcdf.inc'
    INTEGER, INTENT(IN)                :: cdfid   ! netCDF file handle
    INTEGER, INTENT(OUT)               :: nzh     ! number of half-levels
    INTEGER, INTENT(OUT)               :: nzf     ! number of full-levels
    REAL,    INTENT(OUT)               :: ptop    ! Top pressure in Pa
    LOGICAL, INTENT(OUT)               :: clwflag ! Cloud liquid fields avail
    LOGICAL, INTENT(OUT)               :: iceflag ! Ice species avail
    LOGICAL, INTENT(OUT)               :: graupelflag  ! Grauple included

    ! Local variables

    INTEGER :: vid, rcode,mp_level
    CHARACTER(LEN=132) :: dum

    ! Get Ptop, which is not really needed for anything in the em version
    vid = NCVID(cdfid, 'P_TOP', rcode)
    rcode  = NF_GET_VAR_REAL(cdfid,vid,ptop)

    ! Get the vertical dimensions
    rcode = NF_GET_ATT_INT(cdfid, nf_global, 'BOTTOM-TOP_GRID_DIMENSION',nzh)
    nzf = nzh + 1

    ! Determine which microphysics package is used and set flags
    ! accordingly

    rcode = NF_GET_ATT_INT(cdfid, nf_global, 'MP_PHYSICS',mp_level)

    SELECT CASE (mp_level)

      CASE(0)  ! No microphysics
        clwflag = .false.
        iceflag = .false.
        graupelflag = .false.

      CASE(1) ! Kessler warm rain
        clwflag = .true.
        iceflag = .false.
        graupelflag = .false.

      CASE(2) ! Lin et al.
        clwflag = .true.
        iceflag = .true.
        graupelflag = .true.

      CASE(3)   ! NCEP 3-class
        clwflag = .true.
        iceflag = .true.
        graupelflag = .false.
      CASE(4)   ! NCEP 5-class
        clwflag = .true.
        iceflag = .true.
        graupelflag = .false.

      CASE(5)   ! Eta Ferrier 2-class
        clwflag = .true.
        iceflag = .false.
        graupelflag = .false.

      CASE DEFAULT
        print *, 'WARNING:  Cannot determine microphysics option!'
        print *, '          Assuming all species present.'
        clwflag = .true.
        iceflag = .true.
        graupelflag = .true.

    END SELECT

    RETURN

  END SUBROUTINE get_wrf_misc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_wrf2_misc(cdfid, nzh, nzf, ptop, clwflag, iceflag, &
                          graupelflag)

    ! Subroutine to get a few key parameters from the WRF model for the 
    ! model post processor.

    IMPLICIT NONE
    include 'netcdf.inc' 
    INTEGER, INTENT(IN)                :: cdfid   ! netCDF file handle
    INTEGER, INTENT(OUT)               :: nzh     ! number of half-levels
    INTEGER, INTENT(OUT)               :: nzf     ! number of full-levels
    REAL,    INTENT(OUT)               :: ptop    ! Top pressure in Pa
    LOGICAL, INTENT(OUT)               :: clwflag ! Cloud liquid fields avail
    LOGICAL, INTENT(OUT)               :: iceflag ! Ice species avail
    LOGICAL, INTENT(OUT)               :: graupelflag  ! Grauple included

    ! Local variables
     
    INTEGER :: vid, rcode,mp_level
    CHARACTER(LEN=132) :: dum

    ! Get Ptop, which is not really needed for anything in the em version
    vid = NCVID(cdfid, 'P_TOP', rcode)
    rcode  = NF_GET_VAR_REAL(cdfid,vid,ptop)

    ! Get the vertical dimensions
    rcode = NF_GET_ATT_INT(cdfid, nf_global, 'BOTTOM-TOP_GRID_DIMENSION',nzf)
    nzh = nzf - 1

    ! Determine which microphysics package is used and set flags 
    ! accordingly
 
    rcode = NF_GET_ATT_INT(cdfid, nf_global, 'MP_PHYSICS',mp_level) 
   
    SELECT CASE (mp_level)

      CASE(0)  ! No microphysics
        clwflag = .false.
        iceflag = .false. 
        graupelflag = .false.

      CASE(1) ! Kessler warm rain
        clwflag = .true.
        iceflag = .false.
        graupelflag = .false.

      CASE(2) ! Lin et al.
        clwflag = .true.
        iceflag = .true.
        graupelflag = .true.
  
      CASE(3)   ! WSM 3-class Simple
        clwflag = .true.
        iceflag = .true.
        graupelflag = .false.
  
      CASE(4)   ! WSM 5-class Simple
        clwflag = .true.
        iceflag = .true.
        graupelflag = .false.

      CASE(5)   ! Eta Ferrier 2-class
        clwflag = .true.
        iceflag = .false.
        graupelflag = .false.
      
      CASE(6)  ! WSM 6-class
        clwflag = .true.
        iceflag = .true.
        graupelflag = .true.
  
      CASE(98) ! NCEP 3-class
        clwflag = .true.
        iceflag = .false.
        graupelflag = .false.
 
      CASE(99) ! NCEP 5-class
        clwflag = .true.
        iceflag = .true.
        graupelflag = .false.

      CASE DEFAULT
        print *, 'WARNING:  Cannot determine microphysics option!'
        print *, '          Assuming all species present.'
        clwflag = .true.
        iceflag = .true.
        graupelflag = .true.

    END SELECT

    RETURN

  END SUBROUTINE get_wrf2_misc 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_wrf_1d(cdfid, varname, data1d, status)

     IMPLICIT NONE
     include 'netcdf.inc'
     INTEGER, INTENT(IN)                :: cdfid
     CHARACTER(LEN=*),INTENT(IN)        :: varname
     REAL, INTENT(INOUT)                :: data1d(:)
     INTEGER, INTENT(OUT)               :: status
     
     INTEGER                            :: vid
    
     status = NF_INQ_VARID(cdfid,TRIM(varname),vid)
     status = NF_GET_VAR_REAL(cdfid,vid,data1d)
     IF (status.ne.nf_noerr) THEN
       print *, 'NetCDF error getting ', TRIM(varname),status
       status = 1
     ELSE  
       status = 0
     ENDIF
     RETURN
  END SUBROUTINE get_wrf_1d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_wrfnc_3d(cdfid, varname, output_stagger, nx,ny,nz,time, &
                          data3d, status)

    ! Subroutine to return a 3D grid from WRF output.

    IMPLICIT NONE
    INCLUDE 'netcdf.inc'
    
    INTEGER, INTENT(IN)                    :: cdfid
    CHARACTER(LEN=*),INTENT(IN)            :: varname
    CHARACTER(LEN=1),INTENT(IN)            :: output_stagger
    INTEGER,INTENT(IN)                     :: nx,ny,nz
    INTEGER,INTENT(IN)                     :: time
    REAL,INTENT(OUT)                       :: data3d(nx,ny,nz)
    INTEGER, INTENT(OUT)                   :: status

    ! Local vars
    INTEGER                                :: vid,rcode,attid,ndims,natts
    INTEGER                                :: dims(10),dimids(10)
    INTEGER                                :: istart(10),iend(10)
    CHARACTER(LEN=10)                      :: varstagger
    CHARACTER(LEN=80)                      :: vname
    REAL, ALLOCATABLE                      :: dum3d(:,:,:)
    INTEGER                                :: i,nx_in,ny_in,nz_in,ivtype
    CHARACTER(LEN=2)                       :: conv
    status = 0 
    rcode = NF_INQ_VARID(cdfid,TRIM(varname),vid)
    IF (rcode .NE. nf_noerr) THEN
       PRINT *, 'Problem getting varid for ', TRIM(varname), ': ',rcode
       status = 1
       RETURN
    ENDIF
   
    ! Get the dimensions of this variable
    rcode = NF_INQ_VAR(cdfid,vid,vname,ivtype,ndims,dimids,natts)
    IF (ndims .NE. 4) THEN
      PRINT *, 'Dimension problem in get_wrfnc_3d.'
      PRINT *, 'Data in file should be 4D, but has ', ndims, ' dims.'
      status = 1
      RETURN
    ENDIF

    DO i = 1,ndims
      rcode = NF_INQ_DIMLEN(cdfid,dimids(i),dims(i) )
    ENDDO
  
    ! Check the dimensions 
    IF (time .GT. dims(4) ) THEN
      PRINT *, 'Requested time index exceeds input times: ', time, dims(4)
      status = 1
      RETURN
    ENDIF
    nx_in = dims(1)
    ny_in = dims(2)
    nz_in = dims(3)
   
    istart(1) = 1
    iend(1)   = nx_in
    istart(2) = 1
    iend(2)   = ny_in
    istart(3) = 1
    iend(3)   = nz_in
    istart(4) = time
    iend(4)   = 1
       

    rcode = NF_GET_ATT_TEXT(cdfid,vid,'stagger',varstagger)
    IF (rcode .NE. nf_noerr) THEN
      PRINT *, 'Problem getting stagger: ',rcode
      status = 1
      RETURN
    ENDIF
    IF (varstagger(1:1).EQ. "X") THEN
      varstagger = "U"
    ELSEIF (varstagger(1:1).EQ. "Y") THEN
      varstagger = "V"
    ELSE
      varstagger = "T"
    ENDIF

    ! Check input vs. output dimensions

    IF (varstagger(1:1) .EQ. output_stagger) THEN
       ! In/out dimensions must match
       IF ( (nx_in .NE. nx) .OR. (ny_in .NE. ny) .OR. &
            (nz_in .NE. nz) ) THEN
          status = 1
       ELSE 
         conv = "  "
       ENDIF
    ELSE
       SELECT CASE (output_stagger)
       
         CASE ('A')
           IF (varstagger(1:1) .EQ. "T") THEN
             IF ( (nx_in+1 .NE. nx) .OR. (ny_in+1 .NE. ny) .OR. &
                  (nz_in .NE. nz) ) THEN
                status = 1
             ELSE 
               conv = 'TA'
             ENDIF
           ELSEIF(varstagger(1:1) .EQ. "U") THEN
             IF ( (nx_in .NE. nx) .OR. (ny_in+1 .NE. ny) .OR. &
                  (nz_in .NE. nz) ) THEN
                status = 1
             ELSE  
               conv = "UA"
             ENDIF
           ELSEIF(varstagger(1:1) .EQ. "V") THEN
             IF ( (nx_in+1 .NE. nx) .OR. (ny_in .NE. ny) .OR. &
                  (nz_in .NE. nz) ) THEN
                status = 1
             ELSE  
                conv = "VA"
             ENDIF
           ELSE 
             PRINT *, 'Stagger conversion not supported.'
             status = 1
           ENDIF
         CASE('T')
           IF(varstagger(1:1) .EQ. "U") THEN
             IF ( (nx_in-1 .NE. nx) .OR. (ny_in .NE. ny) .OR. &
                  (nz_in .NE. nz) ) THEN
                status = 1
             ELSE  
                conv = "UT"
             ENDIF
           ELSEIF(varstagger(1:1) .EQ. "V") THEN
             IF ( (nx_in .NE. nx) .OR. (ny_in-1 .NE. ny) .OR. &
                  (nz_in .NE. nz) ) THEN
                status = 1
             ELSE 
                conv = "VT"
             ENDIF
           ELSE
             PRINT *, 'Stagger conversion not supported.'
           ENDIF
         CASE DEFAULT
           PRINT *, 'Requested stagger conversion not supported.'
           status = 1
       END SELECT
    ENDIF
    IF (status .NE. 0) THEN
      PRINT *, 'STAGGER_IN/STAGGER_REQ: ',varstagger(1:1),output_stagger 
      PRINT *, 'Mismatch in dimensions: '
      PRINT *, 'NX_IN / NX_REQ : ', nx_in, nx
      PRINT *, 'NY_IN / NY_REQ : ', ny_in, ny
      PRINT *, 'NZ_IN / NZ_REQ : ', nz_in, nz
      RETURN
    ENDIF

    ! If no destaggering is required, we can just populate the output
    ! array directly
  
    IF (conv .EQ. "  " ) THEN
      CALL NCVGT( cdfid, vid, istart, iend, data3d, rcode)
      IF (rcode .NE. NF_NOERR) THEN
        PRINT *, 'Problem getting ', TRIM(varname)
        status = 1 
        RETURN
      ENDIF
    ELSE
      ALLOCATE(dum3d(nx_in,ny_in,nz_in)) 
      CALL NCVGT( cdfid, vid, istart, iend, dum3d, rcode)
      IF (rcode .NE. NF_NOERR) THEN
        PRINT *, 'Problem getting ', TRIM(varname)
        status = 1
        RETURN
      ENDIF
      SELECT CASE (conv)
        CASE ('TA') 
          CALL arakawa_c_t2n(dum3d,nx_in,ny_in,nz_in,data3d)
        CASE ('UA')
          CALL arakawa_c_u2n(dum3d,nx_in,ny_in,nz_in,data3d) 
        CASE ('VA')
          CALL arakawa_c_v2n(dum3d,nx_in,ny_in,nz_in,data3d)
        CASE ('UT')
          CALL arakawa_c_u2t(dum3d,nx_in,ny_in,nz_in,data3d) 
        CASE ('VT')
          CALL arakawa_c_v2t(dum3d,nx_in,ny_in,nz_in,data3d)
        CASE DEFAULT
          PRINT *, 'No stagger case found.  Should not be here.'
          status = 1
      END SELECT
      DEALLOCATE(dum3d)
    ENDIF 

    RETURN 
  END SUBROUTINE get_wrfnc_3d  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_wrfnc_2d(cdfid, varname, output_stagger, nx,ny,time, &
                          data2d, status)

    ! Subroutine to return a 2D grid from WRF output.

    IMPLICIT NONE
    INCLUDE 'netcdf.inc'

    INTEGER, INTENT(IN)                    :: cdfid
    CHARACTER(LEN=*),INTENT(IN)            :: varname
    CHARACTER(LEN=1),INTENT(IN)            :: output_stagger
    INTEGER,INTENT(IN)                     :: nx,ny
    INTEGER,INTENT(IN)                     :: time
    REAL,INTENT(OUT)                       :: data2d(nx,ny)
    INTEGER, INTENT(OUT)                   :: status

    ! Local vars
    INTEGER                                :: vid,rcode,attid,ndims,natts
    INTEGER                                :: dims(10),dimids(10)
    INTEGER                                :: istart(10),iend(10)
    CHARACTER(LEN=10)                      :: varstagger
    CHARACTER(LEN=80)                      :: vname
    REAL, ALLOCATABLE                      :: dum2d(:,:)
    INTEGER                                :: i,nx_in,ny_in,ivtype
    CHARACTER(LEN=2)                       :: conv
    status = 0
    rcode = NF_INQ_VARID(cdfid,TRIM(varname),vid)
    IF (rcode .NE. nf_noerr) THEN
       PRINT *, 'Problem getting varid for ', TRIM(varname), ': ',rcode
       status = 1
       RETURN
    ENDIF

    ! Get the dimensions of this variable
    rcode = NF_INQ_VAR(cdfid,vid,vname,ivtype,ndims,dimids,natts)
    IF (ndims .NE. 3) THEN
      PRINT *, 'Dimension problem in get_wrfnc_2d.'
      PRINT *, 'Data in file should be 3D, but has ', ndims, ' dims.'
      status = 1
      RETURN
    ENDIF

    DO i = 1,ndims
      rcode = NF_INQ_DIMLEN(cdfid,dimids(i),dims(i) )
    ENDDO

    ! Check the dimensions
    IF (time .GT. dims(3) ) THEN
      PRINT *, 'Requested time index exceeds input times: ', time, dims(3)
      status = 1
      RETURN
    ENDIF
    nx_in = dims(1)
    ny_in = dims(2)

    istart(1) = 1
    iend(1)   = nx_in
    istart(2) = 1
    iend(2)   = ny_in
    istart(3) = time
    iend(3)   = 1


    rcode = NF_GET_ATT_TEXT(cdfid,vid,'stagger',varstagger)
    IF (rcode .NE. nf_noerr) THEN
      PRINT *, 'Problem getting stagger: ',rcode
      status = 1
      RETURN
    ENDIF
    IF (varstagger(1:1).EQ. "X") THEN
      varstagger = "U"
    ELSEIF (varstagger(1:1).EQ. "Y") THEN
      varstagger = "V"
    ELSE
      varstagger = "T"
    ENDIF

    ! Check input vs. output dimensions

    IF (varstagger(1:1) .EQ. output_stagger) THEN
       ! In/out dimensions must match
       IF ( (nx_in .NE. nx) .OR. (ny_in .NE. ny)) THEN
          status = 1
       ELSE
         conv = "  "
       ENDIF
    ELSE
       SELECT CASE (output_stagger)

         CASE ('A')
           IF (varstagger(1:1) .EQ. "T") THEN
             IF ( (nx_in+1 .NE. nx) .OR. (ny_in+1 .NE. ny))THEN
                status = 1
             ELSE
               conv = 'TA'
             ENDIF
           ELSEIF(varstagger(1:1) .EQ. "U") THEN
             IF ( (nx_in .NE. nx) .OR. (ny_in+1 .NE. ny)) THEN  
                status = 1
             ELSE
               conv = "UA"
             ENDIF
           ELSEIF(varstagger(1:1) .EQ. "V") THEN
             IF ( (nx_in+1 .NE. nx) .OR. (ny_in .NE. ny)) THEN  
                status = 1
             ELSE
                conv = "VA"
             ENDIF
           ELSE
             PRINT *, 'Stagger conversion not supported.'
             status = 1
           ENDIF
         CASE('T')
           IF(varstagger(1:1) .EQ. "U") THEN
             IF ( (nx_in-1 .NE. nx) .OR. (ny_in .NE. ny) ) THEN 
                status = 1
             ELSE
                conv = "UT"
             ENDIF
           ELSEIF(varstagger(1:1) .EQ. "V") THEN
             IF ( (nx_in .NE. nx) .OR. (ny_in-1 .NE. ny)) THEN  
                status = 1
             ELSE
                conv = "VT"
             ENDIF
           ELSE
             PRINT *, 'Stagger conversion not supported.'
           ENDIF
         CASE DEFAULT
           PRINT *, 'Requested stagger conversion not supported.'
           status = 1
       END SELECT
    ENDIF
    IF (status .NE. 0) THEN
      PRINT *, 'STAGGER_IN/STAGGER_REQ: ',varstagger(1:1),output_stagger 
      PRINT *, 'Mismatch in dimensions: '
      PRINT *, 'NX_IN / NX_REQ : ', nx_in, nx
      PRINT *, 'NY_IN / NY_REQ : ', ny_in, ny
      RETURN
    ENDIF

    ! If no destaggering is required, we can just populate the output
    ! array directly

    IF (conv .EQ. "  " ) THEN
      CALL NCVGT( cdfid, vid, istart, iend, data2d, rcode)
      IF (rcode .NE. NF_NOERR) THEN
        PRINT *, 'Problem getting ', TRIM(varname)
        status = 1
        RETURN
      ENDIF
    ELSE
      ALLOCATE(dum2d(nx_in,ny_in))
      CALL NCVGT( cdfid, vid, istart, iend, dum2d, rcode)
 print *, 'VAR:',varname,' Pre-stagger min/max:',minval(dum2d),maxval(dum2d)
      IF (rcode .NE. NF_NOERR) THEN
        PRINT *, 'Problem getting ', TRIM(varname)
        status = 1
        RETURN
      ENDIF
      SELECT CASE (conv)
        CASE ('TA')
          CALL arakawa_c_t2n(dum2d,nx_in,ny_in,1,data2d)
        CASE ('UA')
          CALL arakawa_c_u2n(dum2d,nx_in,ny_in,1,data2d)
        CASE ('VA')
          CALL arakawa_c_v2n(dum2d,nx_in,ny_in,1,data2d)
        CASE ('UT')
          CALL arakawa_c_u2t(dum2d,nx_in,ny_in,1,data2d)
        CASE ('VT') 
          CALL arakawa_c_v2t(dum2d,nx_in,ny_in,1,data2d)
        CASE DEFAULT
          PRINT *, 'No stagger case found.  Should not be here.'
          status = 1
      END SELECT
      DEALLOCATE(dum2d) 
    IF (varname .EQ. 'GSW') THEN
      WHERE(data2d .LT. 0) data2d = 0.
    ENDIF
    ENDIF

    RETURN
  END SUBROUTINE get_wrfnc_2d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE make_wrf_file_name(prddir,domnum,time_min,fname)

! Creates the WRF output file name based on the working directory,
! domain number, and number of minutes into simulation

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)               :: prddir
    INTEGER, INTENT(IN)                        :: domnum
    INTEGER, INTENT(IN)                        :: time_min
    CHARACTER(LEN=255), INTENT(OUT)            :: fname

    CHARACTER(LEN=2)                           :: domstr
    CHARACTER(LEN=6)                           :: timestr

    WRITE(domstr, '(I2.2)') domnum
    WRITE(timestr, '(I6.6)') time_min

    fname = TRIM(prddir) // '/wrfout_d' // domstr // '_' // timestr
    RETURN
  END SUBROUTINE make_wrf_file_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE make_wrf2_file_name(prddir,domnum,timestr,fname)

    ! Creates the WRF output file name based on the working directory, 
    ! domain number, and number of minutes into simulation

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)               :: prddir
    INTEGER, INTENT(IN)                        :: domnum
    CHARACTER(LEN=24), INTENT(IN)              :: timestr
    CHARACTER(LEN=255), INTENT(OUT)            :: fname

    CHARACTER(LEN=2)                           :: domstr

    WRITE(domstr, '(I2.2)') domnum

    fname = TRIM(prddir) // '/wrfout_d' // domstr // '_' // timestr(1:19)
    RETURN
  END SUBROUTINE make_wrf2_file_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE wrfio_wait(filename,max_wait_sec)

    IMPLICIT NONE
    CHARACTER(LEN=255)  :: filename
    CHARACTER(LEN=8)    :: date_ready
    CHARACTER(LEN=10)   :: time_ready
    LOGICAL             :: file_exists
    LOGICAL             :: file_ready
    INTEGER             :: num_checks
    INTEGER             :: max_wait_sec
    INTEGER, PARAMETER  :: pause_sec = 30
    INTEGER             :: secs_waited
    INTEGER             :: cdf, rcode, dimid, ntimes,status,vid,itimestep
    REAL                :: dt,xtime
    file_ready = .false.
    file_exists = .false.
    num_checks = 0
    print *, "Checking status of ",trim(filename)
    DO WHILE (.NOT.file_ready)
      INQUIRE(FILE=TRIM(filename), EXIST=file_exists)
      IF (.NOT. file_exists) THEN
        print *, 'File not ready: ', TRIM(filename)
        print '(A,I3,A)', 'Sleeping for ', pause_sec, ' seconds'
        CALL sleep(pause_sec)
        num_checks = num_checks + 1
        secs_waited = num_checks * pause_sec
        print '(A,I5,A)', 'Total sleep time now ', secs_waited, ' seconds'
        IF (secs_waited .GE. max_wait_sec) THEN
          PRINT *, 'IO_WAIT:  Timeout waiting for file: ', TRIM(filename)
          PRINT '(A,I5,A)', '    Maximum wait time set to ', max_wait_sec, 's'
          STOP 'io_wait'
        ENDIF
      ELSE 
        num_checks = 0
        ! Make sure it has been populated with at least 1 time period    
        DO WHILE (.NOT. file_ready)
       
          CALL open_wrfnc(filename,cdf,status)  
          rcode = NF_INQ_DIMID(cdf, "Time", dimid)
          rcode = NF_INQ_DIMLEN(cdf, dimid, ntimes)
           
          print *,"rcode/ntimes =", rcode,ntimes
          IF (ntimes.GT.0) THEN
            rcode = NF_INQ_VARID(cdf,"ITIMESTEP",vid)
            rcode = NF_GET_VAR_INT(cdf,vid,itimestep)
            
            !!!!! Modified by Wei-Ting (20130307) :            !!!!!
            !!!!!     Use XTIME instead of ITIMESTEP FOR WRFV3 !!!!!
            if (rcode .lt. 0) then
              write(6,*)' Warning: ITIMESTEP not found, looking for XTIME in wrfio_wait'
              rcode = NF_GET_ATT_REAL(cdf,0,"DT",dt)
              rcode = NF_INQ_VARID(cdf,"XTIME",vid)
              rcode = NF_GET_VAR_REAL(cdf,vid,xtime)
              itimestep = INT(xtime*60/dt)
            endif
            !!!!! End of Modifying !!!!!

            print *,"rcode/itimestep ", rcode,itimestep
            IF ((rcode .EQ. NF_NOERR) .AND. (itimestep .GE. 0))THEN
              CALL date_and_time(date_ready,time_ready)
              PRINT *, TRIM(filename), ' ready at ', date_ready, '/',time_ready
              file_ready = .true.
              !CALL sleep(pause_sec)
            ELSE 
              CALL sleep(pause_sec)
            ENDIF
          ELSE 
            CALL sleep(pause_sec)
            num_checks = num_checks + 1
            secs_waited = num_checks * pause_sec
            IF (secs_waited .GE. max_wait_sec) THEN
              PRINT *, 'IO_WAIT:  Timeout waiting for file: ', TRIM(filename)
              PRINT '(A,I5,A)', '    Maximum wait time set to ', max_wait_sec, 's'
              STOP 'io_wait'
            ENDIF
          ENDIF
          CALL close_wrfnc(cdf)
        ENDDO 
      ENDIF
    ENDDO 
    RETURN
  END SUBROUTINE wrfio_wait

END MODULE wrf_netcdf
