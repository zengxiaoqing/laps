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

MODULE wrfv1_netcdf

! F90 module to deal with reading WRF output files (WRFv1 netcdf)

  USE grid_utils
  PRIVATE
    
  PUBLIC open_wrfnc
  PUBLIC close_wrfnc
  PUBLIC get_wrfnc_map
  PUBLIC get_wrf_time_info
  PUBLIC get_wrf_misc
  PUBLIC get_wrf_scalar
  PUBLIC get_wrf_1d
  PUBLIC get_wrfnc_2d
  PUBLIC get_wrfnc_3d
  PUBLIC make_wrf_file_name 
  PUBLIC make_wrfstatic_file_name
  PUBLIC wrfio_wait

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE open_wrfnc(fname, lun, status)

    ! Opens a WRF v1 netCDF file, returning the CDF id of the file

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
  SUBROUTINE get_wrf_map(lun, projection, lat1, lon1, stdlon, &
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
    vid = NCVID(lun, 'La1', rcode)
    CALL NCVGT1(lun, vid, 1, lat1, rcode)
    vid = NCVID(lun, 'Lo1', rcode)
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

  END SUBROUTINE get_wrf_map 
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
        graupelflag = .true.

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
      PRINT *, 'NZ_IN / NX_REQ : ', nz_in, nz
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

    fname = TRIM(prddir) // '/wrfout_' // domstr // '_' // timestr
    RETURN
  END SUBROUTINE make_wrf_file_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE wrfv1_netcdf
