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

MODULE mm5v3_io

! F90 module to deal with reading data from MM5 V3 data files.

  PRIVATE
  ! MM5v3 Header fields

  TYPE mm5_big_header
    INTEGER               :: bhi (50,20)
    REAL                  :: bhr (20,20)
    CHARACTER (LEN=80)    :: bhic(50,20)
    CHARACTER (LEN=80)    :: bhrc(20,20)
  END TYPE mm5_big_header

  TYPE mm5_sub_header
    INTEGER               :: ndim
    INTEGER               :: start_index  (4)
    INTEGER               :: end_index    (4)
    REAL                  :: xtime
    CHARACTER(LEN=4)      :: staggering
    CHARACTER(LEN=4)      :: ordering
    CHARACTER(LEN=24)     :: current_date
    CHARACTER(LEN=9)      :: name
    CHARACTER(LEN=25)     :: units
    CHARACTER(LEN=46)     :: description
  END TYPE mm5_sub_header

  INTEGER, PARAMETER      :: bh_flag = 0
  INTEGER, PARAMETER      :: sh_flag = 1
  INTEGER, PARAMETER      :: eot_flag = 2
  CHARACTER(LEN=4),PARAMETER       :: crs_test = 'C   '
  CHARACTER(LEN=4),PARAMETER       :: dot_text = 'D   '
  CHARACTER(LEN=4),PARAMETER       :: press3d_text = 'YXP '
  CHARACTER(LEN=4),PARAMETER       :: sigma3d_text = 'YXS '
  CHARACTER(LEN=4),PARAMETER       :: wsigma_text = 'YXW '
  CHARACTER(LEN=4),PARAMETER       :: field2d_text = 'YX  '
  CHARACTER(LEN=4),PARAMETER       :: landuse_text = 'CA  '
  CHARACTER(LEN=4),PARAMETER       :: boundns_text = 'XSB '
  CHARACTER(LEN=4),PARAMETER       :: boundwe_text = 'YSB '
  CHARACTER(LEN=4),PARAMETER       :: boundns_w_text = 'XWB '
  CHARACTER(LEN=4),PARAMETER       :: boundwe_w_text = 'YWB '
  CHARACTER(LEN=4),PARAMETER       :: presslev_text  = 'P   '
  CHARACTER(LEN=4),PARAMETER       :: sigmalev_text  = 'S   '
  INTEGER, PARAMETER               :: lamcon_flag = 1
  INTEGER, PARAMETER               :: polstr_flag = 2
  INTEGER, PARAMETER               :: merctr_flag = 3

  PUBLIC open_mm5v3
  PUBLIC get_mm5_map
  PUBLIC get_mm5_time_info
  PUBLIC get_mm5_misc
  PUBLIC get_mm5_scalar
  PUBLIC get_mm5_1d
  PUBLIC get_mm5_2d
  PUBLIC get_mm5_3d
  PUBLIC make_data_file_name
  PUBLIC make_terrain_file_name
  PUBLIC io_wait

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE open_mm5v3 (fname, lun, status)
    
    ! Opens an MM5V3 file

    IMPLICIT NONE
   
    ! Arguments

    CHARACTER(LEN=255),INTENT(IN)          :: fname
    INTEGER, INTENT(OUT)                   :: lun
    INTEGER, INTENT(OUT)                   :: status

    ! Local variables

    INTEGER                                :: i
    LOGICAL                                :: used

    status = 0
    
    lun = -1
    ! Find an available unit number
    find_lun: DO i = 7,1023
      INQUIRE (UNIT=i, OPENED=used)
      IF(used) THEN
        CYCLE find_lun
      ELSE
        lun = i
        EXIT find_lun
      ENDIF
    ENDDO find_lun
    IF (lun .LT. 0) THEN
      PRINT '(A)', 'No available unit numbers!'
      status = 1
    ELSE
      OPEN ( FILE=TRIM(fname), UNIT=lun, ACCESS='SEQUENTIAL', &
             FORM='UNFORMATTED', STATUS='UNKNOWN', IOSTAT=status, &
             POSITION='REWIND')
    ENDIF
  END SUBROUTINE open_mm5v3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_mm5_big_header(lun, bh)
    IMPLICIT NONE
    INTEGER, INTENT(IN)               :: lun
    TYPE(mm5_big_header), INTENT(OUT) :: bh 
    LOGICAL                           :: opened
    INTEGER                           :: hdr_flag
    INQUIRE (UNIT=lun, OPENED=opened)
    IF (.NOT.opened) THEN
      PRINT *, 'GET_MM5_BIG_HEADER: File not opened yet!'
    ELSE
      REWIND (lun)
      READ (lun) hdr_flag
      IF (hdr_flag .EQ. bh_flag) THEN
        READ(lun) bh%bhi, bh%bhr, bh%bhic, bh%bhrc
      ELSE
        PRINT *, 'GET_MM5_BIG_HEADER: BH not found!'
      ENDIF
    ENDIF
    RETURN
  END SUBROUTINE get_mm5_big_header
         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_mm5_map(lun, projection, proj_cent_lat, proj_cent_lon, &
                         truelat1, truelat2, cone_factor, pole_point, &
                         grid_spacing, nx, ny, &
                         status)

    IMPLICIT NONE
    
    ! Arguments
    INTEGER, INTENT(IN)                    :: lun  ! Logical unit number
    CHARACTER (LEN=32), INTENT(OUT)        :: projection
    REAL, INTENT(OUT)                      :: proj_cent_lat
    REAL, INTENT(OUT)                      :: proj_cent_lon
    REAL, INTENT(OUT)                      :: truelat1
    REAL, INTENT(OUT)                      :: truelat2
    REAL, INTENT(OUT)                      :: cone_factor
    REAL, INTENT(OUT)                      :: pole_point
    REAL, INTENT(OUT)                      :: grid_spacing
    INTEGER, INTENT(OUT)                   :: nx
    INTEGER, INTENT(OUT)                   :: ny
    INTEGER, INTENT(OUT)                   :: status
    ! Local Variables

    LOGICAL                                :: file_opened
    TYPE (mm5_big_header)                  :: bh
    TYPE (mm5_sub_header)                  :: sh
    INTEGER                                :: header_flag

    status = 0
    ! Inquire to see if file is already open. 
    INQUIRE (UNIT=lun, OPENED=file_opened)
    IF (.NOT.file_opened) THEN
      PRINT '(A,I4)', 'File unit not opened: ', lun
      PRINT '(A)', 'Call OPEN_MM5V3 first.'
      status=1
    ENDIF
    IF (status .EQ. 0) THEN
      CALL get_mm5_big_header(lun,bh)
      ! Set projection
      SELECT CASE (bh%bhi(7,1))
        CASE (lamcon_flag)
          projection = 'LAMBERT CONFORMAL               '
        CASE (polstr_flag) 
          projection = 'POLAR STEREOGRAPHIC             '
        CASE (merctr_flag) 
          projection = 'MERCATOR                        '
      END SELECT
      ! Get projection center lat/lon
      proj_cent_lat = bh%bhr(2,1)
      proj_cent_lon = bh%bhr(3,1)
      ! Get cone factor and true latitudes
      truelat1 = bh%bhr(5,1)
      truelat2 = bh%bhr(6,1)
      cone_factor = bh%bhr(4,1)
      pole_point = bh%bhr(7,1)
      ! Get grid spacing and dimensions
      grid_spacing = bh%bhr(9,1)
      nx = bh%bhi(17,1)
      ny = bh%bhi(16,1)

      ! Print out diagnostics
      PRINT '(A)', 'MM5 Map Projection Parameters'
      PRINT '(A)', '-----------------------------'
      PRINT '(2A)','Proj: ', TRIM(projection)
      PRINT '(A,F10.3)', 'Proj Center Lat: ', proj_cent_lat
      PRINT '(A,F10.3)', 'Proj Center Lon: ', proj_cent_lon
      PRINT '(A,F10.3)', 'True Latitude 1: ', truelat1
      PRINT '(A,F10.3)', 'True Latitude 2: ', truelat2
      PRINT '(A,F6.3)', 'Cone Factor: ', cone_factor
      PRINT '(A,F10.3)', 'Pole Point: ', pole_point
      PRINT '(A,F10.1)', 'Grid Spacing: ', grid_spacing
      PRINT '(A,I5)', 'X Dimension: ', nx
      PRINT '(A,I5)', 'Y Dimension: ', ny
    ENDIF
  END SUBROUTINE get_mm5_map   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE find_mm5_var(lun, varname, valtime, sh, status)
 
    ! Searches an MM5V3 file for a variable and returns its subheader

    IMPLICIT NONE
    
    ! Arguments
    INTEGER, INTENT(IN)                 :: lun
    CHARACTER(LEN=9), INTENT(IN)        :: varname
    TYPE(mm5_sub_header), INTENT(OUT)   :: sh
    CHARACTER(LEN=24), INTENT(IN)       :: valtime
    INTEGER, INTENT(OUT)                 :: status

    ! Local variables
    LOGICAL                             :: file_opened
    INTEGER                             :: nx,ny,nz
    INTEGER                             :: header_flag
    LOGICAL                             :: var_found
  
    status = 0
    var_found = .false.
    INQUIRE(UNIT=lun, OPENED=file_opened)
    IF (.NOT.file_opened) THEN 
      PRINT *, 'FIND_MM5_VAR: File not open in search for: ', TRIM(varname)
      status = 1
    ELSE
      REWIND (lun)
      search_loop: DO WHILE (.NOT. var_found)
        READ (lun, END=999, IOSTAT=status) header_flag
        SELECT CASE (header_flag)
          CASE(bh_flag) 
            READ(lun)
          CASE(eot_flag) 
            CYCLE search_loop
          CASE(sh_flag)
            READ (lun) sh%ndim, sh%start_index, sh%end_index, sh%xtime, &
                       sh%staggering, sh%ordering, sh%current_date, &
                       sh%name, sh%units, sh%description
            IF (sh%name .EQ. varname) THEN
              IF (valtime(1:10) .LT. '1960-01-01') THEN
                var_found = .true.
              ELSE
                IF (valtime(1:16).EQ.sh%current_date(1:16))THEN
                  var_found = .true.
                ENDIF
              ENDIF
            ELSE
              READ (lun)      ! skip the next data item
            ENDIF
          CASE DEFAULT
            PRINT '(A,I4)', 'Bad header flag:', header_flag
        END SELECT
      ENDDO search_loop
 999  IF (.NOT. var_found) THEN
        PRINT *, 'FIND_MM5_VAR: Unable to find ', TRIM(varname), &
                 ' at time ', valtime(1:16)
        status = 1
      ENDIF
    ENDIF
  END SUBROUTINE find_mm5_var
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_mm5_scalar(lun, varname, valtime, vardata, status)
  
    ! Gets a scalar data value from an MM5 data file
 
    IMPLICIT NONE

    ! Arguments
    INTEGER, INTENT(IN)                :: lun
    CHARACTER(LEN=9), INTENT(IN)       :: varname
    CHARACTER(LEN=24),INTENT(IN)       :: valtime
    REAL, INTENT(OUT)                :: vardata
    INTEGER, INTENT(OUT)               :: status

    ! Local Variables
    TYPE (mm5_sub_header)              :: sh
    LOGICAL                            :: file_opened

    status = 0
    INQUIRE(UNIT=lun, OPENED=file_opened)
    IF (.NOT. file_opened) THEN
      PRINT '(A)', 'You need to call open_mm5v3 first!'
      status = 1
    ELSE
      CALL find_mm5_var(lun, varname, valtime, sh, status)
      IF (status .NE. 0) THEN 
        PRINT '(3A,I4)', 'Variable ', TRIM(varname), ' not found in unit ',lun
      ELSE
        IF (sh%ndim .EQ. 0) THEN
          READ (lun) vardata
        ELSE
          PRINT '(3A)', 'Variable ', TRIM(varname),' is not a scalar!'
          PRINT '(A,I4)', 'Number of dimensions found: ', sh%ndim
          status = 1
        ENDIF
      ENDIF
    ENDIF
  END SUBROUTINE get_mm5_scalar
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_mm5_1d(lun,varname,valtime,vardata,status)
    
    ! Gets a 1D array variable from the MM5V3 file
   IMPLICIT NONE

    ! Arguments
    INTEGER, INTENT(IN)                :: lun
    CHARACTER(LEN=9), INTENT(IN)       :: varname
    CHARACTER(LEN=24),INTENT(IN)       :: valtime
    REAL, INTENT(OUT)                  :: vardata(:)
    INTEGER, INTENT(OUT)               :: status

    ! Local Variables
    TYPE (mm5_sub_header)              :: sh
    LOGICAL                            :: file_opened

    status = 0
    INQUIRE(UNIT=lun, OPENED=file_opened)
    IF (.NOT. file_opened) THEN
      PRINT '(A)', 'You need to call open_mm5v3 first!'
      status = 1
    ELSE
      CALL find_mm5_var(lun, varname, valtime, sh, status)
      IF (status .NE. 0) THEN 
        PRINT '(3A,I4)', 'Variable ', TRIM(varname), ' not found in unit ',lun
      ELSE
        IF (sh%ndim .EQ. 1) THEN
          READ (lun) vardata
        ELSE
          PRINT '(3A)', 'Variable ', TRIM(varname),' is not 1-D!'
          PRINT '(A,I4)', 'Number of dimensions found: ', sh%ndim
          status = 1
        ENDIF
      ENDIF
    ENDIF
  END SUBROUTINE get_mm5_1d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_mm5_2d(lun,varname,valtime,vardata,output_stagger,status)
    
    ! Gets a 2D array variable from the MM5V3 file, reforms the array
    ! to be (nx,ny) (MM5 uses ny,nx) and optionally transforms the 
    ! data to the desired output grid (cross or dot points) if needed

    IMPLICIT NONE

    ! Arguments
    INTEGER, INTENT(IN)                :: lun
    CHARACTER(LEN=9), INTENT(IN)       :: varname
    CHARACTER(LEN=24),INTENT(IN)       :: valtime
    REAL, INTENT(OUT)                  :: vardata( : , : )
    CHARACTER(LEN=4), INTENT(IN)       :: output_stagger
    INTEGER, INTENT(OUT)               :: status

    ! Local Variables
    TYPE (mm5_sub_header)              :: sh
    LOGICAL                            :: file_opened
    INTEGER                            :: xdim, ydim, i, j
    REAL, ALLOCATABLE                  :: temp_array ( : , : )

    status = 0
    INQUIRE(UNIT=lun, OPENED=file_opened)
    IF (.NOT. file_opened) THEN
      PRINT '(A)', 'You need to call open_mm5v3 first!'
      status = 1
    ELSE
      CALL find_mm5_var(lun, varname, valtime, sh, status)
      IF (status .NE. 0) THEN 
        PRINT '(3A,I4)', 'Variable ', TRIM(varname), ' not found in unit ',lun
      ELSE
        IF (sh%ndim .EQ. 2) THEN
          xdim = sh%end_index(2) - sh%start_index(2) + 1
          ydim = sh%end_index(1) - sh%start_index(1) + 1
          ALLOCATE(temp_array (ydim,xdim))
          READ (lun) temp_array
          ! Re-order the array to by (nx,ny)
          DO j = 1,ydim
            DO i = 1,xdim
              vardata(i,j) = temp_array(j,i)
            ENDDO
          ENDDO
          DEALLOCATE (temp_array)
          ! Determine if we need to do crs2dot or dot2crs
          IF ((sh%staggering .EQ. 'C   ').AND.&
              (output_stagger .EQ. 'D   ') ) THEN
            CALL crs2dot(vardata,xdim,ydim)
          ELSE IF ((sh%staggering .EQ. 'D   ').AND.&
              (output_stagger .EQ. 'C   ') ) THEN
            CALL dot2crs(vardata,xdim,ydim)
          ENDIF
        ELSE
          PRINT '(3A)', 'Variable ', TRIM(varname),' is not 2-D!'
          PRINT '(A,I4)', 'Number of dimensions found: ', sh%ndim
          status = 1
        ENDIF
      ENDIF
    ENDIF
  END SUBROUTINE get_mm5_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_mm5_3d(lun,varname,valtime,vardata,output_stagger,status)
    
    ! Gets a 3D array variable from the MM5V3 file, reforms the array
    ! to be (nx,ny,nz) (MM5 uses ny,nx) and optionally transforms the 
    ! data to the desired output grid (cross or dot points) if needed

    IMPLICIT NONE

    ! Arguments
    INTEGER, INTENT(IN)                :: lun
    CHARACTER(LEN=9), INTENT(IN)       :: varname
    CHARACTER(LEN=24),INTENT(IN)       :: valtime
    REAL, INTENT(OUT)                  :: vardata( : , : , : )
    CHARACTER(LEN=4), INTENT(IN)       :: output_stagger
    INTEGER, INTENT(OUT)               :: status

    ! Local Variables
    TYPE (mm5_sub_header)              :: sh
    LOGICAL                            :: file_opened
    INTEGER                            :: xdim, ydim, zdim, i, j, k
    REAL, ALLOCATABLE                  :: temp_array ( : , : , : )

    status = 0
    INQUIRE(UNIT=lun, OPENED=file_opened)
    IF (.NOT. file_opened) THEN
      PRINT '(A)', 'You need to call open_mm5v3 first!'
      status = 1
    ELSE
      CALL find_mm5_var(lun, varname, valtime, sh, status)
      IF (status .NE. 0) THEN 
        PRINT '(3A,I4)', 'Variable ', TRIM(varname), ' not found in unit ',lun
      ELSE
        IF (sh%ndim .EQ. 3) THEN
          xdim = sh%end_index(2) - sh%start_index(2) + 1
          ydim = sh%end_index(1) - sh%start_index(1) + 1
          zdim = sh%end_index(3) - sh%start_index(3) + 1
          ALLOCATE(temp_array (ydim,xdim, zdim))
          READ (lun) temp_array
          ! Re-order the array to be (nx,ny,nz) and make it go from bottom of
          ! atmosphere to top of atmosphere
          DO k = 1,zdim
            DO j = 1,ydim
              DO i = 1,xdim
                vardata(i,j,k) = temp_array(j,i,zdim-k+1)
              ENDDO
            ENDDO
          ENDDO
          DEALLOCATE(temp_array)

          ! Determine if we need to do crs2dot or dot2crs
          IF ((sh%staggering .EQ. 'C   ').AND.&
              (output_stagger .EQ. 'D   ') ) THEN
            DO k=1,zdim
              CALL crs2dot(vardata(:,:,k),xdim,ydim)
            ENDDO
          ELSE IF ((sh%staggering .EQ. 'D   ').AND.&
                   (output_stagger .EQ. 'C   ') ) THEN
            DO k=1,zdim
              CALL dot2crs(vardata(:,:,k),xdim,ydim)
            ENDDO
          ENDIF
        ELSE
          PRINT '(3A)', 'Variable ', TRIM(varname),' is not 3-D!'
          PRINT '(A,I4)', 'Number of dimensions found: ', sh%ndim
          status = 1
        ENDIF
      ENDIF
    ENDIF
  END SUBROUTINE get_mm5_3d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_mm5_time_info(lun, sim_start_year, sim_start_month, &
                               sim_start_day, sim_start_hour, sim_start_min, &
                               sim_start_sec, sim_start_frac, sim_stop_min, &
                               dom_start_min, dom_stop_min, &
                               fdda_on, fdda_start_min, fdda_stop_min, &
                               tapfrq, buffrq, status)

  ! Gets all of the time variables related to domain start/stop, simulation
  ! start/stop, output frequency, etc.
  
    IMPLICIT NONE
   
    ! Arguments
  
    INTEGER, INTENT(IN)                :: lun
    INTEGER, INTENT(OUT)               :: sim_start_year
    INTEGER, INTENT(OUT)               :: sim_start_month
    INTEGER, INTENT(OUT)               :: sim_start_day
    INTEGER, INTENT(OUT)               :: sim_start_hour
    INTEGER, INTENT(OUT)               :: sim_start_min
    INTEGER, INTENT(OUT)               :: sim_start_sec
    INTEGER, INTENT(OUT)               :: sim_start_frac
    REAL, INTENT(OUT)                  :: sim_stop_min
    REAL, INTENT(OUT)                  :: dom_start_min
    REAL, INTENT(OUT)                  :: dom_stop_min
    LOGICAL,INTENT(OUT)                :: fdda_on
    REAL, INTENT(OUT)                  :: fdda_start_min
    REAL, INTENT(OUT)                  :: fdda_stop_min
    REAL, INTENT(OUT)                  :: tapfrq
    REAL, INTENT(OUT)                  :: buffrq
    INTEGER, INTENT(OUT)               :: status

    ! Locals
    LOGICAL                            :: file_opened
    TYPE(mm5_big_header)               :: bh

    INQUIRE(UNIT=lun, OPENED=file_opened)
    IF (.NOT.file_opened) THEN
      PRINT *, 'GET_MM5_TIME_INFO: Data file not open!'
      status = 1
    ELSE
      ! Find the big header by rewinding to the beginning
      CALL get_mm5_big_header(lun,bh)
      sim_start_year = bh%bhi(5,11)
      sim_start_month = bh%bhi(6,11)
      sim_start_day = bh%bhi(7,11)
      sim_start_hour = bh%bhi(8,11)
      sim_start_min = bh%bhi(9,11)
      sim_start_sec = bh%bhi(10,11)
      sim_start_frac = bh%bhi(11,11)
      sim_stop_min = bh%bhr(1,12)
      tapfrq = bh%bhr(4,12)
      buffrq = bh%bhr(5,12)
      dom_start_min = bh%bhr(1,14)
      dom_stop_min = MIN(bh%bhr(2,14),sim_stop_min)
      print *, 'dom_start_min/dom_stop_min=',dom_start_min,dom_stop_min
      IF ((bh%bhi(1,16).EQ.1).OR.(bh%bhi(6,16).EQ.1).OR. &
         (bh%bhi(14,16).EQ.1)) THEN
         fdda_on = .true.
         fdda_start_min = bh%bhr(1,16)
         fdda_stop_min = bh%bhr(2,16)
      ELSE
         fdda_on = .false.
         fdda_start_min = 0.
         fdda_stop_min = 0.
      ENDIF
    ENDIF
  END SUBROUTINE get_mm5_time_info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE crs2dot(field,dim1,dim2)
  
    ! Transforms grid on cross-points to grid on dot points
    
    IMPLICIT NONE
    INTEGER :: dim1 , dim2
    REAL , DIMENSION(dim1,dim2) :: field,dummy
    INTEGER :: i , j 
   
    dummy(2:dim1-1,2:dim2-1)  = ( field(1:dim1-2,1:dim2-2) +  & 
                                  field(1:dim1-2,2:dim2-1) +  &  
                                  field(2:dim1-1,1:dim2-2) +  &
                                  field(2:dim1-1,2:dim2-1))*0.25
 
    dummy(2:dim1-1,1:dim2:dim2-1)=( field(1:dim1-2,1:dim2-1:dim2-2) + & 
                                    field(2:dim1-1,1:dim2-1:dim2-2)) * 0.5

    dummy(1:dim1:dim1-1,2:dim2-1)=(field(1:dim1-1:dim1-2,1:dim2-2) + & 
                                   field(1:dim1-1:dim1-2,2:dim2-1) )*0.5

    dummy(1:dim1:dim1-1,1:dim2:dim2-1) = field(1:dim1-1:dim1-2,1:dim2-1:dim2-2)

    field = dummy

  END SUBROUTINE crs2dot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE dot2crs(vardot, dim1, dim2)

    ! Transforms variable on dot points to cross points
    IMPLICIT NONE

    INTEGER, INTENT(IN)               :: dim1
    INTEGER, INTENT(IN)               :: dim2
    REAL, INTENT(INOUT)               :: vardot(dim1,dim2)
    REAL                              :: dummy(dim1,dim2)  
    INTEGER                           :: i
    INTEGER                           :: j
  
!-----------------------------------------------------------------------
!--   INTERPOLATE DOT POINT VALUES TO CROSS POINTS USING 
!--   FOUR-POINT INTERPOLATION.
!-----------------------------------------------------------------------

    DO j = 1, dim2-1
      DO i = 1, dim1-1 

        dummy(i,j) = 0.25 * ( vardot(i,j) + vardot(i+1,j) + &
                              vardot(i,j+1) + vardot(i+1,j+1) )

      ENDDO
    ENDDO
    dummy(dim1,:) = vardot(dim1,:)
    dummy(:,dim2) = vardot(:,dim2)
    vardot = dummy

  END SUBROUTINE dot2crs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_mm5_misc(lun, ksig, ptop, PmslBase, TmslBase, &
                          dTdlnPBase, TisoBase, &
                          clwflag, iceflag,graupelflag)

    IMPLICIT NONE
    INTEGER, INTENT(IN)       :: lun
    INTEGER, INTENT(OUT)      :: ksig
    REAL,    INTENT(OUT)      :: ptop
    REAL,    INTENT(OUT)      :: PmslBase
    REAL,    INTENT(OUT)      :: TmslBase
    REAL,    INTENT(OUT)      :: dTdlnPBase
    REAL,    INTENT(OUT)      :: TisoBase
    LOGICAL, INTENT(OUT)      :: clwflag
    LOGICAL, INTENT(OUT)      :: iceflag
    LOGICAL, INTENT(OUT)      :: graupelflag
    LOGICAL                   :: opened
    TYPE (mm5_big_header)     :: bh
    INQUIRE(UNIT=lun, OPENED=opened)

    clwflag = .false.
    iceflag = .false.
    graupelflag = .false.
    IF (.NOT.opened) THEN 
      PRINT *, 'GET_MM5_MISC: File not opened.  Call mm5v3_open first.'
    ELSE
      CALL get_mm5_big_header(lun,bh)
      ksig = bh%bhi(12,11)
      Ptop = bh%bhr(2,2)
      PmslBase = bh%bhr(2,5)
      TmslBase = bh%bhr(3,5)
      dTdlnPbase = bh%bhr(4,5)
      TisoBase = bh%bhr(5,5)
      IF (bh%bhi(20,11).EQ.1) clwflag = .true.
      IF (bh%bhi(18,11).EQ.1) iceflag = .true.
      IF (bh%bhi(19,11).EQ.1) graupelflag = .true.
    ENDIF
    RETURN
  END SUBROUTINE get_mm5_misc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE make_data_file_name(mm5_data_root, domain_num_str, &
                                 split_output, time_period,file_name, &
                                 file_num3)
    IMPLICIT NONE
    CHARACTER(LEN=255),INTENT(IN)    :: mm5_data_root
    CHARACTER(LEN=1),INTENT(IN)      :: domain_num_str
    LOGICAL, INTENT(IN)              :: split_output
    INTEGER, INTENT(IN)              :: time_period
    CHARACTER(LEN=255),INTENT(OUT)   :: file_name
    LOGICAL, INTENT(IN)              :: file_num3
    CHARACTER(LEN=2)                 :: time_period_str
    CHARACTER(LEN=3)                 :: time_period_str3

    file_name = TRIM(mm5_data_root)//'/mm5prd/raw/MMOUT_DOMAIN'//domain_num_str
    IF (split_output) THEN
      IF (.NOT.file_num3) THEN      
        WRITE(time_period_str, '(I2.2)') time_period
        file_name = TRIM(file_name) // '_' // time_period_str
      ELSE
        WRITE(time_period_str3, '(I3.3)') time_period
        file_name = TRIM(file_name) // '_' // time_period_str3
      ENDIF
    ENDIF
    RETURN

  END SUBROUTINE make_data_file_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE make_terrain_file_name(mm5_data_root, domain_num_str,file_name)
    IMPLICIT NONE
    CHARACTER(LEN=255),INTENT(IN)    :: mm5_data_root
    CHARACTER(LEN=1),INTENT(IN)      :: domain_num_str
    CHARACTER(LEN=255),INTENT(OUT)   :: file_name

    file_name = TRIM(mm5_data_root)//'/static/TERRAIN_DOMAIN'//domain_num_str
    RETURN
  END SUBROUTINE make_terrain_file_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE io_wait(filename,max_wait_sec)

    IMPLICIT NONE
    CHARACTER(LEN=255)  :: filename
    CHARACTER(LEN=255)  :: flagfile
    CHARACTER(LEN=8)    :: date_ready
    CHARACTER(LEN=10)    :: time_ready
    LOGICAL             :: file_ready
    INTEGER             :: num_checks
    INTEGER             :: max_wait_sec
    INTEGER, PARAMETER  :: pause_sec = 30
    INTEGER             :: secs_waited
    file_ready = .false.
    num_checks = 0
    flagfile = TRIM(filename) // '.done'
    DO WHILE (.NOT.file_ready)
      INQUIRE(FILE=TRIM(flagfile), EXIST=file_ready)
      ! In case this file was just created, wait to 
      ! give the file a chance to be completely written.  Also, this
      ! keeps us from banging on the disk unnecessarily.
      IF (.NOT. file_ready) THEN
        print *, 'File not ready: ', TRIM(filename)
        print '(A,I3,A)', 'Sleeping for ', pause_sec, ' seconds'
        CALL sleep(pause_sec)
        num_checks = num_checks + 1
        secs_waited = num_checks * pause_sec
        print '(A,I5,A)', 'Total sleep time now ', secs_waited, ' seconds'
        IF (secs_waited .GE. max_wait_sec) THEN
          PRINT *, 'IO_WAIT:  Timeout waiting for file: ', TRIM(filename)
          PRINT '(A,I5,A)', '    Maximum wait time set to ', max_wait_sec, 's'
          STOP 'io_timeout'
        ENDIF
      ELSE 
        CALL date_and_time(date_ready,time_ready) 
        PRINT *, TRIM(filename), ' ready at ', date_ready, '/',time_ready
      ENDIF
    ENDDO 
    RETURN
  END SUBROUTINE io_wait
  END MODULE mm5v3_io
