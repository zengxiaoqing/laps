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

MODULE lapsprep_wrf

! PURPOSE
! =======
! Module to contain the various output routines needed for lapsprep
! to support initializition of the WRF model.
!
! SUBROUTINES CONTAINED
! =====================
! output_gribprep_format  - Used to support WRF initializations
! output_gribprep_header  - Writes the grib prep headers
! output_metgrid_format  - Used to support WRF initializations using WPS
! output_metgrid_header  - Writes the metgrid headers
! REMARKS
! =======
! 
!
! HISTORY
! =======
! 4 Dec 2000 -- Original -- Brent Shaw

  USE setup
  USE laps_static
  USE date_pack
! USE lapsprep_constants
  IMPLICIT NONE

  PRIVATE
  INTEGER, PARAMETER :: gp_version = 4
  INTEGER, PARAMETER :: gp_version_wps = 5
  REAL, PARAMETER    :: xfcst = 0.
  CHARACTER(LEN=32),PARAMETER :: source = &
     'LAPS ANALYSIS                   '
  CHARACTER (LEN=8), PARAMETER:: knownloc='SWCORNER'
  CHARACTER (LEN=24) :: hdate
  INTEGER            :: llflag
  CHARACTER (LEN=9)  :: field
  CHARACTER (LEN=25) :: units
  CHARACTER (LEN=46) :: desc 
  INTEGER, PARAMETER :: output_unit = 78
  REAL, PARAMETER    :: slp_level = 201300.0

  PUBLIC output_gribprep_format
  PUBLIC output_metgrid_format
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE output_gribprep_format(p, t, ht, u, v, rh, slp, psfc, &
                               lwc, rai, sno, ice, pic, snocov,tskin)

  !  Subroutine of lapsprep that will build a file the
  !  WRFSI "gribprep" format that can be read by hinterp

  IMPLICIT NONE

  ! Arguments
 
  REAL, INTENT(IN)                   :: p(:)        ! Pressure (hPa)
  REAL, INTENT(IN)                   :: t(:,:,:)    ! Temperature (K)
  REAL, INTENT(IN)                   :: ht(:,:,:)   ! Height (m)
  REAL, INTENT(IN)                   :: u(:,:,:)    ! U-wind (m s{-1})
  REAL, INTENT(IN)                   :: v(:,:,:)    ! V-wind (m s{-1})
  REAL, INTENT(IN)                   :: rh(:,:,:)   ! Relative Humidity (%)
  REAL, INTENT(IN)                   :: slp(:,:)    ! Sea-level Pressure (Pa)
  REAL, INTENT(IN)                   :: psfc(:,:)   ! Surface Pressure (Pa)
  REAL, INTENT(IN)                   :: lwc(:,:,:)  ! Cloud liquid (kg/kg)
  REAL, INTENT(IN)                   :: rai(:,:,:)  ! Rain (kg/kg)
  REAL, INTENT(IN)                   :: sno(:,:,:)  ! Snow (kg/kg)
  REAL, INTENT(IN)                   :: ice(:,:,:)  ! Ice (kg/kg)
  REAL, INTENT(IN)                   :: pic(:,:,:)  ! Graupel (kg/kg)
  REAL, INTENT(IN)                   :: snocov(:,:) ! Snow cover (fract)
  REAL, INTENT(IN)                   :: tskin(:,:)  ! Skin temperature
  
  ! Local Variables
  
  INTEGER            :: valid_mm, valid_dd
  CHARACTER (LEN=256):: output_file_name
  REAL, ALLOCATABLE  :: d2d(:,:)
  REAL, ALLOCATABLE  :: p_pa(:)
  INTEGER            :: k,yyyyddd
  INTEGER            :: istatus
  REAL               :: r_missing_data
 
  ! Allocate a scratch 2d array
  ALLOCATE (d2d (x,y) )
  ALLOCATE (p_pa (z3+1))
  ! Build the output file name
 
  output_prefix = TRIM(laps_data_root)// '/lapsprd/lapsprep/wrf/LAPS'
  yyyyddd = valid_yyyy*1000 + valid_jjj
  CALL wrf_date_to_ymd(yyyyddd, valid_yyyy, valid_mm, valid_dd) 
  WRITE(hdate, '(I4.4,"-",I2.2,"-",I2.2,"_",I2.2,":",I2.2,":00.0000")') &
          valid_yyyy, valid_mm, valid_dd, valid_hh, valid_min
  IF (valid_min .EQ. 0) THEN
    output_file_name = TRIM(output_prefix) // ':' // hdate(1:13)
  ELSE
    output_file_name = TRIM(output_prefix) // ':' // hdate(1:16)
  ENDIF
  !  Open the file for sequential, unformatted output
  OPEN ( FILE   = TRIM(output_file_name)    , &
         UNIT   = output_unit        , &
         FORM   = 'UNFORMATTED' , &
         STATUS = 'UNKNOWN'     , &
         ACCESS = 'SEQUENTIAL'    )

  ! Convert p levels from mb to Pascals

  p_pa = p * 100.

  ! Set llflag based on grid type

  IF      ( grid_type(1:8)  .EQ. 'mercator'                 ) THEN
    llflag = 1
  ELSE IF ( ( grid_type(1:24) .EQ. 'secant lambert conformal' ) .or. &
           ( grid_type(1:28) .EQ.  'tangential lambert conformal' ) )THEN
    llflag = 3
  ELSE IF ( grid_type(1:19) .EQ. 'polar stereographic'      ) THEN
    llflag = 5
  ELSE
    PRINT '(A,A,A)','Unknown map projection: ',TRIM(grid_type),'.  I quit.'
    STOP 'unknown_projection'
  END IF

  PRINT *, 'GRIBPREP VERSION =', gp_version
  PRINT *, 'SOURCE = ', source
  PRINT *, 'HDATE = ',hdate
  PRINT *, 'XFCST = ', xfcst 
  PRINT *, 'NX = ', X
  PRINT *, 'NY = ', Y
  PRINT *, 'IPROJ = ', LLFLAG
  PRINT *, 'KNOWNLOC = ', knownloc
  PRINT *, 'STARTLAT = ',LA1
  PRINT *, 'STARTLON = ',LO1
  PRINT *, 'DX = ', DX
  PRINT *, 'DY = ', DY
  PRINT *, 'XLONC = ', LOV
  PRINT *, 'TRUELAT1 = ', LATIN1
  PRINT *, 'TRUELAT2 = ', LATIN2  

 ! Output temperature
  field = 'T        '
  units = 'K                        '
  desc  = 'Temperature                                   '                
  PRINT *, 'FIELD = ', field
  PRINT *, 'UNITS = ', units
  PRINT *, 'DESC =  ',desc
  var_t : DO k = 1 , z3 + 1
    IF (( p_pa(k) .GT. 100100).AND.(p_pa(k).LT.200000) ) THEN
      CYCLE var_t
    ENDIF
    CALL write_gribprep_header(field,units,desc,p_pa(k))
    d2d = t(:,:,k)
    WRITE ( output_unit ) d2d
    PRINT '(A,F9.1,A,F5.1,A,F5.1)','Level (Pa):',p_pa(k),' Min: ', &
           MINVAL(d2d),' Max: ', MAXVAL(d2d)
  ENDDO var_t

  ! Do u-component of wind
  field = 'U        '
  units = 'm s{-1}                  '
  desc = 'u-component of velocity, rotated to grid      '
  PRINT *, 'FIELD = ', field
  PRINT *, 'UNITS = ', units
  PRINT *, 'DESC =  ',desc
  var_u : DO k = 1 , z3 + 1
    IF ( ( p_pa(k) .GT. 100100 ) .AND. ( p_pa(k) .LT. 200000 ) ) THEN
      CYCLE var_u
    END IF
    CALL write_gribprep_header(field,units,desc,p_pa(k))
    d2d = u(:,:,k)
    WRITE ( output_unit ) d2d
    PRINT '(A,F9.1,A,F5.1,A,F5.1)', 'Level (Pa):', p_pa(k), ' Min: ', MINVAL(d2d),&
            ' Max: ', MAXVAL(d2d)
  ENDDO var_u

  ! Do v-component of wind
  field = 'V        '
  units = 'm s{-1}                  '
  desc = 'v-component of velocity, rotated to grid      '
  PRINT *, 'FIELD = ', field
  PRINT *, 'UNITS = ', units
  PRINT *, 'DESC =  ',desc
  var_v : DO k = 1 , z3 + 1
    IF ( ( p_pa(k) .GT. 100100 ) .AND. ( p_pa(k) .LT. 200000 ) ) THEN
      CYCLE var_v
    END IF
    CALL write_gribprep_header(field,units,desc,p_pa(k))
    d2d = v(:,:,k)
    WRITE ( output_unit ) d2d
    PRINT '(A,F9.1,A,F5.1,A,F5.1)', 'Level (Pa):', p_pa(k), ' Min: ', MINVAL(d2d),&
            ' Max: ', MAXVAL(d2d)
  ENDDO var_v

  ! Relative Humidity
  field = 'RH       '
  units = '%                        '
  desc  = 'Relative humidity                             '
  PRINT *, 'FIELD = ', field
  PRINT *, 'UNITS = ', units
  PRINT *, 'DESC =  ',desc
  var_rh : DO k = 1 , z3 + 1
    IF ( ( p_pa(k) .GT. 100100 ) .AND. ( p_pa(k) .LT. 200000 ) ) THEN
      CYCLE var_rh
    END IF
    CALL write_gribprep_header(field,units,desc,p_pa(k))
    d2d = rh(:,:,k)
    WRITE ( output_unit ) d2d
    PRINT '(A,F9.1,A,F5.1,A,F5.1)', 'Level (Pa):', p_pa(k), ' Min: ', MINVAL(d2d),&
            ' Max: ', MAXVAL(d2d)
  ENDDO var_rh

  ! Do the heights
  field = 'HGT      '
  units = 'm                        '
  desc  = 'Geopotential height                           '
  PRINT *, 'FIELD = ', field
  PRINT *, 'UNITS = ', units
  PRINT *, 'DESC =  ',desc
  var_ht : DO k = 1 , z3 + 1
    IF ( ( p_pa(k) .GT. 100100 ) .AND. ( p_pa(k) .LT. 200000 ) ) THEN
      CYCLE var_ht
    END IF
    CALL write_gribprep_header(field,units,desc,p_pa(k))
    d2d = ht(:,:,k)
    WRITE ( output_unit ) d2d
    PRINT '(A,F9.1,A,F8.1,A,F8.1)', 'Level (Pa):', p_pa(k), ' Min: ', MINVAL(d2d),&
            ' Max: ', MAXVAL(d2d)
  ENDDO var_ht

  ! Terrain height
  field = 'SOILHGT '
  units = 'm                        '
  desc  = 'Height of topography                          '
  PRINT *, 'FIELD = ', field
  PRINT *, 'UNITS = ', units
  PRINT *, 'DESC =  ',desc
  CALL write_gribprep_header(field,units,desc,p_pa(z3+1))
  d2d = ht(:,:,z3+1)
  WRITE ( output_unit ) d2d
  PRINT '(A,F9.1,A,F9.1,A,F9.1)', 'Level (Pa):', p_pa(z3+1), ' Min: ', MINVAL(d2d),&
            ' Max: ', MAXVAL(d2d)

  ! Skin temperature
  field = 'SKINTEMP '
  units = 'K                        '
  desc  = 'Skin temperature                              '
  PRINT *, 'FIELD = ', field
  PRINT *, 'UNITS = ', units
  PRINT *, 'DESC =  ',desc
  CALL write_gribprep_header(field,units,desc,p_pa(z3+1))
  WRITE ( output_unit ) tskin
  PRINT '(A,F9.1,A,F9.1,A,F9.1)', 'Level (Pa):', p_pa(z3+1), &
       ' Min: ', MINVAL(tskin), ' Max: ', MAXVAL(tskin)

  ! Sea-level Pressure field
  field = 'PMSL     '
  units = 'Pa                       '
  desc  = 'Sea-level pressure                            '
  PRINT *, 'FIELD = ', field
  PRINT *, 'UNITS = ', units
  PRINT *, 'DESC =  ',desc
  CALL write_gribprep_header(field,units,desc,slp_level)
  WRITE ( output_unit ) slp
  PRINT '(A,F9.1,A,F9.1,A,F9.1)', 'Level (Pa):', slp_level, ' Min: ', MINVAL(slp),&
            ' Max: ', MAXVAL(slp)

  ! Surface Pressure field
  field = 'PSFC     '
  units = 'Pa                       '
  desc  = 'Surface pressure                              '
  PRINT *, 'FIELD = ', field
  PRINT *, 'UNITS = ', units
  PRINT *, 'DESC =  ',desc
  CALL write_gribprep_header(field,units,desc,p_pa(z3+1))
  WRITE ( output_unit ) psfc
  PRINT '(A,F9.1,A,F9.1,A,F9.1)', 'Level (Pa):', p_pa(z3+1), ' Min: ', MINVAL(psfc),&
            ' Max: ', MAXVAL(psfc)

  call get_r_missing_data(r_missing_data,istatus)
	print*,'output_gribprep_format r_missing_data ',r_missing_data
  if(istatus .ne. 1)then
      write(6,*)' Bad status for r_missing_data'
      stop
  endif

  IF (MINVAL(snocov).GE.0) THEN
    ! Water equivalent snow depth
    field = 'SNOWCOVR '
    units = '(DIMENSIONLESS)          '
    desc  = 'Snow cover flag                               '
    PRINT *, 'FIELD = ', field
    PRINT *, 'UNITS = ', units
    PRINT *, 'DESC =  ',desc
    CALL write_gribprep_header(field,units,desc,p_pa(z3+1))


!   Initialize output snow cover field to missing value
    d2d=-999.

!   Convert from fraction to mask using namelist entry snow_thresh
    if(snow_thresh .le. 1.0)then
        WHERE(snocov .ge. snow_thresh .AND. snocov .ne. r_missing_data) d2d = 1.0
        WHERE(snocov .lt. snow_thresh .AND. snocov .ne. r_missing_data) d2d = 0.0
    endif

    WRITE ( output_unit ) d2d
    PRINT '(A,F9.1,A,F9.2,A,F9.2)', 'Level (Pa):', p_pa(z3+1), &
        ' Min: ', MINVAL(d2d),&
        ' Max: ', MAXVAL(d2d) 

  ENDIF

  ! Get cloud species if this is a hot start
  IF (hotstart) THEN
    field = 'QLIQUID  '
    units = 'kg kg{-1}               '
    desc  = 'Cloud liquid water mixing ratio             '
    PRINT *, 'FIELD = ', field
    PRINT *, 'UNITS = ', units
    PRINT *, 'DESC =  ',desc    
    var_lwc : DO k = 1 , z3
      IF ( ( p_pa(k) .GT. 100100 ) .AND. ( p_pa(k) .LT. 200000 ) ) THEN
        CYCLE var_lwc
      END IF
      CALL write_gribprep_header(field,units,desc,p_pa(k))
      d2d = lwc(:,:,k)
      WRITE ( output_unit ) d2d
!HJ: W>=D+3. from F8.6 to F9.6 10/14/2013
      PRINT '(A,F9.1,A,F9.6,A,F9.6)', 'Level (Pa):', p_pa(k), ' Min: ', MINVAL(d2d),&
            ' Max: ', MAXVAL(d2d)
    END DO var_lwc 

    ! Cloud ice   
    field = 'QICE     '
    units = 'kg kg{-1}               '
    desc  = 'Cloud ice mixing ratio                      '
    PRINT *, 'FIELD = ', field
    PRINT *, 'UNITS = ', units
    PRINT *, 'DESC =  ',desc    
    var_ice: DO k = 1 , z3
      IF ( ( p_pa(k) .GT. 100100 ) .AND. ( p_pa(k) .LT. 200000 ) ) THEN
        CYCLE var_ice
      END IF
      CALL write_gribprep_header(field,units,desc,p_pa(k))
      d2d = ice(:,:,k)
      WRITE ( output_unit ) d2d
!HJ: W>=D+3. from F8.6 to F9.6 10/14/2013
      PRINT '(A,F9.1,A,F9.6,A,F9.6)', 'Level (Pa):', p_pa(k), ' Min: ', MINVAL(d2d),&
            ' Max: ', MAXVAL(d2d)
    END DO var_ice

    ! Cloud rain
    field = 'QRAIN    '
    units = 'kg kg{-1}               '
    desc  = 'Rain water mixing ratio                     '
    PRINT *, 'FIELD = ', field
    PRINT *, 'UNITS = ', units
    PRINT *, 'DESC =  ',desc    
    var_rai: DO k = 1 , z3
      IF ( ( p_pa(k) .GT. 100100 ) .AND. ( p_pa(k) .LT. 200000 ) ) THEN
        CYCLE var_rai
      END IF
      CALL write_gribprep_header(field,units,desc,p_pa(k))
      d2d = rai(:,:,k)
      WRITE ( output_unit ) d2d
!HJ: W>=D+3. from F8.6 to F9.6 10/14/2013
      PRINT '(A,F9.1,A,F9.6,A,F9.6)', 'Level (Pa):', p_pa(k), ' Min: ', MINVAL(d2d),&
            ' Max: ', MAXVAL(d2d)
    END DO var_rai

   ! Snow
    field = 'QSNOW    '
    units = 'kg kg{-1}               '
    desc  = 'Snow mixing ratio                           '
    PRINT *, 'FIELD = ', field
    PRINT *, 'UNITS = ', units
    PRINT *, 'DESC =  ',desc    
    var_sno: DO k = 1 , z3
      IF ( ( p_pa(k) .GT. 100100 ) .AND. ( p_pa(k) .LT. 200000 ) ) THEN
        CYCLE var_sno
      END IF
      CALL write_gribprep_header(field,units,desc,p_pa(k))
      d2d = sno(:,:,k)
      WRITE ( output_unit ) d2d
!HJ: W>=D+3. from F8.6 to F9.6 10/14/2013
      PRINT '(A,F9.1,A,F9.6,A,F9.6)', 'Level (Pa):', p_pa(k), ' Min: ', MINVAL(d2d),&
            ' Max: ', MAXVAL(d2d)
    END DO var_sno

    ! Graupel
    field = 'QGRAUPEL '
    units = 'kg kg{-1}               '
    desc  = 'Graupel mixing ratio                        '
    PRINT *, 'FIELD = ', field
    PRINT *, 'UNITS = ', units
    PRINT *, 'DESC =  ',desc    
    var_pic: DO k = 1 , z3
      IF ( ( p_pa(k) .GT. 100100 ) .AND. ( p_pa(k) .LT. 200000 ) ) THEN
        CYCLE var_pic
      END IF
      CALL write_gribprep_header(field,units,desc,p_pa(k))
      d2d = pic(:,:,k)
      WRITE ( output_unit ) d2d
!HJ: W>=D+3. from F8.6 to F9.6 10/14/2013
      PRINT '(A,F9.1,A,F9.6,A,F9.6)', 'Level (Pa):', p_pa(k), ' Min: ', MINVAL(d2d),&
            ' Max: ', MAXVAL(d2d)
    END DO var_pic

  ENDIF

  CLOSE (output_unit)
  DEALLOCATE (d2d)
  DEALLOCATE (p_pa)
  RETURN
  END SUBROUTINE output_gribprep_format
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE write_gribprep_header(field,units,desc,level)
 
  ! Writes the gribprep header given the filed, units, description, and level

  IMPLICIT NONE
  CHARACTER(LEN=9), INTENT(IN)  :: field
  CHARACTER(LEN=25),INTENT(IN)  :: units
  CHARACTER(LEN=46),INTENT(IN)  :: desc
  REAL, INTENT(IN)              :: level
  
  WRITE ( output_unit ) gp_version
  WRITE ( output_unit ) hdate,xfcst,source,field,units,desc,level,x,y,llflag
  SELECT CASE (llflag)
    CASE(1)
      WRITE ( output_unit ) knownloc,la1,lo1,dx,dy,latin1
    CASE(3)
      WRITE ( output_unit ) knownloc,la1,lo1,dx,dy,lov,latin1,latin2
    CASE(5)
      WRITE ( output_unit ) knownloc,la1,lo1,dx,dy,lov,latin1
  END SELECT

  END SUBROUTINE write_gribprep_header
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE output_metgrid_format(p, t, ht, u, v, w, rh, slp, psfc, &
                               lwc, rai, sno, ice, pic, snocov,tskin)

  !  Subroutine of lapsprep that will build a file the
  !  WPS format that can be read by metgrid

  IMPLICIT NONE

  ! Arguments

  REAL, INTENT(IN)                   :: p(:)        ! Pressure (hPa)
  REAL, INTENT(IN)                   :: t(:,:,:)    ! Temperature (K)
  REAL, INTENT(IN)                   :: ht(:,:,:)   ! Height (m)
  REAL, INTENT(IN)                   :: u(:,:,:)    ! U-wind (m s{-1})
  REAL, INTENT(IN)                   :: v(:,:,:)    ! V-wind (m s{-1})
  REAL, INTENT(IN)                   :: w(:,:,:)    ! W-wind (m s{-1})
  REAL, INTENT(IN)                   :: rh(:,:,:)   ! Relative Humidity (%)
  REAL, INTENT(IN)                   :: slp(:,:)    ! Sea-level Pressure (Pa)
  REAL, INTENT(IN)                   :: psfc(:,:)   ! Surface Pressure (Pa)
  REAL, INTENT(IN)                   :: lwc(:,:,:)  ! Cloud liquid (kg/kg)
  REAL, INTENT(IN)                   :: rai(:,:,:)  ! Rain (kg/kg)
  REAL, INTENT(IN)                   :: sno(:,:,:)  ! Snow (kg/kg)
  REAL, INTENT(IN)                   :: ice(:,:,:)  ! Ice (kg/kg)
  REAL, INTENT(IN)                   :: pic(:,:,:)  ! Graupel (kg/kg)
  REAL, INTENT(IN)                   :: snocov(:,:) ! Snow cover (fract)
  REAL, INTENT(IN)                   :: tskin(:,:)  ! Skin temperature

  ! Local Variables

  INTEGER            :: valid_mm, valid_dd
  CHARACTER (LEN=256):: output_file_name
  REAL, ALLOCATABLE  :: d2d(:,:)
  REAL, ALLOCATABLE  :: p_pa(:)
  INTEGER            :: k,yyyyddd
  INTEGER            :: istatus
  REAL               :: r_missing_data

  ! Allocate a scratch 2d array
  ALLOCATE (d2d (x,y) )
  ALLOCATE (p_pa (z3+1))
  ! Build the output file name

  output_prefix = TRIM(laps_data_root)// '/lapsprd/lapsprep/wps/LAPS'
  yyyyddd = valid_yyyy*1000 + valid_jjj
  CALL wrf_date_to_ymd(yyyyddd, valid_yyyy, valid_mm, valid_dd)
  WRITE(hdate, '(I4.4,"-",I2.2,"-",I2.2,"_",I2.2,":",I2.2,":00.0000")') &
          valid_yyyy, valid_mm, valid_dd, valid_hh, valid_min
  IF (valid_min .EQ. 0) THEN
    output_file_name = TRIM(output_prefix) // ':' // hdate(1:13)
  ELSE
    output_file_name = TRIM(output_prefix) // ':' // hdate(1:16)
  ENDIF
  !  Open the file for sequential, unformatted output
  OPEN ( FILE   = TRIM(output_file_name)    , &
         UNIT   = output_unit        , &
         FORM   = 'UNFORMATTED' , &
         STATUS = 'UNKNOWN'     , &
         ACCESS = 'SEQUENTIAL'    )

  ! Convert p levels from mb to Pascals

  p_pa = p * 100.

  ! Set llflag based on grid type

  IF      ( grid_type(1:8)  .EQ. 'mercator'                 ) THEN
    llflag = 1
  ELSE IF ( ( grid_type(1:24) .EQ. 'secant lambert conformal' ) .or. &
           ( grid_type(1:28) .EQ.  'tangential lambert conformal' ) )THEN
    llflag = 3
  ELSE IF ( grid_type(1:19) .EQ. 'polar stereographic'      ) THEN
    llflag = 5
  ELSE
    PRINT '(A,A,A)','Unknown map projection: ',TRIM(grid_type),'.  I quit.'
    STOP 'unknown_projection'
  END IF

  PRINT *, 'METGRID VERSION =', gp_version_wps
  PRINT *, 'SOURCE = ', source
  PRINT *, 'HDATE = ',hdate
  PRINT *, 'XFCST = ', xfcst
  PRINT *, 'NX = ', X
  PRINT *, 'NY = ', Y
  PRINT *, 'IPROJ = ', LLFLAG
  PRINT *, 'KNOWNLOC = ', knownloc
  PRINT *, 'STARTLAT = ',LA1
  PRINT *, 'STARTLON = ',LO1
  PRINT *, 'DX = ', DX
  PRINT *, 'DY = ', DY
  PRINT *, 'XLONC = ', LOV
  PRINT *, 'TRUELAT1 = ', LATIN1
  PRINT *, 'TRUELAT2 = ', LATIN2

 ! Output temperature
  field = 'TT       '
  units = 'K                        '
  desc  = 'Temperature                                   '
  PRINT *, 'FIELD = ', field
  PRINT *, 'UNITS = ', units
  PRINT *, 'DESC =  ',desc
  var_t : DO k = 1 , z3 + 1
    IF (( p_pa(k) .GT. 100100).AND.(p_pa(k).LT.200000) ) THEN
      CYCLE var_t
    ENDIF
    CALL write_metgrid_header(field,units,desc,p_pa(k))
    d2d = t(:,:,k)
    WRITE ( output_unit ) d2d
    PRINT '(A,F9.1,A,F5.1,A,F5.1)','Level (Pa):',p_pa(k),' Min: ', &
           MINVAL(d2d),' Max: ', MAXVAL(d2d)
  ENDDO var_t

  ! Do u-component of wind
  field = 'UU       '
  units = 'm s-1                    '
  desc = 'U                                             '
  PRINT *, 'FIELD = ', field
  PRINT *, 'UNITS = ', units
  PRINT *, 'DESC =  ',desc
  var_u : DO k = 1 , z3 + 1
    IF ( ( p_pa(k) .GT. 100100 ) .AND. ( p_pa(k) .LT. 200000 ) ) THEN
      CYCLE var_u
    END IF
    CALL write_metgrid_header(field,units,desc,p_pa(k))
    d2d = u(:,:,k)
    WRITE ( output_unit ) d2d
    PRINT '(A,F9.1,A,F5.1,A,F5.1)', 'Level (Pa):', p_pa(k), ' Min: ', MINVAL(d2d),&
            ' Max: ', MAXVAL(d2d)
  ENDDO var_u

  ! Do v-component of wind
  field = 'VV       '
  units = 'm s-1                    '
  desc = 'V                                             '
  PRINT *, 'FIELD = ', field
  PRINT *, 'UNITS = ', units
  PRINT *, 'DESC =  ',desc
  var_v : DO k = 1 , z3 + 1
    IF ( ( p_pa(k) .GT. 100100 ) .AND. ( p_pa(k) .LT. 200000 ) ) THEN
      CYCLE var_v
    END IF
    CALL write_metgrid_header(field,units,desc,p_pa(k))
    d2d = v(:,:,k)
    WRITE ( output_unit ) d2d
    PRINT '(A,F9.1,A,F5.1,A,F5.1)', 'Level (Pa):', p_pa(k), ' Min: ', MINVAL(d2d),&
            ' Max: ', MAXVAL(d2d)
  ENDDO var_v

  ! Do w-component of wind
  IF (use_laps_vv) THEN
    field = 'VVEL     '
    units = 'm s-1                   '
    desc = 'W                                             '
    PRINT *, 'FIELD = ', field
    PRINT *, 'UNITS = ', units
    PRINT *, 'DESC =  ',desc
    var_vv : DO k = 1 , z3 + 1
      IF ( ( p_pa(k) .GT. 100100 ) .AND. ( p_pa(k) .LT. 200000 ) ) THEN
        CYCLE var_vv
      END IF
      CALL write_metgrid_header(field,units,desc,p_pa(k))
      d2d = w(:,:,k)
      WRITE ( output_unit ) d2d
      PRINT '(A,F9.1,A,F9.1,A,F9.1)', 'Level (Pa):', p_pa(k), ' Min: ', MINVAL(d2d),&
              ' Max: ', MAXVAL(d2d)
  ENDDO var_vv
  ENDIF

  ! Relative Humidity
  field = 'RH       '
  units = '%                        '
  desc  = 'Relative Humidity                             '
  PRINT *, 'FIELD = ', field
  PRINT *, 'UNITS = ', units
  PRINT *, 'DESC =  ',desc
  var_rh : DO k = 1 , z3 + 1
    IF ( ( p_pa(k) .GT. 100100 ) .AND. ( p_pa(k) .LT. 200000 ) ) THEN
      CYCLE var_rh
    END IF
    CALL write_metgrid_header(field,units,desc,p_pa(k))
    d2d = rh(:,:,k)
    WRITE ( output_unit ) d2d
    PRINT '(A,F9.1,A,F5.1,A,F5.1)', 'Level (Pa):', p_pa(k), ' Min: ', MINVAL(d2d),&
            ' Max: ', MAXVAL(d2d)
  ENDDO var_rh

  ! Do the heights
  field = 'HGT      '
  units = 'm                        '
  desc  = 'Geopotential height                           '
  PRINT *, 'FIELD = ', field
  PRINT *, 'UNITS = ', units
  PRINT *, 'DESC =  ',desc
  var_ht : DO k = 1 , z3 + 1
    IF ( ( p_pa(k) .GT. 100100 ) .AND. ( p_pa(k) .LT. 200000 ) ) THEN
      CYCLE var_ht
    END IF
    CALL write_metgrid_header(field,units,desc,p_pa(k))
    d2d = ht(:,:,k)
    WRITE ( output_unit ) d2d
    PRINT '(A,F9.1,A,F8.1,A,F8.1)', 'Level (Pa):', p_pa(k), ' Min: ', MINVAL(d2d),&
            ' Max: ', MAXVAL(d2d)
  ENDDO var_ht

  ! Terrain height
  field = 'SOILHGT '
  units = 'm                        '
  desc  = 'Terrain field of source analysis              '
  PRINT *, 'FIELD = ', field
  PRINT *, 'UNITS = ', units
  PRINT *, 'DESC =  ',desc
  CALL write_metgrid_header(field,units,desc,p_pa(z3+1))
  d2d = ht(:,:,z3+1)
  WRITE ( output_unit ) d2d
  PRINT '(A,F9.1,A,F9.1,A,F9.1)', 'Level (Pa):', p_pa(z3+1), ' Min: ', MINVAL(d2d),&
            ' Max: ', MAXVAL(d2d)

  ! Skin temperature
  IF (use_laps_skintemp) THEN
      field = 'SKINTEMP '
      units = 'K                        '
      desc  = 'Skin temperature                              '
      PRINT *, 'FIELD = ', field
      PRINT *, 'UNITS = ', units
      PRINT *, 'DESC =  ',desc
      CALL write_metgrid_header(field,units,desc,p_pa(z3+1))
      WRITE ( output_unit ) tskin
      PRINT '(A,F9.1,A,F9.1,A,F9.1)', 'Level (Pa):', p_pa(z3+1), &
         ' Min: ', MINVAL(tskin), ' Max: ', MAXVAL(tskin)
  ENDIF

  ! Sea-level Pressure field
  field = 'PMSL     '
  units = 'Pa                       '
  desc  = 'Sea-level Pressure                            '
  PRINT *, 'FIELD = ', field
  PRINT *, 'UNITS = ', units
  PRINT *, 'DESC =  ',desc
  CALL write_metgrid_header(field,units,desc,slp_level)
  WRITE ( output_unit ) slp
  PRINT '(A,F9.1,A,F9.1,A,F9.1)', 'Level (Pa):', slp_level, ' Min: ', MINVAL(slp),&
            ' Max: ', MAXVAL(slp)

  ! Surface Pressure field
  field = 'PSFC     '
  units = 'Pa                       '
  desc  = 'Surface Pressure                              '
  PRINT *, 'FIELD = ', field
  PRINT *, 'UNITS = ', units
  PRINT *, 'DESC =  ',desc
  CALL write_metgrid_header(field,units,desc,p_pa(z3+1))
  WRITE ( output_unit ) psfc
  PRINT '(A,F9.1,A,F9.1,A,F9.1)', 'Level (Pa):', p_pa(z3+1), ' Min: ', MINVAL(psfc),&
            ' Max: ', MAXVAL(psfc)

  call get_r_missing_data(r_missing_data,istatus)
	print*,'output_metgrid_format r_missing_data ',r_missing_data
  if(istatus .ne. 1)then
      write(6,*)' Bad status for r_missing_data'
      stop
  endif

  ! Get cloud species if this is a hot start
  IF (hotstart) THEN
    field = 'QC       '
    units = 'kg kg-1                 '
    desc  = 'Cloud liquid water mixing ratio             '
    PRINT *, 'FIELD = ', field
    PRINT *, 'UNITS = ', units
    PRINT *, 'DESC =  ',desc
    var_lwc : DO k = 1 , z3
      IF ( ( p_pa(k) .GT. 100100 ) .AND. ( p_pa(k) .LT. 200000 ) ) THEN
        CYCLE var_lwc
      END IF
      CALL write_metgrid_header(field,units,desc,p_pa(k))
      d2d = lwc(:,:,k)
      WRITE ( output_unit ) d2d
!HJ: W>=D+3. 10/14/2013
      PRINT '(A,F9.1,A,F9.6,A,F9.6)', 'Level (Pa):', p_pa(k), ' Min: ', MINVAL(d2d),&
            ' Max: ', MAXVAL(d2d)
    END DO var_lwc

    ! Cloud ice
    field = 'QI       '
    units = 'kg kg-1                 '
    desc  = 'Cloud ice mixing ratio                      '
    PRINT *, 'FIELD = ', field
    PRINT *, 'UNITS = ', units
    PRINT *, 'DESC =  ',desc
    var_ice: DO k = 1 , z3
      IF ( ( p_pa(k) .GT. 100100 ) .AND. ( p_pa(k) .LT. 200000 ) ) THEN
        CYCLE var_ice
      END IF
      CALL write_metgrid_header(field,units,desc,p_pa(k))
      d2d = ice(:,:,k)
      WRITE ( output_unit ) d2d
!HJ: W>=D+3. 10/14/2013
      PRINT '(A,F9.1,A,F9.6,A,F9.6)', 'Level (Pa):', p_pa(k), ' Min: ', MINVAL(d2d),&
            ' Max: ', MAXVAL(d2d)
    END DO var_ice


    ! Cloud rain
    field = 'QR       '
    units = 'kg kg-1                 '
    desc  = 'Rain water mixing ratio                     '
    PRINT *, 'FIELD = ', field
    PRINT *, 'UNITS = ', units
    PRINT *, 'DESC =  ',desc
    var_rai: DO k = 1 , z3
      IF ( ( p_pa(k) .GT. 100100 ) .AND. ( p_pa(k) .LT. 200000 ) ) THEN
        CYCLE var_rai
      END IF
      CALL write_metgrid_header(field,units,desc,p_pa(k))
      d2d = rai(:,:,k)
      WRITE ( output_unit ) d2d
!HJ: W>=D+3. 10/14/2013
      PRINT '(A,F9.1,A,F9.6,A,F9.6)', 'Level (Pa):', p_pa(k), ' Min: ', MINVAL(d2d),&
            ' Max: ', MAXVAL(d2d)
    END DO var_rai

   ! Snow
    field = 'QS       '
    units = 'kg kg-1                 '
    desc  = 'Snow mixing ratio                           '
    PRINT *, 'FIELD = ', field
    PRINT *, 'UNITS = ', units
    PRINT *, 'DESC =  ',desc
    var_sno: DO k = 1 , z3
      IF ( ( p_pa(k) .GT. 100100 ) .AND. ( p_pa(k) .LT. 200000 ) ) THEN
        CYCLE var_sno
      END IF
      CALL write_metgrid_header(field,units,desc,p_pa(k))
      d2d = sno(:,:,k)
      WRITE ( output_unit ) d2d
!HJ: W>=D+3. 10/14/2013
      PRINT '(A,F9.1,A,F9.6,A,F9.6)', 'Level (Pa):', p_pa(k), ' Min: ', MINVAL(d2d),&
            ' Max: ', MAXVAL(d2d)
    END DO var_sno

    ! Graupel
    field = 'QG       '
    units = 'kg kg-1                 '
    desc  = 'Graupel mixing ratio                        '
    PRINT *, 'FIELD = ', field
    PRINT *, 'UNITS = ', units
    PRINT *, 'DESC =  ',desc
    var_pic: DO k = 1 , z3
      IF ( ( p_pa(k) .GT. 100100 ) .AND. ( p_pa(k) .LT. 200000 ) ) THEN
        CYCLE var_pic
      END IF
      CALL write_metgrid_header(field,units,desc,p_pa(k))
      d2d = pic(:,:,k)
      WRITE ( output_unit ) d2d
!HJ: W>=D+3. 10/14/2013
      PRINT '(A,F9.1,A,F9.6,A,F9.6)', 'Level (Pa):', p_pa(k), ' Min: ', MINVAL(d2d),&
            ' Max: ', MAXVAL(d2d)
    END DO var_pic

  ENDIF

  CLOSE (output_unit)
  DEALLOCATE (d2d)
  DEALLOCATE (p_pa)
  RETURN
  END SUBROUTINE output_metgrid_format
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE write_metgrid_header(field,units,desc,level)

  ! Writes the gribprep header given the filed, units, description, and level

  IMPLICIT NONE
  CHARACTER(LEN=9), INTENT(IN)  :: field
  CHARACTER(LEN=25),INTENT(IN)  :: units
  CHARACTER(LEN=46),INTENT(IN)  :: desc
  REAL, INTENT(IN)              :: level
  LOGICAL, PARAMETER            :: is_wind_grid_rel = .true.
  REAL                          :: radius_of_earth_m,radius_of_earth_km
  INTEGER                       :: istatus
 
  call get_earth_radius(radius_of_earth_m,istatus)
  radius_of_earth_km = radius_of_earth_m / 1000.
  print*,'write_metgrid_header radius_of_earth_km ',radius_of_earth_km
  if(istatus .ne. 1)then
      write(6,*)' Bad status for radius_of_earth_km'
      stop
  endif

  WRITE ( output_unit ) gp_version_wps
  WRITE ( output_unit ) hdate,xfcst,source,field,units,desc,level,x,y,llflag
  SELECT CASE (llflag)
    CASE(1)
      WRITE ( output_unit ) knownloc,la1,lo1,dx,dy,latin1,radius_of_earth_km
    CASE(3)
      WRITE ( output_unit ) knownloc,la1,lo1,dx,dy,lov,latin1,latin2,radius_of_earth_km
    CASE(5)
      WRITE ( output_unit ) knownloc,la1,lo1,dx,dy,lov,latin1,radius_of_earth_km
  END SELECT
  WRITE ( output_unit) is_wind_grid_rel

  END SUBROUTINE write_metgrid_header

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE lapsprep_wrf
  
