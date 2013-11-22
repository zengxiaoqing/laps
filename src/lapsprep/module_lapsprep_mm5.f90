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



MODULE lapsprep_mm5

! PURPOSE
! =======
! Module to contain the various output routines needed for lapsprep
! to support initializition of the MM5v3 model.
!
! SUBROUTINES CONTAINED
! =====================
! output_pregrid_format  - Used to support MM5 initializations
! write_pregrid_header   - Writes pregrid headers
!
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
  IMPLICIT NONE

  PRIVATE
  INTEGER, PARAMETER :: pg_version = 3
  REAL, PARAMETER    :: xfcst = 0.
  CHARACTER (LEN=24) :: hdate
  INTEGER            :: llflag
  CHARACTER (LEN=9)  :: field
  CHARACTER (LEN=25) :: units
  CHARACTER (LEN=46) :: desc 
  INTEGER, PARAMETER :: output_unit = 78
  REAL, PARAMETER    :: slp_level = 201300.0

  PUBLIC output_pregrid_format
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE output_pregrid_format(p, t, ht, u, v, rh, slp, &
                               lwc, rai, sno, ice, pic, snocov, tskin)

  !  Subroutine of lapsprep that will build a file in the
  !  MM5v3 pregrid format that can be read by REGRID

  IMPLICIT NONE

  ! Arguments
 
  REAL, INTENT(IN)                   :: p(:)        ! Pressure (hPa)
  REAL, INTENT(IN)                   :: t(:,:,:)    ! Temperature (K)
  REAL, INTENT(IN)                   :: ht(:,:,:)   ! Height (m)
  REAL, INTENT(IN)                   :: u(:,:,:)    ! U-wind (m s{-1})
  REAL, INTENT(IN)                   :: v(:,:,:)    ! V-wind (m s{-1})
  REAL, INTENT(IN)                   :: rh(:,:,:)   ! Relative Humidity (%)
  REAL, INTENT(IN)                   :: slp(:,:)    ! Sea-level Pressure (Pa)
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
  CHARACTER(LEN=256)  :: sstfile
  REAL, ALLOCATABLE  :: d2d(:,:)
  REAL, ALLOCATABLE  :: p_pa(:)
  INTEGER            :: k,yyyyddd
 
  ! Allocate a scratch 2d array
  ALLOCATE (d2d (x,y) )
  ALLOCATE (p_pa (z3+1))
  ! Build the output file name
 
  output_prefix = TRIM(laps_data_root)// '/lapsprd/lapsprep/mm5/LAPS'
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
            ( grid_type(1:28) .EQ. 'tangential lambert conformal') ) THEN
    llflag = 3
  ELSE IF ( grid_type(1:19) .EQ. 'polar stereographic'      ) THEN
    llflag = 5
  ELSE
    PRINT '(A,A,A)','Unknown map projection: ',TRIM(grid_type),'.  I quit.'
    STOP 'unknown_projection'
  END IF

  PRINT *, 'PREGRID VERSION =', pg_version
  PRINT *, 'HDATE = ',hdate
  PRINT *, 'XFCST = ', xfcst 
  PRINT *, 'NX = ', X
  PRINT *, 'NY = ', Y
  PRINT *, 'IPROJ = ', LLFLAG
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
    CALL write_pregrid_header(field,units,desc,p_pa(k))
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
    CALL write_pregrid_header(field,units,desc,p_pa(k))
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
    CALL write_pregrid_header(field,units,desc,p_pa(k))
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
    CALL write_pregrid_header(field,units,desc,p_pa(k))
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
    CALL write_pregrid_header(field,units,desc,p_pa(k))
    d2d = ht(:,:,k)
    WRITE ( output_unit ) d2d
    PRINT '(A,F9.1,A,F8.1,A,F8.1)', 'Level (Pa):', p_pa(k), ' Min: ', MINVAL(d2d),&
            ' Max: ', MAXVAL(d2d)
  ENDDO var_ht

  ! Terrain height
  field = 'HGT      '
  units = 'm                        '
  desc  = 'Height of topography                          '
  PRINT *, 'FIELD = ', field
  PRINT *, 'UNITS = ', units
  PRINT *, 'DESC =  ',desc
  CALL write_pregrid_header(field,units,desc,p_pa(z3+1))
  d2d = ht(:,:,z3+1)
  WRITE ( output_unit ) d2d
  PRINT '(A,F9.1,A,F9.1,A,F9.1)', 'Level (Pa):', p_pa(z3+1), ' Min: ', MINVAL(d2d),&
            ' Max: ', MAXVAL(d2d)

  ! Sea-level Pressure field
  field = 'PMSL     '
  units = 'Pa                       '
  desc  = 'Sea-level pressure                            '
  PRINT *, 'FIELD = ', field
  PRINT *, 'UNITS = ', units
  PRINT *, 'DESC =  ',desc
  CALL write_pregrid_header(field,units,desc,slp_level)
  WRITE ( output_unit ) slp
  PRINT '(A,F9.1,A,F9.1,A,F9.1)', 'Level (Pa):', slp_level, ' Min: ', MINVAL(slp),&
            ' Max: ', MAXVAL(slp)


  ! Snow cover
  IF ((MINVAL(snocov) .GE. 0.).AND.(MAXVAL(snocov) .LT. 1.1)) THEN
    ! Water equivalent snow depth
    field = 'SNOWCOVR '
    units = '(DIMENSIONLESS)          '
    desc  = 'Snow cover flag                               '
    PRINT *, 'FIELD = ', field
    PRINT *, 'UNITS = ', units
    PRINT *, 'DESC =  ',desc
    CALL write_pregrid_header(field,units,desc,p_pa(z3+1))

    ! Convert the snow cover fraction (ranges from 0 -> 1.) into a flag
    ! value using the snow_thresh parameter set in the namelist.
    d2d(:,:) = 0.
    WHERE(snocov .GE. snow_thresh) d2d = 1.

    WRITE ( output_unit ) d2d
    PRINT '(A,F9.1,A,F9.2,A,F9.2)', 'Level (Pa):',p_pa(z3+1), &
       ' Min: ', MINVAL(d2d),&
       ' Max: ', MAXVAL(d2d) 

  ENDIF
  
  ! Get cloud species if this is a hot start
  IF (hotstart) THEN
    field = 'CLW      '
    units = 'kg kg{-1}               '
    desc  = 'Cloud liquid water mixing ratio             '
    PRINT *, 'FIELD = ', field
    PRINT *, 'UNITS = ', units
    PRINT *, 'DESC =  ',desc    
    var_lwc : DO k = 1 , z3
      IF ( ( p_pa(k) .GT. 100100 ) .AND. ( p_pa(k) .LT. 200000 ) ) THEN
        CYCLE var_lwc
      END IF
      CALL write_pregrid_header(field,units,desc,p_pa(k))

      d2d = lwc(:,:,k) 
      
      WRITE ( output_unit ) d2d
!HJ: W>=D+3. from F8.6 to F9.6 10/14/2013
      PRINT '(A,F9.1,A,F9.6,A,F9.6)', 'Level (Pa):', p_pa(k), ' Min: ', MINVAL(d2d),&
            ' Max: ', MAXVAL(d2d)
    END DO var_lwc 

    ! Cloud ice   
    field = 'ICE      '
    units = 'kg kg{-1}               '
    desc  = 'Cloud ice mixing ratio                      '
    PRINT *, 'FIELD = ', field
    PRINT *, 'UNITS = ', units
    PRINT *, 'DESC =  ',desc    
    var_ice: DO k = 1 , z3
      IF ( ( p_pa(k) .GT. 100100 ) .AND. ( p_pa(k) .LT. 200000 ) ) THEN
        CYCLE var_ice
      END IF
      CALL write_pregrid_header(field,units,desc,p_pa(k))
      d2d = ice(:,:,k)
      WRITE ( output_unit ) d2d
!HJ: W>=D+3. from F8.6 to F9.6 10/14/2013
      PRINT '(A,F9.1,A,F9.6,A,F9.6)', 'Level (Pa):', p_pa(k), ' Min: ', MINVAL(d2d),&
            ' Max: ', MAXVAL(d2d)
    END DO var_ice

    ! Cloud rain
    field = 'RNW      '
    units = 'kg kg{-1}               '
    desc  = 'Rain water mixing ratio                     '
    PRINT *, 'FIELD = ', field
    PRINT *, 'UNITS = ', units
    PRINT *, 'DESC =  ',desc    
    var_rai: DO k = 1 , z3
      IF ( ( p_pa(k) .GT. 100100 ) .AND. ( p_pa(k) .LT. 200000 ) ) THEN
        CYCLE var_rai
      END IF
      CALL write_pregrid_header(field,units,desc,p_pa(k))
      d2d = rai(:,:,k)
      WRITE ( output_unit ) d2d
!HJ: W>=D+3. from F8.6 to F9.6 10/14/2013
      PRINT '(A,F9.1,A,F9.6,A,F9.6)', 'Level (Pa):', p_pa(k), ' Min: ', MINVAL(d2d),&
            ' Max: ', MAXVAL(d2d)
    END DO var_rai

   ! Snow
    field = 'SNOW     '
    units = 'kg kg{-1}               '
    desc  = 'Snow mixing ratio                           '
    PRINT *, 'FIELD = ', field
    PRINT *, 'UNITS = ', units
    PRINT *, 'DESC =  ',desc    
    var_sno: DO k = 1 , z3
      IF ( ( p_pa(k) .GT. 100100 ) .AND. ( p_pa(k) .LT. 200000 ) ) THEN
        CYCLE var_sno
      END IF
      CALL write_pregrid_header(field,units,desc,p_pa(k))
      d2d = sno(:,:,k)
      WRITE ( output_unit ) d2d
!HJ: W>=D+3. from F8.6 to F9.6 10/14/2013
      PRINT '(A,F9.1,A,F9.6,A,F9.6)', 'Level (Pa):', p_pa(k), ' Min: ', MINVAL(d2d),&
            ' Max: ', MAXVAL(d2d)
    END DO var_sno

    ! Graupel
    field = 'GRAUPEL  '
    units = 'kg kg{-1}               '
    desc  = 'Graupel mixing ratio                        '
    PRINT *, 'FIELD = ', field
    PRINT *, 'UNITS = ', units
    PRINT *, 'DESC =  ',desc    
    var_pic: DO k = 1 , z3
      IF ( ( p_pa(k) .GT. 100100 ) .AND. ( p_pa(k) .LT. 200000 ) ) THEN
        CYCLE var_pic
      END IF
      CALL write_pregrid_header(field,units,desc,p_pa(k))
      d2d = pic(:,:,k)
      WRITE ( output_unit ) d2d
!HJ: W>=D+3. from F8.6 to F9.6 10/14/2013
      PRINT '(A,F9.1,A,F9.6,A,F9.6)', 'Level (Pa):', p_pa(k), ' Min: ', MINVAL(d2d),&
            ' Max: ', MAXVAL(d2d)
    END DO var_pic

  ENDIF

  CLOSE (output_unit)

  ! Change for RSA to write skin temp out as a LAPS:TSKIN file
  ! Skin temperature field
  sstfile = TRIM(output_prefix) // ':TSKIN'
  OPEN(UNIT=output_unit, FILE=sstfile,FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
  field = 'SKINTEMP '
  units = 'K                        '
  desc  = 'Skin temperature                              '
  PRINT *, 'FIELD = ', field
  PRINT *, 'UNITS = ', units
  PRINT *, 'DESC =  ',desc
  CALL write_pregrid_header(field,units,desc,p_pa(z3+1))
  WRITE ( output_unit ) tskin
  PRINT '(A,F9.1,A,F9.1,A,F9.1)', 'Level (Pa):', p_pa(z3+1), ' Min: ', &
      MINVAL(tskin), ' Max: ', MAXVAL(tskin)
  CLOSE(output_unit)

  DEALLOCATE (d2d)
  DEALLOCATE (p_pa)
  RETURN
  END SUBROUTINE output_pregrid_format
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE write_pregrid_header(field,units,desc,level)
 
  ! Writes the gribprep header given the filed, units, description, and level

  IMPLICIT NONE
  CHARACTER(LEN=9), INTENT(IN)  :: field
  CHARACTER(LEN=25),INTENT(IN)  :: units
  CHARACTER(LEN=46),INTENT(IN)  :: desc
  REAL, INTENT(IN)              :: level
  
  WRITE ( output_unit ) pg_version
  WRITE ( output_unit ) hdate,xfcst,field,units,desc,level,x,y,llflag
  SELECT CASE (llflag)
    CASE(1)
      WRITE ( output_unit ) la1,lo1,dx,dy,latin1
    CASE(3)
      WRITE ( output_unit ) la1,lo1,dx,dy,lov,latin1,latin2
    CASE(5)
      WRITE ( output_unit ) la1,lo1,dx,dy,lov,latin1
  END SELECT

  END SUBROUTINE write_pregrid_header
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE lapsprep_mm5
  
