!dis    Forecast Systems Laboratory
!dis    NOAA/OAR/ERL/FSL
!dis    325 Broadway
!dis    Boulder, CO     80303
!dis
!dis    Forecast Research Division
!dis    Local Analysis and Prediction Branch
!dis    LAPS
!dis
!dis    This software and its documentation are in the public domain and
!dis    are furnished "as is."  The United States government, its
!dis    instrumentalities, officers, employees, and agents make no
!dis    warranty, express or implied, as to the usefulness of the software
!dis    and documentation for any purpose.  They assume no responsibility
!dis    (1) for the use of the software and documentation; or (2) to provide
!dis     technical support to users.
!dis
!dis    Permission to use, copy, modify, and distribute this software is
!dis    hereby granted, provided that the entire disclaimer notice appears
!dis    in all copies.  All modifications to this software must be clearly
!dis    documented, and are solely the responsibility of the agent making
!dis    the modifications.  If significant modifications or enhancements
!dis    are made to this software, the FSL Software Policy Manager
!dis    (softwaremgr@fsl.noaa.gov) should be notified.
!dis
!dis
!dis
!dis
!dis
!dis
!dis
!
MODULE lapsprep_output

! PURPOSE
! =======
! Module to contain the various output routines needed for lapsprep.
!
! SUBROUTINES CONTAINED
! =====================
! output_pregrid_v3  - Used to support MM5v3 initializations
!
! REMARKS
! =======
! For new models, add output routines to this module.
!
! HISTORY
! =======
! 28 Nov 2000 -- Original -- Brent Shaw

  USE setup
  USE laps_static
  USE date_pack

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE output_pregrid_v3(p, t, ht, u, v, rh, slp, &
                               lwc, rai, sno, ice, pic)

  !  Subroutine of lapsprep that will build a file in MM5 Pregrid v3 format
  !  to allow a LAPS initialization of MM5.  

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

  ! Local Variables

  CHARACTER (LEN=24)  :: hdate
  INTEGER :: llflag
  CHARACTER (LEN=9) :: field
  CHARACTER (LEN=25) :: units
  CHARACTER (LEN=46) :: desc
  INTEGER            :: valid_mm, valid_dd
  CHARACTER (LEN=256):: output_file_name
  REAL, ALLOCATABLE  :: d2d(:,:)
  REAL, ALLOCATABLE  :: p_pa(:)
  INTEGER, PARAMETER :: output_unit = 78
  INTEGER            :: k,yyyyddd
 
  ! Allocate a scratch 2d array
  ALLOCATE (d2d (x,y) )
  ALLOCATE (p_pa (z3+1))
  ! Build the output file name
 
  output_prefix = TRIM(laps_data_root)// '/lapsprd/dprep/' // &
                  TRIM(output_format) // '_INIT'
  yyyyddd = valid_yyyy*1000 + valid_jjj
  CALL wrf_date_to_ymd(yyyyddd, valid_yyyy, valid_mm, valid_dd) 
  WRITE(hdate, '(I4.4,"-",I2.2,"-",I2.2,"_",I2.2,":",I2.2,":00.0000")') &
          valid_yyyy, valid_mm, valid_dd, valid_hh, valid_min
  output_file_name = TRIM(output_prefix) // ':' // hdate(1:13)

  !  Open the file for sequential, unformatted output
  OPEN ( FILE   = TRIM(output_file_name)    , &
         UNIT   = output_unit        , &
         FORM   = 'UNFORMATTED' , &
         STATUS = 'UNKNOWN'     , &
         ACCESS = 'SEQUENTIAL'    )

  ! Convert p levels from mb to Pascals

  p_pa = p * 1000.

  ! Set llflag based on grid type

  IF      ( grid_type(1:8)  .EQ. 'mercator'                 ) THEN
    llflag = 1
  ELSE IF ( grid_type(1:24) .EQ. 'secant lambert conformal' ) THEN
    llflag = 3
  ELSE IF ( grid_type(1:19) .EQ. 'polar stereographic'      ) THEN
    llflag = 5
  ELSE
    PRINT '(A,A,A)','Unknown map projection: ',TRIM(grid_type),'.  I quit.'
    STOP 'unknown_projection'
  END IF

  PRINT *, 'HDATE = ',hdate
  PRINT *, 'XFCST = 0.' 
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
    WRITE (output_unit ) 3
    WRITE ( output_unit ) hdate , 0. ,  field , units , desc , &
                          p_pa(k), x, y, llflag
    IF (llflag .EQ. 1) THEN
      WRITE ( output_unit ) la1 , lo1 , dx , dy , latin1
    ELSE IF ( llflag .EQ. 3 ) THEN
      WRITE ( output_unit ) la1 , lo1 , dx , dy , lov , latin1 , latin2
    ELSE IF ( llflag .EQ. 5 ) THEN
      WRITE ( output_unit ) la1 , lo1 , dx , dy , lov , latin1
    END IF
    d2d = t(:,:,k)
    WRITE ( output_unit ) d2d
    PRINT '(A,F9.1,A,F5.1,A,F5.1)', 'Level (Pa):', p_pa(k), ' Min: ', MINVAL(d2d),&
            ' Max: ', MAXVAL(d2d)
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
    WRITE ( output_unit ) 3
    WRITE ( output_unit ) hdate , 0. ,  field , units , desc , &
                          p_pa(k) , x , y , llflag
    IF      ( llflag .EQ. 1 ) THEN
      WRITE ( output_unit ) la1 , lo1 , dx , dy , latin1
    ELSE IF ( llflag .EQ. 3 ) THEN
      WRITE ( output_unit ) la1 , lo1 , dx , dy , lov , latin1 , latin2
    ELSE IF ( llflag .EQ. 5 ) THEN
      WRITE ( output_unit ) la1 , lo1 , dx , dy , lov , latin1
    END IF
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
    WRITE ( output_unit ) 3
    WRITE ( output_unit ) hdate , 0. ,  field , units , desc , &
                          p_pa(k) , x , y , llflag
    IF      ( llflag .EQ. 1 ) THEN
      WRITE ( output_unit ) la1 , lo1 , dx , dy , latin1
    ELSE IF ( llflag .EQ. 3 ) THEN
      WRITE ( output_unit ) la1 , lo1 , dx , dy , lov , latin1 , latin2
    ELSE IF ( llflag .EQ. 5 ) THEN
      WRITE ( output_unit ) la1 , lo1 , dx , dy , lov , latin1
    END IF
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
    WRITE ( output_unit ) 3
    WRITE ( output_unit ) hdate , 0. ,  field , units , desc , &
                          p_pa(k) , x , y , llflag
    IF      ( llflag .EQ. 1 ) THEN
      WRITE ( output_unit ) la1 , lo1 , dx , dy , latin1
    ELSE IF ( llflag .EQ. 3 ) THEN
      WRITE ( output_unit ) la1 , lo1 , dx , dy , lov , latin1 , latin2
    ELSE IF ( llflag .EQ. 5 ) THEN
      WRITE ( output_unit ) la1 , lo1 , dx , dy , lov , latin1
    END IF
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
    WRITE ( output_unit ) 3
    WRITE ( output_unit ) hdate , 0. ,  field , units , desc , &
                          p_pa(k) , x , y , llflag
    IF      ( llflag .EQ. 1 ) THEN
      WRITE ( output_unit ) la1 , lo1 , dx , dy , latin1
    ELSE IF ( llflag .EQ. 3 ) THEN
      WRITE ( output_unit ) la1 , lo1 , dx , dy , lov , latin1 , latin2
    ELSE IF ( llflag .EQ. 5 ) THEN
      WRITE ( output_unit ) la1 , lo1 , dx , dy , lov , latin1
    END IF
    d2d = ht(:,:,k)
    WRITE ( output_unit ) d2d
    PRINT '(A,F9.1,A,F8.1,A,F8.1)', 'Level (Pa):', p_pa(k), ' Min: ', MINVAL(d2d),&
            ' Max: ', MAXVAL(d2d)
  ENDDO var_ht

  ! Sea-level Pressure field
  field = 'PMSL     '
  units = 'Pa                       '
  desc  = 'Sea-level pressure                            '
  PRINT *, 'FIELD = ', field
  PRINT *, 'UNITS = ', units
  PRINT *, 'DESC =  ',desc
  WRITE ( output_unit ) 3
  WRITE ( output_unit ) hdate , 0. ,  field , units , desc , &
          201300., x , y , llflag
  IF      ( llflag .EQ. 1 ) THEN
    WRITE ( output_unit ) la1 , lo1 , dx , dy , latin1
  ELSE IF ( llflag .EQ. 3 ) THEN
    WRITE ( output_unit ) la1 , lo1 , dx , dy , lov , latin1 , latin2
  ELSE IF ( llflag .EQ. 5 ) THEN 
    WRITE ( output_unit ) la1 , lo1 , dx , dy , lov , latin1
  END IF
  WRITE ( output_unit ) slp
  PRINT '(A,F9.1,A,F9.1,A,F9.1)', 'Level (Pa):', p_pa(k), ' Min: ', MINVAL(slp),&
            ' Max: ', MAXVAL(slp)

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
      WRITE ( output_unit ) 3
      WRITE ( output_unit ) hdate , 0. ,  field , units , desc , p_pa(k) , x , y , llflag
      IF      ( llflag .EQ. 1 ) THEN
        WRITE ( output_unit ) la1 , lo1 , dx , dy , latin1
      ELSE IF ( llflag .EQ. 3 ) THEN
        WRITE ( output_unit ) la1 , lo1 , dx , dy , lov , latin1 , latin2
      ELSE IF ( llflag .EQ. 5 ) THEN
        WRITE ( output_unit ) la1 , lo1 , dx , dy , lov , latin1
      END IF
      d2d = lwc(:,:,k)
      WRITE ( output ) d2d
      PRINT '(A,F9.1,A,F8.6,A,F8.6)', 'Level (Pa):', p_pa(k), ' Min: ', MINVAL(d2d),&
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
      WRITE ( output_unit ) 3
      WRITE ( output_unit ) hdate , 0. ,  field , units , desc , p_pa(k) , x , y , llflag
      IF      ( llflag .EQ. 1 ) THEN
        WRITE ( output_unit ) la1 , lo1 , dx , dy , latin1
      ELSE IF ( llflag .EQ. 3 ) THEN
        WRITE ( output_unit ) la1 , lo1 , dx , dy , lov , latin1 , latin2
      ELSE IF ( llflag .EQ. 5 ) THEN
        WRITE ( output_unit ) la1 , lo1 , dx , dy , lov , latin1
      END IF
      d2d = ice(:,:,k)
      WRITE ( output ) d2d
      PRINT '(A,F9.1,A,F8.6,A,F8.6)', 'Level (Pa):', p_pa(k), ' Min: ', MINVAL(d2d),&
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
      WRITE ( output_unit ) 3
      WRITE ( output_unit ) hdate , 0. ,  field , units , desc , p_pa(k) , x , y , llflag
      IF      ( llflag .EQ. 1 ) THEN
        WRITE ( output_unit ) la1 , lo1 , dx , dy , latin1
      ELSE IF ( llflag .EQ. 3 ) THEN
        WRITE ( output_unit ) la1 , lo1 , dx , dy , lov , latin1 , latin2
      ELSE IF ( llflag .EQ. 5 ) THEN
        WRITE ( output_unit ) la1 , lo1 , dx , dy , lov , latin1
      END IF
      d2d = rai(:,:,k)
      WRITE ( output_unit ) d2d
      PRINT '(A,F9.1,A,F8.6,A,F8.6)', 'Level (Pa):', p_pa(k), ' Min: ', MINVAL(d2d),&
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
      WRITE ( output_unit ) 3
      WRITE ( output_unit ) hdate , 0. ,  field , units , desc , p_pa(k) , x , y , llflag
      IF      ( llflag .EQ. 1 ) THEN
        WRITE ( output_unit ) la1 , lo1 , dx , dy , latin1
      ELSE IF ( llflag .EQ. 3 ) THEN
        WRITE ( output_unit ) la1 , lo1 , dx , dy , lov , latin1 , latin2
      ELSE IF ( llflag .EQ. 5 ) THEN
        WRITE ( output_unit ) la1 , lo1 , dx , dy , lov , latin1
      END IF
      d2d = sno(:,:,k)
      WRITE ( output_unit ) d2d
      PRINT '(A,F9.1,A,F8.6,A,F8.6)', 'Level (Pa):', p_pa(k), ' Min: ', MINVAL(d2d),&
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
      WRITE ( output_unit ) 3
      WRITE ( output_unit ) hdate , 0. ,  field , units , desc , p_pa(k) , x , y , llflag
      IF      ( llflag .EQ. 1 ) THEN
        WRITE ( output_unit ) la1 , lo1 , dx , dy , latin1
      ELSE IF ( llflag .EQ. 3 ) THEN
        WRITE ( output_unit ) la1 , lo1 , dx , dy , lov , latin1 , latin2
      ELSE IF ( llflag .EQ. 5 ) THEN
        WRITE ( output_unit ) la1 , lo1 , dx , dy , lov , latin1
      END IF
      d2d = pic(:,:,k)
      WRITE ( output_unit ) d2d
      PRINT '(A,F9.1,A,F8.6,A,F8.6)', 'Level (Pa):', p_pa(k), ' Min: ', MINVAL(d2d),&
            ' Max: ', MAXVAL(d2d)
    END DO var_pic

  ENDIF

  CLOSE (output_unit)
  DEALLOCATE (d2d)
  DEALLOCATE (p_pa)
  RETURN
  END SUBROUTINE output_pregrid_v3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE lapsprep_output
