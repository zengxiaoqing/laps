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
MODULE lapsprep_rams

! PURPOSE
! =======
! Module to contain the various output routines needed for lapsprep.
!
! SUBROUTINES CONTAINED
! =====================
! output_ralph2_format  - Used to support RAMS 4.x initializations
!
! REMARKS
! =======

!
! HISTORY
! =======
! 28 Nov 2000 -- Original -- Brent Shaw

  USE setup
  USE laps_static
  USE date_pack

  PRIVATE
  PUBLIC output_ralph2_format
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE output_ralph2_format(p, u, v, t, ht, rh, slp, psfc, snodep)

  ! Subroutine to output data in RAMS 4.x "RALPH 2" format. We do not have 
  ! SST data in LAPS yet, so that is not included for now.  
 
  ! Note that the P array has a bogus 2001 mb value as the last entry
  ! to designate the surface for MM5 applications.  The surface values
  ! of the state variables are contained in the last layer of their
  ! respective 3D arrays. 

  ! The u/v winds are grid-relative, and for now Ralph 2 expects 
  ! true winds, so this routine must rotate them to true.

  ! Surface and MSL Pressure must be converted to mb for RAMS.

  ! P levels must be written as integer values

  ! RH must be written as a fraction 

  ! 3d Variables are stored top down, but RAMS wants them 
  ! bottom up.

  ! Note that RALPH only supports Lat/Lon or Lambert projections.
  ! However, Lambert can be used to fake a polar stereographic by
  ! setting both standard lats equal to the pole. However, if the LAPS
  ! grid is a mercator, then you are out of luck.

  IMPLICIT NONE
  REAL, INTENT(IN)               :: p(:)      !Pressure levels in mb
  REAL, INTENT(IN)               :: u(:,:,:)  !U-component of wind wrt grid in m/s
  REAL, INTENT(IN)               :: v(:,:,:)  !V-component of wind wrf grid in m/s 
  REAL, INTENT(IN)               :: t(:,:,:)  !Temperature in K
  REAL, INTENT(IN)               :: ht(:,:,:) !Geopotential Height in m
  REAL, INTENT(IN)               :: rh(:,:,:) !Relative Humidity in %
  REAL, INTENT(IN)               :: slp(:,:)  !MSL Pressure in Pa
  REAL, INTENT(IN)               :: psfc(:,:) !Surface Pressure in Pa
  REAL, INTENT(IN)               :: snodep(:,:)! Snow depth in m

  ! Local Variables

  INTEGER, PARAMETER             :: marker = 999999
  INTEGER, PARAMETER             :: version= 2
  INTEGER, PARAMETER             :: vt_increment = 0
  INTEGER, PARAMETER             :: latlon_proj = 1
  INTEGER, PARAMETER             :: lambert_proj = 2
  INTEGER, PARAMETER             :: level_flag = 1
  INTEGER, PARAMETER             :: output_unit = 78
  INTEGER                        :: proj_flag
  INTEGER                        :: yyyyddd, valid_mm, valid_dd
  INTEGER                        :: i,j,k
  INTEGER,ALLOCATABLE            :: p_int(:)
  REAL                           :: latin1_out
  REAL                           :: latin2_out
  CHARACTER (LEN=256)            :: output_file_name
  CHARACTER (LEN=15)             :: date_string

  ! Build the output file name
 
  output_prefix = TRIM(laps_data_root)// '/lapsprd/dprep/' // &
                  TRIM(output_format(1:3)) // '_init'
  yyyyddd = valid_yyyy*1000 + valid_jjj
  CALL wrf_date_to_ymd(yyyyddd, valid_yyyy, valid_mm, valid_dd) 
  WRITE(date_string,'(I4.4,"-",I2.2,"-",I2.2,"-",I2.2,I2.2)') &
          valid_yyyy, valid_mm, valid_dd, valid_hh, valid_min
  output_file_name = TRIM(output_prefix) // ':' // date_string

  !  Open the file for sequential, unformatted output
  OPEN ( FILE   = TRIM(output_file_name)    , &
         UNIT   = output_unit        , &
         FORM   = 'FORMATTED' , &
         STATUS = 'UNKNOWN'     , &
         ACCESS = 'SEQUENTIAL'    )
  REWIND(output_unit)
  ! Write the header record as described in the Ralph 2 documentation

  ! Record 1 - Header
  WRITE(output_unit, '(I6.6,I4)') marker, version

  ! Record 2 - Time and dimension information

  WRITE(output_unit, '(I4.4,2I4,I6,I4,2I4)') valid_yyyy, valid_mm, &
                                               valid_hh*100+valid_min, &
                                               vt_increment, z3, x, y

  ! Record 3 - Projection info
  IF      ( grid_type(1:8)  .EQ. 'latlon'                 ) THEN
    proj_flag = 1
    latin1_out = latin1
    latin2_out = latin2
  ELSE IF ( grid_type(1:24) .EQ. 'secant lambert conformal' ) THEN
    proj_flag = 2
    latin1_out = latin1
    latin2_out = latin2
  ELSE IF ( grid_type(1:19) .EQ. 'polar stereographic'      ) THEN
    PRINT '(A)', 'RALPH2 does not support polar stereographic.'
    PRINT '(A)', ' Attempting to fake it with lambert conformal.'
    ! This fake only works if this is a polar stereographic, so
    ! we need to make sure this is actually not a local 
    ! stereographic.  
    IF (ABS(latin2).NE.90.) THEN
       PRINT '(A)', 'This is a local stereographic, so I quit.'
       stop 'unsupported projection'
    ELSE
      If (latin2 .EQ. -90) THEN
        latin1_out = -89.99999
      ELSE
        latin1_out = 89.999999
      ENDIF
    ENDIF
  ELSE
    PRINT '(A,A,A)','RALPH2 unsupported map projection: ', &
        TRIM(grid_type),'.  I quit.'
    STOP 'unsupported projection'
  END IF
  WRITE(output_unit,'(I1,2F10.1,7F10.3)') proj_flag, dx, dy, &
                                          la1, lo1, la2, lo2, &
                                          latin1_out, lov, &
                                          latin2_out

  ! Record 4 - Vertical level specification
  ALLOCATE (p_int (z3+1) )
  p_int = NINT(p)
  WRITE(output_unit, *) level_flag, p_int(z3:1:-1)
  DEALLOCATE(p_int)
  ! Now write out the pressure level data, looping by level

  level_loop: DO k = z3,1,-1

    WRITE(output_unit,'(8F10.3)') (( u(i,j,k),i=1,x),j=1,y)
    WRITE(output_unit,'(8F10.3)') (( v(i,j,k),i=1,x),j=1,y)
    WRITE(output_unit,'(8F10.3)') (( t(i,j,k),i=1,x),j=1,y)
    WRITE(output_unit,'(8F10.2)') (( ht(i,j,k),i=1,x),j=1,y)
    WRITE(output_unit,'(8F10.5)') (( rh(i,j,k)*.01,i=1,x),j=1,y)
 
  ENDDO level_loop

  ! Write out the surface fields
   
   WRITE(output_unit,'(8F10.3)') (( slp(i,j)*.01,i=1,x),j=1,y)
   WRITE(output_unit,'(8F10.3)') (( psfc(i,j)*.01,i=1,x),j=1,y)
   WRITE(output_unit,'(8F10.3)') (( t(i,j,z3+1),i=1,x),j=1,y)
   WRITE(output_unit,'(8F10.4)') (( snodep(i,j),i=1,x),j=1,y)

  ! Close the file
  
  CLOSE(output_unit)
  RETURN
  END SUBROUTINE output_ralph2_format
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    
END MODULE lapsprep_rams
