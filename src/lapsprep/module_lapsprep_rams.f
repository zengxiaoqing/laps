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

  SUBROUTINE output_ralph2_format(p, u, v, t, ht, rh, slp, psfc, snocov,tskin)

  ! Subroutine to output data in RAMS 4.x "RALPH 2" format. The skin temp from
  ! LAPS (lsx/tgd) is used for SST.
 
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
  REAL, INTENT(IN)               :: snocov(:,:)! Snow cover (fract)
  REAL, INTENT(IN)               :: tskin(:,:) ! Skin temp
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
  REAL, ALLOCATABLE              :: ut ( : , : )
  REAL, ALLOCATABLE              :: vt ( : , : )
  CHARACTER (LEN=256)            :: output_file_name
  CHARACTER (LEN=15)             :: date_string

  ! Build the output file name
 
  output_prefix = TRIM(laps_data_root)// '/lapsprd/lapsprep/rams/LAPS'
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

  WRITE(output_unit, '(I4.4,2I4,I6,4I4)') valid_yyyy, valid_mm, &
                                      valid_dd,valid_hh*100+valid_min, &
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
    IF (ABS(latin2).NE.90.) THEN
       PRINT '(A)', 'This is a local stereographic, so I quit.'
       stop 'unsupported projection'
    ELSE
      proj_flag = 3
      latin1_out = latin1
      latin2_out = latin2
    ENDIF
  ELSE
    PRINT '(A,A,A)','RALPH2 unsupported map projection: ', &
        TRIM(grid_type),'.  I quit.'
    STOP 'unsupported projection'
  END IF
  WRITE(output_unit,'(I1,2F10.1,7F10.3)') proj_flag, dx*1000, dy*1000, &
                                          la1, lo1, la2, lo2, &
                                          latin1_out, lov, &
                                          latin2_out

  ! Record 4 - Vertical level specification
  ALLOCATE (p_int (z3+1) )
  p_int = NINT(p)
  WRITE(output_unit, *) level_flag, p_int(z3:1:-1)
  DEALLOCATE(p_int)

  ALLOCATE (ut(x,y)) ! Array for true u-winds
  ALLOCATE (vt(x,y)) ! Array for true v-winds
  ! Now write out the pressure level data, looping by level

  level_loop: DO k = z3,1,-1

    ! Rotate the u and v winds to true if proj_flag EQ 2 
    IF (proj_flag .EQ. 2) THEN
      rotate_winds_j: DO j=1,y
        rotate_winds_i: DO i=1,x
          CALL uvgrid_to_uvtrue(u(i,j,k),v(i,j,k), &
                                ut(i,j),vt(i,j),&
                                lons(i,j))
        ENDDO rotate_winds_i
      ENDDO rotate_winds_j
    ELSE
      ut(:,:) = u(:,:,k)
      vt(:,:) = v(:,:,k)
    ENDIF
     
    WRITE(output_unit,'(8F10.3)') (( ut(i,j),i=1,x),j=1,y)
    WRITE(output_unit,'(8F10.3)') (( vt(i,j),i=1,x),j=1,y)
    WRITE(output_unit,'(8F10.3)') (( t(i,j,k),i=1,x),j=1,y)
    WRITE(output_unit,'(8F10.2)') (( ht(i,j,k),i=1,x),j=1,y)
    WRITE(output_unit,'(8F10.5)') (( rh(i,j,k)*.01,i=1,x),j=1,y)
 
  ENDDO level_loop
  DEALLOCATE (ut)
  DEALLOCATE (vt)

  ! Write out the surface fields
   
   WRITE(output_unit,'(8F10.3)') (( slp(i,j)*.01,i=1,x),j=1,y)
   WRITE(output_unit,'(8F10.3)') (( psfc(i,j)*.01,i=1,x),j=1,y)
   WRITE(output_unit,'(8F10.3)') (( t(i,j,z3+1),i=1,x),j=1,y)
   WRITE(output_unit,'(8F10.4)') (( snocov(i,j),i=1,x),j=1,y)
   WRITE(output_unit,'(8F10.4)') (( tskin(i,j),i=1,x),j=1,y)
  ! Close the file
  
  CLOSE(output_unit)
  RETURN
  END SUBROUTINE output_ralph2_format
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    
END MODULE lapsprep_rams
