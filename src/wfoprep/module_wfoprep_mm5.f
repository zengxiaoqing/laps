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



MODULE wfoprep_mm5

! PURPOSE
! =======
! Module to contain the various output routines needed for lapsprep
! to support initializition of the MM5v3 model.
!
! SUBROUTINES CONTAINED
! =====================
! output_mm5v3_basic     - Outputs state variables on pressure surfaces
!                          plus MSLP and topography
! output_mm5v3_sfc       - Outputs the surface fields
! write_pregrid_header   - Writes pregrid headers
!
! REMARKS
! =======
! 
!
! HISTORY
! =======
! 1 Oct 2001 -- Original -- Brent Shaw

  USE map_utils
  IMPLICIT NONE

  INTEGER, PARAMETER :: pg_version = 3
  INTEGER, PARAMETER :: output_unit = 78
  PUBLIC output_mm5v3_basic, output_mm5v3_sfc
  INTEGER, PARAMETER, PUBLIC   :: MM5MODE_NEW = 1
  INTEGER, PARAMETER, PUBLIC   :: MM5MODE_APPEND = 2
  LOGICAL, PARAMETER           :: verbose = .true.
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE output_mm5v3_basic(i4time_cycle, i4time_valid, proj, &
          np_ht, np_t, np_u, np_v, np_rh, &
          ht_plevels, t_plevels, u_plevels, v_plevels, rh_plevels, &
          z, t, u, v, rh, mslp, topo, &
          ext_data_path, output_name, mode, istatus)


     IMPLICIT NONE
     INTEGER, INTENT(IN)         :: i4time_cycle
     INTEGER, INTENT(IN)         :: i4time_valid
     TYPE(proj_info),INTENT(IN)  :: proj
     INTEGER, INTENT(IN)         :: np_ht
     INTEGER, INTENT(IN)         :: np_t
     INTEGER, INTENT(IN)         :: np_u
     INTEGER, INTENT(IN)         :: np_v
     INTEGER, INTENT(IN)         :: np_rh
     REAL, INTENT(IN)            :: ht_plevels(np_ht)
     REAL, INTENT(IN)            :: t_plevels(np_t)
     REAL, INTENT(IN)            :: u_plevels(np_u)
     REAL, INTENT(IN)            :: v_plevels(np_v)
     REAL, INTENT(IN)            :: rh_plevels(np_rh)
     REAL, INTENT(IN)            :: z ( : , : , : )
     REAL, INTENT(IN)            :: t ( : , : , : )
     REAL, INTENT(IN)            :: u ( : , : , : )
     REAL, INTENT(IN)            :: v ( : , : , : )
     REAL, INTENT(IN)            :: rh ( : , : , : )
     REAL, INTENT(IN)            :: mslp ( : , : )
     REAL, INTENT(IN)            :: topo ( : , : )
     CHARACTER(LEN=256),INTENT(IN):: ext_data_path
     CHARACTER(LEN=32 ),INTENT(IN):: output_name
     INTEGER, INTENT(IN)         :: mode
     INTEGER, INTENT(OUT)        :: istatus
   
     INTEGER                     :: k
     CHARACTER(LEN=256)          :: outfile
     CHARACTER(LEN=24)           :: atime
     CHARACTER(LEN=3)            :: amonth
     CHARACTER(LEN=2)            :: amonth_num
     INTEGER                     :: m
     CHARACTER (LEN=24) :: hdate
     CHARACTER (LEN=9)  :: field
     CHARACTER (LEN=25) :: units
     CHARACTER (LEN=46) :: desc     
     REAL               :: xfcst 
     REAL, PARAMETER    :: slp_level = 201300.0
     REAL, PARAMETER    :: sfc_level = 200100.0  
     istatus = 1
  
     CALL make_pregrid_filename(ext_data_path,output_name, &
        i4time_valid,outfile) 

     ! Compute xfcst
     xfcst = FLOAT(i4time_valid - i4time_cycle)/3600.
 
     ! Compute hdate
     CALL make_hdate_from_i4time(i4time_valid,hdate)

     ! Open the file, using the mode dependency
     PRINT *, 'Opening file: ', TRIM(outfile)
     IF (mode .EQ. MM5MODE_NEW) THEN
       OPEN ( FILE   = TRIM(outfile)    , &
         UNIT   = output_unit        , &
         FORM   = 'UNFORMATTED' , &
         STATUS = 'REPLACE'     , &
         ACCESS = 'SEQUENTIAL'    )
     ELSE IF (mode .EQ. MM5MODE_APPEND) THEN
       OPEN ( FILE   = TRIM(outfile)    , &
         UNIT   = output_unit        , &
         FORM   = 'UNFORMATTED' , &
         STATUS = 'UNKNOWN'     , &
         ACCESS = 'SEQUENTIAL', &
         POSITION = 'APPEND'    ) 
     ELSE 
       PRINT *, 'Uknown open mode for MM5v3 Output: ',mode
       istatus = 0
       RETURN
     ENDIF

     IF (verbose) THEN
       PRINT *, 'PREGRID VERSION =', pg_version
       PRINT *, 'HDATE = ',hdate
       PRINT *, 'XFCST = ', xfcst 
       PRINT *, 'NX = ', proj%nx
       PRINT *, 'NY = ', proj%ny
       PRINT *, 'IPROJ = ', proj%code
       PRINT *, 'STARTLAT = ',proj%lat1
       PRINT *, 'STARTLON = ',proj%lon1
       PRINT *, 'DX = ', proj%dx*0.001
       PRINT *, 'DY = ', proj%dx*0.001
       PRINT *, 'XLONC = ', proj%stdlon
       PRINT *, 'TRUELAT1 = ', proj%truelat1
       PRINT *, 'TRUELAT2 = ', proj%truelat2
     ENDIF

     ! Output temperature
     field = 'T        '
     units = 'K                        '
     desc  = 'Temperature                                   '      
     var_t : DO k = 1 , np_t
       CALL write_pregrid_header(proj,field,units,desc,t_plevels(k),hdate,xfcst)
       WRITE ( output_unit ) t(:,:,k)
       IF (verbose) THEN
         PRINT '(A,1x,A,A,F9.1,A,F5.1,A,F5.1)',field,units, &
               'Level (Pa):',t_plevels(k),' Min: ', &
                MINVAL(t(:,:,k)),' Max: ', MAXVAL(t(:,:,k))
       ENDIF
     ENDDO var_t

     ! Do u-component of wind
     field = 'U        '
     units = 'm s{-1}                  '
     desc = 'u-component of velocity, rotated to grid      '
     var_u : DO k = 1 , np_u
       CALL write_pregrid_header(proj,field,units,desc,u_plevels(k),hdate,xfcst)
       WRITE ( output_unit ) u(:,:,k)
       IF (verbose) THEN
         PRINT '(A,1x,A,A,F9.1,A,F5.1,A,F5.1)', field,units, &
                'Level (Pa):', u_plevels(k),  &
                ' Min: ', MINVAL(u(:,:,k)),&
                ' Max: ', MAXVAL(u(:,:,k))
       ENDIF
     ENDDO var_u

     ! Do v-component of wind
     field = 'V        '
     units = 'm s{-1}                  '
     desc = 'v-component of velocity, rotated to grid      '
     var_v : DO k = 1 , np_v
       CALL write_pregrid_header(proj,field,units,desc,v_plevels(k),hdate,xfcst)
       WRITE ( output_unit ) v(:,:,k)
       If (verbose) THEN
         PRINT '(A,1x,A,A,F9.1,A,F5.1,A,F5.1)',field,units, &
            'Level (Pa):', v_plevels(k), &
            ' Min: ', MINVAL(v(:,:,k)),&
            ' Max: ', MAXVAL(v(:,:,k))
       ENDIF
     ENDDO var_v

     ! Relative Humidity
     field = 'RH       '
     units = '%                        '
     desc  = 'Relative humidity                             '
     var_rh : DO k = 1 , np_rh
       CALL write_pregrid_header(proj,field,units,desc,rh_plevels(k),hdate,xfcst)
       WRITE ( output_unit ) rh(:,:,k)
       IF (verbose) THEN
         PRINT '(A,1x,A,A,F9.1,A,F5.1,A,F5.1)',field,units, &
            'Level (Pa):', rh_plevels(k),&
            ' Min: ', MINVAL(rh(:,:,k)),&
            ' Max: ', MAXVAL(rh(:,:,k))
       ENDIF
     ENDDO var_rh

     ! Do the heights
     field = 'HGT      '
     units = 'm                        '
     desc  = 'Geopotential height                           '
     var_ht : DO k = 1 , np_ht
       CALL write_pregrid_header(proj,field,units,desc,ht_plevels(k),hdate,xfcst)
       WRITE ( output_unit ) z(:,:,k)
       IF (verbose) THEN
         PRINT '(A,x,A,A,F9.1,A,F8.1,A,F8.1)', field,units, &
            'Level (Pa):', ht_plevels(k), &
            ' Min: ', MINVAL(z(:,:,k)),&
            ' Max: ', MAXVAL(z(:,:,k))
       ENDIF
     ENDDO var_ht

     ! Terrain height
     field = 'SOILHGT  '
     units = 'm                        '
     desc  = 'Height of topography                          '
     CALL write_pregrid_header(proj,field,units,desc,sfc_level,hdate,xfcst)
     WRITE ( output_unit ) topo
     IF (verbose) THEN
       PRINT '(a,x,a,A,F9.1,A,F9.1,A,F9.1)', field,units, & 
            'Level (Pa):',sfc_level, &
            ' Min: ', MINVAL(topo),&
            ' Max: ', MAXVAL(topo)
     ENDIF
     ! Sea-level Pressure field
     field = 'PMSL     '
     units = 'Pa                       '
     desc  = 'Sea-level pressure                            '
     CALL write_pregrid_header(proj,field,units,desc,slp_level,hdate,xfcst)
     WRITE ( output_unit ) mslp
     IF (verbose) THEN
       PRINT '(A,x,A,A,F9.1,A,F9.1,A,F9.1)',field,units, &
            'Level (Pa):', slp_level, &
            ' Min: ', MINVAL(mslp),&
            ' Max: ', MAXVAL(mslp)
     ENDIF
     CLOSE(output_unit)
    RETURN
  END SUBROUTINE output_mm5v3_basic
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE output_mm5v3_sfc(i4time_cycle,i4time_valid,proj, &
                              tsf,usf,vsf,rhsf, &
                              ext_data_path,output_name, mode, istatus)

    IMPLICIT NONE

    INTEGER,INTENT(IN)                     :: i4time_cycle
    INTEGER,INTENT(IN)                     :: i4time_valid
    TYPE(proj_info),INTENT(IN)             :: proj
    REAL,INTENT(IN)                        :: tsf(:,:)
    REAL,INTENT(IN)                        :: usf(:,:)
    REAL,INTENT(IN)                        :: vsf(:,:)
    REAL,INTENT(IN)                        :: rhsf(:,:)
    CHARACTER(LEN=256),INTENT(IN)          :: ext_data_path
    CHARACTER(LEN=32),INTENT(IN)           :: output_name
    INTEGER,INTENT(IN)                     :: mode
    INTEGER,INTENT(OUT)                    :: istatus

    CHARACTER(LEN=256)                     :: outfile

    CHARACTER (LEN=24) :: hdate
    CHARACTER (LEN=9)  :: field
    CHARACTER (LEN=25) :: units
    CHARACTER (LEN=46) :: desc  
    REAL               :: xfcst   
    REAL, PARAMETER    :: slp_level = 201300.0
    REAL, PARAMETER    :: sfc_level = 200100.0  
    istatus = 1

    CALL make_pregrid_filename(ext_data_path,output_name, &
                               i4time_valid,outfile)       
    ! Compute xfcst
    xfcst = FLOAT(i4time_valid - i4time_cycle)/3600.

    ! Compute hdate
    CALL make_hdate_from_i4time(i4time_valid,hdate)  
   
    ! Open the file, using the mode dependency
    PRINT *, 'Opening file: ', TRIM(outfile)
    IF (mode .EQ. MM5MODE_NEW) THEN
      OPEN ( FILE   = TRIM(outfile)    , &
             UNIT   = output_unit        , &
             FORM   = 'UNFORMATTED' , &
             STATUS = 'REPLACE'     , &
             ACCESS = 'SEQUENTIAL'    )
    ELSE IF (mode .EQ. MM5MODE_APPEND) THEN
      OPEN ( FILE   = TRIM(outfile)    , &
             UNIT   = output_unit        , &
             FORM   = 'UNFORMATTED' , &
             STATUS = 'UNKNOWN'     , &
             ACCESS = 'SEQUENTIAL', &
             POSITION = 'APPEND'    )
    ELSE
      PRINT *, 'Uknown open mode for MM5v3 Output: ',mode
      istatus = 0
      RETURN
    ENDIF                  

    ! Output temperature
    field = 'T        '
    units = 'K                        '
    desc  = 'Temperature                                   '
    CALL write_pregrid_header(proj,field,units,desc,sfc_level,hdate,xfcst)
    WRITE ( output_unit ) tsf
    IF (verbose) THEN
      PRINT '(A,x,A,A,F9.1,A,F5.1,A,F5.1)',field,units, &
           'Level (Pa):',sfc_level,' Min: ', &
           MINVAL(tsf),' Max: ', MAXVAL(tsf)
    ENDIF

    ! Output u-wind component
    field = 'U        '
    units = 'm s{-1}                  '
    desc  = 'u-component of wind, rotated to grid          '
    CALL write_pregrid_header(proj,field,units,desc,sfc_level,hdate,xfcst)
    WRITE ( output_unit ) usf
    IF (verbose) THEN
      PRINT '(A,x,A,A,F9.1,A,F5.1,A,F5.1)',field,units, &
          'Level (Pa):',sfc_level,' Min: ', &
           MINVAL(usf),' Max: ', MAXVAL(usf)   
    ENDIF

    ! Output v-wind component
    field = 'V        '
    units = 'm s{-1}                  '
    desc  = 'v-component of wind, rotated to grid          '
    CALL write_pregrid_header(proj,field,units,desc,sfc_level,hdate,xfcst)
    WRITE ( output_unit ) vsf
    IF (verbose) THEN
      PRINT '(A,X,A,A,F9.1,A,F5.1,A,F5.1)',field,units, &
            'Level (Pa):',sfc_level,' Min: ', &
           MINVAL(vsf),' Max: ', MAXVAL(vsf) 
    ENDIF

    ! Output Relative Humidity
    field = 'RH       '
    units = '%                        '
    desc  = 'Relative Humidity                             '
    CALL write_pregrid_header(proj,field,units,desc,sfc_level,hdate,xfcst)
    WRITE ( output_unit ) rhsf
    IF (verbose) THEN
      PRINT '(A,X,A,A,F9.1,A,F5.1,A,F5.1)', field,units, &
           'Level (Pa):',sfc_level,' Min: ', &
           MINVAL(rhsf),' Max: ', MAXVAL(rhsf)            
    ENDIF
               
    CLOSE (output_unit) 
    RETURN
  END SUBROUTINE output_mm5v3_sfc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE output_mm5v3_fsst(i4time_cycle,i4time_valid,proj, &
                              fsst,ext_data_path,output_name, mode, istatus)

    IMPLICIT NONE

    INTEGER,INTENT(IN)                     :: i4time_cycle
    INTEGER,INTENT(IN)                     :: i4time_valid
    TYPE(proj_info),INTENT(IN)             :: proj
    REAL,INTENT(IN)                        :: fsst(:,:)
    CHARACTER(LEN=256),INTENT(IN)          :: ext_data_path
    CHARACTER(LEN=32),INTENT(IN)           :: output_name
    INTEGER,INTENT(IN)                     :: mode
    INTEGER,INTENT(OUT)                    :: istatus

    CHARACTER(LEN=256)                     :: outfile

    CHARACTER (LEN=24) :: hdate
    CHARACTER (LEN=9)  :: field
    CHARACTER (LEN=25) :: units
    CHARACTER (LEN=46) :: desc  
    REAL               :: xfcst   
    REAL, PARAMETER    :: sfc_level = 200100.0
    istatus = 1

    outfile = TRIM(ext_data_path) // TRIM(output_name) // &
               ':FSST' 

    ! Compute xfcst
    xfcst = FLOAT(i4time_valid - i4time_cycle)/3600.

    ! Compute hdate
    CALL make_hdate_from_i4time(i4time_valid,hdate)
    ! Open the file, using the mode dependency

    PRINT *, 'Opening file: ', TRIM(outfile)
    IF (mode .EQ. MM5MODE_NEW) THEN
      OPEN ( FILE   = TRIM(outfile)    , &
             UNIT   = output_unit        , &
             FORM   = 'UNFORMATTED' , &
             STATUS = 'REPLACE'     , &
             ACCESS = 'SEQUENTIAL'    )
    ELSE IF (mode .EQ. MM5MODE_APPEND) THEN
      OPEN ( FILE   = TRIM(outfile)    , &
             UNIT   = output_unit        , &
             FORM   = 'UNFORMATTED' , &
             STATUS = 'UNKNOWN'     , &
             ACCESS = 'SEQUENTIAL', &
             POSITION = 'APPEND'    )
    ELSE
      PRINT *, 'Uknown open mode for MM5v3 Output: ',mode
      istatus = 0
      RETURN
    ENDIF

    ! Output temperature
    field = 'SST      '
    units = 'K                        '
    desc  = 'Sea-surface Temperature (faked)               '
    CALL write_pregrid_header(proj,field,units,desc,sfc_level,hdate,xfcst)
    WRITE ( output_unit ) fsst
    IF (verbose) THEN
      PRINT '(A,x,A,A,F9.1,A,F5.1,A,F5.1)',field,units, &
           'Level (Pa):',sfc_level,' Min: ', &
           MINVAL(fsst),' Max: ', MAXVAL(fsst)
    ENDIF
    CLOSE (output_unit)
    RETURN
  END SUBROUTINE output_mm5v3_fsst
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE make_pregrid_filename(ext_data_path,output_name,i4time,filename)

    IMPLICIT NONE
    CHARACTER(LEN=256),INTENT(IN)           :: ext_data_path
    CHARACTER(LEN=32), INTENT(IN)           :: output_name
    INTEGER, INTENT(IN)                     :: i4time
    CHARACTER(LEN=256),INTENT(OUT)          :: filename
    CHARACTER(LEN=24)                       :: hdate
    
    CALL make_hdate_from_i4time(i4time,hdate)
 
    ! Build the output file name

    filename = TRIM(ext_data_path) // TRIM(output_name) // &
               ':' // hdate(1:13)

    RETURN

  END SUBROUTINE make_pregrid_filename
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE make_hdate_from_i4time(i4time,hdate)

    IMPLICIT NONE
    INTEGER, INTENT(IN)                    :: i4time
    CHARACTER(LEN=24),INTENT(OUT)          :: hdate

    CHARACTER(LEN=3)                       :: amonth(12)
    CHARACTER(LEN=24)                      :: atime
    INTEGER                                :: istatus
    INTEGER                                :: m
    INTEGER                                :: dom
    CHARACTER(LEN=2)                       :: amonth_num

    DATA amonth/'JAN','FEB','MAR','APR','MAY','JUN',  &
                 'JUL','AUG','SEP','OCT','NOV','DEC'/   

    ! Build the file name by converting the i4time_valid
    ! to YYYY-MM-DD_HH format
   
    CALL cv_i4tim_asc_lp(i4time,atime,istatus)
    IF (istatus .NE. 1) THEN
      PRINT *, 'Problem converting i4time!'
      RETURN
    ENDIF
    findmonth: DO m = 1, 12
      IF (atime(4:6).EQ.amonth(m)) EXIT findmonth
    ENDDO findmonth

    ! Ensure we make day into a 2-digit value
    READ(atime(1:2),'(I2)') dom
    WRITE(atime(1:2),'(I2.2)') dom
    WRITE(amonth_num, '(I2.2)') m
    hdate  = atime(8:11) // '-' // &
             amonth_num // '-' // &
             atime(1:2) // '_' // &
             atime(13:23) // '00'  
    RETURN 
  END SUBROUTINE make_hdate_from_i4time                 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE write_pregrid_header(proj,field,units,desc,level,hdate,xfcst)
 
  ! Writes the gribprep header given the filed, units, description, and level

  IMPLICIT NONE
  TYPE(proj_info),INTENT(IN)    :: proj
  CHARACTER(LEN=9), INTENT(IN)  :: field
  CHARACTER(LEN=25),INTENT(IN)  :: units
  CHARACTER(LEN=46),INTENT(IN)  :: desc
  REAL, INTENT(IN)              :: level
  CHARACTER(LEN=24),INTENT(IN)  :: hdate
  REAL,INTENT(IN)               :: xfcst
 
  WRITE ( output_unit ) pg_version
  WRITE ( output_unit ) hdate,xfcst,field,units,desc,level,proj%nx, &
                        proj%ny,proj%code
  SELECT CASE (proj%code)
    CASE(0)
      WRITE ( output_unit) proj%lat1, proj%lon1,proj%dlat,proj%dlon 
    CASE(1)
      WRITE ( output_unit ) proj%lat1,proj%lon1,proj%dx*0.001,proj%dx*0.001, &
                            proj%truelat1
    CASE(3)
      WRITE ( output_unit ) proj%lat1,proj%lon1,proj%dx*0.001,proj%dx*0.001, &
                            proj%stdlon,proj%truelat1,proj%truelat2
    CASE(5)
      WRITE ( output_unit ) proj%lat1,proj%lon1,proj%dx*0.001,proj%dx*0.001, &
                            proj%stdlon,proj%truelat1
  END SELECT

  END SUBROUTINE write_pregrid_header
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE wfoprep_mm5
  
