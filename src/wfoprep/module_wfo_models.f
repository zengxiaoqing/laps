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

MODULE wfo_models

  ! Module that contains routines and variable pertaining to the setup
  ! of the background model grid being processed by sbnprep.

  USE map_utils
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'
  INTEGER                :: nfstatus
  INTEGER                :: maxlevels = 100
  INTEGER                :: maxtimes  = 100 
  
  PUBLIC open_wfofile, close_wfofile,get_wfomodel_var_levels, &
         get_wfomodel_fcsttimes, get_wfomodel_proj,get_wfomodel_var_inv,&
         read_wfomodel_data,get_wfomodel_topo
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE open_wfofile(nfname,nfid,istatus)

    ! Given a netCDF file name, this routine opens it and returns the
    ! integer netCDF file handle.

    IMPLICIT NONE
    CHARACTER(LEN=256), INTENT(IN)  :: nfname  ! Input file name
    INTEGER,INTENT(OUT)             :: nfid    ! Output unit number
    INTEGER,INTENT(OUT)             :: istatus ! Status flag (1=succes)
    
    istatus = 1
    nfstatus = NF_OPEN(nfname, NF_NOWRITE, nfid)
    IF (nfstatus .NE. NF_NOERR) THEN
      PRINT *, 'Problem with netCDF file:',TRIM(nfname)
      PRINT *, 'NetCDF Error = ', nfstatus
      istatus = 0
    ENDIF        
    RETURN 
  END SUBROUTINE open_wfofile  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE close_wfofile(nfid,istatus)

    ! Closes an open NetCDF file given the integer netCDF file handle.

    IMPLICIT NONE
    INTEGER, INTENT(IN)            :: nfid
    INTEGER, INTENT(OUT)           :: istatus
    
    istatus = 1
    nfstatus = NF_CLOSE(nfid)
    IF (nfstatus .NE. NF_NOERR) THEN
      PRINT *, 'Problem closing netCDF file.'
      PRINT *, 'NetCDF Error = ', nfstatus
      istatus = 0
    ENDIF
    RETURN
  END SUBROUTINE close_wfofile    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_wfomodel_proj(nfid,mname,proj,istatus)
   
    ! Populates a projection information data structure defined
    ! by module_map_utils.F for an already open WFO netCDF model
    ! file.  Use open_wfofile to open the file before calling this
    ! routine.  Getting projection info is a bit tricky, because the AWIPS
    ! model files have often been clipped and do not give the grid spacing
    ! at the actual true latitude of the projection.  Hence, the model name, 
    ! which corresponds to the last two elements of the directory the model
    ! file is located in (e.g., "CONUS212/MesoEta"), is used to make some
    ! educated guesses. 

    IMPLICIT NONE
    INTEGER, INTENT(IN)           :: nfid
    CHARACTER(LEN=132),INTENT(IN) :: mname
    TYPE(proj_info),INTENT(OUT)  :: proj  ! Declared via "USE map_utils" 
    INTEGER, INTENT(OUT)          :: istatus

    INTEGER                       :: nx,ny
    INTEGER                       :: vid,attnum,attid
    REAL                          :: truelat1,truelat2
    REAL                          :: stdlon
    REAL                          :: latsw,lonsw
    REAL                          :: latne,lonne
    REAL                          :: latne_c,lonne_c
    REAL                          :: difflat, difflon
    REAL                          :: dx,dy
    REAL                          :: truedx
    REAL                          :: latdxdy,londxdy
    CHARACTER(LEN=132)            :: gproj
    

    istatus = 1
    gproj(1:132) = ' ' 
    ! Get the horizontal dimensions
    nfstatus = NF_INQ_DIMID(nfid, 'x', vid)
    IF (nfstatus .NE. NF_NOERR) THEN
      PRINT *, 'Problem getting x dimension ID'
      istatus = 0
      RETURN
    ENDIF

    nfstatus = NF_INQ_DIMLEN(nfid, vid, nx)
    IF (nfstatus .NE. NF_NOERR) THEN
      PRINT *, 'Problem getting x dimension'
      istatus = 0
      RETURN
    ENDIF

    nfstatus = NF_INQ_DIMID(nfid, 'y', vid)
    IF (nfstatus .NE. NF_NOERR) THEN
      PRINT *, 'Problem getting y dimension ID'
      istatus = 0
      RETURN
    ENDIF

    nfstatus = NF_INQ_DIMLEN(nfid,vid,ny)
    IF (nfstatus .NE. NF_NOERR) THEN
      PRINT *, 'Problem getting y dimension'
      istatus = 0
      RETURN
    ENDIF                                           

    ! Get Projection Info
    nfstatus = NF_INQ_ATTID(nfid,attid,'projName',attnum)
    nfstatus = NF_GET_ATT_TEXT(nfid,attid,'projName',gproj)
    nfstatus = NF_INQ_ATTID(nfid,attid,'centralLat',attnum)
    nfstatus = NF_GET_ATT_REAL(nfid,attid,'centralLat',truelat1)
    nfstatus = NF_INQ_ATTID(nfid,attid,'centralLon',attnum)
    nfstatus = NF_GET_ATT_REAL(nfid,attid,'centralLon',stdlon)
    nfstatus = NF_INQ_ATTID(nfid,attid,'rotation',attnum)
    nfstatus = NF_GET_ATT_REAL(nfid,attid,'rotation',truelat2)
    nfstatus = NF_INQ_ATTID(nfid,attid,'lat00',attnum)
    nfstatus = NF_GET_ATT_REAL(nfid,attid,'lat00',latsw)
    nfstatus = NF_INQ_ATTID(nfid,attid,'lon00',attnum)
    nfstatus = NF_GET_ATT_REAL(nfid,attid,'lon00',lonsw)
    nfstatus = NF_INQ_ATTID(nfid,attid,'latNxNy',attnum)
    nfstatus = NF_GET_ATT_REAL(nfid,attid,'latNxNy',latne)
    nfstatus = NF_INQ_ATTID(nfid,attid,'lonNxNy',attnum)
    nfstatus = NF_GET_ATT_REAL(nfid,attid,'lonNxNy',lonne)
    nfstatus = NF_INQ_ATTID(nfid,attid,'dxKm',attnum)
    nfstatus = NF_GET_ATT_REAL(nfid,attid,'dxKm',dx)
    nfstatus = NF_INQ_ATTID(nfid,attid,'dyKm',attnum)
    nfstatus = NF_GET_ATT_REAL(nfid,attid,'dyKm',dy)
    nfstatus = NF_INQ_ATTID(nfid,attid,'latDxDy',attnum)
    nfstatus = NF_GET_ATT_REAL(nfid,attid,'latDxDy',latdxdy)
    nfstatus = NF_INQ_ATTID(nfid,attid,'lonDxDy',attnum)
    nfstatus = NF_GET_ATT_REAL(nfid,attid,'lonDxDy',londxdy)               

    ! Determine grid spacing at true latitude.  We have to do some
    ! level of hard coding for efficiency here, because the SBN data feed
    ! does not provide grid spacing at the true latitude.  Rather, it provides
    ! a "self-computed" grid spacing at the approximate center of the domain.  So,
    ! we will make a guess at the true dx, then use the mapping routines to
    ! verify that we get the correct coordinate conversion when using the "guessed"
    ! value. Down the road, we could put in some iterative method to do this using
    ! the provided dx/dy values as a starting point.


    IF (mname(1:8).EQ.'CONUS211') THEN
      !  This is the NCEP 211 grid, which is not really 48km
      truedx = 81270.50  ! meters at 25.0 N
      CALL map_set(PROJ_LC,latsw,lonsw,truedx,stdlon,truelat1,truelat2, &
                   nx,ny,proj)
    ELSE IF (mname(1:8).EQ.'CONUS212') THEN
      truedx = 40635.25
      CALL map_set(PROJ_LC,latsw,lonsw,truedx,stdlon,truelat1,truelat2, &
                   nx,ny,proj)
    ELSE IF (mname(1:8).EQ.'CONUS215') THEN
      truedx = 20317.625
      CALL map_set(PROJ_LC,latsw,lonsw,truedx,stdlon,truelat1,truelat2, &
                   nx,ny,proj)
    ELSE IF (mname(1:6).EQ.'LATLON') THEN
      ! True dx for this grid is in lat/lon increment
      truedx = 1.25  ! degrees
      ! For LATLON projections, the latitude increment is stored in truelat1 and
      ! the longitudinal increment is put into stdlon.  For the AVN data on SBN,
      ! The origin in this data set is listed as the NW corner, but the
      ! array is actually set up with the origin at the southwest.
      CALL map_set(PROJ_LATLON,latsw,lonsw,0.,truedx,truedx,0.,nx,ny,proj)
    ELSE
      PRINT *, 'Model type not yet supported: ', TRIM(mname)
      istatus = 0      
      RETURN
    ENDIF

    ! Verify that the computed upper-right corner matches the specified upper right
    ! corner
    CALL ij_to_latlon(proj, FLOAT(nx),FLOAT(ny),latne_c,lonne_c)
    print *, 'Comp latne/lonne = ', latne_c,lonne_c
    print *, 'Spec latne/lonne = ', latne,lonne
    difflat = ABS(latne_c - latne)
    difflon = ABS(lonne_c - lonne)
    IF ((difflat .GT. 0.001).OR.(difflon .GT.0.001)) THEN
      PRINT *, 'Problem with projection information:'
      PRINT *, 'Specified lat/lon at nx/ny = ', latne, lonne 
      PRINT *, 'Computed lat/lon at nx/ny  = ', latne_c, lonne_c
      istatus = 0
      RETURN
    ENDIF                                                                        

    ! If we made it this far, print out some diagnostic information.
    PRINT *, 'WFO model map projection info for ', TRIM(mname)
    PRINT *, '------------------------------------------------------------'
    PRINT *, 'Projection Type:              ', TRIM(gproj)
    PRINT *, 'Southwest Lat/Lon:            ', proj%lat1, proj%lon1
    PRINT *, 'dLat/dLon (LATLON Proj only): ', proj%dlat,proj%dlon
    PRINT *, 'Standard Lon:                 ', proj%stdlon
    PRINT *, 'Truelat1/Truelat2:            ', proj%truelat1,proj%truelat2
    PRINT *, 'Hemi (1 for NH, -1 for SH):   ', proj%hemi
    PRINT *, 'Cone factor:                  ', proj%cone
    PRINT *, 'Pole point (i/j):             ', proj%polei,proj%polej
    RETURN 
  END SUBROUTINE get_wfomodel_proj
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_wfomodel_var_levels(nfid,varname,varid,n_levels,levels_c, &
                                   n_plevels,plevels_r,pbot_ind, ptop_ind, &
                                   havesfc, sfc_ind, &
                                   istatus)
    
    ! Given a file handle for a previously opened netCDF WFO
    ! AWIPS-format model file and a character variable name,
    ! this routine will return various pieces of information
    ! about the variable, including number of total levels,
    ! an array of level IDs, number of pressure levels, array
    ! of pressure levels in Pa, and a flag as to whether or
    ! not there is a surface value.  It also returns the netCDF
    ! integer variable ID.

    IMPLICIT NONE
    INTEGER, INTENT(IN)            :: nfid
    CHARACTER(LEN=10),INTENT(IN)   :: varname
    INTEGER, INTENT(OUT)           :: varid
    INTEGER, INTENT(OUT)           :: n_levels
    CHARACTER(LEN=10),INTENT(OUT)  :: levels_c(maxlevels)
    INTEGER, INTENT(OUT)           :: n_plevels
    REAL, INTENT(OUT)              :: plevels_r(maxlevels)
    INTEGER, INTENT(OUT)           :: ptop_ind
    INTEGER, INTENT(OUT)           :: pbot_ind
    INTEGER, INTENT(OUT)           :: sfc_ind
    LOGICAL, INTENT(OUT)           :: havesfc
    INTEGER, INTENT(OUT)           :: istatus

    CHARACTER(LEN=32)              :: levname
    INTEGER                        :: levid
    INTEGER                        :: dimid(4)
    CHARACTER(LEN=10)              :: sfc_level
    CHARACTER(LEN=10)              :: level_txt
    INTEGER                        :: k,kp
    REAL                           :: press_mb
    CHARACTER(len=2)               :: dummy2
    istatus = 1
    n_levels = 0
    n_plevels = 0
    havesfc = .false.
    varid = -1
    pbot_ind = -1
    ptop_ind = -1
    sfc_ind = -1

    IF ( (TRIM(varname).EQ.'t').OR. &
        ( TRIM(varname).EQ.'rh')) THEN
       sfc_level = 'FHAG 2    '
    ELSE IF ((TRIM(varname).EQ.'uw').OR. &
             (TRIM(varname).EQ.'vw'))THEN
       sfc_level = 'FHAG 10   '
    ELSE IF ((TRIM(varname).EQ.'emsp').OR.&
             (TRIM(varname).EQ.'pmsl').OR.&
             (TRIM(varname).EQ.'mmsp'))THEN
       sfc_level = 'MSL       '
    ELSE
       sfc_level = 'UNKNOWN   '
    ENDIF

    ! Get the variable ID for the variable requested
    nfstatus = NF_INQ_VARID(nfid,varname,varid)
    IF (nfstatus .NE. NF_NOERR) THEN
      PRINT *, 'Variable not found: ', varname
      istatus = 0
      RETURN
    ENDIF

    ! Build the name of the level variable   
    levname = TRIM(varname) //'Levels'
    nfstatus = NF_INQ_VARID(nfid,levname,levid)
    IF (nfstatus .NE. NF_NOERR) THEN
      PRINT *, 'No Levels variable found for ',TRIM(levname)
      istatus = 0
      RETURN
    ENDIF
    nfstatus = NF_INQ_VARDIMID(nfid,levid,dimid)
    nfstatus = NF_INQ_DIMLEN(nfid,dimid(2),n_levels)
    IF (n_levels .GT. 0) THEN
      nfstatus = NF_GET_VAR_TEXT(nfid,levid,levels_c)
      ! Scan for number of levels on pressure surfaces
      DO k = 1, n_levels
        level_txt = levels_c(k)
        IF (level_txt(1:2) .EQ. 'MB') THEN
          n_plevels = n_plevels + 1

          ! Assume the pressure level data is contiguous from
          ! the lower atmosphere upward in the array 
          IF (pbot_ind .LE. 0) pbot_ind = k
        ELSE IF (level_txt .EQ. sfc_level)THEN
            havesfc = .true.
            IF (sfc_ind .LE. 0) sfc_ind = k
        ENDIF
      ENDDO
      IF (n_plevels .GT. 0) THEN
        ptop_ind = pbot_ind + n_plevels-1
        kp = 1
        DO k=1,n_levels
          level_txt = levels_c(k)
          IF (level_txt(1:2).EQ.'MB')THEN
            READ(level_txt,'(A2,1x,F7.0)') dummy2, press_mb
            plevels_r(kp) = press_mb * 100.
            kp = kp + 1
          ENDIF
        ENDDO
      ENDIF
    ENDIF  
    PRINT *, 'Level info for variable ', TRIM(varname)
    PRINT *, '-------------------------------------------------'
    PRINT *, 'Num Total Levels:        ', n_levels
    PRINT *, 'Level IDs: ', levels_c(1:n_levels)
    PRINT *, 'Num Pressure Levels:     ', n_plevels
    PRINT *, 'Pressure levels (Pa):', plevels_r(1:n_plevels)
    PRINT *, 'Found surface value:     ', havesfc
    PRINT *, 'KBOTP/KTOPP/KSFC:   ', pbot_ind,ptop_ind,sfc_ind
    RETURN
  END SUBROUTINE get_wfomodel_var_levels
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_wfomodel_fcsttimes(nfid,ntimes,fcstsec,istatus)
  
    ! Given the integer file handle of a previously opened
    ! AWIPS netCDF model file, this routine returns the  number 
    ! of output forecast times and each valtime-reftime in seconds.

    IMPLICIT NONE
    
    INTEGER, INTENT(IN)                :: nfid
    INTEGER, INTENT(OUT)               :: ntimes
    INTEGER, INTENT(OUT)               :: fcstsec(maxtimes)
    INTEGER, INTENT(OUT)               :: istatus

    INTEGER                            :: timeid
    istatus = 1

    nfstatus = NF_INQ_DIMID(nfid,'n_valtimes',timeid)
    IF (nfstatus .NE. NF_NOERR)  THEN
      PRINT *, 'Could not obtain n_valtimes variable ID'
      istatus = 0
      RETURN
    ENDIF
    nfstatus = NF_INQ_DIMLEN(nfid,timeid,ntimes)
    nfstatus = NF_INQ_VARID(nfid,'valtimeMINUSreftime',timeid)
    IF (nfstatus .NE. NF_NOERR) THEN
      PRINT *, 'Could not obtain valtimeMINUSreftime'
      istatus = 0
      RETURN
    ENDIF
    nfstatus = NF_GET_VARA_INT(nfid,timeid,1,ntimes,fcstsec(1:ntimes))
    IF (nfstatus .NE. NF_NOERR) THEN
      PRINT *, 'Problem obtaining the valid times'
      istatus = 0 
    ENDIF 
    PRINT *, 'Forecast hours found: ', fcstsec(1:ntimes)/3600
    RETURN
  END SUBROUTINE get_wfomodel_fcsttimes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE read_wfomodel_data(nfid,vid,proj,time_ind, &
            start_lev,stop_lev,data,istatus)

    ! Subroutine to read a specific variable (given by integer
    ! variable ID) from a specific file (given by integer file ID)
    ! from the start/stop values for each dimension of the array.   The
    ! routine presumes you have already called open_wfofile and 
    ! get_wfomodel_var_levels for the variable you want to obtain.

    IMPLICIT NONE
    INTEGER, INTENT(IN)                   :: nfid
    INTEGER, INTENT(IN)                   :: vid
    TYPE(proj_info),INTENT(IN)            :: proj
    INTEGER,INTENT(IN)                    :: time_ind
    INTEGER,INTENT(IN)                    :: start_lev
    INTEGER,INTENT(IN)                    :: stop_lev
    REAL,INTENT(OUT)                      :: data(proj%nx,proj%ny, &
                                              stop_lev-start_lev+1)
    INTEGER, INTENT(OUT)                  :: istatus

    INTEGER                               :: start_ind(4)
    INTEGER                               :: count_ind(4)

    istatus = 1
    start_ind(4) = time_ind
    count_ind(4) = 1
    start_ind(3) = start_lev
    count_ind(3) = stop_lev - start_lev + 1
    start_ind(2) = 1
    count_ind(2) = proj%ny
    start_ind(1) = 1
    count_ind(1) = proj%nx

    nfstatus = NF_GET_VARA_REAL(nfid,vid,start_ind,count_ind,data)
    IF (nfstatus .NE. NF_NOERR) THEN
      PRINT *, 'Problem obtaining variable for VID = ', vid
      istatus = 0
    ENDIF
    RETURN
  END SUBROUTINE read_wfomodel_data 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_wfomodel_topo(nfid,proj,topo,istatus)
  
    ! This routine will obtain the topography data from
    ! an AWIPS model file.  It assumes you have already opened
    ! the file to get nfid and called the get_wfomodel_proj routine
    ! to get the proj structure and that topo is already allocated

    IMPLICIT NONE
    INTEGER, INTENT(IN)           :: nfid
    TYPE(proj_info),INTENT(IN)    :: proj
    REAL, INTENT(OUT)             :: topo(proj%nx,proj%ny)
    INTEGER, INTENT(OUT)          :: istatus
    INTEGER                       :: topoid 
    istatus = 1

    nfstatus = NF_INQ_VARID(nfid,'staticTopo',topoid)
    IF(nfstatus.ne.NF_NOERR) then
      print *, NF_STRERROR(nfstatus)
      print *,'NF_INQ_VARID staticTopo'
      istatus = 0
      RETURN                  
    ENDIF

    nfstatus = NF_GET_VAR_REAL(nfid,topoid,topo)
        IF(nfstatus.ne.NF_NOERR) then
      print *, NF_STRERROR(nfstatus)
      print *,'NF_GET_VAR_REAL staticTopo'
      istatus = 0
      RETURN
    ENDIF                   

    RETURN
  END SUBROUTINE get_wfomodel_topo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_wfomodel_var_inv(nfid,varname,nlevs,ntimes,inv,istatus)
  
    ! Subroutine to return the inventory of a specific variable
    ! from an AWIPS netCDF model bigfile.  The inventory is returned
    ! as a 2D array of logical flags dimensioned (n_levels,ntimes)

    IMPLICIT NONE
    INTEGER, INTENT(IN)              :: nfid
    CHARACTER(LEN=10),INTENT(IN)     :: varname
    INTEGER,INTENT(IN)               :: nlevs
    INTEGER,INTENT(IN)               :: ntimes
    LOGICAL,INTENT(OUT)              :: inv(nlevs,ntimes)
    INTEGER,INTENT(OUT)              :: istatus

    CHARACTER(LEN=1),ALLOCATABLE     :: invchar(:,:)
    CHARACTER(LEN=16)                :: invname
    INTEGER                          :: invid

    istatus = 1
    invname = TRIM(varname) //'Inventory'
    nfstatus = NF_INQ_VARID(nfid,invname,invid)
    IF (nfstatus .NE. NF_NOERR) THEN
      PRINT *, 'No variable found: ', TRIM(invname)
      istatus = 0
      RETURN
    ENDIF
  
    ALLOCATE(invchar(nlevs,ntimes)) 
    nfstatus = NF_GET_VAR_TEXT(nfid,invid,invchar)
    IF (nfstatus .NE. NF_NOERR) THEN
      PRINT *, 'Problem getting inventory.'
      istatus = 0
      RETURN
    ENDIF

    inv(:,:) = .false.
    WHERE (invchar .EQ. '1') inv = .true.
    DEALLOCATE(invchar)
    RETURN
  END SUBROUTINE get_wfomodel_var_inv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE wfo_models

                                                                      
