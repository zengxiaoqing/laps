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

MODULE wrfsi_static

  USE time_utils
  ! F90 module to interact with the WRFSI static file, which is found
  ! in MOAD_DATAROOT/static.  
  !
  IMPLICIT NONE

  ! Declare an allocatable array to use for fields which may have a variable
  ! 3rd dimension (such as categorical land use fraction, soil category fraction, etc.) 
  ! that can be allocated and returned to the calling routine
 
   REAL, ALLOCATABLE                 :: static_landusef(:,:,:)
   REAL, ALLOCATABLE                 :: static_soilbot(:,:,:)
   REAL, ALLOCATABLE                 :: static_soiltop(:,:,:)
  


  ! Delcare some variable name IDs.  These are variables that
  ! map the netCDF names based on the CDL.  The naming convention is
  ! as follows:
  !
  ! name_var_N:  Non-staggered variable
  ! name_var_M:  Variable on the mass points of the staggered WRF grid
  ! name_var_U:  Variable on the same points as the u-wind 
  ! name_var_V:     "     "   "   "     "     "  "  v-wind 


  ! Grid spacing stuff
  CHARACTER(LEN=2), PARAMETER  :: NAME_DX = 'Dx'
  CHARACTER(LEN=2), PARAMETER  :: NAME_DY = 'Dy'
  
  ! Projection info stuff
  CHARACTER(LEN=6), PARAMETER  :: NAME_TRUELAT1 = 'Latin1'
  CHARACTER(LEN=6), PARAMETER  :: NAME_TRUELAT2 = 'Latin2'
  CHARACTER(LEN=3), PARAMETER  :: NAME_STDLON = 'LoV'
  CHARACTER(LEN=3), PARAMETER  :: NAME_LAT1 = 'La1'
  CHARACTER(LEN=3), PARAMETER  :: NAME_LON1 = 'Lo1'
  CHARACTER(LEN=9), PARAMETER  :: NAME_TYPE = 'grid_type'
  CHARACTER(LEN=2), PARAMETER  :: NAME_NX = 'Nx'
  CHARACTER(LEN=2), PARAMETER  :: NAME_NY = 'Ny'
  ! Latitude and longitude arrays

  CHARACTER(LEN=3), PARAMETER  :: NAME_LAT_N = 'lat'
  CHARACTER(LEN=3), PARAMETER  :: NAME_LON_N = 'lon'
  CHARACTER(LEN=3), PARAMETER  :: NAME_LAT_T = 'lac'
  CHARACTER(LEN=3), PARAMETER  :: NAME_LON_T = 'loc'
  CHARACTER(LEN=3), PARAMETER  :: NAME_LAT_U = 'lab'
  CHARACTER(LEN=3), PARAMETER  :: NAME_LON_U = 'lob'
  CHARACTER(LEN=3), PARAMETER  :: NAME_LAT_V = 'laa'
  CHARACTER(LEN=3), PARAMETER  :: NAME_LON_V = 'loa'
!mp-BLS
  CHARACTER(LEN=3), PARAMETER  :: NAME_LAT_H= 'lah'
  CHARACTER(LEN=3), PARAMETER  :: NAME_LON_H = 'loh'
  CHARACTER(LEN=3), PARAMETER  :: NAME_LAT_W = 'lav'
  CHARACTER(LEN=3), PARAMETER  :: NAME_LON_W = 'lov'
!mp-BLS

  ! Stuff related to topography
  CHARACTER(LEN=3), PARAMETER  :: NAME_TER_N = 'avg'
  CHARACTER(LEN=3), PARAMETER  :: NAME_TER_T = 'avc'
  CHARACTER(LEN=3), PARAMETER  :: NAME_TERSD_N = 'std'
  CHARACTER(LEN=3), PARAMETER  :: NAME_TERENV_N = 'env'
  CHARACTER(LEN=3), PARAMETER  :: NAME_TERGRADLN_N = 'sln'
  CHARACTER(LEN=3), PARAMETER  :: NAME_TERGRADLT_N = 'slt'

  ! Coriolis
  CHARACTER(LEN=3), PARAMETER  :: NAME_HCOR_T = 'cph'
  CHARACTER(LEN=3), PARAMETER  :: NAME_VCOR_T = 'cpv'

  ! Map factors
  CHARACTER(LEN=3), PARAMETER  :: NAME_MAPFAC_N = 'mfl'
  CHARACTER(LEN=3), PARAMETER  :: NAME_MAPFAC_T = 'mfc'
  CHARACTER(LEN=3), PARAMETER  :: NAME_MAPFAC_U = 'mfb'
  CHARACTER(LEN=3), PARAMETER  :: NAME_MAPFAC_V = 'mfa'

  ! SIN of Alpha angles (angle between longitude and stdlon)
  CHARACTER(LEN=3), PARAMETER  :: NAME_SINALPHA_T = 'spr'
  CHARACTER(LEN=3), PARAMETER  :: NAME_COSALPHA_T = 'cpr'

  ! Land use categories
  CHARACTER(LEN=3), PARAMETER  :: NAME_LANDUSE_T = 'use'

  ! Land mask field
  CHARACTER(LEN=3), PARAMETER  :: NAME_LWMASK_T = 'lnd'
  ! Albedo climatology
  CHARACTER(LEN=3), PARAMETER  :: NAME_ALBEDO_N = 'alb'

  ! Mean annual deep soil temperature on mass grid.
  CHARACTER(LEN=3), PARAMETER  :: NAME_AMT_T = 'tmp'

  INCLUDE "netcdf.inc"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE open_wrfsi_static(dataroot,nestid,cdfid)
  
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN)   :: dataroot
    INTEGER, INTENT(IN)            :: nestid
    INTEGER, INTENT(OUT)           :: cdfid
    CHARACTER(LEN=255)            :: staticfile
    LOGICAL                       :: static_exists
    INTEGER                       :: status
    CHARACTER(LEN=2)              :: nestid_str

    WRITE(nestid_str, '(I2.2)') nestid
    staticfile = TRIM(dataroot) // '/static/static.wrfsi.d' // nestid_str
    INQUIRE(FILE=staticfile, EXIST=static_exists)
    IF (static_exists) THEN
      status = NF_OPEN(TRIM(staticfile),NF_NOWRITE,cdfid)
      IF (status .NE. NF_NOERR) THEN
        PRINT '(A,I5)', 'NetCDF error opening WRF static file: ',status
        STOP 'open_wrfsi_static'
      END IF 
    ELSE
!mp-BLS
!       search for rotlat version??
!      PRINT '(A)', 'Static file not found ', staticfile
!      PRINT '(A)', 'Look for NMM version'
      staticfile = TRIM(dataroot) // '/static/static.wrfsi.rotlat'
      INQUIRE(FILE=staticfile, EXIST=static_exists)
      IF (static_exists) THEN
        status = NF_OPEN(TRIM(staticfile),NF_NOWRITE,cdfid)
        IF (status .NE. NF_NOERR) THEN
          PRINT '(A,I5)', 'NetCDF error opening WRF static file: ',status
          STOP 'open_wrfsi_static'
        END IF
      ELSE

        PRINT '(A)', 'rotlat Static file not found, either: ', staticfile
        STOP 'open_wrfsi_static'
      ENDIF

    ENDIF
    RETURN
  END SUBROUTINE open_wrfsi_static      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_wrfsi_static_dims(dataroot, nestid, nx, ny)
  
    ! Subroutine to return the horizontal dimensions of WRF static file
    ! contained in the input dataroot

    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN)  :: dataroot
    INTEGER         , INTENT(IN)  :: nestid
    INTEGER         , INTENT(OUT) :: nx
    INTEGER         , INTENT(OUT) :: ny

    INTEGER                       :: cdfid,vid, status

    CALL open_wrfsi_static(dataroot,nestid,cdfid)
    status = NF_INQ_DIMID(cdfid, 'x', vid)
    status = NF_INQ_DIMLEN(cdfid, vid, nx)
    status = NF_INQ_DIMID(cdfid, 'y', vid)
    status = NF_INQ_DIMLEN(cdfid, vid, ny) 
    status = NF_CLOSE(cdfid)  
    RETURN
  END SUBROUTINE get_wrfsi_static_dims     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_wrfsi_static_proj(dataroot,nestid,proj_type, lat1, lon1, dx, dy, &
                                   stdlon, truelat1, truelat2)

    ! Returns basic projection information from the WRF static file found
    ! in dataroot

    IMPLICIT NONE
    CHARACTER(LEN=*) , INTENT(IN)      :: dataroot
    INTEGER, INTENT(IN)                :: nestid
    CHARACTER(LEN=32), INTENT(OUT)     :: proj_type
    REAL             , INTENT(OUT)     :: lat1
    REAL             , INTENT(OUT)     :: lon1
    REAL             , INTENT(OUT)     :: dx
    REAL             , INTENT(OUT)     :: dy
    REAL             , INTENT(OUT)     :: stdlon
    REAL             , INTENT(OUT)     :: truelat1
    REAL             , INTENT(OUT)     :: truelat2
   
    INTEGER                            :: cdfid, vid,status
    INTEGER                            :: nx,ny
    REAL,ALLOCATABLE                   :: lats(:,:),lons(:,:)
    CHARACTER(LEN=132)                 :: grid_type
    
    CALL get_wrfsi_static_dims(dataroot,nestid,nx,ny)
    ALLOCATE(lats(nx,ny))
    ALLOCATE(lons(nx,ny))
    CALL open_wrfsi_static(dataroot,nestid,cdfid)   
    status = NF_INQ_VARID( cdfid , NAME_TYPE, vid )   
    status = NF_GET_VAR_TEXT( cdfid , vid , grid_type )  
    IF (grid_type(1:19) .EQ. 'polar stereographic') THEN 
      proj_type = 'POLAR STEREOGRAPHIC             ' 
    ELSE IF ( grid_type(1:24) .EQ. 'secant lambert conformal' ) THEN
      proj_type = 'LAMBERT CONFORMAL               '   
    ELSE IF ( grid_type(1:28) .EQ. 'tangential lambert conformal') THEN
      proj_type = 'LAMBERT CONFORMAL               '  
    ELSE IF (grid_type(1:8)  .EQ. 'mercator'                 ) THEN
      proj_type = 'MERCATOR                        '    
!mp-BLS
    ELSE IF (grid_type(1:15)  .EQ. 'rotated lat-lon'                 ) THEN
      write(6,*) 'setting proj_type to ROTATED LATLON'
      proj_type = 'ROTATED LATLON                  '    
!mp-BLS
    ELSE
      PRINT '(A,A)', 'Unrecognized Projection:', proj_type
      STOP 'GET_WRFSI_STATIC_PROJ'
    END IF                                                       

    ! Get SW Corner lat/lon of projection
    CALL get_wrfsi_static_2d(dataroot,nestid,NAME_LAT_N,lats)
    CALL get_wrfsi_static_2d(dataroot,nestid,NAME_LON_N,lons)
    lat1 = lats(1,1)
    lon1 = lons(1,1)
    IF (lon1 .LT. -180.) lon1 = lon1 + 360.
    IF (lon1 .GT. +180.) lon1 = lon1 - 360.
    PRINT '(A,F10.2,A,F10.2)', 'WRF Lat1 = ',lat1, &
        ' WRF lon1 = ', lon1
    ! Get dx and dy, convert to meters from kilometers

    status =  NF_INQ_VARID( cdfid, NAME_DX, vid )
    status = NF_GET_VAR_REAL(cdfid,vid,dx)
    status = NF_INQ_VARID( cdfid, NAME_DY, vid )
    status = NF_GET_VAR_REAL(cdfid,vid,dy)
    PRINT '(A,F10.2,A,F10.2)', 'WRF Delta-x = ',dx, &
        ' WRF Delta-y = ', dy  
    ! Get standard longitude
    status = NF_INQ_VARID ( cdfid , NAME_STDLON, vid )
    status = NF_GET_VAR_REAL(cdfid , vid , stdlon)    
    IF (stdlon .LT. -180.) stdlon = stdlon + 360.
    IF (stdlon .GT. 180.) stdlon = stdlon - 360.
    PRINT '(A,F10.3)', 'WRF Standard Lon = ', stdlon
    ! Get true latitudes
    status = NF_INQ_VARID(cdfid , NAME_TRUELAT1 , vid)
    status = NF_GET_VAR_REAL( cdfid , vid , truelat1)
    status = NF_INQ_VARID(cdfid , NAME_TRUELAT2, vid)
    status = NF_GET_VAR_REAL( cdfid , vid , truelat2 )
    PRINT '(A,2F10.3)', 'WRF Standard Lats = ', truelat1, truelat2
    status = NF_CLOSE(cdfid)
    DEALLOCATE(lats)
    DEALLOCATE(lons)
    RETURN
  END SUBROUTINE get_wrfsi_static_proj
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_wrfsi_static_latlon(dataroot, nestid,stagger, lat, lon)

    ! Subroutine to get lat/lon arrays for desired grid stagger

    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN)        :: dataroot
    INTEGER , INTENT(IN)                :: nestid
    CHARACTER(LEN=1), INTENT(IN)        :: stagger
    REAL                                :: lat(:,:)
    REAL                                :: lon(:,:)

    INTEGER                             :: status, cdfid, vid
    CHARACTER(LEN=3)                    :: varname_lat, varname_lon

    CALL open_wrfsi_static(dataroot, nestid,cdfid)    

    IF ( (stagger .EQ. 'N').OR.(stagger .EQ. 'n').OR.(stagger.EQ. ' '))THEN
      varname_lat = NAME_LAT_N
      varname_lon = NAME_LON_N
    ELSE IF ( (stagger .EQ. 'T').OR.(stagger.EQ.'t'))THEN
      varname_lat = NAME_LAT_T
      varname_lon = NAME_LON_T  
    ELSE IF ( (stagger .EQ. 'U').OR.(stagger.EQ.'u'))THEN
      varname_lat = NAME_LAT_U
      varname_lon = NAME_LON_U
    ELSE IF ( (stagger .EQ. 'V').OR.(stagger.EQ.'v'))THEN 
      varname_lat = NAME_LAT_V
      varname_lon = NAME_LON_V
!mp-BLS
    ELSE IF ( (stagger .EQ. 'H').OR.(stagger.EQ.'h'))THEN 
      varname_lat = NAME_LAT_H
      varname_lon = NAME_LON_H
    ELSE IF ( (stagger .EQ. 'W').OR.(stagger.EQ.'w'))THEN 
      varname_lat = NAME_LAT_W
      varname_lon = NAME_LON_W
!mp-BLS
    ELSE
      PRINT '(2A)', 'Unrecongized stagger code: ', stagger
      STOP 'GET_WRFSI_STATIC_LATLON'
    ENDIF
  
    status = NF_INQ_VARID(cdfid, varname_lat, vid)
    status = NF_GET_VAR_REAL(cdfid,vid,lat)    
    status = NF_INQ_VARID(cdfid, varname_lon, vid)
    status = NF_GET_VAR_REAL(cdfid,vid,lon) 
  
    status = NF_CLOSE(cdfid) 
    RETURN

  END SUBROUTINE get_wrfsi_static_latlon 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_wrfsi_static_mapfac(dataroot, nestid,stagger, mapfac) 

    ! Subroutine to get alpha arrays for desired grid stagger

    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN)        :: dataroot
    INTEGER, INTENT(IN)                 :: nestid
    CHARACTER(LEN=1), INTENT(IN)        :: stagger
    REAL, INTENT(OUT)                   :: mapfac(:,:)

    INTEGER                             :: cdfid, vid, status
    CHARACTER(LEN=3)                    :: varname_map

    CALL open_wrfsi_static(dataroot,nestid,cdfid)
    IF ( (stagger .EQ. 'N').OR.(stagger .EQ. 'n').OR.(stagger.EQ. ' '))THEN
      varname_map = NAME_MAPFAC_N
    ELSE IF ( (stagger .EQ. 'U').OR.(stagger.EQ.'u'))THEN
      varname_map = NAME_MAPFAC_U
    ELSE IF ( (stagger .EQ. 'V').OR.(stagger.EQ.'v'))THEN
      varname_map = NAME_MAPFAC_V
    ELSE IF ( (stagger .EQ. 'T').OR.(stagger.EQ.'t')) THEN
      varname_map = NAME_MAPFAC_T
    ELSE
      PRINT '(2A)', 'Unrecongized stagger code: ', stagger
      STOP 'GET_WRFSI_STATIC_MAPFAC'
    ENDIF

    status = NF_INQ_VARID(cdfid, varname_map, vid)
    status = NF_GET_VAR_REAL(cdfid,vid,mapfac)

    status = NF_CLOSE(cdfid)
    RETURN

  END SUBROUTINE get_wrfsi_static_mapfac
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_wrfsi_static_landmask(dataroot, nestid,landmask)
 
    ! Subroutine to get the 2D land mask field from the WRFSI static
    ! file
   
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN)        :: dataroot
    INTEGER, INTENT(IN)                 :: nestid
    REAL            , INTENT(OUT)       :: landmask(:,:)

    INTEGER                             :: cdfid, vid, status
    INTEGER                             :: lu_water
    INTEGER                             :: lu_ice
    CHARACTER(LEN=4)                    :: lu_source
    INTEGER                             :: n_cat, dim3,nx,ny
    REAL, ALLOCATABLE                   :: landuse(:,:,:)
    INTEGER                             :: i,j

    CALL open_wrfsi_static(dataroot,nestid,cdfid)  
    status = NF_INQ_VARID(cdfid, NAME_LWMASK_T, vid)
    print '(a,2i4)','landmask inq_varid status and vid = ' ,status,vid
    IF (status .NE. NF_ENOTVAR) THEN
      status = NF_GET_VAR_REAL(cdfid,vid,landmask)
      landmask = FLOAT(NINT(landmask))
      IF ( (MINVAL(landmask) .LT. 0.) .OR. &
           (MAXVAL(landmask) .GT. 1.) ) THEN
        PRINT '(A,2F12.8)', 'Landmask values LT 0 or GT 1 found!:', &
          minval(landmask), maxval(landmask)
          STOP 'GET_WRFSI_STATIC_LANDMASK'
      ENDIF
      landmask = FLOAT(NINT(landmask)) ! Prevents precision problems
    ELSE
      PRINT '(A)','Problem getting WRF landmask from static.'
      PRINT '(A)', '(Old version of static file?)'
      PRINT '(A)', 'Attempting to use land use...'
      CALL get_wrfsi_static_landuse_info(dataroot,nestid,lu_source,n_cat,dim3, &
                                         lu_water,lu_ice)
      CALL get_wrfsi_static_dims(dataroot,nestid, nx,ny)
      ALLOCATE (landuse(nx,ny,dim3))
      CALL get_wrfsi_static_landuse(dataroot, nestid,landuse)
      IF (dim3 .EQ. 1) THEN
        PRINT '(A)', 'Using dominant category version...'
        WHERE(NINT(landuse(:,:,1)) .NE. lu_water) landmask = 1.
        WHERE(NINT(landuse(:,:,1)) .EQ. lu_water) landmask = 0.
      ELSE IF (dim3 .EQ. n_cat) THEN
        PRINT '(A)','Using fractional category version with a 50% water threshold...'
        landmask = 0.
        DO j = 1,ny
          DO i=1,nx
            IF (landuse(i,j,lu_water).LT. 0.5) landmask(i,j) = 1.
          ENDDO
        ENDDO
      ENDIF
      DEALLOCATE(landuse)
    ENDIF
    status = NF_CLOSE(cdfid)
    RETURN
  END SUBROUTINE get_wrfsi_static_landmask
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_wrfsi_static_landuse_info(dataroot, nestid, source, n_categories, dim3, &
                                           water_ind,ice_ind )
 
    ! Subroutine to query what type of land use we have in the static file
    !  Source is the source of data (only 'USGS' is supported for now)
    !  n_categories tells how many different categories
    !  dim3 tells the third dimension of the landuse data array:
    !     landuse(nx,ny,dim3), where dim3 = 1 if it is a dominant
    !      category value or dim3 = n_categories if it is 
    !      fraction by category
    !  water_ind is the index value for the water type
    !  ice_ind is the index value for the ice type
   
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN)        :: dataroot
    INTEGER, INTENT(IN)                 :: nestid
    CHARACTER(LEN=4), INTENT(OUT)       :: source
    INTEGER         , INTENT(OUT)       :: n_categories
    INTEGER         , INTENT(OUT)       :: dim3
    INTEGER         , INTENT(OUT)       :: water_ind
    INTEGER         , INTENT(OUT)       :: ice_ind

    INTEGER                             :: cdfid, vid, status
    INTEGER                             :: dimid(4)
   

    CALL open_wrfsi_static(dataroot,nestid,cdfid)
    ! Get land use source (hard coded for USGS 24 cat)
    source = 'USGS'

    ! Get number of categories
    n_categories = 24
    
    ! Get value that represents water (hard coded for USGS 24 cat)
    water_ind = 16
    
    ! Get value that represents ice (hard coded for USGS 24 cat)
    ice_ind = 24

    ! Query for dimension info
    status = NF_INQ_VARID(cdfid,NAME_LANDUSE_T,vid)
    status = NF_INQ_VARDIMID(cdfid,vid,dimid)
    status = NF_INQ_DIMLEN(cdfid,dimid(3),dim3)
    
    PRINT '(A)', 'STATIC FILE LANDUSE INFO:'
    PRINT '(2A)', 'SOURCE: ', source
    PRINT '(A,I4)', 'NUMBER OF CATEGORIES: ',n_categories
    PRINT '(A,I4)', 'NUMBER OF LEVELS IN ARRAY: ', dim3
    PRINT '(A,I4)', 'WATER VALUE: ', water_ind
    PRINT '(A,I4)', 'ICE VALUE: ', ice_ind
    status = NF_CLOSE(cdfid)
    RETURN
  END SUBROUTINE get_wrfsi_static_landuse_info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_wrfsi_static_landuse(dataroot,nestid,landuse)
    
    ! Gets the landuse data
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN)  :: dataroot
    INTEGER, INTENT(IN)           :: nestid
    REAL, INTENT(INOUT)             :: landuse(:,:,:)
 
    INTEGER                             :: cdfid, vid, status
   
    CALL open_wrfsi_static(dataroot,nestid,cdfid)
    status = NF_INQ_VARID(cdfid,NAME_LANDUSE_T,vid)
    status = NF_GET_VAR_REAL(cdfid,vid,landuse)
    IF (status .NE. NF_NOERR) THEN
      PRINT '(A)', 'Problem getting land use data.'
    ENDIF  
    status = NF_CLOSE(cdfid)
    RETURN
  END SUBROUTINE get_wrfsi_static_landuse  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_wrfsi_static_landusef(dataroot, nestid,source,ncat, water, ice)
    
    ! Reads the individual 2D categorical landuse fraction arrays from
    ! the static file and populates the 3D variable delcared at the top of this
    ! module (static_landusef)

    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN)         :: dataroot
    INTEGER, INTENT(IN)                  :: nestid
    CHARACTER(LEN=4), INTENT(OUT)        :: source
    INTEGER         , INTENT(OUT)        :: ncat
    INTEGER         , INTENT(OUT)        :: water
    INTEGER         , INTENT(OUT)        :: ice
    INTEGER                              :: dom_nx,dom_ny
    INTEGER                              :: cdfid,vid,status
    INTEGER                              :: cat
    CHARACTER(LEN=3)                     :: vname
 
 
    ! For now, we only support 16-category WMO/FAO data set, so hard code the
    ! info
    source = 'USGS'
    ncat = 24
    water = 16
    ice = 24
    PRINT '(A)', 'Attempting to read categorical landuse fractions...'
    PRINT '(2A)','Source = ', source
    PRINT '(A,I2)','Number of categories = ', ncat
    PRINT '(A,I2)','Index used for water = ', water
    PRINT '(A,I2)','Index used for ice = ', ice
    CALL get_wrfsi_static_dims(dataroot, nestid,dom_nx, dom_ny)
    IF (ALLOCATED(static_landusef)) DEALLOCATE(static_landusef)
    ALLOCATE(static_landusef(dom_nx,dom_ny,ncat))
    CALL open_wrfsi_static(dataroot,nestid,cdfid)      
    DO cat = 1, ncat
      WRITE(vname,'("u",I2.2)') cat 
      status = NF_INQ_VARID(cdfid, vname, vid) 
      status = NF_GET_VAR_REAL(cdfid,vid,static_landusef(:,:,cat))
      IF (status .NE. NF_NOERR) THEN
        PRINT '(A)', 'Problem getting landuse category #',cat
      ENDIF
    ENDDO
    status = NF_CLOSE(cdfid)
    RETURN
  END SUBROUTINE get_wrfsi_static_landusef
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  SUBROUTINE get_wrfsi_static_soil(dataroot, nestid,source,ncat, st_water)
    
    ! Reads the 2-layer soil categorical fractions and populates the arrays
    ! static_soiltop and static_soilbot

    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN)         :: dataroot 
    INTEGER, INTENT(IN)                  :: nestid
    CHARACTER(LEN=4), INTENT(OUT)        :: source
    INTEGER         , INTENT(OUT)        :: ncat
    INTEGER         , INTENT(OUT)        :: st_water

    INTEGER                              :: dom_nx,dom_ny
    INTEGER                              :: cdfid,vid,status
    INTEGER                              :: cat
    CHARACTER(LEN=3)                     :: vname
 
 
    ! For now, we only support 16-category WMO/FAO data set, so hard code the
    ! info
    source = 'FAO '
    ncat = 16
    st_water = 14
    PRINT '(A)', 'Attempting to read categorical soil type fractions...'
    PRINT '(2A)','Source = ', source
    PRINT '(A,I2)','Number of categories = ', ncat
    CALL get_wrfsi_static_dims(dataroot, nestid,dom_nx, dom_ny)
    IF (ALLOCATED(static_soiltop)) DEALLOCATE(static_soiltop)
    IF (ALLOCATED(static_soilbot)) DEALLOCATE(static_soilbot)
    ALLOCATE(static_soiltop(dom_nx,dom_ny,ncat))
    ALLOCATE(static_soilbot(dom_nx,dom_ny,ncat))
    CALL open_wrfsi_static(dataroot,nestid,cdfid)      
    DO cat = 1, ncat
      WRITE(vname,'("b",I2.2)') cat 
      status = NF_INQ_VARID(cdfid, vname, vid) 
      status = NF_GET_VAR_REAL(cdfid,vid,static_soilbot(:,:,cat))
      IF (status .NE. NF_NOERR) THEN
        PRINT '(A)', 'Problem getting bottom soil category #',cat
      ENDIF
      WRITE(vname,'("t",I2.2)') cat 
      status = NF_INQ_VARID(cdfid, vname, vid) 
      status = NF_GET_VAR_REAL(cdfid,vid,static_soiltop(:,:,cat))
      IF (status .NE. NF_NOERR) THEN
        PRINT '(A)', 'Problem getting top soil category #',cat
      ENDIF
    ENDDO
    status = NF_CLOSE(cdfid)
    RETURN
  END SUBROUTINE get_wrfsi_static_soil  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_wrfsi_static_monthly(dataroot,nestid, dtypechar, time, data)
   
    ! Returns a time-interpolated (valid for time) 2D array for either
    ! greenness ("g") or albedo ("a"), based on user-supplied type character.
    ! The monthly values are 
    ! valid on the 15th day of each month.  This routine only interpolates to the
    ! nearest day and does not account for leap years, but this should not be any
    ! big deal.  

    IMPLICIT NONE
    CHARACTER(LEN=*) , INTENT(IN)               :: dataroot
    INTEGER, INTENT(IN)                         :: nestid
    CHARACTER(LEN=19), INTENT(IN)               :: time    ! yyyy-mm-dd_hh:mm:ss
    CHARACTER(LEN=1) , INTENT(IN)               :: dtypechar
    REAL             , INTENT(OUT)              :: data(:,:)
    
    INTEGER                                     :: midmonth_day(12)
    INTEGER                                     :: valid_day
    INTEGER                                     :: yyyyjjj
    REAL                                        :: sss
    INTEGER                                     :: m, d1, d2, m1, m2
    INTEGER                                     :: nx,ny
    CHARACTER(LEN=3)                            :: varname
    REAL, ALLOCATABLE                           :: data1(:,:), data2(:,:)
    REAL                                        :: w1, w2

    ! midmonth_day is the julian day of the year corresponding to the 15th day
    ! of each month for a standard (non-leap) year
    DATA midmonth_day / 15, 43, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349 /

    ! Check data type character to make sure it is either greenness or albedo.
    IF ((dtypechar .NE. "a").AND.(dtypechar.NE."g")) THEN
      PRINT *, 'Unknown data type character passed into get_wrfsi_static_monthly:', &
               dtypechar
      PRINT *,'Current supported values are a (albedo) and g (greenness fraction).'
      STOP 'get_Wrfsi_monthly'
    ELSE IF (dtypechar .EQ. 'a') THEN
      PRINT *, 'Getting time-interpolated albedo.'
    ELSE IF (dtypechar .EQ. 'g') THEN
      PRINT *, 'Getting time-interpolated greenness.'
    ENDIF

    ! Convert date string into integer yyyyjjj and sss
    CALL mm5_to_wrf_date(time, yyyyjjj, sss)
    valid_day = MOD(yyyyjjj,1000)
    PRINT *, 'Time-interpolating to day # ', valid_day
    ! Find bounding months
    IF ((valid_day .LT. midmonth_day(1)) .OR. (valid_day .GT. midmonth_day(12))) THEN
      ! December and January are bounding months
      d1 = midmonth_day(12)
      d2 = midmonth_day(1)
      m1 = 12
      m2 = 1
    ELSE
      find_bounds: DO m = 1, 11
        d1 = midmonth_day(m)
        d2 = midmonth_day(m+1)
        IF (valid_day .EQ. d1) THEN
           d2 = d1
           m1 = m
           m2 = m1
           EXIT find_bounds
        ELSE IF (valid_day .EQ. d2) THEN
           d1 = d2
           m1 = m + 1
           m2 = m1
           EXIT find_bounds
        ELSE IF ((valid_day .GT. d1).AND.(valid_day .LT. d2)) THEN
           m1 = m
           m2 = m + 1
           EXIT find_bounds
        ENDIF
      ENDDO find_bounds
    ENDIF

    ! If d1 = d2, then we don't need any interpolation, just get that month's 
    ! data values
    IF ( d1 .EQ. d2) THEN
      WRITE(varname, '(A1,I2.2)') dtypechar,m1
      CALL get_wrfsi_static_2d(dataroot,nestid,varname,data)
    ELSE
      ! We need to get the two months of bounding data and time interpolate
      CALL get_wrfsi_static_dims(dataroot,nestid, nx, ny)
      ALLOCATE(data1 (nx,ny))
      ALLOCATE(data2 (nx,ny))
      WRITE(varname, '(A1,I2.2)') dtypechar,m1
      CALL get_wrfsi_static_2d(dataroot,nestid,varname,data1)
      WRITE(varname, '(A1,I2.2)') dtypechar,m2
      CALL get_wrfsi_static_2d(dataroot,nestid,varname,data2)

      ! Compute weights
      IF (d2 .GT. d1) THEN
        w1 = ( FLOAT(d2) - FLOAT(valid_day) ) / FLOAT(d2-d1)
      ELSE ! We must be between Dec 15 and Jan 15
        IF (valid_day .LT. midmonth_day(1)) THEN ! We are in January
           w1 = ( FLOAT(d2) - FLOAT(valid_day) ) / 31.
        ELSE ! We are in December
           w1 = ( 366. - FLOAT(valid_day) + FLOAT(midmonth_day(1)) ) / 31.
        ENDIF
      ENDIF
      w2 = 1. - w1
      data = w1*data1 + w2*data2
      DEALLOCATE(data1)
      DEALLOCATE(data2)
    ENDIF
    RETURN
  END SUBROUTINE get_wrfsi_static_monthly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_wrfsi_static_monthly_all(dataroot, nestid, dtypechar, data)
    IMPLICIT NONE
    CHARACTER(LEN=*) , INTENT(IN)               :: dataroot
    INTEGER, INTENT(IN)                         :: nestid
    CHARACTER(LEN=1) , INTENT(IN)               :: dtypechar
    REAL             , INTENT(OUT)              :: data(:,:,:)
    INTEGER                                     :: m,nx,ny
    REAL, ALLOCATABLE                           :: month_slab(:,:)
    CHARACTER(LEN=3)                            :: varname

    ! Allocate the dummy array to hold each months data
    CALL get_wrfsi_static_dims(dataroot,nestid, nx, ny)
    ALLOCATE (month_slab(nx,ny))
    month_slab(:,:) = 0.
    month_loop: DO m = 1, 12  
      WRITE (varname, '(a1,I2.2)') dtypechar, m
      CALL get_wrfsi_static_2d(dataroot,nestid,varname,month_slab)
      data(:,:,m) = month_slab(:,:)
      month_slab(:,:) = 0.
    ENDDO month_loop
    DEALLOCATE(month_slab)
    RETURN
  END SUBROUTINE get_wrfsi_static_monthly_all
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_wrfsi_static_2d(dataroot, nestid,varname, data)

    ! Gets any 2D variable from the static file
    CHARACTER(LEN=*), INTENT(IN)  :: dataroot
    INTEGER, INTENT(IN)           :: nestid
    CHARACTER(LEN=*), INTENT(IN)  :: varname
    REAL, INTENT(OUT)             :: data(:,:)
 
    INTEGER                             :: cdfid, vid, status
   
    CALL open_wrfsi_static(dataroot,nestid,cdfid)
    status = NF_INQ_VARID(cdfid,varname,vid)
    status = NF_GET_VAR_REAL(cdfid,vid,data)
    IF (status .NE. NF_NOERR) THEN
      PRINT '(A)', 'Problem getting 2D data.'
    ENDIF 
    status = NF_CLOSE(cdfid) 
    RETURN
  END SUBROUTINE get_wrfsi_static_2d    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
END MODULE wrfsi_static
