MODULE laps_static

  IMPLICIT NONE
  INTEGER                     :: x              ! X-dimension of grid
  INTEGER                     :: y              ! Y-dimension of grid
  INTEGER                     :: z2             ! Z-dimension of 2-d grids
  INTEGER                     :: z3             ! Z-dimension of 3-d press
  REAL                        :: la1            ! SW Corner Lat
  REAL                        :: lo1            ! SE Corner Lon
  REAL                        :: la2            ! NE Corner Lat
  REAL                        :: lo2            ! NE Corner Lon
  REAL                        :: dx             ! X-direction grid spacing
  REAL                        :: dy             ! Y-direction grid spacing
  REAL                        :: lov            ! Orientation Longitude
  REAL                        :: latin1         ! Standard Lat 1
  REAL                        :: latin2         ! Standard Lat 2
  REAL, ALLOCATABLE           :: topo(:,:)      ! LAPS Topographic height
  REAL, ALLOCATABLE           :: lats(:,:)      ! LAPS Latitude Array
  REAL, ALLOCATABLE           :: lons(:,:)      ! LAPS Longitude Array
  CHARACTER (LEN=132)         :: grid_type      ! Map projection type

CONTAINS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE get_horiz_grid_spec(laps_data_root)

    IMPLICIT NONE
    CHARACTER (LEN=*), INTENT(IN) :: laps_data_root
    CHARACTER (LEN=256)           :: static_file
    CHARACTER (LEN=132)           :: dum
    INTEGER :: cdfid, rcode,xid,yid,vid
    INTEGER, DIMENSION(2) :: startc, countc
    INTEGER, DIMENSION(4) :: start, count
    INCLUDE "netcdf.inc"

    static_file = TRIM(laps_data_root)//'/static/static.nest7grid'

    ! Open the static file which is a netCDF file

    cdfid = NCOPN(TRIM(static_file),NCNOWRIT, rcode)

    ! Get x/y dimensions.  
 
    xid = NCDID (cdfid, 'x', rcode)
    yid = NCDID (cdfid, 'y', rcode)
    CALL NCDINQ ( cdfid, xid, dum, x, rcode )
    CALL NCDINQ ( cdfid, yid, dum, y, rcode )
    
    ! Get the grid projection type

    startc = (/ 1 , 1 /)
    countc = (/ 132 , 1 /)
    vid = NCVID ( cdfid , 'grid_type' , rcode )
    CALL NCVGTC ( cdfid , vid , startc , countc , grid_type , 132 , rcode ) 

    ! Get the corner points by reading lat/lon array
    ALLOCATE ( lats (x,y) )
    ALLOCATE ( lons (x,y) )
    vid = NCVID (cdfid, 'lat', rcode)
    start = (/ 1 , 1 , 1 , 1 /)
    count = (/ x , y , 1 , 1 /)
    CALL NCVGT ( cdfid , vid , start , count , lats, rcode ) 
    vid = NCVID (cdfid, 'lon', rcode)
    start = (/ 1 , 1 , 1 , 1 /)
    count = (/ x , y , 1 , 1 /)
    CALL NCVGT ( cdfid , vid , start , count , lons, rcode ) 
    la1 = lats(1,1)
    lo1 = lons(1,1)
    la2 = lats(x,y)
    lo2 = lons(x,y)
    IF (lo1 .gt. 180) lo1 = lo1 - 360.
    IF (lo2 .gt. 180) lo2 = lo2 - 360.
    IF (la1 .gt. 270) la1 = la1 - 360.
    IF (la2 .gt. 270) la2 = la2 - 360.
    print '(A,2F10.2)', 'SW Corner Lat/Lon = ', la1, lo1
    print '(A,2F10.2)', 'NE Corner Lat/Lon = ', la2, lo2

    ! Get other projection parameters

    vid = NCVID ( cdfid , 'Dx' , rcode )
    CALL NCVGT1 ( cdfid , vid , 1 , dx , rcode )
    dx = dx / 1000

    vid = NCVID ( cdfid , 'Dy' , rcode )
    CALL NCVGT1 ( cdfid , vid , 1 , dy , rcode )
    dy = dy / 1000

    vid = NCVID ( cdfid , 'LoV' , rcode )
    CALL NCVGT1 ( cdfid , vid , 1 , lov , rcode )
    IF ( lov .GT. 180 ) THEN
      lov = lov - 360
    END IF

    vid = NCVID ( cdfid , 'Latin1' , rcode )
    CALL NCVGT1 ( cdfid , vid , 1 , latin1 , rcode )

    vid = NCVID ( cdfid , 'Latin2' , rcode )
    CALL NCVGT1 ( cdfid , vid , 1 , latin2 , rcode )

    ! Get the topography

    ALLOCATE(topo(x,y))
    vid = NCVID (cdfid, 'avg', rcode)
    start = (/ 1 , 1 , 1 , 1 /)
    count = (/ x , y , 1 , 1 /)
    CALL NCVGT ( cdfid , vid , start , count , topo, rcode ) 
  
  END SUBROUTINE get_horiz_grid_spec
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE laps_static
