MODULE laps_static

  IMPLICIT NONE
  INTEGER :: x, y, z2, z3
  REAL    :: la1, lo1, dx, dy, lov, latin1, latin2
  REAL, ALLOCATABLE :: topo(:,:)
  CHARACTER (LEN=132) grid_type

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

    ! Get the various projection parameters

    vid = NCVID ( cdfid , 'La1' , rcode )
    CALL NCVGT1 ( cdfid , vid , 1 , la1 , rcode )
    IF ( la1 .GT. 270 ) THEN
      PRINT '(A,F9.4,A,F9.4,A)', &
         'Weird lower left corner latitude in LAPS data: ',la1, &
         '.  Replacing it with ',la1-360,'.'
      la1=la1-360
    END IF  

    vid = NCVID ( cdfid , 'Lo1' , rcode )
    CALL NCVGT1 ( cdfid , vid , 1 , lo1 , rcode )
    IF ( lo1 .GT. 180 ) THEN
      lo1 = lo1 - 360
    END IF

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
