MODULE grid_utils
  
! This module contains utilities related to grid manipulation
! (e.g., creating data for the B and C grids on an Arakawa-C
!  type stagger when only the A grid is available)


CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE arakawa_c_n2t(datain, nx, ny, nz, dataout)
    
    ! Staggers a 3D array of data from the non-staggered points
    ! to the mass grid of an Arakawa C stagger.

    IMPLICIT NONE
    INTEGER, INTENT(IN)                :: nx
    INTEGER, INTENT(IN)                :: ny
    INTEGER, INTENT(IN)                :: nz
    REAL, INTENT(IN)                   :: datain(nx,ny,nz)
    REAL, INTENT(OUT)                  :: dataout(nx,ny,nz)
   
    INTEGER                            :: i,j,k
    PRINT *, 'Staggering to T grid (Arakawa C)'
    DO k = 1, nz
      DO j = 1, ny-1
        DO i = 1, nx-1
          dataout(i,j,k) = 0.25 * ( datain(i,j,k)    + &
                                    datain(i+1,j,k)  + &
                                    datain(i+1,j+1,k)+ &
                                    datain(i,j+1,k) )
        ENDDO
        ! Fill unused rightmost column
        dataout(nx,j,k) = dataout(nx-1,j,k)
      ENDDO
      ! Fill unused uppermost row
      dataout(:,ny,k) = dataout(:,ny-1,k)
    ENDDO
    RETURN
  END SUBROUTINE arakawa_c_n2t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE arakawa_c_n2u(datain, nx, ny, nz, dataout)

    ! Staggers a 3D array of data from the non-staggered points
    ! to the U grid of an Arakawa C stagger.

    IMPLICIT NONE
    INTEGER, INTENT(IN)                :: nx
    INTEGER, INTENT(IN)                :: ny
    INTEGER, INTENT(IN)                :: nz
    REAL, INTENT(IN)                   :: datain(nx,ny,nz)
    REAL, INTENT(OUT)                  :: dataout(nx,ny,nz)

    INTEGER                            :: i,j,k
    PRINT *, 'Staggering to U grid (Arakawa C)'
    DO k = 1, nz
      DO j = 1, ny-1
        DO i = 1, nx
          dataout(i,j,k) = 0.50 * ( datain(i,j,k)    + &
                                    datain(i,j+1,k) )
        ENDDO
      ENDDO
      ! Fill unused uppermost row
      dataout(:,ny,k) = dataout(:,ny-1,k)
    ENDDO
    RETURN
  END SUBROUTINE arakawa_c_n2u                    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE arakawa_c_n2v(datain, nx, ny, nz, dataout)

    ! Staggers a 3D array of data from the non-staggered points
    ! to the V grid of an Arakawa C stagger.

    IMPLICIT NONE
    INTEGER, INTENT(IN)                :: nx
    INTEGER, INTENT(IN)                :: ny
    INTEGER, INTENT(IN)                :: nz
    REAL, INTENT(IN)                   :: datain(nx,ny,nz)
    REAL, INTENT(OUT)                  :: dataout(nx,ny,nz)

    INTEGER                            :: i,j,k
    PRINT *, 'Staggering to V Grid (Arakawa C)'
    DO k = 1, nz
      DO j = 1, ny
        DO i = 1, nx-1
          dataout(i,j,k) = 0.50 * ( datain(i,j,k)    + &
                                    datain(i+1,j,k) )
        ENDDO
        ! Fill unused right column
        dataout(nx,j,k) = dataout(nx-1,j,k)
      ENDDO
    ENDDO
    RETURN
  END SUBROUTINE arakawa_c_n2v            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE arakawa_c_t2n(datain,nx_t,ny_t,nz,dataout)

    ! Destaggers a WRF Arakawa C from the staggered thermodynamic
    ! points to the non-staggered points.   Note that this routine
    ! returns an array that is one element larger in each direction
    ! than the input array.
    !
    !  Example:  Input 3x3 "T" points, return 4x4 "N" points
    !
    !           N   N   N   N
    !             T   T   T
    !           N   N   N   N 
    !             T   T   T
    !           N   N   N   N
    !             T   T   T 
    !           N   N   N   N

    IMPLICIT NONE

    INTEGER, INTENT(IN)       :: nx_t
    INTEGER, INTENT(IN)       :: ny_t    
    INTEGER, INTENT(IN)       :: nz
    REAL, INTENT(IN)          :: datain(nx_t,ny_t,nz)       
    REAL, INTENT(OUT)         :: dataout(nx_t+1,ny_t+1,nz)   

    INTEGER                   :: i,j,k

    vertical_loop:  DO k=1,nz

      ! First, compute all of the interior points

      DO j = 2, ny_t
        DO i = 1, nx_t

           dataout(i,j,k) = 0.25*( datain(i-1,j-1,k) + datain(i-1,j,k) + &
                                  datain(i,j,k) + datain(i,j-1,k) )

        ENDDO
      ENDDO

      ! Now, extrapolate upper and lower rows, except corner points

      DO i = 2, nx_t
 
        dataout(i,1,k) = 2.0* dataout(i,2,k)-dataout(i,3,k)
        dataout(i,ny_t+1,k) = 2.0*dataout(i,ny_t,k)-dataout(i,ny_t-1,k)

      ENDDO

      ! Extrapolate left and right columns, except corner points

      DO j = 2, ny_t
    
        dataout(1,j,k) = 2.0*dataout(2,j,k)-dataout(3,j,k)
        dataout(nx_t+1,j,k) = 2.0*dataout(nx_t,j,k)-dataout(nx_t-1,j,k)

      ENDDO

      ! Compute corner point values by solving for 4 point average

      dataout(1,1,k) = 4.0 * datain(1,1,k) - &
                             dataout(1,2,k) - &
                             dataout(2,2,k) - &
                             dataout(2,1,k)
      dataout(1,ny_t+1,k) = 4.0 * datain(1,ny_t,k) - &
                                  dataout(1,ny_t,k) - &
                                  dataout(2,ny_t+1,k) - &
                                  dataout(2,ny_t,k)
      dataout(nx_t+1,ny_t+1,k) = 4.0 * datain(nx_t,ny_t,k) - &
                                       dataout(nx_t,ny_t,k) - &
                                       dataout(nx_t,ny_t+1,k) - &
                                       dataout(nx_t+1,ny_t,k)
      dataout(nx_t+1,1,k) = 4.0 * datain(nx_t,1,k) - &
                                dataout(nx_t,1,k) - &
                                dataout(nx_t,2,k) - &
                                dataout(nx_t+1,2,k)

    ENDDO vertical_loop       
    RETURN
  END SUBROUTINE arakawa_c_t2n
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE arakawa_c_u2n(datain,nx_u,ny_u,nz,dataout)

    ! Destaggers from the WRF "u" grid to the non-staggered grid.  The
    ! return array will be 1 element larger in the Y direction than
    ! the input array.  Example:
    !
    !   Input 4x3 "U" grid, return 4x4 "N" grid:
    !
    !            N   N   N   N
    !            U   U   U   U
    !            N   N   N   N
    !            U   U   U   U
    !            N   N   N   N
    !            U   U   U   U
    !            N   N   N   N

    IMPLICIT NONE

    INTEGER,  INTENT(IN)   :: nx_u
    INTEGER,  INTENT(IN)   :: ny_u
    INTEGER,  INTENT(IN)   :: nz
    REAL,     INTENT(IN)   :: datain(nx_u,ny_u,nz)
    REAL,     INTENT(OUT)  :: dataout(nx_u,ny_u+1,nz)

    INTEGER                :: i,k
    INTEGER                :: nx,ny

    nx = nx_u
    ny = ny_u+1

    DO k = 1, nz
   
      ! Linear interpolation along each column, except top/bottom rows
      DO i = 1, nx
   
        ! Average of points above and below to fill interior rows 
        dataout(i,2:ny_u,k) = 0.5*(datain(i,1:ny_u-1,k)+datain(i,2:ny_u,k))
 
        ! Fill bottom row
        dataout(i,1,k) = 2.0 * dataout(i,2,k) - dataout(i,3,k)

        ! Fill top row
        dataout(i,ny,k) = 2.0 * dataout(i,ny-1,k) - dataout(i,ny-2,k)

      ENDDO
    ENDDO 
    RETURN
  END SUBROUTINE arakawa_c_u2n
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE arakawa_c_v2n(datain,nx_v,ny_v,nz,dataout)

    ! Destaggers from the WRF "v" grid to the non-staggered grid.  The
    ! return array will be 1 element larger in the X direction than
    ! the input array.  Example:
    !
    !   Input 3x4 "V" grid, return 4x4 "N" grid:
    !
    !            N V N V N V N
    !           
    !            N V N V N V N
    !       
    !            N V N V N V N
    !        
    !            N V N V N V N

    IMPLICIT NONE
    INTEGER,  INTENT(IN)   :: nx_v
    INTEGER,  INTENT(IN)   :: ny_v
    INTEGER,  INTENT(IN)   :: nz
    REAL,     INTENT(IN)   :: datain(nx_v,ny_v,nz)
    REAL,     INTENT(OUT)  :: dataout(nx_v+1,ny_v,nz)

    INTEGER                :: j,k
    INTEGER                :: nx,ny

    nx = nx_v + 1
    ny = ny_v

    DO k = 1, nz

      ! Linear interpolation along each row, except left/right columns
      DO j = 1, ny

        ! Average of points above and below to fill interior rows 
        dataout(2:nx_v,j,k) = 0.5*(datain(1:nx_v-1,j,k)+datain(2:nx_v,j,k))

        ! Fill left column
        dataout(1,j,k) = 2.0 * dataout(2,j,k) - dataout(3,j,k)

        ! Fill top row
        dataout(nx,j,k) = 2.0 * dataout(nx-1,j,k) - dataout(nx-2,j,k)

      ENDDO
    ENDDO
    RETURN
  END SUBROUTINE arakawa_c_v2n
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE arakawa_c_v2t(datain,nx_v,ny_v,nz,dataout)

    ! Destaggers from the WRF "v" grid to the mass grid.  The
    ! return array will be 1 element smaller in the Y direction than
    ! the input array.  Example:
    !
    !   Input 3x4 "V" grid, return 3x3 "T" grid:
    !
    !              V   V   V  
    !              T   T   T
    !              V   V   V  
    !              T   T   T
    !              V   V   V  
    !              T   T   T
    !              V   V   V  

    IMPLICIT NONE
    INTEGER,  INTENT(IN)   :: nx_v
    INTEGER,  INTENT(IN)   :: ny_v
    INTEGER,  INTENT(IN)   :: nz
    REAL,     INTENT(IN)   :: datain(nx_v,ny_v,nz)
    REAL,     INTENT(OUT)  :: dataout(nx_v,ny_v-1,nz)

    INTEGER                :: j,k
    INTEGER                :: nx,ny

    nx = nx_v 
    ny = ny_v - 1

    DO k = 1, nz

      ! Linear interpolation along each column
      DO j = 1, ny

        dataout(:,j,k) = 0.5*(datain(:,j,k) + datain(:,j+1,k))

      ENDDO
    ENDDO
    RETURN
  END SUBROUTINE arakawa_c_v2t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE arakawa_c_u2t(datain,nx_u,ny_u,nz,dataout)

    ! Destaggers from the WRF "u" grid to the mass grid.  The
    ! return array will be 1 element smaller in the X direction than
    ! the input array.  Example:
    !
    !   Input 4x3 "U" grid, return 3x3 "T" grid:
    !
    !                           
    !            U T U T U T U
    !                             
    !            U T U T U T U
    !                             
    !            U T U T U T U
    !                           
    IMPLICIT NONE
    INTEGER,  INTENT(IN)   :: nx_u
    INTEGER,  INTENT(IN)   :: ny_u
    INTEGER,  INTENT(IN)   :: nz
    REAL,     INTENT(IN)   :: datain(nx_u,ny_u,nz)
    REAL,     INTENT(OUT)  :: dataout(nx_u-1,ny_u,nz)

    INTEGER                :: i,k
    INTEGER                :: nx,ny

    nx = nx_u - 1
    ny = ny_u

    DO k = 1, nz

      ! Linear interpolation along each row     
      DO i = 1, nx

        dataout(i,:,k) = 0.5*(datain(i,:,k) + datain(i+1,:,k))

      ENDDO
    ENDDO
    RETURN
  END SUBROUTINE arakawa_c_u2t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE wlevs2hlevs(datain,nx,ny,nz_in,dataout)

    ! Vertically destaggers an array from the full W levels to
    ! the half levels.  Output array is one less in the z dimension
    ! than the input array

    IMPLICIT NONE
    INTEGER, INTENT(IN)   :: nx,ny,nz_in
    REAL, INTENT(IN)      :: datain(nx,ny,nz_in)
    REAL, INTENT(OUT)     :: dataout(nx,ny,nz_in-1)

    INTEGER :: k

    DO k = 1,nz_in-1

      dataout(:,:,k) = 0.5*(datain(:,:,k)+datain(:,:,k+1))
  
    ENDDO
    RETURN
  END SUBROUTINE wlevs2hlevs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE grid_utils
