
subroutine kessler_mr2z(nx,ny,nz &
                       ,rho &                                ! kg/m^3
                       ,rainmr,icemr,snowmr,graupelmr  &
                       ,refl)

! Subroutine to compute estimated radar reflectivity (Z) from
! the precipitation mixing ratios.  The estimation
! is done using formulae from Kessler (1969) and 
! Rogers and Yau (1989).  

! Adapted from USAF Weather Agency routine.  
! Brent Shaw, NOAA Forecast System Lab, Dec 2000

! Put into subroutine by Steve Albers
  
implicit none

integer :: nx,ny,nz,i,j,k
real, parameter :: svnfrth=7.0/4.0
real, dimension(nx,ny,nz) :: rho,rainmr,icemr,snowmr,graupelmr,refl

refl=0.0

do j=1,ny
do i=1,nx
   do k=1,nz

!     Compute the basic reflectivity 
      refl(i,j,k) =17300.0 * &
                  (rho(i,j,k) * 1000.0 * &
                   MAX(0.0,rainmr(i,j,k)))**svnfrth

!     Add the ice component
      refl(i,j,k)=refl(i,j,k) + &
                  38000.0*(rho(i,j,k) * 1000.0 * &
                  MAX(0.0,icemr(i,j,k)+snowmr(i,j,k)+graupelmr(i,j,k)))**2.2

   enddo

enddo
enddo

return
end

