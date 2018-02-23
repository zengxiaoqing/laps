
MODULE mem_allsky

!     Input arrays on model grid
      real, allocatable, dimension(:,:,:) :: pres_3d
      real, allocatable, dimension(:,:,:) :: heights_3d
      real, allocatable, dimension(:,:,:) :: clwc_3d    ! kg/m**3
      real, allocatable, dimension(:,:,:) :: cice_3d    ! kg/m**3
      real, allocatable, dimension(:,:,:) :: rain_3d    ! kg/m**3
      real, allocatable, dimension(:,:,:) :: snow_3d    ! kg/m**3
      real, allocatable, dimension(:,:,:) :: aod_3d

!     Local arrays on model grid
      real, allocatable, dimension(:,:,:) :: transm_3d
      real, allocatable, dimension(:,:,:,:) :: transm_4d
      real, allocatable, dimension(:,:,:,:) :: uprad_4d 
      real, allocatable, dimension(:,:,:) :: upxrad_3d 
      real, allocatable, dimension(:,:,:) :: upyrad_3d 

!     2D arrays on sky grid
      real, allocatable, dimension(:,:) :: aod_ill_opac
      real, allocatable, dimension(:,:) :: aod_ill_opac_potl

!     Various non-gridded variables
      real ghi_sim
      real alpha_ha
      real ext_g(3)
      integer mode_aero_cld /1/ ! treat aerosols more as clouds [1,2,3]

      PUBLIC alloc_allsky, dealloc_allsky

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     subroutine alloc_allsky(ni,nj,nk,nc,istatus)   ! I/O

!       Allocate some though not all arrays mentioned above

        allocate(pres_3d(ni,nj,nk))
        allocate(heights_3d(ni,nj,nk))
        allocate(clwc_3d(ni,nj,nk))
        allocate(cice_3d(ni,nj,nk))
        allocate(rain_3d(ni,nj,nk))
        allocate(snow_3d(ni,nj,nk))
        allocate(aod_3d(ni,nj,nk))
        allocate(transm_3d(ni,nj,nk))
        allocate(transm_4d(ni,nj,nk,nc))
        allocate(uprad_4d(ni,nj,nk,nc))
        allocate(upxrad_3d(ni,nj,nk))
        allocate(upyrad_3d(ni,nj,nk))

        write(6,*)' allsky successfully allocated'

        return
   
     END subroutine alloc_allsky

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     subroutine dealloc_allsky()             

!       Deallocate some though not all arrays mentioned above

        deallocate(pres_3d)
        deallocate(heights_3d)
        deallocate(clwc_3d)
        deallocate(cice_3d)
        deallocate(rain_3d)
        deallocate(snow_3d)
        deallocate(aod_3d)
        deallocate(transm_3d)
        deallocate(transm_4d)
        deallocate(uprad_4d)
        deallocate(upxrad_3d)
        deallocate(upyrad_3d)

        write(6,*)' allsky successfully deallocated'

        return
   
     END subroutine dealloc_allsky

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE mem_allsky
