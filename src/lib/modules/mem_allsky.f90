
MODULE mem_allsky

!     Input arrays on model grid
      real, allocatable, dimension(:,:,:) :: pres_3d
      real, allocatable, dimension(:,:,:) :: heights_3d
      real, allocatable, dimension(:,:,:) :: clwc_3d
      real, allocatable, dimension(:,:,:) :: cice_3d
      real, allocatable, dimension(:,:,:) :: rain_3d
      real, allocatable, dimension(:,:,:) :: snow_3d
      real, allocatable, dimension(:,:,:) :: aod_3d

!     Local arrays on model grid
      real, allocatable, dimension(:,:,:) :: transm_3d
      real, allocatable, dimension(:,:,:,:) :: transm_4d
      real, allocatable, dimension(:,:,:,:) :: uprad_4d 

!     2D arrays on sky grid
      real, allocatable, dimension(:,:) :: aod_ill_opac
      real, allocatable, dimension(:,:) :: aod_ill_opac_potl

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

        write(6,*)' allsky successfully deallocated'

        return
   
     END subroutine dealloc_allsky

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE mem_allsky
