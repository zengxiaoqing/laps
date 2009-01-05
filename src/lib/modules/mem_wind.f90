
Module mem_wind

type wind_fields
real, pointer, dimension(:,:,:)  :: uanl, vanl, wanl
real, pointer, dimension(:,:)    :: uanl_sfcitrp, vanl_sfcitrp
end type

type(wind_fields) :: wind

integer num_wind_obs

Contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine alloc_wind_arrays (nxl, nyl, nzl)

implicit none
integer :: nxl, nyl, nzl

integer :: nt

   allocate(wind%uanl (nxl,nyl,nzl))
   allocate(wind%vanl (nxl,nyl,nzl))
   allocate(wind%wanl (nxl,nyl,nzl))
   allocate(wind%uanl_sfcitrp (nxl,nyl))
   allocate(wind%vanl_sfcitrp (nxl,nyl))

  
return
end subroutine alloc_wind_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine deallocate_wind_arrays()

implicit none


integer :: nt

   if (associated(wind%uanl)) deallocate(wind%uanl )
   if (associated(wind%vanl)) deallocate(wind%vanl )
   if (associated(wind%wanl)) deallocate(wind%wanl )
   if (associated(wind%uanl_sfcitrp)) deallocate(wind%uanl_sfcitrp )
   if (associated(wind%vanl_sfcitrp)) deallocate(wind%vanl_sfcitrp )

return
end subroutine deallocate_wind_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine nullify_wind_arrays ()

implicit none

integer :: nt

   if (associated(wind%uanl)) nullify(wind%uanl )
   if (associated(wind%vanl)) nullify(wind%vanl )
   if (associated(wind%wanl)) nullify(wind%wanl )
   if (associated(wind%uanl_sfcitrp)) nullify(wind%uanl_sfcitrp )
   if (associated(wind%vanl_sfcitrp)) nullify(wind%vanl_sfcitrp )

return
end subroutine nullify_wind_arrays
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end Module
