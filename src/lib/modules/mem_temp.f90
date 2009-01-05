!
! Copyright (C) 2007  ; All Rights Reserved ; ATMET, LLC
! 
! This file is free software; you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software 
! Foundation; either version 2 of the License, or (at your option) any later version.
! 
! This software is distributed in the hope that it will be useful, but WITHOUT ANY 
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
! PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with this 
! program; if not, write to the Free Software Foundation, Inc., 
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!======================================================================================

Module mem_temp



type temp_fields
real, pointer, dimension(:,:,:)  :: t3, ht, p3
end type

type(temp_fields) :: temp

! Pointers for renaming arrays
real, pointer, dimension(:,:,:) :: &
         temp_3d,heights_3d,pres_3d_pa

integer num_temp_obs

Contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine alloc_temp_arrays (nxl, nyl, nzl)

implicit none
integer :: nxl, nyl, nzl

integer :: nt

   allocate(temp%t3 (nxl,nyl,nzl))
   allocate(temp%ht (nxl,nyl,nzl))
   allocate(temp%p3 (nxl,nyl,nzl))

  
return
end subroutine alloc_temp_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine point_temp_arrays ()

implicit none

temp_3d =>  temp%t3   
heights_3d =>  temp%ht  
pres_3d_pa  =>  temp%p3  
  
return
end subroutine point_temp_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine deallocate_temp_arrays()

implicit none


integer :: nt

   if (associated(temp%t3)) deallocate(temp%t3 )
   if (associated(temp%ht)) deallocate(temp%ht )
   if (associated(temp%p3)) deallocate(temp%p3 )

return
end subroutine deallocate_temp_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine nullify_temp_arrays ()

implicit none

integer :: nt

   if (associated(temp%t3)) nullify(temp%t3 )
   if (associated(temp%ht)) nullify(temp%ht )
   if (associated(temp%p3)) nullify(temp%p3 )

return
end subroutine nullify_temp_arrays
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end Module
