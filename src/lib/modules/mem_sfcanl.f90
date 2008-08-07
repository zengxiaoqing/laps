
Module mem_sfcanl



type sfcanl_fields
real, pointer, dimension(:,:)    ::  u,v,pr,t,td,vv,rh,pm,tad  &
                                    ,th,the,hi,ps,vor,qm,qcon,div,thad  &
                                    ,qad,spd,css,vis,fwx,tgd
end type

type(sfcanl_fields) :: sfcanl

! Pointers for renaming arrays in lapsvanl 
real, pointer, dimension(:,:) :: &
         u_a,v_a,p_a,t,td,vv,rh,hi,mslp,tadv,theta,thetae,psfc  &
        ,vort,q,qcon,div,thadv,qadv,spd,cssi,vis,fire,tgd_k


Contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine alloc_sfcanl_arrays (nxl, nyl)

implicit none
integer :: nxl, nyl

   allocate(sfcanl%u    (nxl,nyl))
   allocate(sfcanl%v    (nxl,nyl))
   allocate(sfcanl%pr   (nxl,nyl))
   allocate(sfcanl%t    (nxl,nyl))
   allocate(sfcanl%td   (nxl,nyl))
   allocate(sfcanl%vv   (nxl,nyl))
   allocate(sfcanl%rh   (nxl,nyl))
   allocate(sfcanl%pm   (nxl,nyl))
   allocate(sfcanl%tad  (nxl,nyl))
   allocate(sfcanl%th   (nxl,nyl))
   allocate(sfcanl%the  (nxl,nyl))
   allocate(sfcanl%hi   (nxl,nyl))
   allocate(sfcanl%ps   (nxl,nyl))
   allocate(sfcanl%vor  (nxl,nyl))
   allocate(sfcanl%qm   (nxl,nyl))
   allocate(sfcanl%qcon (nxl,nyl))
   allocate(sfcanl%div  (nxl,nyl))
   allocate(sfcanl%thad (nxl,nyl))
   allocate(sfcanl%qad  (nxl,nyl))
   allocate(sfcanl%spd  (nxl,nyl))
   allocate(sfcanl%css  (nxl,nyl))
   allocate(sfcanl%vis  (nxl,nyl))
   allocate(sfcanl%fwx  (nxl,nyl))
   allocate(sfcanl%tgd  (nxl,nyl))

  
return
end subroutine alloc_sfcanl_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine point_sfcanl_arrays ()

implicit none

u_a   =>  sfcanl%u   
v_a   =>  sfcanl%v   
p_a   =>  sfcanl%pr  
t     =>  sfcanl%t   
td    =>  sfcanl%td  
vv    =>  sfcanl%vv  
rh    =>  sfcanl%rh  
hi    =>  sfcanl%pm  
mslp  =>  sfcanl%tad 
tadv  =>  sfcanl%th  
theta =>  sfcanl%the 
thetae=>  sfcanl%hi  
psfc  =>  sfcanl%ps  
vort  =>  sfcanl%vor 
q     =>  sfcanl%qm  
qcon  =>  sfcanl%qcon
div   =>  sfcanl%div 
thadv =>  sfcanl%thad
qadv  =>  sfcanl%qad 
spd   =>  sfcanl%spd 
cssi  =>  sfcanl%css 
vis   =>  sfcanl%vis 
fire  =>  sfcanl%fwx 
tgd_k =>  sfcanl%tgd 

  
return
end subroutine point_sfcanl_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine deallocate_sfcanl_arrays()

implicit none

   if (associated(sfcanl%u    )) deallocate(sfcanl%u    )
   if (associated(sfcanl%v    )) deallocate(sfcanl%v    )
   if (associated(sfcanl%pr   )) deallocate(sfcanl%pr   )
   if (associated(sfcanl%t    )) deallocate(sfcanl%t    )
   if (associated(sfcanl%td   )) deallocate(sfcanl%td   )
   if (associated(sfcanl%vv   )) deallocate(sfcanl%vv   )
   if (associated(sfcanl%rh   )) deallocate(sfcanl%rh   )
   if (associated(sfcanl%pm   )) deallocate(sfcanl%pm   )
   if (associated(sfcanl%tad  )) deallocate(sfcanl%tad  )
   if (associated(sfcanl%th   )) deallocate(sfcanl%th   )
   if (associated(sfcanl%the  )) deallocate(sfcanl%the  )
   if (associated(sfcanl%hi   )) deallocate(sfcanl%hi   )
   if (associated(sfcanl%ps   )) deallocate(sfcanl%ps   )
   if (associated(sfcanl%vor  )) deallocate(sfcanl%vor  )
   if (associated(sfcanl%qm   )) deallocate(sfcanl%qm   )
   if (associated(sfcanl%qcon )) deallocate(sfcanl%qcon )
   if (associated(sfcanl%div  )) deallocate(sfcanl%div  )
   if (associated(sfcanl%thad )) deallocate(sfcanl%thad )
   if (associated(sfcanl%qad  )) deallocate(sfcanl%qad  )
   if (associated(sfcanl%spd  )) deallocate(sfcanl%spd  )
   if (associated(sfcanl%css  )) deallocate(sfcanl%css  )
   if (associated(sfcanl%vis  )) deallocate(sfcanl%vis  )
   if (associated(sfcanl%fwx  )) deallocate(sfcanl%fwx  )
   if (associated(sfcanl%tgd  )) deallocate(sfcanl%tgd  )

return
end subroutine deallocate_sfcanl_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine nullify_sfcanl_arrays ()

implicit none

   if (associated(sfcanl%u    )) nullify(sfcanl%u    )
   if (associated(sfcanl%v    )) nullify(sfcanl%v    )
   if (associated(sfcanl%pr   )) nullify(sfcanl%pr   )
   if (associated(sfcanl%t    )) nullify(sfcanl%t    )
   if (associated(sfcanl%td   )) nullify(sfcanl%td   )
   if (associated(sfcanl%vv   )) nullify(sfcanl%vv   )
   if (associated(sfcanl%rh   )) nullify(sfcanl%rh   )
   if (associated(sfcanl%pm   )) nullify(sfcanl%pm   )
   if (associated(sfcanl%tad  )) nullify(sfcanl%tad  )
   if (associated(sfcanl%th   )) nullify(sfcanl%th   )
   if (associated(sfcanl%the  )) nullify(sfcanl%the  )
   if (associated(sfcanl%hi   )) nullify(sfcanl%hi   )
   if (associated(sfcanl%ps   )) nullify(sfcanl%ps   )
   if (associated(sfcanl%vor  )) nullify(sfcanl%vor  )
   if (associated(sfcanl%qm   )) nullify(sfcanl%qm   )
   if (associated(sfcanl%qcon )) nullify(sfcanl%qcon )
   if (associated(sfcanl%div  )) nullify(sfcanl%div  )
   if (associated(sfcanl%thad )) nullify(sfcanl%thad )
   if (associated(sfcanl%qad  )) nullify(sfcanl%qad  )
   if (associated(sfcanl%spd  )) nullify(sfcanl%spd  )
   if (associated(sfcanl%css  )) nullify(sfcanl%css  )
   if (associated(sfcanl%vis  )) nullify(sfcanl%vis  )
   if (associated(sfcanl%fwx  )) nullify(sfcanl%fwx  )
   if (associated(sfcanl%tgd  )) nullify(sfcanl%tgd  )

return
end subroutine nullify_sfcanl_arrays
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end Module
