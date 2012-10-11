Module mem_vol  

integer, allocatable, dimension(:,:,:) :: Reflectivity, Reflectivity_HI, RadialVelocity, RadialVelocity_HI

real, allocatable, dimension(:,:) :: elevationR, elevationR_HI, elevationV, elevationV_HI          

real, allocatable, dimension(:,:) :: azimuthR, azimuthR_HI, azimuthV, azimuthV_HI          

real, allocatable, dimension(:) :: distanceR, distanceR_HI, distanceV, distanceV_HI, nyquistVelocityV, nyquistVelocityV_HI          

integer :: gateR, gateR_HI, gateV, gateV_HI, radialR, radialR_HI, &
           radialV, radialV_HI, scanR, scanR_HI, scanV, &
           scanV_HI,nf_fid, nf_vid, nf_status  

character*8 :: c8_fname_format

end Module
