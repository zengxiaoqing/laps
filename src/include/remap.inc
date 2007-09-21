c       File: remap.inc

c
c       Velocity Obs
c
        real grid_rvel(NX_L,NY_L,NZ_L)  !  Radial radar velocities
        real grid_rvel_sq(NX_L,NY_L,NZ_L)
        real grid_nyq(NX_L,NY_L,NZ_L)
        integer ngrids_vel(NX_L,NY_L,NZ_L)
        integer n_pot_vel(NX_L,NY_L,NZ_L)
c
c       Reflectivity Obs
c
        real grid_ref (NX_L,NY_L,NZ_L)  !  Radar reflectivities
        integer ngrids_ref (NX_L,NY_L,NZ_L)
        integer n_pot_ref (NX_L,NY_L,NZ_L)

c
!       save
!    :    grid_rvel,grid_rvel_sq,grid_nyq,ngrids_vel,n_pot_vel,
!    :    grid_ref,ngrids_ref,n_pot_ref,


!       real out_array_4d(NX_L,NY_L,NZ_L,3)
!       equivalence (out_array_4d(1,1,1,1),grid_ref(1,1,1))
!       equivalence (out_array_4d(1,1,1,2),grid_rvel(1,1,1))
!       equivalence (out_array_4d(1,1,1,3),grid_nyq(1,1,1))
