
MODULE cloud_rad    

!     Cloud Radiation and Microphysics Parameters

!     Backscattering fractions
!     http://www.ugr.es/~aquiran/ciencia/arti27.pdf (aerosols)
      real, parameter :: bksct_eff_clwc    = .063
      real, parameter :: bksct_eff_cice    = .14
      real, parameter :: bksct_eff_rain    = .063
      real, parameter :: bksct_eff_snow    = .14
      real, parameter :: bksct_eff_graupel = .30
      real, parameter :: bksct_eff_aero    = .125 

!     Scattering efficiencies
      real, parameter :: q_clwc    = 2.0
      real, parameter :: q_cice    = 2.0
      real, parameter :: q_rain    = 1.0
      real, parameter :: q_snow    = 1.0
      real, parameter :: q_graupel = 1.0

!     Densities
      real, parameter :: rholiq     =   1e3 ! kilograms per cubic meter
      real, parameter :: rhosnow    = .07e3 ! kilograms per cubic meter
      real, parameter :: rhograupel = .50e3 ! kilograms per cubic meter

!     Effective radii
      real, parameter :: reff_clwc    = .000007 ! m
      real, parameter :: reff_cice    = .000034 ! m
      real, parameter :: reff_rain    = .000750 ! m
      real, parameter :: reff_snow    = .004000 ! m
      real, parameter :: reff_graupel = .010000 ! m

!     GHI related
      real, parameter :: ghi_zen_toa = 1361.5 ! solar const, W/m**2 at 1AU
      real, parameter :: zen_kt = 0.815       ! zenithal attenuation of GHI

      PUBLIC albedo_to_clouds, albedo_to_clouds2

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     subroutine albedo_to_clouds(albedo                                  & ! I
                                ,cloud_rad_trans,cloud_od,cloud_opacity)   ! O

        use mem_namelist, ONLY: r_missing_data

        clear_albedo = .2097063
        cloud_albedo = .4485300
!       cloud_albedo = .40

        if(albedo .eq. r_missing_data)then
            cloud_rad_trans = r_missing_data
            cloud_od = r_missing_data
            cloud_opacity = r_missing_data
            return
        endif

!       Uniform clouds are presently assumed in the grid box or cloud layer
        cloud_frac = 1.0

        arg = albedo

        call stretch2(clear_albedo,cloud_albedo,0.,1.,arg)

        cloud_albedo = min(arg,0.99999)      ! Back Scattered

        cloud_rad_trans = 1.0 - cloud_albedo ! Fwd Scattered + Direct Transmission

        bksc_eff_od = -log(cloud_rad_trans)  ! Tau * Back Scat Efficiency

        cloud_od = bksc_eff_od / 0.10        ! Tau

        cloud_opacity = 1.0 - exp(-cloud_od) ! 1 - Direct Transmission

        return
   
     END subroutine albedo_to_clouds    

     subroutine albedo_to_clouds2(albedo                                 & ! I
                                 ,cloud_trans_l,cloud_trans_i            & ! O
                                 ,cloud_od_l,cloud_od_i                  & ! O
                                 ,cloud_opac_l,cloud_opac_i)               ! O

        use mem_namelist, ONLY: r_missing_data

        if(albedo .eq. r_missing_data)then
            cloud_trans_l = r_missing_data
            cloud_trans_i = r_missing_data
            cloud_od_l = r_missing_data
            cloud_od_i = r_missing_data
            cloud_opac_l = r_missing_data
            cloud_opac_i = r_missing_data
            return
        endif

        cloud_albedo = min(albedo,0.99999)     ! Back Scattered

        cloud_trans_l = 1.0 - cloud_albedo     ! Fwd Scattered + Direct Transmission
        cloud_trans_i = 1.0 - cloud_albedo     ! Fwd Scattered + Direct Transmission

        cloud_od_l = cloud_albedo / (bksct_eff_clwc * (1.-cloud_albedo)) ! Tau
        cloud_od_i = cloud_albedo / (bksct_eff_cice * (1.-cloud_albedo)) ! Tau

        cloud_opac_l = 1.0 - exp(-cloud_od_l) ! 1 - Direct Transmission
        cloud_opac_i = 1.0 - exp(-cloud_od_i) ! 1 - Direct Transmission

        return
   
     END subroutine albedo_to_clouds2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE cloud_rad   
