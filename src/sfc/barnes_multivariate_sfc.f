
        subroutine barnes_multivariate_sfc(to_2d_in,ni,nj        ! Inputs
     1                   ,r_missing_data                         ! Input
     1                   ,max_snd                                ! Input
     1                   ,rms_thresh_norm                        ! Input
     1                   ,rinst_err                              ! Input
     1                   ,weight_bkg_const                       ! Input
     1                   ,n_fnorm                                ! Input
     1                   ,t_2d                                   ! Output
     1                   ,istatus)                               ! Output


        integer nk
        parameter (nk=1)

        real*4 t_2d(ni,nj)                              ! Analyzed field
        real*4 to_2d_in(ni,nj)                          ! Observations
        real*4 to_2d(ni,nj)                             ! Observations

        logical l_struct

        real*4 wt_2d(ni,nj)
        integer*4 n_obs_lvl

        logical l_analyze(nk),l_3d

        integer*4  n_fnorm

        data l_struct /.true./

        dimension fnorm(0:n_fnorm)

!       include 'tempobs.inc'

        write(6,*)' Subroutine Barnes_univariate_sfc'

        ISTAT = INIT_TIMER()

        call get_grid_spacing_cen(grid_spacing_cen_m,istatus)
        if(istatus .ne. 1)stop

        l_analyze = .false.             ! array set
        wt_2d = 0.                      ! array set

        write(6,*)' setting observation to array weights'

        sumsq_inst = 0.
        n_obs_valid = 0
            
        do i = 1,ni
        do j = 1,nj
            if(to_2d_in(i,j) .ne. 0.0)then
                wt_2d(i,j) = 1.0 / rinst_err**2
                sumsq_inst = sumsq_inst + 1./wt_2d(i,j)
                n_obs_valid = n_obs_valid + 1
                to_2d(i,j) = to_2d_in(i,j)
            else
                to_2d(i,j) = r_missing_data
            endif
        enddo ! j
        enddo ! i

        if(n_obs_valid .gt. 0)then
            rms_inst = sqrt(sumsq_inst/float(n_obs_valid))
        else
            rms_inst = 0.
        endif

        write(6,*)' n_obs_valid,rms_inst = ',n_obs_valid,rms_inst       

!       Set the rms threshold for iteration cutoff, based on instrument error
        rms_thresh = rms_thresh_norm * rms_inst

        write(6,*)'rms_thresh_norm,rms_thresh'
     1            ,rms_thresh_norm,rms_thresh      

        l_3d = .true.          ! Use successive correction even though it's 2D
        n_var = 1
        rep_pres_intvl = 5000. ! Hardwire should work for a 2-D analysis

        call barnes_multivariate(
     1                      t_2d                            ! Outputs
     1                     ,n_var                           ! Input
     1                     ,ni,nj,nk,grid_spacing_cen_m     ! Inputs
     1                     ,rep_pres_intvl                  ! Input
     1                     ,to_2d,wt_2d,fnorm,n_fnorm       ! Inputs
     1                     ,l_analyze,l_3d,rms_thresh       ! Input
     1                     ,weight_bkg_const                ! Input
     1                     ,n_obs_lvl,istatus)              ! Outputs

        return
        end

