
        subroutine barnes_multivariate_sfc(to_2d_in,ni,nj        ! Inputs
     1                   ,r_missing_data                         ! Input
     1                   ,max_snd                                ! Input
     1                   ,rms_thresh_norm                        ! Input
     1                   ,rinst_err                              ! Input
     1                   ,weight_bkg_const                       ! Input
     1                   ,n_fnorm                                ! Input
     1                   ,l_boundary                             ! Input
     1                   ,n_obs,obs_barnes                       ! Input
     1                   ,t_2d                                   ! Output
     1                   ,istatus)                               ! Output

        integer*4 max_obs
        parameter (max_obs = 40000)       

        include 'barnesob.inc'
        type (barnesob_qc) obs_barnes(max_obs)
        type (barnesob)    obs_barnes_valid(max_obs)

        integer nk
        parameter (nk=1)

        real*4 t_2d(ni,nj)                              ! Analyzed field
        real*4 to_2d_in(ni,nj)                          ! Observations
        real*4 to_2d(ni,nj)                             ! Observations

        real*4 wt_2d(ni,nj)
        integer*4 n_obs_lvl

        logical l_analyze(nk), l_3d, l_use_ob, l_boundary(ni,nj)

        logical l_barnes_wide 
        logical l_struct_in, l_struct_out, l_not_struct_out
        data l_struct_in  /.false./    ! Using data structures?
        data l_struct_out /.true./     ! Using data structures?
        data l_barnes_wide /.true./    ! Using barnes_wide routine on boundary?

        integer*4  n_fnorm
        dimension fnorm(0:n_fnorm)

        write(6,*)' Subroutine Barnes_univariate_sfc'

        ISTAT = INIT_TIMER()

        call get_grid_spacing_cen(grid_spacing_cen_m,istatus)
        if(istatus .ne. 1)return

        l_analyze = .false.             ! array set
        wt_2d = 0.                      ! array set

        write(6,*)' setting observation to array weights'

        sumsq_inst = 0.
        n_obs_valid = 0
        ibound = 0
        iskip_bound = 10 ! How far to skip between boundary obs that get used.
                         ! 1 means use all boundary obs, 2 means use every
                         ! other boundary ob, 3 means use every 3rd boundary ob.

        boundary_err = rinst_err
            
        if(l_struct_in)then
          if(l_barnes_wide)then
            do iob = 1,n_obs
                if(obs_barnes(iob)%qc)then
                    i = obs_barnes(iob)%i
                    j = obs_barnes(iob)%j              
                    if(l_boundary(i,j))then
                        ibound = ibound + 1
                        if(ibound .ne. (ibound/iskip_bound)*iskip_bound
     1                                                            )then        
                            l_use_ob = .false.
                        endif
                    endif
                endif
            enddo ! iob
          endif ! l_barnes_wide

!         Place valid obs from obs_barnes into obs_barnes_valid
          do iob = 1,n_obs
            if(obs_barnes(iob)%qc)then
                n_obs_valid = n_obs_valid + 1
                obs_barnes_valid(n_obs_valid)%i = obs_barnes(iob)%i
                obs_barnes_valid(n_obs_valid)%j = obs_barnes(iob)%j
                obs_barnes_valid(n_obs_valid)%k = obs_barnes(iob)%k
                obs_barnes_valid(n_obs_valid)%weight 
     1                                        = obs_barnes(iob)%weight       
                obs_barnes_valid(n_obs_valid)%value(1) 
     1                                        = obs_barnes(iob)%value(1)       
            endif ! valid ob
          enddo ! iob

        else ! .not. l_struct_in
          do i = 1,ni
          do j = 1,nj
            l_use_ob = .true.

!           Determine if point would represent a boundary ob that gets skipped
            if(l_barnes_wide)then
                if(l_boundary(i,j))then
                    ibound = ibound + 1
                    if(ibound .ne. (ibound/iskip_bound)*iskip_bound)then       
                        l_use_ob = .false.
                    endif
                endif
            endif

            if(to_2d_in(i,j) .ne. 0.0 .and. l_use_ob)then
                wt_2d(i,j) = 1.0 / rinst_err**2 ! Set differently for boundary?
                sumsq_inst = sumsq_inst + 1./wt_2d(i,j)
                n_obs_valid = n_obs_valid + 1
!               to_2d(i,j) = to_2d_in(i,j)

!               Fill data structure element
                obs_barnes_valid(n_obs_valid)%i = i
                obs_barnes_valid(n_obs_valid)%j = j
                obs_barnes_valid(n_obs_valid)%k = 1
                obs_barnes_valid(n_obs_valid)%weight=1.0 / rinst_err**2       
                obs_barnes_valid(n_obs_valid)%value(1) = to_2d_in(i,j)

            else
!               to_2d(i,j) = r_missing_data
            endif

          enddo ! j
          enddo ! i
        endif

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

        l_not_struct_out = .not. l_struct_out
        n_var = 1
        rep_pres_intvl = 5000. ! Hardwire should work for a 2-D analysis

        call barnes_multivariate(
     1                      t_2d                                  ! Outputs
     1                     ,n_var,n_obs_valid,obs_barnes_valid    ! Input
     1                     ,ni,nj,nk,grid_spacing_cen_m           ! Inputs
     1                     ,rep_pres_intvl                        ! Input
     1                     ,to_2d,wt_2d,fnorm,n_fnorm             ! Inputs
     1                     ,l_analyze,l_not_struct_out,rms_thresh ! Input
     1                     ,weight_bkg_const                      ! Input
     1                     ,n_obs_lvl,istatus)                    ! Outputs

        return
        end

