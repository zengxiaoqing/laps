
        subroutine barnes_multivariate_sfc_jacket(c_field,obs,mxstn   ! I
     1                                           ,tb                  ! I
     1                                           ,badflag,ni,nj       ! I
     1                                           ,rms_thresh_norm     ! I
     1                                           ,bad_mult            ! I
     1                                           ,t_2d,istatus)       ! O

!       This subroutine can be a substitute to the call to 'spline' and 
!       interfaces more directly with 'barnes_multivariate_sfc'. It also
!       makes more complete use the data structures for storing the obs.

        character*(*) c_field

!       Input observation data structure (all variables)
        include 'sfcob.inc'
        type (sfcob) obs(mxstn)

!       Second observation data structure (just one variable)
        include 'barnesob.inc'
        type (barnesob_qc) obs_barnes(mxstn)
        type (barnesob)    obs_barnes_valid(mxstn)

        real*4 tb(ni,nj)                                ! Background field
        real*4 t_2d(ni,nj)                              ! Analyzed field
        real*4 ob_diff(mxstn)
        real*4 ob_full(mxstn)
        real*4 ob_bkg(mxstn)

        write(6,*)' Subroutine barnes_multivariate_sfc_jacket for...'
     1           ,c_field      

!       Transfer from 'obs' data structure to 'obs_barnes' data structure
        n_obs = 0
        if(c_field .eq. 'tgd')then
            rinst_err = 1.5
            do i = 1,mxstn
                n_obs = n_obs + 1

                ista = obs(i)%i
                jsta = obs(i)%j
                obs_barnes(i)%i = ista
                obs_barnes(i)%j = jsta

                obs_barnes(i)%k = 1
                obs_barnes(i)%weight = 1. / rinst_err**2

                if(obs(i)%sfct_f .ne. badflag)then
                    obs_barnes(i)%qc = .true.
                else
                    obs_barnes(i)%qc = .false.
                endif

!               Reject obs (for now) if they are outside of the domain
                if(ista .lt. 1 .or. ista .gt. ni 
     1        .or. jsta .lt. 1 .or. jsta .gt. nj )then
                    obs_barnes(i)%qc = .false.
                endif

!               Calculate difference of ob from the background
                if(obs_barnes(i)%qc)then
                    ob_full(i) = obs(i)%sfct_f
                    ob_bkg(i)  = tb(ista,jsta)
                    ob_diff(i) = ob_full(i) - ob_bkg(i)
                    obs_barnes(i)%value(1) = ob_diff(i)
                endif

            enddo ! i
        endif

!       Calculate rms (stdev) of the ob-background differences
        cnt = 0
        sumsq = 0.
        do i = 1,mxstn
            if(obs_barnes(i)%qc)then
                cnt = cnt + 1.
                sumsq = sumsq + ob_diff(i)**2
            endif
        enddo ! i

        if(cnt .gt. 0.)then
            std = sqrt(sumsq / cnt)
        else
            std = 0.
        endif

        bad = bad_mult * std
	print *,' std dev: ',std,', bad value: ',bad

!       Eliminate bad data from the data structure
        do i = 1,mxstn
            if(obs_barnes(i)%qc)then
                if(ob_diff(i) .ge. bad)then
                    obs_barnes(i)%qc = .false.
                    write(6,1099) i, j, ob_full(i), ob_diff(i)
1099	            format(1x,'bad data at i,j ',2i5,': value ',e12.4       
     1                    ,', diff ',e12.4)
                endif
            endif
        enddo ! i

        call get_fnorm_max(imax,jmax,r0_norm,r0_value_min,fnorm_max)   
        n_fnorm = int(fnorm_max) + 1

        weight_bkg_const = 5e28 ! a la wind.nl

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)return

        if(.true.)return ! Temporary

!       Note that we do not need to fill to_2d_dum when l_struct = .true.
!       Note that we don't need to fill l_boundary when l_struct = .true.
!                                                  and  l_barnes_wide = .false.

        call barnes_multivariate_sfc(to_2d_dum,ni,nj             ! Inputs
     1                   ,r_missing_data                         ! Input
     1                   ,rms_thresh_norm                        ! Input
     1                   ,rinst_err                              ! Input
     1                   ,weight_bkg_const                       ! Input
     1                   ,n_fnorm                                ! Input
     1                   ,l_boundary,.false.,.true.              ! Input
     1                   ,n_obs_valid,obs_barnes_valid           ! Input
     1                   ,t_2d                                   ! Output
     1                   ,istatus)                               ! Output

        return
        end

        subroutine barnes_multivariate_sfc(to_2d_in,ni,nj        ! Inputs
     1                   ,r_missing_data                         ! Input
     1                   ,rms_thresh_norm                        ! Input
     1                   ,rinst_err                              ! Input
     1                   ,weight_bkg_const                       ! Input
     1                   ,n_fnorm                                ! Input
     1                   ,l_boundary,l_barnes_wide,l_struct_in   ! Input
     1                   ,n_obs,obs_barnes                       ! Input
     1                   ,t_2d                                   ! Output
     1                   ,istatus)                               ! Output

        include 'barnesob.inc'
        type (barnesob_qc) obs_barnes(n_obs)
        type (barnesob)    obs_barnes_valid(n_obs)

        integer nk
        parameter (nk=1)

        real*4 t_2d(ni,nj)                              ! Analyzed field
        real*4 to_2d_in(ni,nj)                          ! Observations
        real*4 to_2d_dum(ni,nj)                         

        real*4 wt_2d(ni,nj)
        integer*4 n_obs_lvl

        logical l_analyze(nk), l_use_ob, l_boundary(ni,nj)

        logical l_barnes_wide          ! Use 'barnes_wide' routine on boundary?
        logical l_struct_in,l_struct_out,l_not_struct_out ! Use data structures?
        data l_struct_out /.true./     ! Use data structures?
!       data l_barnes_wide /.true./    

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

!               Fill data structure element
                obs_barnes_valid(n_obs_valid)%i = i
                obs_barnes_valid(n_obs_valid)%j = j
                obs_barnes_valid(n_obs_valid)%k = 1
                obs_barnes_valid(n_obs_valid)%weight=1.0 / rinst_err**2       
                obs_barnes_valid(n_obs_valid)%value(1) = to_2d_in(i,j)

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
     1                     ,to_2d_dum,wt_2d,fnorm,n_fnorm         ! Inputs
     1                     ,l_analyze,l_not_struct_out,rms_thresh ! Input
     1                     ,weight_bkg_const                      ! Input
     1                     ,n_obs_lvl,istatus)                    ! Outputs

        return
        end

