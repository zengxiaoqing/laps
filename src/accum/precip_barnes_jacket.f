
        subroutine precip_barnes_jacket(          c_field             ! I
     1                                           ,ilaps,jlaps         ! I
     1                                           ,pcp_gauge           ! I
     1                                           ,maxsta              ! I
     1                                           ,pcp_bkg_in          ! I
     1                                           ,badflag,ni,nj       ! I
     1                                           ,topo,ldf            ! I
     1                                           ,wt_bkg_a            ! I
     1                                           ,pcp_2d_in,istatus)  ! O

        include 'constants.inc'
        character*(*) c_field

        include 'barnesob.inc' ! Observation data structure to be passed onward
        type (barnesob_qc) obs_barnes(maxsta)

        integer ilaps(maxsta),jlaps(maxsta)

        real pcp_bkg_in(ni,nj)                        ! Background field
        real pcp_2d_in(ni,nj)                         ! Analyzed field (IN)
        real topo(ni,nj),ldf(ni,nj)                   ! Topo & Landfrac
        real pcp_gauge(maxsta)
        real ob_bkg(maxsta)
        real ob_full(maxsta)
        real ob_diff(maxsta)
        real wt_bkg_a(ni,nj)                         

        logical l_boundary(ni,nj)

        write(6,*)' Subroutine precip_barnes_jacket for...'
     1           ,c_field      

        call get_systime_i4(i4time_sys,istatus) ! temporary while no time wt

        n_obs = 0
        do i = 1,maxsta
            n_obs = n_obs + 1

            ista = ilaps(i)
            jsta = jlaps(i)
            obs_barnes(i)%i = ista
            obs_barnes(i)%j = jsta

            ob_full(i) = pcp_gauge(i) ! inches
            if(pcp_gauge(i) .ge. 0.)then
                obs_barnes(i)%qc = .true.
            else
                obs_barnes(i)%qc = .false.
            endif

            obs_barnes(i)%k = 1
            obs_barnes(i)%weight = 1. ! / rinst_err**2
            obs_barnes(i)%i4time = i4time_sys ! no time weight for now

!           Reject obs (for now) if they are outside of the domain
            if(ista .lt. 1 .or. ista .gt. ni 
     1    .or. jsta .lt. 1 .or. jsta .gt. nj )then
                obs_barnes(i)%qc = .false.
            endif

!           Calculate difference of ob from the background
            if(obs_barnes(i)%qc)then
                ob_bkg(i)  = pcp_bkg_in(ista,jsta)
                ob_diff(i) = ob_full(i) - ob_bkg(i)
                obs_barnes(i)%value(1) = ob_diff(i)
            endif

        enddo ! i

!       Calculate rms (stdev) of the ob-background differences
        cnt = 0
        sumsq = 0.
        do i = 1,maxsta
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


        call get_fnorm_max(ni,nj,r0_norm,r0_value_min,fnorm_max)   
        n_fnorm = int(fnorm_max) + 1

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)return

!       if(.true.)return ! Temporary

!       Note that we don't need to fill l_boundary when l_struct = .true.
!                                                  and  l_barnes_wide = .false.

        call precip_barnes_sfc(ni,nj                             ! Inputs
     1                   ,r_missing_data                         ! Input
     1                   ,rinst_err                              ! Input
     1                   ,bad                                    ! Input
     1                   ,wt_bkg_a                               ! Input
     1                   ,n_fnorm                                ! Input
     1                   ,l_boundary,.false.,.true.              ! Input
     1                   ,maxsta,obs_barnes                      ! Input
     1                   ,topo,ldf                               ! Input
     1                   ,pcp_2d_in                              ! Output
     1                   ,istatus)                               ! Output

        write(6,*)' Adding incremental analysis to background'       

        call array_range(pcp_2d_in,ni,nj,pmin,pmax,r_missing_data)
        write(6,*)' pcp_2d_in range',pmin,pmax

        call array_range(pcp_bkg_in,ni,nj,pmin,pmax,r_missing_data)
        write(6,*)' pcp_bkg_in range',pmin,pmax
        
        call add(pcp_2d_in,pcp_bkg_in,pcp_2d_in,ni,nj)

        write(6,*)' Returning from precip_barnes_jacket'
        write(6,*)

        return
        end

        subroutine precip_barnes_sfc(ni,nj                       ! Inputs
     1                   ,r_missing_data                         ! Input
     1                   ,rinst_err                              ! Input
     1                   ,bad                                    ! Input
     1                   ,wt_bkg_a                               ! Input
     1                   ,n_fnorm                                ! Input
     1                   ,l_boundary,l_barnes_wide,l_struct_in   ! Input
     1                   ,n_obs,obs_barnes                       ! Input
     1                   ,topo,ldf                               ! Input
     1                   ,pcp_2d_in                              ! Output
     1                   ,istatus)                               ! Output

        include 'barnesob.inc'
        type (barnesob_qc) obs_barnes(n_obs)
        type (barnesob)    obs_barnes_valid(n_obs)

        integer nk
        parameter (nk=1)

        real pcp_2d_in(ni,nj)                         ! Analyzed field
        real wt_bkg_a(ni,nj)                         
        real topo(ni,nj),ldf(ni,nj)                   ! Topo & Landfrac

        real wt_2d(ni,nj)
        integer n_obs_lvl

        logical l_analyze(nk), l_use_ob, l_boundary(ni,nj)

        logical l_barnes_wide          ! Use 'barnes_wide' routine on boundary?

        logical l_struct_in,l_struct_out,l_not_struct_out ! Use data structures?
        data l_struct_out /.true./     ! Use data structures?

        integer  n_fnorm
        dimension fnorm(0:n_fnorm)

        write(6,*)' Subroutine precip_barnes_sfc'

        I4_elapsed = ishow_timer()

        call get_grid_spacing_cen(grid_spacing_cen_m,istatus)
        if(istatus .ne. 1)return

        call get_systime_i4(i4time_sys,istatus) ! temporary while no time wt

        l_analyze = .false.             ! array set
        wt_2d = 0.                      ! array set

        write(6,*)' setting observation to array weights'

        sumsq_inst = 0.
        n_obs_valid = 0
        ibound = 0
        iskip_bound = 10 ! How far to skip between boundary obs that get used.
                         ! 1 means use all boundary obs, 2 means use every
                         ! other boundary ob, 3 means use every 3rd boundary ob.

        if(l_struct_in)then

!         Place valid obs from obs_barnes into obs_barnes_valid
          do iob = 1,n_obs
            if(obs_barnes(iob)%qc)then
                n_obs_valid = n_obs_valid + 1
                obs_barnes_valid(n_obs_valid)%i = obs_barnes(iob)%i
                obs_barnes_valid(n_obs_valid)%j = obs_barnes(iob)%j
                obs_barnes_valid(n_obs_valid)%k = obs_barnes(iob)%k
                obs_barnes_valid(n_obs_valid)%weight 
     1                                        = obs_barnes(iob)%weight
                obs_barnes_valid(n_obs_valid)%vert_rad_rat = 1.
                obs_barnes_valid(n_obs_valid)%i4time
     1                                        = obs_barnes(iob)%i4time
                obs_barnes_valid(n_obs_valid)%value(1) 
     1                                        = obs_barnes(iob)%value(1)
                sumsq_inst = sumsq_inst 
     1                     + 1. / obs_barnes_valid(n_obs_valid)%weight

                obs_barnes_valid(n_obs_valid)%elev 
     1                                        = obs_barnes(iob)%elev
                obs_barnes_valid(n_obs_valid)%ldf = obs_barnes(iob)%ldf
C
C TH: 29 November 2002 Begin Hack
C
                obs_barnes_valid(n_obs_valid)%mask_sea =
     1             obs_barnes(iob)%mask_sea
C
C TH: End hack.
C

            endif ! valid ob
          enddo ! iob

        endif

        if(n_obs_valid .gt. 0)then
            rms_inst = sqrt(sumsq_inst/float(n_obs_valid))
        else
            rms_inst = 0.
        endif

        write(6,*)' n_obs_valid,rms_inst = ',n_obs_valid,rms_inst       

!       Set the rms threshold for iteration cutoff, based on instrument error
        rms_thresh = 0.01 ! rms_inst (inches)

        write(6,*)'rms_thresh',rms_thresh      

        l_not_struct_out = .not. l_struct_out
        n_var = 1
        rep_pres_intvl = 5000. ! Hardwire should work for a 2-D analysis

        r0_barnes_max_m = 140000.
        barnes_conv_rate = 0.5

!       Perform gauge analysis of incremental precip in inches
        call barnes_multivariate(
     1                      pcp_2d_in                             ! Outputs
     1                     ,n_var,n_obs_valid,obs_barnes_valid    ! Input
     1                     ,ni,nj,nk,grid_spacing_cen_m           ! Inputs
     1                     ,rep_pres_intvl                        ! Input
     1                     ,wt_bkg_a                              ! Input
     1                     ,i4_loop_total                         ! Output
     1                     ,wt_2d                                 ! Not used
     1                     ,fnorm,n_fnorm                         ! Inputs
     1                     ,l_analyze,l_not_struct_out,rms_thresh ! Input
     1                     ,r0_barnes_max_m,barnes_conv_rate      ! Input
     1                     ,topo,ldf,ni,nj                        ! Input
     1                     ,n_obs_lvl,istatus)                    ! Outputs

        call stats(pcp_2d_in,ni,nj)

        return
        end


