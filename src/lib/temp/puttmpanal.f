cdis    Forecast Systems Laboratory
cdis    NOAA/OAR/ERL/FSL
cdis    325 Broadway
cdis    Boulder, CO     80303
cdis 
cdis    Forecast Research Division
cdis    Local Analysis and Prediction Branch
cdis    LAPS 
cdis 
cdis    This software and its documentation are in the public domain and 
cdis    are furnished "as is."  The United States government, its 
cdis    instrumentalities, officers, employees, and agents make no 
cdis    warranty, express or implied, as to the usefulness of the software 
cdis    and documentation for any purpose.  They assume no responsibility 
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis    
cdis    Permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  All modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making 
cdis    the modifications.  If significant modifications or enhancements 
cdis    are made to this software, the FSL Software Policy Manager  
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 

        subroutine put_temp_anal(i4time_needed
     1          ,ni,nj,nk                        ! Input
     1          ,heights_3d                      ! Output
     1          ,lat,lon,topo                    ! Input
     1          ,temp_sfc_k                      ! Input
     1          ,pres_sfc_pa                     ! Input
     1          ,pres_msl_pa                     ! Input
     1          ,ilaps_cycle_time                ! Input
     1          ,grid_spacing_m                  ! Input
     1          ,comment_2d                      ! Output
     1          ,temp_3d,pres_3d,istatus)        ! Output

!              1991     Steve Albers    Original Version
!          Oct 1991     Steve Albers    Add Sfc pres as input to more accurately
!                                       place sfc temp analysis in 3D domain
!          Apr 1992     Steve Albers    Reduced QC check from 200K to 180K
!       17 Dec 1992     Steve Albers    Add call for RASS data
!       10 Jun 1993     Steve Albers    Use RAMS for background if available
!          Oct 1993     Steve Albers    Add arrays for multiple RASSes
!          Oct 1993     Steve Albers    Add arrays for outputting height field
!          Dec 1993     Steve Albers    Add call to get_heights_hydrostatic
!          Aug 1994     Steve Albers    Pass in rh_3d
!        9 Sep 1994     Steve Albers    Remove explicit RAMS T call
!          Sep 1994     Steve Albers    Convert from RH to SH model input
!          Dec 1994     S.A.            Option to use sfc data only below 750mb
!          1995 Dec 19  Steve Albers    Added timer calls
!          1997 Jun 16  Ken Dritz       Changed NZ_L_MAX to nk, making theta
!                                       and height arrays automatic.

!       In this version, sfc pressure is used to correctly place the surface
!       in the pressure coordinate grid. However, DIFFERENCES in heights for
!       dealing with lapse rates are calculated using the standard atmosphere.
!       This should be an acceptable approximation for this application.

        use mem_namelist, ONLY: l_read_raob=>l_read_raob_t 
     1                         ,l_use_raob=>l_use_raob_t
     1                         ,mode_adjust_heights
     1                         ,weight_bkg_const=>weight_bkg_const_temp
     1                         ,pres_mix_thresh
     1                         ,rms_thresh=>rms_thresh_temp
     1                         ,max_snd_grid,max_snd_levels,max_obs

        character*125 comment_2d
        character*3 var_2d

        real temp_3d(ni,nj,nk) ! Output
        real heights_3d(ni,nj,nk) ! Output

        real sh_3d(ni,nj,nk)     ! Local
        real pres_3d(ni,nj,nk)   ! Output
        real temp_sfc_k(ni,nj)   ! Input
        real pres_sfc_pa(ni,nj)  ! Input
        real pres_msl_pa(ni,nj)  ! Input
        real theta(nk)

        real bkg_500(ni,nj)

        real lat(ni,nj),lon(ni,nj),topo(ni,nj)

        real ht_ref(ni,nj),pres_ref(ni,nj)

        character*9 asc9_tim

        real diff_tol,cold_thresh
        parameter (diff_tol = 25.)
        parameter (cold_thresh = 170.)

        logical l_fill

        O_K(T_K,P_PA)   =   O( T_K-273.15 , P_PA/100. )  + 273.15
        TDA_K(T_K,P_PA) = TDA( T_K-273.15 , P_PA/100. )  + 273.15

        ISTAT = INIT_TIMER()

        write(6,*)' Welcome to subroutine put_temp_anal'

        icen = ni/2
        jcen = nj/2

!       Initialize diagnostic variables
        diff_thmax = 0.
        i_diff_thmax = 0
        j_diff_thmax = 0
        k_diff_thmax = 0
        theta_diff_thmax = 0.

        diff_max = 0.
        i_diff_max = 0
        j_diff_max = 0
        k_diff_max = 0
        t_diff_max = 0.

        diff_min = 0.
        i_diff_min = 0
        j_diff_min = 0
        k_diff_min = 0
        t_diff_min = 0.

!       Get BACKGROUND MODEL Data

        write(6,*)' Getting BACKGROUND MODEL temperatures'

        var_2d = 'T3'
        l_fill = .true.

        call get_modelfg_3d(i4time_needed,var_2d,ni,nj,nk
     1                                   ,temp_3d,istatus)

        if(istatus .ne. 1)then
            write(6,*)' ERROR: No Model First Guess for ',var_2d
            write(6,*)
     1     ' Returning from PUT_TEMP_ANAL without writing LT1 file'
            return
        endif

        write(6,*)' Getting BACKGROUND MODEL heights'

        var_2d = 'HT'
        l_fill = .true.

        call get_modelfg_3d(i4time_needed,var_2d,ni,nj,nk
     1                                   ,heights_3d,istatus)

        if(istatus .ne. 1)then
            write(6,*)' ERROR: No Model First Guess for ',var_2d
            write(6,*)
     1     ' Returning from PUT_TEMP_ANAL without writing LT1 file'
            return
        endif

!       QC check the background model heights against the topo field
        do i = 1,ni
        do j = 1,nj
            if(heights_3d(i,j,1) .gt. topo(i,j))then
                write(6,*)' WARNING: QC check failed, lowest background'
     1                   ,' height level extends above terrain field '  
                write(6,*)'i,j,height,topo'
     1                    ,i,j,heights_3d(i,j,1),topo(i,j)
                go to 150
            endif
        enddo ! j
        enddo ! i
 150    continue

        write(6,*)' Getting BACKGROUND MODEL SH'

        var_2d = 'SH'
        l_fill = .true.

        call get_modelfg_3d(i4time_needed,var_2d,ni,nj,nk
     1                               ,sh_3d,istatus)

        if(istatus .ne. 1)then
            write(6,*)' ERROR: No Model First Guess for ',var_2d
            write(6,*)
     1     ' Returning from PUT_TEMP_ANAL without writing LT1 file'
            return
        endif
!!
!!      added by pas, 28 jun 91...i4_filename not defined, so getting 0 for
!!      i4time and 600010000 for filename...not gt.
!!

        i4_filename = i4time_needed
        call make_fnam_lp(i4_filename,asc9_tim,istatus)

        write(6,*) ' i4time_needed = ', i4time_needed
        write(6,*) ' i4_filename = ', i4_filename
        write(6,*) ' asc9_tim = ', asc9_tim

        call get_pres_3d(i4time_needed,ni,nj,nk,pres_3d,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error: Bad status returned from get_pres_3d'       
            return
        endif

        call insert_tobs(i4time_needed                   ! Input
     1                  ,lat,lon                         ! Input
     1                  ,heights_3d                      ! Input
     1                  ,sh_3d                           ! Input
     1                  ,pres_3d                         ! Input
     1                  ,temp_3d                         ! Input/Output
     1                  ,ilaps_cycle_time                ! Input
     1                  ,l_read_raob                     ! Input
     1                  ,l_use_raob                      ! Input
     1                  ,weight_bkg_const                ! Input
     1                  ,rms_thresh                      ! Input
     1                  ,ni,nj,nk                        ! Input
     1                  ,max_snd_grid,max_snd_levels     ! Input
     1                  ,max_obs                         ! Input
     1                  ,grid_spacing_m                  ! Input
     1                  ,num_temp_obs                    ! Output
     1                  ,istatus)                        ! Output
        if(istatus .ne. 1)then
            write(6,*)' Error: Bad status returned from insert_tobs'       
            return
        endif

        write(comment_2d,1)asc9_tim,num_temp_obs
 1      format(a9,' num_temp_obs = ',i9)

!       Height Analysis (done prior to surface temp analysis insertion)
        write(6,*)
        if(mode_adjust_heights .eq. 1)then ! Store model fg 500 heights
            arg = rlevel_of_field(50000.
     1                           ,pres_3d,ni,nj,nk,icen,jcen,istatus)
            if(istatus .ne. 1)return
            k_ref = nint(arg)

            write(6,*)' Storing bkg ht level ',k_ref
            do i = 1,ni
            do j = 1,nj
                bkg_500(i,j) = heights_3d(i,j,k_ref) 
                if(bkg_500(i,j) .lt. 3000. .or. 
     1             bkg_500(i,j) .gt. 6500.       )then
                    write(6,*)' Error: bkg ht failed QC check at level'       
     1                       ,' nearest 500mb ',i,j,k_ref,bkg_500(i,j)
                    istatus = 0
                    return
                endif ! QC check
            enddo ! j
            enddo ! i
        endif

!       Hydrostatically integrate the temps to get heights using the surface
!       pressure as a reference. The heights_3d array will now contain the
!       integrated heights instead of the model background heights.

        if(mode_adjust_heights .eq. 0 .OR. 
     1     mode_adjust_heights .eq. 1     )then 
            pres_ref = pres_sfc_pa
            ht_ref = topo
        elseif(mode_adjust_heights .eq. 2)then ! use MSLP as a reference
            write(6,*)' Hydrostatic integration - MSLP reference'
            pres_ref = pres_msl_pa
            ht_ref = 0.
        else
            write(6,*)' ERROR: mode_adjust_heights = '
     1                        ,mode_adjust_heights
            istatus = 0
            return
        endif

        write(6,*)' Calling get_heights_hydrostatic'
        call get_heights_hydrostatic(temp_3d,pres_ref,pres_3d,sh_3d
     1                              ,ht_ref,ni,nj,nk,heights_3d,istatus)       

        if(mode_adjust_heights .eq. 1)then 

            write(6,*)' Hydrostatic integration - 500mb ht reference'

!           Adjust height field to model fg 500 heights
            call adjust_heights(temp_3d,heights_3d,bkg_500
     1                         ,ni,nj,nk,k_ref,istatus)

        endif

        I4_elapsed = ishow_timer()

!       Insert Surface Temp at Lowest Levels
        blayer_thk_pres = 2500. 

        write(6,*)' Inserting Surface Data in Lower Levels'
     1             ,blayer_thk_pres

!       Set up stuff for mixed layer top
        iwarn = 0
        diff_temp_max = 0.
        i_max_diff = 0
        j_max_diff = 0
        sum_pres = 0.

        do i = 1,ni
        do j = 1,nj
            sum_pres = sum_pres + pres_sfc_pa(i,j)
        enddo ! j
        enddo ! i

        pres_ave_domain = sum_pres / float(ni*nj)
        pres_mix_top = pres_ave_domain - pres_mix_thresh
        write(6,*)' Average sfc pressure over the domain = '
     1            ,pres_ave_domain
        write(6,*)' pres_mix_thresh = ',pres_mix_thresh
        write(6,*)' pres_mix_top = ',pres_mix_top

        do i = 1,ni
        do j = 1,nj
            if(pres_mix_top .lt. pres_3d(i,j,nk))then ! QC check
                pres_mix_top = pres_3d(i,j,nk)
            endif

            rk_mix = rlevel_of_field(pres_mix_top
     1                              ,pres_3d,ni,nj,nk,i,j,istatus)
            if(istatus .ne. 1)return

            k_mix  = int(rk_mix)
            frac_k_mix = rk_mix - k_mix

!           Find Temp at Top of Boundary Lyr According to Upper Level Anal
            rk_sfc = rlevel_of_field(pres_sfc_pa(i,j)
     1                              ,pres_3d,ni,nj,nk,i,j,istatus)
            if(istatus .ne. 1)return

            k_sfc = int(rk_sfc)
            if(k_sfc .lt. 1 .or. k_sfc .ge. nk)then
                write(6,*)' Error, k_sfc/rk_sfc/p_sfc = '
     1                    ,k_sfc,rk_sfc,pres_sfc_pa(i,j)
                istatus = 0
                return
            endif

            pres_intvl = pres_3d(i,j,k_sfc) - pres_3d(i,j,k_sfc+1)

!           This quantity limits the correction from LAPS surface data when we
!           are below the ground so things don't get too out of hand.
            frac_bias_max = (pres_intvl + blayer_thk_pres) 
     1                     / blayer_thk_pres
            frac_bias_max = 1.0 ! Potential mod for AFWA concerns

!           QC Check (Compare MODEL Temps to LAPS Sfc Temp)
!           For this check, interpolation in P space is sufficient
            frac_k_sfc = rk_sfc - k_sfc
            temp_sfc_intrpl = temp_3d(i,j,k_sfc) * (1.0 - frac_k_sfc)       
     1                      + temp_3d(i,j,k_sfc+1)  *     frac_k_sfc

!           Store the sfc temperature that is finally used in a local variable
            temp_sfc_eff = temp_sfc_k(i,j)

            diff_intrpl = temp_sfc_eff - temp_sfc_intrpl ! LAPS - MODEL
            if(diff_intrpl .gt. diff_tol)then
                write(6,111)i,j,temp_sfc_eff,temp_sfc_intrpl,diff_intrpl
     1                     ,k_sfc,temp_3d(i,j,k_sfc)
     1                           ,temp_3d(i,j,k_sfc+1)
111             format('  LAPS Sfc Temps disagree with MDL',2i4,3f8.1
     1                ,i4,2f8.1)
                temp_sfc_eff = temp_sfc_intrpl + diff_tol

            elseif(diff_intrpl .lt. -25.)then
                write(6,111)i,j,temp_sfc_eff,temp_sfc_intrpl,diff_intrpl
     1                     ,k_sfc,temp_3d(i,j,k_sfc)
     1                           ,temp_3d(i,j,k_sfc+1)
                temp_sfc_eff = temp_sfc_intrpl - diff_tol

            endif

!           Apply "theta check" upper bound to 'temp_sfc_eff'
            temp_mix_intrpl = temp_3d(i,j,k_mix) * (1.0 - frac_k_mix)       
     1                      + temp_3d(i,j,k_mix+1)  *     frac_k_mix

            theta_max_sfc = O_K(temp_mix_intrpl,pres_mix_top)            
            theta_now_sfc = O_K(temp_sfc_eff,pres_sfc_pa(i,j))
            temp_sfc_now = temp_sfc_eff

            if(pres_sfc_pa(i,j) .lt. pres_mix_top)then    ! sfc above mix top

!               Constrain the sfc temp to <= 3-d interpolated temp
                if(temp_sfc_now .gt. temp_sfc_intrpl)then 
                    temp_sfc_eff = temp_sfc_intrpl
                    diff_temp = temp_sfc_now - temp_sfc_eff
                    if(diff_temp .gt. diff_temp_max)then
                        diff_temp_max = diff_temp
                        i_max_diff = i
                        j_max_diff = j
                    endif
                    if(diff_temp .gt. 2.0)iwarn = 1
                endif                        

            else                                         ! sfc below mix top

!               Constrain the sfc theta to <= theta at top of the mixing layer
                if(theta_now_sfc .gt. theta_max_sfc)then 
                    temp_sfc_eff = TDA_K(theta_max_sfc,pres_sfc_pa(i,j))
                    diff_temp = temp_sfc_now - temp_sfc_eff
                    if(diff_temp .gt. diff_temp_max)then
                        diff_temp_max = diff_temp
                        i_max_diff = i
                        j_max_diff = j
                    endif
                    if(diff_temp .gt. 2.0)iwarn = 1
                endif                        

            endif                                        ! sfc above mix top?

!           Fill in from level 1 (even if below the terrain) up through top of
!           boundary layer using the surface temperature and the first guess
!           (MODEL + RASS). This is phased in by calculating a residual of LAPS
!           sfc minus the first guess. The residual is ramped in by adding
!           half the bias in the middle of the boundary layer, the full bias
!           at the surface, etc. This helps preserve the vertical temperature
!           structure of the first guess within the boundary layer while
!           insuring consistency between the final LT1 temps and the LAPS sfc
!           temp analysis.

            pres_top_pa = pres_sfc_pa(i,j) - blayer_thk_pres

            height_top = psatoz(pres_top_pa      * .01)   ! Standard Atmosphere
            height_sfc = psatoz(pres_sfc_pa(i,j) * .01)   ! Standard Atmosphere

            rk_top = rlevel_of_field(pres_top_pa
     1                              ,pres_3d,ni,nj,nk,i,j,istatus)
            if(istatus .ne. 1)return

            k_top = int(rk_top)

            sfc_bias = temp_sfc_eff - temp_sfc_intrpl

            do k = 1,k_top
                height_k = psatoz(pres_3d(i,j,k) * .01)   ! Standard Atmosphere

!               A Standard Atmosphere calculation is probably sufficient
                frac_bias = (height_top - height_k) /    
     1                      (height_top - height_sfc)

!               Potential mod for AFWA concerns - apply reverse ramp underground
                if(frac_bias .gt. frac_bias_max)then
                    frac_bias = (2. * frac_bias_max) - frac_bias
                    frac_bias = max(frac_bias,0.)
                endif

                frac_bias = min(frac_bias, frac_bias_max) ! Prevent overshooting below sfc

                temp_ref = temp_3d(i,j,k)
                temp_3d(i,j,k) = temp_3d(i,j,k) + sfc_bias * frac_bias

                if(k .gt. k_sfc)then
                    if(temp_3d(i,j,k) - temp_ref .lt. diff_min)then
                        diff_min = temp_3d(i,j,k) - temp_ref
                        i_diff_min = i
                        j_diff_min = j
                        k_diff_min = k
                        t_diff_min = temp_3d(i,j,k)
c                       write(6,211)diff_min,t_diff_min,i_diff_min,
c       1                               j_diff_min,k_diff_min
                    endif

                    if(temp_3d(i,j,k) - temp_ref .gt. diff_max)then
                        diff_max = temp_3d(i,j,k) - temp_ref
                        i_diff_max = i
                        j_diff_max = j
                        k_diff_max = k
                        t_diff_max = temp_3d(i,j,k)
c                       write(6,221)diff_max,t_diff_max,i_diff_max,
c       1                               j_diff_max,k_diff_max
                    endif

                endif
            enddo ! k

!           QC check of temp_sfc_eff - note this might be different from LSX T
            if(temp_sfc_eff .lt. 200. .or. temp_sfc_eff .gt. 400.)then
                write(6,*)' Error: Bad Sfc or MDL Temp'
     1                    ,i,j,temp_sfc_eff,temp_sfc_k(i,j)
     1                    ,temp_sfc_intrpl
                istatus = 0
                return
            endif

!           Ensure that dtheta/dz > 0; Adjust Temps if Necessary to adiabatic
            theta_sfc_k = O_K(temp_sfc_eff,pres_sfc_pa(i,j))

!           ichk = 27
!           jchk = 40

            ichk = 17
            jchk = 12

            if(i .eq. ichk .and. j .eq. jchk)then
                write(6,*)' t_sfc/p_sfc/theta_sfc_k',temp_sfc_eff
     1                             ,pres_sfc_pa(i,j),theta_sfc_k
            endif


!           QC the 3D temps and calculate theta in the column
            do k = 1,nk
                if(temp_3d(i,j,k) .lt. cold_thresh
     1                            .or. temp_3d(i,j,k) .gt. 400.)then
                    write(6,*)' Error: Bad 3D/sfc Temp',i,j,k
     1                        ,temp_3d(i,j,k),temp_sfc_eff
!                   write(6,111)i,j,temp_sfc_eff,temp_sfc_intrpl,diff_intrpl
                    if(k .ge. k_sfc)then
                        write(6,*)' k >= k_sfc --> return'
                        istatus = 0
                        return
                    endif
                endif

                theta(k) = O_K(temp_3d(i,j,k),pres_3d(i,j,k))

            enddo ! k

!           Adiabatically Adjust Downward from the Surface
!           The layer between each successive pair of grid points is examined
!           to see if it is superadiabatic.
!           Note that an earlier adjustment produced a uniform lapse rate
!           from level 1 through the top of the boundary layer.
            k = k_sfc
            theta_ref = theta_sfc_k
            do while (k .ge. 1)
                if(theta(k) .gt. theta_ref)then ! Superadiabatic
c                   if(i .eq. 1)
c       1               write(6,101)i,j,k,theta(k),theta_ref,rk_sfc
101                 format(' Adiabatic Lapse Rt Adj',3i3,1x,f8.1,4x,f8.1
     1                    ,f8.2)
                    theta(k) = theta_ref
                endif

                theta_ref = theta(k) ! Reset the reference theta
                k = k-1
            enddo ! while k >= 1

!           Adiabatically Adjust Upward from the Surface
!           The layer between each successive pair of grid points is examined
!           to see if it is superadiabatic.
!           Note that an earlier adjustment produced a uniform lapse rate
!           from level 1 through the top of the boundary layer.

            k = k_sfc + 1
            theta_ref = theta_sfc_k
            do while (k .le. nk)
                if(theta(k) .lt. theta_ref)then ! Superadiabatic
                    if(i .eq. ichk .and. j .eq. jchk)
     1          write(6,101)i,j,k,theta(k),theta_ref,rk_sfc

                    if(theta_ref - theta(k) .gt. diff_thmax)then
                        diff_thmax = theta_ref - theta(k)
                        i_diff_thmax = i
                        j_diff_thmax = j
                        k_diff_thmax = k
                        theta_diff_thmax = theta_ref
c                       write(6,201)diff_thmax,theta_diff_thmax,i_diff_thmax,
c       1                               j_diff_thmax,k_diff_thmax
                    endif

                    theta(k) = theta_ref
                endif

                theta_ref = theta(k) ! Reset the reference theta
                k = k+1
            enddo ! while k <= nk

            do k = 1,nk
                temp_3d(i,j,k) = TDA_K(theta(k),pres_3d(i,j,k))
                temp_3d(i,j,k) = max(temp_3d(i,j,k),cold_thresh)
            enddo ! k

        enddo ! j
        enddo ! i

        write(6,*)
        write(6,*)' Temperature Analysis Complete'

        I4_elapsed = ishow_timer()

        write(6,*)' Checking consistency of sfc and 3-D temps'

        if(iwarn .eq. 1)then
            write(6,*)' WARNING: sfc temps had to be reduced to be'
     1               ,' within adiabatic/3d tolerances,'
     1               ,' max reduction = ',diff_temp_max
     1               ,' at ',i_max_diff,j_max_diff
        else
            write(6,*)' Input sfc temps were adiabatically/3D'
     1               ,' consistent with other data within tolerances,'
     1               ,' max reduction = ',diff_temp_max
     1               ,' at ',i_max_diff,j_max_diff
        endif

!       Double Check 3D Temps against Sfc Temps
        write(6,*)' Double checking consistency of sfc and 3-D temps'
        diffmax = 0.
        do i = 1,ni
        do j = 1,nj
            rk_sfc = rlevel_of_field(pres_sfc_pa(i,j)
     1                              ,pres_3d,ni,nj,nk,i,j,istatus)
            if(istatus .ne. 1)return

            k_sfc = int(rk_sfc)
            frac_k_sfc = rk_sfc - k_sfc
            temp_sfc_intrpl = 
     1                temp_3d(i,j,k_sfc  ) * (1.0 - frac_k_sfc)
     1              + temp_3d(i,j,k_sfc+1) *        frac_k_sfc

            diff = abs(temp_sfc_k(i,j) - temp_sfc_intrpl)

            if(diff .gt. diffmax)then
                diffmax = diff
                d_diff = temp_sfc_k(i,j) - temp_sfc_intrpl
                i_diff = i
                j_diff = j
                rm_diff = temp_sfc_intrpl
                t_diff  = temp_sfc_k(i,j)
                p_diff = pres_sfc_pa(i,j)
            endif

        enddo ! j
        enddo ! i

!       This should be fairly close to diff_temp_max 
!                                      (blayer/adiabatic sfc T reduction)
        write(6,*)' Max diff of input surface T - interpolated 3D T = '       
     1           ,d_diff,i_diff,j_diff,t_diff,rm_diff,p_diff

        write(6,*)' Checking adjustments to 3-D temps from adiabatic'
     1           ,' and surface considerations'

        write(6,201)diff_thmax,theta_diff_thmax,i_diff_thmax,
     1                          j_diff_thmax,k_diff_thmax
201     format('  Maximum Adiabatic Adjustment of Theta',f8.1
     1        ,' to ',f8.1,' at ',i3,i4,i3)

        write(6,211)diff_min,t_diff_min,i_diff_min,
     1                       j_diff_min,k_diff_min
211     format('  Largest Cold Adjustment of Temp      ',f8.1
     1        ,' to ',f8.1,' at ',i3,i4,i3)

        write(6,221)diff_max,t_diff_max,i_diff_max,
     1                       j_diff_max,k_diff_max
221     format('  Largest Warm Adjustment of Temp      ',f8.1
     1        ,' to ',f8.1,' at ',i3,i4,i3)

!       Former location of Height Analysis

!       Write column of output
        write(6,*)
     1    ' Final temp/ht column at (ichk,jchk,k) with LAPS sfc t/p in'       
     1                                  ,temp_sfc_k(ichk,jchk)
     1                                  ,pres_sfc_pa(ichk,jchk)
        do k = 1,nk
            write(6,1001)k,temp_3d(ichk,jchk,k),heights_3d(ichk,jchk,k)
1001        format(i4,f7.1,f7.0)
        enddo ! k

        I4_elapsed = ishow_timer()

!       QC check the heights against the topo field
        do i = 1,ni
        do j = 1,nj
            if(heights_3d(i,j,1) .gt. topo(i,j))then
                write(6,*)' ERROR: QC check failed, lowest height'
     1                   ,' level extends above terrain field '  
                write(6,*)'i,j,height,topo'
     1                    ,i,j,heights_3d(i,j,1),topo(i,j)
                write(6,*)
     1          ' Returning from PUT_TEMP_ANAL without writing LT1 file'       
                istatus = 0
                return
            endif
        enddo ! j
        enddo ! i
        write(6,*)' Passed QC check of heights against the topo field'

        return
        end

        subroutine write_temp_anal(i4time,imax,jmax,kmax,temp_3d
     1                  ,heights_3d,comment_2d,istatus)

        integer nf
        parameter (nf = 2)

        character*31 EXT

        character*125 comment_a(nf),comment_2d
        character*10 units_a(nf)
        character*3 var_a(nf)

        real temp_3d(imax,jmax,kmax)
        real heights_3d(imax,jmax,kmax)
        real output_4d(imax,jmax,kmax,2) ! Local

        EXT = 'lt1'

        write(6,*)
     1 ' Writing out Temp/Height Analyses to LT1 (or equiv) Directory'

        var_a(1) = 'T3' ! newvar = 'T3', oldvar = 'T'
        var_a(2) = 'HT'

        units_a(1) = 'K'
        units_a(2) = 'M'

        do i = 1,nf
            comment_a(i) = comment_2d
        enddo ! i

        call move_3d(temp_3d(1,1,1)   ,output_4d(1,1,1,1)
     1                                                ,imax,jmax,kmax)
        call move_3d(heights_3d(1,1,1),output_4d(1,1,1,2)
     1                                                ,imax,jmax,kmax)      

        call put_laps_multi_3d(i4time,EXT,var_a,units_a,
     1          comment_a,output_4d,imax,jmax,kmax,nf,istatus)

        if(istatus .ne. 1)then
            write(6,*)' Error in put_laps_multi_3d for LT1 file'
        endif

        return
        end
! Hongli Jiang added to deal with sigma_ht. 11/2/2011
        subroutine write_temp_anal_ht(i4time,imax,jmax,kmax,temp_3d
     1                  ,pressure_3d,comment_2d,istatus)

        integer nf
        parameter (nf = 2)

        character*31 EXT

        character*125 comment_a(nf),comment_2d
        character*10 units_a(nf)
        character*3 var_a(nf)

        real temp_3d(imax,jmax,kmax)
        real pressure_3d(imax,jmax,kmax)
        real output_4d(imax,jmax,kmax,2) ! Local

        EXT = 'lt1'

        write(6,*)
     1 ' Writing out Temp/Height Analyses to LT1 (or equiv) Directory'

        var_a(1) = 'T3' ! newvar = 'T3', oldvar = 'T'
        var_a(2) = 'P3'

        units_a(1) = 'K'
        units_a(2) = 'Pa'

        do i = 1,nf
            comment_a(i) = comment_2d
        enddo ! i

        call move_3d(temp_3d(1,1,1)   ,output_4d(1,1,1,1)
     1                                                ,imax,jmax,kmax)
        call move_3d(pressure_3d(1,1,1),output_4d(1,1,1,2)
     1                                                ,imax,jmax,kmax)      

        call put_laps_multi_3d(i4time,EXT,var_a,units_a,
     1          comment_a,output_4d,imax,jmax,kmax,nf,istatus)

        if(istatus .ne. 1)then
            write(6,*)' Error in put_laps_multi_3d for LT1 file'
        endif

        return
        end

