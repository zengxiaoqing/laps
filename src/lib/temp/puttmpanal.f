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

        subroutine put_temp_anal(i4time_needed,ni,nj,nk  ! Input
     1          ,ni_maps,nj_maps                 ! Input  (no longer needed)
     1          ,dum1_3d                         ! Dummy
     1          ,dum2_3d                         ! Dummy
     1          ,field_3d_maps_2                 ! Dummy  (no longer needed)
     1          ,local_nx_m,local_ny_m,iden_ratio,bias_3d ! Dummy
     1          ,array_density_box,n_density_box ! Dummy
     1          ,r0_array_in,r0_array_out        ! Dummy
     1          ,dum1_array,dum2_array           ! Dummy
     1          ,dum3_array,dum4_array           ! Dummy
     1          ,heights_3d                      ! Output
     1          ,output_4d                       ! Dummy
     1          ,lat,lon,topo                    ! Input
     1          ,temp_sfc_k                      ! Input
     1          ,pres_sfc_pa                     ! Input
     1          ,iflag_write                     ! Input
     1          ,ilaps_cycle_time                ! Input
     1          ,grid_spacing_m                  ! Input
     1          ,sh_3d                           ! Local
     1          ,temp_3d,istatus)                ! Output

        entry put_temp_anal_97(i4time_needed
     1          ,ni,nj,nk                        ! Input
     1          ,dum1_3d                         ! Dummy
     1          ,dum2_3d                         ! Dummy
     1          ,bias_3d                         ! Dummy
     1          ,r0_array_in,r0_array_out        ! Dummy
     1          ,dum1_array,dum2_array           ! Dummy
     1          ,dum3_array,dum4_array           ! Dummy
     1          ,heights_3d                      ! Output
     1          ,output_4d                       ! Dummy
     1          ,lat,lon,topo                    ! Input
     1          ,temp_sfc_k                      ! Input
     1          ,pres_sfc_pa                     ! Input
     1          ,iflag_write                     ! Input
     1          ,ilaps_cycle_time                ! Input
     1          ,grid_spacing_m                  ! Input
     1          ,sh_3d                           ! Local
     1          ,temp_3d,istatus)                ! Output

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

        character*50 DIRECTORY
        character*31 ext

        character*10  units_2d
        character*125 comment_2d
        character*3 var_2d

        real*4 temp_3d(ni,nj,nk) ! Output
        real*4 heights_3d(ni,nj,nk) ! Output
        real*4 output_4d(ni,nj,nk,2)

        real*4 sh_3d(ni,nj,nk)     ! Local
        real*4 temp_sfc_k(ni,nj)   ! Input
        real*4 pres_sfc_pa(ni,nj)  ! Input
        real*4 theta(nk)
        real*4 height(nk)

        real*4 dum1_3d(ni,nj,nk)          ! Local
        real*4 dum2_3d(ni,nj,nk)          ! Local

!       These Dummy arrays are passed in
        real*4 bias_3d(ni,nj,nk)
        real*4    r0_array_in(ni,nj)
        real*4    r0_array_out(ni,nj)

        real*4 dum1_array(ni,nj)
        real*4 dum2_array(ni,nj)
        real*4 dum3_array(ni,nj)
        real*4 dum4_array(ni,nj)

        real*4 lat(ni,nj),lon(ni,nj),topo(ni,nj)

        character*9 asc9_tim

        real*4 diff_tol,cold_thresh
        parameter (diff_tol = 25.)
        parameter (cold_thresh = 170.)

        logical l_fill,l_adjust_heights

        O_K(T_K,P_PA)   =   O( T_K-273.15 , P_PA/100. )  + 273.15
        TDA_K(T_K,P_PA) = TDA( T_K-273.15 , P_PA/100. )  + 273.15

        ISTAT = INIT_TIMER()

        write(6,*)' Welcome to subroutine put_temp_anal'

        i4time_raob_window = 0 ! 43200

        l_adjust_heights = .true.

        do k = 1,nk
            height(k) = height_of_level(k)
        enddo

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

!       Get RAMS/MODEL Data
c  following line changed from T to T3 for v3 readlapsdata LW 9/97
        var_2d = 'T3'
        l_fill = .true.

        call get_modelfg_3d(i4time_needed,var_2d,ni,nj,nk
     1                                   ,temp_3d,istatus)

        if(istatus .ne. 1)then
            write(6,*)
     1     ' Returning from PUT_TEMP_ANAL without writing LT1 (or equiva
     1lent) file'
            return
        endif

        write(6,*)' Getting RAMS/MODEL heights'

        var_2d = 'HT'
        l_fill = .true.

        call get_modelfg_3d(i4time_needed,var_2d,ni,nj,nk
     1                                   ,heights_3d,istatus)

        if(istatus .ne. 1)then
            write(6,*)
     1     ' Returning from PUT_TEMP_ANAL without writing LT1 (or equiva
     1lent) file'
            return
        endif

        write(6,*)' Getting RAMS/MODEL SH'

        var_2d = 'SH'
        l_fill = .true.

        call get_modelfg_3d(i4time_needed,var_2d,ni,nj,nk
     1                               ,sh_3d,istatus)

        if(istatus .ne. 1)then
            write(6,*)
     1     ' Returning from PUT_TEMP_ANAL without writing LT1 (or equiva
     1lent) file'
            return
        endif
!!
!!      added by pas, 28 jun 91...i4_filename not defined, so getting 0 for
!!      i4time and 600010000 for filename...not gt.
!!

        i4_filename = i4time_needed
        call make_fnam_lp(i4_filename,asc9_tim,istatus)
        comment_2d(1:9) = asc9_tim
        write(6,*) ' i4time_needed = ', i4time_needed
        write(6,*) ' i4_filename = ', i4_filename
        write(6,*) ' asc9_tim = ', asc9_tim

!       Test for bad temperatures
        if(temp_3d(29,29,13) .lt. 173.)then
            istatus = 0
            write(6,*)
     1   ' Bad data from interpolated Background Temps LGA/LGF/RAM, '
     1          ,'no LT1 (or equiv) written'
            return
        endif

        call insert_tsnd      (i4time_needed        ! Input
     1               ,lat,lon            ! Input
     1               ,heights_3d         ! Input
     1               ,sh_3d              ! Input
     1               ,dum1_3d            ! Used as wt_3d dummy
     1               ,dum2_3d            ! Used as bias_obs_3d dummy
     1               ,temp_3d            ! Input/Output
     1               ,bias_3d            ! Dummy
     1               ,r0_array_in,r0_array_out                 ! Dummy
     1               ,ilaps_cycle_time   ! Input
     1               ,i4time_raob_window ! Input
     1               ,ni,nj,nk           ! Input
     1               ,grid_spacing_m     ! Input
     1               ,istatus)           ! Output

        if(istatus .ne. 1)then
            write(6,*)' Bad status returned from insert_tsnd'
            return
        endif

!       Insert Surface Temp at Lowest Levels
        blayer_thk_pres = 5000. ! 7500.
        pres_intvl = pressure_of_level(1) - pressure_of_level(2)

!       This quantity limits the correction from LAPS surface data when we
!       are below the ground so things don't get too out of hand.
        frac_bias_max = (pres_intvl + blayer_thk_pres) / blayer_thk_pres

        write(6,*)' Inserting Surface Data in Lower Levels',blayer_thk_p
     1res
     1                                             ,frac_bias_max

        do i = 1,ni
        do j = 1,nj

!           Find Temperature at Top of Boundary Layer According to
!                                                             Upper Level Anal
            rk_sfc = zcoord_of_pressure(pres_sfc_pa(i,j))
            k_sfc = int(rk_sfc)

!           QC Check (Compare MODEL Temps to LAPS Sfc Temp)
!           For this check, interpolation in P space is sufficient
            k_sfc_qc = max(k_sfc,1)
            frac_k_sfc = rk_sfc - k_sfc_qc
            temp_sfc_intrpl = temp_3d(i,j,k_sfc_qc) * (1.0 - frac_k_sfc)       
     1                      + temp_3d(i,j,k_sfc_qc+1)  *     frac_k_sfc

!           Store the sfc temperature in a local variable
            temp_sfc_eff = temp_sfc_k(i,j)

            if(.false.)then ! Use sfc only if when terrain is below about 750mb
                frac_sfc = (pres_sfc_pa(i,j) - 72500.) / 5000.
                frac_sfc = min(max(frac_sfc,0.0),1.0)
                temp_sfc_eff = frac_sfc * temp_sfc_eff
     1                       + (1.0 - frac_sfc) * temp_sfc_intrpl
            endif

            diff_intrpl = temp_sfc_eff - temp_sfc_intrpl ! LAPS - MODEL
            if(diff_intrpl .gt. diff_tol)then
                write(6,111)i,j,temp_sfc_eff,temp_sfc_intrpl,diff_intrpl
     1                     ,k_sfc_qc,temp_3d(i,j,k_sfc_qc)
     1                              ,temp_3d(i,j,k_sfc_qc+1)
111             format('  LAPS Sfc Temps disagree with MDL',2i4,3f8.1
     1                ,i4,2f8.1)
                temp_sfc_eff = temp_sfc_intrpl + diff_tol
            elseif(diff_intrpl .lt. -25.)then
                write(6,111)i,j,temp_sfc_eff,temp_sfc_intrpl,diff_intrpl
     1                     ,k_sfc_qc,temp_3d(i,j,k_sfc_qc)
     1                              ,temp_3d(i,j,k_sfc_qc+1)
                temp_sfc_eff = temp_sfc_intrpl - diff_tol
            endif

            pres_top_pa = pres_sfc_pa(i,j) - blayer_thk_pres

            height_top = psatoz(pres_top_pa      * .01)
            height_sfc = psatoz(pres_sfc_pa(i,j) * .01)

            rk_top = zcoord_of_pressure(pres_top_pa)
            k_top = int(rk_top)

!           Fill in from level 1 (even if below the terrain) up through top of
!           boundary layer using the surface temperature and the first guess
!           (MODEL + RASS). This is phased in by calculating a residual of LAPS
!           sfc minus the first guess. The residual is ramped in by adding
!           half the bias in the middle of the boundary layer, the full bias
!           at the surface, etc. This helps preserve the vertical temperature
!           structure of the first guess within the boundary layer while
!           insuring consistency between the final LT1 temps and the LAPS sfc
!           temp analysis.

            rk_sfc = zcoord_of_pressure(pres_sfc_pa(i,j))
            k_sfc = int(rk_sfc)
            k_sfc = max(k_sfc,1)
            frac_k_sfc = rk_sfc - k_sfc
            temp_sfc_intrpl = temp_3d(i,j,k_sfc  ) * (1.0 - frac_k_sfc)
     1              + temp_3d(i,j,k_sfc+1) *        frac_k_sfc

            sfc_bias = temp_sfc_eff - temp_sfc_intrpl

            do k = 1,k_top

                frac_bias = (height_top - height(k)) / (height_top - hei
     1ght_sfc)
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

!           Insure that dtheta/dz > 0; Adjust Temps if Necessary to adiabatic

            if(temp_sfc_eff .lt. 200. .or. temp_sfc_eff .gt. 400.)then
                write(6,*)' Bad Sfc Temp',i,j,temp_sfc_eff
                istatus = 0
                return
            endif

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
     1                          .or. temp_3d(i,j,k) .gt. 400.)then
                    write(6,*)' Bad 3D/sfc Temp',i,j,k,temp_3d(i,j,k)
     1                                          ,temp_sfc_eff
!                   write(6,111)i,j,temp_sfc_eff,temp_sfc_intrpl,diff_intrpl
                    if(k .ge. k_sfc)then
                        istatus = 0
                        return
                    endif
                endif

                theta(k) = O_K(temp_3d(i,j,k),zcoord_of_level(k))

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
     1,f8.2)
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
                temp_3d(i,j,k) = TDA_K(theta(k),zcoord_of_level(k))
                temp_3d(i,j,k) = max(temp_3d(i,j,k),cold_thresh)
            enddo ! k

        enddo ! j
        enddo ! i

        write(6,*)' Temperature Analyses Read In and Adjusted'

        write(6,201)diff_thmax,theta_diff_thmax,i_diff_thmax,
     1                          j_diff_thmax,k_diff_thmax
201     format('  Maximum Adiabatic Adjustment of ',f8.1,' to ',f8.1,' a
     1t ',
     1                     i3,i4,i3)

        write(6,211)diff_min,t_diff_min,i_diff_min,
     1                          j_diff_min,k_diff_min
211     format('  Largest Cold Adjustment of ',5x,f8.1,' to ',f8.1,' at 
     1',
     1                     i3,i4,i3)

        write(6,221)diff_max,t_diff_max,i_diff_max,
     1                          j_diff_max,k_diff_max
221     format('  Largest Warm Adjustment of ',5x,f8.1,' to ',f8.1,' at 
     1',
     1                     i3,i4,i3)

!       Double Check 3D Temps against Sfc Temps
!       Here, interpolation in standard atmosphere height space is done
        diffmax = 0.
        do i = 1,ni
        do j = 1,nj
            rk_sfc = zcoord_of_pressure(pres_sfc_pa(i,j))
            k_sfc = int(rk_sfc)
            k_sfc_qc = max(k_sfc,1)
            frac_k_sfc = rk_sfc - k_sfc_qc
            temp_sfc_intrpl = temp_3d(i,j,k_sfc_qc  ) * (1.0 - frac_k_sf
     1c)
     1              + temp_3d(i,j,k_sfc_qc+1) *        frac_k_sfc

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

        write(6,*)' Max difference of sfc temps - interpolated 3D temps 
     1= '
     1  ,d_diff,i_diff,j_diff,rm_diff,t_diff,p_diff

        if(l_adjust_heights)then ! Store model fg 500 heights
            k_ref = zcoord_of_pressure(50000.)
            do i = 1,ni
            do j = 1,nj
                r0_array_in(i,j) = heights_3d(i,j,k_ref) ! Using up a dummy array
            enddo ! j
            enddo ! i
        endif

        write(6,*)' Calling get_heights_hydrostatic'
        call get_heights_hydrostatic(temp_3d,pres_sfc_pa,topo,
     1          dum1_array,dum2_array,dum3_array,dum4_array,
     1                                  ni,nj,nk,heights_3d)

        if(l_adjust_heights)then ! Adjust height field to model fg 500 heights
            call adjust_heights(temp_3d,heights_3d,r0_array_in
     1                         ,ni,nj,nk,k_ref,istatus)

        endif


        write(6,*)' Final temp/ht column at (ichk,jchk,k) with LAPS sfc 
     1t/p in'
     1                                  ,temp_sfc_k(ichk,jchk)
     1                                  ,pres_sfc_pa(ichk,jchk)
        do k = 1,nk
            write(6,1001)k,temp_3d(ichk,jchk,k),heights_3d(ichk,jchk,k)
1001        format(i4,f7.1,f7.0)
        enddo ! k

        I4_elapsed = ishow_timer()

        if(iflag_write .eq. 1)then
            call write_temp_anal(i4time_needed,ni,nj,nk,temp_3d
     1                  ,heights_3d,output_4d,comment_2d,istatus)
        endif

        I4_elapsed = ishow_timer()

        return
        end

        subroutine write_temp_anal(i4time,imax,jmax,kmax,temp_3d
     1                  ,heights_3d,output_4d,comment_2d,istatus)

        integer*4 nf
        parameter (nf = 2)

        character*31 EXT

        character*125 comment_a(nf),comment_2d
        character*10 units_a(nf)
        character*3 var_a(nf)

        real*4 temp_3d(imax,jmax,kmax)
        real*4 heights_3d(imax,jmax,kmax)
        real*4 output_4d(imax,jmax,kmax,2)

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
