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

        subroutine get_vis(i4time,solar_alt,l_use_vis              ! Input
     1                    ,i4_sat_window,i4_sat_window_offset      ! Input
     1                    ,cloud_frac_vis_a,albedo,ihist_alb       ! Output
     1                    ,ni,nj,nk,r_missing_data                 ! Input
     1                    ,istatus)                                ! Output

!       Steve Albers 1997

        integer*4 ihist_alb(-10:20)
        integer*4 ihist_frac_sat(-10:20)
        real*4 albedo(ni,nj)

!       This stuff is for reading VIS data from LVD file
        real*4 solar_alt(ni,nj)
        real*4 cloud_frac_vis_a(ni,nj)
        integer*4 mxstn
        parameter (mxstn = 100)       ! max number of "stations" in data file
        character*9 filename
        character*31 ext

        character*3 lvd_ext
        data lvd_ext /'lvd'/

        character var*3,comment*125,units*10
        logical l_use_vis

!       Initialize histograms
        do i = -10,20
            ihist_alb(i) = 0
            ihist_frac_sat(i) = 0
        enddo ! i

!       Initialize
        do i = 1,ni
        do j = 1,nj
            albedo(i,j) = r_missing_data
            cloud_frac_vis_a(i,j) = r_missing_data
        enddo ! j
        enddo ! i

        n_missing_albedo = ni*nj

!       Determine whether to use VIS / ALBEDO data
        if(.not. l_use_vis)then
            write(6,*)' Warning: l_use_vis set to not use vis data'
            istatus = 0
            return
        endif

!       Read in albedo data
        write(6,*)' Getting the VIS data from LVD file'
        call make_fnam_lp(i4time,filename,istatus)

        ext = lvd_ext
        var = 'ALB'
        ilevel = 0
        call get_laps_2dvar(i4time+i4_sat_window_offset,i4_sat_window
     1                     ,i4time_nearest,EXT,var,units
     1                     ,comment,ni,nj,albedo,ilevel,istatus)
        write(6,*)' istatus from albedo data = ',istatus
        if(istatus .ne. 1 .and. istatus .ne. -1)then
            write(6,*)' Warning: could not read albedo - '
     1               ,'return from get_vis'
            istatus = 0
            return
        endif

        n_missing_albedo = 0

!       Horizontal array loop
        do i = 1,ni
        do j = 1,nj

!         We will now only use the VIS data if the solar alt exceeds 15 deg
!         11 degrees now used to allow 15 min slack in data availability
          if(solar_alt(i,j) .lt. 11.0)then
              if(albedo(i,j) .ne. r_missing_data)then
                  write(6,*)' Error -  albedo not missing:',solar_alt(i,
     1j)
                  stop
              endif
!             albedo(i,j) = r_missing_data
          endif

          if(albedo(i,j) .ne. r_missing_data)then

!           Translate the albedo into cloud fraction

!           Store histogram information for satellite data
            iscr_alb  = nint(albedo(i,j)*10.)
            iscr_alb  = min(max(iscr_alb,-10),20)
            ihist_alb(iscr_alb) = ihist_alb(iscr_alb) + 1

            cloud_frac_vis = albedo_to_cloudfrac(albedo(i,j))

!           Fudge the frac at low solar elevation angles
!           Note the ramp extrapolates down to 9 deg to account for slight
!           errors in determining the solar elevation
            if(.false.)then
!           if(solar_alt(i,j) .lt. 20. .and. solar_alt(i,j) .ge. 9.)then
                frac = (20. - solar_alt(i,j)) / 10.
                term1 = .13 * frac
                term2 = 1. + term1
                cloud_frac_vis = (cloud_frac_vis + term1) * term2
            endif

            iscr_frac_sat = nint(cloud_frac_vis*10.)
            iscr_frac_sat = min(max(iscr_frac_sat,-10),20)
            ihist_frac_sat(iscr_frac_sat) = 
     1      ihist_frac_sat(iscr_frac_sat) + 1

!           Make sure satellite cloud fraction is between 0 and 1
            if(cloud_frac_vis .le. 0.0)cloud_frac_vis = 0.0
            if(cloud_frac_vis .ge. 1.0)cloud_frac_vis = 1.0

            cloud_frac_vis_a(i,j) = cloud_frac_vis

          else
            n_missing_albedo =  n_missing_albedo + 1
            cloud_frac_vis_a(i,j) = r_missing_data

          endif

        enddo ! i
        enddo ! j

        write(6,*)
        write(6,*)' N_MISSING_ALBEDO = ',n_missing_albedo
        write(6,*)

        if(n_missing_albedo .eq. ni*nj)then ! Return with status = 0
            write(6,*)' All albedos were missing - return from get_vis'
            istatus = 0
            return
        endif

        write(6,*)'              HISTOGRAMS'
        write(6,*)' I          ',
     1  ' Albedo  Cld Frac Sat'
        do i = -5,15
            write(6,11)i,ihist_alb(i),ihist_frac_sat(i)
11          format(i4,i12,i12,i12,i12)
        enddo ! i

        write(6,*)
        istatus = 1

        return
        end

        subroutine insert_vis(i4time,clouds_3d,cld_hts
     1      ,topo,cloud_frac_vis_a,albedo,ihist_alb                   ! I
     1      ,istat_39_a                                               ! I
     1      ,ni,nj,nk,r_missing_data                                  ! I
     1      ,vis_radar_thresh_cvr,vis_radar_thresh_dbz                ! I
     1      ,istat_radar,radar_ref_3d,klaps,ref_base
     1      ,dbz_max_2d,surface_sao_buffer,istatus)

!       Steve Albers 1999 Added 3.9u cloud clearing step in NULL Mode
!                         This can be turned on if desired or another cloud
!                         mask input(s) can be used for cloud clearing,
!                         especially if it applies throughout the column.

        integer*4 ihist_alb(-10:20)
        integer*4 ihist_frac_sat(-10:20)
        integer*4 ihist_frac_in(-10:20)
        integer*4 ihist_frac_out(-10:20)
        integer*4 ihist_frac_in_out(-10:20,-10:20)
        integer*4 ihist_frac_in_sat(-10:20,-10:20)
        integer*4 ihist_colmaxin_sat(-10:20,-10:20)
        integer*4 ihist_colmaxout_sat(-10:20,-10:20)
        integer*4 istat_39_a(ni,nj)
        real*4 albedo(ni,nj)
        real*4 topo(ni,nj)
        real*4 dbz_max_2d(ni,nj)
        real*4 radar_ref_3d(ni,nj,klaps)
        real*4 clouds_3d(ni,nj,nk)
        real*4 cld_hts(nk)

!       This stuff is for reading VIS data from LVD file
        real*4 cloud_frac_vis_a(ni,nj)
        integer*4 mxstn
        parameter (mxstn = 100)       ! max number of "stations" in data file

        write(6,*)' subroutine insert_vis...'

!       Initialize histograms
        do i = -10,20
            ihist_frac_sat(i) = 0
            ihist_frac_in(i) = 0
            ihist_frac_out(i) = 0

            do j = -10,20
                ihist_frac_in_out(i,j) = 0
                ihist_frac_in_sat(i,j) = 0
                ihist_colmaxin_sat(i,j) = 0
                ihist_colmaxout_sat(i,j) = 0
            enddo ! j

        enddo ! i

        n_missing_albedo = 0
        n_missing_uprb = 0
        n_vis_mod = 0
        n_39_mod = 0

        diffin_sum  = 0.
        diffout_sum = 0.
        diffin_sumsq  = 0.
        diffout_sumsq = 0.

!       Horizontal array loop
        do i = 1,ni
        do j = 1,nj

          if(cloud_frac_vis_a(i,j) .ne. r_missing_data)then
              cloud_frac_uprb = cloud_frac_vis_a(i,j)
          else
              n_missing_albedo =  n_missing_albedo + 1
  
              if(istat_39_a(i,j) .eq. -1)then
!                 cloud_frac_uprb = 0.
                  cloud_frac_uprb = r_missing_data
              else
                  cloud_frac_uprb = r_missing_data
              endif

          endif

          if(cloud_frac_uprb .ne. r_missing_data)then

            iscr_frac_sat = nint(cloud_frac_uprb*10.)
            iscr_frac_sat = min(max(iscr_frac_sat,-10),20)
            ihist_frac_sat(iscr_frac_sat) = 
     1      ihist_frac_sat(iscr_frac_sat) + 1

!           Make sure satellite cloud fraction is between 0 and 1
            if(cloud_frac_uprb .le. 0.0)cloud_frac_uprb = 0.0
            if(cloud_frac_uprb .ge. 1.0)cloud_frac_uprb = 1.0

            iscr_frac_sat = nint(cloud_frac_uprb*10.)

            colmaxin = 0.
            colmaxout = 0.

            iblank_radar = 0
            iset_vis = 0

            do k = 1,nk
                if(clouds_3d(i,j,k) .gt. 1.0)then
!                   cloud_frac_in = 1.0
                    cloud_frac_in = clouds_3d(i,j,k)
                    write(6,*)' Warning, clouds_3d > 1'
     1                       ,i,j,k,clouds_3d(i,j,k)
!                   stop
                else
                    cloud_frac_in = clouds_3d(i,j,k)
                endif

!               Modify the cloud field with the vis input - allow .3 vis err?
                if(cld_hts(k) .gt. topo(i,j) + surface_sao_buffer)then
                    cushion = 0.3
                else
                    cushion = 0.0
                endif

!               if(cloud_frac_in - cloud_frac_uprb .gt. cushion)then
                if(clouds_3d(i,j,k) - cloud_frac_uprb .gt. cushion)then     
                   cloud_frac_out = cloud_frac_uprb

!                  Determine if we need to reconcile VIS with radar
                   if(      istat_radar .eq. 1
     1                .and. dbz_max_2d(i,j) .ne. r_missing_data
     1                .and. dbz_max_2d(i,j) .gt. ref_base
     1                                               )then ! Valid radar echo
                       if(cloud_frac_out .lt. vis_radar_thresh_cvr)then
                           if(dbz_max_2d(i,j) .lt. 
     1                                        vis_radar_thresh_dbz)then

!                              Blank out Radar, Normal VIS Clearing
                               iblank_radar = 1

                           else ! Set cf = VIS THRESH
                               if(.true.)then ! Take action
                                   cloud_frac_out = vis_radar_thresh_cvr
                               endif
                               iset_vis = 1
                           endif
                       endif
                   endif

                   if(cloud_frac_in - cloud_frac_out .gt. .01)then
                       if(cloud_frac_vis_a(i,j) .ne. r_missing_data)then
                           n_vis_mod = n_vis_mod + 1
                       else
                           n_39_mod = n_39_mod + 1
                       endif
                   endif

                   clouds_3d(i,j,k) = cloud_frac_out   ! Modify the output
                else
                   cloud_frac_out = cloud_frac_in
                endif

!               Update Histograms
                iscr_frac_in = nint(cloud_frac_in*10.)
                ihist_frac_in(iscr_frac_in) = 
     1          ihist_frac_in(iscr_frac_in) + 1

                iscr_frac_out = nint(cloud_frac_out*10.)
                ihist_frac_out(iscr_frac_out) =
     1          ihist_frac_out(iscr_frac_out) + 1

                ihist_frac_in_sat(iscr_frac_in,iscr_frac_sat) =
     1          ihist_frac_in_sat(iscr_frac_in,iscr_frac_sat) + 1

                ihist_frac_in_out(iscr_frac_in,iscr_frac_out) =
     1          ihist_frac_in_out(iscr_frac_in,iscr_frac_out) + 1

                colmaxin  = max(colmaxin,cloud_frac_in)
                colmaxout = max(colmaxout,cloud_frac_out)

            enddo ! k

!           Reconcile VIS with radar
            if(iblank_radar .eq. 1)then ! NO VIS / WEAK ECHO

!               Blank out radar column for this grid point

                if(colmaxout .le. vis_radar_thresh_cvr)then
                    write(6,1)i,j,colmaxout,dbz_max_2d(i,j)
     1                       ,cloud_frac_uprb
1                   format(
     1              ' VIS_RDR - Blank out radar: cvr/dbz/vis      < '       
     1                    ,2i4,f8.2,f8.1,f8.2)

                else ! Some of the column remains above the VIS threshold
                     ! We are "saved" by the VIS cushion
                     ! May not show up in comparisons

                    write(6,2)i,j,colmaxout,dbz_max_2d(i,j)
     1                       ,cloud_frac_uprb
2                   format(
     1              ' VIS_RDR - Blank out radar: cvr/dbz/vis-s    < '
     1                    ,2i4,f8.2,f8.1,f8.2)

                endif

                if(.true.)then        ! Take action
                    dbz_max_2d(i,j) = ref_base
                    do kl = 1,klaps
                        radar_ref_3d(i,j,kl) = ref_base
                    enddo ! kl
                endif

            elseif(iset_vis .eq. 1)then ! NO VIS / STRONG ECHO
!               Cloud cvr has been reset to threshold value above VIS

                if(colmaxout .le. vis_radar_thresh_cvr)then
                    write(6,3)i,j,colmaxout,dbz_max_2d(i,j)
     1                       ,cloud_frac_uprb,vis_radar_thresh_cvr
3                   format(
     1              ' VIS_RDR - Reset vis:       cvr/dbz/vis/thr* > '
     1                    ,2i4,f8.2,f8.1,2f8.2)

                else ! Some of the column remains above the VIS threshold
                     ! We are "saved" by the VIS cushion
                     ! May not show up in comparisons
                     ! Is resetting the VIS perhaps not necessary?

                    write(6,4)i,j,colmaxout,dbz_max_2d(i,j)
     1                       ,cloud_frac_uprb,vis_radar_thresh_cvr
4                   format(
     1              ' VIS_RDR - Reset vis:       cvr/dbz/vis/thr-s> '
     1                    ,2i4,f8.2,f8.1,2f8.2)


                endif

            endif

            if(j .eq. (j-9/10)*10+9
!           if(j .gt. 28 .and. j .lt. 38
     1                  .and. colmaxin-colmaxout .gt. 0.1)then
                write(6,12)i,j,colmaxin,colmaxout
12              format(1x,'Vismod',2i4,2f5.2)
            endif

            iscr_colmaxin  = nint(colmaxin*10.)
            iscr_colmaxout = nint(colmaxout*10.)
            ihist_colmaxin_sat(iscr_colmaxin,iscr_frac_sat)
     1    = ihist_colmaxin_sat(iscr_colmaxin,iscr_frac_sat) + 1
            ihist_colmaxout_sat(iscr_colmaxout,iscr_frac_sat)
     1    = ihist_colmaxout_sat(iscr_colmaxout,iscr_frac_sat) + 1

            diffin  = colmaxin  - cloud_frac_uprb
            diffout = colmaxout - cloud_frac_uprb
            diffin_sum  = diffin_sum  + diffin
            diffout_sum = diffout_sum + diffout
            diffin_sumsq  = diffin_sumsq  + diffin**2
            diffout_sumsq = diffout_sumsq + diffout**2

          else ! missing upper bound data
            n_missing_uprb =  n_missing_uprb + 1

          endif ! if upper bound value is missing

        enddo ! i
        enddo ! j

        write(6,*)
        write(6,*)' N_MISSING_ALBEDO = ',n_missing_albedo
        write(6,*)' N_MISSING_UPRB = ',n_missing_uprb
        write(6,*)' N_VIS_MOD = ',n_vis_mod
        write(6,*)' N_39_MOD = ',n_39_mod
        write(6,*)

        write(6,*)'              HISTOGRAMS'
        write(6,*)' I          ',
     1  ' Albedo  Cld Frac Sat  Cld Frac In  Cld Frac Out'
        do i = -5,15
            write(6,11)i,ihist_alb(i),ihist_frac_sat(i),ihist_frac_in(i)
     1                          ,ihist_frac_out(i)
11          format(i4,i12,i12,i12,i12)
        enddo ! i

        write(6,*)
        write(6,*)
     1  '               Input vs. Satellite Cloud Fraction Histogram'       
        write(6,*)
        write(6,*)'                           SATELLITE CLOUD FRACTION'
        do i = 0,10
            write(6,21)(ihist_frac_in_sat(i,j),j=0,10)
21          format(1x,11i7)
        enddo ! i

        write(6,*)
        write(6,*)
     1  '                  Input vs. Output Cloud Fraction Histogram'
        write(6,*)
        write(6,*)'                             OUTPUT CLOUD FRACTION'
        do i = 0,10
            write(6,21)(ihist_frac_in_out(i,j),j=0,10)
        enddo ! i

        write(6,*)
        write(6,*)
     1  '          Column Max Input vs. Satellite Fraction Histogram'
        write(6,*)
        write(6,*)'                           SATELLITE CLOUD FRACTION'
        do i = 0,10
            write(6,21)(ihist_colmaxin_sat(i,j),j=0,10)
        enddo ! i

        write(6,*)
        write(6,*)
     1  '          Column Max Output vs. Satellite Fraction Histogram'
        write(6,*)
        write(6,*)'                           SATELLITE CLOUD FRACTION'
        do i = 0,10
            write(6,21)(ihist_colmaxout_sat(i,j),j=0,10)
        enddo ! i

        r_present = ni*nj - n_missing_uprb
        if(r_present .gt. 0.)then ! write stats
            write(6,31)diffin_sum/r_present,sqrt(diffin_sumsq/r_present)
31          format(' VIS STATS: Mean/RMS input residual  = ',2f8.3)
            write(6,41)
     1      diffout_sum/r_present,sqrt(diffout_sumsq/r_present)
41          format(' VIS STATS: Mean/RMS output residual = ',2f8.3)
        endif

        write(6,*)
        istatus = 1

        return
        end

