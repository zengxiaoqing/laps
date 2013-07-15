cdis   
cdis    Open Source License/Disclaimer, Forecast Systems Laboratory
cdis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
cdis    
cdis    This software is distributed under the Open Source Definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    In particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - Redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - Redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - All modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - If significant modifications or enhancements are made to this
cdis    software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
cdis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
cdis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
cdis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
cdis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
cdis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
cdis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
cdis   
cdis
cdis
cdis   
cdis

        subroutine insert_vis(i4time,clouds_3d,cld_hts
     1      ,topo,cloud_frac_vis_a,albedo,ihist_alb                   ! I
     1      ,istat_39_a,l_use_39                                      ! I
     1      ,ni,nj,nk,r_missing_data                                  ! I
     1      ,vis_radar_thresh_cvr,vis_radar_thresh_dbz                ! I
     1      ,istat_radar,radar_ref_3d,klaps,ref_base
     1      ,solar_alt,solar_az
     1      ,di_dh,dj_dh                                              ! I
     1      ,dbz_max_2d,surface_sao_buffer,istatus)

!       Steve Albers 1999 Added 3.9u cloud clearing step in NULL Mode
!                         This can be turned on if desired or another cloud
!                         mask input(s) can be used for cloud clearing,
!                         especially if it applies throughout the column.

!       Visible satellite implied cloud fraction is used as an upper bound 
!       to clear out clouds.

!       Visible satellite is reconciled with the radar data in the course 
!       of using satellite to clear out clouds. If the visible satellite 
!       says no cloud and the radar echo is non-existent or <30 dBZ, then the 
!       echo (along with any pre-existing cloud) is blanked out. If the radar 
!       echo is strong then visible is overruled and a minimal value of cloud 
!       fraction of 0.2 is left in.

        use cloud_rad

        integer ihist_alb(-10:20)
        integer ihist_frac_sat(-10:20)
        integer ihist_frac_in(-10:20)
        integer ihist_frac_out(-10:20)
        integer ihist_frac_in_out(-10:20,-10:20)
        integer ihist_frac_in_sat(-10:20,-10:20)
        integer ihist_colmaxin_sat(-10:20,-10:20)
        integer ihist_colmaxout_sat(-10:20,-10:20)
        integer istat_39_a(ni,nj)
        logical l_use_39, l_39_clr_2d
        real albedo(ni,nj)
        real topo(ni,nj)
        real dbz_max_2d(ni,nj)
        real radar_ref_3d(ni,nj,klaps)
        real clouds_3d(ni,nj,nk)
        real cld_hts(nk)
        real solar_alt(ni,nj)
        real solar_az(ni,nj)
        real r_shadow_3d(ni,nj,nk)
        real di_dh(ni,nj)                      
        real dj_dh(ni,nj)                      

!       This stuff is for reading VIS data from LVD file
        real cloud_frac_vis_a(ni,nj)
        integer mxstn
        parameter (mxstn = 100)       ! max number of "stations" in data file

        if(l_use_39)then
            write(6,*)' subroutine insert_vis (with 3.9u)...'
        else
            write(6,*)' subroutine insert_vis...'
        endif

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
        n_39_clr = 0
        n_39_clr_2d = 0

        diffin_sum  = 0.
        diffout_sum = 0.
        diffin_sumsq  = 0.
        diffout_sumsq = 0.

!       Horizontal array loop
        do i = 1,ni
        do j = 1,nj    

          if(j .eq. nj/2)then
              idebug = 1
          else
              idebug = 0
          endif

!         Calculate upper bound to cloud cover (through the column)
          if(cloud_frac_vis_a(i,j) .ne. r_missing_data)then ! VIS (Daytime)
              call albedo_to_clouds(albedo(i,j),trans,cloud_od,cloud_op)
              cloud_frac_uprb = cloud_frac_vis_a(i,j)

              if(idebug .eq. 1)then
                  write(6,51)cloud_frac_vis_a(i,j),albedo(i,j)
     1                       ,trans,cloud_od,cloud_op
51                format('cf_vis/albedo/trans/od/op',5f9.3)
              endif

          else                                              ! 3.9u (Nighttime)
              n_missing_albedo =  n_missing_albedo + 1
  
              if(istat_39_a(i,j) .eq. -1 .and. l_use_39)then
                  cloud_frac_uprb = 0.
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

            l_39_clr_2d = .false.

            do k = 1,nk
                it = i - nint(di_dh(i,j) * cld_hts(k))
                jt = j - nint(dj_dh(i,j) * cld_hts(k))
                it = max(min(it,ni),1)
                jt = max(min(jt,nj),1)

                call qc_clouds_0d(i,j,k,clouds_3d(it,jt,k)
     1                           ,ni,nj,.false.)       

                cloud_frac_in = clouds_3d(it,jt,k)

                if(idebug .eq. 1 .and. k .eq. nk/2)then
                    write(6,201)i,j,it,jt,di_dh(i,j),dj_dh(i,j)
     1                         ,cloud_frac_in,cld_hts(k)
 201                format('i,j,it,jt,didh,djdh,cloud_frac_in,cldht'
     1                    ,4i6,2f9.5,f9.3,f9.1)
                endif

!               Modify the cloud field with the vis input - allow .3 vis err?
                if(cld_hts(k) .gt. topo(i,j) + surface_sao_buffer)then
                    cushion = 0.0 ! 0.3
                else
                    cushion = 0.0
                endif

                if(cloud_frac_in - cloud_frac_uprb .gt. cushion)then
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
                           n_39_clr = n_39_clr + 1
                           l_39_clr_2d = .true.
                       endif
                   endif

                   clouds_3d(it,jt,k) = cloud_frac_out   ! Modify the output
                else
                   cloud_frac_out = cloud_frac_in
                endif

!               Update Histograms
                iscr_frac_in = nint(cloud_frac_in*10.)
                if(iscr_frac_in .lt. -10)then
                    write(6,*)' Bounds error for iscr_frac_in'
     1                       ,i,j,k,cloud_frac_in
                    iscr_frac_in = -10
                endif
                ihist_frac_in(iscr_frac_in) = 
     1          ihist_frac_in(iscr_frac_in) + 1

                iscr_frac_out = nint(cloud_frac_out*10.)
                if(iscr_frac_out .lt. -10)then
                    write(6,*)' Bounds error for iscr_frac_out'
     1                       ,i,j,k,cloud_frac_out
                    iscr_frac_out = -10               
                endif
                ihist_frac_out(iscr_frac_out) =
     1          ihist_frac_out(iscr_frac_out) + 1

                ihist_frac_in_sat(iscr_frac_in,iscr_frac_sat) =
     1          ihist_frac_in_sat(iscr_frac_in,iscr_frac_sat) + 1

                ihist_frac_in_out(iscr_frac_in,iscr_frac_out) =
     1          ihist_frac_in_out(iscr_frac_in,iscr_frac_out) + 1

                colmaxin  = max(colmaxin,cloud_frac_in)
                colmaxout = max(colmaxout,cloud_frac_out)

            enddo ! k

            if(l_39_clr_2d)n_39_clr_2d = n_39_clr_2d + 1

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

        pct_use_vis = (1.0 - float(n_missing_albedo)/float(ni*nj))*100.       

        write(6,*)
        write(6,101)pct_use_vis
 101    format(' Visible albedo data used over ',f6.2,'% of domain')       
        write(6,*)' N_MISSING_ALBEDO = ',n_missing_albedo
        write(6,*)' N_MISSING_UPRB = ',n_missing_uprb
        write(6,*)' N_VIS_MOD = ',n_vis_mod
        write(6,*)' N_39_CLR = ',n_39_clr
        write(6,*)' N_39_CLR_2D = ',n_39_clr_2d
        write(6,*)

        call get_grid_spacing_cen(grid_spacing_m,istatus)
        call get_cloud_shadow(i4time,clouds_3d,cld_hts,r_shadow_3d
     1                             ,ni,nj,nk
     1                             ,solar_alt,solar_az 
     1                             ,grid_spacing_m,r_missing_data)

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


        subroutine get_cloud_shadow(i4time,clouds_3d,cld_hts,r_shadow_3d
     1                             ,ni,nj,nk
     1                             ,solar_alt,solar_az 
     1                             ,grid_spacing_m,r_missing_data)

        include 'trigd.inc'

!       Determine shadow regions of the input 3-D cloud array

        real r_shadow_3d(ni,nj,nk)
        real clouds_3d(ni,nj,nk)
        real cld_hts(nk)
        real solar_alt(ni,nj)
        real solar_az(ni,nj)

        write(6,*)' Subroutine get_cloud_shadow'

        ishadow_tot = 0                          

!       First set shadow based on approximate ceiling cover (zero order)
!       do j = 1,nj
!       do i = 1,ni
!         cvr_max = 0.
!         do k = nk,1,-1
!           cvr_max = max(cvr_max,clouds_3d(i,j,k))
!           r_shadow_3d(i,j,k) = cvr_max             
!         enddo ! k
!       enddo ! i
!       enddo ! j

        do j = 1,nj
        do i = 1,ni

!         Trace towards sun from each grid point
          solar_altitude_deg = solar_alt(i,j) 
          solar_azi_deg = solar_az(i,j)

!         Get direction cosines based on azimuth
          xcos = sind(solar_azi_deg)
          ycos = cosd(solar_azi_deg)

          ishadow = 0
          if(solar_altitude_deg .gt. 20.)then
            do k = 1,nk-1  
              r_shadow_3d(i,j,k) = 0.

              do kk = k+1,nk
                dk = kk-k
                dz = cld_hts(kk) - cld_hts(k)
                dxy = dz * tand(solar_altitude_deg)

                dx = dxy * xcos
                dy = dxy * ycos

                idelt = nint(dx / grid_spacing_m)
                jdelt = nint(dy / grid_spacing_m)
              
                if(idelt .ne. 0 .or. jdelt .ne. 0)then
                  inew = i + idelt
                  jnew = j + jdelt
                  if(inew .gt. 1 .and. inew .le. ni .AND.
     1               jnew .gt. 1 .and. jnew .le. nj      )then
                    cvr_path = clouds_3d(i+idelt,j+jdelt,kk)
                    r_shadow_3d(i,j,k)=max(r_shadow_3d(i,j,k),cvr_path)      
                  endif
                endif
              enddo ! kk

              if(r_shadow_3d(i,j,k) .gt. 0.5)then
                ishadow = 1                
                if(ishadow_tot .le. 10)then
                  write(6,*)' Found shadow ',i,j
                endif
              endif

            enddo ! k

          else
            r_shadow_3d(i,j,:) = r_missing_data

          endif ! high enough sun to use vis   

          ishadow_tot = ishadow_tot + ishadow

        enddo ! i
        enddo ! j

        write(6,*)' Number of columns with shadow = ',ishadow_tot

        return
        end
