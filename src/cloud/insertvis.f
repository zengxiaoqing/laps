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
     1      ,topo,cloud_frac_vis_a,albedo,mode_refl,ihist_alb         ! I
     1      ,istat_39_a,l_use_39,idb,jdb                              ! I
     1      ,ni,nj,nk,r_missing_data                                  ! I
     1      ,vis_radar_thresh_cvr,vis_radar_thresh_dbz                ! I
     1      ,istat_radar,radar_ref_3d,klaps,ref_base
     1      ,solar_alt,solar_az
     1      ,di_dh,dj_dh,i_fill_seams                                 ! I
     1      ,cldtop_tb8_m,cldht_prlx_top                              ! I
     1      ,cloud_albedo,cloud_od,cloud_op                           ! O
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
        real albedo(ni,nj)            ! Measured from satellite
        real cloud_albedo(ni,nj)      ! Cloud albedo
        real cloud_od(ni,nj)          ! Cloud optical depth
        real cloud_op(ni,nj)          ! Cloud opacity
        real topo(ni,nj)
        real dbz_max_2d(ni,nj)
        real radar_ref_3d(ni,nj,klaps)
        real clouds_3d(ni,nj,nk)
        real clouds_3d_buf(ni,nj,nk)
        real cld_hts(nk)
        real solar_alt(ni,nj)
        real solar_az(ni,nj)
        real r_shadow_3d(ni,nj,nk)
        real di_dh(ni,nj)                      
        real dj_dh(ni,nj)                      
        real cldtop_tb8_m(ni,nj)     ! Input
        real cldht_prlx_top(ni,nj)
        integer i_fill_seams(ni,nj)

!       This stuff is for reading VIS data from LVD file
        real cloud_frac_vis_a(ni,nj)  ! presently used rather than 'albedo'
        integer mxstn
        parameter (mxstn = 100)       ! max number of "stations" in data file

        call get_grid_spacing_cen(grid_spacing_m,istatus)

!       Initialize
        cloud_albedo = r_missing_data
        cloud_od = r_missing_data
        cloud_op = r_missing_data

        mode_prlx = 3 ! 2 (fixed) or 3 (variable)
        cldht_prlx_fixed = 10000.

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
        n_ir = 0
        n_no_ir = 0

        diffin_sum  = 0.
        diffout_sum = 0.
        diffin_sumsq  = 0.
        diffout_sumsq = 0.

        clouds_3d_buf(:,:,:) = clouds_3d(:,:,:)

!       Horizontal array loop
        do i = 1,ni
        do j = 1,nj    

!         This can be set up for mode_prlx
!         if(cldtop_tb8_m(i,j) .ne. r_missing_data)then
!             cldht_prlx_top(i,j) = cldtop_tb8_m(i,j)
!             cldht_prlx_top(i,j) = cldht_prlx_fixed
!         else
!             cldht_prlx_top(i,j) = cldht_prlx_fixed
!         endif
!         cldht_prlx_top(i,j) = min(max(cldtop_tb8_m(i,j),8000.),12000.)

          if(i .eq. idb .AND. j .eq. jdb)then
              idebug = 1
          else
              idebug = 0
          endif

!         Calculate upper bound to cloud cover (through the column)
          if(cloud_frac_vis_a(i,j) .ne. r_missing_data)then ! VIS (Daytime)
              albedo_eff = cloud_frac_vis_a(i,j)
              call albedo_to_clouds2(albedo_eff,trans,trans_i
     1                              ,cloud_od(i,j),cloud_od_i
     1                              ,cloud_op(i,j),cloud_op_i)
              cloud_frac_uprb = cloud_frac_vis_a(i,j)

              if(idebug .eq. 1)then
                  write(6,51)cloud_frac_vis_a(i,j),albedo_eff
     1                      ,albedo(i,j),trans,cloud_od(i,j)
     1                      ,cloud_op(i,j),cldht_prlx_top(i,j)
51                format(' CTR cf_vis/albeff/albsat/trans/od/op/prlx/'
     1                  ,7f9.3)
              endif

              cloud_albedo(i,j) = albedo_eff

          else                                              ! 3.9u (Nighttime)
              n_missing_albedo =  n_missing_albedo + 1
  
              if(istat_39_a(i,j) .eq. -1 .and. l_use_39)then
                  cloud_frac_uprb = 0.
              else
!                 cloud_frac_uprb = r_missing_data          ! prior strategy
                  cloud_frac_uprb = 0.5                     ! limit IR clouds
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
                if(mode_prlx .eq. 3)then
                    cldht_prlx = cldht_prlx_top(i,j) ! * 2.
                elseif(mode_prlx .eq. 2)then
                    cldht_prlx = cldht_prlx_fixed
                else
                    cldht_prlx = cld_hts(k)
                endif 

                it = i - nint(di_dh(i,j) * cldht_prlx)
                jt = j - nint(dj_dh(i,j) * cldht_prlx)
                it = max(min(it,ni),1)
                jt = max(min(jt,nj),1)

                if(it .eq. idb .AND. jt .eq. jdb)then
                  idebug = 1
                else
                  idebug = 0
                endif

                if(i_fill_seams(i,j) .ne. 0)then
                  itn = min(i,it)
                  itx = max(i,it)
!                 itn = it
!                 itx = it
                else ! consider gradients
                  if(mode_prlx .eq. 3)then
                    ip = min(i+1,ni)
                    itp = ip - nint(di_dh(ip,j)*cldht_prlx_top(ip,j))
                    itp = min(max(itp,1),ni)
                  else
                    itp = it
                  endif

!                 if(itp-it .ge. 2 .or. .true.)then
                  if(itp-it .ge. 2)then
                    itn = it
                    itx = min(it+1,ni)
                  else
                    itn = it
                    itx = it
                  endif

                  if(mode_prlx .eq. 3)then
                    jp = min(j+1,nj)
                    jtp = jp - nint(di_dh(i,jp)*cldht_prlx_top(i,jp))
                    jtp = min(max(jtp,1),nj)
                  else
                    jtp = jt
                  endif

!                 if(jtp-jt .ge. 2 .or. .true.)then
                  if(jtp-jt .ge. 2)then
                    jtn = jt
                    jtx = min(jt+1,nj)
                  else
                    jtn = jt
                    jtx = jt
                  endif
                endif

                call qc_clouds_0d(i,j,k,clouds_3d(it,jt,k)
     1                           ,ni,nj,.false.)       

                cloud_frac_in = clouds_3d(it,jt,k)

!               Consider that if 'ir_present_here' is true then cushion
!               should be non-zero (for high clouds).

!               Clear stronger with low clouds, weaker with high clouds
                thr_vis_clr = 0.10
                if(cld_hts(k) .gt. topo(i,j) + 5000.)then
                    cushion = thr_vis_clr * (1.0 - cloud_frac_uprb) ! 999. ! reduce clearing
!               Modify the cloud field with the vis input - allow .3 vis err?
                elseif(cld_hts(k).gt.topo(i,j) + surface_sao_buffer)then
                    cushion = 0.0 ! 0.3
                else
                    cushion = 0.0
                endif

!               Test input clouds against upper bound + cushion
                if(cloud_frac_in - cloud_frac_uprb .gt. cushion)then
                   cloud_frac_out = max(cloud_frac_uprb,cushion)

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

!                  Check for here & near IR cloud signature
                   nsmooth = (int(8000. / grid_spacing_m) / 2) + 1
                   il = max(i-nsmooth,1)
                   ih = min(i+nsmooth,ni)
                   jl = max(j-nsmooth,1)
                   jh = min(j+nsmooth,nj)

                   ir_present_here = 1
                   if(cldtop_tb8_m(i,j) .eq. r_missing_data)then
                       ir_present_here = 0
                   endif

                   vis_present_near = 1
                   
                   ir_present_near = 1
                   do ii = il,ih
                   do jj = jl,jh
                       if(cldtop_tb8_m(ii,jj) .eq. r_missing_data)then
                           ir_present_near = 0
                       endif
                   enddo ! jj
                   enddo ! ii

                   if(ir_present_near .eq. 1)then
                       n_ir = n_ir + 1
                   else
                       n_no_ir = n_no_ir + 1
                   endif

                   clouds_3d_buf(itn:itx,jtn:jtx,k) = cloud_frac_out ! Modify output
                   if(idebug .eq. 1)then
                       write(6,202)k,itn,jtn,itx,jtx,ir_present_here
     1                            ,cloud_frac_out,cloud_frac_in,cushion
202                    format(' cloud_frac_out modified out/in/cushion '
     1                       ,9x,6i5,3f8.2,' CTR')
                   endif

                else ! keep cloud_frac the same
                   cloud_frac_out = cloud_frac_in
                   if(idebug .eq. 1)then
                       write(6,203)k,itn,itx,jtn,jtx,ir_present_here
     1                            ,cloud_frac_out
203                    format(' cloud_frac_out kept at                 '
     1                       ,9x,6i5,f8.2,' CTR')
                   endif

                endif ! modify cloud_frac

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

                if(idebug .eq. 1)then
                    write(6,201)k,i,j,it,jt,di_dh(i,j),dj_dh(i,j)
     1                         ,cld_hts(k),cloud_frac_in,cloud_frac_out
     1                         ,cloud_frac_uprb,cushion
201                 format(
     1              ' ijk,it,jt,didh,djdh,cldht,cldfracin/out,uprb,cush'
     1                    ,i4,4i5,2f9.5,f9.1,2f7.3,2f7.3,' CTR')
                endif

            enddo ! k

            if(l_39_clr_2d)n_39_clr_2d = n_39_clr_2d + 1

!           Reconcile VIS with radar
            if(iblank_radar .eq. 1)then ! NO VIS / WEAK ECHO

!               Blank out radar column for this grid point

                if(colmaxout .le. vis_radar_thresh_cvr)then
                    if(idebug .eq. 1)then
                        write(6,1)i,j,colmaxout,dbz_max_2d(i,j)
     1                           ,cloud_frac_uprb
                    endif
1                   format(
     1              ' VIS_RDR - Blank out radar: cvr/dbz/vis      < '       
     1                    ,2i4,f8.2,f8.1,f8.2)

                else ! Some of the column remains above the VIS threshold
                     ! We are "saved" by the VIS cushion
                     ! May not show up in comparisons

                    if(idebug .eq. 1)then
                        write(6,2)i,j,colmaxout,dbz_max_2d(i,j)
     1                           ,cloud_frac_uprb
                    endif
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
                    if(idebug .eq. 1)then
                        write(6,3)i,j,colmaxout,dbz_max_2d(i,j)
     1                       ,cloud_frac_uprb,vis_radar_thresh_cvr
                    endif
3                   format(
     1              ' VIS_RDR - Reset vis:       cvr/dbz/vis/thr* > '
     1                    ,2i4,f8.2,f8.1,2f8.2)

                else ! Some of the column remains above the VIS threshold
                     ! We are "saved" by the VIS cushion
                     ! May not show up in comparisons
                     ! Is resetting the VIS perhaps not necessary?

                    if(idebug .eq. 1)then
                        write(6,4)i,j,colmaxout,dbz_max_2d(i,j)
     1                       ,cloud_frac_uprb,vis_radar_thresh_cvr
                    endif
4                   format(
     1              ' VIS_RDR - Reset vis:       cvr/dbz/vis/thr-s> '
     1                    ,2i4,f8.2,f8.1,2f8.2)


                endif

            endif

!           if(j .eq. (j-9/10)*10+9
            if(idebug .eq. 1                 
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

        clouds_3d(:,:,:) = clouds_3d_buf(:,:,:)

        pct_use_vis = (1.0 - float(n_missing_albedo)/float(ni*nj))*100.       

        write(6,*)
        write(6,101)pct_use_vis
 101    format(' Visible albedo data used over ',f6.2,'% of domain')       
        write(6,*)' N_MISSING_ALBEDO = ',n_missing_albedo
        write(6,*)' N_MISSING_UPRB = ',n_missing_uprb
        write(6,*)' N_VIS_MOD = ',n_vis_mod
        write(6,*)' N_39_CLR = ',n_39_clr
        write(6,*)' N_39_CLR_2D = ',n_39_clr_2d
        write(6,*)' N_IR = ',n_ir
        write(6,*)' N_NO_IR = ',n_no_ir
        write(6,*)

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
