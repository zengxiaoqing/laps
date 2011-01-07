
        subroutine compare_analysis_to_rad(i4time,ni,nj
     1  ,cvr_sao_max,solar_alt
     1  ,cloud_frac_vis_a,tb8_k,t_gnd_k,t_sfc_k,cvr_max,r_missing_data
     1  ,dbz_max_2d,cld_snd,ista_snd,max_cld_snd,cld_hts,KCLOUD
     1  ,rad_s,n_cld_snd,c_stations,lat_s,lon_s,elev_s,maxstns
     1  ,swi_2d)                                                 ! O

        real cloud_frac_vis_a(ni,nj),tb8_k(ni,nj),t_gnd_k(ni,nj)
     1        ,t_sfc_k(ni,nj),cvr_max(ni,nj),cvr_sao_max(ni,nj)
     1        ,dbz_max_2d(ni,nj),solar_alt(ni,nj),swi_2d(ni,nj)

!       How much the solar radiation varies with changes in cloud fraction
        real cvr_scl_a(ni,nj) 

        real rad_clr(ni,nj)

        real airmass_2d(ni,nj)
        real pw_2d(ni,nj)
        real trans_h2o_2d(ni,nj)

        real dum_2d(ni,nj)

        real cld_snd(max_cld_snd,KCLOUD)
        integer ista_snd(max_cld_snd)
        real cld_hts(KCLOUD)

        character c_stations(maxstns)*(*)
        character stn(maxstns)*20

        real lat_s(maxstns), lon_s(maxstns), elev_s(maxstns)
        real rad_s(maxstns)
        real rad2_s(maxstns)
        real swi_s(maxstns)
        real resid_s(maxstns)
        real cvr_s(maxstns)

        real dum_s(maxstns)

        real ea(maxstns)

        integer ii_s(maxstns)
        integer jj_s(maxstns)
        integer doy

        character*3 c3_discrep
        character*1 c1_c
        character title*60, ver_file*256, filename*9

        write(6,*)' Comparing cloud/sat/sfc data to solar radiation'

        do j = 1,nj
        do i = 1,ni
            cvr_max(i,j) = min(cvr_max(i,j),1.00)

            if(solar_alt(i,j) .gt. 0.)then
                airmass_2d(i,j) = min(1./cosd(solar_alt(i,j)),40.)
            else
                airmass_2d(i,j) = 40.
            endif

        enddo ! i
        enddo ! j

        ea = 100. ! default
        resid_s = r_missing_data
        call make_fnam_lp(i4time,filename,istatus)

!       Get zenithal solar radiation
        read(filename(3:5),*)doy
        rad_dist_factor = radnorm(doy)
        write(6,*)' distance rad_factor = ',rad_dist_factor

        solar_constant = 1353.
        solar_irradiance = solar_constant * rad_dist_factor
        transmittance = 0.850 ! catchall for various types of absorption / scattering

        model = 1
        if(model .eq. 1)then
            rad_zenith_clr = solar_irradiance * transmittance
            write(6,*)' rad_zenith_clr = ',rad_zenith_clr  
        elseif(model .eq. 2)then ! Laue formula
            a = 0.14
            h = 1600. / 1000. ! elevation in km
!           direct = solar_irradiance
!    1             * ( (1. - a*h) * 0.7**(airmass_2d**0.678) + a*h )
            global = 1.1 * direct
        elseif(model .eq. 3)then ! Drexel model
            pw_2d = 1.0 ! cm
            trans_h2o_2d = 1.0 - 0.077 * (pw_2d * airmass_2d)**0.3
        endif

!       Calculate solar radiation field
!       Some type of regression to the solar radiation obs can be considered later
        do j = 1,nj
        do i = 1,ni
!           Set how much solar radiation varies with changes in cloud fraction
            if(cloud_frac_vis_a(i,j) .ne. r_missing_data)then
                cvr_scl_a(i,j) = 1.0 ! scaling where we have VIS data
            else
                cvr_scl_a(i,j) = 0.6 ! scaling without VIS data
            endif
            cvr_max(i,j) = min(cvr_max(i,j),1.00)

            if(model .eq. 1)then ! simple formula (radiation on horizontal)
                rad_clr(i,j) = max(rad_zenith_clr 
     1                       * sind(solar_alt(i,j)),0.)

            elseif(model .eq. 2)then ! start to using Laue formula 
                                     ! (normal to sun's rays)
                a = 0.14
                h = 1600. / 1000. ! elevation in km
                if(solar_alt(i_i,i_j) .gt. 0.)then
!                   direct = solar_irradiance
!    1                    * ( (1. - a*h) * 0.7**(airmass**0.678) + a*h )
                    global = 1.1 * direct
                endif

            endif

!           100% cloud cover has (1-cvr_scl) of possible solar radiation
            swi_2d(i,j) = rad_clr(i,j) 
     1                  * (1.0- (cvr_scl_a(i,j)*cvr_max(i,j)))
        enddo ! i
        enddo ! j

        if(.true.)then
            write(6,*)' cvr_scl_a(1,1) = ',cvr_scl_a(1,1)

            iwrite = 0

            n_ovc = 0
            n_tovc = 0
            n_sct = 0
            n_bkn = 0
            n_tsct = 0
            n_tbkn = 0
            ovc_sum = 0.
            tovc_sum = 0.
            sct_sum = 0.
            bkn_sum = 0.
            tsct_sum = 0.
            tbkn_sum = 0.

            sumobs = 0.
            sumanl = 0.
            sumalt = 0.
            sumresid = 0.
            sumscl = 0.
            cnt = 0.

            sumcld = 0.

            do ista = 1,maxstns  
              stn(ista) = c_stations(ista)(1:3)

!             write(6,*)ista,rad_s(ista)
              swi_s(ista) = r_missing_data 
              rad2_s(ista) = r_missing_data 

              if(rad_s(ista) .gt. 0.)then
                call latlon_to_rlapsgrid(lat_s(ista),lon_s(ista),lat,lon
     1                          ,ni,nj,ri,rj,istatus)

                i_i = nint(ri)
                i_j = nint(rj)

                ii_s(ista) = i_i
                jj_s(ista) = i_j

                if(i_i .ge. 3 .and. i_i .le. ni-2 .and.
     1             i_j .ge. 3 .and. i_j .le. nj-2            )then

                    if(iwrite .eq. iwrite/20*20)then
                        write(6,*)'sv '
                        write(6,*)'sv Sta   i    j   VIS frac tb8_k  '
     1                  //'t_gnd_k t_sfc_k cv_s_mx cvr_mx '
     1                  //'solalt 9pt  rad_an '
     1                  //'rad_ob rad_th ratio cv_sol  df'
                    endif

!                   Calculate 9pt cover
                    cvr_9pt = 0.
                    do ii = -1,1
                    do jj = -1,1
                        cvr_9pt = cvr_9pt + cvr_max(i_i+ii,i_j+jj)
                    enddo ! jj
                    enddo ! ii
                    cvr_9pt = cvr_9pt / 9.

!                   Calculate 25pt cover
                    cvr_25pt = 0.
                    do ii = -2,2
                    do jj = -2,2
                        cvr_25pt = cvr_25pt + cvr_max(i_i+ii,i_j+jj)
                    enddo ! jj
                    enddo ! ii
                    cvr_25pt = cvr_25pt / 25.

                    iwrite = iwrite + 1

                    c1_c = ' '
                    c3_discrep = '   '

                    if(solar_alt(i_i,i_j) .gt. 0.)then
                        rad_ratio = rad_s(ista) / rad_clr(i_i,i_j)
                        cv_solar = (1.0 - rad_ratio)       ! 100% cloud cover 
     1                           / (cvr_scl_a(i_i,i_j))    ! has (1-cvr_scl) of possible
                                                           ! solar radiation
                        cv_diff = cv_solar - cvr_max(i_i,i_j)

!                       Determine residual of clear sky vs observed radiation
                        resid_s(ista) = 1.0 - rad_ratio

                    else
                        rad_ratio = 0.
                        cv_solar = 0.
                        cv_diff = 0.
                        resid_s(ista) = 0.0

                    endif

                    rad2_s(ista) = rad_s(ista)            
                    swi_s(ista) = swi_2d(i_i,i_j)

                    if(swi_s(ista) .gt. 0.)then
                        radob_ratio = rad_s(ista) / swi_s(ista)         
                    else
                        radob_ratio = 1.0
                    endif

!                   QC checks
                    if(cvr_max(i_i,i_j) .le. .10 .and. 
     1                 radob_ratio .lt. 0.3      .and.
     1                 swi_s(ista) .ge. 100.           )then
                        c1_c = '-' ! Suspected low
                    endif

                    if(radob_ratio .gt. 3.0      .and.      
     1                 rad_s(ista) .ge. 400.           )then
                        c1_c = '+' ! Suspected high
                    endif

                    if(radob_ratio .lt. 0.1 .and. 
     1                 swi_s(ista) .ge. 100.      )then
                        c1_c = '*' ! QC'd out
                        rad2_s(ista) = r_missing_data
                    endif

                    if(rad_s(ista) - rad_clr(i_i,i_j) .gt. 500.)then
                       if(rad_clr(i_i,i_j) .gt. 100.)then
                          if(rad_s(ista)/rad_clr(i_i,i_j) .gt. 2.5)then
                             c1_c = '*'
                             rad2_s(ista) = r_missing_data
                          endif
                       endif
                    endif

                    write(6,1111,err=1112)'sv',c_stations(ista)(1:3)
     1                           ,i_i,i_j
     1                           ,cloud_frac_vis_a(i_i,i_j)
     1                           ,tb8_k(i_i,i_j)
     1                           ,t_gnd_k(i_i,i_j)
     1                           ,t_sfc_k(i_i,i_j)
     1                           ,cvr_sao_max(i_i,i_j)
     1                           ,cvr_max(i_i,i_j)
     1                           ,solar_alt(i_i,i_j)
     1                           ,cvr_9pt
!    1                           ,cvr_25pt
     1                           ,swi_2d(i_i,i_j)
     1                           ,rad_s(ista)
     1                           ,rad_clr(i_i,i_j)
     1                           ,rad_ratio
     1                           ,cv_solar
     1                           ,cv_diff
     1                           ,c1_c
1111                format(1x,a2,1x,a3,2i5,f8.2,3f8.1,f7.2,f8.2,f6.1       
     1                    ,f6.2,f8.1,2f7.1,3f6.2,1x,a,1x)

                    sumobs = sumobs + rad_s(ista)
                    sumanl = sumanl + swi_2d(i_i,i_j)
                    sumcld = sumcld + cvr_max(i_i,i_j)
                    sumalt = sumalt + solar_alt(i_i,i_j)
                    sumresid = sumresid + resid_s(ista) 
                    sumscl = sumscl + cvr_scl_a(i_i,i_j)
                    cnt = cnt + 1.

                    cvr_s(ista) = cvr_max(i_i,i_j)

1112            endif ! ob is in domain
              endif ! ista .ne. 0 (valid value)
            enddo ! isnd

            write(6,*)
            write(6,*)' Generic stats:'
            call stats_1d(maxstns,swi_s,rad2_s
     1                   ,'Solar Radiation (QCed): '
     1                   ,a_t,b_t,xbar,ybar
     1                   ,bias,std,r_missing_data,istatus)

!           Calculate other stats
            if(cnt .gt. 0.)then
                write(6,*)' sw radiation comparison stats'
                write(6,*)' obs / anl ratio = ',sumobs/sumanl
                write(6,801)sumcld/cnt,sumalt/cnt,xbar,ybar
801             format(
     1          '  means: cloud frac, solar alt = ',f7.2,f8.1,3x,
     1          '  analyzed, observed radiation = ',2f9.2)

                write(6,802)sumresid/sumcld, sumscl/cnt
802             format('  sensitivity of radiation to cloud '
     1                ,'fraction - measured/used: ',2f8.3)
            endif

            write(6,*)
            write(6,*)' Do regression on residuals vs cloud fraction'
            call regression(maxstns,cvr_s,resid_s
     1                   ,'Residuals vs cloud fraction: '
     1                   ,a_t,b_t,xbar,ybar
     1                   ,bias,std,r_missing_data,istatus)

!           Open file and call verification routine
            iunit = 11
            call get_directory('log', ver_file, len)
            ver_file = ver_file(1:len)//'qc/laps_sol.ver.'
     1                                //filename(6:9)
            call s_len(ver_file, len)
            open(iunit,file=ver_file(1:len),status='unknown',err=990)

            title = 'Solar Radiation (before qc): '

            call verify(swi_2d,rad_s,stn,maxstns,title,iunit 
     &                 ,ni,nj,maxstns,dum_s,dum_s,dum_2d
     &                 ,ii_s,jj_s,ea,badflag)

            close(iunit)

990         continue

        endif ! do radiation comparison to analysis

        return
        end

