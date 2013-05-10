
        subroutine compare_analysis_to_radiometer(i4time,ni,nj
     1  ,slwc_int                               
     1  ,r_missing_data
     1  ,c_stations,lat_s,lon_s,elev_s,cloud_liquid_int,maxstns)

        real slwc_int(ni,nj)                                

        character c_stations(maxstns)*(*)
        character stn(maxstns)*20
        character a24time*24

        real lat_s(maxstns), lon_s(maxstns), elev_s(maxstns)
        real cloud_liquid_int(maxstns)
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

        write(6,*)' Comparing cloud/sat/sfc data to radiometer'         

        if(.true.)then
            iwrite = 0

            sumobs = 0.
            sumanl = 0.
            sumalt = 0.
            sumresid = 0.
            sumscl = 0.
            cnt = 0.

            sumcld = 0.
            sumsnow = 0.
            sumclr = 0.

            do ista = 1,maxstns  
              stn(ista) = c_stations(ista)(1:3)

!             write(6,*)ista,cloud_liquid_int(ista)
              swi_s(ista) = r_missing_data 
              rad2_s(ista) = r_missing_data 

              if(cloud_liquid_int(ista) .ge. 0.)then ! valid value
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
     1                  //'solalt cv_r rad_an '
     1                  //'rad_ob rad_th ratio cv_sol  df'
                    endif

                    iwrite = iwrite + 1

                    c1_c = ' '
                    c3_discrep = '   '

!                   QC checks
                    if(radob_ratio .lt. 0.5      .and.
     1                 cloud_liquid_int(ista) .ge. 100.           )then
                        c1_c = '+' ! Suspected high analysis
                    endif

                    if(radob_ratio .gt. 2.0      .and.      
     1                 cloud_liquid_int(ista) .ge. 100.           )then
                        c1_c = '-' ! Suspected low analysis
                    endif

                    if(radob_ratio .lt. 0.1 .and. 
     1                 swi_s(ista) .ge. 100.      )then
                        c1_c = '*' ! QC'd out
                        rad2_s(ista) = r_missing_data
                    endif

                    if(cloud_liquid_int(ista) - rad_clr(i_i,i_j) 
     1                                              .gt. 500.)then
                       if(rad_clr(i_i,i_j) .gt. 100.)then
                          if(cloud_liquid_int(ista)/rad_clr(i_i,i_j) 
     1                                              .gt. 2.5)then
                             c1_c = '*' ! QC'd out
                             rad2_s(ista) = r_missing_data
                          endif
                       endif
                    endif

                    call qc_solar_ob(cloud_liquid_int(ista)
     1                              ,solar_alt(i_i,i_j)
     1                              ,r_missing_data,iqc
     1                              ,rad_clr(i_i,i_j))      

                    if(iqc .ne. 0)then
                       c1_c = '*' ! QC'd out
                       rad2_s(ista) = r_missing_data
                    endif

                    write(6,1111,err=1112)'sv',c_stations(ista)(1:3)
     1                           ,i_i,i_j
     1                           ,cloud_frac_vis_a(i_i,i_j)
     1                           ,tb8_k(i_i,i_j)
     1                           ,t_gnd_k(i_i,i_j)
     1                           ,t_sfc_k(i_i,i_j)
     1                           ,cvr_snow(i_i,i_j)
     1                           ,cvr_max(i_i,i_j)
     1                           ,solar_alt(i_i,i_j)
     1                           ,cvr_rad(i_i,i_j)
!    1                           ,cvr_25pt
     1                           ,swi_2d(i_i,i_j)
     1                           ,cloud_liquid_int(ista)
     1                           ,rad_clr(i_i,i_j)
     1                           ,rad_ratio
     1                           ,cv_solar
     1                           ,cv_diff
     1                           ,c1_c
1111                format(1x,a2,1x,a3,2i5,f8.2,3f8.1,f7.2,f8.2,f6.1       
     1                    ,f6.2,f8.1,2f7.1,3f6.2,1x,a,1x)

                    sumobs = sumobs + cloud_liquid_int(ista)
                    sumanl = sumanl + swi_2d(i_i,i_j)
                    sumcld = sumcld + cvr_rad(i_i,i_j)
                    sumsnow = sumsnow + cvr_snow(i_i,i_j)
                    sumalt = sumalt + solar_alt(i_i,i_j)
                    sumresid = sumresid + resid_s(ista) 
                    sumclr = sumclr + rad_clr(i_i,i_j)
                    cnt = cnt + 1.

                    cvr_s(ista) = cvr_rad(i_i,i_j)

1112            endif ! ob is in domain
              endif ! valid value                   
            enddo ! isnd

            write(6,*)
            write(6,*)' Generic stats:'
            call stats_1d(maxstns,swi_s,rad2_s
     1                   ,'Solar Radiation (QCed): '
     1                   ,a_t,b_t,xbar,ybar
     1                   ,bias,std,r_missing_data,istatus)

!           Write out line of stats for gnuplot (if data are valid)
            call cv_i4tim_asc_lp(i4time,a24time,istatus)
            if(cnt .gt. 0)then
                write(6,710)a24time,xbar,ybar,std,sumclr/cnt
            else
                write(6,*)' Write missing data due to zero count'
                pmsg = -99.9
                write(6,710)a24time,pmsg,pmsg,pmsg,pmsg       
            endif
710         format(1x,a24,4f10.3,' gnuplot')

!           Calculate other stats
            if(cnt .gt. 0.)then
                write(6,*)' radiometer comparison stats'
                if(sumanl .gt. 0.)then
                    write(6,*)' obs / anl ratio = ',sumobs/sumanl
                else
                    write(6,*)' obs / anl ratio not computed'
                endif

                write(6,801)sumcld/cnt,sumsnow/cnt,sumalt/cnt,xbar,ybar
     1                     ,sumclr/cnt
801             format(
     1          '  means: cloud frac, snow cover, solar alt = '
     1          ,2f7.2,f8.1,2x,
     1          '  analyzed, observed, clr_sky rad =',3f8.1)

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

            call verify(swi_2d,cloud_liquid_int,stn,maxstns,title,iunit 
     &                 ,ni,nj,maxstns,dum_s,dum_s,dum_2d
     &                 ,ii_s,jj_s,ea,badflag)

            close(iunit)

990         continue

        endif ! do radiometer comparison to analysis

        return
        end

