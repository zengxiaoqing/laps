
        subroutine compare_analysis_to_rad(ni,nj,cvr_sao_max,solar_alt
     1  ,cloud_frac_vis_a,tb8_k,t_gnd_k,t_sfc_k,cvr_max,r_missing_data
     1  ,dbz_max_2d,cld_snd,ista_snd,max_cld_snd,cld_hts,KCLOUD
     1  ,rad_s,n_cld_snd,c_stations,lat_s,lon_s,elev_s,maxstns)

        real cloud_frac_vis_a(ni,nj),tb8_k(ni,nj),t_gnd_k(ni,nj)
     1        ,t_sfc_k(ni,nj),cvr_max(ni,nj),cvr_sao_max(ni,nj)
     1        ,dbz_max_2d(ni,nj),solar_alt(ni,nj)

        real cld_snd(max_cld_snd,KCLOUD)
        integer ista_snd(max_cld_snd)
        real cld_hts(KCLOUD)

        character c_stations(maxstns)*(*)
        real lat_s(maxstns), lon_s(maxstns), elev_s(maxstns)
        real rad_s(maxstns)

        character*3 c3_discrep
        character*1 c1_c

        do j = 1,nj
        do i = 1,ni
            cvr_max(i,j) = min(cvr_max(i,j),1.00)
        enddo ! i
        enddo ! j

        if(.true.)then
            write(6,*)
     1      ' Comparing cloud/sat/sfc data to solar radiation'
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


            do ista = 1,maxstns  
!             write(6,*)ista,rad_s(ista)
              if(rad_s(ista) .gt. 0.)then
                call latlon_to_rlapsgrid(lat_s(ista),lon_s(ista),lat,lon
     1                          ,ni,nj,ri,rj,istatus)

                i_i = nint(ri)
                i_j = nint(rj)

                if(i_i .ge. 3 .and. i_i .le. ni-2 .and.
     1             i_j .ge. 3 .and. i_j .le. nj-2            )then

                    if(iwrite .eq. iwrite/20*20)then
                        write(6,*)
                        write(6,*)'Sta   i    j   VIS frac tb8_k  '
     1                  //'t_gnd_k t_sfc_k cldsnd cv_sa_mx cvr_mx '
     1                  //'solalt 9pt  25pt  '
     1                  //' rad   rad_th ratio cv_sol  df'
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
                        rad_clr = 1150. * sind(solar_alt(i_i,i_j))
                        rad_ratio = rad_s(ista) / rad_clr
                        cv_solar = (1.0 - rad_ratio) / 0.8 ! 100% cloud cover 
                                                           ! has 20% of possible
                                                           ! solar radiation
                        cv_diff = cv_solar - cvr_max(i_i,i_j)
                    else
                        rad_clr = 0.
                        rad_ratio = 0.
                        cv_solar = 0.
                        cv_diff = 0.
                    endif

                    write(6,1111,err=1112)c_stations(ista)(1:3)
     1                           ,i_i,i_j
     1                           ,cloud_frac_vis_a(i_i,i_j)
     1                           ,tb8_k(i_i,i_j)
     1                           ,t_gnd_k(i_i,i_j)
     1                           ,t_sfc_k(i_i,i_j)
     1                           ,cld_snd_max
     1                           ,cvr_sao_max(i_i,i_j)
     1                           ,cvr_max(i_i,i_j)
     1                           ,solar_alt(i_i,i_j)
     1                           ,cvr_9pt
     1                           ,cvr_25pt
     1                           ,rad_s(ista)
     1                           ,rad_clr      
     1                           ,rad_ratio
     1                           ,cv_solar
     1                           ,cv_diff
1111                format(1x,a3,2i5,f8.2,3f8.1,f7.2,2f8.2,f6.1       
     1                    ,f6.2,f7.2,2f7.1,3f6.2)

1112            endif ! ob is in domain
              endif ! ista .ne. 0 (valid value)
            enddo ! isnd

            write(6,*)

        endif ! do radiation comparison to analysis

        return
        end

