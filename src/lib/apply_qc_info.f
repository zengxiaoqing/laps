
	subroutine apply_qc_info(n_obs_b,r_missing_data
     1                       ,mo,mf,mt
     1                       ,stn_a,bkg_a,obs_a,diff_a,nsta
     1                       ,bias_a,obs_mean,obs_std
     1                       ,mxstn,obs,stations
     1                       ,t_s,t_ea,td_s,td_ea
     1                       ,dd_s,dd_ea,ff_s,ff_ea)

        character*5 stn_a(mo)
        real bkg_a(mo,mf,mt)
        real obs_a(mo,mf,mt)
        real diff_a(mo,mf,mt)
        real bias_a(mo,mf) 
        real obs_mean(mo,mf)      
        real obs_std(mo,mf) 
c
c..... Stuff for the sfc data and other station info (LSO +)
c
        include 'sfcob.inc'
        type (sfcob) obs(mxstn)

	real lat_s(mxstn), lon_s(mxstn), elev_s(mxstn)
	real t_s(mxstn), t_ea(mxstn), max24t(mxstn), min24t(mxstn)
        real td_s(mxstn), td_ea(mxstn), rh_s(mxstn), rh_ea(mxstn)

        real dd_s(mxstn), ddg_s(mxstn), dd_ea(mxstn)
        real ff_s(mxstn), ffg_s(mxstn), ff_ea(mxstn)

        real alt_s(mxstn), alt_ea(mxstn), delp(mxstn)
	real pstn_s(mxstn), pmsl_s(mxstn), p_ea(mxstn), pred_s(mxstn)

	real store_hgt(mxstn,5) 

        real vis_s(mxstn), vis_ea(mxstn)
        real solar_s(mxstn), solar_ea(mxstn)

        real sfct(mxstn), sfct_ea(mxstn)
        real sfcm(mxstn), sfcm_ea(mxstn)
        real pcp1(mxstn), pcp3(mxstn), pcp6(mxstn), pcp24(mxstn)
        real snow(mxstn), snow_ea(mxstn), pcp_ea(mxstn)

	real rii(mxstn), rjj(mxstn)
c
	integer kloud_s(mxstn), obstime(mxstn)
        integer wmoid(mxstn), delpch(mxstn)
	integer ii(mxstn), jj(mxstn)
c
	character atime_s*24
	character store_amt(mxstn,5)*4
        character stations(mxstn)*20, provider(mxstn)*11,stn20*20      
        character reptype(mxstn)*6, autostntype(mxstn)*6
        character wx_s(mxstn)*25 

        real lapse_t, lapse_td

        if(n_obs_b .gt. 0)then

!           Apply QC information                            
            iv_t   = 1
            iv_td  = 2
            iv_spd = 5
            iv_dir = 9

            lapse_t = -.01167
            lapse_td = -.007

            do mm = 1,n_obs_b
              do ista = 1,nsta
                if(stations(mm)(1:5) .eq. stn_a(ista))then
                   goto 130
                endif
              enddo 
              goto 150
130           continue
              if(mm .le. 50)then
                write(6,*)' Found a QC match',mm,ista,stn_a(ista)
              endif

!             Test biases and flag ob using expected accuracy info

!             Temperature bias
              if(bias_a(ista,iv_t) .ne. r_missing_data)then
                biast_corr  = bias_a(ista,iv_t)  
     1                      - lapse_t *obs(mm)%elev_diff

                if(abs(biast_corr) .gt. 8.0)then
                  if(t_s(mm) .ne. badflag)then
                    obs(mm)%t_ea_f = abs(biast_corr)
                    t_ea(mm)       = abs(biast_corr)
                    write(6,135)mm,ista,obs(mm)%i,obs(mm)%j
     1                       ,stn_a(ista),biast_corr
     1                       ,bias_a(ista,iv_t),obs(mm)%elev_diff
     1                       ,obs(mm)%ldf
135		    format(' setting t_ea based on bias '
     1                    ,i6,i6,4x,2i5,1x,a5,1x,2f8.2,f8.1,f7.2)
                  endif                
                endif
              endif

!             Stuck temperature
              if(obs_std(ista,iv_t) .ne. r_missing_data)then
                if(obs_std(ista,iv_t)  .lt. 0.1)then
                  if(t_s(mm) .ne. badflag)then
                    obs(mm)%t_ea_f = 50.                  
                    t_ea(mm)       = 50.                 
                    write(6,136)mm,ista,obs(mm)%i,obs(mm)%j
     1                       ,stn_a(ista),obs_std(ista,iv_t) 
136		    format(' stuck temperature for   '
     1                    ,i6,i6,4x,2i5,1x,a5,1x,f8.1,'deg')
                    do it = 1,mt              
                      write(6,*)' stuck temp time ',it 
     1                         ,obs_a(ista,iv_t,it)
     1                         ,obs_a(ista,iv_t,it)
                    enddo
                  endif
                endif
              endif

!             Dewpoint bias
              if(bias_a(ista,iv_td) .ne. r_missing_data)then
                biastd_corr = bias_a(ista,iv_td) 
     1                      - lapse_td*obs(mm)%elev_diff
      
                if(abs(biastd_corr) .gt. 8.0)then
                  if(td_s(mm) .ne. badflag)then
                    obs(mm)%td_ea_f = abs(biastd_corr)
                    td_ea(mm)       = abs(biastd_corr)
                    write(6,137)mm,ista,obs(mm)%i,obs(mm)%j
     1                       ,stn_a(ista),biastd_corr
     1                       ,bias_a(ista,iv_td),obs(mm)%elev_diff
137		    format(' setting td_ea based on bias'
     1                    ,i6,i6,4x,2i5,1x,a5,1x,2f8.2,f8.1)
                  endif                
!                 if(abs(biastd_corr) .gt. 1000.)then
!                   do it = 1,mt              
!                     write(6,*)' td processing error',it 
!    1                         ,bkg_a(ista,iv_td,it)
!    1                         ,obs_a(ista,iv_td,it)
!    1                         ,diff_a(ista,iv_td,it)
!                   enddo
!                 endif
                endif
              endif

!             Stuck wind direction
              if(obs_std(ista,iv_dir) .ne. r_missing_data)then
                if(obs_std(ista,iv_dir)  .lt. 3.0 .and. 
     1             obs_mean(ista,iv_spd) .gt. 0.)then
                  if(dd_s(mm) .ne. badflag)then
                    obs(mm)%dd_ea_deg = 180.                  
                    dd_ea(mm)         = 180.                 
                    write(6,138)mm,ista,obs(mm)%i,obs(mm)%j
     1                       ,stn_a(ista),obs_mean(ista,iv_dir) 
138		    format(' stuck wind direction for   '
     1                    ,i6,i6,4x,2i5,1x,a5,1x,f8.1,'deg')
                    do it = 1,mt              
                      write(6,*)' stuck time ',it 
     1                         ,obs_a(ista,iv_dir,it)
     1                         ,obs_a(ista,iv_spd,it)
                    enddo
                  endif
                endif
              endif

!             Stuck wind speed
              if(obs_mean(ista,iv_spd) .ne. r_missing_data)then
                if(obs_mean(ista,iv_spd) .le. 0.01)then   
                  if(ff_s(mm) .ne. badflag)then
                    obs(mm)%ff_ea_kt = 50.                  
                    ff_ea(mm)        = 50.                 
                    write(6,139)mm,ista,obs(mm)%i,obs(mm)%j
     1                       ,stn_a(ista),obs_mean(ista,iv_spd)               
139		    format(' stuck wind speed for       '
     1                    ,i6,i6,4x,2i5,1x,a5,1x,f8.1,'kt')
                    do it = 1,mt              
                      write(6,*)' stuck time ',it 
     1                         ,obs_a(ista,iv_dir,it)
     1                         ,obs_a(ista,iv_spd,it)
                    enddo
                  endif
                endif
              endif

150         enddo ! mm

            I4_elapsed = ishow_timer()

        endif

        return
        end
