


       subroutine put_sfc_bal(i4time,t_bal,ht_bal,u_bal,v_bal         ! Input
     1                       ,topo,ni,nj,nk                           ! Input
!    1                       ,rlat,rlon                               ! Input
     1                       ,istatus                              )  ! Output

       real ht_bal(ni,nj,nk)
       real u_bal(ni,nj,nk)
       real v_bal(ni,nj,nk)
       real t_bal(ni,nj,nk)

       real topo(ni,nj),rk_terrain(ni,nj)
       real rlat(ni,nj), rlon(ni,nj)

       real u_bal_sfc(ni,nj)
       real v_bal_sfc(ni,nj)
       real t_bal_sfc(ni,nj)

       real div_bal_sfc(ni,nj)
       real vort_bal_sfc(ni,nj)
       real dx(ni,nj), dy(ni,nj)

       real mslp_bal_sfc(ni,nj)
       real redp_bal_sfc(ni,nj)
       real stnp_bal_sfc(ni,nj)

       real pres_3d(ni,nj,nk)
       real pres_1d(nk)

       real td_bal_sfc(ni,nj)
       real rh_bal_sfc(ni,nj) ! 0-1
       real mr_bal_sfc(ni,nj)
       real sat_t(ni,nj)
       real sat_td(ni,nj)

       real k_to_c

       character*150 dir_in,dir_out

       integer n_sfc
       parameter (n_sfc=24)

       character*3 var(n_sfc),ext
       character*10 lvl_coord_req(n_sfc)
       character*10 units_req(n_sfc)
       character*125 comment_req(n_sfc)
       integer lvl_req(n_sfc)

       real data(ni,nj,n_sfc)

       write(6,*)' subroutine put_sfc_bal...'

       ext = 'LSX'
       
       var(1) = 'U'		! u-wind (m/s)
       var(2) = 'V'		! v-wind (m/s)
       var(3) = 'P'		! reduced press (Pa)
       var(4) = 'T'		! temp (K)
       var(5) = 'TD'		! dew point (K)
       var(6) = 'VV'		! vert. vel (m/s)
       var(7) = 'RH'		! relative humidity (%)
       var(8) = 'HI'		! Heat Index (K)
       var(9) = 'MSL'		! MSL pressure (Pa)
       var(10) = 'TAD'		! temperature advection (K/sec)
       var(11) = 'TH'		! potential temp (K)
       var(12) = 'THE'		! equiv pot temp (K)
       var(13) = 'PS'		! surface press (Pa)
       var(14) = 'VOR'		! sfc vorticity (/s)
       var(15) = 'MR'		! mixing ratio (g/kg)
       var(16) = 'MRC'		! moisture convergence (g/kg/s)
       var(17) = 'DIV'		! sfc divergence (/s)
       var(18) = 'THA'		! pot temp adv (K/s)
       var(19) = 'MRA'		! moisture adv (g/kg/s)
       var(20) = 'SPD'		! wind speed (m/s)
       var(21) = 'CSS'		! CSSI 
       var(22) = 'VIS'		! Visibility (m)
       var(23) = 'FWX'		! Fire threat index (integer)
       var(24) = 'TGD'		! Ground Temperature

       lvl_req = 0 ! set entire array to 0
 
!      Read Pre-balanced LSX surface analysis
       call get_directory(ext,dir_in,len_dir_in)
       call read_laps_data(i4time,dir_in,ext,ni,nj,n_sfc,n_sfc,     ! I
     1                     var,lvl_req,                             ! I
     1                     lvl_coord_req,                           ! O
     1                     units_req,comment_req,data,istatus)      ! O

!      Determine terrain location in 3-D grid
       do j = 1,nj
       do i = 1,ni

!          Interpolate from three dimensional grid to terrain surface
           zlow = height_to_zcoord2(topo(i,j),ht_bal,ni,nj,nk
     1                                                  ,i,j,istatus)
           if(istatus .ne. 1)then
               write(6,*)' lapswind_anal: error in height_to_zcoord2'
     1                  ,' in sfc wind interpolation',istatus
               write(6,*)i,j,zlow,topo(i,j),
     1                   (ht_bal(i,j,k),k=1,nk)
               istatus = 0
               return
           endif

           rk_terrain(i,j) = zlow

       enddo ! i
       enddo ! j

!      Perform interpolation for U and V wind components
       do j = 1,nj
       do i = 1,ni
           klow = max(rk_terrain(i,j),1.)
           khigh = klow + 1
           fraclow  = 0. ! float(khigh) - rk_terrain(i,j)
           frachigh = 1. ! 1.0 - fraclow

           if(u_bal(i,j,khigh) .eq. 1e-30)then
               khigh = khigh + 1
           endif

           if( u_bal(i,j,klow)  .eq. r_missing_data
     1    .or. v_bal(i,j,klow)  .eq. r_missing_data
     1    .or. u_bal(i,j,khigh) .eq. r_missing_data
     1    .or. v_bal(i,j,khigh) .eq. r_missing_data        )then

               write(6,3333)i,j
3333           format(' Warning: cannot interpolate to sfc at ',2i3)
               u_bal_sfc(i,j) = r_missing_data
               v_bal_sfc(i,j) = r_missing_data

           else
               u_bal_sfc(i,j) = u_bal(i,j,klow ) * fraclow
     1                        + u_bal(i,j,khigh) * frachigh

               v_bal_sfc(i,j) = v_bal(i,j,klow ) * fraclow
     1                        + v_bal(i,j,khigh) * frachigh

           endif

       enddo ! j
       enddo ! i

!      Recalculate Vorticity and Divergence
!      call get_grid_spacing_array(rlat,rlon,,ni,nj,dx,dy)
!      call vortdiv(u_bal_sfc,v_bal_sfc,vort_bal_sfc,div_bal_sfc
!    1             ,ni,nj,dx,dy)

!      Perform interpolation for Temperature
       do j = 1,nj
       do i = 1,ni
           klow = max(rk_terrain(i,j),1.)
           khigh = klow + 1
           fraclow = float(khigh) - rk_terrain(i,j)
           frachigh = 1.0 - fraclow

           if( t_bal(i,j,klow)  .eq. r_missing_data
     1    .or. t_bal(i,j,khigh) .eq. r_missing_data        )then

               write(6,3333)i,j
               t_bal_sfc(i,j) = r_missing_data

           else
               t_bal_sfc(i,j) = t_bal(i,j,klow ) * fraclow
     1                        + t_bal(i,j,khigh) * frachigh

           endif

       enddo ! j
       enddo ! i

!      Perform interpolation for height / surface pressure / MSLP/ reduced P
!      possibly use height to pressure routine in conversions.f

       call get_laps_redp(redp_lvl,istatus)
       if(istatus .ne. 1) return

       call get_pres_3d(i4time,ni,nj,nk,pres_3d,istatus)
       if(istatus .ne. 1) return

       do j = 1,nj
       do i = 1,ni
           ht_mslp = 0.
           ht_redp = redp_lvl
           ht_stnp = topo(i,j)

           pres_1d(:) = pres_3d(i,j,:)

           mslp_bal_sfc(i,j) = height_to_pressure(ht_mslp,ht_bal
     1                                         ,pres_1d,ni,nj,nk,i,j)

           redp_bal_sfc(i,j) = height_to_pressure(ht_redp,ht_bal
     1                                         ,pres_1d,ni,nj,nk,i,j)

           stnp_bal_sfc(i,j) = height_to_pressure(ht_stnp,ht_bal
     1                                         ,pres_1d,ni,nj,nk,i,j)


       enddo ! i
       enddo ! j

!      Recalculate moisture variables (td,rh,mr)
       td_bal_sfc = data(:,:,5) 

!      Constrain dewpoint
       do i = 1,ni
       do j = 1,nj

!          Calculate RH from Td and T
           if(td_bal_sfc(i,j) .ne. r_missing_data .and.
     1        t_bal_sfc(i,j)  .ne. r_missing_data        )then

!             Avoid supersaturation
              td_bal_sfc(i,j) = min(td_bal_sfc(i,j),t_bal_sfc(i,j)) 

           endif

       enddo ! j
       enddo ! i

!      Calculate humidity
       call hum(t_bal_sfc,td_bal_sfc,rh_bal_sfc,ni,nj,sat_t,sat_td)

!      Note W = 622. * x / (p-x), where x is ambient vapor pressure in mb
       mr_bal_sfc = (622. * sat_td) / (stnp_bal_sfc/100. - sat_td)

!      Substitute balanced fields in data array
       data(:,:,1) = u_bal_sfc
       data(:,:,2) = v_bal_sfc
       data(:,:,4) = t_bal_sfc
       data(:,:,5) = td_bal_sfc
       data(:,:,7) = rh_bal_sfc * 100. ! convert to percent
       data(:,:,9) = mslp_bal_sfc
       data(:,:,3) = redp_bal_sfc
       data(:,:,13) = stnp_bal_sfc
!      data(:,:,14) = vort_bal_sfc
       data(:,:,15) = mr_bal_sfc
!      data(:,:,17) = div_bal_sfc

!      Write balanced LSX file
       call get_directory('balance/lsx',dir_out,len_dir_out)
       print*,'Writing out balanced lsx...'
       call write_laps_data(i4time,dir_out,ext,ni,nj,n_sfc,n_sfc,    ! I
     1                      var,lvl_req,                             ! I
     1                      lvl_coord_req,                           ! O
     1                      units_req,comment_req,data,istatus)      ! I/O
       if (istatus .NE. 1) then
         print*,'Problem writing out lsx!'
       endif

       istatus = 1

       return
       end
 
       
   
