      Program lsr_driver
c
      include      'lapsparms.cmn'
      include      'lsr_dims.inc'

      integer     n_elems(max_sat)
      integer     n_lines(max_sat)
      integer     ismsng
      
      real*4        r_channel_wavelengths(max_ch,max_sat)
      character     c_sat_id(max_sat)*6
      character     c_sounding_path(max_sat)*200

      Call get_laps_config('nest7grid',IStatus)
      if(istatus.eq.1)then
        write(*,*)'LAPS Parameters obtained'
      else
         write(*,*)'IStatus = ',IStatus,'Error - Get_LAPS_Config'
         write(*,*)'Terminating LAPS-lsr. Sat Sounder remapping'
         stop
      end if

c
c get the number of satellites and channels data/static/sat_sounder.nl
c
      call get_sat_sounder_info(n_sat,c_sat_id,
     &n_channels,c_sounding_path,n_elems,n_lines,r_channel_wave
     &lengths,ismsng,pct_req_lsr,istatus)

      if(istatus.ne.1)then
         write(6,*)'Error returned from get_sat_sounder_info'
         goto 1000
      endif

      do i=1,n_sat

         call lsr_driver_sub(NX_L_CMN,NY_L_CMN,n_channels,
     &     n_elems(i),n_lines(i),ismsng,c_sat_id(i),c_sounding_path(i),
     &     r_channel_wavelengths(i,1),pct_req_lsr,istatus)
         if(istatus.ne.1)then
            if(i.eq.n_sat)then
               write(6,*)'No data processed in lsr_driver_sub'
               write(6,*)'Finished in lsr_driver'
            else
               write(6,*)'Try for another satellite'
            endif
         else
            write(6,*)'Data processed in lsr_driver_sub'
         endif

      enddo

1000  stop
      end
c
c=======================================================
c
      subroutine lsr_driver_sub(nx_l,ny_l,n_channels,
     &nelems,nlines,ismsng,c_sat_id,c_sounding_path,rch_wvlngth,
     &pct_req_lsr,istatus)
c
      implicit none
c
      integer     i_sat
      integer     nlines
      integer     nelems
      integer     nx_l,ny_l
      integer     n_channels

      integer     icnt(n_channels)
      integer     jcnt(n_channels)
      integer     icount
      integer     ismsng

      integer     ndimx(nlines)
      integer     ndimy,ndimch

      Real*8        orbitAttitude(336)
      Real*8        lineTimeBeg(nlines,n_channels)
      Real*8        lineTimeEnd(nlines,n_channels)
      Real*8        t, f_time

      real*4        lat(nx_l,ny_l)
      real*4        lon(nx_l,ny_l)
      real*4        r_llij_lut_ri(nx_l,ny_l)
      real*4        r_llij_lut_rj(nx_l,ny_l)
      real*4        sa(nx_l,ny_l)
      real*4        sc(nx_l,ny_l)
      real*4        st(nx_l,ny_l)
      real*4        laps_data(nx_l,ny_l,n_channels)
      real*4        grid_spacing
      real*4        grid_spacing_km
      real*4        r_sndr_res_km
      real*4        r_grid_ratio
      real*4        rcount,rsb,rsg
c
      real*4        data(nx_l,ny_l,2)
      real*4        rline(nx_l,ny_l)
      real*4        rpix(nx_l,ny_l)
      real*4        sndr_rad(nelems,nlines,n_channels)
      real*4        result
      real*4        xconv,yconv
      real*4        rch_wvlngth(n_channels)
      real*4        r_missing_data
      real*4        rmintime,rmaxtime
      real*4        rltb,rlte
      real*4        pct_req_lsr

      real*4        scalingBias(nlines,n_channels)
      real*4        scalingGain(nlines,n_channels)

      Integer     ewCycles,ewIncs
      Integer     nsCycles,nsIncs
      Integer     nw_pix,nw_line
      Integer     se_pix,se_line
      Integer     istatus
      Integer     iostatus
      Integer     mstatus
      Integer     i,j,k,n
      Integer     lend
      Integer     time_spec(2)
      Integer     imcI4
      Integer     int4(2)
      Integer     imax,jmax
      Integer     i4time_data
      Integer     i4time_data_orig
      Integer     ires_x,ires_y
      Integer     imaximum(n_channels)
      Integer     iminimum(n_channels)
      Integer     i2_missing_data
      Integer     ltindex
c
      Integer     instr
      Integer     isndrdata(nelems,nlines,n_channels)
      REAL*8        wavelength(n_channels)
      Character*1   imc(4)
      Character*100 filename
      Character*255 c_filename_sat
      Character*200 c_sounding_path
      Character*255 datapath
      character*125 comment_ll(2)
      character*10  units_ll(2)
      character*3   var_ll(2)
      character*2   cch
      character*150 dir_static
      Character     f9time*9
      Character     c_sat_id*6
      Character     cid*2
c
c =============================================================
c =============================================================
c                get sndr data
c =============================================================
c
      n=index(c_sounding_path,' ')-1
      write(*,*)'Data pathname: ',c_sounding_path(1:n)

      call get_sounding_data_cdf(c_sat_id,
     &                         c_sounding_path,
     &                         i4time_data,
     &                         c_filename_sat,
     &                         ires_x,ires_y, 
     &                         nelems,
     &                         nlines,
     &                         n_channels,
     &                         isndrdata,
     &                         wavelength,
     &                         scalingBias,
     &                         scalingGain,
     &                         nw_pix,nw_line,
     &                         se_pix,se_line,
     &                         ewCycles,
     &                         ewIncs,
     &                         nsCycles,
     &                         nsIncs,
     &                         f_time,
     &                         lineTimeBeg,lineTimeEnd,
     &                         imcI4,
     &                         orbitAttitude,
     &                         ndimx,ndimy,ndimch,
     &                         istatus)
      if(istatus.eq.1)then
         write(6,*)'Satellite data obtained'
         write(6,*)
      else
         write(6,*)'Data NOT obtained for ',c_sounding_path(1:n)
         goto 1000
      endif

      call get_i2_missing_data(i2_missing_data,iostatus)
      if(iostatus.ne.1)then
         write(6,*)'Error getting i2_missing_data'
         goto 1000
      endif

      write(6,*)'Set missing sat-sndr to i2_missing'
      write(6,*)'----------------------------------'
      do i=1,n_channels
         write(6,*)'Channel # ',i
         call set_missing_sndr(isndrdata(1,1,i),
     &               nelems,nlines,
     &               ndimx,ndimy,
     &               ismsng,
     &               i2_missing_data,
     &               mstatus)

      end do

      if(mstatus .ne. 0)then
         write(6,*)'Missing data found in isndrdata'
         write(6,*)'Number found: ',mstatus
      else
         write(6,*)'Found all isndrdata good in set_missing_sndr'
      endif
c
c =============================================================
c                   get laps domain lat/lon
c =============================================================
c
      call get_directory('static',dir_static,lend)
      var_ll(1) = 'LAT'
      var_ll(2) = 'LON'

      write(*,*)'Get LAPS lat/lon grid'
      call rd_laps_static(dir_static,'nest7grid',nx_l,ny_l,2,
     &var_ll,units_ll,comment_ll,data,grid_spacing,istatus)

      if(istatus.eq.1)then
         write(*,*)'LAPS lat/lon grid obtained'
         write(*,*)
         do j=1,ny_l
            do i=1,nx_l
               lat(i,j)=data(i,j,1)
               lon(i,j)=data(i,j,2)
            end do
         end do
      else
         write(*,*)'Unable to get lat/lon data'
         write(*,*)'sounding process terminating'
         stop
      end if
      grid_spacing_km=grid_spacing/1000.
c
c =================================================================
c       generate satellite-to-laps remapping table
c =================================================================
c
      write(6,*)'Compute sat-2-laps look-up-table'
      call gen_gvrsndr_lut_lsr(c_filename_sat,nlines,nelems,wavelength,
     &ires_x,ires_y,r_sndr_res_km,nw_pix,nw_line,se_pix,se_line,
     &ewCycles,ewIncs,nsCycles,nsIncs,f_time,orbitattitude,n_channels,
     &nx_l,ny_l,lat,lon,r_llij_lut_ri,r_llij_lut_rj,pct_req_lsr,istatus)

      if(istatus.eq.0)then
         write(6,*)'Sounder Nav computed'
      else
         write(6,*)'Found too many points outside of domain '
         goto 1000
      endif
c
c ===================================================================
c       compute count to radiance look up table
c ===================================================================
c
      call count_range(ndimx,ndimy,ndimch,nelems,nlines,n_channels,
     &isndrdata,imaximum,iminimum,istatus)

c     write(6,*)'Compute count2radiance lut'
c     call count2radiance_lut(n_channels,nlines,maxcnt,scalingBias,
c    &scalingGain,imaximum,iminimum,cnt2rad)
c
c
      write(6,*)'Image motion Comp ',imc
      write(6,*)'framsStartTime: ',f_time
      write(6,*)'ew cycle/inc: ',ewCycles,ewIncs
      write(6,*)'ns cycle/inc: ',nsCycles,nsIncs

      write(6,*)'netCDF file properly read'
      write(6,*)'nw pix/nw line :',nw_pix,nw_line
      write(6,*)
c
c ===================================================================
c               scale sounder counts to radiances.
c ===================================================================
c
      call get_r_missing_data(r_missing_data,iostatus)

      do k=1,n_channels

        do j=1,ndimy
        do i=1,ndimx(j)
           sndr_rad(i,j,k)=r_missing_data
        enddo
        enddo

        icnt(k) = 0
        jcnt(k) = 0
        rsb=scalingBias(1,k)
        rsg=scalingGain(1,k)

c       if(rsg.gt.1.0e-10)then                     !vis channel rsg is small. no rad conversion.
c                                                   this test doesn't work on linux-alpha (jet).

        if(k.lt.n_channels)then                    !the last channel is visible, so no radiance calc.

           write(6,*)'Scaling Bias/Gain ',k,rsb,rsg
           do j=1,ndimy
           do i=1,ndimx(j)

              rcount=float( isndrdata(i,j,k) )
              if(rcount.ge.imaximum(k).or.rcount.le.0.0)then  !this would include i2_missing_data (=-99)
                 icnt(k) = icnt(k) + 1
              else
                 sndr_rad(i,j,k)=(rcount-rsb)/rsg
                 jcnt(k) = jcnt(k) + 1
              endif

           enddo
           enddo

           print*,'Ch ',k,' # not used (gt imax or < 0): ',icnt(k)
c       print*,'Points Used: Sndrdata < imax and > 0: ',jcnt(k)
c       write(6,*)

        endif
      enddo
c 
c ================================================================
c          Determine representative time
c ================================================================
c
      rmintime=9999999999.
      rmaxtime=0.
      do j=1,ny_l
      do i=1,nx_l
         if(r_llij_lut_ri(i,j).ne.r_missing_data.and.
     &      r_llij_lut_rj(i,j).ne.r_missing_data)then
            ltindex=int(r_llij_lut_rj(i,j)+0.5)
            rltb=lineTimeBeg(ltindex,1)
            rlte=lineTimeEnd(ltindex,1)
            rmintime=min(rmintime,rltb)
            rmaxtime=max(rmaxtime,rlte)
         endif
      enddo
      enddo
      write(6,*)'max/min line times ',rmaxtime,rmintime
      i4time_data_orig=i4time_data
      if(rmaxtime.gt.0.0.and.rmintime.gt.0.0)then
         i4time_data=nint((rmaxtime+rmintime)/2.)+315619200
      elseif(rmaxtime.gt.0.0)then
         i4time_data=int(rmaxtime)+315619200
      else
         i4time_data=int(rmintime)+315619200
      endif
c
c ================================================================
c        Remap to laps domain
c ================================================================
c
      r_grid_ratio = r_sndr_res_km/grid_spacing_km
      write(6,*)'r grid ratio: ',r_grid_ratio

      do i = 1,n_channels

         call satdat2laps_sndr(nx_l,ny_l,
     &                  r_grid_ratio,
     &                  r_missing_data,
     &                  sndr_rad(1,1,i),
     &                  r_llij_lut_ri,
     &                  r_llij_lut_rj,
     &                  nlines, ! sndr array dimensions
     &                  nelems, !       "
     &                  sa,sc,st,
     &                  istatus)
     
         call move(sa,laps_data(1,1,i),nx_l,ny_l)

      enddo

      call write_lsr(nx_l,ny_l,n_channels,c_sat_id,rch_wvlngth,
     &laps_data,i4time_data,istatus)
      if(istatus.ne.1)then
         write(6,*)'Failed to write lsr'
      else
         call make_fnam_lp(i4time_data_orig,f9time,istatus)
         write(6,*)'Original Data Filetime: ',f9time
         call make_fnam_lp(i4time_data,f9time,istatus)
         write(6,*)'Computed Data Filetime: ',f9time
      endif

      goto 1000

901   write(6,*)'Error opening output file'
      goto 1000
902   write(6,*)'Error writing to output file'
      goto 1000
998   write(6,*)'Problem reading systime.dat'

1000  return
      end
