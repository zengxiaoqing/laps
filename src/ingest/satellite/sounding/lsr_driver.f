      Program lsr_driver
c
      include      'lapsparms.cmn'

      integer*4     max_sat
      integer*4     max_ch
      parameter     (max_sat=2,max_ch=19)
      integer*4     n_elems(max_sat)
      integer*4     n_lines(max_sat)
      integer*4     ismsng
      
      real*4        r_channel_wavelengths(max_ch,max_sat)
      character     c_sat_id(max_sat)*6
      character     c_sounding_path(max_sat)*200
      character     cmode*10

      Call get_laps_config('nest7grid',IStatus)
      if(istatus.eq.1)then
        write(*,*)'LAPS Parameters obtained'
      else
         write(*,*)'IStatus = ',IStatus,'Error - Get_LAPS_Config'
         write(*,*)'Terminating LAPS-lsr. Sat Sounder remapping'
         stop
      end if

      cmode='noinstall'
c
c Subroutine reads static/lsr/lsr.parms files
c
c     call readwrite_lsr_parms(cmode,nlines_max,nelems_max,n_channels,
c    &rdum,rdum,idum,idum,idum,idum,rdum,idum,idum,idum,
c    &rdum,idum,idum,idum,idum,istatus)
c     if(istatus.ne.1)then
c        write(6,*)'Error returned from readwrite_ln3_parms'
c        write(6,*)'Terminating'
c        goto 1000
c     endif
c
c get the number of satellites and channels from nest7grid.parms
c
      call get_sat_sounder_info(n_sat,c_sat_id,n_channels,
     &c_sounding_path,n_elems,n_lines,r_channel_wavelengths,
     &ismsng,istatus)

      if(istatus.ne.1)then
         write(6,*)'Error returned from get_sat_sounder_info'
         goto 1000
      endif

      do i=1,n_sat

         call lsr_driver_sub(NX_L_CMN,NY_L_CMN,n_channels,
     &     n_elems(i),n_lines(i),ismsng,c_sat_id(i),c_sounding_path(i),
     &     r_channel_wavelengths(i,1),istatus)
         if(istatus.ne.1)then
            if(i.eq.n_sat)then
               write(6,*)'No data processed in lsr_driver_sub'
               write(6,*)'Finished in lsr_driver'
            else
               write(6,*)'Try for another satellite'
            endif
         else
            write(6,*)'Data processed in lsr_driver_sub'
            goto 1000
         endif

      enddo

1000  stop
      end
c
c=======================================================
c
      subroutine lsr_driver_sub(nx_l,ny_l,n_channels,
     &nelems,nlines,ismsng,c_sat_id,c_sounding_path,rch_wvlngth,
     &istatus)
c
      implicit none
c
      integer*4     i_sat
      integer*4     nlines
      integer*4     nelems
      integer*4     nx_l,ny_l
      integer*4     n_channels

      integer*4     icnt(n_channels)
      integer*4     jcnt(n_channels)
      integer*4     icount
      integer*4     ismsng

      integer*4     ndimx(nlines)
      integer*4     ndimy,ndimch

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

      real*4        scalingBias(nlines,n_channels)
      real*4        scalingGain(nlines,n_channels)

      Integer*2     ewCycles,ewIncs
      Integer*2     nsCycles,nsIncs
      Integer*2     nw_pix,nw_line
      Integer*4     istatus
      Integer*4     iostatus
      Integer*4     mstatus
      Integer*4     i,j,k,n
      Integer*4     lend
      Integer*4     time_spec(2)
      Integer*4     imcI4
      Integer*4     int4(2)
      Integer*4     imax,jmax
      Integer*4     i4time_data
      Integer*4     i4time_data_orig
      Integer*4     imaximum(n_channels)
      Integer*4     iminimum(n_channels)
      Integer*4     i2_missing_data
      Integer*4     ltindex
c
      Integer*4     instr
      Integer*4     isndrdata(nelems,nlines,n_channels)
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
      Character     cmode*10
      Character     f9time*9
      Character     c_sat_id*6
      Character     cid*2
c =============================================================
c
c
c get sndr data
c
      n=index(c_sounding_path,' ')-1
      write(*,*)'Data pathname: ',c_sounding_path(1:n)

      call get_sounding_data_cdf(c_sat_id,
     &                         c_sounding_path,
     &                         i4time_data,
     &                         c_filename_sat,
     &                         nelems,
     &                         nlines,
     &                         n_channels,
     &                         isndrdata,
     &                         wavelength,
     &                         scalingBias,
     &                         scalingGain,
     &                         nw_pix,nw_line,
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
c compute count to radiance look up table
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
c get laps domain lat/lon
c
c Definitions needed for acquiring LAPS latitude and longitude arrays.
c -------------------------------------------------------------------
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
c generate sounding db to laps remapping lut
c
      write(6,*)'Compute sat-2-laps look-up-table'
      call gen_gv8sndr_lut_lsr(c_filename_sat,nlines,nelems,r_sndr_res_k
     &m,n_channels,nx_l,ny_l,lat,lon,r_llij_lut_ri,r_llij_lut_rj,
     &istatus)

      if(istatus.eq.0)then
         write(6,*)'Sounder Nav computed'
      else
         write(6,*)'Found points outside of domain '
      endif
c
c Next is to scale the sounder counts to radiances.
c
      call get_r_missing_data(r_missing_data,iostatus)

      do k=1,n_channels

        icnt(k) = 0
        jcnt(k) = 0
        rsb=scalingBias(1,k)
        rsg=scalingGain(1,k)

        write(6,*)'Scaling Bias/Gain ',k,rsb,rsg
        if(rsb.eq.0.or.rsg.eq.0.0)then
           write(*,*)'WARNING: rsg or rsg = 0.0'
           goto 875
        endif

        do j=1,ndimy
        do i=1,ndimx(j)

           rcount=float(isndrdata(i,j,k))
           if(rcount.ge.imaximum(k).or.rcount.le.0.0)then  !this would include i2_missing_data (=-99)
              icnt(k) = icnt(k) + 1
              sndr_rad(i,j,k)=r_missing_data
           else
              sndr_rad(i,j,k)=(rcount-rsb)/rsg
              jcnt(k) = jcnt(k) + 1
           endif

        enddo
        enddo

        write(6,*)'Channel ',k
        write(6,*)'Points Not Used: gt imax or < 0: ',icnt(k)
        write(6,*)'Points Used: Sndrdata < imax and > 0: ',jcnt(k)
        write(6,*)

875     enddo
c 
c Determine representative time 11-14-97 (J.Smart)
c
      i4time_data_orig=i4time_data
      rmintime=9999999999.
      rmaxtime=0.
      do j=1,ny_l
      do i=1,nx_l
         if(r_llij_lut_ri(i,j).ne.r_missing_data.and.
     &      r_llij_lut_rj(i,j).ne.r_missing_data)then
            ltindex=int(r_llij_lut_rj(i,j))
            rltb=lineTimeBeg(ltindex,1)
            rlte=lineTimeEnd(ltindex,1)
            rmintime=min(rmintime,rltb)
            rmaxtime=max(rmaxtime,rlte)
         endif
      enddo
      enddo

      write(6,*)'max/min line times ',rmaxtime,rmintime

      if(rmintime.le.0.or.rmaxtime.le.0)then
         if(rmintime.le.0.and.rmaxtime.gt.0)then
            i4time_data=((nint(rmaxtime)+315619200)+i4time_data_orig)/2.0
         elseif(rmintime.gt.0.and.rmaxtime.le.0)then
            i4time_data=((nint(rmintime)+315619200)+i4time_data_orig)/2.0
         endif
      else
         i4time_data=nint((rmaxtime+rmintime)/2.)+315619200
      endif
c
c remap to laps domain
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
     
         call load_array(sa,nx_l,ny_l,laps_data(1,1,i))
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
