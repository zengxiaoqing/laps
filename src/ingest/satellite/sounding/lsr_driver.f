      Program lsr_driver
c
      include      'lsr_dims.inc'

      
      real        r_channel_wavelengths(max_ch,max_sat)
      character     c_sat_id(max_sat)*6
      character     c_sounding_path(max_sat)*200

      integer       nx,ny
      integer       ismsng

      call get_grid_dim_xy(nx,ny,istatus)
      if (istatus .ne. 1) then
         write (6,*) 'Error getting horizontal domain dimensions'
         stop
      endif
c
c get the number of satellites and channels data/static/sat_sounder.nl
c
      call get_sat_sounder_info(n_sat,c_sat_id
     &,n_channels,c_sounding_path,r_channel_wavelengths
     &,ismsng,pct_req_lsr,istatus)

      if(istatus.ne.1)then
         print*,'Error returned from get_sat_sounder_info'
         goto 1000
      endif

      do i=1,n_sat

         call lsr_driver_sub(nx,ny,n_channels,
     &ismsng,c_sat_id(i),c_sounding_path(i),
     &r_channel_wavelengths(i,1),pct_req_lsr,istatus)

         print*
         if(istatus.ne.1)then
            if(i.eq.n_sat)then
               print*,'Data was not processed in lsr_driver_sub'
               print*,'Finished in lsr_driver'
            else
               print*,'Data was not processed in lsr_driver_sub'
               print*,'Try for another satellite'
            endif
         else
            print*,'Data was processed in lsr_driver_sub'
         endif
         print*

      enddo

      print*,' *** Finished in lsr_driver ***'

1000  stop
      end
c
c=======================================================
c
      subroutine lsr_driver_sub(nx_l,ny_l,n_channels,
     &ismsng,c_sat_id,c_sounding_path,rch_wvlngth,
     &pct_req_lsr,istatus)
c
      implicit none
c
      integer     i_sat
      integer     nlines
      integer     nelems
      integer     n_channels
      integer     nx_l,ny_l

      integer     icnt(n_channels)
      integer     jcnt(n_channels)
      integer     nradcnt(n_channels)
      integer     icount
      integer     ismsng

      integer     ndimx,ndimy,ndimch

      Real*8        orbitAttitude(336)
      Real*8        t, f_time

      Real*8,       allocatable :: lineTimeBeg(:,:)
      Real*8,       allocatable :: lineTimeEnd(:,:)

      real        lat(nx_l,ny_l)
      real        lon(nx_l,ny_l)
      real        r_llij_lut_ri(nx_l,ny_l)
      real        r_llij_lut_rj(nx_l,ny_l)
      real        sa(nx_l,ny_l)
      real        sc(nx_l,ny_l)
      real        st(nx_l,ny_l)
      real        laps_data(nx_l,ny_l,n_channels)
      real        grid_spacing
      real        grid_spacing_km
      real        r_sndr_res_km
      real        r_grid_ratio
      real        rcount,rsb,rsg
c
      real        data(nx_l,ny_l,2)
      real        rline(nx_l,ny_l)
      real        rpix(nx_l,ny_l)
      real        result
      real        xconv,yconv
      real        rch_wvlngth(n_channels)
      real        r_missing_data
      real        rmintime,rmaxtime
      real        rltb,rlte
      real        pct_req_lsr

      real,       allocatable :: sndr_rad(:,:,:)
      real,       allocatable :: scalingBias(:,:)
      real,       allocatable :: scalingGain(:,:)
      Integer,    allocatable :: isndrdata(:,:,:)

      Integer     ewCycles,ewIncs
      Integer     nsCycles,nsIncs
      Integer     nw_pix,nw_line
      Integer     se_pix,se_line
      Integer     istatus
      Integer     iostatus
      Integer     mstatus
      Integer     i,j,k,n,lf
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
      Integer     itstatus
      Integer     init_timer
      Integer     ishow_timer
c
      Integer       instr
      REAL*8        wavelength(n_channels)
      Character*1   imc(4)
      Character*255 c_filename_sat
      Character*255 datapath
      character*200 c_dataroot
      Character*200 c_sounding_path
      character*125 comment_ll(2)
      Character*100 filename
      character*10  c10_grid_fname
      character*10  units_ll(2)
      character*9   c_filetime_sat
      character*3   var_ll(2)
      character*2   cch
      character*150 dir_static
      Character     f9time*9
      Character     c_sat_id*6
      Character     cid*2
c
c =============================================================
c =============================================================
c      get sounder dimensions data for current sounder file
c =============================================================
c
      call find_domain_name(c_dataroot,c10_grid_fname,istatus)
      call s_len(c10_grid_fname,lf)

      n=index(c_sounding_path,' ')-1
      write(*,*)'Data pathname: ',c_sounding_path(1:n)

      call get_sounding_info_cdf(c_sat_id,
     &                         c_sounding_path,
     &                         i4time_data,
     &                         c_filename_sat,
     &                         ndimx,ndimy,ndimch,
     &                         istatus)
      if(istatus.eq.1)then
         write(6,*)'Satellite data obtained'
         write(6,*)
      else
         write(6,*)'Data NOT obtained for ',c_sounding_path(1:n)
         goto 1000
      endif

      call get_r_missing_data(r_missing_data,iostatus)
      call get_i2_missing_data(i2_missing_data,iostatus)
      if(iostatus.ne.1)then
         write(6,*)'Error getting missing_data flag'
         goto 1000
      endif
c
c =============================================================
c read sounder data
c =============================================================
c
      write(6,*)'Read sounder database '

      allocate (isndrdata(ndimx,ndimy,ndimch)
     &         ,linetimebeg(ndimy,ndimch)
     &         ,linetimeend(ndimy,ndimch)
     &         ,scalingBias(ndimy,ndimch)
     &         ,scalingGain(ndimy,ndimch) )

      itstatus=init_timer()
      itstatus=ishow_timer()

      isndrdata=0

      Call Read_sounder_db_cdf(c_filename_sat,
     &                         ndimx,ndimy,ndimch,
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
     &                         imc,ires_x,ires_y,
     &                         orbitAttitude,
     &                         istatus)

      if(istatus.ne.1)then
         print*,'Failed to read sounder data '
         return
      endif

      where (isndrdata .lt. 0) isndrdata=isndrdata+65535

      itstatus=ishow_timer()
      write(6,*)'Elapsed time (sec): ',itstatus
c
c !for sounder data the image motion compensation is off (ie.,=1)
c
      write(6,*)'Set missing sat-sndr to i2_missing'
      write(6,*)'----------------------------------'
      do i=1,n_channels
         write(6,*)'Channel # ',i
         call set_missing_sndr(isndrdata(1,1,i),
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
      call rd_laps_static(dir_static,c10_grid_fname,nx_l,ny_l,2,
     &var_ll,units_ll,comment_ll,data,grid_spacing,istatus)

      if(istatus.eq.1)then
         print*,'LAPS lat/lon grid obtained'
         lat(:,:)=data(:,:,1)
         lon(:,:)=data(:,:,2)
      else
         print*,'Unable to get lat/lon data'
         print*,'sounding process terminating'
         stop
      end if

      grid_spacing_km=grid_spacing/1000.
c
c =================================================================
c       generate satellite-to-laps remapping table
c =================================================================
c
      write(6,*)'Compute sat-2-laps look-up-table'
      call gen_gvrsndr_lut_lsr(c_filename_sat,ndimy,ndimx,wavelength,
     &ires_x,ires_y,r_sndr_res_km,nw_pix,nw_line,se_pix,se_line,
     &ewCycles,ewIncs,nsCycles,nsIncs,f_time,orbitattitude,ndimch,
     &nx_l,ny_l,lat,lon,r_llij_lut_ri,r_llij_lut_rj,pct_req_lsr,istatus)

      if(istatus.eq.1)then

         write(6,*)'Sounder Nav computed'

      else

         print*,'Found too many points outside of domain '

         deallocate (isndrdata
     &              ,lineTimeBeg
     &              ,lineTimeEnd
     &              ,scalingBias
     &              ,scalingGain)
         goto 1000

      endif
c
c ===================================================================
c       compute count to radiance look up table
c ===================================================================
c
      print*,'Allocate real sndr_rad array'
      allocate (sndr_rad(ndimx,ndimy,n_channels))

      call count_range(ndimx,ndimy,ndimch,isndrdata
     &,imaximum,iminimum,istatus)
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
c scale sounder counts to radiances. load channels desired for proc.
c ===================================================================
c
      print*,'Scale counts to radiance with rsb/rsg'
      print*,'---------------------------------------'

      do k=1,n_channels

       do j=1,ndimy
       do i=1,ndimx
          sndr_rad(i,j,k)=r_missing_data
       enddo
       enddo

       icnt(k) = 0
       jcnt(k) = 0
       nradcnt(k)=0
       rsb=scalingBias(1,k)
       rsg=scalingGain(1,k)

       if(rsg.gt.0.0)then 

        if(k.lt.n_channels)then          !the last channel is visible, so no radiance calc.

           print*,'Scaling Bias/Gain ',k,rsb,rsg

           do j=1,ndimy
           do i=1,ndimx

              rcount=float( isndrdata(i,j,k) )
              if(rcount.ge.imaximum(k).or.rcount.le.0.0)then  !this would include i2_missing_data (=-99)
                 icnt(k) = icnt(k) + 1
              else
                 sndr_rad(i,j,k)=(rcount-rsb)/rsg
                 jcnt(k) = jcnt(k) + 1
              endif

              if(sndr_rad(i,j,k).le. 0.0)then
c                print*,'rcount/rsb/rsg: ',rcount,rsb,rsg
                 nradcnt(k)=nradcnt(k)+1
              endif

           enddo
           enddo

           print*
           print*,'Ch ',k,' # not used (gt imax or < 0): ',icnt(k)
           if(nradcnt(k).gt.0)then
              print*,' Neg rad found: channel/# ',k,nradcnt(k)
           endif
c       write(6,*)

        endif

       else

        print*,' Scaling Gain = 0.'
        print*,' Radiance not computed for channel ',k

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
     &      r_llij_lut_rj(i,j).ne.r_missing_data.and.
     &      r_llij_lut_ri(i,j).gt.0.            .and.
     &      r_llij_lut_rj(i,j).gt.0.                )then
            ltindex=int(r_llij_lut_rj(i,j)+0.5)
            rltb=lineTimeBeg(ltindex,1)
            rlte=lineTimeEnd(ltindex,1)
            rmintime=min(rmintime,rltb)
            rmaxtime=max(rmaxtime,rlte)
         endif
      enddo
      enddo
      print*,' max/min line times ',rmaxtime,rmintime
      i4time_data_orig=i4time_data
      if(rmaxtime.gt.0.0.and.rmintime.gt.0.0)then
         print*,' Using max/min average for filetime'
         i4time_data=nint((rmaxtime+rmintime)/2.)+315619200
      elseif(rmaxtime.gt.0.0)then
         print*,' Using max line time for filetime'
         i4time_data=int(rmaxtime)+315619200
      elseif(rmintime.gt.0.0)then
         print*,' Using min line time for filetime'
         i4time_data=int(rmintime)+315619200
      else
         print*,' Using original i4time for filetime'
         i4time_data=i4time_data_orig
      endif

      deallocate (isndrdata
     &           ,lineTimeBeg
     &           ,lineTimeEnd
     &           ,scalingBias
     &           ,scalingGain)

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
     &                  ndimy,        ! sndr array dimensions
     &                  ndimx,        !       "
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
         if(istatus.ne.1)then
            print*,'Failed to make f9time original'
         else
            print*,'Original Data Filetime: ',f9time
         endif
         call make_fnam_lp(i4time_data,f9time,istatus)
         if(istatus.ne.1)then
            print*,'Failed to make f9time computed'
         else
            print*,'Computed Data Filetime: ',f9time
         endif
      endif

      deallocate (sndr_rad)

      goto 1000

901   write(6,*)'Error opening output file'
      goto 1000
902   write(6,*)'Error writing to output file'
      goto 1000
998   write(6,*)'Problem reading systime.dat'

1000  return
      end
