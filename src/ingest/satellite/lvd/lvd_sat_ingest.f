      Program lvd_sat_ingest
c
c	3-6-97	J. Smart	New main driver for laps satellite image ingest process.
c				Purpose for having this new top-level driver is to accomdate
c				dynamic memory requirements.
c				1. Include remapping LUT as subroutine.
c				2. Include acquisition of domain parameters from
c				   static/ nest7grid.parms.
c				3. lvd_driver now a subroutine called by this main routine.
c
c       9-12-97 J. Smart        Dynamic array development. Renamed this to sub(routine)
c       2-19-98 J. Smart        Incorporate satellite_master.nl. This eliminates the need for
c                               all the separate files containing nav info for each sat and
c                               each format type.
c                               Made this the main program once again.
c       12-28-98   "            Added 'include lapsparms.cmn' and call get_laps_config
c
      Implicit None

      Integer nx_l
      Integer ny_l
      Integer i,j,k
      Integer ispec
      Integer i4time_now_gg
      Integer i4time_cur
      Integer i4time_sys
      Integer istatus
      Integer laps_cycle_time
      Integer nav_status
      Integer nchannels

      include 'satellite_dims_lvd.inc'
      include 'satellite_common_lvd.inc'
      include 'satellite_namelist_lvd.cmn'

      character*3 chtype(maxchannel)
      character*9 cfname_cur
      character*9 cfname_sys
      character   cgeneric_dataroot*255
      character   c_gridfname*50

      real, allocatable :: gri(:,:,:)
      real, allocatable :: grj(:,:,:)
c
c ========================== START ==============================
c 
      call get_grid_dim_xy(nx_l,ny_l,istatus)
      if(istatus.ne.1)then
         print*,'Error returned from get_config'
         goto 1000
      endif

      call find_domain_name(cgeneric_dataroot,c_gridfname,istatus)
      write(6,*)'namelist parameters obtained: ',c_gridfname
c
      call config_satellite_lvd(istatus)
      if(istatus.ne.1)then
         write(*,*)'Error config_satellite - Cannot continue'
         stop
      endif

      i4time_cur = i4time_now_gg()
      call get_systime(i4time_sys,cfname_sys,istatus)
      call make_fnam_lp(i4time_cur,cfname_cur,istatus)

c this is designed to allow archive data runs!
      call get_laps_cycle_time(laps_cycle_time,istatus)
      if(i4time_cur-i4time_sys .gt. 2*laps_cycle_time)then
         print*,'Set current time to contents of systime.dat'
         cfname_cur=cfname_sys
         i4time_cur=i4time_sys
      endif

      allocate (gri(nx_l,ny_l,maxchannel),
     &          grj(nx_l,ny_l,maxchannel))

c---------------------------------------------------------------
c Compute array dimensions for ir, vis, and wv.
c
      do k=1,maxsat
       if(isats(k).eq.1)then

       do j=1,maxtype
        if(itypes(j,k).eq.1)then

         nav_status=0

50       nchannels=0
         do 4 i=1,maxchannel
          if(ichannels(i,j,k).eq.1)then
           nchannels=nchannels+1
           chtype(nchannels)=c_channel_types(i,j,k)
c--------------------------------------------------------------------
c         call lvd_file_specifier(c_channel_types(i,j,k),ispec,istatus)
c         if(istatus.ne.0)then
c            write(6,*)'Error status returned from lvd_file_specifier'
c            goto 1000
c         endif

c         if(ispec.eq.1)then
c            nelemvis=i_end_vis(j,k)-i_start_vis(j,k)+1
c            nlinesvis=j_end_vis(j,k)-j_start_vis(j,k)+1
c         elseif(ispec.eq.2.or.ispec.eq.4.or.ispec.eq.5)then
c            nelemir=i_end_ir(j,k)-i_start_ir(j,k)+1
c            nlinesir=j_end_ir(j,k)-j_start_ir(j,k)+1
c         elseif(ispec.eq.3)then
c            nelemwv=i_end_wv(j,k)-i_start_wv(j,k)+1
c            nlineswv=j_end_wv(j,k)-j_start_wv(j,k)+1
c         endif
c--------------------------------------------------------------------

         endif
4        enddo
 
         print*,'=================================================='
         print*,'          lvd process information'
         print*,'--------------------------------------------------'
         print*,'Satellite ID: ',c_sat_id(k)
         print*,'Satellite TYPE: ',c_sat_types(j,k)
         write(6,40)(chtype(i),i=1,nchannels)
40       format(1x,'Satellite CHANNELS:',5(1x,a3))
         print*,'Path-to-raw-satellite'
         print*,'--------------------------------------------------'
         do i=1,maxchannel
           print*,' ',i,' ',TRIM(path_to_raw_sat(i,j,k))
         enddo
 
c       if( (nlinesvis.eq.0.and.nelemvis.eq.0).and.
c    +      (nlinesir .eq.0.and.nelemir .eq.0).and.
c    +      (nlineswv .eq.0.and.nelemwv .eq.0) ) then
c            print*, 'All satellite array dimensions = 0 '
c            print*, 'Abort. Check static/satellite_lvd.nl'
c            stop
c       endif
c       if(nlinesir.eq.0.or.nelemir.eq.0)then
c            nlinesir=1
c            nelemir =1
c       endif
c       if(nlinesvis.eq.0.or.nelemvis.eq.0)then
c            nlinesvis=1
c            nelemvis =1
c       endif
c       if(nlineswv.eq.0.or.nelemwv.eq.0)then
c            nlineswv=1
c            nelemwv =1
c       endif

         call compute_nav_llij(nx_l,ny_l,maxchannel,nchannels,
     &c_sat_id(k),c_sat_types(j,k),chtype,k,j,cfname_cur,gri,grj,
     &nav_status)

         if(nav_status.eq.1)then
             print*,'Success computing mapping arrays ',c_sat_id(k)
     +,'/',c_sat_types(j,k)

c             call config_satellite_lvd(istatus)
c             goto 50

         elseif(nav_status.lt.0)then
              print*,'ERROR returned from compute_nav_llij - stop'
              goto 1000
         endif
c
c ================================================================
c
         call lvd_driver_sub(nx_l,ny_l,k,j,n_images,
     &                      chtype,i4time_cur,i_delta_sat_t_sec,
     &                      gri,grj,istatus)

         if(istatus.ne.1)then
            write(6,*)'NO data was processed by lvd_driver_sub'
         else
            write(6,*)'Data was processed by lvd_driver_sub'
         endif

c =================================================================
        endif 
       enddo
       endif
      enddo

1000  stop
      end
