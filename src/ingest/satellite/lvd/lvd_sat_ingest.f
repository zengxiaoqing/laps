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
      Integer nsat,ntype

      include 'satellite_dims_lvd.inc'
      include 'satellite_common_lvd.inc'
      include 'satellite_namelist_lvd.cmn'

      character*3   chtype(maxchannel)
      character*9   cfname_cur
      character*9   cfname_sys
      character*6   csatid(maxsat)
      character*3   csattypes(maxtype*maxsat)
      character*3   cchanneltypes(maxchannel*maxtype*maxsat)
      character*200 cpath2sat(maxtype*maxsat)
      real solalt_thr_vis

      character   cgeneric_dataroot*255
      character   c_gridfname*50

      real, allocatable :: sri(:,:,:)
      real, allocatable :: srj(:,:,:)
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

      allocate (sri(nx_l,ny_l,maxchannel),
     &          srj(nx_l,ny_l,maxchannel))

c---------------------------------------------------------------
c Compute array dimensions for ir, vis, and wv.
c
      nsat=0
      do k=1,maxsat

       if(isats(k).eq.1)then
        nsat=nsat+1
        ntype=0

        write(6,*)
        write(6,*)' Searching for valid types for selected satellite '
     1            ,k,c_sat_id(k)

        do j=1,maxtype

         write(6,*)' itypes(j,k) = ',j,k,itypes(j,k)

         if(itypes(j,k).eq.1)then
          ntype=ntype+1
          nchannels=0
          nav_status=0
          do i=1,maxchannel
           if(ichannels(i,j,k).eq.1)then
            nchannels=nchannels+1
            chtype(nchannels)=c_channel_types(i,j,k)
           endif
          enddo

          if(nchannels.eq.0)then
            print*,'!!Error: No channels specified',
     &' for this satellite and type: ',c_sat_id(k),c_sat_types(j,k)
            print*,'!!Terminating!!'
            goto 1000
          endif
 
          print*,'=================================================='
          print*,'          lvd process information'
          print*,'--------------------------------------------------'
          print*,'Satellite ID: ',c_sat_id(k)
          print*,'Satellite TYPE: ',c_sat_types(j,k)
          write(6,40)(chtype(i),i=1,nchannels)
40        format(1x,'Satellite CHANNELS:',6(1x,a3))
          print*,'Path-to-raw-satellite'
          print*,'--------------------------------------------------'
          do i=1,maxchannel
            print*,' ',i,' ',TRIM(path_to_raw_sat(i,j,k))
          enddo

          if(trim(c_sat_types(j,k)) .ne. 'rll' .and.
     &       trim(c_sat_types(j,k)) .ne. 'cms'      )then
 
!           Obtain satellite sri/srj pixel coordinates for each model gridpoint
            call compute_nav_llij(nx_l,ny_l,maxchannel,nchannels,
     &c_sat_id(k),c_sat_types(j,k),chtype,k,j,cfname_cur,sri,srj,
     &nav_status)

            if(nav_status.eq.1)then
               print*,'Success computing mapping arrays ',c_sat_id(k)
     +               ,'/',c_sat_types(j,k)

            elseif(nav_status.lt.0)then
             print*,'ERROR returned from compute_nav_llij - stop'
             print*,'!!Terminating!!'
             goto 1000
            endif

          else
             print*,' Skip compute_nav_llij for type '
     &             ,trim(c_sat_types(j,k))
             print*,' Try reading lat/lon arrays from sat data file'
             print*,' for computing gri/j. For now sri and srj will'            
             print*,' be blank'            

          endif
c
c ================================================================
c
          call lvd_driver_sub(nx_l,ny_l,k,j,n_images,
     &                      chtype,i4time_cur,i_delta_sat_t_sec,
     &                      sri,srj,istatus)

          if(istatus.ne.1)then
            print*,'NO data processed by lvd_driver_sub: ',
     &c_sat_id(k),'/',c_sat_types(j,k)
          else
            print*,'Data processed by lvd_driver_sub: ',
     &c_sat_id(k),'/',c_sat_types(j,k)
          endif

c =================================================================
         endif 
        enddo
        if(ntype.eq.0)then
         print*,'!!Error: No type specified for this satellite:'
         print*,'   k,c_sat_id(k)',k,c_sat_id(k)
         print*,'!!Terminating!!'
         goto 1000
        endif
       endif ! isats(k) = 1
      enddo

      if(nsat.eq.0)then
       print*,'!!Error: No satellites specified for this lvd run'
       print*,'!!Check static/satellite_lvd.nl:  NSATS or CSATID'
       print*,'!!Terminating!!'
      endif
1000  deallocate (sri,srj)

      write(6,*)' Program lvd_sat_ingest complete...'

      stop
      end
