cdis    Forecast Systems Laboratory
cdis    NOAA/OAR/ERL/FSL
cdis    325 Broadway
cdis    Boulder, CO     80303
cdis 
cdis    Forecast Research Division
cdis    Local Analysis and Prediction Branch
cdis    LAPS 
cdis 
cdis    This software and its documentation are in the public domain and 
cdis    are furnished "as is."  The United States government, its 
cdis    instrumentalities, officers, employees, and agents make no 
cdis    warranty, express or implied, as to the usefulness of the software 
cdis    and documentation for any purpose.  They assume no responsibility 
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis    
cdis    Permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  All modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making 
cdis    the modifications.  If significant modifications or enhancements 
cdis    are made to this software, the FSL Software Policy Manager  
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
      subroutine lvd_driver_sub(nx_l,ny_l,
     &               isat,jtype,nimages,
     &               n_ir_lines,n_ir_elem,
     &               n_wv_lines,n_wv_elem,
     &               n_vis_lines,n_vis_elem,
     &               chtype,maxchannels,nchannels,
     &               i4time_cur,
     &               lvd_status)
c
c Program drives generation of LAPS lvd.  Processes satellite data.
c
c The satillite data are acquired from several sources which are defined
c in /static/satellite_lvd.nl.
c
c More information on the lvd process is found in the README in this directory.
c
c Note: the channel types for this process (chtype) and the number of channels
c       nchannels will be less than or equal to (but never greater than) the
c       maxchannel variable within the include.
c
c Variables computed in subroutine that processes IR satellite data
c
c       ta4                               Channel 1 (3.9)  B-temps
c       tb4                                      "         B-temps averaged
c       ta6                               Channel 2 (6.7)  B-temps
c       tb6                                      "         B-temps averaged       
c       ta8              RA       O       Channel 4 (11.2) Brightness temps (averaged)
c       tb8              RA       O       Channel 4 (11.2) Brightness temps (warm pixel)
c       tc8              RA       O           "            Brightness temps (filtered)
c       ta12                              Channel 5 (12.0) B-temps
c       tb12                                     "          "       averaged
c       laps_vis         RA       O       Visible (raw)
c       vis_norm         RA       O       Visible (normalized)
c       albedo           RA       O       Albedo (0.0 -- 1.0)
c
c
      implicit none

      include 'satellite_dims_lvd.inc'  ! dimensions for maxsat, maxtype, maxchannel
c                                     used for variable in include 'satellite_common_lvd.inc'

      integer   nx_l,ny_l
      integer   n_ir_elem,n_ir_lines
      integer   n_vis_elem,n_vis_lines
      integer   n_wv_elem,n_wv_lines
      integer   maxchannels
      integer   nchannels
      integer   nimages
      integer   isat,jtype

      integer   n_lvd_fields_max
      parameter (n_lvd_fields_max = 14)
      integer   max_files
      parameter (max_files=1500)
 
      integer   n1,n2,nn,ns
      integer   lend,len_lvd

      character*3 csattype
      character*3 c_type(maxchannels,max_files)
      character*3 csat_type
      character*3 chtype(maxchannels)
      character*6 csatid
      character*9 c_fname_cur
      character*9 c_fname
      character*10 cmode
      character*255 cname

      real image_vis (n_vis_elem,n_vis_lines,nimages) 
      real image_ir  (n_ir_elem,n_ir_lines,nimages)
      real image_12  (n_ir_elem,n_ir_lines,nimages)
      real image_39  (n_ir_elem,n_ir_lines,nimages)
      real image_67  (n_wv_elem,n_wv_lines,nimages)

      real laps_data(nx_l,ny_l,n_lvd_fields_max)
      real visnorm(nx_l,ny_l)
      real visraw(nx_l,ny_l)
      real albedo(nx_l,ny_l)
      real ta8(nx_l,ny_l)
      real tb8(nx_l,ny_l)
      real tc8(nx_l,ny_l)    ! 8 is the 11.2 (referred to as ir here) micron data
      real ta12(nx_l,ny_l)
      real tb12(nx_l,ny_l)   ! 12 is the 12.0 (referred to as 12 here) micron data
      real ta4(nx_l,ny_l)
      real tb4(nx_l,ny_l)    ! 4 is the 3.9 (referred to as 39 here) micron data
      real ta6(nx_l,ny_l)
      real tb6(nx_l,ny_l)    ! 6 is the 6.7 (referred to as 67 here) micron data
      real lat(nx_l,ny_l)
      real lon(nx_l,ny_l)
      real topo(nx_l,ny_l)
      real r_grid_ratio(maxchannels,nimages)

      real	range_m
      real      sublon_d
      real      sublat_d
      real      smsng
      real      r_missing_data
      real      radtodeg

      logical   lvis_flag
      logical   lsatqc
      logical   l_lut_flag

      integer	iskip_bilin
      integer   i,j,k,l
      integer   n,nf,np
      integer   ispec
      integer   nlf
      integer   nlf_prev
      integer   in
c
      real*4      vis_cnt_to_cnt_lut(0:1023)
      real*4      ir_cnt_to_btemp_lut(0:1023) !this one is 11u
      real*4      r12_cnt_to_btemp_lut(0:1023)
      real*4      r39_cnt_to_btemp_lut(0:1023)
      real*4      r67_cnt_to_btemp_lut(0:1023)
      real*4      r_llij_lut_ri(nx_l,ny_l,maxchannels)
      real*4      r_llij_lut_rj(nx_l,ny_l,maxchannels)
      real*4      good_vis_data_thresh

      character*100 path
c
c dimensions for lat/lon
c
      character*125 comment_ll(2)
      character*10 units_ll(2)
      character*3 var_ll(2)
c
c dimensions for lvd
c
      character*125 c_lvd(n_lvd_fields_max)
      character*10 units_lvd(n_lvd_fields_max)
      character*3 var_lvd(n_lvd_fields_max)
      character*4 lvl_coord_lvd(n_lvd_fields_max)
      integer lvl_lvd(n_lvd_fields_max)
      character*150 dir_lvd
      character*31 ext_lvd

      real*4 grid_spacing_laps_m
      real*4 r_image_res_m(maxchannels,nimages)
      real*4 r_image_status(maxchannels,max_files)

      integer ishow_timer
      integer init_timer
      integer i4time_cur
      integer i4time_data(max_files)
      integer istat
      integer gstatus
      integer istatus
      integer istatus_vis(3)
      integer itstatus
      integer lvd_status
      integer nft,ntm(max_files)

      include 'satellite_common_lvd.inc'
c =========================================================================
c ----------------------------- START -------------------------------------
c
      itstatus=init_timer()
      itstatus=ishow_timer()

      csatid = c_sat_id(isat)
      csattype = c_sat_types(jtype,isat)
      smsng = float(i_msng_sat_flag(jtype,isat))
c ----------------------------------------------------------------------
c if current time is at beginning of new day, then adjust time back just
c a few seconds to allow any data just before top of hour to have a chance
c at being processed now.
c ---------------------------------------------------------------------
      call make_fnam_lp(i4time_cur,c_fname_cur,istatus)
      if(c_fname_cur(6:9).eq.'0000')then
         i4time_cur=i4time_cur-15
         call make_fnam_lp(i4time_cur,c_fname_cur,istatus)
      endif

      write(6,*)'Current LVD process time: ',
     &'filename: ',c_fname_cur,' i4time: ',i4time_cur
c ---------------------------------------------
c acquiring LAPS latitude and longitude arrays.
c ---------------------------------------------
      call get_domain_laps(nx_l,ny_l,'nest7grid',lat,lon,topo,
     &grid_spacing_laps_m,istatus)
      if(istatus.eq.1)then
         write(6,*)'LAPS lat/lon/grid_spacing obtained'
         write(6,*)
      else
         write(6,*)'Error getting LAPS lat/lon data'
         stop
      end if
c --------------------------------------------------------------------
c Read look-up table for mapping lat/lon data pixels to real i/j pairs
c --------------------------------------------------------------------
      if(csatid.ne.'gmssat')then

         call readlut(csatid,csattype,maxchannels,nchannels,
     &chtype,nx_l,ny_l,r_llij_lut_ri,r_llij_lut_rj,istatus)

      if(istatus.eq.1)then

         write(6,*)'LUT not obtained: ',csatid,'/',csattype

c           write(6,*)'Computing lut using genlvdlut_lvd'
c        call genlvdlut_sub(nx_l,ny_l,gstatus)
c           call genlvdlut_lvd(nx_l,ny_l,lat,lon,jtype,isat,
c    +gstatus)
c           if(gstatus.lt.0)then
c              write(6,*)'Error generating LUT - terminating'
c              goto 910
c           else
c              write(6,*)'rewrite satellite_lvd.nl'
c              write(6,*)
c              call rewrite_satellite_lvd_nl(istatus)
c              call readlut(csatid,csattype,maxchannels,nchannels,
c    &chtype,nx_l,ny_l,r_llij_lut_ri,r_llij_lut_rj,istatus)
c              if(istatus.eq.1)then
c                 write(6,*)'Error reading new luts - terminating'
c                 goto 909
c              endif
c           endif

      else
         write(6,*)'Got the mapping look-up-tables '

c           write(6,*)'Check if luts are up-to-date'

c           call check_luts(c_fname_cur,isat,jtype,
c    &chtype,maxchannels,nchannels,l_lut_flag,istatus)

c           if(l_lut_flag.and.istatus.eq.0)then
c              print*,'*******************************************'
c              write(6,*)'Found difference in nav parms',
c    +' - rebuild the lut'
c              print*,'*******************************************'
c              call genlvdlut_lvd(nx_l,ny_l,lat,lon,jtype,isat,
c    +gstatus)
c              if(gstatus.lt.0)then
c                 write(6,*)'Error generating LUT - terminating'
c                 goto 910
c              else
c                 write(6,*)'**********************************'
c                 write(6,*)
c                 call rewrite_satellite_lvd_nl(istatus)
c                 call readlut(csatid,csattype,maxchannels,nchannels,
c    &chtype,nx_l,ny_l,r_llij_lut_ri,r_llij_lut_rj,istatus)
c                 if(istatus.eq.1)then
c                    print*,'Error reading new luts - terminating'
c                    goto 909
c                 endif
c              endif
c           elseif(istatus.eq.0)then
c              write(6,*)'Lut checked out ok'
c              write(6,*)
c           else
c              write(6,*)'Error status returned from check_lut'
c              goto 910
c           endif

      endif
      endif
c
c -------------------------------------------------------------------------
c Determine solar-altitude and set flag for visible sat data availability.
c This flag will force  this process to not wait for vis data at the end of
c the day when the sun is setting.
c
      lvis_flag = .false.   !assume that it is available
      call set_vis_flag(i4time_cur,lat,lon,nx_l,ny_l,lvis_flag)
      if(lvis_flag)then
         write(6,*)'lvis_Flag set indicating no visible data'
      endif
c
c --------------------------------------------------------------------------
c Find and read current satellite files... as many as 4 ir channels and vis.
c --------------------------------------------------------------------------
      if(csattype.eq.'cdf'.or.csattype.eq.'gvr'.or.
     &   csattype.eq.'wfo')then

       write(6,*)'Using getcdf_satdat routine'

         call getcdf_satdat(csatid,
     &                      csattype,
     &                      nchannels,chtype,
     &                      path_to_raw_sat(1,jtype,isat),
     &                      c_fname_cur,lvis_flag,
     &                      i_delta_sat_t_sec,
     &                      n_ir_lines, n_ir_elem,
     &                      n_vis_lines,n_vis_elem,
     &                      n_wv_lines,n_wv_elem,
     &                      maxchannels,nimages,
     &                      nft,ntm,c_type,max_files,
     &                      image_ir,image_vis,
     &                      image_12,image_39,image_67,
     &                      i4time_data,
     &                      istatus)

         if(istatus .ne. 1)then
            write(6,*)'Did not get data for ',c_fname_cur
            goto 998
         endif

      elseif(csattype.eq.'asc')then   !then we are using ascii files for raw ingest sat data

         write(6,*)'datapath: ',path_to_raw_sat(1,jtype,isat)(1:in)
         write(6,*)'Using getascii_satdat routine'

c  only possible to have one time for ascii files (nft=1); however, the number of
c  matches for this time (ntm) >= 0 depending on the result in getascii_satdat.

         nft=1
         write(6,*)'ascii satellite data ingest currently disabled'

c        call getascii_satdat(i4time_cur,lvis_flag,
c    &                        nchannels,chtype,
c    &                        n_ir_lines, n_ir_elem,
c    &                        n_vis_lines,n_vis_elem,
c    &                        n_wv_lines,n_wv_elem,
c    &                        path_to_raw_sat(1,jtype,isat), 
c    &                        ntm(nft),c_type(1,1),maxchannel,
c    &                        image_ir(1,1,1),
c    &                        image_vis(1,1,1),
c    &                        image_12(1,1,1),
c    &                        image_39(1,1,1),
c    &                        image_67(1,1,1),
c    &                        i4time_data(nft),
c    &                        r_image_res_m(1,1),
c    &                        istatus)

c        if(istatus .eq. 1)then
c           write(6,*)'Did not get data for ',c_fname_cur
c           goto 998
c        end if
c
c 5-15-97: JSmart added gwc satdat switch
c
      elseif(csattype.eq.'gwc')then

         write(6,*)'Search for data: ',csattype
         write(6,*)'Using getafgwc_satdat routine'

         nft=1
         call getafgwc_satdat(isat,jtype,
     &                        maxchannels,nchannels,chtype,
     &                        i4time_cur,lvis_flag,
     &                        n_ir_lines,n_ir_elem,
     &                        n_vis_lines,n_vis_elem,
     &                        n_wv_lines,n_wv_elem,
     &                        ntm(nft),max_files,c_type,
     &                        image_ir,image_vis,
     &                        image_12,image_39,image_67,
     &                        i4time_data(nft),
     &                        istatus)

         if(istatus.ne.0)then
            write(6,*)'Did not get data for ',c_fname_cur
            goto 998
         endif 

      endif
c --------------------------------------------------------------------
c Get image resolution information
c --------------------------------
      do j = 1,nft
      do i = 1,ntm(j)

         call lvd_file_specifier(c_type(i,j),ispec,istatus)
         goto(301,302,303,302,302)ispec

301           r_image_res_m(i,j)=r_resolution_x_vis(jtype,isat)
              goto 304
302           r_image_res_m(i,j)=r_resolution_x_ir(jtype,isat)
              goto 304
303           r_image_res_m(i,j)=r_resolution_x_wv(jtype,isat)

304      continue

      enddo
      enddo
c
c -----------------------------------------------------------------------
c Compute or read ir/vis count to brightness temp (Tb)/vis count-to-count.
c -----------------------------------------------------------------------
c
      if(csattype.eq.'cdf'.or.csattype.eq.'wfo')then

         write(6,*)'Compute ',csatid,' cnt-to-btemp lookup tables'
         call genbtemplut(2,r39_cnt_to_btemp_lut,istatus)
         if(istatus.ne.1)then
            write(6,*)'Error computing 39 lut'
         endif
         call genbtemplut(3,r67_cnt_to_btemp_lut,istatus)
         if(istatus.ne.1)then
            write(6,*)'Error computing 67 lut'
         endif
         call genbtemplut(4,ir_cnt_to_btemp_lut,istatus)
         if(istatus.ne.1)then
            write(6,*)'Error computing ir lut'
         endif
         call genbtemplut(5,r12_cnt_to_btemp_lut,istatus)
         if(istatus.ne.1)then
            write(6,*)'Error computing 12 lut'
         endif

      elseif(csattype.eq.'gvr'
     &.or.csattype.eq.'gwc')then

         write(6,*)'get gvarimage count to btemp lut'
         write(6,*)'and gvarimage vis count-to-count lut'

         do j = 1,nft
         do i = 1,ntm(j)

         call lvd_file_specifier(c_type(i,j),ispec,istat)

         goto(1,2,3,4,5)ispec     !ispec = 1 is visible data

1           if(csatid.eq.'gmssat')goto 6    

            call read_gvarimg_cnt2btemp_lut(csatid,
     &c_type(i,j),vis_cnt_to_cnt_lut,istatus)
            goto 6

2           call read_gvarimg_cnt2btemp_lut(csatid,
     &c_type(i,j),r39_cnt_to_btemp_lut,istatus)
            goto 6

3           call read_gvarimg_cnt2btemp_lut(csatid,
     &c_type(i,j),r67_cnt_to_btemp_lut,istatus)
            goto 6

4           call read_gvarimg_cnt2btemp_lut(csatid,
     &c_type(i,j),ir_cnt_to_btemp_lut,istatus)
            goto 6

5           call read_gvarimg_cnt2btemp_lut(csatid,
     &c_type(i,j),r12_cnt_to_btemp_lut,istatus)

6        enddo
         enddo
c
c           call count2radiance_lut(n_vis_lines,scalingBias,
c    &scalingGain,cnt2rad(1,i))

      endif
c
c ------------------------------------------------------------
c check for and fill-in for any missing data in current images
c ------------------------------------------------------------
c
      lsatqc=.true.
      if(csattype.eq.'asc'.or.csattype.eq.'gwc')lsatqc=.false.
      
      if(lsatqc)then

         write(6,*)'Entering satdatfill routine'

         call satdatfill(csatid,csattype,nft,ntm,
     &   n_ir_elem,n_ir_lines,n_vis_elem,n_vis_lines,
     &   n_wv_elem,n_wv_lines,c_type,smsng,
     &   maxchannels,nimages,max_files,
     &   i4time_data,path_to_raw_sat(1,jtype,isat),
     &   image_ir,image_39,image_12,image_67,image_vis,
     &   r_image_status)

      else         !  if(csatid.ne.'gmssat')then

         write(6,*)'TEST = skipping satfill section'
         write(6,*)'==============================='

         write(6,*)'Not using fill routine '
         write(6,*)'Only use set_missing_flag '

         call set_missing_flag(csatid,csattype,n_ir_elem,n_ir_lines,
     &             n_vis_elem,n_vis_lines,n_wv_elem,n_wv_lines,
     &             nft,ntm,c_type,smsng,maxchannels,nimages,
     &             image_ir,image_39,image_12,image_67,image_vis,
     &             r_image_status)
         write(6,*)

      endif

      iskip_bilin = 1
c
c ---------------------------------------------------
c satellite range and sub-longitude (namelist items).
c ---------------------------------------------------
      radtodeg=180.0/acos(-1.)
      range_m = sat_range_m(isat)
      sublon_d = r_sat_sub_lon(isat)*radtodeg
      sublat_d = r_sat_sub_lat(isat)*radtodeg
      write(6,*)'range_m = ',range_m
      write(6,*)'sublat_d = ',sublat_d
      write(6,*)'sublon_d = ',sublon_d
      write(6,*)'LAPS grid spacing (m): ',grid_spacing_laps_m
      write(6,*)
c
c -------------------------------------------------
c output info concerning the state of these images.
c -------------------------------------------------
      do i = 1,nft
      do j = 1,ntm(i)

         write(6,*)'Satellite Channel: ',c_type(j,i)

         if(r_image_status(j,i) .gt. 0.0)then
            write(6,*)'Some Bad data found = ', c_type(j,i)
     &,r_image_status(j,i)
         else
            write(6,*)'Good data type/status: ',c_type(j,i)
     &,r_image_status(j,i)
         endif
c
c compute grid ratio input/output resolutions
c
      r_grid_ratio(j,i)=r_image_res_m(j,i)/grid_spacing_laps_m

         write(6,*)'Image resolution (m): ',r_image_res_m(j,i)
         write(6,*)'I/O Grid ratio: ',r_grid_ratio(j,i)
         write(6,*)

      enddo
      enddo
c
      call get_r_missing_data(r_missing_data,istatus)
      if(istatus.ne.1)then
         write(6,*)'error getting r_missing_data'
         goto 16
      endif
c
c ------------------------------------------------------------
c convert from counts to brightness temps for CDF data use the
c pre-generated lut's. For ascii data divide all by 10.
c ------------------------------------------------------------
       if(csattype.ne.'asc')then
          write(6,*)
          write(6,*)'Convert counts to brightness temps'
          do i = 1,nft
          do j = 1,ntm(i)
             if(r_image_status(j,i).le.0.333)then

             call lvd_file_specifier(c_type(j,i),ispec,istat)
             goto(15,11,12,13,14)ispec

13           write(6,*)'For ',c_type(j,i)
             call btemp_convert(n_ir_elem,n_ir_lines,
     &                        ir_cnt_to_btemp_lut,
     &                        r_missing_data,
     &                        image_ir(1,1,i))
             goto 17
 
11           write(6,*)'For ',c_type(j,i)
             call btemp_convert(n_ir_elem,n_ir_lines,
     &                        r39_cnt_to_btemp_lut,
     &                        r_missing_data,
     &                        image_39(1,1,i))
             goto 17

12          write(6,*)'For ',c_type(j,i)
             call btemp_convert(n_wv_elem,n_wv_lines,
     &                        r67_cnt_to_btemp_lut,
     &                        r_missing_data,
     &                        image_67(1,1,i))
             goto 17

14           write(6,*)'For ',c_type(j,i)
             call btemp_convert(n_ir_elem,n_ir_lines,
     &                        r12_cnt_to_btemp_lut,
     &                        r_missing_data,
     &                        image_12(1,1,i))
             endif

             goto 17

15           if(csattype.eq.'gvr'.or.csattype.eq.'gwc')then
                if(csatid.eq.'gmssat')goto 17
                call btemp_convert(n_vis_elem,n_vis_lines,
     &                          vis_cnt_to_cnt_lut,
     &                          r_missing_data,
     &                          image_vis(1,1,i))
              write(*,*)'Converted GVAR vis data using cnt-to-cnt lut'
             else
              write(*,*)'Not converting ',csattype,' vis data'
             endif
17        enddo
          enddo
          write(6,*)'Done with conversions '
          write(6,*)

       elseif(csattype.eq.'asc')then
c
c note that brightness temps in the ascii file have 1 significant
c decimal digit and have been acquired as integers so convert here.
c 
          write(6,*)'Convert ascii btemps from Integer to Floating pt.'
          do i = 1,nft
          do j = 1,ntm(i)
             call lvd_file_specifier(c_type(j,i),ispec,istat)
             goto(26,22,23,24,25)ispec

22           call btemp_convert_asc(n_ir_lines,n_ir_elem,
     &               r_missing_data,image_39(1,1,i),istatus)
             goto 26

23           call btemp_convert_asc(n_wv_lines,n_wv_elem,
     &               r_missing_data,image_67(1,1,i),istatus)
             goto 26

24           call btemp_convert_asc(n_ir_lines,n_ir_elem,
     &               r_missing_data,image_ir(1,1,i),istatus)
             goto 26

25           call btemp_convert_asc(n_ir_lines,n_ir_elem,
     &               r_missing_data,image_12(1,1,i),istatus)

26        enddo
          enddo
          write(6,*)'Done with conversion'
          write(6,*)
       else
          write(6,*)'Cannot convert image to btemps yet'
       endif
c
c ------------------------------------------------------------------------------------------------
c Start processing satellite image data found for the current run time.
c nft = number of file times. ntm(i) is the number of time matches for each file time
c nft can be > 1 if, 1) the code has been idle; 2) if there is more
c than one time due to rapid scan; 3) there are file times within threshold i_sat_delta_t_sec
c (found in static/nest7grid.parms). nft must never exceed parameter max_images (satellite_lvd.nl).
c This is insured within the getcdf, getafgwc, etc code.
c ------------------------------------------------------------------------------------------------
      write(6,*)
      write(6,*)'Ready to remap satellite data'
      write(6,*)'-----------------------------'
c
c This for output.  LAPS LVD files as indicated.
c
      call get_directory('lvd',dir_lvd,len_lvd)
      dir_lvd=dir_lvd(1:len_lvd)//csatid//'/'
      len_lvd=index(dir_lvd,' ')-1
      ext_lvd = 'lvd'
c
c initialize output array
c
      do i=1,n_lvd_fields_max
         lvl_lvd(i) = 0
         lvl_coord_lvd(i) = 'AGL'
         units_lvd(i)='K'
         do k=1,ny_l
         do j=1,nx_l
            laps_data(j,k,i)=r_missing_data
         enddo
         enddo
      enddo

      do i = 1,nft

         nlf = 0
         nlf_prev = 1

         call make_fnam_lp(i4time_data(i),c_fname,istatus)

c
c ----------  GMS SATELLITE SWITCH -------
         if(csatid.eq.'gmssat')goto 310


         do j = 1,ntm(i)

            n=index(c_type(j,i),' ')-1
            if(n.le.0)n=3
            call lvd_file_specifier(c_type(j,i),ispec,istat)

            if(ispec.eq.4)then
            if(r_image_status(j,i).le.0.3333)then

               call process_ir_satellite(i4time_data(i),
     &                      nx_l,ny_l,lat,lon,
     &                      n_ir_lines,n_ir_elem,
     &                      r_grid_ratio(j,i),
     &                      image_ir(1,1,i),
     &                      r_llij_lut_ri(1,1,ispec),
     &                      r_llij_lut_rj(1,1,ispec),
     &                      c_type(j,i),
     &                      ta8,tb8,tc8,
     &                      istatus)

               if(istatus .ne. 1)then
                  write(*,*)'Error processing IR Satillite Data'
               else
                  nlf=nlf+1
                  call move(ta8,laps_data(1,1,nlf),nx_l,ny_l)
                  var_lvd(nlf)  = 'S8A'       ! satellite, channel-4, averaged
                  c_lvd(nlf)=csatid//' (11.2u) IR B-TEMPS - AVERAGED'
                  nlf=nlf+1
                  call move(tb8,laps_data(1,1,nlf),nx_l,ny_l)
                  var_lvd(nlf)='S8W'       ! satellite, channel-4, warm pixel
                  c_lvd(nlf)=csatid//' (11.2u) IR B-TEMPS; WARM PIX'
                  nlf=nlf+1
                  call move(tc8,laps_data(1,1,nlf),nx_l,ny_l)
                  var_lvd(nlf)='S8C'       ! satellite, channel-4, warm pixel
                  c_lvd(nlf)=csatid//' (11.2u) IR B-TEMPS - FILTERED'
               endif

            else
               write(6,*)'IR image not processed: missing ir data'
            endif

            elseif(ispec.eq.2)then
            if(r_image_status(j,i).lt.0.3333)then

               call process_ir_satellite(i4time_data(i),
     &                      nx_l,ny_l,lat,lon,
     &                      n_ir_lines,n_ir_elem,
     &                      r_grid_ratio(j,i),
     &                      image_39(1,1,i),
     &                      r_llij_lut_ri(1,1,ispec),
     &                      r_llij_lut_rj(1,1,ispec),
     &                      c_type(j,i),
     &                      ta4,tb4,tb4,
     &                      istatus)

               if(istatus .ne. 1)then
                  write(*,*)'Error processing IR Satillite Data'
               else
                  nlf=nlf+1
                  call move(ta4,laps_data(1,1,nlf),nx_l,ny_l)
                  var_lvd(nlf) = 'S3A'       ! satellite, , averaged
                  c_lvd(nlf)=csatid//' (3.9u) IR B-TEMPS - AVERAGED'
                  nlf=nlf+1
                  call move(tb4,laps_data(1,1,nlf),nx_l,ny_l)
                  var_lvd(nlf)  = 'S3C'       ! satellite, , filtered
                  c_lvd(nlf)=csatid//' (3.9u) IR B-TEMPS - FILTERED'
               end if
            else
                   write(6,*)'39u image not processed: missing data'
            endif
            elseif(ispec.eq.5)then
            if(r_image_status(j,i).lt.0.3333)then

               call process_ir_satellite(i4time_data(i),
     &                      nx_l,ny_l,lat,lon,
     &                      n_ir_lines,n_ir_elem,
     &                      r_grid_ratio(j,i),
     &                      image_12(1,1,i),
     &                      r_llij_lut_ri(1,1,ispec),
     &                      r_llij_lut_rj(1,1,ispec),
     &                      c_type(j,i),
     &                      ta12,tb12,tb12,
     &                      istatus)

               if(istatus .ne. 1)then
                  write(*,*)'Error processing IR Satillite Data'
               else
                  nlf=nlf+1
                  call move(ta12,laps_data(1,1,nlf),nx_l,ny_l)
                  var_lvd(nlf) = 'SCA'       ! satellite, averaged
                  c_lvd(nlf)=csatid//' (12.0u) IR B-TEMPS - AVERAGED'
                  nlf=nlf+1
                  call move(tb12,laps_data(1,1,nlf),nx_l,ny_l)
                  var_lvd(nlf) = 'SCC'       ! satellite, averaged
                  c_lvd(nlf)=csatid//' (12.0u) IR B-TEMPS - FILTERED'
               endif
            else
                   write(6,*)'12u image not processed: missing data'
            endif
            elseif(ispec.eq.3)then
            if(r_image_status(j,i).lt.0.3333)then

               call process_ir_satellite(i4time_data(i),
     &                      nx_l,ny_l,lat,lon,
     &                      n_wv_lines,n_wv_elem,
     &                      r_grid_ratio(j,i),
     &                      image_67(1,1,i),
     &                      r_llij_lut_ri(1,1,ispec),
     &                      r_llij_lut_rj(1,1,ispec),
     &                      c_type(j,i),
     &                      ta6,tb6,tb6,
     &                      istatus)

               if(istatus .ne. 1)then
                  write(*,*)'Error processing wv Satillite Data'
               else
                  nlf=nlf+1
                  call move(ta6,laps_data(1,1,nlf),nx_l,ny_l)
                  var_lvd(nlf) = 'S4A'       ! satellite, averaged
                  c_lvd(nlf)=csatid//' (6.7u) IR B-TEMPS - AVERAGED'
                  nlf=nlf+1
                  call move(tb6,laps_data(1,1,nlf),nx_l,ny_l)
                  var_lvd(nlf) = 'S4C'       ! satellite, filtered
                  c_lvd(nlf)=csatid//' (6.7u) IR B-TEMPS - FILTERED'
               end if
            else
                  write(6,*)'wv image not processed: missing data'
            endif
            elseif(ispec.eq.1)then
            if(r_image_status(j,i).lt.0.3333)then

               call process_vis_satellite(csatid,
     &                      csattype,
     &                      i4time_data(i),
     &                      nx_l,ny_l,lat,lon,
     &                      n_vis_lines,n_vis_elem,    !array dimensions
     &                      r_grid_ratio(j,i),
     &                      image_vis(1,1,i),
     &                      r_llij_lut_ri(1,1,ispec),
     &                      r_llij_lut_rj(1,1,ispec),
     &                      sublat_d,sublon_d,range_m,
     &                      l_national,iskip_bilin,
     &                      visraw,visnorm,albedo,
     &                      istatus_vis)

               good_vis_data_thresh=(nx_l*ny_l)*0.3333

               if(abs(istatus_vis(1)).lt.good_vis_data_thresh)then
                  nlf=nlf+1
                  call move(visraw,laps_data(1,1,nlf),nx_l,ny_l)
                  var_lvd(nlf) = 'SVS'       ! satellite, visible
                  c_lvd(nlf)=csatid//' (VISIBLE) SATELLITE - RAW'
                  units_lvd(nlf) = 'COUNTS'
               endif

               if(abs(istatus_vis(2)).lt.good_vis_data_thresh)then
                  nlf=nlf+1
                  call move(visnorm,laps_data(1,1,nlf),nx_l,ny_l)
                  var_lvd(nlf) = 'SVN'       ! satellite, visible, normalized
                  c_lvd(nlf)=csatid//' (VISIBLE) SATELLITE - NORM'
                  units_lvd(nlf) = 'COUNTS'
               endif

               if(abs(istatus_vis(3)).lt.good_vis_data_thresh)then
                  nlf=nlf+1
                  call move(albedo,laps_data(1,1,nlf),nx_l,ny_l)
                  var_lvd(nlf) = 'ALB'       ! albedo
                  c_lvd(nlf)= csatid//' (VISIBLE) ALBEDO'
               endif

            else
               write(6,*)'vis not processed: too much missing data'
            endif

            endif

            write(6,*)'Number of lvd fields so far: ',nlf
            write(6,*)'     New Fields Written:'
            do k=nlf_prev,nlf
               write(6,132)var_lvd(k)
132            format(8x,a3)
            enddo
            nlf_prev = nlf+1
            write(6,*)

         enddo
         goto 311
c
c following routine handles the case for which the data have already
c been mapped to the laps domain. AFWA's GMS so far.

310      call loadlapsdata(nx_l,ny_l,maxchannels,n_lvd_fields_max,
     &                     ntm(i),c_type(1,i),r_image_status(1,i),
     &                     csatid,image_vis(1,1,i),
     &                     image_39(1,1,i),image_67(1,1,i),
     &                     image_ir(1,1,i),image_12(1,1,i),
     &                     var_lvd,c_lvd,units_lvd,
     &                     nlf,laps_data,istatus)

311      write(6,*)' Writing lvd. Total # of fields: ',nlf
         write(6,*)'    to ',dir_lvd(1:len_lvd)
         call write_laps_data(i4time_data(i),
     &                      dir_lvd,
     &                      ext_lvd,
     &                      nx_l,ny_l,
     &                      nlf,
     &                      nlf,
     &                      var_lvd,
     &                      lvl_lvd,
     &                      lvl_coord_lvd,
     &                      units_lvd,
     &                      c_lvd,
     &                      laps_data,
     &                      istatus)
         if(istatus.eq.1)then
            write(6,*)'*****************************'
            write(*,*)'lvd file successfully written'
            write(*,*)'for: ',c_fname
            write(*,*)'i4 time: ',i4time_data(i)
            write(6,*)'*****************************'
         else
            write(*,*)' Error writing lvd file for this time'
            write(*,*)' i4Time: ',i4time_data(i)
            write(*,*)' File Time: ',c_fname
         endif

997   enddo

      lvd_status = 1

      goto 16

 99   write(6,*)'Error opening count LUT: terminating. NO LVD'
      goto 16

909   write(6,*)'Error opening ll/ij look up table'
      goto 16

910   write(6,*)'Error getting mapping lut'
      goto 16

998   write(*,*)'No ',c_sat_id(isat),"/",c_sat_types(jtype,isat),
     &' satellite data.'

 16   write(*,*)' lvd driver sub completed'
      itstatus=ishow_timer()
      write(6,*)'Elapsed time (sec): ',itstatus
      return
      end
