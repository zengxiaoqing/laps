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
     &               i4time_cur,i_delta_t,
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
      parameter (max_files=500)
 
      integer   len_lvd

      character*3 csattype
      character*3 c_type(maxchannels,max_files)
      character*3 chtype(maxchannels)
      character*6 csatid
      character*9 c_fname_cur
      character*9 c_fname
      character*9 fname_ctp
      character*50 c_gridfname
      character*255 c_generic_dataroot

      real image_vis (n_vis_elem,n_vis_lines,nimages) 
      real image_ir  (n_ir_elem,n_ir_lines,nimages)
      real image_12  (n_ir_elem,n_ir_lines,nimages)
      real image_39  (n_ir_elem,n_ir_lines,nimages)
      real image_67  (n_wv_elem,n_wv_lines,nimages)
c
c this stuff for cloud top pressure
c
      real,     allocatable  :: rlctp(:,:)
      real,     allocatable  :: rlca(:,:)
      real,     allocatable  :: rlct(:,:)
      real,     allocatable  :: ri4time_ob(:,:)
      real,     allocatable  :: ctp_data(:,:,:)
c
c netCDF for cloud top pressure
c
      character*125 c_ctp(4)
      character*10 units_ctp(4)
      character*3 var_ctp(4)
      character*4 lvl_coord_ctp(4)
      integer lvl_ctp(4)
      integer ldctp
      character*150 dir_ctp
      character*31 ext_ctp

      character path_to_ctp*256
      character c8_project*8

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
      real      rcal

      logical   lvis_flag
      logical   lsatqc
      logical   l_lut_flag
      logical   objects_allocated

      integer   i,j,k,l,n
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
      integer i_delta_t
      integer i4time_data(max_files)
      integer i4time_ctp_data
      integer iwindow_ctp
      integer istat
      integer istatus
      integer istatus_vis(3)
      integer istatus_ctp
      integer itstatus
      integer lvd_status
      integer nft,ntm(max_files),nft_prior

      real*4 favgth39u
      real*4 favgth67u
      real*4 favgth11u
      real*4 favgth12u
c
c these should be seasonally dependent. 6-13-99. 
c used for "bad" meteosat (11u) data. 
      data favgth39u /210.0/
      data favgth67u /209.0/
      data favgth11u /239.0/
      data favgth12u /237.0/

      data rcal/0.106178/  !as specified by EUMETSAT User Services 5-May-99.

      include 'satellite_common_lvd.inc'
c =========================================================================
c ----------------------------- START -------------------------------------
c
      itstatus=init_timer()
      itstatus=ishow_timer()

      objects_allocated = .false.  !used for NIMBUS cloud_top_pressure processing

      do i=1,nimages
         call zero(image_ir(1,1,i),n_ir_elem,n_ir_lines)
         call zero(image_vis(1,1,i),n_vis_elem,n_vis_lines)
         call zero(image_67(1,1,i),n_wv_elem,n_wv_lines)
         call zero(image_12(1,1,i),n_ir_elem,n_ir_lines)
      enddo

      call find_domain_name(c_generic_dataroot,c_gridfname,istatus)

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
      call get_domain_laps(nx_l,ny_l,c_gridfname,lat,lon,topo,
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
      if(csatid.ne.'gmssat'.or.csattype.eq.'twn')then

         call readlut(csatid,csattype,maxchannels,nchannels,
     &chtype,nx_l,ny_l,r_llij_lut_ri,r_llij_lut_rj,istatus)

         if(istatus.eq.1)then
            write(6,*)'LUT not obtained: ',csatid,'/',csattype
            istatus = -1
            return
         else
            write(6,*)'Got the mapping look-up-tables '
         endif
         if(csattype.eq.'twn')then
            where (r_llij_lut_ri .lt. 0.5)r_llij_lut_ri=1.0
            where (r_llij_lut_rj .lt. 0.5)r_llij_lut_rj=1.0
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
         write(6,*)'lvis_flag set: do not wait for vis data'
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
     &                      i_delta_t,
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

         call s_len(path_to_raw_sat(1,jtype,isat),in)
         write(6,*)'datapath: ',path_to_raw_sat(1,jtype,isat)(1:in)
         write(6,*)'Using getascii_satdat routine'

c  only possible to have one time for ascii files (nft=1); however, the number of
c  matches for this time (ntm) >= 0 depending on the result in getascii_satdat.

         nft=1
         write(6,*)'ascii satellite data ingest currently disabled'

c        call getascii_satdat(i4time_cur,lvis_flag,i_delta_t,
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

         nft=1    !default to start. can change within this routine depending on AFWA file times.
         nft_prior=nft
         call getafgwc_satdat(isat,jtype,
     &                        maxchannels,nchannels,chtype,
     &                        i4time_cur,lvis_flag,
     &                        n_ir_lines,n_ir_elem,
     &                        n_vis_lines,n_vis_elem,
     &                        n_wv_lines,n_wv_elem,
     &                        nft,ntm,max_files,c_type,
     &                        image_ir,image_vis,
     &                        image_12,image_39,image_67,
     &                        i4time_data,
     &                        istatus)

         if(nft.gt.nft_prior.and.nft.lt.3)then
            print*,'nft incremented in getafgwc_satdat.'
     &,' Move image_ir(1,1,1) to image_ir(1,1,2)'
         call move(image_ir(1,1,1),image_ir(1,1,2),n_ir_elem,n_ir_lines)
         elseif(nft.ge.3)then
            print*,'nft > 2. lvd for AFWA cannot handle this. terminate'
            stop
         endif

         if(istatus.ne.0)then
            write(6,*)'Did not get data for ',c_fname_cur
            goto 998
         endif 

c June 2001 added Taiwan (gms) sat ingest

      elseif(csattype.eq.'twn')then

         call read_gms_taiwan(path_to_raw_sat(1,jtype,isat)
     &,n_lines_ir(jtype,isat),n_pixels_ir(jtype,isat)        !<-- full array size raw data
     &,maxchannels,max_files,nchannels,csatid,csattype
     &,chtype,i4time_cur,n_ir_elem,n_ir_lines,n_vis_elem
     &,n_vis_lines,n_wv_elem,n_wv_lines,image_ir,image_vis
     &,image_67,image_12,nimages,nft,ntm,c_type,i4time_data
     &,istatus)

         if(istatus.ne.1)then
            print*,'Error returned: read_gms_taiwan'
            return
         endif

      endif
c --------------------------------------------------------------------
c Get image resolution information
c --------------------------------
      do j = 1,nft
      do i = 1,ntm(j)

         call lvd_file_specifier(c_type(i,j),ispec,istatus)
         if(ispec.eq.1)then
            r_image_res_m(i,j)=r_resolution_x_vis(jtype,isat)
         elseif(ispec.eq.2.or.ispec.eq.4.or.ispec.eq.5)then
            r_image_res_m(i,j)=r_resolution_x_ir(jtype,isat)
         elseif(ispec.eq.3)then
            r_image_res_m(i,j)=r_resolution_x_wv(jtype,isat)
         endif

      enddo
      enddo
c
c -----------------------------------------------------------------------
c Compute or read ir/vis count to brightness temp (Tb)/vis count-to-count.
c -----------------------------------------------------------------------
c
      if(csattype.eq.'cdf'.or.csattype.eq.'wfo'.or.
     &   csatid.eq.'meteos')then

         write(6,*)'Compute ',csatid,' cnt-to-btemp lookup tables'

         do j = 1,nft
         do i = 1,ntm(j)

            call lvd_file_specifier(c_type(i,j),ispec,istat)

            if(ispec.eq.2)then
               call genbtemplut(csatid,csattype,ispec,rcal,
     &                          r39_cnt_to_btemp_lut,istatus)
               if(istatus.ne.1)then
                  write(6,*)'Error computing 39 lut'
               endif
            elseif(ispec.eq.3)then
               call genbtemplut(csatid,csattype,ispec,rcal,
     &                          r67_cnt_to_btemp_lut,istatus)
               if(istatus.ne.1)then
                  write(6,*)'Error computing 67 lut'
               endif
            elseif(ispec.eq.4)then
               call genbtemplut(csatid,csattype,ispec,rcal,
     &                          ir_cnt_to_btemp_lut,istatus)
               if(istatus.ne.1)then
                  write(6,*)'Error computing ir lut'
               endif
            elseif(ispec.eq.5)then
               call genbtemplut(csatid,csattype,ispec,rcal,
     &                          r12_cnt_to_btemp_lut,istatus)
               if(istatus.ne.1)then
                  write(6,*)'Error computing 12 lut'
               endif
            endif

         enddo
         enddo

      elseif(csattype.eq.'gvr'
     &.or.csattype.eq.'gwc'.or.csattype.eq.'twn')then

         write(6,*)'get gvarimage count to btemp lut'
         write(6,*)'and gvarimage vis count-to-count lut'

         do j = 1,nft
         do i = 1,ntm(j)

            call lvd_file_specifier(c_type(i,j),ispec,istat)

            if( ispec.eq.1.and.
     &        ( csatid.ne.'gmssat'.or.csatid.ne.'meteos'))then
               call read_gvarimg_cnt2btemp_lut(csatid,
     &c_type(i,j),vis_cnt_to_cnt_lut,istatus)
            elseif(ispec.eq.2)then
               call read_gvarimg_cnt2btemp_lut(csatid,
     &c_type(i,j),r39_cnt_to_btemp_lut,istatus)
            elseif(ispec.eq.3)then
               call read_gvarimg_cnt2btemp_lut(csatid,
     &c_type(i,j),r67_cnt_to_btemp_lut,istatus)
            elseif(ispec.eq.4)then
               call read_gvarimg_cnt2btemp_lut(csatid,
     &c_type(i,j),ir_cnt_to_btemp_lut,istatus)
            elseif(ispec.eq.5)then
               call read_gvarimg_cnt2btemp_lut(csatid,
     &c_type(i,j),r12_cnt_to_btemp_lut,istatus)

            endif

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
      if(csattype.eq.'asc'.or.csattype.eq.'gwc'.or.
     &csattype.eq.'twn')lsatqc=.false.
      
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
                if(ispec.eq.4)then

                   write(6,*)'For ',c_type(j,i)
                   call btemp_convert(n_ir_elem,n_ir_lines,
     &                        ir_cnt_to_btemp_lut,
     &                        r_missing_data,
     &                        image_ir(1,1,i))
                elseif(ispec.eq.2)then 
 
                   write(6,*)'For ',c_type(j,i)
                   call btemp_convert(n_ir_elem,n_ir_lines,
     &                        r39_cnt_to_btemp_lut,
     &                        r_missing_data,
     &                        image_39(1,1,i))
                elseif(ispec.eq.3)then

                   write(6,*)'For ',c_type(j,i)
                   call btemp_convert(n_wv_elem,n_wv_lines,
     &                        r67_cnt_to_btemp_lut,
     &                        r_missing_data,
     &                        image_67(1,1,i))

                elseif(ispec.eq.5)then

                   write(6,*)'For ',c_type(j,i)
                   call btemp_convert(n_ir_elem,n_ir_lines,
     &                        r12_cnt_to_btemp_lut,
     &                        r_missing_data,
     &                        image_12(1,1,i))

                elseif(ispec.eq.1)then

                   if( (csatid.ne.'gmssat'.and.csatid.ne.'meteos').and.
     &                 (csattype.eq.'gvr'.or.csattype.eq.'gwc') )then

                       call btemp_convert(n_vis_elem,n_vis_lines,
     &                          vis_cnt_to_cnt_lut,
     &                          r_missing_data,
     &                          image_vis(1,1,i))
                       write(*,*)'VIS data converted - cnt-to-cnt lut'
                   else
                       write(*,*)'Not converting ',csattype,' vis data'
                   endif
                endif
             endif
          enddo
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
             if(ispec.eq.2)then
                call btemp_convert_asc(n_ir_lines,n_ir_elem,
     &               r_missing_data,image_39(1,1,i),istatus)
             elseif(ispec.eq.3)then
                call btemp_convert_asc(n_wv_lines,n_wv_elem,
     &               r_missing_data,image_67(1,1,i),istatus)
             elseif(ispec.eq.4)then
                call btemp_convert_asc(n_ir_lines,n_ir_elem,
     &               r_missing_data,image_ir(1,1,i),istatus)
             elseif(ispec.eq.5)then
                call btemp_convert_asc(n_ir_lines,n_ir_elem,
     &               r_missing_data,image_12(1,1,i),istatus)
             endif

          enddo
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
         if(csatid.eq.'gmssat'.and.csattype.ne.'twn')goto 310


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
                  if(csatid.eq.'meteos')then
                     call check_field_ave(nx_l,ny_l,ta8,favgth11u
     &                                   ,istatus)
                  endif
                  if(istatus.eq.1)then

c                    if(csattype.eq.'twn')then
                        call filter_2dx(ta8,nx_l,ny_l,1, 0.5)
                        call filter_2dx(ta8,nx_l,ny_l,1,-0.5)
c                    endif

                     nlf=nlf+1
                     call move(ta8,laps_data(1,1,nlf),nx_l,ny_l)
                     var_lvd(nlf)  = 'S8A'    ! satellite, channel-4, averaged
                     c_lvd(nlf)=csatid//' (11.2u) IR B-TEMPS - AVERAGED'
                     nlf=nlf+1
                     call move(tb8,laps_data(1,1,nlf),nx_l,ny_l)
                     var_lvd(nlf)='S8W'       ! satellite, channel-4, warm pixel
                     c_lvd(nlf)=csatid//' (11.2u) IR B-TEMPS; WARM PIX'
                     nlf=nlf+1
                     call move(tc8,laps_data(1,1,nlf),nx_l,ny_l)
                     var_lvd(nlf)='S8C'       ! satellite, channel-4, warm pixel
                     c_lvd(nlf)=csatid//' (11.2u) IR B-TEMPS - FILTERED'
                   else
                     print*,'No output for this channel: ',c_type(j,k)
                   endif
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
                  if(csatid.eq.'meteos')then
                     call check_field_ave(nx_l,ny_l,ta4,favgth39u
     &                                   ,istatus)
                  endif
                  if(istatus.eq.1)then
                     nlf=nlf+1
                     call move(ta4,laps_data(1,1,nlf),nx_l,ny_l)
                     var_lvd(nlf) = 'S3A'       ! satellite, , averaged
                     c_lvd(nlf)=csatid//' (3.9u) IR B-TEMPS - AVERAGED'
                     nlf=nlf+1
                     call move(tb4,laps_data(1,1,nlf),nx_l,ny_l)
                     var_lvd(nlf)  = 'S3C'       ! satellite, , filtered
                     c_lvd(nlf)=csatid//' (3.9u) IR B-TEMPS - FILTERED'
                  else
                     print*,'No output for this channel: ',c_type(j,k)
                  end if
               endif
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
                  if(csatid.eq.'meteos')then
                     call check_field_ave(nx_l,ny_l,ta12,favgth12u
     &                                   ,istatus)
                  endif
                  if(istatus.eq.1)then
                     nlf=nlf+1
                     call move(ta12,laps_data(1,1,nlf),nx_l,ny_l)
                     var_lvd(nlf) = 'SCA'       ! satellite, averaged
                     c_lvd(nlf)=csatid//' (12.0u) IR B-TEMPS - AVERAGED'
                     nlf=nlf+1
                     call move(tb12,laps_data(1,1,nlf),nx_l,ny_l)
                     var_lvd(nlf) = 'SCC'       ! satellite, averaged
                     c_lvd(nlf)=csatid//' (12.0u) IR B-TEMPS - FILTERED'
                  else
                     print*,'No output for this channel: ',c_type(j,k)
                  end if
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
                  if(csatid.eq.'meteos')then
                     call check_field_ave(nx_l,ny_l,ta6,favgth67u
     &                                   ,istatus)
                  endif
                  if(istatus.eq.1)then
                     nlf=nlf+1
                     call move(ta6,laps_data(1,1,nlf),nx_l,ny_l)
                     var_lvd(nlf) = 'S4A'       ! satellite, averaged
                     c_lvd(nlf)=csatid//' (6.7u) IR B-TEMPS - AVERAGED'
                     nlf=nlf+1
                     call move(tb6,laps_data(1,1,nlf),nx_l,ny_l)
                     var_lvd(nlf) = 'S4C'       ! satellite, filtered
                     c_lvd(nlf)=csatid//' (6.7u) IR B-TEMPS - FILTERED'
                  else
                     print*,'No output for this channel: ',c_type(j,k)
                  end if
               endif
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
c
c *** istatus_v() is < 0. Determine if we have enough vis data.
c
               good_vis_data_thresh=(nx_l*ny_l)*0.1000
               if( (nx_l*ny_l+istatus_vis(1)).gt.
     &              good_vis_data_thresh)then

                  nlf=nlf+1
                  call move(visraw,laps_data(1,1,nlf),nx_l,ny_l)
                  var_lvd(nlf) = 'SVS'       ! satellite, visible
                  c_lvd(nlf)=csatid//' (VISIBLE) SATELLITE - RAW'
                  units_lvd(nlf) = 'COUNTS'
               endif

               if( (nx_l*ny_l+istatus_vis(2)).gt.
     &              good_vis_data_thresh)then

                  nlf=nlf+1
                  call move(visnorm,laps_data(1,1,nlf),nx_l,ny_l)
                  var_lvd(nlf) = 'SVN'       ! satellite, visible, normalized
                  c_lvd(nlf)=csatid//' (VISIBLE) SATELLITE - NORM'
                  units_lvd(nlf) = 'COUNTS'
               endif

               if( (nx_l*ny_l+istatus_vis(3)) .gt.
     &             good_vis_data_thresh)then

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

311      if(nlf .gt. 0)then
            write(6,*)' Writing lvd. Total # of fields: ',nlf
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

         else

            print*,'No fields processed. No lvd written ',c_fname

         endif

997   enddo

      goto 17

998   write(*,*)'No ',c_sat_id(isat),"/",c_sat_types(jtype,isat),
     &' satellite image data.'


17    call get_c8_project(c8_project,istatus)
      if(c8_project.eq.'NIMBUS')then

         print*,'check for new cloud top pressure (C02) files'
         iwindow_ctp=4000
         print*,'ctp time window (sec) = ',iwindow_ctp
         path_to_ctp='/public/data/sat/nesdis/goes8/cloudtop/ascii'
         call check_for_new_ctp(iwindow_ctp,istatus_ctp)

         if(istatus_ctp.eq.1)then
            ext_ctp='ctp'
            if(.not. objects_allocated)then
               allocate (rlctp(nx_l,ny_l),rlca(nx_l,ny_l)
     &                  ,rlct(nx_l,ny_l),ctp_data(nx_l,ny_l,4))
               allocate (ri4time_ob(nx_l,ny_l))
               objects_allocated = .true.
            endif

            print*,'calling read_cld_top_p'

            call read_cld_top_p(nx_l,ny_l,path_to_ctp
     &,rlctp,rlca,rlct,ri4time_ob,iwindow_ctp,i4time_ctp_data
     &,fname_ctp,istatus)
            if(istatus.ne. 1)then
               print*,'returned from read_cld_top_p'
               print*,'no data returned from read_cld_top_p'
               print*,'no ctp made for this time'
               print*
            else
c
c setup for writing ctp
c
             call get_directory('ctp',dir_ctp,ldctp)
c            call move(rlctp,ctp_data(1,1,1),nx_l,ny_l)
             ctp_data(:,:,1)=rlctp
             var_ctp(1) = 'PCT'
             c_ctp(1)=csatid//' NESDIS derived cloud top pressure'
             units_ctp(1)='PA'
             lvl_ctp(1) = 0
             lvl_coord_ctp(1)='AGL'

             ctp_data(:,:,2)=rlca
             var_ctp(2) = 'LCA'
             c_ctp(2)=csatid//' NESDIS derived cloud amount'
             units_ctp(2)='%'
             lvl_ctp(2) = 0
             lvl_coord_ctp(2)='AGL'

c            call move(rlct,ctp_data(1,1,3),nx_l,ny_l)
             ctp_data(:,:,3)=rlct
             var_ctp(3)='CTT'
             c_ctp(3)=csatid//' NESDIS derived cloud top temperature'
             units_ctp(3)='K'
             lvl_ctp(3) = 0
             lvl_coord_ctp(3)='AGL'

c            call move(ri4time_ob,ctp_data(1,1,4),nx_l,ny_l)
             ctp_data(:,:,4)=ri4time_ob
             var_ctp(4)='I4T'
             c_ctp(4)=csatid//' i4time of obs'
             units_ctp(4)='sec'
             lvl_ctp(4) = 0
             lvl_coord_ctp(4)='AGL'

             print*,'writing cld top pressure file ',ext_ctp
             print*,'dir = ',dir_ctp(1:ldctp)
c            dir_ctp=dir_ctp(1:ldctp)//'/'//fname_ctp
c            print*,'path/filename out: ',dir_ctp

             call write_laps_data(i4time_ctp_data
     &,dir_ctp,ext_ctp,nx_l,ny_l,4,4,var_ctp,lvl_ctp
     &,lvl_coord_ctp,units_ctp,c_ctp,ctp_data,istatus)
             if(istatus.eq.1)then
                write(6,*)'*****************************'
                write(*,*)'ctp file successfully written'
                write(*,*)'for: ',fname_ctp
                write(*,*)'i4 time: ',i4time_ctp_data
                write(6,*)'*****************************'
             else
                write(*,*)' Error writing ctp file for this time'
                write(*,*)' i4Time: ',i4time_ctp_data
                write(*,*)' File Time: ',fname_ctp
             endif

            endif

         else
            print*,'No new cld top pressure files to process'
            print*
         endif

      endif

      lvd_status = 1

      goto 16

 99   write(6,*)'Error opening count LUT: terminating. NO LVD'
      goto 16

909   write(6,*)'Error opening ll/ij look up table'
      goto 16

910   write(6,*)'Error getting mapping lut'

 16   write(*,*)' lvd driver sub completed'
      itstatus=ishow_timer()
      write(6,*)'Elapsed time (sec): ',itstatus
      return
      end
