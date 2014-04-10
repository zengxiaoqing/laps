cdis   
cdis    Open Source License/Disclaimer, Forecast Systems Laboratory
cdis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
cdis    
cdis    This software is distributed under the Open Source Definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    In particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - Redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - Redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - All modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - If significant modifications or enhancements are made to this
cdis    software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
cdis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
cdis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
cdis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
cdis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
cdis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
cdis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
cdis   
cdis
cdis
cdis   
cdis

        subroutine laps_deriv_sub(i4time,
     1                        NX_L,NY_L,
     1                        NZ_L,
     1                        N_PIREP,
     1                        maxstns,
     1                        max_snd_grid,max_snd_levels,
     1                        n_prods,
     1                        iprod_number,
     1                        temp_3d,
     1                        heights_3d,
     1                        rh_3d_pct,
     1                        pres_sfc_pa,
     1                        t_sfc_k,
     1                        dbz_max_2d,istatus_lps,    ! O
     1                        twet_snow,                 ! O
     1                        j_status,istatus)

        use cloud_rad ! Cloud Radiation and Microphysics Parameters

        integer       ss_normal,sys_bad_prod,sys_no_data,
     1                  sys_abort_prod

        parameter (ss_normal      =1, ! success
     1             sys_bad_prod   =2, ! inappropriate data, insufficient data
     1             sys_no_data    =3, ! no data
     1             sys_abort_prod =4) ! failed to make a prod

!       1991     Steve Albers - Original Version
!       1993 Mar Steve Albers - Add LMT product
!       1995 Jul Steve Albers - Fix coverage threshold on cloud base
!                               Added cloud ceiling
!                               Added sfc or 2-D cloud type.
!       1995 Nov 1  S. Albers - Add diagnostic output of cloud ceiling
!                               in precip type comparisons
!       1995 Nov 2  S. Albers - Use SAO's to add drizzle in sfc "thresholded"
!                               precip type calculation (lct-PTT)
!       1995 Nov 10 S. Albers - Better handling of cloud tops at top of domain
!                               when "three-dimensionalizing" radar data
!       1995 Nov 29 S. Albers - Improve use of SAO's to add snow in PTT field.
!                               Cloud ceiling threshold replaced with
!                               thresholds on cloud cover and sfc dewpoint
!                               depression.
!       1995 Dec 4  S. Albers - Use SAO's to add rain in sfc "thresholded"
!                               precip type calculation (lct-PTT)
!       1995 Dec 13 S. Albers - Now calls get_radar_ref
!       1996 Aug 22 S. Albers - Now calls read_radar_3dref
!       1996 Oct 10 S. Albers - Max SAO cloud cover is now 1.00 + some other
!                               misc cleanup.
!       1997 Jul 31 K. Dritz  - Removed include of lapsparms.for.
!       1997 Jul 31 K. Dritz  - Added call to get_i_perimeter.
!       1997 Jul 31 K. Dritz  - Removed PARAMETER statements for IX_LOW,
!                               IX_HIGH, IY_LOW, and IY_HIGH, and instead
!                               compute them dynamically (they are not used
!                               as array bounds, only passed in a call).
!       1997 Jul 31 K. Dritz  - Added NX_L, NY_L as dummy arguments.
!       1997 Jul 31 K. Dritz  - Added call to get_r_missing_data.
!       1997 Jul 31 K. Dritz  - Removed PARAMETER statements for default_base,
!                               default_top, and default_ceiling, and instead
!                               compute them dynamically.
!       1997 Jul 31 K. Dritz  - Removed PARAMETER statement for Nhor.  Now
!                               initialize c1_name_array dynamically instead
!                               of with a DATA statement.
!       1997 Jul 31 K. Dritz  - Added NZ_L as dummy argument.
!       1997 Jul 31 K. Dritz  - Added call to get_ref_base.
!       1997 Aug 01 K. Dritz  - Added maxstns, IX_LOW, IX_HIGH, IY_LOW, and
!                               IY_HIGH as arguments in call to insert_sao.
!       1997 Aug 01 K. Dritz  - Also now pass r_missing_data to barnes_r5.
!       1997 Aug 01 K. Dritz  - Pass r_missing_data to insert_sat.
!       1997 Aug 01 K. Dritz  - Pass ref_base to rfill_evap.

!       Prevents clearing out using satellite (hence letting SAOs dominate)
!       below this altitude (M AGL)
        real surface_sao_buffer
        parameter (surface_sao_buffer = 800.)

        real thresh_cvr,default_top,default_base,default_clear_cover
     1                   ,default_ceiling

        parameter       (thresh_cvr = 0.65) ! Used to "binaryize" cloud cover

        parameter       (default_clear_cover = .01)

        real thresh_cvr_base,thresh_cvr_top,thresh_cvr_ceiling
        parameter (thresh_cvr_base = 0.1)
        parameter (thresh_cvr_top  = 0.1)
        parameter (thresh_cvr_ceiling = thresh_cvr)

        real thresh_thin_lwc_ice     ! Threshold cover for thin cloud LWC/ICE
        parameter (thresh_thin_lwc_ice = 0.1)

        real vis_radar_thresh_cvr,vis_radar_thresh_dbz
        parameter (vis_radar_thresh_cvr = 0.2)  ! 0.2, 0.0
        parameter (vis_radar_thresh_dbz = 10.)  ! 5. , -99.

        real lat(NX_L,NY_L),lon(NX_L,NY_L)
        real topo(NX_L,NY_L)
        real rlaps_land_frac(NX_L,NY_L)
        real solar_alt(NX_L,NY_L)
        real solar_ha(NX_L,NY_L)

        real k_to_c

        logical l_packed_output
        logical l_evap_radar

        data l_packed_output /.false./
        data l_evap_radar /.false./

        logical l_fill
        logical l_flag_mvd
        logical l_flag_cloud_type
        logical l_flag_icing_index
        logical l_flag_bogus_w, l_bogus_radar_w, l_deep_vv
        logical l_flag_pcp_type        
        logical l_parse

        logical l_sao_lso
        data l_sao_lso /.true./ ! Do things the new way?

        logical l_perimeter
        data l_perimeter /.true./ ! Use SAOs just outside domain?

        include 'laps_cloud.inc'

!       Nominal cloud heights. Actual ones used are fitted to the terrain.
        real cld_hts_new(KCLOUD)

        data cld_hts_new/1200.,1300.,1400.,1500.,1600.,1700.,1800.,
     11900.,2000.,2100.,2200.,2400.,2600.,2800.,3000.,3200.,
     23400.,3600.,3800.,4000.,4500.,5000.,5500.,6000.,6500.,
     37000.,7500.,8000.,8500.,9000.,9500.,10000.,11000.,12000.,
     413000.,14000.,15000.,16000.,17000.,18000.,19000.,20000./

        integer iarg

        equivalence (cld_hts,cld_hts_new)

        real clouds_3d(NX_L,NY_L,KCLOUD)

        integer ihist_alb(-10:20)

        real cloud_top(NX_L,NY_L)
        real cloud_base(NX_L,NY_L)
        real cloud_ceiling(NX_L,NY_L)

        real cldtop_m(NX_L,NY_L)
        real cldtop_m_co2(NX_L,NY_L)
        real cldtop_m_tb8(NX_L,NY_L)

        real cld_pres_1d(KCLOUD)
        real pres_3d(NX_L,NY_L,NZ_L)
        real clouds_3d_pres(NX_L,NY_L,NZ_L)

        real CVHZ(NX_L,NY_L)
        real CVHZ1(NX_L,NY_L),CVEW1(NX_L,KCLOUD)
        real cvr_max(NX_L,NY_L),CVEW2(NX_L,KCLOUD)
        real cvr_sao_max(NX_L,NY_L)
        real cvr_snow_cycle(NX_L,NY_L)
        real cvr_water_temp(NX_L,NY_L)
        real cvr_snow(NX_L,NY_L)
        real band8_mask(NX_L,NY_L)

        character*4 radar_name
        character*31 radarext_3d_cloud
        real radar_ref_3d(NX_L,NY_L,NZ_L)
        real heights_3d(NX_L,NY_L,NZ_L)

        Real vv_to_height_ratio_Cu
        Real vv_to_height_ratio_Sc
        Real vv_for_St

        real mvd_3d(NX_L,NY_L,NZ_L)
!       real lwc_res_3d(NX_L,NY_L,NZ_L)
        real w_3d(NX_L,NY_L,NZ_L)

        integer icing_index_3d(NX_L,NY_L,NZ_L)

        integer cldpcp_type_3d(NX_L,NY_L,NZ_L) ! Also contains 3D precip type

!       Output array declarations
        real out_array_3d(NX_L,NY_L,NZ_L)

        real, allocatable, dimension(:,:,:) :: slwc
        real, allocatable, dimension(:,:,:) :: cice
        real slwc_int(NX_L,NY_L)
        real cice_int(NX_L,NY_L)
        real rain_int(NX_L,NY_L)
        real snow_int(NX_L,NY_L)
        real pice_int(NX_L,NY_L)

        real pcpcnc(NX_L,NY_L,NZ_L)
        real raicnc(NX_L,NY_L,NZ_L)
        real snocnc(NX_L,NY_L,NZ_L)
        real piccnc(NX_L,NY_L,NZ_L)

        real cldamt(NX_L,NY_L)
        real cldalb_in(NX_L,NY_L)
        real cldalb_out(NX_L,NY_L)
        real cldod_out(NX_L,NY_L)

!       real snow_2d(NX_L,NY_L)

        character*2 c2_precip_types(0:10)

        character*20 c_z2m

        data c2_precip_types
     1  /'  ','Rn','Sn','Zr','Sl','Ha','L ','ZL','  ','  ','  '/

        character*3 c3_pt_flag
        character*1 c1_r,c1_s
        character*8 c8_project

        integer i2_pcp_type_2d(NX_L,NY_L)
        real r_pcp_type_2d(NX_L,NY_L)

        real dum1_array(NX_L,NY_L)
        real dum2_array(NX_L,NY_L)
        real dum3_array(NX_L,NY_L)
        real dum4_array(NX_L,NY_L)

      ! Used for "Potential" Precip Type
        logical l_mask_pcptype(NX_L,NY_L)
        integer ibase_array(NX_L,NY_L)
        integer itop_array(NX_L,NY_L)

        logical l_unresolved(NX_L,NY_L)

        character*1 c1_name_array(NX_L,NY_L)
        character*9 filename

        character*35 TIME
        character*13 filename13

        integer MAX_FIELDS
        parameter (MAX_FIELDS = 10)

        character*255 c_filespec
        character var*3,comment*125,directory*150,ext*31,units*10
        character*3 exts(20)
        character*3 var_a(MAX_FIELDS)
        character*125 comment_a(MAX_FIELDS)
        character*10  units_a(MAX_FIELDS)

!       Arrays used to read in satellite data
        real tb8_k(NX_L,NY_L)
        real tb8_cold_k(NX_L,NY_L)
        real albedo(NX_L,NY_L)
        real cloud_frac_vis_a(NX_L,NY_L)
        real cloud_frac_co2_a(NX_L,NY_L)

        real temp_3d(NX_L,NY_L,NZ_L)
        real rh_3d_pct(NX_L,NY_L,NZ_L)
        real model_3d(NX_L,NY_L,NZ_L)

        real t_sfc_k(NX_L,NY_L)
        real t_gnd_k(NX_L,NY_L)
        real sst_k(NX_L,NY_L)
        real td_sfc_k(NX_L,NY_L)
        real pres_sfc_pa(NX_L,NY_L)

!       Declarations for LSO file stuff
        real lat_s(maxstns), lon_s(maxstns), elev_s(maxstns)
        real cover_s(maxstns)
        real t_s(maxstns), td_s(maxstns), pr_s(maxstns), sr_s(maxstns)
        real dd_s(maxstns), ff_s(maxstns)
        real ddg_s(maxstns), ffg_s(maxstns)
        real vis_s(maxstns)
        real pstn_s(maxstns),pmsl_s(maxstns),alt_s(maxstns)
        real store_hgt(maxstns,5),ceil(maxstns),lowcld(maxstns)
        real cover_a(maxstns),rad_s(maxstns)
        integer obstime(maxstns),kloud(maxstns),idp3(maxstns)
        character store_emv(maxstns,5)*1,store_amt(maxstns,5)*4
        character wx_s(maxstns)*8, obstype(maxstns)*8
        character atime*24, infile*256

        integer STATION_NAME_LEN
        parameter (STATION_NAME_LEN = 3)                   
        character c_stations(maxstns)*(STATION_NAME_LEN)    

        character asc_tim_9*9

        real ri_s(maxstns), rj_s(maxstns)

!       Product # notification declarations
        integer j_status(20),iprod_number(20)

!       Stuff for 2d fields
        real ref_mt_2d(NX_L,NY_L)
        real dbz_low_2d(NX_L,NY_L)
        real dbz_max_2d(NX_L,NY_L)

!       SFC precip and cloud type (LCT file)
        real r_pcp_type_thresh_2d(NX_L,NY_L)
        real r_cld_type_2d(NX_L,NY_L)

        character*40 c_vars_req
        character*180 c_values_req

        character*3 lso_ext
        data lso_ext /'lso'/

        ISTAT = INIT_TIMER()

        write(6,*)' Welcome to the laps_deriv_sub (derived cloud prods)'       

        allocate( slwc(NX_L,NY_L,NZ_L), STAT=istat_alloc )
        if(istat_alloc .ne. 0)then
            write(6,*)' ERROR: Could not allocate slwc'
            stop
        endif

        allocate( cice(NX_L,NY_L,NZ_L), STAT=istat_alloc )
        if(istat_alloc .ne. 0)then
            write(6,*)' ERROR: Could not allocate cice'
            stop
        endif

        call get_r_missing_data(r_missing_data,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'Error calling get_r_missing_data'
           stop
        endif

        call get_deriv_parms(mode_evap,l_bogus_radar_w,              ! O
     1                       l_deep_vv,                              ! O
     1                       vv_to_height_ratio_Cu,                  ! O
     1                       vv_to_height_ratio_Sc,                  ! O
     1                       vv_for_St,                              ! O
     1                       c_z2m,                                  ! O
     1                       thresh_cvr_cty_vv,thresh_cvr_lwc,       ! O
     1                       twet_snow,                              ! O
     1                       hydrometeor_scale_cldliq,               ! O
     1                       hydrometeor_scale_cldice,               ! O
     1                       hydrometeor_scale_pcp,                  ! O
     1                       istatus)                                ! O
        if (istatus .ne. 1) then
           write (6,*) 'Error calling get_deriv_parms'
           stop
        endif

        if(mode_evap .gt. 0)l_evap_radar = .true.

        default_base     = r_missing_data
        default_top      = r_missing_data
        default_ceiling  = r_missing_data

        do j = 1,NY_L
           do i = 1,NX_L
              c1_name_array(i,j) = ' '
           enddo
        enddo

        call get_ref_base(ref_base,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'Error getting ref_base'
           stop
        endif

c Determine the source of the radar data
        c_vars_req = 'radarext_3d'

        call get_static_info(c_vars_req,c_values_req,1,istatus)

        if(istatus .eq. 1)then
            radarext_3d_cloud = c_values_req(1:3)
        else
            write(6,*)' Error getting radarext_3d'
            goto 9999
        endif

!       radarext_3d_cloud = radarext_3d
!       radarext_3d_cloud = 'v02'

        write(6,*)' radarext_3d_cloud = ',radarext_3d_cloud

c read in laps lat/lon and topo
        call get_laps_domain_95(NX_L,NY_L,lat,lon,topo
     1           ,rlaps_land_frac,grid_spacing_cen_m,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error getting LAPS domain'
            goto 9999
        endif

        call get_laps_cycle_time(ilaps_cycle_time,istatus)
        if(istatus .eq. 1)then
            write(6,*)' ilaps_cycle_time = ',ilaps_cycle_time
        else
            write(6,*)' Error getting laps_cycle_time'
            goto 9999
        endif

        call get_pres_3d(i4time,NX_L,NY_L,NZ_L,pres_3d,istatus)
        if(istatus .ne. 1)goto 9999

        n_lc3 = 1
        n_lps = 2
        n_lcb = 3
        n_lcv = 4
        n_lcp = 5
        n_lwc = 6
        n_lil = 7
        n_lct = 8
        n_lmd = 9
        n_lco = 10
        n_lrp = 11
        n_lty = 12
        n_lmt = 13

        exts(n_lc3) = 'lc3'
        exts(n_lcp) = 'lcp'
        exts(n_lwc) = 'lwc'
        exts(n_lil) = 'lil'
        exts(n_lcb) = 'lcb'
        exts(n_lct) = 'lct'
        exts(n_lcv) = 'lcv'
        exts(n_lmd) = 'lmd'
        exts(n_lco) = 'lco'
        exts(n_lps) = 'lps'
        exts(n_lrp) = 'lrp'
        exts(n_lty) = 'lty'
        exts(n_lmt) = 'lmt'

        do i = 1,13
            j_status(i) = sys_no_data
        enddo

        n_prods = 9
        iprod_start = 5
        iprod_end = 13

        if(.true.)then                    ! Read data, then calc derived fields
            I4_elapsed = ishow_timer()

            write(6,*)
            write(6,*)'Reading lc3,lt1,lps,lsx,lcv,lcb '
     1               ,'to calculate derived fields'

!           Read in data (lc3 - clouds_3d)
            ext = 'lc3'
            call get_clouds_3dgrid(i4time,i4time_lc3
     1                         ,NX_L,NY_L,KCLOUD,ext
     1                         ,clouds_3d,cld_hts,cld_pres_1d,istatus)
            if(istatus .ne. 1 .or. i4time .ne. i4time_lc3)then
                write(6,*)' Error reading 3D Clouds'
                goto 999
            endif

!           Read in data (lps - radar_ref_3d)
            var = 'REF'
            ext = 'lps'
            call get_laps_3d(i4time,NX_L,NY_L,NZ_L
     1       ,ext,var,units,comment,radar_ref_3d,istatus_lps)

            if(istatus_lps .ne. 1)then
                write(6,*)' Warning: could not read lps 3d ref, filling'
     1                   ,' array with r_missing_data'
                call constant_3d(radar_ref_3d,r_missing_data
     1                          ,NX_L,NY_L,NZ_L)           
                istat_radar_2dref = 0
                istat_radar_3dref = 0
                istat_radar_3dref_orig = 0

            else  ! istatus_lps = 1
                read(comment,510)istat_radar_2dref,istat_radar_3dref
     1                          ,istat_radar_3dref_orig
 510            format(23x,3i3)

!               Obtain column max ref
                call get_max_reflect(radar_ref_3d,NX_L,NY_L,NZ_L
     1                              ,ref_base,dbz_max_2d)

            endif ! istatus_lps

            write(6,*)' istatus_lps = ',istatus_lps

            var = 'LCV'
            ext = 'lcv'
            call get_laps_2d(i4time,ext,var,units,comment
     1                  ,NX_L,NY_L,cvr_max,istatus)

            if(istatus .ne. 1)THEN
                write(6,*)' Error Reading Cvr_max Analysis - abort'
                goto999
            endif

            var = 'CLA'
            ext = 'lcv'
            call get_laps_2d(i4time,ext,var,units,comment
     1                  ,NX_L,NY_L,cldalb_in,istatus)

            if(istatus .ne. 1 .and. istatus .ne. -1)THEN
                write(6,*)' Error Reading Cloud Albedo Analysis - abort'
                goto999
            endif

            var = 'CCE'
            ext = 'lcb'
            call get_laps_2d(i4time,ext,var,units,comment
     1                   ,NX_L,NY_L,cloud_ceiling,istatus)
            if(istatus .ne. 1 .and. istatus .ne. -1)THEN
                write(6,*)' Error Reading cld_ceiling Analysis - abort'       
                goto999
            endif

            var = 'LCB'
            ext = 'lcb'
            call get_laps_2d(i4time,ext,var,units,comment
     1                   ,NX_L,NY_L,cloud_base,istatus)
            if(istatus .ne. 1 .and. istatus .ne. -1)THEN
                write(6,*)' Error Reading cld_base Analysis - abort'       
                goto999
            endif

!           Access SAO data from LSO files
            ext = 'lso'
            call get_directory(ext,directory,len_dir) ! Returns directory
            infile = directory(1:len_dir)//filename13(i4time,ext(1:3))
            call read_surface_old(infile,maxstns,atime,n_meso_g,
     1           n_meso_pos,
     1           n_sao_g,n_sao_pos_g,n_sao_b,n_sao_pos_b,n_obs_g,
     1           n_obs_pos_g,n_obs_b,
     1           n_obs_pos_b,                   ! We use this as an obs count
     1           c_stations,obstype,        
     1           lat_s,lon_s,elev_s,wx_s,t_s,td_s,dd_s,ff_s,ddg_s,
     1           ffg_s,pstn_s,pmsl_s,alt_s,kloud,ceil,lowcld,cover_a,
     1           rad_s,idp3,store_emv,
     1           store_amt,store_hgt,vis_s,obstime,istat_sfc)

            I4_elapsed = ishow_timer()

            call get_c8_project(c8_project,istatus)
            if(istatus .ne. 1)goto 999

            do i = 1,n_obs_pos_b 

!               Does this station not report precip?
                if(obstype(i)(1:4) .eq. 'MESO'     .OR.
     1             obstype(i)(1:4) .eq. 'CDOT'     .OR.
     1             obstype(i)(7:8) .eq. '1A'       .OR.
     1             wx_s(i)(1:7)    .eq. 'UNKNOWN'  .OR.       ! Already in LSO
     1             l_parse(c8_project,'AFGWC')     .OR.
     1             l_parse(c8_project,'AFWA')             )then
                    wx_s(i) = 'UNKNOWN'

                endif

            enddo ! i

        endif

!       Calculate derived fields
        write(6,*)
        write(6,*)' Calculating Derived Fields'

!       Write out cloud grid in pressure coordinates
        write(6,*)' Writing out grid in pressure coordinates'

        call interp_height_pres_fast(NX_L,NY_L,NZ_L,kcloud
     1  ,clouds_3d_pres,clouds_3d,heights_3d,cld_hts,istatus)

        var = 'LCP'
        ext = 'lcp'
        units = 'Fractional'
        comment = 'LAPS Cloud Cover'
        call put_laps_3d(i4time,ext,var,units,comment,clouds_3d_pres
     1                                          ,NX_L,NY_L,NZ_L)
        j_status(n_lcp) = ss_normal
        I4_elapsed = ishow_timer()


!       Calculate and write out LWC, MVD, and Icing Index
        write(6,*)
        write(6,*)' Calling LWC etc. routine (get_cloud_deriv)'
        iflag_slwc = 13 ! New SMF LWC
        l_flag_mvd = .true.
        l_flag_cloud_type = .true.
        l_flag_icing_index = .true.
        l_flag_bogus_w = .true.
        l_flag_pcp_type = .true.

!       if(.not. l_flag_cloud_type)then ! Read in instead
!       endif

        call get_cloud_deriv(
     1                NX_L,NY_L,NZ_L,clouds_3d,cld_hts,
     1                temp_3d,rh_3d_pct,heights_3d,pres_3d,
     1                istat_radar_3dref,radar_ref_3d,grid_spacing_cen_m,       
     1                l_mask_pcptype,ibase_array,itop_array,
     1                iflag_slwc,slwc,cice,
     1                thresh_cvr_cty_vv,thresh_cvr_lwc,
     1                l_flag_cloud_type,cldpcp_type_3d,
     1                l_flag_mvd,mvd_3d,
     1                l_flag_icing_index,icing_index_3d,
     1                vv_to_height_ratio_Cu,                               ! I
     1                vv_to_height_ratio_Sc,                               ! I
     1                vv_for_St,                                           ! I
     1                l_flag_bogus_w,w_3d,l_bogus_radar_w,                 ! I
     1                l_deep_vv,                                           ! I
     1                twet_snow,                                           ! I
     1                l_flag_pcp_type,                                     ! I
     1                istatus)                                             ! O
        if(istatus .ne. 1)then
            write(6,*)' Bad status return from get_cloud_deriv'
            goto 999
        endif

        I4_elapsed = ishow_timer()

        if(l_flag_cloud_type)then

!           Write 3D Cloud Type
!           4 most significant bits are precip type, other 4 are cloud type
            do k = 1,NZ_L
            do j = 1,NY_L
            do i = 1,NX_L
                iarg = cldpcp_type_3d(i,j,k)
                out_array_3d(i,j,k) = iarg - iarg/16*16         ! 'CTY'
            enddo
            enddo
            enddo

            ext = 'cty'
            var = 'CTY'
            units = 'NONE'
            comment = 
     1         'Cloud Type: (1-10) - St,Sc,Cu,Ns,Ac,As,Cs,Ci,Cc,Cb'
            call put_laps_3d(i4time,ext,var,units
     1                      ,comment,out_array_3d             
     1                      ,NX_L,NY_L,NZ_L)

            I4_elapsed = ishow_timer()

        endif ! l_flag_cloud_type

!       Calculate "SFC" cloud type
!       Now this is simply the cloud type of the lowest significant
!       (> 0.65 cvr) layer present. If a CB is present, it is used instead.
        do i = 1,NX_L
        do j = 1,NY_L
            r_cld_type_2d(i,j) = 0

!           Pick lowest "significant" layer
            do k = NZ_L,1,-1
                cld_type_3d = out_array_3d(i,j,k)
                if(cld_type_3d .gt. 0)then
                    r_cld_type_2d(i,j) = cld_type_3d
                endif
            enddo ! k

!           Pick a CB if present
            do k = NZ_L,1,-1
                cld_type_3d = out_array_3d(i,j,k)
                if(cld_type_3d .eq. 10)then
                    r_cld_type_2d(i,j) = cld_type_3d
                endif
            enddo ! k

        enddo ! j
        enddo ! i


!       Convert SLWC and CICE by applying scale factor to parcel method values
        if(hydrometeor_scale_cldliq .ge. 0.)then
            ratio_cldliq =  hydrometeor_scale_cldliq
        else
            ratio_cldliq = -hydrometeor_scale_cldliq / 
     1                  (grid_spacing_cen_m/1000.)
        endif

        if(hydrometeor_scale_cldice .ge. 0.)then
            ratio_cldice =  hydrometeor_scale_cldice
        else
            ratio_cldice = -hydrometeor_scale_cldice / 
     1                  (grid_spacing_cen_m/1000.)
        endif

        cld_ice_ub_gpm3 = .03

        do k = 1,NZ_L
        do j = 1,NY_L
        do i = 1,NX_L
            if(slwc(i,j,k) .ne. r_missing_data)then
                slwc(i,j,k) = (slwc(i,j,k) * ratio_cldliq)
            endif
            if(cice(i,j,k) .ne. r_missing_data)then
!               cice(i,j,k) = (cice(i,j,k) * ratio_cldice)
                cice(i,j,k) = min(cice(i,j,k),cld_ice_ub_gpm3)                
        
            endif
        enddo 
        enddo
        enddo

        write(6,*)
        write(6,*)' Inserting thin clouds into LWC/ICE fields'
        call insert_thin_lwc_ice(clouds_3d,clouds_3d_pres,heights_3d
     1       ,temp_3d,cldalb_in,cld_hts,NX_L,NY_L,NZ_L,KCLOUD
     1       ,thresh_thin_lwc_ice       
     1       ,pres_3d,slwc,cice,istatus)

        if(istatus .ne. 1)then
            write(6,*)' Bad status return from insert_thin_lwc_ice'
            goto 999
        endif

        where(cice(:,:,:) .ne. r_missing_data .AND.
     1        cice(:,:,:) .gt. cld_ice_ub_gpm3      )    
              cice(:,:,:) = cld_ice_ub_gpm3
        end where

        I4_elapsed = ishow_timer()

!       Calculate and Write Integrated LWC
        write(6,*)
        write(6,*)' Calculating Integrated LWC and CICE'
        call integrate_slwc(slwc,heights_3d,NX_L,NY_L,NZ_L,slwc_int)
        call integrate_slwc(cice,heights_3d,NX_L,NY_L,NZ_L,cice_int)

!       Calculate cloud optical depth and cloud albedo
        const_lwp = (1.5 * rholiq) / (rholiq     * reff_clwc)
        const_lwp_bks = const_lwp * bksct_eff_clwc

        const_iwp = (1.5 * rholiq) / (rholiq     * reff_cice)
        const_iwp_bks = const_iwp * bksct_eff_cice

        const_rwp = (1.5 * rholiq) / (rholiq     * reff_rain)
        const_rwp_bks = const_rwp * bksct_eff_rain
  
        const_swp = (1.5 * rholiq) / (rhosnow    * reff_snow)
        const_swp_bks = const_swp * bksct_eff_snow

        const_gwp = (1.5 * rholiq) / (rhograupel * reff_graupel)
        const_gwp_bks = const_gwp * bksct_eff_graupel

!       Cloud amount is opacity of cloud liquid and cloud ice hydrometeors
        do j = 1,NY_L
        do i = 1,NX_L
            cldamt(i,j) = 1. - (exp( -(const_lwp * slwc_int(i,j) 
     1                               + const_iwp * cice_int(i,j))) ) 
        enddo ! i
        enddo ! j

        I4_elapsed = ishow_timer()

!       DERIVED RADAR/PRECIP STUFF
        if(istat_radar_3dref .eq. 1)then ! LMT

            if(l_evap_radar)then 

                write(6,*)' Calling rfill_evap: mode_evap = ',mode_evap       

!               Use LPS reflectivity field

!               Use Cloud Base field

                I4_elapsed = ishow_timer()

!               Apply evaporation subroutine
                call rfill_evap(radar_ref_3d,NX_L,NY_L,NZ_L
     1          ,cloud_base,lat,lon,topo,mode_evap
     1          ,temp_3d,rh_3d_pct,cldpcp_type_3d,heights_3d,istatus
     1          ,ref_base)

                I4_elapsed = ishow_timer()

            endif ! l_evap_radar, etc.

          ! Do Max Tops

            I4_elapsed = ishow_timer()

            write(6,*)' Getting Max Tops'

            call get_maxtops(radar_ref_3d,heights_3d,NX_L,NY_L,NZ_L
     1                      ,ref_mt_2d)

!           Get LAPS reflectivities at the surface (or immediately above it)
            write(6,*)' Getting Low Level Reflectivity'
            call get_low_ref(radar_ref_3d,pres_sfc_pa,NX_L,NY_L,NZ_L
     1                      ,dbz_low_2d)

            istat_pty = 0

!           Do SFC precip type
            var = 'TD'
            ext = 'lsx'
            call get_laps_2d(i4time,ext,var,units,comment
     1                  ,NX_L,NY_L,td_sfc_k,istatus)
            if(istatus .ne. 1)then
                write(6,*)
     1          ' Error reading SFC TD - SFC Precip Type not computed'       
                goto700
            endif

            I4_elapsed = ishow_timer()

!           Note that pres_sfc_pa was already read in above
            call get_sfc_preciptype(pres_sfc_pa,t_sfc_k,td_sfc_k
     1                             ,cldpcp_type_3d,twet_snow
     1                             ,dbz_low_2d,i2_pcp_type_2d
     1                             ,NX_L,NY_L,NZ_L)

!           Compute thresholded precip type
            do i = 1,NX_L
            do j = 1,NY_L
                iarg = i2_pcp_type_2d(i,j)
                r_pcp_type_2d(i,j) = iarg/16

!               Apply a threshold to the sfc precip type that depends
!               on both the low level reflectivity and the "potential"
!               sfc precip type.

                if(r_pcp_type_2d(i,j) .eq. 1.0        ! Rain
     1        .or. r_pcp_type_2d(i,j) .eq. 3.0        ! ZR
     1        .or. r_pcp_type_2d(i,j) .eq. 4.0        ! IP
     1        .or. r_pcp_type_2d(i,j) .eq. 5.0        ! Hail
     1                                                      )then

                    r_pcp_type_thresh_2d(i,j) = r_pcp_type_2d(i,j)

                    if(dbz_low_2d(i,j) .lt. 13.0)then
                        r_pcp_type_thresh_2d(i,j) = 0. ! No Precip
                    endif
                else ! Apply dry threshold to snow

                    call nowrad_virga_correction(r_pcp_type_2d(i,j), ! I
     1                                    r_pcp_type_thresh_2d(i,j), ! O
     1                                    t_sfc_k(i,j),              ! I
     1                                    td_sfc_k(i,j),             ! I
     1                                    istat_radar_3dref_orig)    ! I

                    if(.true.)then ! testcode
                        if(r_pcp_type_2d(i,j) .eq. 2.0 .and.
     1                      r_pcp_type_thresh_2d(i,j) .eq. 0.0)then
                            dbz_low_2d(i,j) = ref_base
                        endif
                    endif

                endif

            enddo
            enddo

            I4_elapsed = ishow_timer()

!           Add SAO Drizzle to the thresh sfc precip type (if applicable out
!           of radar range)
            call sao_drizzle_correction(r_pcp_type_thresh_2d,NX_L,NY_L
     1              ,n_obs_pos_b,obstype,wx_s,lat_s,lon_s,maxstns
     1              ,ri_s,rj_s,lat,lon
     1              ,t_sfc_k
     1              ,cloud_ceiling,r_missing_data)

            I4_elapsed = ishow_timer()

!           Add SAO Rain to the thresh sfc precip type (if applicable out
!           of radar range)
            call sao_rain_correction(r_pcp_type_thresh_2d,NX_L,NY_L
     1              ,n_obs_pos_b,obstype,wx_s,lat_s,lon_s,maxstns
     1              ,ri_s,rj_s,lat,lon
     1              ,t_sfc_k,td_sfc_k,twet_snow
     1              ,dbz_low_2d
     1              ,cvr_max,r_missing_data)

            I4_elapsed = ishow_timer()

!           Add SAO Snow to the thresh sfc precip type (if applicable out
!           of radar range)
            call sao_snow_correction(r_pcp_type_thresh_2d,NX_L,NY_L
     1              ,n_obs_pos_b,obstype,wx_s,lat_s,lon_s,maxstns
     1              ,ri_s,rj_s,lat,lon
     1              ,t_sfc_k,td_sfc_k,twet_snow
     1              ,dbz_low_2d
     1              ,cvr_max,r_missing_data)

            I4_elapsed = ishow_timer()

!           Increase radar reflectivity threshold for precip in those areas
!           where the radar has echo, but the SAO says no precip.
            call sao_precip_correction(r_pcp_type_thresh_2d,NX_L,NY_L
     1              ,n_obs_pos_b,obstype,wx_s,lat_s,lon_s,maxstns
     1              ,ri_s,rj_s,lat,lon
     1              ,t_sfc_k,td_sfc_k
     1              ,dbz_low_2d
     1              ,cvr_max,r_missing_data)

            I4_elapsed = ishow_timer()

!           Write LMT/LLR
!           Note that these arrays start off with 1 as the first index
            var_a(1) = 'LMT'
            var_a(2) = 'LLR'
            ext = 'lmt'
            units_a(1) = 'M'
            units_a(2) = 'DBZ'
            comment_a(1) = 'LAPS Maximum Tops'
            comment_a(2) = 'LAPS Low Level Reflectivity'

            call move(ref_mt_2d, out_array_3d(1,1,1),NX_L,NY_L)
            call move(dbz_low_2d,out_array_3d(1,1,2),NX_L,NY_L)

            call put_laps_multi_2d(i4time,ext,var_a,units_a,
     1          comment_a,out_array_3d,NX_L,NY_L,2,istatus)

            if(istatus .eq. 1)then
                j_status(n_lmt) = ss_normal
                write(6,*)' Success in writing out LMT'
            else
                write(6,*)' Error detected writing out LMT'
            endif


!           Compare precip type to the obs
            iarg = 0
            if(istat_radar_3dref .eq. 1 .and. istat_sfc .eq. 1)then
              call make_fnam_lp(i4time,asc_tim_9,istatus)
              write(6,*)' ',asc_tim_9,' __'
              write(6,*)' Comparing precip type to the obs:',n_obs_pos_b
              write(6,*)' Sta PTY PTT dbz   T(anl)  Td(anl)'//
     1                  '  T(ob)  Td(ob)  PTY(ob)'//
     1                  '        Tw(anl) Tw(ob) Elev(ob) Celg'
              do i = 1,n_obs_pos_b
                call latlon_to_rlapsgrid(lat_s(i),lon_s(i),lat,lon
     1                     ,NX_L,NY_L,ri,rj,istatus)

                i_i = nint(ri)
                i_j = nint(rj)

                if(i_i .ge. 1 .and. i_i .le. NX_L .and.
     1             i_j .ge. 1 .and. i_j .le. NY_L            )then
                    i_pty = i2_pcp_type_2d(i_i,i_j) / 16
                    i_ptt = r_pcp_type_thresh_2d(i_i,i_j)

!                   Does LAPS analyze or does the SAO report wx here
                    if(  (i_pty .gt. 0 .or. wx_s(i) .ne. '        ')

!                        and is this an SAO station that reports precip?
!    1                  .AND.  obstype(i)(1:4) .ne. 'MESO'
!    1                  .AND.  obstype(i)(1:4) .ne. 'CDOT'
!    1                  .AND.  obstype(i)(7:8) .ne. '1A'
     1                  .AND.  wx_s(i)(1:7)    .ne. 'UNKNOWN'
     1                                                         )then
                        c3_pt_flag = ' __'
                    else
                        c3_pt_flag = '   '
                    endif

                    t_sfc_c  = k_to_c(t_sfc_k(i_i,i_j))
                    td_sfc_c = k_to_c(td_sfc_k(i_i,i_j))
                    p_sfc_mb = pres_sfc_pa(i_i,i_j) / 100.

                    tw_sfc_c = tw(t_sfc_c,td_sfc_c,p_sfc_mb)

                    t_s_c  = f_to_c(t_s(i))
                    td_s_c = f_to_c(td_s(i))

                    if(t_s_c .gt. -50. .and. td_s_c .gt. -50.)then
                        tw_s_c = tw(t_s_c,td_s_c,p_sfc_mb)
                    else
                        tw_s_c = -99.
                    endif

                    if(cloud_ceiling(i_i,i_j) .ne. r_missing_data)
     1                                                     then
                        iceil = nint(cloud_ceiling(i_i,i_j))
                    else
                        iceil = 99999
                    endif

!                   SNOW
                    call parse_wx_pcp(wx_s(i),'S',ipresent,istatus)       
                    if(ipresent .eq. 1)then
                        c1_s = 'S'
                    else
                        c1_s = ' '
                    endif

!                   RAIN
                    call parse_wx_pcp(wx_s(i),'R',ipresent,istatus)       
                    if(ipresent .eq. 1)then
                        c1_r = 'R'
                    else
                        c1_r = ' '
                    endif

                    call filter_string(wx_s(i))

                    write(6,1101,err=1102)c_stations(i)(1:3)
     1                          ,c2_precip_types(i_pty)
     1                          ,c2_precip_types(i_ptt)
     1                          ,int(dbz_low_2d(i_i,i_j))
     1                          ,t_sfc_c
     1                          ,td_sfc_c
     1                          ,t_s_c
     1                          ,td_s_c
     1                          ,wx_s(i)
     1                          ,c3_pt_flag
     1                          ,tw_sfc_c
     1                          ,tw_s_c
     1                          ,elev_s(i)
     1                          ,iceil
     1                          ,cvr_max(i_i,i_j)
     1                          ,c1_r,c1_s,obstype(i)(7:8)
1101                format(1x,a3,2x,a2,2x,a2,i4,4f8.1,3x,a8,2x,a3
     1                                ,3f8.1,i7,f5.2,1x,2a1,1x,a2
     1                                ,' ptvrf')
1102            endif ! ob is in domain
              enddo ! i
            endif

            istat_pty = 1

        else ! set SFC preciptype to missing
            write(6,*)
     1      ' No SFC preciptype calculated due to lack of radar data'       
            do i = 1,NX_L
            do j = 1,NY_L
                r_pcp_type_thresh_2d(i,j) = r_missing_data
                r_pcp_type_2d(i,j) = r_missing_data
            enddo ! j
            enddo ! i

        endif ! istat_radar_3dref

 700    continue

        if(istat_radar_3dref .eq. 1)then
            if(l_flag_bogus_w .and. l_bogus_radar_w) then    ! Adan add
!             Re-calculate cloud bogus omega within radar echo area
!             Add by Adan
              call get_radar_deriv(NX_L,NY_L,NZ_L,grid_spacing_cen_m,       
     1                           r_missing_data,
     1                           radar_ref_3d,clouds_3d,cld_hts,
     1                           temp_3d,heights_3d,pres_3d,
     1                           ibase_array,itop_array,thresh_cvr,
     1                           vv_to_height_ratio_Cu,                 ! I
     1                           cldpcp_type_3d,w_3d,istat_radar_deriv)       
              if(istat_radar_deriv .ne. 1)then
                write(6,*)' Bad status return from get_radar_deriv'
              endif
            endif                           ! l_flag_bogus_w (Adan add)

            I4_elapsed = ishow_timer()

            write(6,*)' Computing Precip Concentration'

!           Calculate 3D Precip Concentration in kg/m**3
            call cpt_pcp_cnc(radar_ref_3d,temp_3d           ! Input
     1                                  ,rh_3d_pct          ! Input
     1                                  ,cldpcp_type_3d     ! Input
     1                                  ,NX_L,NY_L,NZ_L     ! Input
     1                                  ,c_z2m              ! Input
     1                                  ,pres_3d            ! Input
     1                                  ,pcpcnc             ! Output
     1                                  ,raicnc             ! Output
     1                                  ,snocnc             ! Output
     1                                  ,piccnc)            ! Output

!           Calculate Integrated Rainwater
            write(6,*)
            write(6,*)' Calculating Integrated Rainwater'
            call integrate_slwc(raicnc,heights_3d,NX_L,NY_L,NZ_L
     1                                                     ,rain_int)
            write(6,*)' Integrated rain range is ',minval(rain_int)
     1                                           ,maxval(rain_int)

            I4_elapsed = ishow_timer()

!           Calculate Integrated Snow
            write(6,*)
            write(6,*)' Calculating Integrated Snow'
            call integrate_slwc(snocnc,heights_3d,NX_L,NY_L,NZ_L
     1                                                     ,snow_int)
            write(6,*)' Integrated snow range is ',minval(snow_int)
     1                                            ,maxval(snow_int)

            I4_elapsed = ishow_timer()

!           Calculate Integrated Precipitating Ice
            write(6,*)
            write(6,*)' Calculating Integrated Precip. Ice'
            call integrate_slwc(piccnc,heights_3d,NX_L,NY_L,NZ_L
     1                                                     ,pice_int)
            write(6,*)' Integrated pice range is ',minval(pice_int)
     1                                            ,maxval(pice_int)

            I4_elapsed = ishow_timer()

            nf = 6

        else
            nf = 2

        endif

!       Rain, Snow, and Graupel are added for the cldalb & simvis computation
        do j = 1,NY_L
        do i = 1,NX_L
            cldalb_out(i,j) 
     1                  = 1. - (exp( -(const_lwp_bks * slwc_int(i,j) + 
     1                                 const_iwp_bks * cice_int(i,j) + 
     1                                 const_rwp_bks * rain_int(i,j) + 
     1                                 const_swp_bks * snow_int(i,j) + 
     1                                 const_gwp_bks * pice_int(i,j)) ))
            cldod_out(i,j) = const_lwp * slwc_int(i,j)
     1                     + const_iwp * cice_int(i,j)  
     1                     + const_rwp * rain_int(i,j)  
     1                     + const_swp * snow_int(i,j)  
     1                     + const_gwp * pice_int(i,j)  
            if(i .eq. NX_L/2 .AND. j .eq. NY_L/2)then
                write(6,*)'i,j,lat,lon,cld_od:'
     1                    ,i,j,lat(i,j),lon(i,j),cldod_out(i,j)
            endif
        enddo ! i
        enddo ! j

!       Write LIL file
!       Note that these arrays start off with 1 as the first index
        ext = 'lil'
        var_a(1) = 'LIL'
        var_a(2) = 'LIC'
        var_a(3) = 'COD'
        var_a(4) = 'CLA'
        units_a(1) = 'M'
        units_a(2) = 'M'
        units_a(3) = ' '
        units_a(4) = ' '
        comment_a(1) = 'Integrated Cloud Liquid'
        comment_a(2) = 'Integrated Cloud Ice'
        comment_a(3) = 'Cloud Optical Depth'
        comment_a(4) = 'Cloud Albedo'

        call move(slwc_int,  out_array_3d(1,1,1),NX_L,NY_L)
        call move(cice_int,  out_array_3d(1,1,2),NX_L,NY_L)
        call move(cldod_out, out_array_3d(1,1,3),NX_L,NY_L)
        call move(cldalb_out,out_array_3d(1,1,4),NX_L,NY_L)

        call put_laps_multi_2d(i4time,ext,var_a,units_a,
     1      comment_a,out_array_3d,NX_L,NY_L,4,istatus)

        if(istatus .eq. 1)then
            j_status(n_lil) = ss_normal
            write(6,*)' Success in writing out LIL'
        else
            write(6,*)' Error detected writing out LIL'
        endif

!       Convert SLWC and CICE from g/m**3 to kg/m**3
        do k = 1,NZ_L
        do j = 1,NY_L
        do i = 1,NX_L
            if(slwc(i,j,k) .ne. r_missing_data)then
                slwc(i,j,k) = (slwc(i,j,k) * 1.0) / 1e3
            endif
            if(cice(i,j,k) .ne. r_missing_data)then
                cice(i,j,k) = (cice(i,j,k) * 1.0) / 1e3
            endif
        enddo 
        enddo
        enddo

!       Apply hydrometeor scale to precip
        if(hydrometeor_scale_pcp .ge. 0.)then
            ratio_pcp =  hydrometeor_scale_pcp
        else
            ratio_pcp = -hydrometeor_scale_pcp / 
     1                  (grid_spacing_cen_m/1000.)
        endif

        pcpcnc = pcpcnc * ratio_pcp
        raicnc = raicnc * ratio_pcp
        snocnc = snocnc * ratio_pcp
        piccnc = piccnc * ratio_pcp

!       Write out Cloud Liquid Water, Cloud Ice and Precip Content Fields
        var_a(1) = 'LWC'
        var_a(2) = 'ICE'
        var_a(3) = 'PCN'
        var_a(4) = 'RAI'
        var_a(5) = 'SNO'
        var_a(6) = 'PIC'

        ext = 'lwc'

        units_a(1) = 'KG/M**3'
        units_a(2) = 'KG/M**3'
        units_a(3) = 'KG/M**3'
        units_a(4) = 'KG/M**3'
        units_a(5) = 'KG/M**3'
        units_a(6) = 'KG/M**3'

        comment_a(1) = 'Cloud Liquid Water Content - LAPS Smith Feddes'       
        comment_a(2) = 'Cloud Ice Content - LAPS Smith Feddes'
        comment_a(3) = 'Precipitate Content'
        comment_a(4) = 'Rain Content'
        comment_a(5) = 'Snow Content'
        comment_a(6) = 'Precipitating Ice Content'

        slwc_max = maxval(slwc)
        cice_max = maxval(cice)
        pcpcnc_max = maxval(pcpcnc)
        raicnc_max = maxval(raicnc)
        snocnc_max = maxval(snocnc)
        piccnc_max = maxval(piccnc)

        write(6,*)' Max slwc = ',slwc_max
        write(6,*)' Max cice = ',cice_max
        write(6,*)' Max pcpcnc = ',pcpcnc_max
        write(6,*)' Max raicnc = ',raicnc_max
        write(6,*)' Max snocnc = ',snocnc_max
        write(6,*)' Max piccnc = ',piccnc_max

        if(slwc_max   .le. 1e6 .AND. cice_max   .le. 1e6 .AND.
     1     pcpcnc_max .le. 1e6 .AND. raicnc_max .le. 1e6 .AND.
     1     snocnc_max .le. 1e6 .AND. piccnc_max .le. 1e6       )then

            call put_laps_3d_multi(i4time,ext,var_a,units_a,comment_a
     1                        ,slwc,cice
     1                        ,pcpcnc,raicnc
     1                        ,snocnc,piccnc
     1                        ,NX_L,NY_L,NZ_L
     1                        ,NX_L,NY_L,NZ_L
     1                        ,NX_L,NY_L,NZ_L
     1                        ,NX_L,NY_L,NZ_L
     1                        ,NX_L,NY_L,NZ_L
     1                        ,NX_L,NY_L,NZ_L
     1                        ,nf,istatus)
            if(istatus .eq. 1)j_status(n_lwc) = ss_normal
        else
            write(6,*)' ERROR: large hydrometeor values LWC not written'
        endif

!       Write out Mean Volume Diameter field (Potentially compressible)
        ext = 'lmd'
        var = 'LMD'
        units = 'M'
        comment = 'Mean Volume Diameter of Cloud Droplets'
        call put_laps_3d(i4time,ext,var,units,comment,mvd_3d
     1                                           ,NX_L,NY_L,NZ_L)
        j_status(n_lmd) = ss_normal

!       Write out Icing Index field  (Potentially compressible)
        ext = 'lrp'
        do k = 1,NZ_L
        do j = 1,NY_L
        do i = 1,NX_L
            out_array_3d(i,j,k) = icing_index_3d(i,j,k)
        enddo
        enddo
        enddo

        var = 'LRP'
        units = 'NONE'
        comment = 'Icing Severity Index '//
     1  '1-LtCnt,2-MdCnt,3-HvCnt,4-LtInt,5-MdInt,6-HvInt'
        call put_laps_3d(i4time,ext,var,units,comment,out_array_3d
     1                                  ,NX_L,NY_L,NZ_L)

        j_status(n_lrp) = ss_normal

!       Write 3D Precip Type
!       4 most significant bits are precip type, other 4 are cloud type
        do k = 1,NZ_L
        do j = 1,NY_L
        do i = 1,NX_L
            iarg = cldpcp_type_3d(i,j,k)
            out_array_3d(i,j,k) = iarg/16                   ! 'PTY'
        enddo
        enddo
        enddo

        ext = 'pty'
        var = 'PTY'
        units = 'NONE'
        comment = 'Precip Type: 0-None,1-Rain,2-Snow,3-ZR,4-IP,5-Hail'
        call put_laps_3d(i4time,ext,var,units
     1                  ,comment,out_array_3d(1,1,1)
     1                  ,NX_L,NY_L,NZ_L)

        I4_elapsed = ishow_timer()

!       Write sfc precip and cloud type
        var_a(1) = 'PTY'
        var_a(2) = 'PTT'
        var_a(3) = 'SCT'
        ext = 'lct'
        units_a(1) = 'UNDIM'
        units_a(2) = 'UNDIM'
        units_a(3) = 'UNDIM'
        comment_a(1) = 'LAPS Precip Type (Unthresholded): '//
     1                 '0:Np 1:Rn 2:Sn 3:Zr 4:Sl 5:Ha 6:L  7:ZL'
        comment_a(2) = 'LAPS Precip Type (Refl Threshold): '//
     1                 '0:Np 1:Rn 2:Sn 3:Zr 4:Sl 5:Ha 6:L  7:ZL'
        comment_a(3) = 'LAPS Cloud Type '//
     1                 '0:Clr 1:St 2:Sc 3:Cu 4:Ns 5:Ac '//
     1                 '6:As 7:Cs 8:Ci 9:Cc 10: Cb'

        call move(r_pcp_type_2d,       out_array_3d(1,1,1),NX_L,NY_L)
        call move(r_pcp_type_thresh_2d,out_array_3d(1,1,2),NX_L,NY_L)
        call move(r_cld_type_2d,       out_array_3d(1,1,3),NX_L,NY_L)        

        call put_laps_multi_2d(i4time,ext,var_a,units_a,
     1                 comment_a,out_array_3d,NX_L,NY_L,3,istatus)

        if(istatus .eq. 1)j_status(n_lct) = ss_normal

!       Write out Cloud derived Omega field
        var = 'COM'
        ext = 'lco'
        units = 'PA/S'
        comment = 'LAPS Cloud Derived Omega'
        call put_laps_3d(i4time,ext,var,units,comment,w_3d
     1                  ,NX_L,NY_L,NZ_L)
        j_status(n_lco) = ss_normal

        I4_elapsed = ishow_timer()

        istatus = 1

!       Read sounding metadata to get integrated cloud liquid obs
!       lun = 89
!       ext = 'snd'
!       call read_snd_metadata(lun,i4time,ext                         ! I
!    1                        ,MAX_SND_GRID,MAX_SND_LEVELS            ! I
!    1                        ,lat,lon,imax,jmax                      ! I
!    1                        ,n_profiles                             ! O
!    1                        ,nlevels_obs_pr,lat_pr,lon_pr,elev_pr   ! O
!    1                        ,c5_name,i4time_ob_pr,obstype           ! O
!    1                        ,cloud_base_temp,cloud_liquid           ! O
!    1                        ,istatus)                               ! O

!       write(6,*)' back from read_snd_metadata, # soundings = '
!    1            ,n_profiles

999     continue

        write(6,*)' Notifications'
        do i = iprod_start,iprod_end
            write(6,*)' ',exts(i),' ',j_status(i),' ',i
        enddo ! i

9999    deallocate( slwc )
        deallocate( cice )

        write(6,*)' End of laps_deriv_sub'

        return
        end


