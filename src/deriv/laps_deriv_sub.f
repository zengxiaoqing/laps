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

        subroutine laps_deriv_sub(i4time,
     1                        NX_L,NY_L,
     1                        NZ_L,
     1                        N_PIREP,
     1                        maxstns,
     1                        max_cld_snd,
     1                        n_prods,
     1                        iprod_number,
     1                        temp_3d,
     1                        heights_3d,
     1                        rh_3d_pct,
     1                        pres_sfc_pa,
     1                        t_sfc_k,
     1                        j_status,istatus)

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
!       1997 Jul 31 K. Dritz  - Added N_PIREP, maxstns, and max_cld_snd as
!                               dummy arguments.  Removed the PARAMETER
!                               statement for max_cld_snd.
!       1997 Jul 31 K. Dritz  - Added call to get_ref_base.
!       1997 Aug 01 K. Dritz  - Added maxstns, IX_LOW, IX_HIGH, IY_LOW, and
!                               IY_HIGH as arguments in call to insert_sao.
!       1997 Aug 01 K. Dritz  - Also now pass r_missing_data to barnes_r5.
!       1997 Aug 01 K. Dritz  - Pass r_missing_data to insert_sat.
!       1997 Aug 01 K. Dritz  - Pass ref_base to rfill_evap.

!       Prevents clearing out using satellite (hence letting SAOs dominate)
!       below this altitude (M AGL)
        real*4 surface_sao_buffer
        parameter (surface_sao_buffer = 800.)

        real*4 thresh_cvr,default_top,default_base,default_clear_cover
     1                   ,default_ceiling

        parameter       (thresh_cvr = 0.65) ! Used to "binaryize" cloud cover

        parameter       (default_clear_cover = .01)

        real*4 thresh_cvr_base,thresh_cvr_top,thresh_cvr_ceiling
        parameter (thresh_cvr_base = 0.1)
        parameter (thresh_cvr_top  = 0.1)
        parameter (thresh_cvr_ceiling = thresh_cvr)

        real*4 thresh_thin_lwc_ice     ! Threshold cover for thin cloud LWC/ICE
        parameter (thresh_thin_lwc_ice = 0.1)

        real*4 vis_radar_thresh_cvr,vis_radar_thresh_dbz
        parameter (vis_radar_thresh_cvr = 0.2)  ! 0.2, 0.0
        parameter (vis_radar_thresh_dbz = 10.)  ! 5. , -99.

        real*4 lat(NX_L,NY_L),lon(NX_L,NY_L)
        real*4 topo(NX_L,NY_L)
        real*4 rlaps_land_frac(NX_L,NY_L)
        real*4 solar_alt(NX_L,NY_L)
        real*4 solar_ha(NX_L,NY_L)

        logical l_packed_output
        logical l_evap_radar

        data l_packed_output /.false./
        data l_evap_radar /.false./

        logical l_fill
        logical l_flag_mvd
        logical l_flag_cloud_type
        logical l_flag_icing_index
        logical l_flag_bogus_w
        logical l_flag_snow_potential
        logical l_parse

        logical l_sao_lso
        data l_sao_lso /.true./ ! Do things the new way?

        logical l_perimeter
        data l_perimeter /.true./ ! Use SAOs just outside domain?

        include 'laps_cloud.inc'

!       Nominal cloud heights. Actual ones used are fitted to the terrain.
        real*4 cld_hts_new(KCLOUD)

        data cld_hts_new/1200.,1300.,1400.,1500.,1600.,1700.,1800.,
     11900.,2000.,2100.,2200.,2400.,2600.,2800.,3000.,3200.,
     23400.,3600.,3800.,4000.,4500.,5000.,5500.,6000.,6500.,
     37000.,7500.,8000.,8500.,9000.,9500.,10000.,11000.,12000.,
     413000.,14000.,15000.,16000.,17000.,18000.,19000.,20000./

        integer iarg

        equivalence (cld_hts,cld_hts_new)

        real*4 clouds_3d(NX_L,NY_L,KCLOUD)

        integer ihist_alb(-10:20)

        real*4 cloud_top(NX_L,NY_L)
        real*4 cloud_base(NX_L,NY_L)
        real*4 cloud_ceiling(NX_L,NY_L)
        real*4 cloud_base_buf(NX_L,NY_L)

        real*4 cldtop_m(NX_L,NY_L)
        real*4 cldtop_m_co2(NX_L,NY_L)
        real*4 cldtop_m_tb8(NX_L,NY_L)

        real*4 cld_pres_1d(KCLOUD)
        real*4 pres_3d(NX_L,NY_L,NZ_L)
        real*4 clouds_3d_pres(NX_L,NY_L,NZ_L)

        real*4 CVHZ(NX_L,NY_L)
        real*4 CVHZ1(NX_L,NY_L),CVEW1(NX_L,KCLOUD)
        real*4 cvr_max(NX_L,NY_L),CVEW2(NX_L,KCLOUD)
        real*4 cvr_sao_max(NX_L,NY_L)
        real*4 cvr_snow_cycle(NX_L,NY_L)
        real*4 cvr_water_temp(NX_L,NY_L)
        real*4 cvr_snow(NX_L,NY_L)
        real*4 band8_mask(NX_L,NY_L)

        character*4 radar_name
        character*31 radarext_3d_cloud
        real*4 radar_ref_3d(NX_L,NY_L,NZ_L)
        real*4 heights_3d(NX_L,NY_L,NZ_L)

        real*4 mvd_3d(NX_L,NY_L,NZ_L)
!       real*4 lwc_res_3d(NX_L,NY_L,NZ_L)
        real*4 w_3d(NX_L,NY_L,NZ_L)

        integer icing_index_3d(NX_L,NY_L,NZ_L)

        integer cldpcp_type_3d(NX_L,NY_L,NZ_L) ! Also contains 3D precip type

!       Output array declarations
!       real*4 out_array_4d(NX_L,NY_L,NZ_L,2)
        real*4, allocatable, dimension(:,:,:,:) :: out_array_4d
        real*4 out_array_3d(NX_L,NY_L,NZ_L)

        real*4 slwc(NX_L,NY_L,NZ_L),slwc_int(NX_L,NY_L)
        real*4 cice(NX_L,NY_L,NZ_L)

        real*4 pcpcnc(NX_L,NY_L,NZ_L)
        real*4 raicnc(NX_L,NY_L,NZ_L)
        real*4 snocnc(NX_L,NY_L,NZ_L)
        real*4 piccnc(NX_L,NY_L,NZ_L)

!       real*4 snow_2d(NX_L,NY_L)

        character*2 c2_precip_types(0:10)

        data c2_precip_types
     1  /'  ','Rn','Sn','Zr','Sl','Ha','L ','ZL','  ','  ','  '/

        character*3 c3_pt_flag
        character*1 c1_r,c1_s
        character*8 c8_project

        integer i2_pcp_type_2d(NX_L,NY_L)
        real*4 r_pcp_type_2d(NX_L,NY_L)

        real*4 dum1_array(NX_L,NY_L)
        real*4 dum2_array(NX_L,NY_L)
        real*4 dum3_array(NX_L,NY_L)
        real*4 dum4_array(NX_L,NY_L)

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
        real*4 tb8_k(NX_L,NY_L)
        real*4 tb8_cold_k(NX_L,NY_L)
        real*4 albedo(NX_L,NY_L)
        real*4 cloud_frac_vis_a(NX_L,NY_L)
        real*4 cloud_frac_co2_a(NX_L,NY_L)

        real*4 temp_3d(NX_L,NY_L,NZ_L)
        real*4 rh_3d_pct(NX_L,NY_L,NZ_L)
        real*4 model_3d(NX_L,NY_L,NZ_L)

        real*4 t_sfc_k(NX_L,NY_L)
        real*4 t_gnd_k(NX_L,NY_L)
        real*4 sst_k(NX_L,NY_L)
        real*4 td_sfc_k(NX_L,NY_L)
        real*4 pres_sfc_pa(NX_L,NY_L)

!       Declarations for LSO file stuff
        real*4 lat_s(maxstns), lon_s(maxstns), elev_s(maxstns)
        real*4 cover_s(maxstns)
        real*4 t_s(maxstns), td_s(maxstns), pr_s(maxstns), sr_s(maxstns)
        real*4 dd_s(maxstns), ff_s(maxstns)
        real*4 ddg_s(maxstns), ffg_s(maxstns)
        real*4 vis_s(maxstns)
        real*4 pstn_s(maxstns),pmsl_s(maxstns),alt_s(maxstns)
        real*4 store_hgt(maxstns,5),ceil(maxstns),lowcld(maxstns)
        real*4 cover_a(maxstns),rad_s(maxstns)
        integer obstime(maxstns),kloud(maxstns),idp3(maxstns)
        character store_emv(maxstns,5)*1,store_amt(maxstns,5)*4
        character wx_s(maxstns)*8, obstype(maxstns)*8
        character atime*24, infile*70

        integer STATION_NAME_LEN
        parameter (STATION_NAME_LEN = 3)                   
        character c_stations(maxstns)*(STATION_NAME_LEN)    

        character asc_tim_9*9

        real*4 ri_s(maxstns), rj_s(maxstns)

!       Product # notification declarations
        integer j_status(20),iprod_number(20)

!       Stuff for 2d fields
        real*4 ref_mt_2d(NX_L,NY_L)
        real*4 dbz_low_2d(NX_L,NY_L)
        real*4 dbz_max_2d(NX_L,NY_L)

!       SFC precip and cloud type (LCT file)
        real*4 r_pcp_type_thresh_2d(NX_L,NY_L)
        real*4 r_cld_type_2d(NX_L,NY_L)

        character*40 c_vars_req
        character*180 c_values_req

        character*3 lso_ext
        data lso_ext /'lso'/

        ISTAT = INIT_TIMER()

        write(6,*)' Welcome to the LAPS gridded cloud analysis'

        call get_r_missing_data(r_missing_data,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'Error calling get_r_missing_data'
           stop
        endif

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
            return
        endif

!       radarext_3d_cloud = radarext_3d
!       radarext_3d_cloud = 'v02'

        write(6,*)' radarext_3d_cloud = ',radarext_3d_cloud

c read in laps lat/lon and topo
        call get_laps_domain_95(NX_L,NY_L,lat,lon,topo
     1           ,rlaps_land_frac,grid_spacing_cen_m,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error getting LAPS domain'
            return
        endif

        call get_laps_cycle_time(ilaps_cycle_time,istatus)
        if(istatus .eq. 1)then
            write(6,*)' ilaps_cycle_time = ',ilaps_cycle_time
        else
            write(6,*)' Error getting laps_cycle_time'
            return
        endif

        call get_pres_3d(i4time,NX_L,NY_L,NZ_L,pres_3d,istatus)
        if(istatus .ne. 1)return

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
                return
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

            endif ! istatus_lps

            var = 'LCV'
            ext = 'lcv'
            call get_laps_2d(i4time,ext,var,units,comment
     1                  ,NX_L,NY_L,cvr_max,istatus)

            if(istatus .ne. 1)THEN
                write(6,*)' Error Reading Cvr_max Analysis - abort'
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
            if(istatus .ne. 1)return

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
        l_flag_snow_potential = .true.

        call get_cloud_deriv(
     1                NX_L,NY_L,NZ_L,clouds_3d,cld_hts,
     1                temp_3d,rh_3d_pct,heights_3d,pres_3d,
     1                istat_radar_3dref,radar_ref_3d,grid_spacing_cen_m,       
     1                l_mask_pcptype,ibase_array,itop_array,
     1                iflag_slwc,slwc,cice,thresh_cvr,
     1                l_flag_cloud_type,cldpcp_type_3d,
     1                l_flag_mvd,mvd_3d,
     1                l_flag_icing_index,icing_index_3d,
     1                l_flag_bogus_w,w_3d,istatus)
!    1                l_flag_snow_potential,snow_2d,lwc_res_3d)

        if(istatus .ne. 1)then
            write(6,*)' Bad status return from get_cloud_deriv'
            return
        endif

        I4_elapsed = ishow_timer()

        write(6,*)
        write(6,*)' Inserting thin clouds into LWC/ICE fields'
        call insert_thin_lwc_ice(clouds_3d,clouds_3d_pres,heights_3d
     1       ,temp_3d,cld_hts,NX_L,NY_L,NZ_L,KCLOUD,thresh_thin_lwc_ice       
     1       ,pres_3d,slwc,cice,istatus)

        if(istatus .ne. 1)then
            write(6,*)' Bad status return from insert_thin_lwc_ice'
            return
        endif

        I4_elapsed = ishow_timer()

!       Calculate and Write Integrated LWC
        write(6,*)
        write(6,*)' Calculating Integrated LWC'
        call integrate_slwc(slwc,heights_3d,NX_L,NY_L,NZ_L,slwc_int)

        var = 'LIL'
        ext = 'lil'
        units = 'M'
        comment = 'no comment'
        call put_laps_2d(i4time,ext,var,units,comment
     1  ,NX_L,NY_L,slwc_int,istatus)
        if(istatus .eq. 1)j_status(n_lil) = ss_normal
        I4_elapsed = ishow_timer()

!       DERIVED RADAR/PRECIP STUFF
        if(istat_radar_3dref .eq. 1)then ! LMT

            if(l_evap_radar 
!    1             .and. istatus_rh .eq. 1
     1             .and. istat_radar_3dref_orig .eq. 1)then ! Reread radar data

                write(6,*)
     1              ' Rereading radar data to apply evaporation/fill'

!               Get time of radar file of the indicated appropriate extension
                call get_filespec(radarext_3d_cloud(1:3),2,c_filespec
     1                                                    ,istatus)
                call get_file_time(c_filespec,i4time,i4time_radar)

                call read_radar_3dref(i4time_radar,
!    1                 1200,i4time_radar,
     1                 .true.,ref_base,                                 ! I
     1                 NX_L,NY_L,NZ_L,radarext_3d_cloud,
     1                 lat,lon,topo,.true.,.true.,
     1                 heights_3d,
     1                 radar_ref_3d,
     1                 rlat_radar,rlon_radar,rheight_radar,radar_name,     
     1                 n_ref_grids,istat_radar_2dref,istat_radar_3dref)       

                I4_elapsed = ishow_timer()

!               Apply evaporation subroutine
                call rfill_evap(radar_ref_3d,NX_L,NY_L,NZ_L
     1          ,.true.,.true.
     1          ,lat,lon,topo,rlat_radar,rlon_radar,rheight_radar
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


!           Note that pres_sfc_pa was already read in above
            call get_sfc_preciptype(pres_sfc_pa,t_sfc_k,td_sfc_k
     1                             ,cldpcp_type_3d
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

!           Add SAO Drizzle to the thresh sfc precip type (if applicable out
!           of radar range)
            call sao_drizzle_correction(r_pcp_type_thresh_2d,NX_L,NY_L
     1              ,n_obs_pos_b,obstype,wx_s,lat_s,lon_s,maxstns
     1              ,ri_s,rj_s,lat,lon
     1              ,t_sfc_k
     1              ,cloud_ceiling,r_missing_data)

!           Add SAO Rain to the thresh sfc precip type (if applicable out
!           of radar range)
            call sao_rain_correction(r_pcp_type_thresh_2d,NX_L,NY_L
     1              ,n_obs_pos_b,obstype,wx_s,lat_s,lon_s,maxstns
     1              ,ri_s,rj_s,lat,lon
     1              ,t_sfc_k,td_sfc_k
     1              ,dbz_low_2d
     1              ,cvr_max,r_missing_data)

!           Add SAO Snow to the thresh sfc precip type (if applicable out
!           of radar range)
            call sao_snow_correction(r_pcp_type_thresh_2d,NX_L,NY_L
     1              ,n_obs_pos_b,obstype,wx_s,lat_s,lon_s,maxstns
     1              ,ri_s,rj_s,lat,lon
     1              ,t_sfc_k,td_sfc_k
     1              ,dbz_low_2d
     1              ,cvr_max,r_missing_data)


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
     1                  '        Tw(anl)  Tw(ob)  Celg'
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

                    t_sfc_c  = t_sfc_k(i_i,i_j)  - 273.15    ! K to C
                    td_sfc_c = td_sfc_k(i_i,i_j) - 273.15    ! K to C
                    p_sfc_mb = pres_sfc_pa(i_i,i_j) / 100.

                    tw_sfc_c = tw(t_sfc_c,td_sfc_c,p_sfc_mb)

                    t_s_c  = (t_s(i)-32.) / 1.8              ! F to C
                    td_s_c = (td_s(i)-32.) / 1.8             ! F to C

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
     1                          ,iceil
     1                          ,cvr_max(i_i,i_j)
     1                          ,c1_r,c1_s,obstype(i)(7:8)
1101                format(1x,a3,2x,a2,2x,a2,i4,4f8.1,3x,a8,2x,a3
     1                                ,2f8.1,i7,f5.2,1x,2a1,1x,a2)
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

        if(istat_radar_3dref .eq. 1)then

            I4_elapsed = ishow_timer()

            write(6,*)' Computing Precip Concentration'

!           Calculate 3D Precip Concentration in kg/m**3
            call cpt_pcp_cnc(radar_ref_3d,temp_3d,cldpcp_type_3d    ! Input
     1                                  ,NX_L,NY_L,NZ_L     ! Input
     1                                  ,pcpcnc             ! Output
     1                                  ,raicnc             ! Output
     1                                  ,snocnc             ! Output
     1                                  ,piccnc)            ! Output

            I4_elapsed = ishow_timer()

            nf = 6

        else
            nf = 2

        endif

!       Convert SLWC and CICE from g/m**3 to kg/m**3
 700    do k = 1,NZ_L
        do j = 1,NY_L
        do i = 1,NX_L
            if(slwc(i,j,k) .ne. r_missing_data)then
                slwc(i,j,k) = slwc(i,j,k) / 1e3
            endif
            if(cice(i,j,k) .ne. r_missing_data)then
                cice(i,j,k) = cice(i,j,k) / 1e3
            endif
        enddo 
        enddo
        enddo

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

!       call move_3d(slwc(1,1,1)  ,out_array_4d(1,1,1,1),NX_L,NY_L,NZ_L)       
!       call move_3d(cice(1,1,1)  ,out_array_4d(1,1,1,2),NX_L,NY_L,NZ_L)
!       call move_3d(pcpcnc(1,1,1),out_array_4d(1,1,1,3),NX_L,NY_L,NZ_L)
!       call move_3d(raicnc(1,1,1),out_array_4d(1,1,1,4),NX_L,NY_L,NZ_L)
!       call move_3d(snocnc(1,1,1),out_array_4d(1,1,1,5),NX_L,NY_L,NZ_L)
!       call move_3d(piccnc(1,1,1),out_array_4d(1,1,1,6),NX_L,NY_L,NZ_L)

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
     1  '1-LtInt,2-MdInt,3-HvInt,4-LtCnt,5-MdCnt,6-HvCnt'
        call put_laps_3d(i4time,ext,var,units,comment,out_array_3d
     1                                  ,NX_L,NY_L,NZ_L)

        j_status(n_lrp) = ss_normal

!       Write out Cloud/3D Precip Type Field
        ext = 'lty'
        var_a(1) = 'PTY'
        var_a(2) = 'CTY'
        units_a(1) = 'NONE'
        units_a(2) = 'NONE'
        comment_a(1) = 'Precip Type: 0-None,1-Rain,2-Snow,3-ZR,4-IP,5-Ha
     1il'
        comment_a(2) = 'Cloud Type: (1-10) - St,Sc,Cu,Ns,Ac,As,Cs,Ci,Cc,
     1Cb'

        allocate( out_array_4d(NX_L,NY_L,NZ_L,2), STAT=istat_alloc )
        if(istat_alloc .ne. 0)then
            write(6,*)' ERROR: Could not allocate out_array_4d'
            stop
        endif

!       4 most significant bits are precip type, other 4 are cloud type
        do k = 1,NZ_L
        do j = 1,NY_L
        do i = 1,NX_L
            iarg = cldpcp_type_3d(i,j,k)
            out_array_4d(i,j,k,1) = iarg/16                   ! 'PTY'
            out_array_4d(i,j,k,2) = iarg - iarg/16*16         ! 'CTY'
        enddo
        enddo
        enddo

        call put_laps_multi_3d(i4time,ext,var_a,units_a,comment_a
     1                  ,out_array_4d,NX_L,NY_L,NZ_L,2,istatus)

        j_status(n_lty) = istatus

        I4_elapsed = ishow_timer()

!       Calculate "SFC" cloud type
!       Now this is simply the cloud type of the lowest significant
!       (> 0.65 cvr) layer present. If a CB is present, it is used instead.
        do i = 1,NX_L
        do j = 1,NY_L
            r_cld_type_2d(i,j) = 0

!           Pick lowest "significant" layer
            do k = NZ_L,1,-1
                cld_type_3d = out_array_4d(i,j,k,2)
                if(cld_type_3d .gt. 0)then
                    r_cld_type_2d(i,j) = cld_type_3d
                endif
            enddo ! k

!           Pick a CB if present
            do k = NZ_L,1,-1
                cld_type_3d = out_array_4d(i,j,k,2)
                if(cld_type_3d .eq. 10)then
                    r_cld_type_2d(i,j) = cld_type_3d
                endif
            enddo ! k

        enddo ! j
        enddo ! i

        deallocate( out_array_4d )

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

999     continue

        write(6,*)' Notifications'
        do i = iprod_start,iprod_end
            write(6,*)' ',exts(i),' ',j_status(i),' ',i
        enddo ! i

        write(6,*)' End of laps_deriv_sub'

        return
        end


