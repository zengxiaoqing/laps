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

        subroutine laps_cloud(i4time,
     1                        NX_L,NY_L,
     1                        NZ_L,
     1                        N_PIREP,
     1                        maxstns,
     1                        max_cld_snd,
     1                        i_diag,
     1                        n_prods,
     1                        iprod_number,
     1                        isplit,
     1                        j_status)

        integer       ss_normal,sys_bad_prod,sys_no_data,
     1                sys_abort_prod

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
!       1997 Jul 31 K. Dritz  - Changed LAPS_DOMAIN_FILE to 'nest7grid'.
!       1997 Jul 31 K. Dritz  - Added call to get_ref_base.
!       1997 Aug 01 K. Dritz  - Compute NX_DIM_LUT, NY_DIM_LUT, and n_fnorm as
!                               they were previously computed in barnes_r5.
!       1997 Aug 01 K. Dritz  - Added NX_DIM_LUT, NY_DIM_LUT, IX_LOW, IX_HIGH,
!                               IY_LOW, IY_HIGH, and n_fnorm as arguments in
!                               call to barnes_r5.
!       1997 Aug 01 K. Dritz  - Added maxstns, IX_LOW, IX_HIGH, IY_LOW, and
!                               IY_HIGH as arguments in call to insert_sao.
!       1997 Aug 01 K. Dritz  - Also now pass r_missing_data to barnes_r5.
!       1997 Aug 01 K. Dritz  - Pass r_missing_data to insert_sat.
!       1997 Aug 01 K. Dritz  - Pass ref_base to rfill_evap.

!       Prevents clearing out using satellite (hence letting SAOs dominate)
!       below this altitude (M AGL)
        real*4 surface_sao_buffer
        parameter (surface_sao_buffer = 800.)

        real*4 thresh_cvr_smf,default_top,default_base
     1                       ,default_clear_cover,default_ceiling

        parameter       (thresh_cvr_smf = 0.65) ! Used to "binaryize" cloud cover

        parameter       (default_clear_cover = .01)

        real*4 thresh_cvr_base,thresh_cvr_top,thresh_cvr_ceiling
        parameter (thresh_cvr_base = 0.1)
        parameter (thresh_cvr_top  = 0.1)
        parameter (thresh_cvr_ceiling = thresh_cvr_smf)

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

        logical l_packed_output, l_use_vis
        logical l_evap_radar

        data l_packed_output /.false./
        data l_evap_radar /.false./

        logical l_fill
        logical l_flag_mvd
        logical l_flag_cloud_type
        logical l_flag_icing_index
        logical l_flag_bogus_w
        logical l_flag_snow_potential

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

        REAL*4 cldcv1(NX_L,NY_L,KCLOUD)
        REAL*4 cf_modelfg(NX_L,NY_L,KCLOUD)
        REAL*4 t_modelfg(NX_L,NY_L,KCLOUD)

        real*4 clouds_3d(NX_L,NY_L,KCLOUD)

        integer ista_snd(max_cld_snd)
        real*4 cld_snd(max_cld_snd,KCLOUD)
        real*4 wt_snd(max_cld_snd,KCLOUD)
        real*4 cvr_snd(max_cld_snd)
        integer i_snd(max_cld_snd)
        integer j_snd(max_cld_snd)

        integer ihist_alb(-10:20)

        real*4 cloud_top(NX_L,NY_L)
        real*4 cloud_base(NX_L,NY_L)
        real*4 cloud_ceiling(NX_L,NY_L)
        real*4 cloud_base_buf(NX_L,NY_L)

        real*4 cldtop_m(NX_L,NY_L)
        real*4 cldtop_m_co2(NX_L,NY_L)
        real*4 cldtop_m_tb8(NX_L,NY_L)

        real*4 cldcv_sao(NX_L,NY_L,KCLOUD)
        real*4 cld_pres_1d(KCLOUD)
        real*4 pressures_pa(NZ_L)
        real*4 clouds_3d_pres(NX_L,NY_L,NZ_L)
        real*4 wtcldcv(NX_L,NY_L,KCLOUD)

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
        real*4 dum_3d(NX_L,NY_L,NZ_L)
       
        real*4 heights_3d(NX_L,NY_L,NZ_L)

        real*4 mvd_3d(NX_L,NY_L,NZ_L)
!       real*4 lwc_res_3d(NX_L,NY_L,NZ_L)
        real*4 w_3d(NX_L,NY_L,NZ_L)

        integer icing_index_3d(NX_L,NY_L,NZ_L)

        integer cldpcp_type_3d(NX_L,NY_L,NZ_L) ! Also contains 3D precip type

!       Output array declarations (equivalences are used to share space)
        real*4 out_array_4d(NX_L,NY_L,NZ_L,3)
        real*4 out_array_3d(NX_L,NY_L,NZ_L)

        real*4 slwc(NX_L,NY_L,NZ_L),slwc_int(NX_L,NY_L)
        real*4 cice(NX_L,NY_L,NZ_L)
        real*4 pcpcnc(NX_L,NY_L,NZ_L)

!       real*4 snow_2d(NX_L,NY_L)

        character*2 c2_precip_types(0:10)

        data c2_precip_types
     1  /'  ','Rn','Sn','Zr','Sl','Ha','L ','ZL','  ','  ','  '/

        character*3 c3_pt_flag
        character*1 c1_r,c1_s

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
        real*4 dd_s(maxstns), ff_s(maxstns), ddg_s(maxstns)
        real*4 ffg_s(maxstns), vis_s(maxstns)
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
        character*100 c_values_req

        character*3 lso_ext
        data lso_ext /'lso'/

        ISTAT = INIT_TIMER()

        write(6,*)' Welcome to the LAPS gridded cloud analysis'

        call get_i_perimeter(I_PERIMETER,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'Error calling get_i_perimeter'
           stop
        endif

        IX_LOW  = 1    - I_PERIMETER
        IX_HIGH = NX_L + I_PERIMETER
        IY_LOW  = 1    - I_PERIMETER
        IY_HIGH = NY_L + I_PERIMETER

        NX_DIM_LUT = NX_L + I_PERIMETER - 1
        NY_DIM_LUT = NY_L + I_PERIMETER - 1

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

!     This should work out to be slightly larger than needed
      n_fnorm =
     1      1.6 * ( (NX_DIM_LUT*NX_DIM_LUT) + (NY_DIM_LUT*NY_DIM_LUT) )


c Determine the source of the radar data
        c_vars_req = 'radarext_3d'

        call get_static_info(c_vars_req,c_values_req,1,istatus)

        if(istatus .eq. 1)then
            radarext_3d_cloud = c_values_req(1:3)
        else
            write(6,*)' Error getting radarext_3d'
            return
        endif

        write(6,*)' radarext_3d_cloud = ',radarext_3d_cloud

c read in laps lat/lon and topo
        call get_laps_domain_95(NX_L,NY_L,lat,lon,topo
     1           ,rlaps_land_frac,grid_spacing_cen_m,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error getting LAPS domain'
            return
        endif

        write(6,*)' Actual grid spacing in domain center = '
     1                              ,grid_spacing_cen_m

        call get_laps_cycle_time(ilaps_cycle_time,istatus)
        if(istatus .eq. 1)then
            write(6,*)' ilaps_cycle_time = ',ilaps_cycle_time
        else
            write(6,*)' Error getting laps_cycle_time'
            return
        endif

        call get_cloud_parms(l_use_vis,pct_req_lvd_s8a,istatus)
        if(istatus .ne. 1)then
            write(6,*)' laps_cloud_sub: Error getting cloud parms'
            stop
        endif

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

        n_prods = 4
        iprod_start = 1
        iprod_end = n_prods

        ext = 'lc3'

c Dynamically adjust the heights of the cloud levels to the terrain
        write(6,*)' Fitting cloud heights to the terrain'
        topo_min = 1e10
        do j = 1,NY_L
        do i = 1,NX_L
            topo_min = min(topo_min,topo(i,j))
        enddo
        enddo
        write(6,*)'    OLD ht   NEW ht       Lowest terrain point = '
     1           ,topo_min

        range_orig = cld_hts(KCLOUD) - cld_hts(1)
        range_new =  cld_hts(KCLOUD) - topo_min
        range_ratio = range_new / range_orig

        do k = 1,KCLOUD
            cld_hts_orig = cld_hts(k)
            cld_hts(k) = cld_hts(KCLOUD)
     1                - (cld_hts(KCLOUD)-cld_hts_orig) * range_ratio
            write(6,21)k,cld_hts_orig,cld_hts(k)
21          format(1x,i3,2f8.1)
        enddo ! k


        write(6,*)' Initializing fields'
c initialize the fields
        if(l_perimeter)then
          do k = 1,KCLOUD
          do n = 1,max_cld_snd
              cld_snd(n,k) = r_missing_data
          enddo
          enddo
        else
          do k = 1,KCLOUD
          do j = 1,NY_L
          do i = 1,NX_L
            cldcv1(i,j,k) = r_missing_data
          enddo
          enddo
          enddo
        endif

C READ IN LAPS TEMPS
        var = 'T3'
        ext = 'lt1'
        call get_laps_3d(i4time,NX_L,NY_L,NZ_L
     1  ,ext,var,units,comment,temp_3d,istatus)

        if(istatus .ne. 1)then
            write(6,*)' Error reading 3D Temps - get_modelfg not called'
            goto999
        endif

C READ IN LAPS HEIGHTS
        var = 'HT'
        ext = 'lt1'
        call get_laps_3d(i4time,NX_L,NY_L,NZ_L
     1  ,ext,var,units,comment,heights_3d,istatus)

        if(istatus .ne. 1)then
            write(6,*)' Error reading 3D Heights - get_modelfg not calle
     1d'
            goto999
        endif

C OBTAIN MODEL FIRST GUESS CLOUD COVER FIELD
        call get_modelfg(cf_modelfg,t_modelfg,default_clear_cover
     1           ,temp_3d,model_3d,heights_3d,cld_hts
     1              ,i4time,ilaps_cycle_time
     1                  ,NX_L,NY_L,NZ_L,KCLOUD
     1                  ,istatus)

C READ IN RADAR DATA
!       Get time of radar file of the indicated appropriate extension
        call get_filespec(radarext_3d_cloud(1:3),2,c_filespec,istatus)
        call get_file_time(c_filespec,i4time,i4time_radar)

        if(abs(i4time - i4time_radar) .le. 1200)then

            call read_radar_3dref(i4time_radar,
!    1                 1200,i4time_radar,
     1                 .true.,
     1                 NX_L,NY_L,NZ_L,radarext_3d_cloud,
     1                 lat,lon,topo,.true.,.true.,
     1                 heights_3d,
     1                 radar_ref_3d,
     1                 rlat_radar,rlon_radar,rheight_radar,radar_name,     
     1                 n_ref_grids,istat_radar_2dref,istat_radar_3dref)       

        else
            write(6,*)'radar data outside time window'
            n_ref_grids = 0
            istat_radar_2dref = 0
            istat_radar_3dref = 0

        endif

C READ IN AND INSERT SAO DATA
!       Read in surface pressure
        var = 'PS'
        ext = 'lsx'
        call get_laps_2d(i4time,ext,var,units,comment
     1                  ,NX_L,NY_L,pres_sfc_pa,istatus)

        IF(istatus .ne. 1)THEN
            write(6,*)
     1  ' Error Reading Surface Pres Analyses - abort cloud analysis'
            goto999
        endif

        var = 'T'
        ext = 'lsx'
        call get_laps_2d(i4time,ext,var,units,comment
     1                  ,NX_L,NY_L,t_sfc_k,istatus)

        if(istatus .ne. 1)then
            write(6,*)' Error reading SFC Temps - abort cloud analysis'
            goto999
        endif


        write(6,*)
        write(6,*)' Call Ingest/Insert SAO routines'
        n_cld_snd = 0
        call insert_sao(i4time,cldcv1,cf_modelfg,t_modelfg
     1  ,cld_hts
     1  ,lat,lon,topo,t_sfc_k,wtcldcv
     1  ,c1_name_array,l_perimeter,ista_snd
     1  ,cvr_snd,cld_snd,wt_snd,i_snd,j_snd,n_cld_snd,max_cld_snd
     1  ,NX_L,NY_L,KCLOUD
     1  ,n_obs_pos_b,lat_s,lon_s,c_stations    ! returned for precip type comp
     1  ,wx_s,t_s,td_s,obstype                 !    "      "    "     "
     1  ,elev_s                                ! ret for comparisons
     1  ,istat_sfc,maxstns,IX_LOW,IX_HIGH,IY_LOW,IY_HIGH)

        if(istat_sfc .ne. 1)then
            write(6,*)' No SAO data inserted: Aborting cloud analysis'
            goto999
        endif


C READ IN PIREPS
        write(6,*)' Using Pireps stored in LAPS realtime system'

        call insert_pireps(i4time,cldcv1,cld_hts,wtcldcv
     1      ,cld_snd,wt_snd,i_snd,j_snd,n_cld_snd,max_cld_snd
     1      ,lat,lon,NX_L,NY_L,KCLOUD,IX_LOW,IX_HIGH,IY_LOW,IY_HIGH
     1      ,N_PIREP,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error: Bad status from insert_pireps,'
     1               ,' aborting cloud analysis'
            goto999
        endif
        I4_elapsed = ishow_timer()


C DO ANALYSIS on SAO and PIREP data
        write(6,*)
        write(6,*)' Analyzing SFC Obs and PIREP data'

        max_obs = n_cld_snd * KCLOUD

!       Set weight for using model background clouds beyond a certain effective
!       radius of influence from the sfc obs/pireps
        weight_modelfg = 0.    ! Model wt inactive, obs used to infinite radius
!       weight_modelfg = 1.    ! Model used beyond ~100km from nearest obs
!       weight_modelfg = .01   ! Model used beyond ~200km from nearest obs
!       weight_modelfg = .0001 ! Model used beyond ~400km from nearest obs

        call barnes_r5(clouds_3d,NX_L,NY_L,KCLOUD,cldcv1,wtcldcv
     1     ,cf_modelfg,l_perimeter,cld_snd,wt_snd,r_missing_data
     1     ,grid_spacing_cen_m,i_snd,j_snd,n_cld_snd,max_cld_snd
     1     ,max_obs,weight_modelfg,NX_DIM_LUT,NY_DIM_LUT
     1     ,IX_LOW,IX_HIGH,IY_LOW,IY_HIGH,n_fnorm,istatus)      
        if(istatus .ne. 1)then
            write(6,*)
     1      ' Error: Bad status from barnes_r5, aborting cloud analysis'
            goto999
        endif

!       Hold current cloud cover from SAO/PIREPS in a buffer array
        do k = 1,KCLOUD
        do j = 1,NY_L
        do i = 1,NX_L
            cldcv_sao(i,j,k) = clouds_3d(i,j,k)
        enddo
        enddo
        enddo
        I4_elapsed = ishow_timer()

        var = 'SC'
        ext = 'lm2'
        call get_laps_2d(i4time-ilaps_cycle_time,ext,var,units,comment
     1                  ,NX_L,NY_L,cvr_snow,istat_cvr_snow)

        if(istat_cvr_snow .ne. 1)then               ! Try using snow cover field
!       if(istat_cvr_snow .ne. 1 .or. .true.)then   ! Don't use snow cover
            write(6,*)' Error reading snow cover'
            do i = 1,NX_L
            do j = 1,NY_L
                cvr_snow(i,j) = r_missing_data
            enddo
            enddo
        endif

C READ IN SATELLITE DATA
!       Calculate solar altitude
        do j = 1,NY_L
        do i = 1,NX_L
!           call solalt(lat(i,j),lon(i,j),i4time,solar_alt(i,j))
            call solar_position(lat(i,j),lon(i,j),i4time,solar_alt(i,j)
     1                                     ,solar_dec,solar_ha(i,j))
        enddo
        enddo

        call get_vis(i4time,solar_alt
     1              ,cloud_frac_vis_a,albedo,ihist_alb
     1              ,NX_L,NY_L,KCLOUD,r_missing_data,istat_vis)

        call insert_sat(i4time,clouds_3d,cldcv_sao,cld_hts,lat,lon,
     1        tb8_cold_k,tb8_k,grid_spacing_cen_m,surface_sao_buffer,
     1        cloud_frac_vis_a,istat_vis,solar_alt,solar_ha,solar_dec,
     1        cloud_frac_co2_a,rlaps_land_frac,
     1        topo,heights_3d,temp_3d,t_sfc_k,t_gnd_k,sst_k,pres_sfc_pa,       
     1        dum_3d,cldtop_m_co2,cldtop_m_tb8,cldtop_m,cvr_snow,
     1        NX_L,NY_L,KCLOUD,NZ_L,istatus,r_missing_data)

        if(istatus .ne. 1)then
            write(6,*)' Error: Bad status returned from insert_sat'
            goto999
        endif

!       write(6,*)' Cloud top (Band 8 vs. Final Analysis)'
        scale = .0001
!       write(6,301)
301     format('  Cloud Top (km)             Band 8                ',
     1            20x,'              Final Analysis')
!       CALL ARRAY_PLOT(cldtop_m_tb8,cldtop_m,NX_L,NY_L,'HORZ CV',c1_name_array
!       1                                       ,KCLOUD,cld_hts,scale)

        I4_elapsed = ishow_timer()

        istat_radar_3dref_orig = istat_radar_3dref

C       THREE DIMENSIONALIZE RADAR DATA IF NECESSARY (E.G. NOWRAD)
        if(istat_radar_2dref .eq. 1 .and. istat_radar_3dref .eq. 0)then
            write(6,*)' Three dimensionalizing radar data'

!           Clear out radar echo above the highest cloud top
            k_ref_def = nint(zcoord_of_pressure(float(700*100)))

            do j = 1,NY_L
            do i = 1,NX_L
                if(radar_ref_3d(i,j,1) .gt. ref_base)then

                    k_topo = int(zcoord_of_pressure(pres_sfc_pa(i,j)))
                    k_ref = max(k_ref_def,k_topo+2)

!                   k_ref = k_ref_def

                    cloud_top_m = default_top

!                   Test for Cloud Top
                    do k = KCLOUD-1,1,-1
                        if(clouds_3d(i,j,k  ) .gt. thresh_cvr_smf .and.       
     1                     clouds_3d(i,j,k+1) .le. thresh_cvr_smf )then        
                            cloud_top_m = 0.5 * (cld_hts(k) 
     1                                         + cld_hts(k+1))
                            if(cloud_top_m .gt. heights_3d(i,j,NZ_L)
     1                                                            )then
                                cloud_top_m = heights_3d(i,j,NZ_L)
                            endif
                           goto150
                        endif

                    enddo ! k

!                   Add third dimension to radar echo
150                 if(abs(cloud_top_m) .le. 1e6)then ! Valid Cloudtop
                        k_cloud_top = ! nint(height_to_zcoord(cloud_top_m))
     1            nint(height_to_zcoord2(cloud_top_m,heights_3d
     1               ,NX_L,NY_L,NZ_L,i,j,istatus)   )
                        if(istatus .ne. 1)then
                            write(6,*)' Error: Bad status returned'
     1                          ,' from height_to_zcoord2'
     1                          ,cloud_top_m,heights_3d(i,j,NZ_L),i,j       
                            goto999
                        endif

                        k_cloud_top = max(k_ref,k_cloud_top)
                    else
                        k_cloud_top = k_ref
                    endif

                    do k = NZ_L,k_cloud_top,-1
                        radar_ref_3d(i,j,k) = -10.
                    enddo ! k

                endif ! Radar echo at this grid point
            enddo ! i
            enddo ! j

            istat_radar_3dref = 1

            I4_elapsed = ishow_timer()

        elseif(istat_radar_2dref .eq. 1 .and. istat_radar_3dref .eq. 1
     1                                                             )then       
            write(6,*)' Radar data is already three dimensional'

        elseif(istat_radar_2dref .eq. 0 .and. istat_radar_3dref .eq. 0
     1                                                             )then
            write(6,*)' Radar data is unavailable for insertion'

        endif ! Do we have 3-d or 2-d radar data?

C INSERT RADAR DATA
        if(istat_radar_3dref .eq. 1)then
            call get_max_ref(radar_ref_3d,NX_L,NY_L,NZ_L,dbz_max_2d)

            call insert_radar(i4time,clouds_3d,cld_hts
     1          ,temp_3d,t_sfc_k,grid_spacing_cen_m,NX_L,NY_L,NZ_L
     1          ,KCLOUD,cloud_base,cloud_base_buf,ref_base
     1          ,radar_ref_3d,dbz_max_2d,vis_radar_thresh_dbz
     1          ,l_unresolved
     1          ,heights_3d,istatus)

            if(istatus .ne. 1)then
                write(6,*)
     1          ' Error: Bad status returned from insert_radar'      
                goto999
            endif

        endif

        I4_elapsed = ishow_timer()

C       INSERT VISIBLE SATELLITE DATA
        if(istat_vis .eq. 1)then
            call insert_vis(i4time,clouds_3d,cld_hts
     1        ,topo,cloud_frac_vis_a,albedo,ihist_alb
     1        ,NX_L,NY_L,KCLOUD,r_missing_data
     1        ,vis_radar_thresh_cvr,vis_radar_thresh_dbz
     1        ,istat_radar_3dref,radar_ref_3d,NZ_L,ref_base
     1        ,dbz_max_2d,surface_sao_buffer,istatus)
        endif

C Clear out stuff below ground
        do k = 1,KCLOUD
        do j = 1,NY_L
        do i = 1,NX_L
            if(cld_hts(k) .lt. topo(i,j))
     1                   clouds_3d(i,j,k) = default_clear_cover
        enddo
        enddo
        enddo

C OUTPUT ARRAY in HORIZONTAL AND VERTICAL SLICES

C       HORIZONTAL SLICES

        DO K=1,KCLOUD,5
            CALL SLICE(cldcv_sao,NX_L,NY_L,KCLOUD,CVHZ1
     1                ,NX_L,NY_L,1,0,0,K,0,0)
            CALL SLICE(clouds_3d,NX_L,NY_L,KCLOUD,cvr_max
     1                ,NX_L,NY_L,1,0,0,K,0,0)
            write(6,401)k,cld_hts(k)
401         format(4x,'Lvl',i4,f8.0,' m     Before Satellite/Radar',
     1            20x,'              After Satellite/Radar')
            scale = 1.
            CALL ARRAY_PLOT(CVHZ1,cvr_max,NX_L,NY_L,'HORZ CV'
     1                     ,c1_name_array,KCLOUD,cld_hts,scale)
        ENDDO ! k

C       EW SLICES
        DO J=9,NY_L,10
!       DO J=NY_L/2+1,NY_L/2+1
            CALL SLICE(cldcv_sao,NX_L,NY_L,KCLOUD,CVEW1
     1                ,NX_L,KCLOUD,0,1,0,0,J,0)
            CALL SLICE(clouds_3d,NX_L,NY_L,KCLOUD,CVEW2
     1                ,NX_L,KCLOUD,0,1,0,0,J,0)
            write(6,501)j
501         format(5x,'  J =',i4,10x,'  Before Satellite/Radar       ',
     1          21x,'           After Satellite/Radar')

            scale = 1.
            CALL ARRAY_PLOT(CVEW1,cvew2,NX_L,KCLOUD,'VERT CV'
     1                     ,c1_name_array,KCLOUD,cld_hts,scale)
        ENDDO ! j

!       Get Max Cloud Cover
        do j = 1,NY_L
        do i = 1,NX_L
            cvr_sao_max(i,j) = 0.
            cvr_max(i,j) = 0.
            do k = 1,KCLOUD
                cvr_sao_max(i,j) = max(cvr_sao_max(i,j)
     1                                ,cldcv_sao(i,j,k))   
                cvr_max(i,j) = max(cvr_max(i,j),clouds_3d(i,j,k))
            enddo ! k
        enddo ! i
        enddo ! j

        if(istat_radar_3dref .eq. 1)then
            call compare_cloud_radar(radar_ref_3d,dbz_max_2d,cvr_max
     1         ,ref_base,cloud_frac_vis_a
     1         ,vis_radar_thresh_cvr,vis_radar_thresh_dbz,r_missing_data       
     1         ,NX_L,NY_L,NZ_L)
        endif

        write(6,601)
601     format('  Max Cloud Cover              SAO/PIREP           ',
     1            20x,'              Final Analysis')
        scale = 1.
        CALL ARRAY_PLOT(cvr_sao_max,cvr_max,NX_L,NY_L,'HORZ CV'
     1                 ,c1_name_array,KCLOUD,cld_hts,scale)

        write(6,701)
701     format('  Max Cloud Cover           VISIBLE SATELLITE     ',
     1            20x,'              Final Analysis')
        scale = 1.
        CALL ARRAY_PLOT(cloud_frac_vis_a,cvr_max,NX_L,NY_L,'HORZ CV'
     1                 ,c1_name_array,KCLOUD,cld_hts,scale)

        I4_elapsed = ishow_timer()

        do k = 1,NZ_L
            write(6,101)k,heights_3d(NX_L/2,NY_L/2,k)
101         format(1x,i3,f10.2)
        enddo ! k

        I4_elapsed = ishow_timer()

!       Write out Cloud Field
        do k=1,KCLOUD
            cld_pres_1d(k) = pressure_of_rlevel(
     1          height_to_zcoord2(cld_hts(k),heights_3d
     1               ,NX_L,NY_L,NZ_L,NX_L/2,NY_L/2,istatus)   )
            write(6,102)k,cld_pres_1d(k),istatus
102         format(1x,i3,f10.1,i2)
        enddo ! k

        ext = 'lc3'

        write(6,*)' Calling put_clouds_3d'
        call put_clouds_3d(i4time,ext,clouds_3d,cld_hts,cld_pres_1d
     1                                       ,NX_L,NY_L,KCLOUD,istatus)

        j_status(n_lc3) = ss_normal
        I4_elapsed = ishow_timer()

!       Get Cloud Bases and Tops
2000    write(6,*)' Calculating Cloud Base and Top'
        do j = 1,NY_L
        do i = 1,NX_L
            cloud_base(i,j) = default_base
            cloud_ceiling(i,j) = default_ceiling
            cloud_top(i,j)  = default_top

            if(cvr_max(i,j) .ge. thresh_cvr_base)then

!             Test for Cloud Base (MSL)
              do k = KCLOUD-1,1,-1
                if(clouds_3d(i,j,k  ) .lt. thresh_cvr_base .and.
     1             clouds_3d(i,j,k+1) .ge. thresh_cvr_base  )then
                    cloud_base(i,j) = 0.5 * (cld_hts(k) + cld_hts(k+1))
                endif
              enddo ! k

            endif ! Clouds exist in this column


            if(cvr_max(i,j) .ge. thresh_cvr_top)then

!             Test for Cloud Top (MSL)
              do k = 1,KCLOUD-1
                if(clouds_3d(i,j,k  ) .gt. thresh_cvr_top .and.
     1             clouds_3d(i,j,k+1) .le. thresh_cvr_top  )then
                    cloud_top(i,j) = 0.5 * (cld_hts(k) + cld_hts(k+1))
                endif
              enddo ! k

            endif ! Clouds exist in this column

            if(cvr_max(i,j) .ge. thresh_cvr_ceiling)then

!             Test for Cloud Ceiling (AGL)
              do k = KCLOUD-1,1,-1
                if(clouds_3d(i,j,k  ) .lt. thresh_cvr_base .and.
     1             clouds_3d(i,j,k+1) .ge. thresh_cvr_base  )then
                    cloud_ceiling(i,j) = 0.5 * 
     1                    (cld_hts(k) + cld_hts(k+1))
     1                                                - topo(i,j)
                endif
              enddo ! k

            endif ! Clouds exist in this column

        enddo ! i
        enddo ! j

!       Calculate cloud analysis implied snow cover
        call cloud_snow_cvr(cvr_max,cloud_frac_vis_a,cldtop_m_tb8
     1          ,tb8_k,NX_L,NY_L,r_missing_data,cvr_snow_cycle)

        write(6,801)
801     format('                            VISIBLE SATELLITE     ',
     1            20x,'                Snow Cover')
        scale = 1.
        CALL ARRAY_PLOT(cloud_frac_vis_a,cvr_snow_cycle,NX_L,NY_L
     1                 ,'HORZ CV',c1_name_array,KCLOUD,cld_hts,scale)       

        write(6,901)
901     format('                     lm2 (overall) Snow Cover      ',
     1            20x,'      csc  (cycle)  Snow Cover')
        scale = 1.
        CALL ARRAY_PLOT(cvr_snow,cvr_snow_cycle,NX_L,NY_L,'HORZ CV'
     1                  ,c1_name_array,KCLOUD,cld_hts,scale)

        do i = 1,NX_L
        do j = 1,NY_L
            if(cldtop_m(i,j)  .ne. r_missing_data .and.
     1       cvr_max(i,j)   .ge. 0.1                            )then
                band8_mask(i,j) = cldtop_m(i,j)
            else
                band8_mask(i,j) = r_missing_data
            endif
        enddo ! j
        enddo ! i

        write(6,*)' Cloud top (Band 8 vs. Final Analysis)'
        scale = .0001
        write(6,1001)
1001    format('  Cloud Top (km)             Band 8                ',
     1            20x,'              Final Analysis')
        CALL ARRAY_PLOT(cldtop_m_tb8,band8_mask,NX_L,NY_L,'HORZ CV'
     1                  ,c1_name_array,KCLOUD,cld_hts,scale)

!       Write out LCB file (Cloud Base, Top, and Ceiling fields)
        ext = 'lcb'
        call get_directory(ext,directory,len_dir)

        call move(cloud_base    ,out_array_3d(1,1,1),NX_L,NY_L)
        call move(cloud_top     ,out_array_3d(1,1,2),NX_L,NY_L)
        call move(cloud_ceiling ,out_array_3d(1,1,3),NX_L,NY_L)

        call put_clouds_2d(i4time,directory,ext,NX_L,NY_L
     1                                  ,out_array_3d,istatus)
        if(istatus .eq. 1)j_status(n_lcb) = ss_normal


!       This is where we will eventually split the routines, additional data
!       is necessary for more derived fields

        if(istat_radar_3dref .eq. 1)then ! Write out data (lps - radar_ref_3d)
            write(6,*)' Writing out 3D Radar Reflectivity field'

            var = 'REF'
            ext = 'lps'
            units = 'dBZ'
            write(comment,490)istat_radar_2dref,istat_radar_3dref
     1                       ,istat_radar_3dref_orig
 490        format('LAPS Radar Reflectivity',3i3)
            call put_laps_3d(i4time,ext,var,units,comment,radar_ref_3d       
     1                                                ,NX_L,NY_L,NZ_L)
            j_status(n_lps) = ss_normal

        endif ! istat_radar_3dref

        call compare_analysis_to_saos(NX_L,NY_L,cvr_sao_max
     1  ,cloud_frac_vis_a,tb8_k,t_gnd_k,t_sfc_k,cvr_max,r_missing_data
     1  ,cld_snd,ista_snd,max_cld_snd,cld_hts,KCLOUD,n_cld_snd
     1  ,c_stations,lat_s,lon_s,elev_s,maxstns)


!       Write LCV file
        do i = 1,NX_L
        do j = 1,NY_L
            cvr_water_temp(i,j) = r_missing_data
        enddo ! j
        enddo ! i

        ext = 'lcv'
        var_a(1) = 'LCV'
        var_a(2) = 'CSC'
        var_a(3) = 'CWT'
        units_a(1) = 'UNDIM'
        units_a(2) = 'UNDIM'
        units_a(3) = 'K'
        comment_a(1) = 'LAPS Cloud Cover'
        comment_a(2) = 'LAPS Cloud Analysis Implied Snow Cover'
        comment_a(3) = 'LAPS Clear Sky Water Temp'

        call move(cvr_max       ,out_array_3d(1,1,1),NX_L,NY_L)
        call move(cvr_snow_cycle,out_array_3d(1,1,2),NX_L,NY_L)
        call move(cvr_water_temp,out_array_3d(1,1,3),NX_L,NY_L)

        call put_laps_multi_2d(i4time,ext,var_a,units_a,
     1          comment_a,out_array_3d,NX_L,NY_L,3,istatus)

        if(istatus .eq. 1)j_status(n_lcv) = ss_normal

500     continue

999     continue

        write(6,*)' Notifications'
        do i = iprod_start,iprod_end
            write(6,*)' ',exts(i),' ',j_status(i),' ',i
        enddo ! i

        write(6,*)' End of Cloud Analysis Package'

        return
        end



        subroutine put_clouds_2d(i4time,DIRECTORY,ext,imax,jmax
     1                          ,field_2dcloud,istatus)

        integer  nfields
        parameter (nfields = 3)

        character*(*) DIRECTORY
        character*31 ext

        character*125 comment_2d(nfields)
        character*10 units_2d(nfields)
        character*3 var_2d(nfields)
        integer LVL,LVL_2d(nfields)
        character*4 LVL_COORD,LVL_COORD_2d(nfields)

        real*4 field_2dcloud(imax,jmax,nfields)

        write(6,11)directory,ext(1:5)
11      format(' Writing 2d Clouds ',a50,1x,a5,1x,a3)

        lvl = 0

        lvl_coord_2d(1) = 'MSL'
        lvl_coord_2d(2) = 'MSL'
        lvl_coord_2d(3) = 'AGL'

        var_2d(1) = 'LCB'
        var_2d(2) = 'LCT'
        var_2d(3) = 'CCE'

        comment_2d(1) = 'LAPS Cloud Base'
        comment_2d(2) = 'LAPS Cloud Top'
        comment_2d(3) = 'LAPS Cloud Ceiling'

        do k = 1,nfields
            LVL_2d(k) = LVL
            units_2d(k) = 'M'
        enddo

        CALL WRITE_LAPS_DATA(I4TIME,DIRECTORY,ext,imax,jmax,
     1  nfields,nfields,VAR_2D,LVL_2D,LVL_COORD_2D,UNITS_2D,
     1                     COMMENT_2D,field_2dcloud,ISTATUS)

        return
        end

        subroutine put_clouds_3d(i4time,ext,clouds_3d,cld_hts,cld_pres_1
     1d
     1                                          ,ni,nj,nk,istatus)

!       1997 Jul 31 K. Dritz  - Removed include of lapsparms.for, which was
!                               not actually needed for anything.

        character*150 DIRECTORY
        character*31 ext

        integer NZ_CLOUD_MAX
        parameter (NZ_CLOUD_MAX = 42)

        character*125 comment_3d(NZ_CLOUD_MAX),comment_2d
        character*10 units_3d(NZ_CLOUD_MAX),units_2d
        character*3 var_3d(NZ_CLOUD_MAX),var_2d
        integer LVL_3d(NZ_CLOUD_MAX)
        character*4 LVL_COORD_3d(NZ_CLOUD_MAX)

        real*4 clouds_3d(ni,nj,nk)
        real*4 cld_hts(nk)
        real*4 cld_pres_1d(nk)

        call get_directory(ext,directory,len_dir)

        var_2d = 'LC3'

        write(6,11)directory,ext(1:5),var_2d
11      format(' Writing 3d ',a50,1x,a5,1x,a3)

        do k = 1,nk
            units_3d(k)   = 'Fractional'
            lvl_3d(k) = k
            lvl_coord_3d(k) = 'MSL'

            var_3d(k) = var_2d

            write(comment_3d(k),1)cld_hts(k),cld_pres_1d(k)
1           format(2e20.8,' Height MSL, Pressure')

        enddo ! k

        CALL WRITE_LAPS_DATA(I4TIME,DIRECTORY,ext,ni,nj,
     1  nk,nk,VAR_3D,LVL_3D,LVL_COORD_3D,UNITS_3D,
     1                     COMMENT_3D,clouds_3d,ISTATUS)

        return
        end


        subroutine cloud_snow_cvr(cvr_max,cloud_frac_vis_a,cldtop_m_tb8
     1             ,tb8_k,ni,nj,r_missing_data,cvr_snow_cycle)

        real*4 cvr_max(ni,nj)          ! Input
        real*4 cloud_frac_vis_a(ni,nj) ! Input
        real*4 cldtop_m_tb8(ni,nj)     ! Input
        real*4 tb8_k(ni,nj)            ! Input
        real*4 cvr_snow_cycle(ni,nj)   ! Output

        logical l_cvr_max              ! Local

        n_csc_pts = 0
        n_no_csc_pts = 0
        n_clear_pts = 0
        n_cld_pts = 0

        n_snow = 0
        n_nosnow = 0
        n_missing = 0

        do i = 1,ni
        do j = 1,nj

            l_cvr_max = .false.

!           Make sure all neighbors are clear in addition to the grid point
            il = max(i-1,1)
            ih = min(i+1,ni)
            do ii = il,ih
                jl = max(j-1,1)
                jh = min(j+1,nj)
                do jj = jl,jh
                    if(cvr_max(ii,jj) .gt. 0.1)l_cvr_max = .true.
                enddo ! jj
            enddo ! ii

            if(                cvr_max(i,j) .le. 0.1        )then! No Cld Cover
                n_clear_pts = n_clear_pts + 1
            else
                n_cld_pts = n_cld_pts + 1
            endif

            if(        .not. l_cvr_max                          ! No Cld Cover
!           if(                cvr_max(i,j) .le. 0.1            ! No Cld Cover
     1          .and. cloud_frac_vis_a(i,j) .ne. r_missing_data
!    1          .and. cldtop_m_tb8(i,j)     .eq. r_missing_data ! No Band 8 clds
     1                                                    )then ! Definitely clear
                cvr_snow_cycle(i,j) = cloud_frac_vis_a(i,j)
                cvr_snow_cycle(i,j) = max(cvr_snow_cycle(i,j),0.)
                cvr_snow_cycle(i,j) = min(cvr_snow_cycle(i,j),1.)

                n_csc_pts = n_csc_pts + 1

            else                               ! Cloud Cover or missing albedo
                cvr_snow_cycle(i,j) = r_missing_data
                n_no_csc_pts = n_no_csc_pts + 1

            endif

!           If ground is warm then no snow is present
            if(tb8_k(i,j) .ne. r_missing_data
     1       .and.        tb8_k(i,j) .gt. 281.15
     1       .and.         cvr_max(i,j) .le. 0.1
     1                                                       )then
                cvr_snow_cycle(i,j) = 0.
            endif

!           Count various categories of CSC
            if(cvr_snow_cycle(i,j) .ne. r_missing_data)then
                if(cvr_snow_cycle(i,j) .gt. 0.1)then
                    n_snow = n_snow + 1
                else
                    n_nosnow = n_nosnow + 1
                endif
            else
                n_missing = n_missing + 1
            endif

        enddo ! j
        enddo ! i

        write(6,*)' # csc/nocsc/clr/cld pts = ',n_csc_pts,n_no_csc_pts
     1                                       ,n_clear_pts,n_cld_pts
        write(6,1)n_snow,n_nosnow,n_missing
 1      format('  # snow/nosnow/missing ',3i7)

        return
        end

        subroutine get_modelfg(cf_modelfg,t_modelfg,default_clear_cover
     1                    ,temp_3d,model_3d,heights_3d,cld_hts
     1              ,i4time_needed,ilaps_cycle_time
     1              ,ni,nj,klaps,KCLOUD
     1  ,istatus)

!       Obtain model first guess cloud cover and temperature fields
!       This should probably be free of very small scale horizontal structures
!       for best results when combining with the satellite data

!       This routine is currently under development, it just assigns the
!       default value to the model first guess for now.

!                Steve Albers   Use model first guess info to guess cloud cover
!       1994     Steve Albers   Read in SH from model background (instead of
!                                                                         RH)
!       1995 Dec Steve Albers   QC check to prevent cf_modelfg > 1.0

        real*4 heights_3d(ni,nj,klaps)   ! Input
        real*4 temp_3d(ni,nj,klaps)      ! Input
        real*4 model_3d(ni,nj,klaps)     ! Local
        real*4 cf_modelfg(ni,nj,KCLOUD)  ! Output
        real*4 t_modelfg(ni,nj,KCLOUD)   ! Output
        real*4 cld_hts(KCLOUD)           ! Input

        real*4 make_rh

        character*31 ext
        character*3 var_2d
        character*125 comment_2d
        character*10 units_2d

        logical l_fill

!       Initialize model first guess cover field with default value
        do k = 1,KCLOUD
        do j = 1,nj
        do i = 1,ni
            cf_modelfg(i,j,k) = default_clear_cover
            t_modelfg(i,j,k) = 0.
        enddo
        enddo
        enddo

        i_hum_high = 0
        i_hum_low = 0
        i_hum_ok = 0

!       Read model sh
        istat_sh = 0

        if(.true.)then
            write(6,*)
            write(6,*)' Getting first guess cloud cover'
            write(6,*)' Getting MODEL SH background'

!           Get MAPS Data
            var_2d = 'SH'
            l_fill = .true.

            call get_modelfg_3d(i4time_needed,var_2d,ni,nj,klaps
     1                                  ,model_3d,istat_sh)

            if(istat_sh .ne. 1)then
                write(6,*)' No MAPS: No first guess available for SH '
                return
            endif

        endif ! Good status for RAMS data


!       Remap to cloud height grid and convert to cloud cover
        t_ref = 0.

        do k = 1,KCLOUD
        do j = 1,nj
        do i = 1,ni

!           Find the model pressure at this location in the cloud height grid
            if(i-1 .eq. (i-1)/10*10)then ! Update every 10th grid point
                z_laps = height_to_zcoord2(cld_hts(k),heights_3d
     1                  ,ni,nj,klaps,i,j,istatus)
                if(istatus .ne. 1)then
!                   Determine if cloud height grid is above pressure grid
                    if(cld_hts(k) .gt. heights_3d(i,j,klaps))then
                        i_grid_high = i_grid_high + 1
                        cf_modelfg(i,j,k) = default_clear_cover
                        go to 1000
                    else
                        write(6,*)' Error: Bad status from '
     1                           ,'height_to_zcoord2'
                        return
                    endif
                endif

                z_laps = max(1.,min(z_laps,float(klaps)-.001))
                iz_laps = int(z_laps)
                frac = z_laps - iz_laps

                p_modelfg =  pressure_of_level(iz_laps) * (1. - frac)
     1                  +  pressure_of_level(iz_laps+1)  * frac

                p_modelfg_mb = p_modelfg * .01

            endif

!           Find the model temp at this location in the cloud height grid
            t_modelfg(i,j,k) =  temp_3d(i,j,iz_laps)    * (1. - frac)
     1                       +  temp_3d(i,j,iz_laps+1)  * frac

            t_modelfg_c = t_modelfg(i,j,k) - 273.15

!           Find the model sh at this location in the cloud height grid
            if(istat_sh .eq. 1)then
                q_modelfg =  model_3d(i,j,iz_laps)    * (1. - frac)
     1                    +  model_3d(i,j,iz_laps+1)  * frac

                q_modelfg_gkg = q_modelfg * 1000.

                rh_modelfg = make_rh(p_modelfg_mb              ! fractional rh
     1                     ,t_modelfg_c,q_modelfg_gkg,t_ref)

            else
                write(6,*)' Code error in get_modelfg, STOP'
                stop
!               istatus = 0
!               return

            endif

!           QC the rh
            rh_qc = rh_modelfg                          ! fractional rh

            if(rh_qc .gt. 1.0)then
                rh_qc = 1.0
                i_hum_high = i_hum_high + 1
            elseif(rh_qc .lt. 0.0)then
                rh_qc = 0.0
                i_hum_low = i_hum_low + 1
            else
                i_hum_ok = i_hum_ok + 1
            endif

            if(cld_hts(k) .gt. 11000.)rh_qc = .01       ! set upper lvls to dry
                                                        ! counters model (ruc) 
                                                        ! moist bias

            cf_modelfg(i,j,k) = rh_to_cldcv(rh_qc)      ! fractional_rh

 1000   enddo ! i
        enddo ! j
        enddo ! k (cloud height array level)

        write(6,*)' # RH values QCed  high/low/ok'
     1            ,i_hum_high,i_hum_low,i_hum_ok
        write(6,*)' # cloud height grids above pressure grid'
     1            ,i_grid_high

        return
        end

        function rh_to_cldcv(rh)

        rh_to_cldcv = rh                                ! fractional rh

        return
        end

        subroutine compare_analysis_to_saos(ni,nj,cvr_sao_max
     1  ,cloud_frac_vis_a,tb8_k,t_gnd_k,t_sfc_k,cvr_max,r_missing_data
     1  ,cld_snd,ista_snd,max_cld_snd,cld_hts,KCLOUD,n_cld_snd
     1  ,c_stations,lat_s,lon_s,elev_s,maxstns)

        real*4 cloud_frac_vis_a(ni,nj),tb8_k(ni,nj),t_gnd_k(ni,nj)
     1        ,t_sfc_k(ni,nj),cvr_max(ni,nj),cvr_sao_max(ni,nj)

        real*4 cld_snd(max_cld_snd,KCLOUD)
        integer ista_snd(max_cld_snd)
        real*4 cld_hts(KCLOUD)

        character c_stations(maxstns)*(*)
        real*4 lat_s(maxstns), lon_s(maxstns), elev_s(maxstns)

        character*3 c3_discrep
        character*1 c1_c

        do j = 1,nj
        do i = 1,ni
            cvr_max(i,j) = min(cvr_max(i,j),1.00)
        enddo ! i
        enddo ! j

        if(.true.)then
            write(6,*)
     1      ' Comparing cloud/sat/sfc data at SAO/pirep locations'
            iwrite = 0

            n_ovc = 0
            n_tovc = 0
            n_sct = 0
            n_bkn = 0
            n_tsct = 0
            n_tbkn = 0
            ovc_sum = 0.
            tovc_sum = 0.
            sct_sum = 0.
            bkn_sum = 0.
            tsct_sum = 0.
            tbkn_sum = 0.


            do isnd = 1,n_cld_snd
              ista = ista_snd(isnd)
              if(ista .ne. 0)then
                call latlon_to_rlapsgrid(lat_s(ista),lon_s(ista),lat,lon
     1                          ,ni,nj,ri,rj,istatus)

                i_i = nint(ri)
                i_j = nint(rj)

                if(i_i .ge. 3 .and. i_i .le. ni-2 .and.
     1             i_j .ge. 3 .and. i_j .le. nj-2            )then

                    if(iwrite .eq. iwrite/20*20)then
                        write(6,*)
                        write(6,*)'Sta  VIS frac tb8_k  '
     1                  //'t_gnd_k t_sfc_k  cldsnd cv_sa_mx cvr_mx '
     1                  //'snd-ht        9pt    25pt'
                    endif

!                   Calculate 9pt cover
                    cvr_9pt = 0.
                    do ii = -1,1
                    do jj = -1,1
                        cvr_9pt = cvr_9pt + cvr_max(i_i+ii,i_j+jj)
                    enddo ! jj
                    enddo ! ii
                    cvr_9pt = cvr_9pt / 9.

!                   Calculate 25pt cover
                    cvr_25pt = 0.
                    do ii = -2,2
                    do jj = -2,2
                        cvr_25pt = cvr_25pt + cvr_max(i_i+ii,i_j+jj)
                    enddo ! jj
                    enddo ! ii
                    cvr_25pt = cvr_25pt / 25.

                    iwrite = iwrite + 1

                    cld_snd_max = -.999
                    height_max = 0.
                    do k = 1,KCLOUD
                        if(cld_snd(isnd,k) .ne. r_missing_data)then
                            cld_snd_max = max(cld_snd_max
     1                                       ,cld_snd(isnd,k))
                            height_max = cld_hts(k)
                        endif
                    enddo ! k

!                   cld_snd_max = cvr_snd(isnd)

!                   Flag significant discrepancies between SAOs and analysis
                    ht_12000 = elev_s(ista) + 12000./3.281 + 1.
                    c3_discrep = '   '

                    if(cvr_max(i_i,i_j) - cld_snd_max .le. -0.25)then
                        if(cvr_max(i_i,i_j) .le. 0.1)then
                            c3_discrep = ' **'
                        else
                            c3_discrep = ' . '
                        endif
                    endif

                    if(cvr_max(i_i,i_j) - cld_snd_max .ge. +0.25)then
                        if(height_max .gt. ht_12000)then
                            if(cld_snd_max .le. 0.1)then
                               c3_discrep = ' **'
                            else
                               c3_discrep = ' . '
                            endif
                        endif
                    endif

                    if(cvr_max(i_i,i_j) - cld_snd_max .ge. +0.50)then
                        if(height_max .gt. ht_12000)c3_discrep = ' **'
                    endif

                    if(cvr_max(i_i,i_j) - cld_snd_max .le. -0.50)then
                        c3_discrep = ' **'
                    endif

                    if(cvr_sao_max(i_i,i_j) .lt. cld_snd_max - 0.1)then       
                        c3_discrep = ' SS'
                    endif

                    icat1 = 1
                    if(cvr_25pt .ge. 0.10)icat1 = 2
                    if(cvr_25pt .ge. 0.50)icat1 = 3
                    if(cvr_25pt .ge. 0.90)icat1 = 4

                    icat2 = 1
                    if(cld_snd_max .ge. 0.10)icat2 = 2
                    if(cld_snd_max .ge. 0.50)icat2 = 3
                    if(cld_snd_max .ge. 0.90)icat2 = 4

                    if(icat1 .ne. icat2
     1            .and. abs(cld_snd_max - cvr_25pt) .gt. .10
     1                              .AND.
     1           (height_max .gt. ht_12000 .or. cld_snd_max .eq. 1.00)
     1                                                             )then       
                        c1_c = 'C'
                    else
                        c1_c = ' '
                    endif

                    if(abs(cld_snd_max - 1.00) .lt. .01)then
                        n_ovc = n_ovc + 1
                        ovc_sum = ovc_sum + cvr_25pt
                    endif

                    if(height_max .gt. ht_12000)then
                        if(abs(cld_snd_max - .60) .lt. .01)then
                            n_tovc = n_tovc + 1
                            tovc_sum = tovc_sum + cvr_25pt
                        endif
                        if(abs(cld_snd_max - .25) .lt. .01)then
                            n_sct = n_sct + 1
                            sct_sum = sct_sum + cvr_25pt
                        endif
                        if(abs(cld_snd_max - .70) .lt. .01)then
                            n_bkn = n_bkn + 1
                            bkn_sum = bkn_sum + cvr_25pt
                        endif
                        if(abs(cld_snd_max - .15) .lt. .01)then
                            n_tsct = n_tsct + 1
                            tsct_sum = tsct_sum + cvr_25pt
                        endif
                        if(abs(cld_snd_max - .40) .lt. .01)then
                            n_tbkn = n_tbkn + 1
                            tbkn_sum = tbkn_sum + cvr_25pt
                        endif
                    endif ! > 12000 AGL type station

                    write(6,1111,err=1112)c_stations(ista)(1:3)
     1                           ,cloud_frac_vis_a(i_i,i_j)
     1                           ,tb8_k(i_i,i_j)
     1                           ,t_gnd_k(i_i,i_j)
     1                           ,t_sfc_k(i_i,i_j)
     1                           ,cld_snd_max
     1                           ,cvr_sao_max(i_i,i_j)
     1                           ,cvr_max(i_i,i_j)
     1                           ,height_max
     1                           ,c3_discrep
     1                           ,cvr_9pt
     1                           ,cvr_25pt
     1                           ,c1_c
1111                format(1x,a3,f8.2,3f8.1,3f8.2,f8.0,a3,f8.2,f7.2
     1                    ,1x,a1)

1112            endif ! ob is in domain
              endif ! ista .ne. 0 (valid value)
            enddo ! isnd

            if(n_sct .gt. 0)write(6,11)n_sct,sct_sum/float(n_sct)
11          format(' Mean analysis value for  SCT (.25) sao  = ',i3,f7.2
     1)

            if(n_bkn .gt. 0)write(6,12)n_bkn,bkn_sum/float(n_bkn)
12          format(' Mean analysis value for  BKN (.70) sao  = ',i3,f7.2
     1)

            if(n_tsct .gt. 0)write(6,13)n_tsct,tsct_sum/float(n_tsct)
13          format(' Mean analysis value for -SCT (.20) sao  = ',i3,f7.2
     1)

            if(n_tbkn .gt. 0)write(6,14)n_tbkn,tbkn_sum/float(n_tbkn)
14          format(' Mean analysis value for -BKN (.40) sao  = ',i3,f7.2
     1)

            if(n_ovc .gt. 0)write(6,15)n_ovc,ovc_sum/float(n_ovc)
15          format(' Mean analysis value for  OVC (1.00) sao = ',i3,f7.2
     1)

            if(n_tovc .gt. 0)write(6,16)n_tovc,tovc_sum/float(n_tovc)
16          format(' Mean analysis value for -OVC  (.60) sao = ',i3,f7.2
     1)

            write(6,*)

        endif ! do SAO comparison to analysis

        return
        end


        subroutine compare_cloud_radar(radar_ref_3d,dbz_max_2d,cvr_max
     1          ,ref_base,cloud_frac_vis_a
     1          ,vis_radar_thresh_cvr,vis_radar_thresh_dbz
     1          ,r_missing_data,ni,nj,nk)

        real*4 radar_ref_3d(ni,nj,nk)
        real*4 dbz_max_2d(ni,nj)
        real*4 cvr_max(ni,nj)
        real*4 cloud_frac_vis_a(ni,nj)

!       This routine compares the cloud and radar fields and flags
!       remaining differences that weren't caught in earlier processing
!       This routine does not alter any arrays or pass anything back,
!       it is diagnostic only.

        call get_max_ref(radar_ref_3d,ni,nj,nk,dbz_max_2d)

        write(6,*)'Comparing clouds and radar (_RDR)'
        write(6,*)'vis_radar_thresh_cvr = ',vis_radar_thresh_cvr
        write(6,*)'vis_radar_thresh_dbz = ',vis_radar_thresh_dbz

        do i = 1,ni
        do j = 1,nj
            if(       cvr_max(i,j)    .lt. 0.2
     1          .and. dbz_max_2d(i,j) .gt. ref_base
     1          .and. dbz_max_2d(i,j) .ne. r_missing_data
     1                                                            )then

!             We have a discrepancy between the VIS and radar

              if(cvr_max(i,j) .eq. cloud_frac_vis_a(i,j))then ! CVR = VIS

                if(dbz_max_2d(i,j) .lt. vis_radar_thresh_dbz)then
                  write(6,1)i,j,cvr_max(i,j),dbz_max_2d(i,j)
     1                     ,cloud_frac_vis_a(i,j)
1                 format(' VIS_RDR: cvr/dbz/vis <',2i4,f8.2,f8.1,f8.2)

                else
                  write(6,11)i,j,cvr_max(i,j),dbz_max_2d(i,j)
     1                      ,cloud_frac_vis_a(i,j)
11                format(' VIS_RDR: cvr/dbz/vis >',2i4,f8.2,f8.1,f8.2)

                endif


              elseif(cvr_max(i,j) .lt. cloud_frac_vis_a(i,j))then

!                 Don't blame the VIS            CVR < VIS

                  if(dbz_max_2d(i,j) .lt. vis_radar_thresh_dbz)then

!                     We can potentially block out the radar
                      write(6,2)i,j,cvr_max(i,j),dbz_max_2d(i,j)
     1                       ,cloud_frac_vis_a(i,j)
2                     format(' CLD_RDR: cvr/dbz/vis <',2i4,f8.2
     1                                                 ,f8.1,f8.2)

                  else ! Radar is too strong to block out

                      write(6,3)i,j,cvr_max(i,j),dbz_max_2d(i,j)
     1                       ,cloud_frac_vis_a(i,j)
3                     format(' CLD_RDR: cvr/dbz/vis >',2i4,f8.2
     1                                                 ,f8.1,f8.2)

                  endif

              elseif(cvr_max(i,j) .gt. cloud_frac_vis_a(i,j))then

!                 Don't know if VIS lowered cloud cover        CVR > VIS
!                 At least if difference is less than VIS "cushion"

                  if(dbz_max_2d(i,j) .lt. vis_radar_thresh_dbz)then

!                     We can potentially block out the radar
                      write(6,4)i,j,cvr_max(i,j),dbz_max_2d(i,j)
     1                       ,cloud_frac_vis_a(i,j)
4                     format(' ???_RDR: cvr/dbz/vis <',2i4,f8.2
     1                                                 ,f8.1,f8.2)

                  else ! Radar is too strong to block out

                      write(6,5)i,j,cvr_max(i,j),dbz_max_2d(i,j)
     1                       ,cloud_frac_vis_a(i,j)
5                     format(' ???_RDR: cvr/dbz/vis >',2i4,f8.2
     1                                                 ,f8.1,f8.2)

                  endif

              endif

            endif
        enddo ! j
        enddo ! i

        return
        end

