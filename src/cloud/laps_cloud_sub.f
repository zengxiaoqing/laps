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

        subroutine laps_cloud_sub(i4time,
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

!       Obtain cloud parameters
        use mem_namelist, ONLY: l_use_vis,l_use_vis_add                
     1                         ,l_use_vis_partial                      
     1                         ,l_use_39,latency_co2                   
     1                         ,pct_req_lvd_s8a                        
     1                         ,cld_weight_modelfg
     1                         ,i4_sat_window,i4_sat_window_offset     
     1                         ,l_use_metars,l_use_radar,iwrite_output
     1                         ,i_varadj

        integer       ss_normal,sys_bad_prod,sys_no_data,
     1                sys_abort_prod

        parameter (ss_normal      =1, ! success
     1             sys_bad_prod   =2, ! inappropriate data, insufficient data
     1             sys_no_data    =3, ! no data
     1             sys_abort_prod =4) ! failed to make a prod

        include 'cloud.inc'

!       Prevents clearing out using satellite (hence letting SAOs dominate)
!       below this altitude (M AGL)
        real surface_sao_buffer
        parameter (surface_sao_buffer = 800.)

        real thresh_cvr_ceiling,default_top,default_base
     1                           ,default_clear_cover,default_ceiling       

        parameter       (thresh_cvr_ceiling = 0.65) ! Used to "binaryize" cloud cover

        parameter       (default_clear_cover = .001)

        real thresh_cvr_base,thresh_cvr_top
        parameter (thresh_cvr_base = 0.1)
        parameter (thresh_cvr_top  = 0.1)

        real thresh_thin_lwc_ice     ! Threshold cover for thin cloud LWC/ICE
        parameter (thresh_thin_lwc_ice = 0.1)

        real vis_radar_thresh_cvr,vis_radar_thresh_dbz
        parameter (vis_radar_thresh_cvr = 0.2)  ! 0.2, 0.0
        parameter (vis_radar_thresh_dbz = 30.)  ! 10., 5. , -99.

        real lat(NX_L,NY_L),lon(NX_L,NY_L)
        real topo(NX_L,NY_L)
        real rlaps_land_frac(NX_L,NY_L)
        real solar_alt(NX_L,NY_L)
        real solar_az(NX_L,NY_L)
        real solar_ha(NX_L,NY_L)

        logical l_packed_output 
        logical l_use_co2_mode1, l_use_co2_mode2
        logical l_evap_radar, l_get_cloudtype
        logical l_trust_narrowband ! should narrowband be trusted (i.e. not 
                                   ! QC'd out) in the absence of other data?

        logical lstat_co2_a(NX_L,NY_L)

        data l_packed_output /.false./
        data l_evap_radar /.false./
        data l_trust_narrowband /.false./ ! overriden to .true. for NOWRAD

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

        equivalence (cld_hts,cld_hts_new)

        real, allocatable, dimension(:,:,:) :: cldcv1
        real, allocatable, dimension(:,:,:) :: cf_modelfg
        real, allocatable, dimension(:,:,:) :: t_modelfg
        real, allocatable, dimension(:,:,:) :: sh_modelfg
        real, allocatable, dimension(:,:,:) :: ref_modelfg
!       REAL cldcv1(NX_L,NY_L,KCLOUD)
!       REAL cf_modelfg(NX_L,NY_L,KCLOUD)
!       REAL t_modelfg(NX_L,NY_L,KCLOUD)

        real clouds_3d(NX_L,NY_L,KCLOUD)

        integer ista_snd(max_cld_snd)

        real, allocatable, dimension(:,:) :: cld_snd
        real, allocatable, dimension(:,:) :: wt_snd
!       real cld_snd(max_cld_snd,KCLOUD)
!       real wt_snd(max_cld_snd,KCLOUD)

        real cvr_snd(max_cld_snd)
        integer i_snd(max_cld_snd)
        integer j_snd(max_cld_snd)

        integer ihist_alb(-10:20)

        real cloud_top(NX_L,NY_L)
        real cloud_base(NX_L,NY_L)
        real cloud_ceiling(NX_L,NY_L)

        real cldtop_m(NX_L,NY_L)
        real cldtop_co2_m(NX_L,NY_L)
        real cldtop_tb8_m(NX_L,NY_L)
        real cldtop_co2_pa_a(NX_L,NY_L)
        real ht_sao_top(NX_L,NY_L)

        real cldcv_sao(NX_L,NY_L,KCLOUD)
        real cld_pres_1d(KCLOUD)
        real pressures_pa(NZ_L)
        real pres_3d(NX_L,NY_L,NZ_L)
        real wtcldcv(NX_L,NY_L,KCLOUD)

        real CVHZ(NX_L,NY_L)
        real CVHZ1(NX_L,NY_L),CVEW1(NX_L,KCLOUD)
        real cvr_max(NX_L,NY_L),CVEW2(NX_L,KCLOUD)
        real cvr_sao_max(NX_L,NY_L)
        real cvr_snow_cycle(NX_L,NY_L)
        real cvr_water_temp(NX_L,NY_L)
        real cvr_snow(NX_L,NY_L)
        real plot_mask(NX_L,NY_L)
        real plot_maskr(NX_L,NY_L)

        character*4 radar_name
        character*31 radarext_3d_cloud
        real radar_ref_3d(NX_L,NY_L,NZ_L)
        real closest_radar(NX_L,NY_L)
        integer istat_radar_2dref_a(NX_L,NY_L)
        integer istat_radar_3dref_a(NX_L,NY_L)
        logical lstat_radar_3dref_orig_a(NX_L,NY_L)
       
        real heights_3d(NX_L,NY_L,NZ_L)

!       Output array declarations
        real out_array_3d(NX_L,NY_L,8)

!       real snow_2d(NX_L,NY_L)

        character*2 c2_precip_types(0:10)

        data c2_precip_types
     1  /'  ','Rn','Sn','Zr','Sl','Ha','L ','ZL','  ','  ','  '/

        integer i2_pcp_type_2d(NX_L,NY_L)
        real r_pcp_type_2d(NX_L,NY_L)

        real dum1_array(NX_L,NY_L)
        real dum2_array(NX_L,NY_L)
        real dum3_array(NX_L,NY_L)
        real dum4_array(NX_L,NY_L)

!       Arguments for calling get_cloud_deriv
        logical l_mask_pcptype(NX_L,NY_L)
        integer ibase_array(NX_L,NY_L)
        integer itop_array(NX_L,NY_L)
        logical   l_flag_cloud_type,l_flag_mvd,l_flag_icing_index
     1           ,l_flag_bogus_w,l_bogus_radar_w,l_deep_vv
     1           ,l_flag_pcp_type

        logical l_unresolved(NX_L,NY_L)

        character*1 c1_name_array(NX_L,NY_L,KCLOUD)

        integer MAX_FIELDS
        parameter (MAX_FIELDS = 10)

        character*255 c_filespec
        character var*3,comment*125,directory*150,ext*31,units*10
        character*125 comment_tb8,comment_t39,comment_sst,comment_alb       
        character*3 exts(20)
        character*3 var_a(MAX_FIELDS)
        character*125 comment_a(MAX_FIELDS)
        character*10  units_a(MAX_FIELDS)

!       Arrays used to read in satellite data
        real tb8_k(NX_L,NY_L)
        real t39_k(NX_L,NY_L)
        real tb8_cold_k(NX_L,NY_L)
        real albedo(NX_L,NY_L)
        real static_albedo(NX_L,NY_L)              ! Static albedo database
        real cloud_frac_vis_a(NX_L,NY_L)
        real cloud_frac_co2_a(NX_L,NY_L)
        real subpoint_lat_clo_vis(NX_L,NY_L)
        real subpoint_lon_clo_vis(NX_L,NY_L)
        real di_dh(NX_L,NY_L)                      
        real dj_dh(NX_L,NY_L)                      

        integer istat_39_a(NX_L,NY_L)
        integer istat_39_add_a(NX_L,NY_L)
        integer istat_vis_potl_a(NX_L,NY_L)
        integer istat_vis_added_a(NX_L,NY_L)

        real temp_3d(NX_L,NY_L,NZ_L)

        real t_sfc_k(NX_L,NY_L)
        real t_gnd_k(NX_L,NY_L)
        real sst_k(NX_L,NY_L)
        real td_sfc_k(NX_L,NY_L)
        real pres_sfc_pa(NX_L,NY_L)

!       Declarations for LSO file stuff
        real lat_s(maxstns), lon_s(maxstns), elev_s(maxstns)
        real cover_s(maxstns)
        real t_s(maxstns), td_s(maxstns), pr_s(maxstns), sr_s(maxstns)
        real dd_s(maxstns), ff_s(maxstns), ddg_s(maxstns)
        real ffg_s(maxstns), vis_s(maxstns)
        real pstn_s(maxstns),pmsl_s(maxstns),alt_s(maxstns)
        real store_hgt(maxstns,5),ceil(maxstns),lowcld(maxstns)
        real cover_a(maxstns),rad_s(maxstns),solar_ea(maxstns)
        integer obstime(maxstns),kloud(maxstns),idp3(maxstns)
        character store_emv(maxstns,5)*1,store_amt(maxstns,5)*4
        character wx_s(maxstns)*8
        character atype(maxstns)*6
        character reptype(maxstns)*6

        integer STATION_NAME_LEN
        parameter (STATION_NAME_LEN = 3)                   
        character c_stations(maxstns)*(STATION_NAME_LEN)    

        real ri_s(maxstns), rj_s(maxstns)

!       Product # notification declarations
        integer j_status(20),iprod_number(20)

!       Stuff for 2d fields
        real ref_mt_modelfg(NX_L,NY_L)
        real dbz_low_2d(NX_L,NY_L)
        real dbz_max_2d(NX_L,NY_L)
        real rqc_2d(NX_L,NY_L)
        real swi_2d(NX_L,NY_L)

!       SFC precip and cloud type (LCT file)
        real r_pcp_type_thresh_2d(NX_L,NY_L)
        real r_cld_type_2d(NX_L,NY_L)

        character*40 c_vars_req
        character*180 c_values_req

        character*3 lso_ext        
        data lso_ext /'lso'/

        allocate( cldcv1(NX_L,NY_L,KCLOUD), STAT=istat_alloc )
        if(istat_alloc .ne. 0)then
            write(6,*)' ERROR: Could not allocate cldcv1'
        endif

        allocate( cf_modelfg(NX_L,NY_L,KCLOUD), STAT=istat_alloc )
        if(istat_alloc .ne. 0)then
            write(6,*)' ERROR: Could not allocate cf_modelfg'
        endif

        allocate( t_modelfg(NX_L,NY_L,KCLOUD), STAT=istat_alloc )
        if(istat_alloc .ne. 0)then
            write(6,*)' ERROR: Could not allocate t_modelfg'
        endif

        allocate( sh_modelfg(NX_L,NY_L,NZ_L), STAT=istat_alloc )
        if(istat_alloc .ne. 0)then
            write(6,*)' ERROR: Could not allocate sh_modelfg'
        endif

        allocate( ref_modelfg(NX_L,NY_L,NZ_L), STAT=istat_alloc )
        if(istat_alloc .ne. 0)then
            write(6,*)' ERROR: Could not allocate ref_modelfg'
        endif

        allocate( cld_snd(max_cld_snd,KCLOUD), STAT=istat_alloc )
        if(istat_alloc .ne. 0)then
            write(6,*)' ERROR: Could not allocate cld_snd'
        endif

        allocate( wt_snd(max_cld_snd,KCLOUD), STAT=istat_alloc )
        if(istat_alloc .ne. 0)then
            write(6,*)' ERROR: Could not allocate wt_snd'
        endif

        ISTAT = INIT_TIMER()

        write(6,*)' Welcome to the LAPS gridded cloud analysis'

        lstat_radar_3dref_orig_a = .false.

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
              c1_name_array(i,j,:) = ' '
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


!       Fill namelist_parms data structure
        namelist_parms%l_use_metars = l_use_metars
        namelist_parms%l_use_radar = l_use_radar

        if(latency_co2 .ge. 0)then
            l_use_co2_mode1 = .true.
        else
            l_use_co2_mode1 = .false.
        endif

        l_use_co2_mode2 = .false.

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

        topo_min_buff = topo_min - .005 ! Put in a slight buffer so that
                                        ! cloud height grid extends below topo

        write(6,*)'      OLD ht     NEW ht      Lowest terrain point = '       
     1           ,topo_min

        range_orig = cld_hts(KCLOUD) - cld_hts(1)
        range_new =  cld_hts(KCLOUD) - topo_min_buff
        range_ratio = range_new / range_orig

        do k = 1,KCLOUD
            cld_hts_orig = cld_hts(k)
            cld_hts(k) = cld_hts(KCLOUD)
     1                - (cld_hts(KCLOUD)-cld_hts_orig) * range_ratio
            write(6,21)k,cld_hts_orig,cld_hts(k)
21          format(1x,i3,2f10.3)
        enddo ! k

        if(cld_hts(1) .ge. topo_min)then
            write(6,*)' Error: topo extends at or below edge of '
     1               ,'cloud height grid',topo_min,cld_hts(1)
            return
        endif

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

C OBTAIN MODEL FIRST GUESS CLOUD COVER FIELD (along with reflectivity, water vapor)
        call get_modelfg(cf_modelfg,t_modelfg,sh_modelfg,ref_modelfg    ! O
     1           ,default_clear_cover,r_missing_data                    ! I
     1           ,temp_3d,heights_3d,cld_hts                            ! I
     1              ,i4time,ilaps_cycle_time                            ! I
     1                  ,NX_L,NY_L,NZ_L,KCLOUD                          ! I
     1                  ,istatus)                                       ! O

C READ IN RADAR DATA
!       Get time of radar file of the indicated appropriate extension
        call get_filespec(radarext_3d_cloud(1:3),2,c_filespec,istatus)
        call get_file_time(c_filespec,i4time,i4time_radar)

        if(namelist_parms%l_use_radar)then
            i4_tol = 1200
        else
            write(6,*)' Withholding radar data from cloud analysis'
            i4_tol = -1
        endif

        call read_multiradar_3dref(i4time,                               ! I
     1                 i4_tol,i4_ret,                                    ! I
     1                 .true.,ref_base,                                  ! I
     1                 NX_L,NY_L,NZ_L,radarext_3d_cloud,                 ! I
     1                 lat,lon,topo,.true.,.true.,                       ! I
     1                 heights_3d,                                       ! I
     1                 radar_ref_3d,                                     ! O
     1                 rlat_radar,rlon_radar,rheight_radar,radar_name,   ! O
     1                 iqc_2dref,closest_radar,                          ! O
     1                 n_ref_grids,n_2dref,n_3dref,istat_radar_2dref_a,  ! O  
     1                 istat_radar_3dref_a)                              ! O

        if(iqc_2dref .eq. 1)then
            write(6,*)' Good quality 2d radar (e.g. NOWRAD)'
            l_trust_narrowband = .true.
        else
            write(6,*)' Lesser quality 2d radar (e.g. narrowband)'      
        endif

        rqc_2d = 0.
        where (istat_radar_2dref_a .eq. 1)rqc_2d = 2.
        where (istat_radar_3dref_a .eq. 1)rqc_2d = 3.

!       else
!           write(6,*)'radar data outside time window'
!           n_ref_grids = 0
!           call constant_i(istat_radar_2dref_a,0,NX_L,NY_L)       
!           call constant_i(istat_radar_3dref_a,0,NX_L,NY_L)

!       endif

C BLEND IN FIRST GUESS RADAR
        n_fg_radar = 0
        n_fg_echoes = 0

        call get_maxtops(ref_modelfg,heights_3d
     1                  ,NX_L,NY_L,NZ_L,ref_mt_modelfg)

        write(6,*)' Modelfg maxtops range is: '
     1      ,minval(ref_mt_modelfg),maxval(ref_mt_modelfg)
        write(6,*)' Closest radar range is: '
     1      ,minval(closest_radar),maxval(closest_radar)

        range_thresh_lo = 200000.
        range_thresh_hi = 300000.

        radius_earth_8_thirds = 6371.e3 * 2.6666666
        aterm = 1. / radius_earth_8_thirds

        elev_ang_thr = 0.5 ! angle of radar horizon 
        bterm = tand(elev_ang_thr)

        do i = 1,NX_L
        do j = 1,NY_L

!           Determine dynamic cutoff distance between radar data and first guess
            if(ref_mt_modelfg(i,j) .gt. 0.)then ! find distance from radar given
                cterm = ref_mt_modelfg(i,j)     ! echo height and elev angle
                hor_dist = (sqrt(4.*aterm*cterm + bterm**2.) - bterm)       
     1                                / (2.*aterm)
                range_thresh = 
     1              min(max(hor_dist,range_thresh_lo),range_thresh_hi)
            else
                range_thresh = range_thresh_hi
            endif

            if(rqc_2d(i,j) .eq. 0. .OR. 
     1         closest_radar(i,j) .gt. range_thresh)then
                do k = 1,NZ_L
                    if(ref_modelfg(i,j,k) .ne. r_missing_data)then
                        radar_ref_3d(i,j,k) = ref_modelfg(i,j,k)
                        rqc_2d(i,j) = 1.
                        n_fg_radar = n_fg_radar + 1
                        if(ref_modelfg(i,j,k) .gt. ref_base)then
                            n_fg_echoes = n_fg_echoes + 1
                        endif
                    endif ! first guess is present
                enddo ! k
            endif ! radar obs data are absent
        enddo ! j
        enddo ! i

        frac_fg_radar = float(n_fg_radar) / float(NX_L*NY_L*NZ_L)

        write(6,*)'First guess radar used over ',frac_fg_radar*100.
     1           ,'% of domain'

        write(6,*)'Number of first guess echoes is ',n_fg_echoes

C READ IN AND INSERT SAO DATA AS CLOUD SOUNDINGS
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

        var = 'TD'
        ext = 'lsx'
        call get_laps_2d(i4time,ext,var,units,comment
     1                  ,NX_L,NY_L,td_sfc_k,istatus)

        if(istatus .ne. 1)then
            write(6,*)' Error reading SFC Td - abort cloud analysis'
            goto999
        endif

        write(6,*)
        write(6,*)' Call Ingest/Insert SAO routines'
        n_cld_snd = 0
        call insert_sao(i4time,cldcv1,cf_modelfg,t_modelfg             ! I
     1  ,cld_hts,default_clear_cover,namelist_parms                    ! I
     1  ,lat,lon,topo,t_sfc_k,wtcldcv                                  ! I
     1  ,c1_name_array,l_perimeter,ista_snd
     1  ,cvr_snd,cld_snd,wt_snd,i_snd,j_snd,n_cld_snd,max_cld_snd
     1  ,NX_L,NY_L,KCLOUD                                              ! I
     1  ,n_obs_pos_b,lat_s,lon_s,c_stations    ! returned for precip type comp
     1  ,wx_s,t_s,td_s                         !    "      "    "     "
     1  ,elev_s                                ! ret for comparisons
     1  ,rad_s,solar_ea                        !  "   "       "   
     1  ,istat_sfc,maxstns,IX_LOW,IX_HIGH,IY_LOW,IY_HIGH)

        if(istat_sfc .ne. 1)then
            write(6,*)' No SAO data inserted: Aborting cloud analysis'
            goto999
        endif


C READ IN AND INSERT PIREP DATA AS CLOUD SOUNDINGS
        write(6,*)' Using Pireps stored in LAPS realtime system'

        call insert_pireps(i4time,cld_hts
     1      ,default_clear_cover
     1      ,cld_snd,wt_snd,i_snd,j_snd,n_cld_snd,max_cld_snd
     1      ,lat,lon,NX_L,NY_L,KCLOUD,IX_LOW,IX_HIGH,IY_LOW,IY_HIGH
     1      ,N_PIREP,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error: Bad status from insert_pireps,'
     1               ,' aborting cloud analysis'
            goto999
        endif
        I4_elapsed = ishow_timer()

C READ IN AND INSERT CO2 SLICING DATA AS CLOUD SOUNDINGS
        if(l_use_co2_mode1)then
            call insert_co2ctp(i4time,cld_hts,heights_3d                  ! I
     1            ,NX_L,NY_L,NZ_L,KCLOUD,r_missing_data                   ! I
     1            ,l_use_co2_mode1,latency_co2                            ! I
     1            ,default_clear_cover                                    ! I
     1            ,lat,lon,ix_low,ix_high,iy_low,iy_high                  ! I
     1            ,cld_snd,wt_snd,i_snd,j_snd,n_cld_snd,max_cld_snd       ! I/O
     1            ,istatus)                                               ! O
        endif

C DO ANALYSIS to horizontally spread SAO, PIREP, and optionally CO2 data
        write(6,*)
        if(l_use_co2_mode1)then
            write(6,*)' Analyzing SFC Obs, PIREP, and CO2-Slicing data'       
        else
            write(6,*)' Analyzing SFC Obs and PIREP data'
        endif

        max_obs = n_cld_snd * KCLOUD

!       Set weight for using model background clouds beyond a certain effective
!       radius of influence from the sfc obs/pireps
!       weight_modelfg = 0.    ! Model wt inactive, obs used to infinite radius
!       weight_modelfg = 1.    ! Model used beyond ~100km from nearest obs
!       weight_modelfg = .01   ! Model used beyond ~250km from nearest obs
!       weight_modelfg = .0001 ! Model used beyond ~630km from nearest obs
 
        weight_modelfg = cld_weight_modelfg 

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

!       Cloud cover QC check
        call qc_clouds_3d(clouds_3d,NX_L,NY_L,KCLOUD)

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
        call get_sat_data(i4time,i4_sat_window,i4_sat_window_offset,     ! I
     1                    NX_L,NY_L,r_missing_data,                      ! I
     1                    l_use_39,l_use_co2_mode2,latency_co2,          ! I
     1                    lat,lon,                                       ! I
     1                    tb8_k,istat_tb8,comment_tb8,                   ! O
     1                    t39_k,istat_t39,comment_t39,                   ! O
     1                    sst_k,istat_sst,comment_sst,                   ! O
     1                    cldtop_co2_pa_a,cloud_frac_co2_a,              ! O
     1                    istat_co2,lstat_co2_a)                         ! O

!       Calculate solar altitude
        do j = 1,NY_L
        do i = 1,NX_L
            call solar_position(lat(i,j),lon(i,j),i4time,solar_alt(i,j)
     1                                     ,solar_dec,solar_ha(i,j))
            call equ_to_altaz_d(solar_dec,solar_ha(i,j),lat(i,j)         ! I
     1                                   ,altdum,solar_az(i,j))          ! O
        enddo
        enddo

        icen = NX_L/2
        jcen = NY_L/2
        write(6,*)'solar dec/ha ',solar_dec,solar_ha(icen,jcen),
     1            ' at lat/lon ',lat(icen,jcen),lon(icen,jcen)
        write(6,*)'solar altitude = ',solar_alt(icen,jcen)

!       Cloud cover QC check
        call qc_clouds_3d(clouds_3d,NX_L,NY_L,KCLOUD)

        call get_vis(i4time,solar_alt,l_use_vis,l_use_vis_add            ! I
     1              ,l_use_vis_partial,lat,lon                           ! I
     1              ,i4_sat_window,i4_sat_window_offset                  ! I
     1              ,rlaps_land_frac,topo                                ! I
     1              ,cloud_frac_vis_a,albedo,ihist_alb                   ! O
     1              ,static_albedo                                       ! O
     1              ,subpoint_lat_clo_vis,subpoint_lon_clo_vis           ! O 
     1              ,comment_alb                                         ! O
     1              ,NX_L,NY_L,KCLOUD,r_missing_data                     ! O
     1              ,istat_vis_potl_a,istat_vis)                         ! O

        call get_istat_39(t39_k,tb8_k,solar_alt,r_missing_data           ! I
     1                   ,rlaps_land_frac,NX_L,NY_L                      ! I
     1                   ,static_albedo                                  ! I
     1                   ,istat_39_a)                                    ! O

!       Cloud cover QC check
        call qc_clouds_3d(clouds_3d,NX_L,NY_L,KCLOUD)

        call get_pres_3d(i4time,NX_L,NY_L,NZ_L,pres_3d,istatus)

        if(i_varadj .eq. 1)then
            write(6,*)' Call cloud_var before insert_sat'
            call cloud_var(i4time,lat,lon
     1                    ,NX_L,NY_L,NZ_L,KCLOUD,heights_3d,temp_3d
     1                    ,t_gnd_k,clouds_3d,cld_hts,tb8_k
     1                    ,cloud_frac_vis_a
     1                    ,subpoint_lat_clo_vis,subpoint_lon_clo_vis  ! I 
     1                    ,r_missing_data                             ! I
     1                    ,di_dh,dj_dh)                               ! O
            di_dh = 0. ! for testing
            dj_dh = 0. ! for testing
        endif

        I4_elapsed = ishow_timer()

        call insert_sat(i4time,clouds_3d,cldcv_sao,cld_hts,lat,lon,
     1       pct_req_lvd_s8a,default_clear_cover,                       ! I
     1       tb8_k,istat_tb8,                                           ! I
     1       sst_k,istat_sst,                                           ! I
     1       istat_39_a, l_use_39,                                      ! I
     1       di_dh,dj_dh,                                               ! I
     1       istat_39_add_a,                                            ! O
     1       tb8_cold_k,                                                ! O
     1       grid_spacing_cen_m,surface_sao_buffer,                     ! I
     1       cloud_frac_vis_a,istat_vis_potl_a,                         ! I
     1       istat_vis_added_a,                                         ! O
     1       solar_alt,solar_ha,solar_dec,                              ! I
     1       lstat_co2_a, cloud_frac_co2_a, cldtop_co2_pa_a,            ! I
     1       rlaps_land_frac,                                           ! I
     1       topo,heights_3d,temp_3d,t_sfc_k,td_sfc_k,pres_sfc_pa,      ! I
     1       t_modelfg,sh_modelfg,pres_3d,                              ! I
     1       cvr_snow,NX_L,NY_L,KCLOUD,NZ_L,r_missing_data,             ! I
     1       t_gnd_k,                                                   ! O
     1       cldtop_co2_m,cldtop_tb8_m,cldtop_m,ht_sao_top,             ! O
     1       istatus)                                                   ! O

        if(istatus .ne. 1)then
            write(6,*)' Error: Bad status returned from insert_sat'
            goto999
        endif

!       Cloud cover QC check
        call qc_clouds_3d(clouds_3d,NX_L,NY_L,KCLOUD)

        write(6,*)' Cloud top (Band 8 vs. ht_sao_top)'
        scale = .0001
        write(6,301)
301     format('  Cloud Top (km)             Band 8                ',
     1                 23x,'               ht_sao_top')
        CALL ARRAY_PLOT(cldtop_tb8_m,ht_sao_top,NX_L,NY_L,'HORZ CV'
     1                 ,c1_name_array(:,:,1),KCLOUD,cld_hts,scale)

        write(6,*)' Cloud top (Band 8 vs. Satellite Analysis)'
        scale = .0001
        write(6,302)
302     format('  Cloud Top (km)             Band 8                ',
     1                 23x,'          Satellite Analysis')
        CALL ARRAY_PLOT(cldtop_tb8_m,cldtop_m,NX_L,NY_L,'HORZ CV'
     1                 ,c1_name_array(:,:,1),KCLOUD,cld_hts,scale)

        I4_elapsed = ishow_timer()

        n_radar_2dref=0
        n_radar_3dref=0
        n_radar_3dref_orig=0

C       THREE DIMENSIONALIZE RADAR DATA IF NECESSARY (E.G. NOWRAD)
        write(6,*)' Three dimensionalizing radar data'

!       Clear out radar echo above the highest cloud top
        k_ref_def = nint(zcoord_of_pressure(float(700*100)))

        do j = 1,NY_L
        do i = 1,NX_L
            
            if(istat_radar_2dref_a(i,j) .eq. 1 .and.
     1         istat_radar_3dref_a(i,j) .eq. 0       )then

                if(radar_ref_3d(i,j,1) .gt. ref_base)then

                    k_topo = int(zcoord_of_pressure(pres_sfc_pa(i,j)))
                    k_ref = max(k_ref_def,k_topo+2)

!                   k_ref = k_ref_def

                    cloud_top_m = default_top

!                   Test for Cloud Top
                    do k = KCLOUD-1,1,-1
                        if(clouds_3d(i,j,k  ) .gt. thresh_cvr_ceiling 
     1                                       .AND.       
     1                     clouds_3d(i,j,k+1) .le. thresh_cvr_ceiling 
     1                                                            )then 
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
                        k_cloud_top = nint(height_to_zcoord2(cloud_top_m
     1                          ,heights_3d,NX_L,NY_L,NZ_L,i,j,istatus))      
                        if(istatus .ne. 1)then
                            write(6,*)' Error: Bad status returned'
     1                          ,' from height_to_zcoord2'
     1                          ,cloud_top_m,heights_3d(i,j,NZ_L),i,j       
                            goto999
                        endif

                        k_cloud_top = max(k_ref,k_cloud_top)
                    else
                        if(l_trust_narrowband)then
                            k_cloud_top = k_ref
                        else
                            k_cloud_top = 1 ! Will set radar to no echo
                        endif
                    endif

                    do k = NZ_L,k_cloud_top,-1
                        radar_ref_3d(i,j,k) = ref_base
                    enddo ! k

                endif ! Radar echo at this grid point

!               if(k_cloud_top .gt. 1)then
!               We have three-dimensionalized this grid point
                n_radar_2dref = n_radar_2dref + 1
                n_radar_3dref = n_radar_3dref + 1
                istat_radar_3dref_a(i,j) = 1
!               endif

            elseif(istat_radar_2dref_a(i,j) .eq. 1 .and.
     1             istat_radar_3dref_a(i,j) .eq. 1       )then

!               Grid point is already fully three dimensional
                n_radar_2dref = n_radar_2dref + 1
                n_radar_3dref = n_radar_3dref + 1
                n_radar_3dref_orig = n_radar_3dref_orig + 1
                lstat_radar_3dref_orig_a(i,j) = .true.

            endif ! Is this grid point 2-d or 3-d?

        enddo ! i
        enddo ! j

        if(n_radar_2dref .ge. 1)then
            istat_radar_2dref = 1
        else
            istat_radar_2dref = 0
        endif

        if(n_radar_3dref .ge. 1)then
            istat_radar_3dref = 1
        else
            istat_radar_3dref = 0
        endif

!       Note, if orig data is a mixture of 2d and 3d, the istat here gets set
!       to 3d.
        if(n_radar_3dref_orig .ge. 1)then
            istat_radar_3dref_orig = 1
        else
            istat_radar_3dref_orig = 0
        endif

        write(6,*)' n_radar 2dref/3dref_orig/3dref = ' 
     1           ,n_radar_2dref,n_radar_3dref_orig,n_radar_3dref

        write(6,*)' istat_radar 2dref/3dref_orig/3dref = ' 
     1           ,istat_radar_2dref,istat_radar_3dref_orig
     1           ,istat_radar_3dref

        I4_elapsed = ishow_timer()

!       Generate field of radar coverage (2D/3D)
        do i = 1,NX_L
        do j = 1,NY_L
            if(rqc_2d(i,j) .eq. 1.)then
                plot_maskr(i,j) = 10.0 ! 3D first guess data
            elseif(lstat_radar_3dref_orig_a(i,j))then
                plot_maskr(i,j) = 30.0 ! 3D original data
            elseif(istat_radar_2dref_a(i,j) .eq. 1)then
                plot_maskr(i,j) = 20.0 ! 2D original data
            else                    
                plot_maskr(i,j) = 0.0  ! no data
            endif
        enddo ! j
        enddo ! i

C INSERT RADAR DATA
        if(n_radar_3dref .gt. 0)then
            call get_max_reflect(radar_ref_3d,NX_L,NY_L,NZ_L,ref_base       
     1                          ,dbz_max_2d)
 
!           Generate ASCII prelim reflectivity plot
            write(6,1201)
1201        format('                  Radar Coverage',50x
     1                              ,'Prelim Radar Max Reflectivity')       
            scale = 0.01
            CALL ARRAY_PLOT(plot_maskr,dbz_max_2d,NX_L,NY_L,'HORZ CV'
     1                     ,c1_name_array(:,:,1),KCLOUD,cld_hts,scale) ! radar 

            call insert_radar(i4time,clouds_3d,cld_hts
     1          ,temp_3d,t_sfc_k,td_sfc_k                            ! I
     1          ,grid_spacing_cen_m,NX_L,NY_L,NZ_L                   ! I
     1          ,KCLOUD,cloud_base,ref_base                          ! I
     1          ,topo,solar_alt,r_missing_data                       ! I
     1          ,radar_ref_3d,dbz_max_2d                             ! I/O
     1          ,vis_radar_thresh_dbz                                ! I
     1          ,l_unresolved                                        ! O
     1          ,heights_3d                                          ! I
     1          ,istatus)                                            ! O

            if(istatus .ne. 1)then
                write(6,*)
     1          ' Error: Bad status returned from insert_radar'      
                goto999
            endif

        endif

        if(i_varadj .eq. 1)then
            write(6,*)' Call cloud_var before insert vis'
            call cloud_var(i4time,lat,lon
     1                    ,NX_L,NY_L,NZ_L,KCLOUD,heights_3d,temp_3d
     1                    ,t_gnd_k,clouds_3d,cld_hts,tb8_k
     1                    ,cloud_frac_vis_a
     1                    ,subpoint_lat_clo_vis,subpoint_lon_clo_vis  ! I 
     1                    ,r_missing_data                             ! I
     1                    ,di_dh,dj_dh)                               ! O
            di_dh = 0. ! for testing
            dj_dh = 0. ! for testing
        endif

        I4_elapsed = ishow_timer()

C       INSERT VISIBLE / 3.9u SATELLITE IN CLEARING STEP
        if(istat_vis .eq. 1 .OR. (istat_t39 .eq. 1 .and. l_use_39) )then
            call insert_vis(i4time,clouds_3d,cld_hts
     1        ,topo,cloud_frac_vis_a,albedo,ihist_alb                 ! I
     1        ,istat_39_a,l_use_39                                    ! I
     1        ,NX_L,NY_L,KCLOUD,r_missing_data                        ! I
     1        ,vis_radar_thresh_cvr,vis_radar_thresh_dbz              ! I
     1        ,istat_radar_3dref,radar_ref_3d,NZ_L,ref_base
     1        ,solar_alt,solar_az                                     ! I
     1        ,di_dh,dj_dh                                            ! I
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

C ASCII PLOTS in HORIZONTAL AND VERTICAL SLICES

C       HORIZONTAL SLICES

        DO K=4,KCLOUD,4
            CALL SLICE(cf_modelfg,NX_L,NY_L,KCLOUD,CVHZ1
     1                ,NX_L,NY_L,1,0,0,K,0,0)
            CALL SLICE(cldcv_sao,NX_L,NY_L,KCLOUD,cvr_max
     1                ,NX_L,NY_L,1,0,0,K,0,0)
            write(6,401)k,cld_hts(k)
401         format(4x,'Lvl',i4,f8.0,' m     Model First Guess Only',
     1            20x,'              With Point Data Added')
            scale = 1.
            CALL ARRAY_PLOT(CVHZ1,cvr_max,NX_L,NY_L,'HORZ CV'
     1                     ,c1_name_array(:,:,K),KCLOUD,cld_hts,scale)
        ENDDO ! k

        DO K=4,KCLOUD,4
            CALL SLICE(cldcv_sao,NX_L,NY_L,KCLOUD,CVHZ1
     1                ,NX_L,NY_L,1,0,0,K,0,0)
            CALL SLICE(clouds_3d,NX_L,NY_L,KCLOUD,cvr_max
     1                ,NX_L,NY_L,1,0,0,K,0,0)
            write(6,402)k,cld_hts(k)
402         format(4x,'Lvl',i4,f8.0,' m     Before Satellite/Radar',
     1            20x,'              After Satellite/Radar')
            scale = 1.
            CALL ARRAY_PLOT(CVHZ1,cvr_max,NX_L,NY_L,'HORZ CV'
     1                     ,c1_name_array(:,:,K),KCLOUD,cld_hts,scale)
        ENDDO ! k

C       EW SLICES
        write(6,*)
        DO J=10,NY_L,20
            CALL SLICE(cf_modelfg,NX_L,NY_L,KCLOUD,CVEW1
     1                ,NX_L,KCLOUD,0,1,0,0,J,0)
            CALL SLICE(cldcv_sao,NX_L,NY_L,KCLOUD,CVEW2
     1                ,NX_L,KCLOUD,0,1,0,0,J,0)
            write(6,501)j
501         format(5x,'  J =',i4,10x,'  Model First Guess Only       ',
     1          21x,'           With Point Data Added')

            scale = 1.
            CALL ARRAY_PLOT(CVEW1,cvew2,NX_L,KCLOUD,'VERT CV'
     1                     ,c1_name_array,KCLOUD,cld_hts,scale)
        ENDDO ! j

C       EW SLICES
        write(6,*)
        DO J=10,NY_L,20
            CALL SLICE(cldcv_sao,NX_L,NY_L,KCLOUD,CVEW1
     1                ,NX_L,KCLOUD,0,1,0,0,J,0)
            CALL SLICE(clouds_3d,NX_L,NY_L,KCLOUD,CVEW2
     1                ,NX_L,KCLOUD,0,1,0,0,J,0)
            write(6,502)j
502         format(5x,'  J =',i4,10x,'  Before Satellite/Radar       ',
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

!       Log info on cloud holes with cover at least .10 less than all neighbors
        do i = 2,NX_L-1
        do j = 2,NY_L-1
            isign_hole = 0
            do ii = i-1,i+1
            do jj = j-1,j+1
                if(ii .ne. i .or. jj .ne. j)then
                    diff = cvr_max(i,j) - cvr_max(ii,jj)
                    isign_hole = isign_hole + int(sign(1.0,diff+.10))

                endif ! not at the center point

            enddo ! jj
            enddo ! ii

            if(isign_hole .eq. -8)then
                write(6,*)' hole detected ',i,j,cvr_max(i,j)

                do jj = j+1,j-1,-1
                    write(6,511,err=512)
     1                    nint(cvr_max(i-1,jj)*100)
     1                   ,nint(cvr_max(i  ,jj)*100)
     1                   ,nint(cvr_max(i+1,jj)*100)
     1                   ,nint2(tb8_k(i-1,jj),1),nint2(tb8_k(i,jj),1)      
     1                   ,nint2(tb8_k(i+1,jj),1)
     1                   ,nint(t_gnd_k(i-1,jj)),nint(t_gnd_k(i,jj))        
     1                   ,nint(t_gnd_k(i+1,jj))
     1                   ,nint2(tb8_k(i-1,jj)-t_gnd_k(i-1,jj),1)
     1                   ,nint2(tb8_k(i  ,jj)-t_gnd_k(i  ,jj),1)
     1                   ,nint2(tb8_k(i+1,jj)-t_gnd_k(i+1,jj),1)
     1                   ,nint(topo(i-1,jj)),nint(topo(i,jj))        
     1                   ,nint(topo(i+1,jj))
     1                   ,nint(cvr_sao_max(i-1,jj)*100)
     1                   ,nint(cvr_sao_max(i  ,jj)*100)
     1                   ,nint(cvr_sao_max(i+1,jj)*100)
     1                   ,nint2(cloud_frac_vis_a(i-1,jj),100)
     1                   ,nint2(cloud_frac_vis_a(i  ,jj),100)
     1                   ,nint2(cloud_frac_vis_a(i+1,jj),100)
 511                format(1x,3i3,4x,3i4,4x,3i4,4x,3i4,4x,3i5,4x,3i3
     1                    ,4x,3i3)
 512            enddo ! jj

            endif ! cloud hole detected

        enddo ! j
        enddo ! i

        if(istat_radar_3dref .eq. 1)then
            call get_max_reflect(radar_ref_3d,NX_L,NY_L,NZ_L
     1                          ,r_missing_data,dbz_max_2d)

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
     1                 ,c1_name_array(:,:,1),KCLOUD,cld_hts,scale)

        write(6,701)
701     format('  Max Cloud Cover           VISIBLE SATELLITE     ',
     1            20x,'              Final Analysis')
        scale = 1.
        CALL ARRAY_PLOT(cloud_frac_vis_a,cvr_max,NX_L,NY_L,'HORZ CV'
     1                 ,c1_name_array(:,:,1),KCLOUD,cld_hts,scale)

        I4_elapsed = ishow_timer()

        do k = 1,NZ_L
            write(6,101)k,heights_3d(NX_L/2,NY_L/2,k)
101         format(1x,i3,f10.2)
        enddo ! k

        I4_elapsed = ishow_timer()

!       Write out Cloud Field
        do k=1,KCLOUD
            rlevel = height_to_zcoord2(cld_hts(k),heights_3d
     1                   ,NX_L,NY_L,NZ_L,NX_L/2,NY_L/2,istatus)   

            if(rlevel .le. float(NZ_L))then
                cld_pres_1d(k) = pressure_of_rlevel(rlevel)
            else
                cld_pres_1d(k) = 0. ! Cloud height level above pressure grid
            endif

            write(6,102)k,cld_hts(k),cld_pres_1d(k),istatus
102         format(1x,i3,2f10.1,i2)

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
                    cloud_base(i,j) = max(cloud_base(i,j),topo(i,j))
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
                    cloud_ceiling(i,j) = max(cloud_ceiling(i,j),0.)
                endif
              enddo ! k

            endif ! Clouds exist in this column

        enddo ! i
        enddo ! j

!       Calculate cloud analysis implied snow cover
        call cloud_snow_cvr(cvr_max,cloud_frac_vis_a,cldtop_tb8_m
     1          ,tb8_k,NX_L,NY_L,r_missing_data,cvr_snow_cycle)

!       MORE ASCII PLOTS
        write(6,801)
801     format('                            VISIBLE SATELLITE     ',
     1            20x,'      csc  (cycle)  Snow Cover')
        scale = 1.
        CALL ARRAY_PLOT(cloud_frac_vis_a,cvr_snow_cycle,NX_L,NY_L
     1                 ,'HORZ CV',c1_name_array(:,:,1),KCLOUD,cld_hts
     1                 ,scale)       

        write(6,901)
901     format('                     lm2 (overall) Snow Cover      ',
     1            20x,'      csc  (cycle)  Snow Cover')
        scale = 1.
        CALL ARRAY_PLOT(cvr_snow,cvr_snow_cycle,NX_L,NY_L,'HORZ CV'
     1                  ,c1_name_array(:,:,1),KCLOUD,cld_hts,scale)

        do i = 1,NX_L
        do j = 1,NY_L
            if(cldtop_m(i,j)  .ne. r_missing_data .and.
     1         cvr_max(i,j)   .ge. 0.1                            )then       
                plot_mask(i,j) = cloud_top(i,j)            ! Set cloud top mask
            else
                plot_mask(i,j) = r_missing_data
            endif
        enddo ! j
        enddo ! i

        write(6,*)' Cloud top (Band 8 vs. Final Analysis)'
        scale = .0001
        write(6,1001)
1001    format('  Cloud Top (km)             Band 8                ',
     1            20x,'              Final Analysis')
        CALL ARRAY_PLOT(cldtop_tb8_m,plot_mask,NX_L,NY_L,'HORZ CV'
     1                  ,c1_name_array(:,:,1),KCLOUD,cld_hts,scale) ! Plot Band 8 mask


        write(6,1101)
1101    format('  Max Cloud Cover   3.9u     (nocld:cld:added=3:7:9)',      
     1       18x,'              Final Analysis')

!       Set 3.9 micron plot mask
        do i = 1,NX_L
        do j = 1,NY_L
            if(istat_39_a(i,j) .eq. -1)then
                plot_mask(i,j) = 0.3 
            elseif(istat_39_a(i,j) .eq. 0)then
                plot_mask(i,j) = 0.0 
            elseif(istat_39_a(i,j) .eq. 1)then
                if(istat_39_add_a(i,j) .eq. 1)then
                    plot_mask(i,j) = 0.9
                else
                    plot_mask(i,j) = 0.7
                endif
            endif
        enddo ! j
        enddo ! i

        scale = 1.
        CALL ARRAY_PLOT(plot_mask,cvr_max,NX_L,NY_L,'HORZ CV'
     1                 ,c1_name_array(:,:,1),KCLOUD,cld_hts,scale) ! 3.9u mask



!       Write out LCB file (Cloud Base, Top, and Ceiling fields)
        if(iwrite_output .ge. 1)then
            ext = 'lcb'
            call get_directory(ext,directory,len_dir)

            call move(cloud_base    ,out_array_3d(1,1,1),NX_L,NY_L)
            call move(cloud_top     ,out_array_3d(1,1,2),NX_L,NY_L)
            call move(cloud_ceiling ,out_array_3d(1,1,3),NX_L,NY_L)

            call put_clouds_2d(i4time,directory,ext,NX_L,NY_L
     1                                  ,out_array_3d,istatus)
            if(istatus .eq. 1)j_status(n_lcb) = ss_normal
        endif ! iwrite_output

!       This is where we will eventually split the routines, additional data
!       is necessary for more derived fields

        if(n_radar_3dref .gt. 0)then ! Write out data (lps - radar_ref_3d)
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

!           Generate ASCII final reflectivity plot
            write(6,1202)
1202        format('                  Radar Coverage',50x
     1                              ,'Final Radar Max Reflectivity')
            scale = 0.01
            CALL ARRAY_PLOT(plot_maskr,dbz_max_2d,NX_L,NY_L,'HORZ CV'
     1                     ,c1_name_array(:,:,1),KCLOUD,cld_hts,scale) ! radar 

        endif ! n_radar_3dref

        call compare_analysis_to_saos(NX_L,NY_L,cvr_sao_max
     1  ,cloud_frac_vis_a,tb8_k,t_gnd_k,t_sfc_k,cvr_max,r_missing_data
     1  ,dbz_max_2d,cld_snd,ista_snd,max_cld_snd,cld_hts,KCLOUD
     1  ,n_cld_snd,c_stations,lat_s,lon_s,elev_s,maxstns)

!       Reread solar data from latest LSO file
        call read_surface_sa(i4time,maxstns,                       ! I
     1   n_obs_b,c_stations,reptype,atype,                         ! O
     1   lat_s,lon_s,elev_s,wx_s,t_s,td_s,                         ! O
     1   kloud,store_amt,store_hgt,                                ! O
     1   rad_s,solar_ea,obstime,istatus)                           ! O

        call compare_analysis_to_rad(i4time,NX_L,NY_L,cvr_sao_max
     1  ,solar_alt,cvr_snow
     1  ,cloud_frac_vis_a,tb8_k,t_gnd_k,td_sfc_k,cvr_max,r_missing_data
     1  ,dbz_max_2d,cld_snd,ista_snd,max_cld_snd,cld_hts,KCLOUD
     1  ,rad_s,n_cld_snd,c_stations,lat_s,lon_s,elev_s
     1  ,maxstns,n_obs_b,swi_2d)

!       Write LCV file
        if(iwrite_output .ge. 1)then
            do i = 1,NX_L
            do j = 1,NY_L
                cvr_water_temp(i,j) = r_missing_data
            enddo ! j
            enddo ! i

            ext = 'lcv'
            var_a(1) = 'LCV'
            var_a(2) = 'CSC'
            var_a(3) = 'CWT'
            var_a(4) = 'S8A'
            var_a(5) = 'S3A'
            var_a(6) = 'ALB'
            var_a(7) = 'RQC'
            var_a(8) = 'SWI'
            units_a(1) = 'UNDIM'
            units_a(2) = 'UNDIM'
            units_a(3) = 'K'
            units_a(4) = 'K'
            units_a(5) = 'K'
            units_a(6) = ' '
            units_a(7) = ' '
            units_a(8) = 'W/M**2'
            comment_a(1) = 'LAPS Cloud Cover'
            comment_a(2) = 'LAPS Cloud Analysis Implied Snow Cover'
            comment_a(3) = 'LAPS Clear Sky Water Temp'
            comment_a(4) = comment_tb8
            comment_a(5) = comment_t39
            comment_a(6) = comment_alb
            comment_a(7) = 'LAPS Radar Quality'
            comment_a(8) = 'Downward Solar Radiation'


            call move(cvr_max       ,out_array_3d(1,1,1),NX_L,NY_L)
            call move(cvr_snow_cycle,out_array_3d(1,1,2),NX_L,NY_L)
            call move(cvr_water_temp,out_array_3d(1,1,3),NX_L,NY_L)
            call move(tb8_k         ,out_array_3d(1,1,4),NX_L,NY_L)
            call move(t39_k         ,out_array_3d(1,1,5),NX_L,NY_L)
            call move(albedo        ,out_array_3d(1,1,6),NX_L,NY_L)
            call move(rqc_2d        ,out_array_3d(1,1,7),NX_L,NY_L)
            call move(swi_2d        ,out_array_3d(1,1,8),NX_L,NY_L)

            call put_laps_multi_2d(i4time,ext,var_a,units_a,
     1              comment_a,out_array_3d,NX_L,NY_L,8,istatus)

            if(istatus .eq. 1)j_status(n_lcv) = ss_normal
        endif ! iwrite_output

500     continue

        l_get_cloudtype = .false.

        if(l_get_cloudtype)then ! Compute just the cloud type
            l_flag_cloud_type = .true.
            l_flag_mvd = .false.
            l_flag_icing_index = .false.
            l_flag_bogus_w = .false.
            l_deep_vv = .false.
!           call get_cloud_deriv(NX_L,NY_L,NZ_L,clouds_3d,cld_hts         ! I
!!   1                          temp_3d,rh_3d_dum,heights_3d,pres_3d,     ! I
!!   1                          istat_radar,radar_3d,grid_spacing_cen_m,  ! I  
!    1                          l_mask_pcptype,                           ! O
!    1                          ibase_array,itop_array,                   ! O
!!   1                          iflag_slwc,slwc_3d,cice_3d,
!!   1                          thresh_cvr_cty_vv,thresh_cvr_lwc,
!!   1                          l_flag_cloud_type,cldpcp_type_3d,         ! I/O
!!   1                          l_flag_mvd,mvd_3d,
!!   1                          l_flag_icing_index,icing_index_3d,
!!   1                          vv_to_height_ratio_Cu,                    ! I
!!   1                          vv_to_height_ratio_Sc,                    ! I
!!   1                          vv_for_St,                                ! I
!!   1                          l_flag_bogus_w,omega_3d,l_bogus_radar_w,
!    1                          l_deep_vv,                                ! I
!    1                          twet_snow_dum,                            ! I
!!   1                          l_flag_pcp_type,                          ! I
!!   1                          istatus)                                  ! O

!           Write CTY field          

        endif

!       If needed 't_modelfg' and 'sh_modelfg' are available.
!       'rh_modelfg' is sometimes calculated in 'get_model_fg'

        if(i_varadj .eq. 1)then
            write(6,*)' Call cloud_var at end of analysis'
            call cloud_var(i4time,lat,lon
     1                    ,NX_L,NY_L,NZ_L,KCLOUD,heights_3d,temp_3d
     1                    ,t_gnd_k,clouds_3d,cld_hts,tb8_k
     1                    ,cloud_frac_vis_a
     1                    ,subpoint_lat_clo_vis,subpoint_lon_clo_vis  ! I 
     1                    ,r_missing_data                             ! I
     1                    ,di_dh,dj_dh)                               ! O
        endif

999     continue

        write(6,*)' Notifications'
        do i = iprod_start,iprod_end
            write(6,*)' ',exts(i),' ',j_status(i),' ',i
        enddo ! i

        write(6,*)' End of Cloud Analysis Package'

        deallocate(cldcv1)
        deallocate(cf_modelfg)
        deallocate(t_modelfg)
        deallocate(sh_modelfg)
        deallocate(ref_modelfg)
        deallocate(cld_snd)
        deallocate(wt_snd)

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
        character*4 LVL_COORD_2d(nfields)

        real field_2dcloud(imax,jmax,nfields)

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

        subroutine put_clouds_3d(i4time,ext,clouds_3d,cld_hts
     1                          ,cld_pres_1d,ni,nj,nk,istatus)

!       1997 Jul 31 K. Dritz  - Removed include of lapsparms.for, which was
!                               not actually needed for anything.

        character*150 DIRECTORY
        character*31 ext

        integer NZ_CLOUD_MAX
        parameter (NZ_CLOUD_MAX = 42)

        character*125 comment_3d(NZ_CLOUD_MAX)
        character*10 units_3d(NZ_CLOUD_MAX)
        character*3 var_3d(NZ_CLOUD_MAX),var_2d
        integer LVL_3d(NZ_CLOUD_MAX)
        character*4 LVL_COORD_3d(NZ_CLOUD_MAX)

        real clouds_3d(ni,nj,nk)
        real cld_hts(nk)
        real cld_pres_1d(nk)

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


        subroutine cloud_snow_cvr(cvr_max,cloud_frac_vis_a,cldtop_tb8_m
     1             ,tb8_k,ni,nj,r_missing_data,cvr_snow_cycle)

        real cvr_max(ni,nj)          ! Input
        real cloud_frac_vis_a(ni,nj) ! Input
        real cldtop_tb8_m(ni,nj)     ! Input
        real tb8_k(ni,nj)            ! Input
        real cvr_snow_cycle(ni,nj)   ! Output

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
!    1          .and. cldtop_tb8_m(i,j)     .eq. r_missing_data ! No Band 8 clds
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


        subroutine compare_analysis_to_saos(ni,nj,cvr_sao_max
     1  ,cloud_frac_vis_a,tb8_k,t_gnd_k,t_sfc_k,cvr_max,r_missing_data
     1  ,dbz_max_2d,cld_snd,ista_snd,max_cld_snd,cld_hts,KCLOUD
     1  ,n_cld_snd,c_stations,lat_s,lon_s,elev_s,maxstns)

        real cloud_frac_vis_a(ni,nj),tb8_k(ni,nj),t_gnd_k(ni,nj)
     1        ,t_sfc_k(ni,nj),cvr_max(ni,nj),cvr_sao_max(ni,nj)
     1        ,dbz_max_2d(ni,nj)

        real cld_snd(max_cld_snd,KCLOUD)
        integer ista_snd(max_cld_snd)
        real cld_hts(KCLOUD)

        character c_stations(maxstns)*(*)
        real lat_s(maxstns), lon_s(maxstns), elev_s(maxstns)

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
                        write(6,*)'Sta   i    j   VIS frac tb8_k  '
     1                  //'t_gnd_k t_sfc_k  cldsnd cv_sa_mx cvr_mx '
     1                  //'snd-ht  dbz        9pt    25pt'
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
     1                           ,i_i,i_j
     1                           ,cloud_frac_vis_a(i_i,i_j)
     1                           ,tb8_k(i_i,i_j)
     1                           ,t_gnd_k(i_i,i_j)
     1                           ,t_sfc_k(i_i,i_j)
     1                           ,cld_snd_max
     1                           ,cvr_sao_max(i_i,i_j)
     1                           ,cvr_max(i_i,i_j)
     1                           ,height_max
     1                           ,dbz_max_2d(i_i,i_j)
     1                           ,c3_discrep
     1                           ,cvr_9pt
     1                           ,cvr_25pt
     1                           ,c1_c
1111                format(1x,a3,2i5,f8.2,3f8.1,3f8.2,f8.0,f5.0       
     1                    ,a3,f8.2,f7.2,1x,a1)

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

        real radar_ref_3d(ni,nj,nk)                   ! I
        real dbz_max_2d(ni,nj)                        ! I
        real cvr_max(ni,nj)                           ! I
        real cloud_frac_vis_a(ni,nj)                  ! I

!       This routine compares the cloud and radar fields and flags
!       remaining differences that weren't caught in earlier processing

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

              iblank_radar = 0

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

!                     We will block out the radar
                      write(6,2)i,j,cvr_max(i,j),dbz_max_2d(i,j)
     1                       ,cloud_frac_vis_a(i,j)
2                     format(' CLD_RDR: cvr/dbz/vis <',2i4,f8.2
     1                      ,f8.1,f8.2,' Blank out radar')

                      iblank_radar = 1

                  else ! Radar is too strong to block out

                      write(6,3)i,j,cvr_max(i,j),dbz_max_2d(i,j)
     1                       ,cloud_frac_vis_a(i,j)
3                     format(' CLD_RDR: cvr/dbz/vis >',2i4,f8.2
     1                                                 ,f8.1,f8.2)

                  endif ! ref < vis

              elseif(cvr_max(i,j) .gt. cloud_frac_vis_a(i,j))then

!                 Don't know if VIS lowered cloud cover        CVR > VIS

                  if(dbz_max_2d(i,j) .lt. vis_radar_thresh_dbz)then

!                     We will block out the radar
                      write(6,4)i,j,cvr_max(i,j),dbz_max_2d(i,j)
     1                       ,cloud_frac_vis_a(i,j)
4                     format(' ???_RDR: cvr/dbz/vis <',2i4,f8.2
     1                      ,f8.1,f8.2,' Blank out radar')

                      iblank_radar = 1

                  else ! Radar is too strong to block out

                      write(6,5)i,j,cvr_max(i,j),dbz_max_2d(i,j)
     1                       ,cloud_frac_vis_a(i,j)
5                     format(' ???_RDR: cvr/dbz/vis >',2i4,f8.2
     1                                                 ,f8.1,f8.2)

                  endif ! ref < vis

              endif ! cover compare to vis

              if(iblank_radar .eq. 1)then        ! Take action
                  dbz_max_2d(i,j) = ref_base
                  do kl = 1,nk
                      radar_ref_3d(i,j,kl) = ref_base
                  enddo ! kl
              endif

            endif ! radar echo with low cloud cover
        enddo ! j
        enddo ! i

        return
        end

        function nint2(x,ifactor)

!       This routine helps scale the arguments for ASCII debug printing

        call get_r_missing_data(r_missing_data,istatus)

        if(x .ne. r_missing_data)then
            nint2 = nint(x*float(ifactor))
        else
            nint2 = 999999
        endif

        return
        end

        subroutine qc_clouds_3d(clouds_3d,NX_L,NY_L,KCLOUD)

        real clouds_3d(NX_L,NY_L,KCLOUD)
        logical l_poss_extrap ! used to allow for edge effects from 'barnes_r5'

        nskip_max = 4 ! 'See barnes_r5'

        I4_elapsed = ishow_timer()

        write(6,*)' Subroutine qc_clouds_3d...'

        do i = 1,NX_L
        do j = 1,NY_L
            if(NX_L-i .le. nskip_max .or. NY_L-j .le. nskip_max)then
                l_poss_extrap = .true. ! Extrapolation edge effects possible
            else
                l_poss_extrap = .false.
            endif

            do k = 1,KCLOUD
!               call qc_clouds_0d(i,j,k,clouds_3d(i,j,k)
!    1                           ,NX_L,NY_L,l_poss_extrap)

!               Subroutine code reproduced in calling routine for efficiency
                clouds_0d = clouds_3d(i,j,k)

                if(clouds_0d .gt. 1.0)then 
                    if(.not. l_poss_extrap)then
                        if(clouds_0d .gt. 1.001)then
                            write(6,*)
     1                          ' Error, clouds_0d > 1',i,j,k,clouds_0d       
                            stop
                        else ! just over 1.0 with no edge effect
                            write(6,*)
     1                      ' Warning, clouds_0d > 1 - reset'
     1                      ,i,j,k,clouds_0d
                            clouds_0d = 1.0
                        endif
                    else
                        write(6,*)
     1                 ' Warning, clouds_0d > 1 - reset for edge effect'       
     1                  ,i,j,k,clouds_0d
                        clouds_0d = 1.0
                    endif

                elseif(clouds_0d .lt. 0.0)then 
                    if(l_poss_extrap)then
                        qc_thr = -200.0
                    else
                        qc_thr = -0.0005
                    endif
                    if(clouds_0d .lt. qc_thr)then
                        write(6,*)' Error, clouds_0d << 0',i,j,k
     1                                    ,clouds_0d   
                        stop
                    else 
                        write(6,*)
     1                 ' Warning, clouds_0d < 0 - reset for edge effect'       
     1                  ,i,j,k,clouds_0d
                        clouds_0d = 0.
                    endif
                endif

                clouds_3d(i,j,k) = clouds_0d

            enddo ! k
        enddo ! j
        enddo ! i

        I4_elapsed = ishow_timer()

        return
        end

        subroutine qc_clouds_0d(i,j,k,clouds_3d
     1                         ,NX_L,NY_L,l_poss_extrap)

        real clouds_3d
        logical l_poss_extrap ! used to allow for edge effects from 'barnes_r5'

        if(clouds_3d .gt. 1.0)then 
            if(.not. l_poss_extrap)then
                if(clouds_3d .gt. 1.001)then
                    write(6,*)' Error, clouds_3d > 1',i,j,k,clouds_3d
                    stop
                else ! just over 1.0 with no edge effect
                    write(6,*)
     1              ' Warning, clouds_3d > 1 - reset'
     1              ,i,j,k,clouds_3d
                    clouds_3d = 1.0
                endif
            else
                write(6,*)
     1          ' Warning, clouds_3d > 1 - reset for edge effect'
     1          ,i,j,k,clouds_3d
                clouds_3d = 1.0
            endif
        endif

        if(clouds_3d .lt. 0.0)then 
            write(6,*)' ERROR: clouds_3d < 0 - reset',i,j,k,clouds_3d       
            clouds_3d = 0.0
        endif

        return
        end
