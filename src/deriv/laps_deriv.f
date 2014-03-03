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
        program laps_deriv_main

        integer j_status(20)

        character*9 a9_time

        call get_systime(i4time,a9_time,istatus)
        if(istatus .ne. 1)go to 999

        write(6,*)' systime = ',a9_time

        call get_grid_dim_xy(NX_L,NY_L,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'Error getting horizontal domain dimensions'
           go to 999
        endif

        call get_laps_dimensions(NZ_L,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'Error getting vertical domain dimension'
           go to 999
        endif

        call get_r_missing_data(r_missing_data,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'Error getting r_missing_data'
           go to 999
        endif
          
        call laps_deriv(i4time,
     1                  NX_L,NY_L,
     1                  NZ_L,
     1                  r_missing_data,
     1                  j_status)

999     continue

        end
          
        subroutine laps_deriv(i4time,
     1                  NX_L,NY_L,
     1                  NZ_L,
     1                  r_missing_data,
     1                  j_status)

        use mem_namelist, ONLY: max_snd_grid, max_snd_levels

        integer j_status(20),iprod_number(20)

        real temp_3d(NX_L,NY_L,NZ_L)
        real rh_3d_pct(NX_L,NY_L,NZ_L)
        real td_3d_k(NX_L,NY_L,NZ_L)
        real heights_3d(NX_L,NY_L,NZ_L)
        real u_3d(NX_L,NY_L,NZ_L)
        real v_3d(NX_L,NY_L,NZ_L)

        real temp_sfc_k(NX_L,NY_L)
        real pres_sfc_pa(NX_L,NY_L)
        real rh_sfc_pct(NX_L,NY_L)
        real tpw_2d(NX_L,NY_L)         ! units are M
        real u_sfc_ms(NX_L,NY_L)
        real v_sfc_ms(NX_L,NY_L)

        real dbz_max_2d(NX_L,NY_L)

        real lat(NX_L,NY_L)
        real lon(NX_L,NY_L)
        real topo(NX_L,NY_L)
        real ldf(NX_L,NY_L)

        character*31 EXT

        character*10  units_2d
        character*125 comment_2d
        character*3 var_2d

!       Get parameters for laps_deriv_sub call
        call get_meso_sao_pirep(N_MESO,N_SAO,N_PIREP,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'Error getting N_PIREP'
           go to 999
        endif

        call get_maxstns(maxstns,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'Error getting maxstns'
           go to 999
        endif

        max_cld_snd = maxstns + N_PIREP

!       Read data for laps_deriv_sub and put_stability calls

!       Read LT1 - temp_3d
        var_2d = 'T3'
        ext = 'lt1'
        call get_laps_3d(i4time,NX_L,NY_L,NZ_L
     1      ,ext,var_2d,units_2d,comment_2d,temp_3d,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading 3D Temp'
            return
        endif
        call qc_field_3d('T3',temp_3d,NX_L,NY_L,NZ_L,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading 3D Temp'
            return
        endif

!       Read LT1 - heights_3d
        var_2d = 'HT'
        ext = 'lt1'
        call get_laps_3d(i4time,NX_L,NY_L,NZ_L
     1      ,ext,var_2d,units_2d,comment_2d,heights_3d,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading 3D Heights'
            return
        endif
        call qc_field_3d('HT',heights_3d,NX_L,NY_L,NZ_L,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading 3D Heights'
            return
        endif

!       Read RH
        var_2d = 'RHL'
        ext = 'lh3'
        call get_laps_3d(i4time,NX_L,NY_L,NZ_L
     1      ,ext,var_2d,units_2d,comment_2d,rh_3d_pct,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading 3D RH (lh3/RHL)'
            return
        endif
        call qc_field_3d('RHL',rh_3d_pct,NX_L,NY_L,NZ_L,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading 3D RH (lh3/RHL)'
            return
        endif

!       Read in surface temp data
        var_2d = 'T'
        ext = 'lsx'
        call get_laps_2d(i4time,ext,var_2d,units_2d,comment_2d
     1                  ,NX_L,NY_L,temp_sfc_k,istatus)

        if(istatus .ne. 1)then
            write(6,*)' LAPS Sfc Temp not available'
            go to 999
        endif

!       Read in surface pressure data
        var_2d = 'PS'
        ext = 'lsx'
        call get_laps_2d(i4time,ext,var_2d,units_2d,comment_2d
     1                  ,NX_L,NY_L,pres_sfc_pa,istatus)

        if(istatus .ne. 1)then
            write(6,*)' LAPS Sfc Pres not available'
            go to 999
        endif

!       Read in surface rh data
        var_2d = 'RH'
        ext = 'lsx'
        call get_laps_2d(i4time,ext,var_2d,units_2d,comment_2d
     1                  ,NX_L,NY_L,rh_sfc_pct,istatus)

        if(istatus .ne. 1)then
            write(6,*)' LAPS Sfc RH not available'
            go to 999
        endif

!       Read in surface u component
        var_2d = 'U'
        ext = 'lsx'
        call get_laps_2d(i4time,ext,var_2d,units_2d,comment_2d
     1                  ,NX_L,NY_L,u_sfc_ms,istatus)

        if(istatus .ne. 1)then
            write(6,*)' LAPS Sfc U not available'
            go to 999
        endif

!       Read in surface v component
        var_2d = 'V'
        ext = 'lsx'
        call get_laps_2d(i4time,ext,var_2d,units_2d,comment_2d
     1                  ,NX_L,NY_L,v_sfc_ms,istatus)

        if(istatus .ne. 1)then
            write(6,*)' LAPS Sfc V not available'
            go to 999
        endif

!       Read in tpw data
        var_2d = 'TPW'
        ext = 'lh4'
        call get_laps_2d(i4time,ext,var_2d,units_2d,comment_2d
     1                  ,NX_L,NY_L,tpw_2d,istatus)

        if(istatus .ne. 1)then
            write(6,*)' LAPS TPW not available'
            go to 999
        endif

        write(6,*)
        write(6,*)' Calling laps_deriv_sub'
        call laps_deriv_sub(i4time,                  ! I
     1                  NX_L,NY_L,                   ! I
     1                  NZ_L,                        ! I
     1                  N_PIREP,                     ! I
     1                  maxstns,                     ! I
     1                  max_snd_grid,max_snd_levels, ! I
     1                  n_prods,
     1                  iprod_number,
     1                  temp_3d,                     ! I
     1                  heights_3d,                  ! I
     1                  rh_3d_pct,                   ! I
     1                  pres_sfc_pa,                 ! I
     1                  temp_sfc_k,                  ! I
     1                  dbz_max_2d,istat_lps,        ! O
     1                  twet_snow,                   ! O
     1                  j_status,                    ! O
     1                  istatus1)                    ! O

!       istat_lps = 0 ! Temporary until this is wired in'
        istat_twet_snow = 1 ! if we get this far the read was successful

        call get_laps_domain_95(NX_L,NY_L,lat,lon,topo
     1                         ,ldf,grid_spacing_m,istatus)       
        if(istatus .ne. 1)then
            write(6,*)' Error getting LAPS domain'
            go to 999
        endif

        if(istat_twet_snow .eq. 1)then
            call get_laps_cycle_time(laps_cycle_time,istatus)
            if (istatus .ne. 1) then
                write (6,*) 'Error getting LAPS cycle time'
                go to 999
            endif

            I4_elapsed = ishow_timer()

            write(6,*)
            write(6,*)' Calling put_stability'
            call put_stability(
     1           i4time                          ! I
     1          ,NX_L,NY_L,NZ_L                  ! I
     1          ,heights_3d                      ! I
     1          ,lat,lon,topo                    ! I
     1          ,laps_cycle_time                 ! I
     1          ,temp_3d                         ! I
     1          ,rh_3d_pct                       ! I
     1          ,temp_sfc_k                      ! I
     1          ,pres_sfc_pa                     ! I
     1          ,twet_snow                       ! I
     1          ,td_3d_k                         ! O
     1          ,istat_lst)                      ! O
        else
            write(6,*)' put_stability not called for LST file'

        endif

        I4_elapsed = ishow_timer()

!       If we need space we can deallocate rh_3d_pct here
!       If we need space we can allocate u_3d, v_3d here

        write(6,*)
        write(6,*)' Calling put_derived_wind_prods'
        call put_derived_wind_prods(NX_L,NY_L,NZ_L           ! Input
     1          ,NX_L,NY_L,NZ_L                              ! Input (sic)
     1          ,max_radars_dum,r_missing_data               ! Input
     1          ,i4time                                      ! Input
     1          ,dbz_max_2d,istat_lps                        ! Input
     1          ,lat,lon,topo,ldf                            ! Input
     1          ,tpw_2d                                      ! Input
     1          ,heights_3d                                  ! Input
     1          ,u_3d,v_3d)                                  ! Output

        write(6,*)
        if(istat_lst .eq. 1)then

            write(6,*)' Calling fire_fields'
            call fire_fields(NX_L,NY_L,NZ_L,temp_3d,td_3d_k          ! I
     1                      ,u_3d,v_3d                               ! I
     1                      ,temp_sfc_k,pres_sfc_pa                  ! I
     1                      ,rh_sfc_pct,u_sfc_ms,v_sfc_ms            ! I
     1                      ,r_missing_data,i4time                   ! I
     1                      ,istatus)                                ! O

        else
            write(6,*)' Skipping call to fire_fields'

        endif

 999    write(6,*)' End of subroutine laps_deriv'

        return

        end

 
       subroutine get_deriv_parms(mode_evap,l_bogus_radar_w,              ! O
     1                            l_deep_vv,                              ! O
     1                            vv_to_height_ratio_Cu,                  ! O
     1                            vv_to_height_ratio_Sc,                  ! O
     1                            vv_for_St,                              ! O
     1                            c_z2m,                                  ! O
     1                            thresh_cvr_cty_vv,thresh_cvr_lwc,       ! O
     1                            twet_snow,                              ! O
     1                            hydrometeor_scale_cldliq,               ! O
     1                            hydrometeor_scale_cldice,               ! O
     1                            hydrometeor_scale_pcp,                  ! O
     1                            istatus)                                ! O

       Real vv_to_height_ratio_Cu
       Real vv_to_height_ratio_Sc
       Real vv_for_St
       Real thresh_cvr_cty_vv,thresh_cvr_lwc,twet_snow
       Real hydrometeor_scale_cldliq,hydrometeor_scale_cldice
       Real hydrometeor_scale_pcp

       character*20 c_z2m

       logical l_bogus_radar_w, l_deep_vv

       namelist /deriv_nl/ mode_evap, l_bogus_radar_w, l_deep_vv,
     1                     vv_to_height_ratio_Cu,
     1                     vv_to_height_ratio_Sc,
     1                     vv_for_St,
     1                     c_z2m,
     1                     thresh_cvr_cty_vv,thresh_cvr_lwc,
     1                     twet_snow,
     1                     hydrometeor_scale_cldliq,
     1                     hydrometeor_scale_cldice,
     1                     hydrometeor_scale_pcp
 
       character*150 static_dir,filename

!      Set default value(s)
       c_z2m = 'albers'
       hydrometeor_scale_cldliq = 0.2
       hydrometeor_scale_cldice = 0.2
       hydrometeor_scale_pcp    = 1.0
 
       call get_directory('static',static_dir,len_dir)

       filename = static_dir(1:len_dir)//'/deriv.nl'
 
       open(1,file=filename,status='old',err=900)
       read(1,deriv_nl,err=901)
       close(1)

       print*,'success reading deriv_nl in ',filename
       write(*,deriv_nl)

       istatus = 1
       return

  900  print*,'error opening file ',filename
       istatus = 0
       return

  901  print*,'error reading deriv_nl in ',filename
       write(*,deriv_nl)
       istatus = 0
       return

       end
