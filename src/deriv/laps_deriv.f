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
          
        call laps_deriv(i4time,
     1                  NX_L,NY_L,
     1                  NZ_L,
     1                  j_status)

999     continue

        end
          
        subroutine laps_deriv(i4time,
     1                  NX_L,NY_L,
     1                  NZ_L,
     1                  j_status)

        integer j_status(20),iprod_number(20)

        real*4 temp_3d(NX_L,NY_L,NZ_L)
        real*4 sh_3d_dum(NX_L,NY_L,NZ_L)
        real*4 heights_3d(NX_L,NY_L,NZ_L)
        real*4 temp_sfc_k(NX_L,NY_L)
        real*4 pres_sfc_pa(NX_L,NY_L)

        real*4 lat(NX_L,NY_L)
        real*4 lon(NX_L,NY_L)
        real*4 topo(NX_L,NY_L)

        character*31 EXT

        character*10  units_2d
        character*125 comment_2d
        character*3 var_2d

!       Get parameters for put_derived_wind_prods call
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

        write(6,*)
        write(6,*)' Calling put_derived_wind_prods'
        call put_derived_wind_prods(NX_L,NY_L,NZ_L           ! Input
     1          ,NX_L,NY_L,NZ_L                              ! Input (sic)
     1          ,max_radars_dum,r_missing_data               ! Input
     1          ,i4time)                                     ! Input

!       Read data for laps_deriv_sub and put_stability calls

!       Read lt1 - temp_3d,heights_3d
        var_2d = 'T3'
        ext = 'lt1'
        call get_laps_3d(i4time,NX_L,NY_L,NZ_L
     1      ,ext,var_2d,units_2d,comment_2d,temp_3d,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading 3D Temp'
            return
        endif

        var_2d = 'HT'
        ext = 'lt1'
        call get_laps_3d(i4time,NX_L,NY_L,NZ_L
     1      ,ext,var_2d,units_2d,comment_2d,heights_3d,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading 3D Heights'
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

        write(6,*)
        write(6,*)' Calling laps_deriv_sub'
        call laps_deriv_sub(i4time,              ! I
     1                  NX_L,NY_L,               ! I
     1                  NZ_L,                    ! I
     1                  N_PIREP,                 ! I
     1                  maxstns,                 ! I
     1                  max_cld_snd,             ! I
     1                  n_prods,
     1                  iprod_number,
     1                  temp_3d,                 ! I
     1                  heights_3d,              ! I
     1                  pres_sfc_pa,             ! I
     1                  temp_sfc_k,              ! I
     1                  j_status,                ! O
     1                  istatus1)                ! O

        if(.false. .and. istatus1 .eq. 1)then
            call get_domain_laps(NX_L,NY_L,LAPS_DOMAIN_FILE,lat,lon,topo       
     1                          ,grid_spacing_m,istatus)
            if(istatus .ne. 1)then
                write(6,*)' Error getting LAPS domain'
                go to 999
            endif

            call get_laps_cycle_time(laps_cycle_time,istatus)
            if (istatus .ne. 1) then
                write (6,*) 'Error getting LAPS cycle time'
                go to 999
            endif

            write(6,*)
            write(6,*)' Calling put_stability'
            call put_stability(
     1           i4time_needed                   ! I
     1          ,NX_L,NY_L,NZ_L                  ! I
     1          ,heights_3d                      ! I
     1          ,topo                            ! I
     1          ,laps_cycle_time                 ! I
     1          ,temp_3d                         ! I
     1          ,sh_3d_dum                       ! I
     1          ,temp_sfc_k                      ! I
     1          ,pres_sfc_pa                     ! I
     1          ,istatus)                        ! O
        else
            write(6,*)' put_stability not called for LST file'

        endif

 999    write(6,*)' End of subroutine laps_deriv'

        return

        end

