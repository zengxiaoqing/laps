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

!       1997 Jun        Ken Dritz     Added calls to get_grid_dim_xy,
!                                     get_laps_dimensions to get values of
!                                     NX_L, NY_L, NZ_L.
!       1997 Jun        Ken Dritz     Now pass NX_L, NY_L, NZ_L to laps_temp.
 
        integer*4 j_status(20),iprod_number(20)
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
 
        call laps_temp      (i4time,
     1                       NX_L,NY_L,NZ_L,
     1                       i_diag,
     1                       n_prods,
     1                       iprod_number,
     1                       j_status)

999     continue

        end


        subroutine laps_temp(i4time_needed,
     1                       NX_L,NY_L,NZ_L,
     1                       i_diag,
     1                       n_prods,
     1                       iprod_number,
     1                       j_status)

!       1997 Jun        Ken Dritz     Made NX_L, NY_L, NZ_L dummy arguments,
!                                     making non-dummy arrays dimensioned
!                                     therewith dynamic (automatic).
!       1997 Jun        Ken Dritz     Changed include to 
!                                     laps_static_parameters.inc.
!       1997 Jun        Ken Dritz     Added call to get_laps_cycle_time.
!       1999 Jan        Steve Albers  Added stability stuff

        include 'laps_static_parameters.inc'

        integer*4 j_status(20),iprod_number(20)


!  ************ DECLARATIONS **************************************************

        real*4 output_4d(NX_L,NY_L,NZ_L,2)

        integer*4 iflag_write

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

c read in LAPS_DOMAIN
        call get_domain_laps(NX_L,NY_L,LAPS_DOMAIN_FILE,lat,lon,topo
     1               ,grid_spacing_m,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error getting LAPS domain'
            go to 999
        endif
        call get_laps_cycle_time(laps_cycle_time,istatus)
        if (istatus .ne. 1) then
            write (6,*) 'Error getting LAPS cycle time'
            go to 999
        endif

!       Read in surface temp data
        var_2d = 'T'
        ext = 'lsx'
        call get_laps_2dgrid(i4time_needed,laps_cycle_time/2
     1                      ,i4time_nearest
     1                      ,ext,var_2d,units_2d,comment_2d,NX_L,NY_L       
     1                      ,temp_sfc_k,0,istatus)

        if(istatus .ne. 1)then
            write(6,*)' LAPS Sfc Temp not available'
            write(6,*)' Not calling put_temp_anal'
            go to 999
        endif

!       Read in surface pressure data
        var_2d = 'PS'
        ext = 'lsx'
        call get_laps_2dgrid(i4time_needed,laps_cycle_time/2
     1                      ,i4time_nearest
     1                      ,ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                      ,pres_sfc_pa,0,istatus)

        if(istatus .ne. 1)then
            write(6,*)' LAPS Sfc Pres not available'
            write(6,*)' Not calling put_temp_anal'
            go to 999
        endif

!  ************ UPDATED ARGUMENT LIST ****************************************

        iflag_write = 1 ! Flag on whether to write out LT1 (0 or 1)

        call put_temp_anal(i4time_needed
     1          ,NX_L,NY_L,NZ_L                  ! Input
     1          ,heights_3d                      ! Output
     1          ,lat,lon,topo                    ! Input
     1          ,temp_sfc_k                      ! Input/Output
     1          ,pres_sfc_pa                     ! Input
     1          ,iflag_write                     ! Input
     1          ,laps_cycle_time                 ! Input
     1          ,grid_spacing_m                  ! Input
     1          ,temp_3d,istatus)                ! Output

!  ************ STABILITY TEST STUFF ******************************************

!       This code is here just in case we need it due to a change in stability
!       index strategy.

        if(.false. .and. istatus .eq. 1)then
            call put_stability(
     1           i4time_needed                   ! Input
     1          ,NX_L,NY_L,NZ_L                  ! Input
     1          ,heights_3d                      ! Input
     1          ,topo                            ! Input
     1          ,laps_cycle_time                 ! Input
     1          ,temp_3d                         ! Input
     1          ,sh_3d_dum                       ! Input
     1          ,temp_sfc_k                      ! Input
     1          ,pres_sfc_pa                     ! Input
     1          ,istatus)                        ! Output
        else
            write(6,*)' No stability calculation done for LST file'

        endif

! ************* NOTIFICATION STUFF *********************************************

        iprod_number(1) = 28261 ! LT1
        n_prods = 1

999     if(istatus .eq. 1)then
            j_status(1) = 1 ! Normal product
        else
            j_status(1) = 4 ! No product
        endif

! ****************************************************************************

        return
        end

