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

!       1997 Jun        Ken Dritz     Added calls to get_grid_dim_xy,
!                                     get_laps_dimensions to get values of
!                                     NX_L, NY_L, NZ_L.
!       1997 Jun        Ken Dritz     Now pass NX_L, NY_L, NZ_L to laps_temp.
 
        integer*4 j_status(20),iprod_number(20)
        character*9 a9_time

        call get_systime(i4time,a9_time,istatus)
        if(istatus .ne. 1)go to 999

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

        if(.true. .and. istatus .eq. 1)then
            call put_stability(
     1           i4time_needed                   ! Input
     1          ,NX_L,NY_L,NZ_L                  ! Input
     1          ,heights_3d                      ! Input
     1          ,topo                            ! Input
     1          ,laps_cycle_time                 ! Input
     1          ,temp_3d                         ! Input
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


        subroutine put_stability(
     1           i4time_needed                   ! Input
     1          ,NX_L,NY_L,NZ_L                  ! Input
     1          ,heights_3d                      ! Input
     1          ,topo                            ! Input
     1          ,laps_cycle_time                 ! Input
     1          ,temp_3d                         ! Input
     1          ,temp_sfc_k                      ! Input
     1          ,pres_sfc_pa                     ! Input
     1          ,istatus)                        ! Output

!       Arrays passed in
        real*4 temp_3d(NX_L,NY_L,NZ_L)
        real*4 heights_3d(NX_L,NY_L,NZ_L)
        real*4 temp_sfc_k(NX_L,NY_L)
        real*4 pres_sfc_pa(NX_L,NY_L)
        real*4 topo(NX_L,NY_L)

!       Local declarations for stability 
        real*4 t_sfc_f(NX_L,NY_L)
        real*4 td_sfc_k(NX_L,NY_L)
        real*4 td_sfc_f(NX_L,NY_L)
        real*4 pbe_2d(NX_L,NY_L)
        real*4 nbe_2d(NX_L,NY_L)
        real*4 pres_sfc_mb(NX_L,NY_L)
        real*4 li(NX_L,NY_L)
        real*4 t500laps(NX_L,NY_L)
        real*4 p_1d_pa(NZ_L) ! This should eventually be pres_3d
        character*10 units_2d_a(3)
        character*125 comment_2d_a(3)
        character*3 var_2d_a(3)
        real*4 out_multi_2d(NX_L,NY_L,3)

        character*31 EXT

        character*3 var_2d
        character*10  units_2d
        character*125 comment_2d

        real*4 k_to_f

!       Read in surface dewpoint data
        var_2d = 'TD'
        ext = 'lsx'
        call get_laps_2dgrid(i4time_needed,laps_cycle_time/2
     1                      ,i4time_nearest
     1                      ,ext,var_2d,units_2d,comment_2d,NX_L,NY_L       
     1                      ,td_sfc_k,0,istatus)
        if(istatus .ne. 1)then
            write(6,*)' LAPS Sfc Dewpoint not available'
            write(6,*)' Abort put_stability routine'
            return
        endif

!       Fill p_1d_pa, later we can change this to pres_3d?
        do k = 1,NZ_L
            p_1d_pa(k) = pressure_of_level(k)
        enddo ! k

        call laps_be(NX_L,NY_L,NZ_L
     1              ,temp_sfc_k,td_sfc_k,pres_sfc_pa
     1              ,temp_3d,heights_3d,p_1d_pa,topo   
     1              ,pbe_2d,nbe_2d)

!       Fill pres_sfc_mb
        call move(pres_sfc_pa,pres_sfc_mb,NX_L,NY_L)
        call multcon(pres_sfc_mb,.01,NX_L,NY_L)

!       Convert T, Td to F
        do i = 1,NX_L
        do j = 1,NY_L
            t_sfc_f(i,j)  = k_to_f(temp_sfc_k(i,j))
            td_sfc_f(i,j) = k_to_f(td_sfc_k(i,j))
        enddo ! j
        enddo ! i

        flag = 0.0
        call li_laps(t_sfc_f,td_sfc_f,pres_sfc_mb,t500laps ! Local?
     1              ,i4time_needed,NX_L,NY_L,li,flag,istatus)

!       call move
        call move(pbe_2d,out_multi_2d(1,1,1),NX_L,NY_L)
        call move(nbe_2d,out_multi_2d(1,1,2),NX_L,NY_L)
        call move(li,out_multi_2d(1,1,3),NX_L,NY_L)

!       add var arrays
        ext = 'lst'

        var_2d_a(1) = 'PBE'
        var_2d_a(2) = 'NBE'
        var_2d_a(3) = 'LI'

        units_2d_a(1) = 'J/KG'
        units_2d_a(2) = 'J/KG'
        units_2d_a(3) = 'K'

        call put_laps_multi_2d(i4time_needed,ext,var_2d_a,units_2d_a
     1                        ,comment_2d_a,out_multi_2d,NX_L,NY_L,3    
     1                        ,istatus)
        if(istatus .ne. 1)then
            write(6,*)' LST output error'
        else
            write(6,*)' Successfully wrote LST'
        endif

        return

        end

