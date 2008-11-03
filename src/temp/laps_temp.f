

        subroutine laps_temp(i4time_needed)

!       1997 Jun        Ken Dritz     Made NX_L, NY_L, NZ_L dummy arguments,
!                                     making non-dummy arrays dimensioned
!                                     therewith dynamic (automatic).
!       1997 Jun        Ken Dritz     Changed include to 
!                                     laps_static_parameters.inc.
!       1997 Jun        Ken Dritz     Added call to get_laps_cycle_time.
!       1999 Jan        Steve Albers  Added stability stuff

        use mem_namelist, ONLY: NX_L,NY_L,NZ_L=>nk_laps
     1                         ,laps_cycle_time,c6_maproj
     1                         ,grid_spacing_m_parm=>grid_spacing_m
     1                         ,iwrite_output

        use mem_grid, ONLY: lat,lon,topo
 
        use mem_temp, ONLY: temp_3d,heights_3d,pres_3d_pa

        include 'laps_static_parameters.inc'

        integer j_status(20),iprod_number(20)

!  ************ DECLARATIONS **************************************************

        real output_4d(NX_L,NY_L,NZ_L,2)

!       real temp_3d(NX_L,NY_L,NZ_L)
!       real heights_3d(NX_L,NY_L,NZ_L)
!       real pres_3d_pa(NX_L,NY_L,NZ_L)
        real pres_3d_mb(NX_L,NY_L,NZ_L)
        real temp_sfc_k(NX_L,NY_L)
        real pres_sfc_pa(NX_L,NY_L), pres_sfc_mb(NX_L,NY_L)
        real pres_msl_pa(NX_L,NY_L)
        real pbl_top_pa(NX_L,NY_L), pbl_top_mb(NX_L,NY_L)
        real pbl_depth_m(NX_L,NY_L)

        character*31 EXT

        character*10  units_2d
        character*125 comment_2d
        character*3 var_2d

        parameter (MAX_FIELDS=2)
        real field_array(NX_L,NY_L,MAX_FIELDS)
        character*125 comment_a(MAX_FIELDS)
        character*10 units_a(MAX_FIELDS)
        character*3 var_a(MAX_FIELDS)

!       Obtain grid spacing at the center
!       Test for a conformal map projection
        if(c6_maproj .ne. 'latlon' .and. c6_maproj .ne. 'icshdr')then 
            icen = NX_L/2 + 1
            jcen = NY_L/2 + 1
            call get_grid_spacing_actual(lat(icen,jcen),lon(icen,jcen)       
     1                                  ,grid_spacing_m,istatus)
            if(istatus .ne. 1)then
                write(6,*)' Error return from get_grid_spacing_actual'       
                return
            endif

        else
            write(6,*)' Non-conformal map projection: ',c6_maproj
            write(6,*)' Set grid_spacing_cen_m to parameter value'
            grid_spacing_m = grid_spacing_m_parm

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

!       Read in MSL pressure data
        var_2d = 'MSL'
        ext = 'lsx'
        call get_laps_2dgrid(i4time_needed,laps_cycle_time/2
     1                      ,i4time_nearest
     1                      ,ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                      ,pres_msl_pa,0,istatus)

        if(istatus .ne. 1)then
            write(6,*)' LAPS MSL Pres not available'
            write(6,*)' Not calling put_temp_anal'
            go to 999
        endif

!  ************ UPDATED ARGUMENT LIST ****************************************

        call put_temp_anal(i4time_needed
     1          ,NX_L,NY_L,NZ_L                  ! Input
     1          ,heights_3d                      ! Output
     1          ,lat,lon,topo                    ! Input
     1          ,temp_sfc_k                      ! Input
     1          ,pres_sfc_pa                     ! Input
     1          ,pres_msl_pa                     ! Input
     1          ,laps_cycle_time                 ! Input
     1          ,grid_spacing_m                  ! Input
     1          ,comment_2d                      ! Output
     1          ,temp_3d,pres_3d_pa,istatus)     ! Output

        if(iwrite_output .ge. 0 .and. istatus .eq. 1)then
            call write_temp_anal(i4time_needed,NX_L,NY_L,NZ_L,temp_3d       
     1                  ,heights_3d,comment_2d,istatus)
        endif

        I4_elapsed = ishow_timer()


!  ******************** PBL SECTION ******************************************

        if(.true. .and. istatus .eq. 1)then
            write(6,*)' Start PBL Section'
            pres_3d_mb  = pres_3d_pa / 100.
            pres_sfc_mb = pres_sfc_pa / 100.

            call ghbry (i4time_needed,pres_3d_mb,pres_sfc_mb          ! I
     1                 ,temp_sfc_k,temp_3d                            ! I
     1                 ,pbl_top_mb                                    ! O
     1                 ,NX_L,NY_L,NZ_L                                ! I
     1                 ,istatus)                                      ! O
            if(istatus .ne. 1)then
                write(6,*)' ERROR: on PBL istatus returned from ghbry'       
                return
            endif

            I4_elapsed = ishow_timer()

            pbl_top_pa = pbl_top_mb * 100.

!           Convert to PBL height AGL

!           Note that the 'pres_to_ht' call uses a linear interpolation that
!           can be upgraded (within the routine) to logp interpolation

!           The 'pressure_to_height' call uses log interpolation but will
!           become an unusable routine if we switch away from a constant
!           pressure vertical grid.

            do i = 1,NX_L
            do j = 1,NY_L
                call pres_to_ht(pbl_top_pa(i,j),pres_3d_pa,heights_3d
     1                         ,NX_L,NY_L,NZ_L,i,j,pbl_top_m,istatus)       

!               call pressure_to_height(pbl_top_pa(i,j),heights_3d
!    1                                 ,NX_L,NY_L,NZ_L,i,j
!    1                                 ,pbl_top_m,istatus)       

                pbl_depth_m(i,j) = max(pbl_top_m - topo(i,j),0.)
            enddo ! j
            enddo ! i

!           Write PBL file
            if(iwrite_output .ge. 0)then
                call move(pbl_top_pa ,field_array(1,1,1),NX_L,NY_L)
                call move(pbl_depth_m,field_array(1,1,2),NX_L,NY_L)

                ext = 'pbl'
                var_a(1) = 'PTP'
                var_a(2) = 'PDM'
                units_a(1) = 'Pa'
                units_a(2) = 'M'
                comment_a(1) = 'PBL Top Pressure'
                comment_a(2) = 'PBL Depth'
                call put_laps_multi_2d(i4time_needed,ext,var_a,units_a
     1                                ,comment_a,field_array,NX_L,NY_L
     1                                ,2,istatus)
            endif    

        else
            write(6,*)' No PBL calculation done for PBL file'

        endif

! ************* NOTIFICATION STUFF *********************************************

        I4_elapsed = ishow_timer()

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

