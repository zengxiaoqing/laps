  
        subroutine put_derived_wind_prods
     1    (imax,jmax,kmax                                         ! Input
     1    ,NX_L,NY_L,NZ_L                                         ! Input
     1    ,max_radars,r_missing_data                              ! Input
     1    ,i4time_sys)                                            ! Input

!       1997 Jun     Ken Dritz      Added NX_L, NY_L, NZ_L, and max_radars
!                                   as dummy arguments, making arrays
!                                   (including those in 'windparms.inc')
!                                   declared with those dimensions automatic.
!                                   Also added r_missing_data as dummy
!                                   argument.  (Resizability change)

        include 'windparms.inc'

!       Housekeeping
        integer*4 ss_normal,rtsys_no_data
        parameter (ss_normal        =1, ! success
     1             rtsys_no_data    =3) ! no data

        integer*4 j_status(20),i4time_array(20)
        character*3 exts(20)

!       Input Variables
        real*4 lat(imax,jmax),lon(imax,jmax),topo(imax,jmax)

        real*4 uanl(imax,jmax,kmax),vanl(imax,jmax,kmax)
     1        ,wanl_2d(imax,jmax)
        
        real*4 grid_ra_ref(imax,jmax,kmax)

!       Stuff for SFC and MEAN winds
        real*4 umean(NX_L,NY_L),vmean(NX_L,NY_L)
        real*4 ustorm(NX_L,NY_L),vstorm(NX_L,NY_L)
        real*4 out_lhe_3d(NX_L,NY_L,3)

!       Stuff for reading radar reflectivity
        character*4 radar_name
        character*255 c_filespec
ccc        character*4 radar_name(max_radars)
        character*31 ext_radar

!       Stuff for helicity
        real*4 helicity(NX_L,NY_L)

!       Dummy arrays
        real*4 dum1_2d(NX_L,NY_L)
        real*4 dum2_2d(NX_L,NY_L)
        real*4 dum3_2d(NX_L,NY_L)
        integer*4 idum1_2d(NX_L,NY_L)

!       Used for Steer Grid
!       real*4 iiilut(-NX_L:NX_L,-NY_L:NY_L)

!       Stuff for 2d fields
        real*4 lifted(NX_L,NY_L),liw(NX_L,NY_L)
        real*4 field_array(NX_L,NY_L,2)

!       Radar Stuff
        integer*4  n_fcst_radar
!       parameter (n_fcst_radar = 7200 / laps_cycle_time) ! Out to 2 hours
        parameter (n_fcst_radar = 0) ! No forecasts for now

        real*4 ref_max(NX_L,NY_L,0:n_fcst_radar)

        real*4 ref_curr_2d(NX_L,NY_L)
        real*4 ref_fcst_2d(NX_L,NY_L)


!       Dummy arrays.f steer_grid
!       real*4 dum4_2d(NX_L,NY_L)
!       real*4 dum5_2d(NX_L,NY_L)
!       real*4 dum6_2d(NX_L,NY_L)

        character*125 comment_2d,comment_a(0:10)
        character*10 units_2d,units_a(0:10)
        character*3 var_2d,var_a(0:10)

        character*50 directory
        character*31 ext

        write(6,*)
        write(6,*)' Entering Derived Wind Fields Subroutine',i4time_sys       
        write(6,*)' Version updated 3/29/96'

!       Housekeeping
        n_prods = 4
        do i = 1,n_prods
            j_status(i) = rtsys_no_data
            i4time_array(i) = i4time_sys ! Default Value
        enddo ! i

        istat_lhe = 0
        istat_lmr = 0
        istat_lf1 = 0
        istat_liw = 0

        n_lhe = 1
        n_liw = 2
        n_lmr = 3
        n_lf1 = 4
!       n_vrc = 5

!       exts(n_lwm) = 'lwm'
        exts(n_lhe) = 'lhe'
        exts(n_liw) = 'liw'
        exts(n_lmr) = 'lmr'
        exts(n_lf1) = 'lf1'
!       exts(n_vrc) = 'vrc'


        I4_elapsed = ishow_timer()

        call get_domain_laps(NX_L,NY_L,'nest7grid',lat,lon,topo
     1                                      ,grid_spacing_m,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error getting LAPS domain'
            return
        endif

!       Read in 3D U/V wind data
        ext = 'lw3'
        var_2d = 'U3'
        call get_laps_3d(i4time_sys,NX_L,NY_L,NZ_L
     1          ,ext,var_2d,units_2d,comment_2d,uanl,istat_lw3u)

        var_2d = 'V3'
        call get_laps_3d(i4time_sys,NX_L,NY_L,NZ_L
     1          ,ext,var_2d,units_2d,comment_2d,vanl,istat_lw3v)

        if(istat_lw3u .ne. 1 .or. istat_lw3v .ne. 1)then
            write(6,*)' Error reading in LW3 U/V fields - Abort'
            return
        endif


        call get_laps_cycle_time(ilaps_cycle_time,istatus)
        if(istatus .eq. 1)then
            write(6,*)' ilaps_cycle_time = ',ilaps_cycle_time
        else
            write(6,*)' Error getting laps_cycle_time'
            return
        endif

        write(6,*)
        write(6,*)' Reading reflectivity data'
        i4_tol = max(ilaps_cycle_time / 2, iradar_cycle_time / 2)

        if(.false.)then
!           call get_multiradar_ref(i4time_sys,i4_tol,i4time_radar_a
!     1    ,max_radars,n_radars,ext_radar,r_missing_data
!     1    ,.true.,NX_L,NY_L,NZ_L,lat,lon,topo,.true.,.false.
!     1    ,grid_ra_ref
!     1    ,n_ref_grids
!     1    ,rlat_radar,rlon_radar,rheight_radar,radar_name
!     1    ,istat_radar_2dref,istat_radar_3dref)

        else
            ext_radar = 'vrc'
            call get_directory(ext_radar,directory,len_dir)
            c_filespec = directory(1:len_dir)//'*.'//ext_radar
            call get_file_time(c_filespec,i4time_sys,i4time_radar)
            if(abs(i4time_radar - i4time_sys) .le. i4_tol)then

                call read_radar_2dref(i4time_radar,radar_name,
     1                  NX_L,NY_L,
     1                  ref_max(1,1,0),istat_radar_2dref)

            endif


        endif


        if(istat_radar_2dref .eq. 1)then
            write(6,*)' Radar 2d ref data successfully read in'
        else
            write(6,*)' Radar 2d ref data NOT successfully read in'
        endif


!       Calculate and write out Storm Steering Wind Field (ustorm, vstorm)
        write(6,*)
        write(6,*)' Calculating Mean Wind'

!       Get layer mean wind
        call mean_wind(uanl,vanl,topo,imax,jmax,kmax
     1        ,dum1_2d,dum2_2d,dum3_2d,idum1_2d ! Local arrays in mean_wind
     1        ,umean,vmean,ustorm,vstorm,istatus)

        if(istatus .ne. 1)then
            write(6,*)' FATAL ERROR in MEAN_WIND ROUTINE'
            return
        endif

        I4_elapsed = ishow_timer()

!       write(6,*)
!       write(6,*)' Calculating Storm Motion Grid'

!       Get storm motion from mean wind and tracking info
!       call steer_grid(i4time_lapswind,imax,jmax,kmax
!    1    ,dum1_2d,dum2_2d,dum3_2d,dum4_2d,grid_ra_veldum,grid_ra_ref,du
!    1m5_2d
!    1    ,dum6_2d
!    1    ,lat,lon,standard_latdum,standard_londum
!    1                  ,iiilut,umean,vmean
!    1                  ,ustorm,vstorm,istatus)


        I4_elapsed = ishow_timer()


        if(istat_radar_2dref.eq.1)then ! Get and advect the Reflectivity Stuff (LMR)
                                       ! Note: Ustorm and Vstorm need to be determined

            I4_elapsed = ishow_timer()

            do ifcst = 1,n_fcst_radar
                write(6,*)' Generating advected max reflectivity for pd 
     1',ifcst

                write(6,*)' Grid spacing (m) = ',grid_spacing_m

                call advect(ustorm,vstorm,ref_max(1,1,0)
     1                  ,dum1_2d,grid_spacing_m,imax,jmax
     1                       ,ref_max(1,1,ifcst)
     1   ,float(ilaps_cycle_time*ifcst),1.0,lon,r_missing_data)

            enddo ! i

            I4_elapsed = ishow_timer()

!           Move writing of LMR product to here with appropriate changes
!           Write out LMR file
            write(6,*)' Writing out LMR (or equiv) file for max reflecti
     1vity'
            do ifcst = 0,n_fcst_radar
                minutes_10 = (ilaps_cycle_time * ifcst) / 600
                if(minutes_10 .lt. 10)then
!                   write(var_a(ifcst),101)minutes_10
!101                format('R0',i1)
                    var_a(ifcst) = 'R'
                    write(comment_a(ifcst),111)minutes_10
111                 format('LAPS Max Reflectivity  ',i1,'0 Min Fcst')
                else
                    write(var_a(ifcst),102)minutes_10
102                 format('R',i2)
                    write(comment_a(ifcst),112)minutes_10
112                 format('LAPS Max Reflectivity ' ,i2,'0 Min Fcst')
                endif
                units_a(ifcst) = 'DBZ'
            enddo

            ext = 'lmr'

            call put_laps_multi_2d(i4time_sys,ext,var_a
     1      ,units_a,comment_a,ref_max,NX_L,NY_L,n_fcst_radar+1
     1                                                     ,istatus)

            if(istatus .ne. 1)then
                write(6,*)' Error writing out LMR file'
            endif

            istat_lmr = istatus

        else
            write(6,*)
     1       ' Not writing out Max Reflectivity LMR (or equiv) file'

        endif

!       Get and advect the Reflectivity History Stuff (LF1)
!       Note: Input grid_ra_ref is wiped out here
        write(6,*)
        write(6,*)' Generating max reflectivity history analysis'
        call get_radar_max_pd(i4time_sys-ilaps_cycle_time
     1          ,i4time_sys,imax,jmax,kmax
     1          ,lat,lon,topo,grid_ra_ref,dum1_2d
     1          ,ref_max(1,1,0),frac_sum,istat_radar_hist)

        if(istat_radar_hist.eq.1)then

            I4_elapsed = ishow_timer()

            do ifcst = 1,n_fcst_radar
                write(6,*)
     1     ' Generating advected max reflectivity history for pd',ifcst

                call advect(ustorm,vstorm,ref_max(1,1,0)
     1                  ,dum1_2d,grid_spacing_m,imax,jmax
     1                       ,ref_max(1,1,ifcst)
     1         ,float(ilaps_cycle_time*ifcst),1.0,lon,r_missing_data)

            enddo ! i

            I4_elapsed = ishow_timer()

!           Write out LF1 file
            write(6,*)' Writing out LF1 (or equiv) file for max reflecti
     1vity history'

            do ifcst = 0,n_fcst_radar
                minutes_10 = (ilaps_cycle_time * ifcst) / 600
                if(minutes_10 .lt. 10)then
!                   write(var_a(ifcst),201)minutes_10
!201                format('H0',i1)
                    var_a(ifcst) = 'H'
                    write(comment_a(ifcst),211)minutes_10
211                 format('LAPS Max Reflectivity History  ',i1,'0 Min F
     1cst')
                else
                    write(var_a(ifcst),202)minutes_10
202                 format('H',i2)
                    write(comment_a(ifcst),212)minutes_10
212                 format('LAPS Max Reflectivity History ' ,i2,'0 Min F
     1cst')
                endif
                units_a(ifcst) = 'DBZ'
            enddo

            ext = 'lf1'

            call put_laps_multi_2d(i4time_sys,ext,var_a,units_a
     1          ,comment_a,ref_max,NX_L,NY_L,n_fcst_radar+1,istatus)

            if(istatus .ne. 1)then
                write(6,*)' Error writing out LF1 file'
            endif

            istat_lf1 = istatus

        else
            write(6,*)' Not writing out Max Reflectivity History LF1 (or
     1 equiv) file'
        endif


        if(.true.)then ! Calculate Helicity
                       ! Note: (u+vstorm must be calculated first)

!           Write out Helicity field
            write(6,*)' Calculating Helicity'
            call helicity_laps(uanl,vanl,ustorm,vstorm,topo
     1  ,dum1_2d,dum2_2d,dum3_2d,idum1_2d ! Local arrays in helicity
     1  ,imax,jmax,kmax,helicity,istatus)

            if(istatus .eq. 1)then

                ext = 'lhe'

                var_a(0) = 'LHE'
                var_a(1) = 'MU'
                var_a(2) = 'MV'

                units_a(0) = 'm/s**2'
                units_a(1) = 'm/s'
                units_a(2) = 'm/s'

                comment_a(0) = 'LAPS Helicity'
                comment_a(1) = 'LAPS Mean Wind SFC - 300mb'
                comment_a(2) = 'LAPS Mean Wind SFC - 300mb'

                call move(helicity,out_lhe_3d(1,1,1),NX_L,NY_L)
                call move(ustorm  ,out_lhe_3d(1,1,2),NX_L,NY_L)
                call move(vstorm  ,out_lhe_3d(1,1,3),NX_L,NY_L)

                call put_laps_multi_2d(i4time_sys,ext,var_a
     1          ,units_a,comment_a,out_lhe_3d,NX_L,NY_L,3,istatus)
                 istat_lhe = istatus

            endif

        endif

        I4_elapsed = ishow_timer()

        if(.true.)then ! LIW

!           Calculate Li * 600mb Omega
            write(6,*)' Generating Li * Omega file'

!           Read in LI data
            var_2d = 'LI'
            ext = 'lsx'
            call get_laps_2dgrid(i4time_sys,0,i4time_nearest,
     1          ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                          ,lifted,0,istatus)

!           Read in omega data

!           Determine k coordinate for 600mb for passing in omega (wanl array)
            lvl_liw = nint(zcoord_of_pressure(60000.))
            ipres = nint( pressure_of_level(lvl_liw) / 100. )
            write(6,*)' LVL for Lifted index * omega = ',lvl_liw,ipres

            var_2d = 'OM'
            ext = 'lw3'
            call get_laps_2dgrid(i4time_sys,0,i4time_nearest,
     1          ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                ,wanl_2d,ipres,istat_omega)

            if(istatus .ne. 1 .or. istat_omega .ne. 1)then
                write(6,*)' Error reading Lifted Index / OM data'

            else ! Good Li Data
                call cpt_liw(lifted,wanl_2d,imax,jmax,liw)

                call move(liw    ,field_array(1,1,1),NX_L,NY_L)
                call move(wanl_2d,field_array(1,1,2),NX_L,NY_L)

!               Write out LIW field
!               Note that these arrays start off with 0 as the first index
                var_a(0) = 'LIW'
                var_a(1) = 'W'
                ext = 'liw'
!               call get_directory(ext,directory,len_dir)
                units_a(0) = 'K-Pa/s'
                units_a(1) = 'Pa/s  '
                comment_a(0) = 'Log LAPS Li * 600mb Omega'
                comment_a(1) = 'LAPS 600mb Omega'

                call put_laps_multi_2d(i4time_sys,ext,var_a,units_a
     1                    ,comment_a,field_array,imax,jmax,2,istatus)

                istat_liw = istatus

            endif

        endif

!       Set notification arrays.f derived wind products
!       if(istat_radar_2dref .eq. 1)then ! Set notification arrays for VRC file
!           i4time_array(n_vrc) = i4time_radar
!           j_status(n_vrc) = ss_normal
!       endif


        if(istat_lhe .eq. 1)then
            i4time_array(n_lhe) = i4time_sys
            j_status(n_lhe) = ss_normal
        endif

        if(istat_lmr .eq. 1)then
            i4time_array(n_lmr) = i4time_sys
            j_status(n_lmr) = ss_normal
        endif

        if(istat_lf1 .eq. 1)then
            i4time_array(n_lf1) = i4time_sys
            j_status(n_lf1) = ss_normal
        endif

        if(istat_liw .eq. 1)then
            i4time_array(n_liw) = i4time_sys
            j_status(n_liw) = ss_normal
        endif


        write(6,*)' Status of output products'
        do i = 1,n_prods
            write(6,*)i,' ',exts(i),' '
     1         ,i4time_array(i),j_status(i)
        enddo ! i

        return
        end
