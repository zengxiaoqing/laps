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
  
        subroutine put_derived_wind_prods
     1    (imax,jmax,kmax                                         ! Input
     1    ,NX_L,NY_L,NZ_L                                         ! Input
     1    ,max_radars,r_missing_data                              ! Input
     1    ,i4time_sys                                             ! Input
     1    ,dbz_max_2d,istat_lps                                   ! Input
     1    ,lat,lon,topo,ldf                                       ! Input
     1    ,tpw_2d                                                 ! Input
     1    ,heights_3d                                             ! Input
     1    ,uanl,vanl)                                             ! Output

!       1997 Jun     Ken Dritz      Added NX_L, NY_L, NZ_L, and max_radars
!                                   as dummy arguments, making arrays
!                                   (including those in 'windparms.inc')
!                                   declared with those dimensions automatic.
!                                   Also added r_missing_data as dummy
!                                   argument.  (Resizability change)

        include 'windparms.inc'

!       Housekeeping
        integer ss_normal,rtsys_no_data
        parameter (ss_normal        =1, ! success
     1             rtsys_no_data    =3) ! no data

        integer j_status(20),i4time_array(20)
        character*3 exts(20)

!       Array Variables
        real lat(imax,jmax),lon(imax,jmax),topo(imax,jmax)
        real ldf(imax,jmax)      
        real dx(NX_L,NY_L)
        real dy(NX_L,NY_L)

        real uanl(imax,jmax,kmax),vanl(imax,jmax,kmax)
     1        ,wanl_2d(imax,jmax)
        
        real grid_ra_ref(imax,jmax,kmax)
        real dbz_max_2d(imax,jmax)
        real heights_3d(imax,jmax,kmax)

!       Stuff for SFC and MEAN winds
        real umean(NX_L,NY_L),vmean(NX_L,NY_L)
        real ushear(NX_L,NY_L),vshear(NX_L,NY_L)
        real ustorm(NX_L,NY_L),vstorm(NX_L,NY_L)
        real out_lhe_3d(NX_L,NY_L,5)

!       Stuff for reading radar reflectivity
        character*4 radar_name
        character*255 c_filespec
        character*31 ext_radar

!       Stuff for helicity
        real helicity(NX_L,NY_L)

        real tpw_2d(NX_L,NY_L)         ! units are M
        real upslope_flux(NX_L,NY_L)         

!       Dummy arrays
        real dum1_2d(NX_L,NY_L)
        real dum2_2d(NX_L,NY_L)
        real dum3_2d(NX_L,NY_L)
        integer idum1_2d(NX_L,NY_L)

!       Used for Steer Grid
!       real iiilut(-NX_L:NX_L,-NY_L:NY_L)

!       Stuff for 2d fields
        real lifted(NX_L,NY_L),liw(NX_L,NY_L)
        real field_array(NX_L,NY_L,2)

!       Radar Stuff
        integer  n_fcst_radar
!       parameter (n_fcst_radar = 7200 / laps_cycle_time) ! Out to 2 hours
        parameter (n_fcst_radar = 0) ! No forecasts for now

        real ref_max(NX_L,NY_L,0:n_fcst_radar)

        real ref_curr_2d(NX_L,NY_L)
        real ref_fcst_2d(NX_L,NY_L)

        character*125 comment_2d,comment_a(0:10)
        character*10 units_2d,units_a(0:10)
        character*3 var_2d,var_a(0:10)

        character*31 ext

        character*40 c_vars_req
        character*180 c_values_req

        write(6,*)
        write(6,*)' Entering Derived Wind Fields Subroutine',i4time_sys       

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

        c_vars_req = 'radarext_3d'

        call get_static_info(c_vars_req,c_values_req,1,istatus)

        if(istatus .eq. 1)then
            ext_radar = c_values_req(1:3)
        else
            write(6,*)' Error getting radarext_3d'
            return
        endif

        write(6,*)' ext_radar = ',ext_radar

!       i4_tol = max(ilaps_cycle_time / 2, iradar_cycle_time / 2)

        call get_filespec(ext_radar,2,c_filespec,istatus)

        write(6,*)' istat_lps = ',istat_lps

        if(istat_lps .eq. 1)then
            write(6,*)' passing in "dbz_max_2d" for "ref_max" array'
            ref_max(:,:,0) = dbz_max_2d
            istat_radar_2dref = 1 ! since we can now use the data

        elseif(ext_radar .ne. 'vrc')then  ! original way to get column max ref
            call get_ref_base(ref_base,istatus)
            if(istatus .ne. 1)then
                write(6,*)' Error getting ref_base'
                return
            endif

            i4_tol = 1200

            call read_radar_3dref(i4time_sys,
     1                 i4_tol,i4_ret,                                   ! I/O
     1                 .true.,ref_base,
     1                 NX_L,NY_L,NZ_L,ext_radar,
     1                 lat,lon,topo,.true.,.true.,
     1                 heights_3d,
     1                 grid_ra_ref,
     1                 rlat_radar,rlon_radar,rheight_radar,radar_name,     
     1                 n_ref_grids,istat_radar_2dref,istat_radar_3dref)       

            if(istat_radar_2dref .eq.  1 .or. 
     1         istat_radar_2dref .eq. -1          )then
                write(6,*)' Call get_max_reflect'

                call get_max_reflect(grid_ra_ref,NX_L,NY_L,NZ_L
     1                              ,ref_base,ref_max(1,1,0))

                istat_radar_2dref = 1 ! since we can now use the data

            endif

        else ! ext_radar = 'vrc' (and not lps)
            call read_radar_2dref(i4time_sys,radar_name,
     1                  NX_L,NY_L,
     1                  ref_max(1,1,0),istat_radar_2dref)

            if(istat_radar_2dref .eq. 1)then
                write(6,*)' Radar 2d ref data successfully read in'
            elseif(istat_radar_2dref .eq. -1)then
                write(6,*)' Radar 2d ref: fill missing data'
                do i = 1,NX_L
                do j = 1,NY_L
                    if(ref_max(i,j,0) .eq. r_missing_data)then
                        ref_max(i,j,0) = ref_base
                    endif
                enddo ! j
                enddo ! i
                istat_radar_2dref = 1 ! since we can now use the data
            else
                write(6,*)' Radar 2d ref data NOT successfully'
     1                   ,' read in'
            endif

        endif ! ext_radar

!       Calculate and write out Storm Steering Wind Field (ustorm, vstorm)

!       Get layer mean wind
        if(.false.)then
            call mean_wind(uanl,vanl,topo,imax,jmax,kmax
     1                    ,umean,vmean,ustorm,vstorm,istatus)
        else
            call mean_wind_bunkers(uanl,vanl,topo,imax,jmax,kmax       ! I
     1                    ,heights_3d                                  ! I
     1                    ,umean,vmean,ushear,vshear,ustorm,vstorm     ! O
     1                    ,istatus)                                    ! O
        endif

        if(istatus .ne. 1)then
            write(6,*)' FATAL ERROR in MEAN_WIND ROUTINE'
            return
        endif

        I4_elapsed = ishow_timer()

!       write(6,*)
!       write(6,*)' Calculating Storm Motion Grid'

!       Get storm motion from mean wind and tracking info
!       call steer_grid(i4time_lapswind,imax,jmax,kmax
!    1    ,dum1_2d,dum2_2d,dum3_2d,dum4_2d,grid_ra_veldum,grid_ra_ref
!    1    ,dum5_2d,dum6_2d
!    1    ,lat,lon,standard_latdum,standard_londum
!    1                  ,iiilut,umean,vmean
!    1                  ,ustorm,vstorm,istatus)

        I4_elapsed = ishow_timer()

        if(istat_radar_2dref.eq.1)then ! Get and advect the Reflectivity Stuff (LMR)
                                       ! Note: Ustorm and Vstorm need to be determined

            I4_elapsed = ishow_timer()

            if(n_fcst_radar .gt. 0)then
              index = n_fcst_radar
              do ifcst = 1,index
                write(6,*)
     1            ' Generating advected max reflectivity for pd ',ifcst       

                write(6,*)' Grid spacing (m) = ',grid_spacing_m

                call advect(ustorm,vstorm,ref_max(1,1,0)
     1                  ,dum1_2d,grid_spacing_m,imax,jmax
     1                  ,ref_max(1,1,ifcst)
     1                  ,float(ilaps_cycle_time*ifcst),1.0,lon
     1                  ,r_missing_data)

              enddo ! i
            endif

            I4_elapsed = ishow_timer()

!           Move writing of LMR product to here with appropriate changes
!           Write out LMR file
            write(6,*)' Writing out LMR file for max reflectivity'
            do ifcst = 0,n_fcst_radar
                minutes_10 = (ilaps_cycle_time * ifcst) / 600
                write(var_a(ifcst),101) 'R',minutes_10 
 101            format(a,i2.2)
                if(minutes_10 .lt. 10)then
!                   write(var_a(ifcst),101)minutes_10
!101                format('R0',i1)
cc                    var_a(ifcst) = 'R'
                    write(comment_a(ifcst),111)minutes_10
111                 format('LAPS Max Reflectivity  ',i1,'0 Min Fcst')
                else
cc                    write(var_a(ifcst),102)minutes_10
cc102                 format('R',i2)
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
     1       ' Not writing out Max Reflectivity LMR file'

        endif

!       Get and advect the Reflectivity History Stuff (LF1)
        if(.false.)then ! 'get_radar_max_pd' not yet modernized for radar input
            write(6,*)
            write(6,*)' Generating max reflectivity history analysis'
            call get_radar_max_pd(i4time_sys-ilaps_cycle_time
     1          ,i4time_sys,imax,jmax,kmax,heights_3d,ext_radar
     1          ,max_radar_files                                         ! I
     1          ,lat,lon,topo
     1          ,ref_max(1,1,0),frac_sum,istat_radar_hist)

        else
            istat_radar_hist = 0

        endif

        if(istat_radar_hist.eq.1)then

            I4_elapsed = ishow_timer()

            if(n_fcst_radar .gt. 0)then
              index = n_fcst_radar
              do ifcst = 1,index
                write(6,*)
     1          ' Generating advected max reflectivity history for pd'       
     1          ,ifcst

                call advect(ustorm,vstorm,ref_max(1,1,0)
     1                  ,dum1_2d,grid_spacing_m,imax,jmax
     1                  ,ref_max(1,1,ifcst)
     1         ,float(ilaps_cycle_time*ifcst),1.0,lon,r_missing_data)

              enddo ! ifcst
            endif

            I4_elapsed = ishow_timer()

!           Write out LF1 file
            write(6,*)
     1      ' Writing out LF1 file for max reflectivity history'

            do ifcst = 0,n_fcst_radar
                minutes_10 = (ilaps_cycle_time * ifcst) / 600
                write(var_a(ifcst),101) 'H',minutes_10
                if(minutes_10 .lt. 10)then
!                   write(var_a(ifcst),201)minutes_10
!201                format('H0',i1)
ccc                 var_a(ifcst) = 'H'
                    write(comment_a(ifcst),211)minutes_10
211                 format('LAPS Max Reflectivity History  ',i1
     1                    ,'0 Min Fcst')
                else
ccc                 write(var_a(ifcst),202)minutes_10
ccc202              format('H',i2)
                    write(comment_a(ifcst),212)minutes_10
212                 format('LAPS Max Reflectivity History ' ,i2
     1                    ,'0 Min Fcst')
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
            write(6,*)
     1      ' Not writing out Max Reflectivity History LF1 file'

        endif


        if(.true.)then ! Calculate Helicity
                       ! Note: (u+vstorm must be calculated first)

!           Write out Helicity field
            write(6,*)' Calculating Helicity'
            call helicity_laps(uanl,vanl,ustorm,vstorm
     1                        ,heights_3d,topo        
     1                        ,imax,jmax,kmax,helicity,istatus)

            if(istatus .eq. 1)then
                ext = 'lhe'

                var_a(0) = 'LHE'
                var_a(1) = 'MU'
                var_a(2) = 'MV'
                var_a(3) = 'SHU'
                var_a(4) = 'SHV'

                units_a(0) = 'm/s**2'
                units_a(1) = 'm/s'
                units_a(2) = 'm/s'
                units_a(3) = 'm/s'
                units_a(4) = 'm/s'

                comment_a(0) = 'LAPS Helicity'
                comment_a(1) = 'LAPS Mean Wind 0-6km AGL'
                comment_a(2) = 'LAPS Mean Wind 0-6km AGL'
                comment_a(3) = 'U Shear Component 0-6km AGL'
                comment_a(4) = 'V Shear Component 0-6km AGL'

                call move(helicity,out_lhe_3d(1,1,1),NX_L,NY_L)
                call move(umean   ,out_lhe_3d(1,1,2),NX_L,NY_L)
                call move(vmean   ,out_lhe_3d(1,1,3),NX_L,NY_L)
                call move(ushear  ,out_lhe_3d(1,1,4),NX_L,NY_L)
                call move(vshear  ,out_lhe_3d(1,1,5),NX_L,NY_L)

                call put_laps_multi_2d(i4time_sys,ext,var_a
     1          ,units_a,comment_a,out_lhe_3d,NX_L,NY_L,5,istatus)
                 istat_lhe = istatus

            else
                write(6,*)' WARNING: LHE not written out'

            endif ! istatus

        endif

        I4_elapsed = ishow_timer()

!       Calculate upslope component of moisture flux (PSD conventions)
        call get_grid_spacing_array(lat,lon,NX_L,NY_L,dx,dy)
	call up_mflux(NX_L,NY_L,NZ_L,topo,ldf,dx,dy
     1                     ,uanl,vanl,tpw_2d,upslope_flux
     1                     ,heights_3d,r_missing_data)

        if(.true.)then ! LIW

!           Calculate Li * 600mb Omega
            write(6,*)' Generating Li * Omega file'

!           Read in LI data
            var_2d = 'LI'
            ext = 'lst'
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
                liw = r_missing_data

            else ! Good Li Data
                call cpt_liw(lifted,wanl_2d,imax,jmax,liw)

            endif

            call move(liw    ,field_array(1,1,1),NX_L,NY_L)
!           call move(wanl_2d,field_array(1,1,2),NX_L,NY_L)
            call move(upslope_flux,field_array(1,1,2),NX_L,NY_L)        

!           Write out LIW field
!           Note that these arrays start off with 0 as the first index
            var_a(0) = 'LIW'
            var_a(1) = 'UMF'
            ext = 'liw'
            units_a(0) = 'K-Pa/s'
!           units_a(1) = 'Pa/s  '
            units_a(1) = 'M**2/s'
            comment_a(0) = 'Log LAPS Li * 600mb Omega'
!           comment_a(1) = 'LAPS 600mb Omega'
            comment_a(1) = 'Upslope component of moisture flux'

            call put_laps_multi_2d(i4time_sys,ext,var_a,units_a
     1                ,comment_a,field_array,imax,jmax,2,istatus)

            istat_liw = istatus

        endif

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
