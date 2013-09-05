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

        subroutine lapsplot_3d()

        include 'lapsplot.inc'

        character*1 c_display
        integer ASKI4T
        character*9 asc9_tim
        character*4 RM
        character*35 TIME
        character*4 czoom
        logical l_atms, l_plotobs
        logical l_parse, l_polar, l_cyl

        character*2 c_section

        Common/RamTDims/NXMin,NXMax,NYMin,NYMax

        common /zoom/ zoom,density

!       integer nk
!       common /get_packed_data2/ nk

        real, allocatable :: dx(:,:)            
        real, allocatable :: dy(:,:)            

!       nk = 17 ! A temporary fix as the blockdata won't work

        call get_grid_dim_xy(NX_L,NY_L,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'Error getting horizontal domain dimensions'
           stop
        endif

        call get_laps_dimensions(NZ_L,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'Error getting vertical domain dimension'
           stop
        endif

        call get_max_radars(MAX_RADARS,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'Error getting value of MAX_RADARS'
           stop
        endif

        call get_r_missing_data(r_missing_data,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'Error getting value of r_missing_data'
           stop
        endif

        call get_standard_longitude(standard_longitude,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'Error getting value of standard_longitude'
           stop
        endif

        call get_maxstns(maxstns,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'Error getting value of maxstns'
           stop
        endif

        dyn_low = r_missing_data
        dyn_high = r_missing_data

        allocate (dx(NX_L,NY_L),dy(NX_L,NY_L))

!       read(5,2)c_display
2       format(a1)
        c_display = 'r'

        if(c_display .eq. 'r')then
            write(6,*)
            write(6,*)
            write(6,*)' This program will produce a gmeta file for displ
     1ay'
            write(6,*)
            write(6,*)
     1    ' It will be located in ./gmeta'
        endif

!       Read namelist parameters
        call get_lapsplot_parms(namelist_parms,istatus)       

!       Set Default plot parameters
        plot_parms%rimage_intensity = 1.0
        plot_parms%zoom = 1.0
        plot_parms%xcen = 0.5
        plot_parms%ycen = 0.5
        plot_parms%zoom_wdw = 1.0
        plot_parms%obs_size = 1.0
        plot_parms%ncols = 0.
        plot_parms%icol_barbs = namelist_parms%icol_barbs
        plot_parms%l_hinterp_zoom = .false.

        if(.true.)then       ! testcode
!           call GOPKS(6,10000)
!           call GOPWK(1,51,3)
!           call GACWK(1)

            CALL GOPKS (6,IDUM)
            CALL GOPWK (1,2,1)
            CALL GACWK (1)

            write(6,*)' Calling COLOR'
            iwhite = namelist_parms%i_background_color
            call color(iwhite)

            write(6,*)' Calling GOPWK'
            CALL GOPWK (9,1,3)

        else
            call OPNGKS
        endif

!       Supply reference time for matching the data
        lun = 5
1090    write(6,*)
        write(6,*)
        I4TIME_ref = ASKI4T()
        i4time_ref = (i4time_ref)/60*60

1100    write(6,1110)
1110    format(/////'     [h/hz]  Horizontal Plan View '
     1        ,'   (Const Pressure Level or SFC)'
     1        /'     [x]  Vertical Cross Section'
     1        /'     [s]  Sounding'
     1        /'     [nt]  New Time'
     1        /' ',60x,'[q] QUIT ? ',$)
        read(lun,1111)c_section
1111    format(a)

        if(c_section .eq. 'nt')then
!           write(6,*)' Setting namelist_parms%iraster to -1'
!           write(6,*)' Raster may not work with multiple frames'
!           namelist_parms%iraster = -1
            goto1090
        endif

        if(c_section .eq. 'q' .or. c_section .eq. 'Q')goto999

        call get_laps_cycle_time(laps_cycle_time, istatus)
        if(istatus .ne. 1)then
            write(6,*)' Bad call to get_laps_cycle_time'
        endif

        if(c_section .eq. 'in')then
            write(6,*)' Input image intensity...'
            read(lun,*)plot_parms%rimage_intensity
            go to 1100
        endif

        ifield_found = 1
        i_overlay = 0

        IF(     c_section .eq. 'h' .or. c_section .eq. 'H' 
     1     .or. c_section .eq. '1'
     1     .or. c_section .eq. 'hz'
     1                                                          )THEN

            if(c_section .eq. 'hz')then
                write(6,101)
 101            format(
     1       '    Zoom,Density (contours/sfc wind barbs), Contour width'       
     1                ,15x,'     ? ',$)
 
                read(lun,*)zoom,density
     1                    ,plot_parms%contour_line_width       
     1                    ,plot_parms%xcen
     1                    ,plot_parms%ycen
     1                    ,plot_parms%zoom_wdw
     1                    ,plot_parms%obs_size

            else
                zoom = 1.0
                density = 1.0
                plot_parms%contour_line_width = 1
            endif

            plot_parms%zoom = zoom

            if(MAX_RADARS .ge. 1)then
                L_RADARS = 1
            else
                L_RADARS = 0
            endif

            call lapswind_plot(c_display,i4time_ref,lun,NX_L,NY_L,NZ_L,
     1                         MAX_RADARS,L_RADARS,r_missing_data,
     1                         laps_cycle_time,zoom,density,
     1                         dyn_low,dyn_high,dx,dy,
     1                         plot_parms,namelist_parms,ifield_found)       

            if(ifield_found .eq. 1)then
                write(6,*)' Field found - calling frame'
                call frame
            else
                write(6,*)' Field not found - not calling frame'
            endif

        elseif(c_section .eq. 'x' .or. c_section .eq. 'X'
     1                            .or. c_section .eq. 'xz'
     1                            .or. c_section .eq. '2')THEN
            zoom = 1.0
            plot_parms%zoom = zoom

            if(c_section .eq. 'xz')then
                write(6,102)
 102            format('    Density (contours), Contour Line Width'       
     1                ,13x,'     ? ',$)
                read(lun,*)density,plot_parms%contour_line_width
!               plot_parms%contour_line_width = 1
            else
                density = 1.0
                plot_parms%contour_line_width = 1
            endif

            l_atms = .false.

            NXSECT = nint(sqrt(float(NX_L**2 + NY_L**2))) ! 541

            write(6,*)' NXSECT (NX_P) = ',NXSECT

            NX_C = max(NX_L,NY_L)

            NX_T = nint(121. * density)                      

            write(6,*)' NX_C/NX_T= ',NX_C,NX_T

            call xsect(c_display,i4time_ref,lun,l_atms
     1                ,standard_longitude,NX_L,NY_L,NZ_L
     1                ,NX_C,NZ_L,NXSECT,NX_T      
     1                ,r_missing_data,laps_cycle_time,maxstns
     1                ,dyn_low,dyn_high,dx,dy
     1                ,density,plot_parms,namelist_parms,ifield_found)       

            if(ifield_found .eq. 1)then
                write(6,*)' Field found - calling frame'
                call frame
            else
                write(6,*)' Field not found - not calling frame'
            endif

        elseif(c_section .eq. 's' .or. c_section .eq. 'sz')THEN
            if(c_section .eq. 'sz')then
                read(lun,*)plot_parms%obs_size           
            endif
            
            call plot_sounding(i4time_ref,lun,NX_L,NY_L,NZ_L
     1                        ,r_missing_data,laps_cycle_time,maxstns
     1                        ,i_overlay,plot_parms,namelist_parms
     1                        ,l_plotobs)       

            if(ifield_found .eq. 1)then
                write(6,*)' Field found - calling frame'
                call frame
            endif

        elseif(c_section .eq. 'a' .or. c_section .eq. 'az')THEN
            if(c_section .eq. 'az')then
                write(6,*)' Enter plot info (e.g. 180p, 90c)'
                read(lun,*)czoom                   
                l_polar = l_parse(czoom,'p')
                l_cyl   = l_parse(czoom,'c')
                if((l_polar .OR. l_cyl) .eqv. .false.)then
                    write(6,*)' Set l_polar to TRUE'
                    l_polar = .true.
                endif
                write(6,*)' l_polar = ',l_polar,' l_cyl = ',l_cyl
!               l_polar = .true.
!               l_cyl = .true.
                write(6,*)' Enter minalt,maxalt (e.g. 0,90 0,180)'
                read(lun,*)minalt,maxalt           
            endif

            minazi = 0
            maxazi = maxalt * 4

            call plot_allsky(i4time_ref,lun,NX_L,NY_L,NZ_L
     1                        ,minalt,maxalt,minazi,maxazi
     1                        ,r_missing_data,laps_cycle_time,maxstns
     1                        ,i_overlay,plot_parms,namelist_parms
     1                        ,l_polar,l_cyl)       

            if(ifield_found .eq. 1)then
                write(6,*)' Field found - calling frame'
                call frame
            endif

        endif ! c_section

        write(6,901)
901     format(////20x,'   ----   ADD NEW FRAME??    ----')

        goto 1100 ! loop back for user input

999     continue ! quit/end program

        call sflush

        write(6,*)
        write(6,*)' closing GKS, look for your ./gmeta file'
!       call CLSGKS

        if(c_display .ne. 'v')then
!           write(6,*)' Saving metacode file in for008.dat'
        endif

        deallocate(dx,dy)

        return

        end
