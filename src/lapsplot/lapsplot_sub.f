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
        integer*4 SYS$TRNLOG,ASKI4T
        character*9 asc9_tim
        character*4 RM
        character*35 TIME
        logical l_atms

        character*2 c_section

        Common/RamTDims/NXMin,NXMax,NYMin,NYMax

        common /zoom/ zoom,density


        include 'lapsparms.cmn'

!       integer*4 nk
!       common /get_packed_data2/ nk

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

        call get_maxstns(maxstns,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'Error getting value of maxstns'
           stop
        endif

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

!       Supply reference time for matching the data
        lun = 5

        write(6,*)
        write(6,*)
        I4TIME_ref = ASKI4T()
        i4time_ref = (i4time_ref)/60*60

        if(.true.)then       ! testcode
!           call GOPKS(6,10000)
!           call GOPWK(1,51,3)
!           call GACWK(1)

            CALL GOPKS (6,IDUM)
            CALL GOPWK (1,2,1)
            CALL GACWK (1)

            write(6,*)' Calling COLOR'
            call color

            write(6,*)' Calling GOPWK'
            CALL GOPWK (9,1,3)

        else
            call OPNGKS
        endif

        call get_lapsplot_parms(namelist_parms,istatus)       

1100    write(6,1110)
1110    format(/////'     [h/hz]  Horizontal Plan View '
     1        ,'   (Const Pressure Level or SFC)'
     1        /'     [x]  Vertical Cross Section'
     1        /'     [s]  Sounding'
     1        /' ',60x,'[q] QUIT ? ',$)
        read(lun,1111)c_section
1111    format(a)

        if(c_section .eq. 'q' .or. c_section .eq. 'Q')goto999

        call get_laps_cycle_time(laps_cycle_time, istatus)
        if(istatus .ne. 1)then
            write(6,*)' Bad call to get_laps_cycle_time'
        endif

        IF(     c_section .eq. 'h' .or. c_section .eq. 'H' 
     1     .or. c_section .eq. '1'
     1     .or. c_section .eq. 'hz'
     1                                                          )THEN

            if(c_section .eq. 'hz')then
                write(6,101)
 101            format(
     1       '    Zoom,Density (contours/sfc wind barbs), Contour width'       
     1                ,15x,'     ? ',$)
                read(lun,*)zoom,density,plot_parms%contour_line_width       
!               plot_parms%contour_line_width = 1
            else
                zoom = 1.0
                density = 1.0
                plot_parms%contour_line_width = 1
            endif

            if(MAX_RADARS .ge. 1)then
                L_RADARS = 1
            else
                L_RADARS = 0
            endif

            call lapswind_plot(c_display,i4time_ref,lun,NX_L,NY_L,NZ_L,
     1                         MAX_RADARS,L_RADARS,r_missing_data,
     1                         laps_cycle_time,zoom,density,
     1                         plot_parms,namelist_parms)
            call frame

        elseif(c_section .eq. 'x' .or. c_section .eq. 'X'
     1                            .or. c_section .eq. 'xz'
     1                            .or. c_section .eq. '2')THEN
            zoom = 1.0

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

            write(6,*)' NXSECT = ',NXSECT

            call xsect(c_display,i4time_ref,lun,l_atms
     1                ,standard_longitude,NX_L,NY_L,NZ_L,121,NZ_L,NXSECT       
     1                ,r_missing_data,laps_cycle_time,maxstns
     1                ,density,plot_parms,namelist_parms)

        elseif(c_section .eq. 's' .or. c_section .eq. 'S')THEN
            call plot_sounding(i4time_ref,lun,NX_L,NY_L,NZ_L
     1                        ,r_missing_data,laps_cycle_time,maxstns
     1                        ,namelist_parms)       

        endif ! c_section

        write(6,901)
901     format(////20x,'   ----   ADD NEW FRAME??    ----')


        goto 1100

999     continue

        call sflush

        write(6,*)
        write(6,*)' closing GKS, look for your ./gmeta file'
!       call CLSGKS

        if(c_display .ne. 'v')then
!           write(6,*)' Saving metacode file in for008.dat'
        endif

        return

        end
