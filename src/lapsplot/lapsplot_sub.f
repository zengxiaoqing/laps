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

        subroutine lapsplot_3d()

        character*1 c_display
        integer*4 SYS$TRNLOG,ASKI4T
        character*9 asc9_tim
        character*4 RM
        character*35 TIME
        logical l_atms

        character*1 c_section

        Common/RamTDims/NXMin,NXMax,NYMin,NYMax


        include 'lapsparms.cmn'

!       integer*4 nk
!       common /get_packed_data2/ nk

!       nk = 17 ! A temporary fix as the blockdata won't work

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
     1    ' It will be located in ~/gmeta'
!       1'    For Ramteks other than RMA0, run trans on metacode file.'
!           write(6,*)
!           write(6,*)
!       1'    For display devices beyond FSL, copy over metacode file located'
!           write(6,*)
!       1'    in HOME$[STORMFEST]for008.dat; then translate onto local device.'
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

1100    write(6,1110)
1110    format(/////'     [h]  Horizontal Plan View '
     1        ,'   (Const Pressure Level or SFC)'
     1        /'     [x]  Vertical Cross Section'
     1        /' ',60x,'[q] QUIT ? '$)
        read(lun,1111)c_section
1111    format(a)

        if(c_section .eq. 'q' .or. c_section .eq. 'Q')goto999

        call get_laps_cycle_time(laps_cycle_time, istatus)
        if(istatus .ne. 1)then
            write(6,*)' Bad call to get_laps_cycle_time'
        endif

        IF(c_section .eq. 'h' .or. c_section .eq. 'H' .or. c_section .eq
     1. '1'
     1                                                          )THEN
            call lapswind_plot(c_display,i4time_ref,lun)
            call frame

        elseif(c_section .eq. 'x' .or. c_section .eq. 'X'
     1                                  .or. c_section .eq. '2')THEN
            l_atms = .false.
            call xsect(c_display,i4time_ref,lun,l_atms,standard_longitud
     1e)

        endif ! c_section

        write(6,901)
901     format(////20x,'   ----   ADD NEW FRAME??    ----')


        goto 1100

999     continue

        write(6,*)
        write(6,*)' closing GKS, look for your ~/gmeta file'
        call CLSGKS

        if(c_display .ne. 'v')then
!           write(6,*)' Saving metacode file in for008.dat'
        endif

        return

        end

