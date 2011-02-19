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
        function filename13(i4time,ext)

!       Steve Albers            1990

!       This routine constructs a filename from the logical DATE_TIME and
!       a passed in extension

        character ext*3,filename13*13,cdum13*13

        character*9 asc9_time

        common /laps_diag/ no_laps_diag

        call make_fnam_lp(i4time,asc9_time,istatus)

        cdum13 = asc9_time//'.'//ext

        call downcase(cdum13,cdum13)

        if(no_laps_diag .eq. 0)then
            write(6,*)' filename13 = ',cdum13
        endif

        filename13 = cdum13

        return

        end

        function filename14(i4time,ext)

!       Steve Albers            1990

!       This routine constructs a filename from the logical DATE_TIME and
!       a passed in extension

        character ext*4,filename14*14,cdum14*14

        character*9 asc9_time

        common /laps_diag/ no_laps_diag

        call make_fnam_lp(i4time,asc9_time,istatus)

        cdum14 = asc9_time//'.'//ext

        call downcase(cdum14,cdum14)

        if(no_laps_diag .eq. 0)then
            write(6,*)' filename14 = ',cdum14
        endif

        filename14 = cdum14

        return

        end
