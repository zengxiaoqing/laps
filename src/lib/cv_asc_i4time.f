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

!        character*9 ascii_time
!        integer*4 i4time

!       read(5,1)ascii_time
!1      format(a9)
!       call cv_asc_i4time(ascii_time,I4time)

!        write(6,*)i4time
!       end

        subroutine cv_asc_i4time(ascii_time,I4time)

        character*9 ascii_time
        integer*4 i4time

        read(ascii_time(1:2),2)iyear
        read(ascii_time(3:5),3)idoy
        read(ascii_time(6:7),2)ihour
        read(ascii_time(8:9),2)imin

2       format(i2)
3       format(i3)

!       Valid for years 1960-2060
        if(iyear .lt. 60)iyear = iyear + 100

        lp = (iyear + 3 - 60) / 4

        i4time =  (iyear-60) * 31536000
     1  + (idoy-1)   * 86400
     1  + ihour      * 3600
     1  + imin       * 60

        i4time = i4time + 86400 * lp

        return
        end
