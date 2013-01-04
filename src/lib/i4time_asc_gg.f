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
C
C Given a ASCII time such as 09-mar-1992 12:00:00.00, return the corresponding
C i4time (seconds since 1-1-1960). The use of the seconds field is a bit of 
c a historical relic. This function is mostly used to parse the types of dates 
c a user might type in.
C
        function        i4time_asc_gg(atime,istatus)

        implicit        none

        integer       i4time_asc_gg,  ! this function
     1          i4time_int_lp,  ! function call
     1          istatus,        ! argument (exported)
     1          year,
     1          month,
     1          day,
     1          hour,
     1          min,
     1          sec,
     1          i

        character*24    atime           ! argument (imported)
        character*4     ayear
        character*3     amonth(12)
        character*2     aday
        character*2     ahour
        character*2     amin
        character*2     asec

        data            amonth/'JAN','FEB','MAR','APR','MAY','JUN',
     1                 'JUL','AUG','SEP','OCT','NOV','DEC'/

        ayear=atime(8:11)
        aday=atime(1:2)
        ahour=atime(13:14)
        amin=atime(16:17)
        asec=atime(19:20)

        read(ayear,900)  year
        read(aday,901)  day
        read(ahour,901)  hour
        read(amin,901)  min
        read(asec,901)  sec
900     format(i4)
901     format(i2)

        do i=1,12
                if (atime(4:6) .eq. amonth(i)) then
                        month=i
                        goto 10
                endif
        enddo

 10     continue

        i4time_asc_gg=i4time_int_lp(year,month,day,hour,min,sec,istatus)

        return

        end
