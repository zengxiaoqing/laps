
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


        function a8_to_a9(a8_time)

        character*9 a8_to_a9
        character*9 a9
        character*8 a8_time

        integer*4 imon_a(12)

        data imon_a/0,31,59,90,120,151,181,212,243,273,304,334/

        read(a8_time,1)iyr,imn,idy,ih
1       format(i2,i2,i2,i2)

        id = imon_a(imn)

        idays = id + idy

!       Decide whether to add a day for leap year.
        if(iyr .eq. (iyr / 4) * 4 )then
            if(imn .ge. 3)then
                idays = idays + 1
            endif
        endif

        write(a9,2)iyr,idays,ih,im
2       format(i2,i3,i2,i2)

        if(a9(3:3) .eq. ' ')a9(3:3) = '0'
        if(a9(4:4) .eq. ' ')a9(4:4) = '0'
        if(a9(6:6) .eq. ' ')a9(6:6) = '0'
        if(a9(8:8) .eq. ' ')a9(8:8) = '0'

        a8_to_a9 = a9

        return
        end


