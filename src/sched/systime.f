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


	character*9 asc_time
        character*24 asc_tim_24
        character*2 utc_hour,utc_min
        character*5 utc_date
        character*200 fname 

        call get_laps_config('nest7grid',istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error calling get_laps_config'
            stop
        endif

        call get_laps_cycle_time(ilaps_cycle_time,istatus)
        if(istatus .eq. 1)then
            write(6,*)' ilaps_cycle_time = ',ilaps_cycle_time
        else
            write(6,*)' Error getting laps_cycle_time'
            stop
        endif

!       Write out systime.dat file that has seconds since 1960

        call get_directory('time',fname,len_fname)
        open(11,file=fname(1:len_fname)//'systime.dat',status='unknown')

	i4time = i4time_now_gg()
        i4time_sys = (i4time - 0) / ilaps_cycle_time * ilaps_cycle_time        

	call make_fnam_lp(i4time_sys,asc_time,istatus)

	write(11,2)i4time_sys
2	format(1x,i11)
	write(11,1)asc_time
1	format(1x,a9)

        utc_hour = asc_time(6:7)
        utc_min  = asc_time(8:9)
        utc_date = asc_time(1:5)

3       format(a2)	

        write(11,3)utc_hour

        write(11,3)utc_min

        call cv_i4tim_asc_lp(i4time_sys,asc_tim_24,istatus)
        write(11,5)asc_tim_24(1:14),asc_tim_24(16:17)
 5	format(a14,a2)

        write(11,10)utc_date
 10	format(a5)

        close(11)

!       Write out c_time.dat file that has seconds since 1970

        open(11,file=fname(1:len_fname)//'c_time.dat',status='unknown')

	write(11,101)asc_time
101	format(1x,a9)
	write(11,102)i4time_sys -315619200
102	format(1x,i11)

        close(11)

	end
