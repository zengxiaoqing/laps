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
	subroutine write_table (table_path,
     &             nx,ny,lat,lon,ri,rj,istatus)

	integer nx,ny
	real lat (nx,ny)
	real lon (nx,ny)
	real ri (nx,ny)
	real rj (nx,ny)
	integer istatus
        character*(*) table_path

	istatus = 0
        n=index(table_path,' ')

	open (unit=12,file=table_path(1:n-1),
     &form='unformatted',status='unknown',err=23)

	write(12,err=24) lat
	write(12,err=24) lon
	write(12,err=24) ri
	write(12,err=24) rj

	close (12)

	istatus = 1

        goto 100

23      write(6,*)'Error openning ',table_path(1:n-1)
        goto 100

24      write(6,*)'Error writing to ',table_path(1:n-1)

100     return
        end
