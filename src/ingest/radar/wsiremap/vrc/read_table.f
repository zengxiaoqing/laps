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
       subroutine read_table(table_path,nx,ny,lat,lon,ri,rj,istatus)

        integer nx,ny
        real lat (nx,ny)
        real lon (nx,ny)
        real ri (nx,ny)
        real rj (nx,ny)
        integer istatus
        character*(*) table_path

	istatus = 0

 	open(unit=12,file=table_path,form='unformatted',status='old'
     +       ,err=23)

	read(12,err=23,end=23) lat
	read(12,err=23,end=23) lon
	read(12,err=23,end=23) ri
	read(12,err=23,end=23) rj

	close (12)

	istatus = 1

	return
 23	print *,'Could not open vrc look up table',table_path
        stop 'read_table.f'
	end
