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
c
c
	subroutine write_surface_obs(btime,outfile,n_meso,n_meso_pos,
     &    n_sao_g,n_sao_pos_g,n_sao_b,n_sao_pos_b,n_obs_g,n_obs_pos_g,
     &    n_obs_b,n_obs_pos_b,stations,store,wx,obstype,store_emv,
     &    store_amt,store_hgt,maxsta,num_var,badflag,jstatus)
c
c*****************************************************************************
c
c	Routine to write the LAPS surface data file.   The data is passed
c       to this routine via the 'store' array.
c
c	Changes:
c		P. Stamus  10-27-94  Original version (from get_surface_obs).
c
c*****************************************************************************
c
	real*4 store(maxsta,num_var), store_hgt(maxsta,5)
c
	integer*4 jstatus
c
	character btime*24,outfile*70, store_amt(maxsta,5)*4,
     &		store_emv(maxsta,5)*1,stations(maxsta)*3,
     &		wx(maxsta)*8,obstype(maxsta)*8
c
c
c.....	Write the file.
c
	open(11,file=outfile,status='unknown')
c
c.....	Write the header.
c
	write(11,900) btime,		! time
     &               n_meso,		! # of mesonet stations
     &               n_meso_pos,	! total # mesonet stations possible
     &               n_sao_g,		! # of saos in the laps grid
     &               n_sao_pos_g,	! total # of saos possible in laps grid
     &               n_sao_b,		! # of saos in the box
     &               n_sao_pos_b,	! total # of saos possible in the box
     &               n_obs_g,		! # of obs in the laps grid
     &               n_obs_pos_g,	! total # of obs psbl in the laps grid
     &               n_obs_b,		! # of obs in the box
     &               n_obs_pos_b 	! total # of obs possible in the box
900	format(1x,a24,10(1x,i4))
c
c.....	Write the station data.
c
	do k=1,n_obs_b
	  write(11,901) stations(k)(1:3),(store(k,i),i=1,3),obstype(k),
     &                 nint(store(k,4)),wx(k)
901	  format(1x,a3,1x,2(f7.2,1x),f5.0,1x,a8,1x,i4,1x,a8)
c
	  write(11,903) (store(k,i),i=5,13)
903	  format(4x,2(f6.1,1x),4(f5.0,1x),3(f6.1,1x))
c
	  kkk_s = nint(store(k,14))
	  write(11,905) kkk_s,(store(k,i),i=15,19),nint(store(k,20))
905	  format(4x,i2,2(1x,f7.1),1x,f5.1,1x,f7.3,1x,f6.1,1x,i4)
c
c.....	Write the cloud data if we have any.
c
	  if(kkk_s .gt. 0) then
	    do ii=1,kkk_s
	   write(11,907) store_emv(k,ii),store_amt(k,ii),store_hgt(k,ii)
907	      format(5x,a1,1x,a4,1x,f7.1)
	    enddo !ii
	  endif
c
	enddo !k
	endfile(11)
	close(11)	
c
c..... End of data writing.  Let's go home...
c
	return
	end
