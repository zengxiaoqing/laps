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
        subroutine degfrom(xd,yd,deg_from)

        real*4 deg_from,xd,yd,xx

             if(yd .eq. 0.0 .and. xd .eq. 0.0)then
             deg_from = 360.
             goto 1019
             endif

c   North of current storm.

             if(yd .eq. 0.0 .and. xd .gt. 0.0)then
             deg_from = 360.
             goto 1019
             endif

c   South of current storm.

             if(yd .eq. 0.0 .and. xd .lt. 0.0)then
             deg_from = 180.
             goto 1019
             endif

             xx = xd/yd
             deg_from = atand(xx)

c   Left half.

             if(yd .lt. 0.0)then
             deg_from = 270. - deg_from
             endif

c   Right half.

             if(yd .gt. 0.0)then
             deg_from = 90. - deg_from
             endif


1019     return
        end
