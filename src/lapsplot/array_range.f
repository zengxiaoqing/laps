cdis   
cdis    Open Source License/Disclaimer, Forecast Systems Laboratory
cdis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
cdis    
cdis    This software is distributed under the Open Source Definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    In particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - Redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - Redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - All modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - If significant modifications or enhancements are made to this
cdis    software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
cdis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
cdis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
cdis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
cdis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
cdis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
cdis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
cdis   
cdis
cdis
cdis   
cdis


       subroutine contour_settings(a,ni,nj,clow,chigh,cint,zoom,scale)       

       real*4 a(ni,nj)

       call get_r_missing_data(r_missing_data,istatus)

       call array_range(a,ni,nj,rmin,rmax,r_missing_data)

       rmax = rmax / scale
       rmin = rmin / scale

       range = (rmax-rmin) / sqrt(zoom)

       if(range .gt. 2000)then
           cint = 400.
       elseif(range .gt. 600)then
           cint = 100.
       elseif(range .gt. 200)then
           cint = 50.
       elseif(range .gt. 60)then ! From 60-200, cint = 10 (6  - 20 contours)
           cint = 10.
       elseif(range .gt. 30)then ! From  30-60, cint = 5  (6  - 12 contours)
           cint = 5.
       elseif(range .gt. 6)then  ! From   6-30, cint = 2  (3  - 15 contours)
           cint = 2.
       elseif(range .ge. 1)then  ! From   1- 6, cint = 1  (1  -  6 contours)
           cint = 1.
       else ! range < 1          
           cint = 0.1
       endif

       cint_2 = cint * 2.

       if(cint .ge. 1.0)then
           clow  = int(rmin)/int(cint_2) * int(cint_2)
           chigh = int(rmax)/int(cint) * int(cint) + int(cint)
       else
           clow =  nint(rmin/cint) * cint - cint
           chigh = nint(rmax/cint) * cint + cint
       endif
 
       write(6,*)' Subroutine contour_settings....',range,zoom,cint,scale

       return
       end


       subroutine array_range(a,ni,nj,rmin,rmax,r_missing_data)

       real*4 a(ni,nj)

       rmin =  abs(r_missing_data)
       rmax = -abs(r_missing_data)

       do i = 1,ni
       do j = 1,nj
           if(a(i,j) .ne. r_missing_data)then
               rmin = min(rmin,a(i,j))
               rmax = max(rmax,a(i,j))
           endif
       enddo ! j
       enddo ! i

       return
       end
