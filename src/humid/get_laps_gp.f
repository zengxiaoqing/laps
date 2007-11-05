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
cdis cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
        subroutine get_laps_gp (rnum,rlat,rlon,ix,jy,glat,glon,
     1                           igrid,jgrid)

c       $log: get_laps_gp.for,v $
c revision 1.2  1995/09/13  21:35:19  birk
c added disclaimer to files
c
c revision 1.1  1994/04/25  15:20:18  birk
c initial revision
c


c this routine does 3 things.
c 1 it reads the laps grid  (this could be any one)
c 2 it determines the 4 grid indexes to use for the stations in question
c 3 it figures out the bilinear weights to apply to each of the 4 points
c modified read-laps-static to rd_laps_static 12/6/93 db

        implicit none


c input parameters

      integer rnum,igrid,jgrid
      real rlat(rnum),rlon(rnum)
      integer ix(rnum),jy(rnum)
      real glat(igrid,jgrid),glon(igrid,jgrid)
      character*31 static_ext
      character*200 static_dir
c regular internal variables



        integer i,j,k  !indexes
        integer len
        data static_ext /'nest7grid'/

        call get_directory('static',static_dir,len)


c   coordinates of the interpolation box

        do k = 1,rnum   ! do for all stations

        do j = 1,jgrid-1
        do i = 1,igrid-1

        if(
     1  glat(i,j).le.rlat(k) .and. glon(i,j).le.rlon(k)
     1  ) then

        if(
     1  glat(i+1,j+1).gt.rlat(k) .and. glon(i+1,j+1).gt.rlon(k)
     1  ) then ! we have define a box

        ix(k) = i
        jy(k) = j

        endif
        endif


        enddo !j
        enddo !i
        enddo !k

        return
        end



