cdis    forecast systems laboratory
cdis    noaa/oar/erl/fsl
cdis    325 broadway
cdis    boulder, co     80303
cdis
cdis    forecast research division
cdis    local analysis and prediction branch
cdis    laps
cdis
cdis    this software and its documentation are in the public domain and
cdis    are furnished "as is."  the united states government, its
cdis    instrumentalities, officers, employees, and agents make no
cdis    warranty, express or implied, as to the usefulness of the software
cdis    and documentation for any purpose.  they assume no responsibility
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis
cdis    permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  all modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making
cdis    the modifications.  if significant modifications or enhancements
cdis    are made to this software, the fsl software policy manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
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
      real*4 glat(igrid,jgrid),glon(igrid,jgrid)
      character*31 static_ext
      character*200 static_dir
c regular internal variables



        integer*4 i,j,k  !indexes
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



