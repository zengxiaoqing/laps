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
      subroutine gen_bl_file (i4time,pb,igrid,jgrid,istatus)

c     this routine generates the boundary layer output file determined
c     for the sh mixing algorithm and subsequently for the retroactive
c     inclusion in the surface analysis.
      
c     inputs:
      
c     i4time
c     pb ! top of boundary layer pressure (mb)
      
c     outputs:   none!  this routine generates a file only where
c     filename (ascii time stamp for first record)
c     pb in a simple file (binary write)
      
      implicit none
      
      integer igrid,jgrid
      real pb(igrid,jgrid)
      integer i4time
      integer istatus, len
      character*9 filename
      character*200 pbfile

c     ----- begin exe  last modified 10/26/99 db
      
c     set istatus failure mode
      istatus = 0

      call get_directory('lpbl',pbfile,len)
      pbfile = pbfile(1:len)//'pbl_depth.lpbl'
      
      call make_fnam_lp (i4time,filename,istatus)
      
      if (istatus.eq.1) then
         istatus = 0
      else
         print*, 'Error in routine make_fnam_lp'
         istatus = 0
         return
      endif
      
      open (unit=22, file=pbfile, form ='formatted')
      
      write (22,5) filename,i4time
      write (22,6) pb
 5    format (1x,a9,i12)
 6    format (1x, 7f10.3)
      
      close (22)
      
      istatus = 1
      
      return
      end
