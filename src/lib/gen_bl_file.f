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
      
      open (unit=22, file=pbfile, status='replace')
      
      write (22,5) filename,i4time
      write (22,6) pb
 5    format (1x,a9,i12)
 6    format (1x, 7f10.3)
      
      close (22)
      
      istatus = 1
      
      return
      end
