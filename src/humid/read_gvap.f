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
      subroutine read_gvap (filename, nstations, path_to_gvap8,
     1     path_to_gvap10,time_diff,
     1     lat,lon, wt,w1,w2,w3,
     1     nn, istatus)

c     routine to read goes watervapor for system testing
c     author Dan Birkenheuer
c     2/8/99
c     
c     paramter list
c     filename, name of the wisconsin file
c     nstations, dimension of variables
c     lat, latitude of each station
c     lon, longitude of each station
c     wt, total water
c     w1, layer 1 water
c     w2, layer 2 water
c     w3, layer 3 water
c     nn, number of stations with real data for consideration
c
      implicit none

c     input variables
      character*9 filename
      integer nstations,nn,istatus,idummy,idummy2
      character*256 path_to_gvap8,path_to_gvap10
      integer time_diff         !time allowed for latency (sec)
      real lat(nstations)
      real lon(nstations)
      real wt(nstations)
      real w1(nstations)
      real w2(nstations)
      real w3(nstations)

c     internal variables

      integer i
      real dummy
      integer ptg_index
      character*512 const_file
      integer cf
      character*9 filefound
      character*120 extension
      integer extension_index
      integer istatus_8,istatus_10
      integer file_name_length
      character*120 desired_ext
      integer de_index


      file_name_length = 9
      desired_ext = 'tpw'
      istatus_8=0
      istatus_10 = 0
      call s_len(path_to_gvap8, ptg_index, istatus)

c     reading goes 8

c     get most recent file in directory

      write (6,*)
      write (6,*) 'Acquiring GOES  8 GVAP info'
      write (6,*)



      call get_newest_file (filename, time_diff,file_name_length,
     1     path_to_gvap8,ptg_index,filefound,desired_ext, de_index,
     1     extension, extension_index, istatus)

      if(istatus.ne.1) then     !failure in getting file
         return
      endif

      const_file = path_to_gvap8(1:ptg_index)//'20'//filefound
     1     //'.'//extension(1:extension_index)

      call s_len(const_file, cf, istatus)

      write(6,*) 'opening file ',const_file(1:cf)
      
      open(22, file=const_file(1:cf),form='formatted',
     1     status='old',err = 668)
      read(22,*,end=668,err=668) ! first header line is ignored
      do i = 1,nstations
         read(22,*,end=665,err=665) idummy,idummy,lat(i),lon(i),
     1        idummy,idummy,wt(i), w1(i),w2(i),w3(i)

      enddo

 665  close (22)
      nn = i-1
      if (nn .eq. 0) go to 666
      write(6,*) nn, ' number of records read GOES 8'
      istatus_8 = 1             !assign success to reading goes 8
c     write(6,*) (wt(i),i=1,nn)

      go to 669

 668  close (22)
      write(6,*) 'failed reading GOES 8'
      istatus_8 = 0
      nn = 1

 669  continue

c     reading goes 10

c     get most recent file in directory

      write (6,*)
      write (6,*) 'Acquiring GOES 10 GVAP info'
      write (6,*)

      istatus_10 = 1

      call s_len(path_to_gvap10, ptg_index, istatus)

      desired_ext = 'tpw'

      call get_newest_file (filename, time_diff,file_name_length,
     1     path_to_gvap10,ptg_index,filefound,desired_ext, de_index,
     1     extension, extension_index, istatus)

      if(istatus.ne.1) then     !failure in getting file
         istatus_10 = 0
         go to 667
      endif

      const_file = path_to_gvap10(1:ptg_index)//'20'//filefound
     1     //'.'//extension(1:extension_index)

      call s_len(const_file, cf, istatus)

      write(6,*) 'opening file ',const_file(1:cf)
      
      open(23, file=const_file(1:cf),form='formatted',
     1     status='old',err = 666)
      read(23,*,end=666,err=666) ! first header line is ignored
      do i = 1,nstations
         read(23,*,end=667,err=667) idummy,idummy,lat(i),lon(i),
     1        idummy,idummy,wt(i), w1(i),w2(i),w3(i)

      enddo


 667  close(23)

      if (istatus_8.eq.0) then
         write(6,*) 'GOES 8 failed'
      endif
      if (istatus_10.eq.0) then
         write(6,*) 'GOES 10 failed'
      endif

      if(istatus_8+istatus_10.eq.0) then
         return
      else
         write(6,*) 'got something,  processing....'
         istatus = 1
      endif
      
      nn = i-1 + nn
      if (nn .eq. 0) go to 666
      write(6,*) i-1, ' number of records read GOES 10'
      istatus = 1
c      write(6,*) (wt(i),i=1,nn)

      write (6,*) 'total GVAP records read is ', nn

      return

 666  if(nn.eq.0) then
         istatus = 0
         write(6,*) 'process failure...  no available gvap data'
      else
         istatus = 1
      endif

      return

      end

